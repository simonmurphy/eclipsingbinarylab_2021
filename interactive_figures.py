import matplotlib.pyplot as plt
import ipywidgets as widgets
import numpy as np
from scipy import optimize
import isochrones
from matplotlib.offsetbox import AnchoredText 

def phase(t,t0,P):
    """ Calculate the orbital phase for a given time and period """
    return np.array((t - t0)/P - np.floor((t - t0)/P))

def solve_kepler_eqn(M,e):
    """ Solve Keplers equation M = E - e*sin(E) for E """
    try:
        M[0]
        res = np.zeros(M.shape)
        for i,Mi in enumerate(M):
            tmp,= optimize.fsolve(lambda x: x-e*np.sin(x) - Mi,Mi)
            res[i] = tmp
    except IndexError:
        res, = optimize.fsolve(lambda x: x - e*np.sin(x)-M,M)
    return res

def rv_keplerian(t,t0,P,k,w,e,vsys):
    """ Generate Keplerian RV for elliptical orbit at time <t> """
    # Simplification for circular orbits:
    if e == 0.0:
        f = (2*np.pi/P)*(t-t0)
        return k*np.cos(f + np.pi/2) + vsys
    # Otherwise solve Kepler's equation
    E = solve_kepler_eqn((2*np.pi/P)*(t-t0),e)
    f = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2),np.sqrt(1-e)*np.cos(E/2))
    return k*(np.cos(f + np.radians(w) + np.pi/2) + e*np.cos(np.radians(w))) + vsys

def string_length(phase, flux):
    """ Calculate the point-to-point distance of the phased light curve """
    idx = np.argsort(phase)
    phase, flux = phase[idx], flux[idx]
    dp = np.diff(phase,append=phase[0]-phase[-1]+1)
    df = np.diff(flux,append=flux[0]-flux[-1])
    return np.sum(np.sqrt(dp**2 + df**2))

def chi2(obs, calc, err, dof=None):
    chisq = np.sum((obs-calc)**2/err**2)
    if not dof:
        return chisq
    return chisq / dof

class RVCurve(widgets.HBox):

    def __init__(self, t1, rv1, e_rv1, t2, rv2, e_rv2):
        super().__init__()
        output = widgets.Output()

        self.t1 = t1
        self.t2 = t2
        self.rv1 = rv1
        self.rv2 = rv2
        self.e_rv1 = e_rv1
        self.e_rv2 = e_rv2

        # Initital parameters
        self._t0 = 8101.0
        self._P = 0.8589
        self._k1 = 50.0
        self._k2 = 50.0
        self._vsys = 0.0
        self._e = 0.0
        self._w = 0.0

        # Intial phases
        phase1 = phase(t1,self._t0,self._P)
        phase2 = phase(t2,self._t0,self._P)

        # Continuous phases for model curves
        self._phase = np.linspace(0, 1, 100)
        self._t = self._phase * 2*np.pi

        # Draw data points and initial fits
        with output:
            self.fig, self.ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

        self.obs1 = plt.errorbar(phase1,rv1,e_rv1,fmt='o',color='tab:red',ms=4,ecolor='#999999')
        self.obs2 = plt.errorbar(phase2,rv2,e_rv2,fmt='o',color='tab:blue',ms=4,ecolor='#999999')
        model = rv_keplerian(self._t,0,2*np.pi,self._k1,self._w,self._e,self._vsys)
        self.model1, = plt.plot(self._phase, model, '-', color='tab:red',alpha=0.5,lw=2)
        model = rv_keplerian(self._t,0,2*np.pi,-self._k2,self._w,self._e,self._vsys)
        self.model2, = plt.plot(self._phase, model, '-', color='tab:blue',alpha=0.5,lw=2)
        self.sys, = plt.plot(self._phase, 0*self._phase+self._vsys,lw=2,alpha=0.5,color='tab:green')
        plt.xlim(0,1)
        plt.xlabel('Phase')
        plt.ylabel('Radial velocity (km s$^{-1}$)')
        plt.title('WiFeS radial velocities')
        self.fig.canvas.toolbar_position = 'left'
        self.fig.set_label(' ')

        # Calculate initial chi-squared for data points
        model = rv_keplerian(t1,self._t0,self._P,self._k1,self._w,self._e,self._vsys)
        chisq = chi2(rv1,model,e_rv1,dof=len(rv1)-6)
        self.chisq1 = AnchoredText('$\chi^{2}_{r}$ = %3.2f'%chisq,loc='upper right',prop=dict(color='w'))
        self.chisq1.patch.set(facecolor='tab:red',alpha=0.5)
        self.ax.add_artist(self.chisq1)
        model = rv_keplerian(t2,self._t0,self._P,-self._k2,self._w,self._e,self._vsys)
        chisq = chi2(rv2,model,e_rv2,dof=len(rv2)-6)
        self.chisq2 = AnchoredText('$\chi^{2}_{r}$ = %3.2f'%chisq,loc='upper left',prop=dict(color='w'))
        self.chisq2.patch.set(facecolor='tab:blue',alpha=0.5)
        self.ax.add_artist(self.chisq2)

        # Define widgets
        self.P = widgets.FloatText(value=self._P,description='Period (d)',step=1e-5,style={'description_width': '10em'})
        self.t0 = widgets.FloatText(value=self._t0,description='$t_{0}$ (BJD $-$ 2450000)',step=1e-3,style={'description_width': '10em'})
        self.k1 = widgets.IntSlider(min=0,max=100,step=1,value=self._k1,description='$K_{1}$')
        self.k2 = widgets.IntSlider(min=0,max=200,step=1,value=self._k2,description='$K_{2}$')
        self.vsys = widgets.IntSlider(min=-50,max=50,step=1,value=self._vsys,description='$v_{sys}$')
        self.e = widgets.FloatSlider(min=0,max=0.99,step=0.01,value=self._e,description='$e$')
        self.w = widgets.IntSlider(min=0,max=360,step=1,value=self._w,description='$\omega$')
        self.show_models = widgets.Checkbox(description='Show models',value=True)
        self.show_vsys = widgets.Checkbox(description='Show $v_{sys}$',value=True)
        self.show_grid = widgets.Checkbox(description='Show grid',value=False)

        # Monitor for updates
        self.P.observe(self.update_points,'value')
        self.t0.observe(self.update_points,'value')
        self.k1.observe(self.update_models,'value')
        self.k2.observe(self.update_models,'value')
        self.vsys.observe(self.update_vsys,'value')
        self.e.observe(self.update_models,'value')
        self.w.observe(self.update_models,'value')
        self.show_models.observe(self.update_show_models,'value')
        self.show_vsys.observe(self.update_show_vsys,'value')
        self.show_grid.observe(self.update_show_grid,'value')

        controls = widgets.VBox([self.P,self.t0,self.k1,self.k2,self.vsys,self.e,self.w,\
            self.show_models,self.show_vsys,self.show_grid])

        # Add to children
        self.children = [output,controls]

    def update_chi2(self):
        model = rv_keplerian(self.t1,self.t0.value,self.P.value,self.k1.value,self.w.value,self.e.value,self.vsys.value)
        chi = chi2(self.rv1,model,self.e_rv1,dof=len(self.rv1)-6)
        self.chisq1.txt.set_text('$\chi^{2}_{r}$ = %3.2f'%chi)
        model = rv_keplerian(self.t2,self.t0.value,self.P.value,-self.k2.value,self.w.value,self.e.value,self.vsys.value)
        chi = chi2(self.rv2,model,self.e_rv2,dof=len(self.rv2)-6)
        self.chisq2.txt.set_text('$\chi^{2}_{r}$ = %3.2f'%chi)

    def update_points(self, change):
        """ Some trickery to move the error bars as well as the points """
        # Primary velocities
        phi = phase(self.t1,self.t0.value,self.P.value)
        ln, caps, bars = self.obs1
        ln.set_xdata(phi)
        bars[0].set_segments([[[x,yt], [x,yb]] for x, yt, yb in zip(phi, self.rv1+self.e_rv1, self.rv1-self.e_rv1)])
        # Secondary velocities
        phi = phase(self.t2,self.t0.value,self.P.value)
        ln, caps, bars = self.obs2
        ln.set_xdata(phi)
        bars[0].set_segments([[[x,yt], [x,yb]] for x, yt, yb in zip(phi, self.rv2+self.e_rv2, self.rv2-self.e_rv2)])
        # Update the chisq value
        self.update_chi2()

    def update_models(self, change):
        model = rv_keplerian(self._t,0,2*np.pi,self.k1.value,self.w.value,self.e.value,self.vsys.value)
        self.model1.set_ydata(model)
        model = rv_keplerian(self._t,0,2*np.pi,-self.k2.value,self.w.value,self.e.value,self.vsys.value)
        self.model2.set_ydata(model)
        self.update_chi2()

    def update_vsys(self, change):
        self.update_models(change)
        self.sys.set_ydata(self._t*0 + change.new)

    def update_show_models(self, change):
        for i in (self.model1,self.model2,self.chisq1,self.chisq2):
            i.set_visible(change.new)

    def update_show_vsys(self, change):
        self.sys.set_visible(change.new)

    def update_show_grid(self,change):
        self.ax.grid(change.new)


class LightCurve(widgets.HBox):
        
    def __init__(self, t, flux, P=1.0, t0=8468.287, dP=1e-3):
        super().__init__()
        output = widgets.Output()

        self.t = t
        self.flux = flux
        
        self._t0 = t0
        self._P = P

        # Display phases between -0.25 and 0.75
        phi = phase(t,self._t0,self._P)
        phi[phi>0.75] -= 1

        # Draw data points and initial fits
        with output:
            self.fig, self.ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

        # Sort the phases so we can draw the lines
        idx = np.argsort(phi)
        self.lc, = plt.plot(phi[idx], self.flux[idx], 'o-', ms=2, color='#AAAAAA',mfc='k',mec='k')
        plt.xlim(-0.25,0.75)
        plt.xlabel('Phase')
        plt.ylabel('Flux (normalised)')
        plt.title('$TESS$ phased light curve')
        self.fig.canvas.toolbar_position = 'left'
        self.fig.set_label(' ')
        # Calculate initial string length
        str_len = string_length(phi,flux)
        self.str_len = AnchoredText('String length = %.2f'%str_len, loc='lower right',prop=dict(color='w'))
        self.str_len.patch.set(facecolor='tab:red',alpha=0.5)
        self.ax.add_artist(self.str_len)
        self.lc.set_linestyle('None') # make invisible initially
        self.str_len.set_visible(False)

        # Define widgets
        self.P = widgets.FloatText(value=self._P,description='Period (d)',step=dP,style={'description_width': '10em'})
        self.t0 = widgets.FloatText(value=self._t0,min=self._t0-1,max=self._t0+1,description='$t_{0}$ (BJD $-$ 2450000)',step=1e-3,style={'description_width': '10em'})
        self.show_string = widgets.Checkbox(description='Show string', value=False,style={'description_width': '10em'})
        self.show_grid = widgets.Checkbox(description='Show grid',value=False,style={'description_width': '10em'})

        # Monitor for updates
        self.P.observe(self.update_points,'value')
        self.t0.observe(self.update_points,'value')
        self.show_string.observe(self.update_show_string,'value')
        self.show_grid.observe(self.update_show_grid,'value')

        controls = widgets.VBox([self.P,self.t0,self.show_string,self.show_grid])

        # Add to children
        self.children = [output,controls]        

    def update_points(self, change):
        phi = phase(self.t,self.t0.value,self.P.value)
        phi[phi>0.75] -= 1
        idx = np.argsort(phi)
        self.lc.set_xdata(phi[idx])
        self.lc.set_ydata(self.flux[idx])
        # Update the string length 
        self.str_len.txt.set_text('String length = %.2f'%string_length(phi, self.flux))
    
    def update_show_string(self, change):
        self.str_len.set_visible(change.new)
        ls = '-' if change.new else 'None'
        self.lc.set_linestyle(ls)

    def update_show_grid(self, change):
        self.ax.grid(change.new)


class MassRadiusDiagram(widgets.HBox):
    
    def __init__(self):
        super().__init__()
        output = widgets.Output()

        self.mass1 = 1 # Msun
        self.mass2 = 1 # Msun
        self.radius1 = 1 # Rsun
        self.radius2 = 1 # Rsun
        self.age = 100.0 # Myr
        self.mass_range = [0.015,0.8] # for isochrone plots

        # Draw data points and initial fits
        with output:
            self.fig, self.ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

        # Plot a few isochrones for guidance
        iso_ages = [1,2,3,5,10,20,50]
        mass,radius,ages = isochrones.BHAC15(iso_ages,self.mass_range)
        self.isochrones = []

        for age in np.unique(ages):
            idx = ages == age
            # Append to isochrones list to make it easier to turn off/on
            self.isochrones.append(plt.plot(mass[idx],radius[idx],'-',color='tab:grey',zorder=0,lw=1)[0])
            self.isochrones.append(plt.text(mass[idx][-1]+0.01,radius[idx][-1],'%.0f Myr'%age,ha='left',va='center',color='k',fontsize=9))
        
        # Plot 1 Gyr isochrone over restricted range of masees
        mass,radius,ages = isochrones.BHAC15([1000],[0.05,0.8])
        self.isochrones.append(plt.plot(mass,radius,'k-',color='k',zorder=0,lw=1,alpha=0.8)[0])
        self.isochrones.append(plt.text(mass[-1]+0.01,radius[-1]-0.02,'MS',ha='left',va='center',color='k',fontsize=9))

        # Draw the adjustable isochrone
        mass,radius,ages = isochrones.BHAC15([self.age],self.mass_range)
        self.isochrone, = plt.plot(mass,radius,'-',color='tab:blue',zorder=0,lw=2)
        self.isochrone_label = plt.text(mass[-1]-0.02,radius[-1],'%i Myr'%self.age,\
            va='center',ha='right',color='w',bbox=dict(facecolor='tab:blue', alpha=0.5))

        # Plot the points
        self.stars, = plt.plot([self.mass1,self.mass2],[self.radius1,self.radius2],'o',ms=7,color='tab:red')
        plt.xlim(0,0.8)
        plt.ylim(0.0,2.2)
        plt.xlabel('Mass ($M_{\odot}$)')
        plt.ylabel('Radius ($R_{\odot}$)')
        plt.title('Mass$-$Radius Diagram')
        # Draw the BD and fully convective boundaries
        self.boundaries = []
        self.boundaries.append(plt.axvspan(0.325,0.375,color='#DDDDDD',zorder=0,alpha=0.8))
        self.boundaries.append(plt.text(0.35,1.3,'Fully convective boundary',fontsize=9,rotation=90,color='#555555',ha='center',va='center'))
        self.boundaries.append(plt.axvspan(-0.01,0.08,color='#DDDDDD',zorder=0,alpha=0.8))
        self.boundaries.append(plt.text(0.04,1.3,'Brown dwarfs',fontsize=9,rotation=90,color='#555555',ha='center',va='center'))
        self.fig.canvas.toolbar_position = 'left'
        self.fig.set_label(' ')

        # Define widgets
        self.m1 = widgets.FloatText(value=self.mass1,description='$M_{1}$ (M$_{\odot}$)')
        self.m2 = widgets.FloatText(value=self.mass2,description='$M_{2}$ (M$_{\odot}$)')
        self.r1 = widgets.FloatText(value=self.radius1,description='$R_{1}$ (R$_{\odot}$)')
        self.r2 = widgets.FloatText(value=self.radius2,description='$R_{2}$ (R$_{\odot}$)')
        self.iso_age = widgets.IntSlider(min=1,max=100,step=1,value=self.age,description='Age (Myr)',continuous_update=False)
        self.show_isochrone = widgets.Checkbox(description='Adjustable isochrone',value=True)
        self.show_isochrones = widgets.Checkbox(description='Show isochrones',value=True)
        self.show_boundaries = widgets.Checkbox(description='Show boundaries',value=True)
        self.show_grid = widgets.Checkbox(description='Show grid',value=False)

        # Monitor for updates
        self.m1.observe(self.update_points,'value')
        self.m2.observe(self.update_points,'value')
        self.r1.observe(self.update_points,'value')
        self.r2.observe(self.update_points,'value')
        self.iso_age.observe(self.update_isochrone,'value')
        self.show_isochrone.observe(self.update_show_isochrone,'value')
        self.show_isochrones.observe(self.update_show_isochrones,'value')
        self.show_boundaries.observe(self.update_show_boundaries,'value')
        self.show_grid.observe(self.update_show_grid,'value')

        controls = widgets.VBox([self.m1,self.m2,self.r1,self.r2,self.iso_age,\
            self.show_isochrone,self.show_isochrones,self.show_boundaries,self.show_grid])

        # Add to children
        self.children = [output,controls]        

    def update_points(self, change):
        self.stars.set_data([self.m1.value,self.m2.value],[self.r1.value,self.r2.value])

    def update_isochrone(self, change):
        if change.new > 100:
            self.iso_age.value = 100
        mass,radius,ages = isochrones.BHAC15([self.iso_age.value],self.mass_range)
        self.isochrone.set_data(mass,radius)
        self.isochrone_label.set_text('%i Myr'%self.iso_age.value)
        self.isochrone_label.set_position((mass[-1]-0.02,radius[-1]))

    def update_show_isochrone(self, change):
        self.isochrone.set_visible(change.new)
        self.isochrone_label.set_visible(change.new)
        self.iso_age.disabled = not change.new

    def update_show_isochrones(self, change):
        for i in self.isochrones:
            i.set_visible(change.new)

    def update_show_boundaries(self, change):
        for i in self.boundaries:
            i.set_visible(change.new)

    def update_show_grid(self, change):
        self.ax.grid(change.new)

class Test(widgets.HBox):
    
    def __init__(self, amplitude=0.5, freq=5.0):
        super().__init__()
        output = widgets.Output()

        # Draw data points and initial fits
        with output:
            self.fig, self.ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

        # Plot a test curve
        self.xx = np.arange(0, 2*np.pi,0.01)
        self.sine, = plt.plot(self.xx, amplitude*np.sin(freq*self.xx))
        plt.ylim(-1.1,1.1)
        plt.title('Things seem to be working')
        self.fig.canvas.toolbar_position = 'left'
        self.fig.set_label(' ')


        # Define widgets
        self.freq = widgets.FloatSlider(value=freq,min=0,max=10,description='Frequency')
        self.amplitude = widgets.FloatSlider(value=amplitude,min=0,max=1,step=0.01,description='Amplitude')
        self.show_curve = widgets.Checkbox(value=True,description='Show curve')

        # Monitor for updates
        self.freq.observe(self.update_curve,'value')
        self.amplitude.observe(self.update_curve,'value')
        self.show_curve.observe(self.update_show_curve,'value')

        controls = widgets.VBox([self.freq,self.amplitude,self.show_curve])

        # Add to children
        self.children = [output,controls]        

    def update_curve(self, change):
        new_curve = self.amplitude.value * np.sin(self.xx * self.freq.value)
        self.sine.set_ydata(new_curve)

    def update_show_curve(self,change):
        self.sine.set_visible(change.new)



