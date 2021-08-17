import os, time, datetime
import multiprocessing

SESSION_LENGTH = 240 # in minutes
USED_TOKEN = multiprocessing.Value('i', 0) # a global variable. makes sure we call the worker only once

def clock():
    current_time = datetime.datetime.utcnow()
    return current_time.strftime('%X') 

def dummy_worker():
    if USED_TOKEN.value:
        print("You can only open the session once.")
    else:
        # makes sure the dummy woker is called only once
        USED_TOKEN.value = 1 
        # warn the user that they cannot run twice
        print("Session of " + str(SESSION_LENGTH) + " minutes started at " + clock() + " UTC.")
        # parameters of the pauses
        pause = 5
        chrono = 0
        nb_sleeps = int(SESSION_LENGTH/pause)
        for i in range(nb_sleeps):
            time.sleep(pause*60)
            chrono = chrono + pause
            print('Session open for '+str(chrono)+' minutes...', end="\r") # perdiodic refreshing
        # final message    
        print("Session finished at " + clock() + " UTC. The Binder server will shut down after 10 minutes of inactivity.")
    return

def start_session():
    multiprocessing.Process(target=dummy_worker).start()
