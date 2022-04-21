#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 10:51:26 2020

@author: Ted
"""

import os 
import time
import datetime
import ee



def maximum_no_of_tasks(MaxNActive, waitingPeriod):
  ##maintain a maximum number of active tasks
  time.sleep(10)
  ## initialize submitting jobs
  ts = list(ee.batch.Task.list())

  NActive = 0
  for task in ts:
       if ('RUNNING' in str(task) or 'READY' in str(task)):
           NActive += 1
  ## wait if the number of current active tasks reach the maximum number
  ## defined in MaxNActive
  while (NActive >= MaxNActive):
      time.sleep(waitingPeriod) # if reach or over maximum no. of active tasks, wait
      ts = list(ee.batch.Task.list())
      NActive = 0
      for task in ts:
        if ('RUNNING' in str(task) or 'READY' in str(task)):
          NActive += 1
  return

def taskMonitor(exportTask, waitingPeriod):
  
  if (waitingPeriod < 10):
    waitingPeriod = 10  
  state = exportTask.status().get('state')
  
  #still running
  while (state=='RUNNING' or state=='READY'):
    time.sleep(waitingPeriod)
    state = exportTask.status().get('state')

  #complete
  # if (state=='COMPLETED'):
    # os.system('afplay Sounds/bikeHorn.wav')
  # else:
    # os.system('afplay Sounds/failTrombone.wav')
  
  status = exportTask.status()  
  dSeconds = (status.get('update_timestamp_ms')-status.get('start_timestamp_ms'))//1000
  dt = datetime.timedelta(0,dSeconds)
  print(state + "\truntime: " + str(dt))
  return


def cancelAllTasks():
  for t in ee.batch.Task.list():
    if ('RUNNING' in str(t) or 'READY' in str(t)):
      #for some reach the cancel method throws a 404 error, but only after cancelling the task. Skipping the exception allows you to keep going and cancel them all.
        try:
            t.cancel()
        except Exception:
            pass
            