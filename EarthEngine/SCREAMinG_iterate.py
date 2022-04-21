#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 12:13:13 2020

Ted Langhorst
tlang@live.unc.edu

"""

import ee
import time
import numpy as np
ee.Initialize()  

pixMap = ee.Image("users/tedlanghorst/SWORD_pixmap_v09")

allReaches = (ee.FeatureCollection('users/tedlanghorst/SWORD_reach_v09')
             .filterMetadata('width','greater_than',150)
             .aggregate_array('reach_id').getInfo())

riverReachFilter = filter(lambda r: r%10==1, allReaches)
reachList = list(riverReachFilter)


#reach limits and step
RStart = 10000000000
REnd = 100000000000
dR = 1000000000
reachBounds = range(RStart,REnd+1,dR)

def groupFilt(reach):
  return reach>=RStart and reach<=REnd

allReaches = list(filter(groupFilt, reachList))
groupList = np.unique([r // dR for r in allReaches])
reachBounds = [r * dR for r in groupList]
len(reachBounds)

#settings
YR1 = 2001
DT = 18
AVG_WINDOW = 0
WIDTH_BUFF = 2
SCALE = 30
CRS = "epsg:4326"


import functions_tasking as ft
from functions_watermask import filterSword
from functions_migration import calcMigration
from functions_nodeSummary import calcStats


def SCREAMinG(reachList_internal):
  def reachCalcs(reach):
    #make watermask from JRC dataset
    swordClip = filterSword(reach,reach, WIDTH_BUFF)
    #migration raster data
    mData = calcMigration(YR1, DT, AVG_WINDOW, swordClip, SCALE, CRS)
    #calculate node-level data
    reachStats = calcStats(mData, pixMap, swordClip)
    return reachStats
  
  swordStats = ee.FeatureCollection(ee.List(reachList_internal).map(reachCalcs)).flatten()
  return swordStats


"""
single map
"""
def inRange(reach):
  return reach>=R1 and reach<=R2

for idx in range(len(reachBounds)):
  R1 = reachBounds[idx]
  R2 = R1 + dR
  
  reachListFilt = list(filter(inRange, reachList))
  
  # reachListFilt = reachListFilt[0:150]
  
  if len(reachListFilt)==0:
    # print('no reaches in range {0}:{1}'.format(R1,R2))
    continue

  # #the business
  swordStats = SCREAMinG(reachListFilt)
  
  #export
  exportTask = ee.batch.Export.table.toDrive(
    collection = swordStats,
    description = "{0}_{1}_{2}_{3:02}_{4:02}".format(R1,R2,YR1,DT,AVG_WINDOW),
    folder = "SCREAMING_uncertainty_Pickens_noAvg_08042021_v09",
    fileFormat = 'CSV')
  exportTask.start()
  now = time.localtime()
  print("{0:02}:{1:02}\t submitting {2} reaches from range {3}:{4}".format(
    now.tm_hour,now.tm_min,len(reachListFilt),R1,R2))
  ft.taskMonitor(exportTask, 10)
