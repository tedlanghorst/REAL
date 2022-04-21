#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 12:13:13 2020

@author: Ted


TODO
- check all angle reference/orientation
- implement circular derivative for curvature
- cutoff detection (still water, lost connection)

"""

import ee
import csv
import time
import numpy as np
import pandas as pd
ee.Initialize()  


pixMap = ee.Image("users/tedlanghorst/SWORD_pixmap_v09")


allReaches = (ee.FeatureCollection('users/tedlanghorst/SWORD_reach_v09')
             .filterMetadata('width','greater_than',150)
             .aggregate_array('reach_id').getInfo())

riverReachFilter = filter(lambda r: r%10==1, allReaches)
reachList = list(riverReachFilter)


# pixMap_reachList = [] 
# with open('/Users/Ted/GDrive/SWORD/SWORD_reachList_river_meanGt150.csv') as csvDataFile:
#     csvReader = csv.reader(csvDataFile)
#     next(csvReader) #skip header
#     for row in csvReader:
#         pixMap_reachList.append(int(row[0]))

# v1v8 = pd.read_csv('/Users/Ted/Documents/SWORD_v09/Version_Differences/diff/reachDifferences.csv')

# reachList_df = v1v8[v1v8['v01'].isin(pixMap_reachList)]
# reachList = list(reachList_df['v08'])

# Mississippi
# R1 = 74000000000
# R2 = 75000000000

## Yukon
# R1 = 81200000000
# R2 = 81300000000




## NA
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
  # ft.taskMonitor(exportTask, 30)
  # now = time.localtime()
  # print("{0:02}:{1:02} finished range {2}:{3}".format(now.tm_hour,now.tm_min,R1,R2))


"""
failed reaches. 2/12/2021
v09 run 1
35000000000_36000000000_2000_08_03
Error: Number of pixels requested from Window.max exceeds the maximum allowed (2^31)



Yukon
81210300021
81210300041
81210300051

5/22/20
81210900401
81210900211
81210700081
81210500121

81210900201
"""

"""
failed tasks using landsatMask method
74225000900_74225001000_2000_15_04
74241301600_74241301700_2000_15_04
74249100200_74249100300_2000_15_04
74252100400_74252100500_2000_15_04

"""

"""

exportTask = ee.batch.Export.image.toDrive(
  image = mData.toFloat(),
  description = "mData",
  scale = SCALE,
  crs = CRS,
  maxPixels = int(1E13)).start()









"""