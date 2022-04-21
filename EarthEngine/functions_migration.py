#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
functions for riverbank migration calculations

Ted Langhorst
April 2020

"""
import ee
import functions_watermask as fw

#rate calculations
#find areas of change, calculate bank to bank distance across areas
def diffMap(wm1, wm2, roi, SCALE, CRS):
  dt = wm1.getNumber('year').subtract(wm2.get('year')).abs()
  
  #pixels that were water in wm1 and land in wm2.
  diff = wm1.select('mask').subtract(wm2.select('watermask')).eq(1).selfMask().rename('diff')
  
  #calculate distance across the changed area to the old banks.
  bankDist = diff.focal_max().cumulativeCost(wm1.select('banks'),1000,False).mask(diff).rename('bankDist')
  
  #islands are disconnected from banks and distance calc fails. Need to go back and fix these.
  islands = bankDist.eq(ee.Image.constant(0)).selfMask().clip(roi).rename('islands')
  
  islandCentroids = islands.reduceToVectors(
    geometry = roi,
    scale = SCALE,
    crs = CRS,
    geometryType = "centroid",
    maxPixels = int(1E10),
    tileScale = 4)
  
  islandCenterPix = islands.Not().paint(islandCentroids,1).rename('islandCenter')
  
  #calculate distance to center of island.
  islandDist = islands.focal_max().cumulativeCost(islandCenterPix,1000,False).mask(islands).rename("islandDist")
  
  dist = bankDist.unmask().add(islandDist.unmask()).mask(diff).rename('dist')
  
  bankRate = dist.unmask().mask(wm2.select('banks').gt(0)).divide(dt).rename('bankRate')
  
  diff = diff.multiply(ee.Image.pixelArea()) #change from binary to area
  
  noisyDiff = (wm1.select('noisyRiverMask').subtract(wm2.select('noisyWaterMask'))
               .eq(1).selfMask().rename('noisyDiff')
               .multiply(ee.Image.pixelArea()))

  
  return diff.addBands([dist, bankRate, noisyDiff]).clip(roi)



def calcMigration(yr1, dt, avgWindow, sword_clip, SCALE, CRS):
  roi = ee.Geometry(sword_clip.get('roi'))
  # width = ee.Number(sword_clip.get('width'))
  
  #water masks
  yearList = ee.List([yr1, yr1+avgWindow, yr1+dt, yr1+dt+avgWindow])
  
  
  # pekel
  # wm1 = fw.pekelMask(yearList.get(0), yearList.get(1), sword_clip)
  # wm2 = fw.pekelMask(yearList.get(2), yearList.get(3), sword_clip)
  # noDataMask = wm1.select("noData").Or(wm2.select("noData")).clip(roi)
  # wm1 = wm1.mask(noDataMask.Not())
  # wm2 = wm2.mask(noDataMask.Not())
  
  #Pickens
  wm1 = fw.pickensMask(yearList.get(0), yearList.get(1), sword_clip)
  wm2 = fw.pickensMask(yearList.get(2), yearList.get(3), sword_clip)
  noDataMask = wm1.select("noData").Or(wm2.select("noData")).clip(roi)
  wm1 = wm1.mask(noDataMask.Not())
  wm2 = wm2.mask(noDataMask.Not())
  
  
  #errosion and accretion
  A = diffMap(wm1, wm2, roi, SCALE, CRS).rename(['ADiff','ADist','ARate','noisyADiff'])
  E = diffMap(wm2 ,wm1, roi, SCALE, CRS).rename(['EDiff','EDist','ERate','noisyEDiff'])
  
  #probably a better way to do this
  wm1 = wm1.rename(["t1_mask","t1_banks","t1_bankAspect","t1_bankCurv","t1_bankLen","t1_watermask","t1_noData","t1_noisyWaterMask","t1_noisyRiverMask"])
  wm2 = wm2.rename(["t2_mask","t2_banks","t2_bankAspect","t2_bankCurv","t2_bankLen","t2_watermask","t2_noData","t2_noisyWaterMask","t2_noisyRiverMask"])
  
  
  wm1 = wm1.addBands(
    srcImg = wm1.select('t1_mask').multiply(ee.Image.pixelArea()),
    overwrite = True)
  wm2 = wm2.addBands(
    srcImg = wm2.select('t2_mask').multiply(ee.Image.pixelArea()),
    overwrite = True)
  
  
  
  imgOut = wm1.addBands(wm2).addBands(A).addBands(E).set(
    'year2', wm2.get('year'),
    'SCALE', SCALE,
    'CRS', CRS)
  return imgOut
