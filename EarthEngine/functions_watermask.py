#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
functions for creating watermasks for riverbank migration

Ted Langhorst
April 2020

"""
import ee
# from random import randint

jrcMetadata = ee.Image("JRC/GSW1_2/Metadata")
jrcYearly = ee.ImageCollection("JRC/GSW1_2/YearlyHistory")
jrcMonthly = ee.ImageCollection("JRC/GSW1_2/MonthlyHistory")
pickensYearly = ee.ImageCollection("projects/glad/water/annual")
sword = ee.FeatureCollection("users/tedlanghorst/SWORD_node_v09")

def filterSword(r1, r2, WIDTH_BUFF):
  sword_clip = sword.filterMetadata('reach_id','equals',ee.Number(r1)) #v08
  
  swordWidth = ee.Number(sword_clip.reduceColumns(ee.Reducer.mean(),['width']).getNumber("mean"))
  #replace missing width values with static.
  width = ee.Number(ee.Algorithms.If(
    condition = swordWidth.eq(1),
    trueCase = 500,
    falseCase = swordWidth))
  
  roi = sword_clip.geometry().buffer(width.multiply(WIDTH_BUFF))
  
  return sword_clip.set('width',width,'widthBuff',WIDTH_BUFF,'roi',roi)


daKernel = ee.Kernel.fixed(3,3,
[[0.5,  0.5, 0.5],
[0.5,  -1,  0.5],
[0.5,  0.5, 0.5]])

#kernel that gives unique value for each bank orientation in 3x3
directionKernel = ee.Kernel.fixed(3,3,
[[128,  1, 2],
 [ 64,  256, 4],
 [ 32, 16, 8]])

#translate kernel values to angles. Found these using a MATLAB script.
kernelValues = ee.List([1, 4, 16, 64, 129, 3, 6, 12, 24, 48, 96, 192, 131, 7, 14, 28, 56, 112, 224, 193, 135, 15, 30, 60, 120, 240, 225, 195, 143, 31, 62, 124, 248, 241, 227, 199, 159, 63, 126, 252, 249, 243, 231, 207])

aspectKey = ee.List([4.712389, 0.000000, 1.570796, 3.141593, 4.467410, 4.957368, 6.038207, 0.244979, 1.325818, 1.815775, 2.896614, 3.386571, 4.712389, 5.497787, 0.000000, 0.785398, 1.570796, 2.356194, 3.141593, 3.926991, 5.092895, 5.902679, 0.380506, 1.190290, 1.951303, 2.761086, 3.522099, 4.331883, 5.497787, 0.000000, 0.785398, 1.570796, 2.356194, 3.141593, 3.926991, 4.712389, 5.695183, 0.588003, 0.982794, 2.158799, 2.553590, 3.729595, 4.124386, 5.300392])


#calculates the derivative of the bank angles
#transform angles to x,y vecs, find derivatives, calculate angle.
def circularDerivative(image):
  xy = image.cos().rename('x').addBands(image.sin().rename('y'))
  dxy = xy.convolve(daKernel)
  dTheta = dxy.select('x').atan2(dxy.select('y'))
  return dTheta


#fill in all non water regions with area smaller than areaMax
def fillBars(image, roi):
  def calcArea(feat):
    a = feat.geometry().area(1000)
    bars = a.lte(areaMax)
    return feat.set("fill", bars)
  
  areaMax = ee.Geometry(roi).area().divide(25)
  barPolys = image.unmask().Not().selfMask().reduceToVectors(
    geometry = roi,
    scale = 30,
    eightConnected = False,
    maxPixels = 1E10,
    bestEffort = True).map(calcArea)
  
  filled = image.selfMask().paint(barPolys,"fill");
  
  return filled

def largestBody(image, roi):
  #water body outlines
  maskVec = image.selfMask().reduceToVectors(
    geometry = roi,
    scale = 30,
    eightConnected = False,
    maxPixels = 1E10)
  
  #pick the largest water body in region, paint on new image
  riverVec = ee.Feature(maskVec.sort('count',False).first())
  riverMask = ee.Image().toByte().paint(riverVec,1)
  
  return riverMask

def bankCalcs(riverMask):
  bankAspect = (riverMask.unmask().convolve(directionKernel)
  .remap(kernelValues, aspectKey).rename("bankAspect"))
  
  banks = bankAspect.add(9999).neq(0).selfMask().rename('banks')
  
  bankCurvature = bankAspect.convolve(daKernel).rename('bankCurv')
  bankLen = banks.convolve(ee.Kernel.euclidean(1)).rename('bankLen')
  
  return riverMask.addBands([banks, bankAspect, bankCurvature, bankLen])

def addNoise(year1, year2, roi, errorRates):
  nobs = jrcMetadata.select('valid_obs').reduceRegion(ee.Reducer.mean(),roi).getNumber('valid_obs')
  nObsYr = nobs.divide(36).toInt()
  
  def addYearNoise(img):
    def addObsNoise(num):
      #create random noise image
      monthNoise = ee.Image.random(img.getNumber('system:time_start').add(num))
       
      noisySeasonal = (img.eq(2).unmask()
      .subtract(monthNoise.lte(errorRates['seasonO'])).eq(1) #remove some water
      .add(monthNoise.gte(ee.Number(1).subtract(errorRates['seasonC']))).gte(1)) #add some water
      
      noisyPerm = (img.eq(3).unmask()
      .subtract(monthNoise.lte(errorRates['permO'])).eq(1) #remove some water
      .add(monthNoise.gte(ee.Number(1).subtract(errorRates['permC']))).gte(1)) #add some water
      
      return noisyPerm.Or(noisySeasonal)
    
    #duplicate the yearly map with new noise for each observation
    monthMasks = ee.List.sequence(0,nObsYr).map(addObsNoise)
    #aggregate to a yearly mask
    yearMask = ee.ImageCollection(monthMasks).reduce(ee.Reducer.mean()).gte(0.5)
    
    return yearMask
    
  yearMasks = jrcYearly.filter(ee.Filter.calendarRange(year1,year2,'year')).map(addYearNoise)
  
  noisyMask = ee.ImageCollection(yearMasks).reduce(ee.Reducer.mean()).gte(0.5)
  return noisyMask

def addNoise_pekel_distance(year1, year2, roi, OErr,CErr):
  bins = ee.List.sequence(1,10,1) 

  def addYearNoise(watermask):
    def addObsNoise(num):
      #create random noise image
      noise = ee.Image.random(seed.add(num))
      
      WDist = watermask.cumulativeCost(watermask.Not(),90).mask(watermask)
      WDist = WDist.add(WDist.eq(0).multiply(90))
      WBinImage = WDist.divide(10).add(1).toInt8()
      OThreshImg = WBinImage.remap(bins,OErr).unmask()
      
      LDist = watermask.Not().cumulativeCost(watermask,90).mask(watermask.Not())
      LDist = LDist.add(LDist.eq(0).multiply(90))
      LBinImage = LDist.divide(10).add(1).toInt8()
      CThreshImg = LBinImage.remap(bins,CErr).unmask()
      
      noisyImg = (watermask
      .subtract(noise.lte(OThreshImg)).eq(1) #remove some water
      .add(noise.gte(ee.Image.constant(1).subtract(CThreshImg))).gte(1)) #add some water
  
      return noisyImg
    
    seed = watermask.getNumber("system:time_start")
    watermask = watermask.gte(2);
    #duplicate the yearly map with new noise for each observation
    monthMasks = ee.List.sequence(0,12).map(addObsNoise)
    #aggregate to a yearly mask
    yearMask = ee.ImageCollection(monthMasks).reduce(ee.Reducer.mean()).gte(0.5)
    return yearMask
  
  yearMasks = jrcYearly.filter(ee.Filter.calendarRange(year1,year2,'year')).map(addYearNoise)
  noisyMask = ee.ImageCollection(yearMasks).reduce(ee.Reducer.mean()).gte(0.5)
  return noisyMask

def addNoise_pickens(year1, year2, roi, OErr,CErr):
  bins = ee.List.sequence(1,10,1) 

  def addYearNoise(watermask):
    def addObsNoise(num):
      #create random noise image
      noise = ee.Image.random(seed.add(num))
      
      WDist = watermask.cumulativeCost(watermask.Not(),90).mask(watermask)
      WDist = WDist.add(WDist.eq(0).multiply(90))
      WBinImage = WDist.divide(10).add(1).toInt8()
      OThreshImg = WBinImage.remap(bins,OErr).unmask()
      
      LDist = watermask.Not().cumulativeCost(watermask,90).mask(watermask.Not())
      LDist = LDist.add(LDist.eq(0).multiply(90))
      LBinImage = LDist.divide(10).add(1).toInt8()
      CThreshImg = LBinImage.remap(bins,CErr).unmask()
      
      noisyImg = (watermask
      .subtract(noise.lte(OThreshImg)).eq(1) #remove some water
      .add(noise.gte(ee.Image.constant(1).subtract(CThreshImg))).gte(1)) #add some water
  
      return noisyImg
    
    seed = watermask.getNumber("system:time_start")
#    seed = ee.Number(42)
    watermask = watermask.gte(50);
    #duplicate the yearly map with new noise for each observation
    monthMasks = ee.List.sequence(0,12).map(addObsNoise)
    #aggregate to a yearly mask
    yearMask = ee.ImageCollection(monthMasks).reduce(ee.Reducer.mean()).gte(0.5)
    return yearMask
  
  yearMasks = pickensYearly.filter(ee.Filter.calendarRange(year1,year2,'year')).map(addYearNoise)
  noisyMask = ee.ImageCollection(yearMasks).reduce(ee.Reducer.mean()).gte(0.5)
  return noisyMask
                                   



#watermask function
#samples JRC masks and does some geomorph calcs for a single image.
def pekelMask(year1, year2, sword_clip):
  roi = sword_clip.get('roi')
  
  watermask = (jrcYearly.filter(ee.Filter.calendarRange(year1,year2,'year'))
               .map(lambda img: img.selfMask())
               .reduce(ee.Reducer.median())
               .clip(roi)
               .gte(2)
               .rename('watermask'))
  
   
  # errorRates = {'permO':0.02,
  #               'permC':0.005,
  #               'seasonO':0.25,
  #               'seasonC':0.02}
  # noisyWaterMask = addNoise(year1, year2, roi, errorRates).rename('noisyWaterMask')
  
  OErr = ee.List([0.55, 0.3, 0.15, .1, .08, .05, .03, .02, .01, 0])
  CErr = ee.List([0.15, 0.05, 0.3, .02, .01, .005, 0, 0, 0, 0])
  noisyWaterMask = addNoise_pekel_distance(year1, year2, roi, OErr,CErr).rename("noisyWaterMask")
  
  # filledMask = fillBars(watermask, roi)
  # filledMask = watermask
  
  # sources = ee.Image().toByte().paint(sword_clip, 1).And(filledMask)
  # ccMask = (filledMask.Not().cumulativeCost(
  #   source = sources, 
  #   maxDistance = width.multiply(sword_clip.get("widthBuff")))
  #   .eq(0).And(watermask).unmask() \
  #   .focal_max().focal_min()) #grwl connected pts and closure

  # riverMask = largestBody(ccMask,roi).rename("mask")
  # riverMask = ccMask.rename("mask")
  
  
  riverMask = (largestBody(watermask,roi)
                .unmask()
                .focal_max().focal_min()
                .rename("mask"))
  
  noisyRiverMask = (largestBody(noisyWaterMask,roi)
                .unmask()
                .focal_max().focal_min()
                .rename("noisyRiverMask"))

  
  noData = (jrcYearly.filter(ee.Filter.calendarRange(year1,year2,'year'))
            .map(lambda img: img.eq(0))  #pixels with no data = 1 (collection)
            .reduce(ee.Reducer.mean()).gte(0.75)      #pixels with no data in most of the images
            .unmask()
            .And(riverMask.Not()).rename("noData"))   #dont mask pixels that are water after closure
  
            
  imgOut = (bankCalcs(riverMask)
            .addBands(watermask)
            .addBands(noData)
            .addBands(noisyWaterMask)
            .addBands(noisyRiverMask).clip(roi)
            .set('year',year1,'nAvg', ee.Number(year2).subtract(year1).add(1)))
  
  return imgOut


def pickensMask(year1, year2, sword_clip):
  roi = sword_clip.get('roi')
  
  watermask = (pickensYearly.filter(ee.Filter.calendarRange(year1,year2,'year'))
               .map(lambda img: img.unmask())
               .reduce(ee.Reducer.mean()).clip(roi).gte(50)
               .rename('watermask'))
  
  OErr = ee.List([0.4, 0.2, 0.1, .08, .05, .03, .02, .02, .01, 0])
  CErr = ee.List([0.2, 0.1, 0.5, .03, .02, .01, .005, 0, 0, 0])
  noisyWaterMask = addNoise_pickens(year1, year2, roi, OErr,CErr).rename("noisyWaterMask")
  
  # filledMask = fillBars(watermask, roi)
  # filledMask = watermask
  
  # sources = ee.Image().toByte().paint(sword_clip, 1).And(filledMask)
  # ccMask = (filledMask.Not().cumulativeCost(
  #   source = sources, 
  #   maxDistance = width.multiply(sword_clip.get("widthBuff")))
  #   .eq(0).And(watermask).unmask() \
  #   .focal_max().focal_min()) #grwl connected pts and closure

  # riverMask = largestBody(ccMask,roi).rename("mask")
  # riverMask = ccMask.rename("mask")
  riverMask = (largestBody(watermask,roi)
             .unmask()
             .focal_max().focal_min()
             .rename("mask"))
  
  noisyRiverMask = (largestBody(noisyWaterMask,roi)
             .unmask()
             .focal_max().focal_min()
             .rename("noisyRiverMask"))
  
  noData = (pickensYearly.filter(ee.Filter.calendarRange(year1,year2,'year'))
            .map(lambda img: img.gte(0).unmask().Not())
            .reduce(ee.Reducer.mean()).gte(0.75)
            .And(riverMask.Not()).rename("noData"))   #dont mask pixels that are water after closure
  

            
  imgOut = (bankCalcs(riverMask).clip(roi)
            .addBands(watermask)
            .addBands(noData) #nodata kluge to get Pickens data to play nice.
            .addBands(noisyWaterMask)
            .addBands(noisyRiverMask).clip(roi)
            .set('year',year1,'nAvg', ee.Number(year2).subtract(year1).add(1)))
  
  return imgOut

def landsatMask(year1, year2, sword_clip):
  import functions_landsat as lsfun
  
  roi = sword_clip.get('roi')
  
  collection = (lsfun.merge_collections()
                .filterBounds(roi)
                .filterMetadata('CLOUD_COVER_LAND','less_than',50))
  
  #classify water and veg
  def internalMap(img):
    img = lsfun.CalculateWaterAddFlagsSR(img)
    return img.updateMask(img.select('fmask').lte(1))
  
  collection = (collection.filter(ee.Filter.calendarRange(year1,year2,'year'))
                .map(internalMap))
  
  watermask = (collection.select('corridor')
               .reduce(ee.Reducer.mean())
               .gte(0.5)
               .focal_max().focal_min()
               .rename('watermask')
               .set('nImages',collection.size()))
  
  riverMask = largestBody(watermask,roi).rename('mask')

  imgOut = (bankCalcs(riverMask).clip(roi)
            .addBands(watermask)
            .addBands(ee.Image.constant(0)) #temporary kluge to get Pickens data to play nice with the Pekel method.
            .set('year',year1,'nAvg', ee.Number(year2).subtract(year1).add(1)))
  
  return imgOut






