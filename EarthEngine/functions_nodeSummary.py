#!/usr/bin/env python3
# -*- coding utf-8 -*-
"""
functions for summarizing data to SWORD

Ted Langhorst
April 2020

"""
import ee


def groupReduce(dataImg, pixMap, redux):
  #when using a weighted reducer (e.g. mean), need to pass this function the ".unweighted()" version. 
  #Not exactly sure why, maybe just how multiple bands works with weighted reducers.
  
  internalReducer = redux.forEach(dataImg.bandNames()).group(
        groupField = dataImg.bandNames().length(),
        groupName = "node_id")

  imgStats = dataImg.addBands(pixMap).reduceRegion(
    reducer = internalReducer,
    scale = dataImg.getNumber("SCALE"),
    crs = dataImg.getString("CRS"),
    bestEffort = False,
    maxPixels = int(1E10),
    tileScale = 4).get("groups")

  featOut = ee.List(imgStats).map(lambda m: ee.Feature(None,m))
  
  return ee.FeatureCollection(featOut)

#Special group reduce that uses x,y and weights
def xyWeightedGroupReduce(dataImg, weights, pixMap, outName):
  xStats = dataImg.select('x').addBands(weights).addBands(pixMap) \
  .reduceRegion(
    reducer = ee.Reducer.mean().splitWeights().group(
      groupField = 2,
      groupName = "node_id"),
    scale = dataImg.getNumber("SCALE"),
    crs = dataImg.getString("CRS"),
    maxPixels = 1E10,
    tileScale = 4).get("groups")
    
  yStats = dataImg.select('y').addBands(weights).addBands(pixMap) \
  .reduceRegion(
    reducer = ee.Reducer.mean().splitWeights().group(
      groupField = 2,
      groupName = "node_id"),
    scale = dataImg.getNumber("SCALE"),
    crs = dataImg.getString("CRS"),
    maxPixels = 1E10,
    tileScale = 4).get("groups")

  #changes from the weird list/object format to feature collection with property "fieldName" instead of "mean"
  def stats2fc(stats, fieldName):
    def internalMap(s):
      tempFeat = ee.Feature(None,s)
      return ee.Feature(None).set(fieldName,tempFeat.get("mean"),"node_id",tempFeat.get("node_id"))
      
    return ee.FeatureCollection(ee.List(stats).map(internalMap))

  
  featX = stats2fc(xStats, 'mx')
  featY = stats2fc(yStats, 'my')
  featXY = nodeJoin(featX,featY)
  ang = xy2ang_fc(featXY, outName)
  
  return ang


def ang2xy(image):
  xy = image.cos().rename("x") \
      .addBands(image.sin().rename("y"))
  return xy

def xy2ang_fc(fc, outName):
  theta = ee.String(outName)
  R = theta.cat("R")
      
  def internalMap(f):
    tempDict = ee.Algorithms.If(
      condition = f.propertyNames().contains("mx"),
      falseCase = ee.Dictionary.fromLists([theta, R],[-9999, -1]),
      trueCase = ee.Dictionary.fromLists([theta, R],
        [ee.Number(f.get("mx")).atan2(f.get("my")),
        ee.Number(f.get("mx")).pow(2) \
        .add(ee.Number(f.get("my")).pow(2)).sqrt()]))
    return ee.Feature(None).set("node_id",f.get("node_id")).set(tempDict)

  return fc.map(internalMap)

def setConstant(fc, propDict):
  def internalMap(feat):
    return ee.Feature(feat).set(propDict)
  return ee.FeatureCollection(fc).map(internalMap)


def nodeJoin(f1,f2):
  def internalMap(fMap):
    match = f2.filter(ee.Filter.eq("node_id",fMap.get("node_id")))
    fOut = ee.Algorithms.If(
      condition = match.size().gt(0),
      trueCase = fMap.copyProperties(match.first()),
      falseCase = fMap)
    return fOut

  return f1.map(internalMap)

def addMissingProperties(feat):
    nProperties = ee.Number(feat.propertyNames().length())
    featOut = ee.Algorithms.If(
                  condition = nProperties.lte(23),
                  falseCase = feat,
                  trueCase = feat.set({
                          'nProperties':nProperties,
                          'ADiff': -1,
                          'EDiff': -1,
                          't1_mask': -1,
                          't2_mask': -1,
                          't1_bankLen': -1,
                          't2_bankLen': -1,
                          'noisyADiff': -1,
                          'noisyEDiff': -1,
                          'ARate': -1,
                          'ERate': -1,
                          'ADir': -9999,
                          'ADirR': -1,
                          'EDir': -9999,
                          'EDirR': -1}))
    return featOut

#
def calcStats(mData, pixMap, swordClip):
  pixMapMask = (pixMap.clip(mData.geometry()) 
  .eq(ee.Image.constant(swordClip.aggregate_array('node_id')))
  .reduce(ee.Reducer.anyNonZero()))
  pixMap = pixMap.updateMask(pixMapMask)
  
  #reduce the data in different ways
  sumData = mData.select(["ADiff","EDiff","t1_mask","t2_mask","t1_bankLen","t2_bankLen","noisyADiff","noisyEDiff"])
  sumFeat = groupReduce(sumData, pixMap, ee.Reducer.sum())
  nodeData = nodeJoin(swordClip,sumFeat)
  
  meanData = mData.select(["ARate","ERate"])
  meanFeat = groupReduce(meanData, pixMap, ee.Reducer.mean().unweighted())
  nodeData = nodeJoin(nodeData,meanFeat)
    
  EDir = xyWeightedGroupReduce(ang2xy(mData.select("t1_bankAspect")), \
  mData.select("EDist").mask(mData.select("t1_banks")),pixMap, "EDir")
  nodeData = nodeJoin(nodeData,EDir)
   
  ADir = xyWeightedGroupReduce(ang2xy(mData.select("t2_bankAspect")),
  mData.select("ADist").mask(mData.select("t2_banks")),pixMap, "ADir")
  nodeData = nodeJoin(nodeData,ADir)
  
  #add null data to any feature that is missing SCREAMinG data.
  nodeData = nodeData.map(addMissingProperties)
  
  dataOut = ee.FeatureCollection(nodeData).copyProperties(mData)
  
  return dataOut






















