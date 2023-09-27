# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 18:26:43 2020

@author: Melanie Martyn Rosco

This script is based on the Classification script from Chris Holden and Florian Beyer
"""



import pandas as pd
import numpy as np
import geopandas as gpd
import math
import sys
from osgeo import gdal, ogr, gdal_array, gdalconst, osr
from osgeo_utils import gdal_sieve
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier



# --------------------------------------------------------
# function to parse feature from GeoDataFrame in such a manner that rasterio wants them
def getFeatures(gdf):
    import json
    return [json.loads(gdf.to_json())['features'][0]['geometry']]


# Conver coordinates from one EPSG to another EPSG
def coord_conversion(inputEPSG, outputEPSG, coords):
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(inputEPSG)
    
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(outputEPSG)
    
    coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    
    x_list = []
    y_list = []
    
    for index, row in coords.iterrows():
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(row['x'],row['y'])
        point.Transform(coordTransform)
        x_list.append(point.GetX())
        y_list.append(point.GetY())
        newcoords = pd.DataFrame({'x':x_list, 'y':y_list})
    
    return newcoords


# Calculate Earth-Sun distance
def Earth_Sun_distance(year, month, day, hh, mm, ss):
    # Universal time
    UT = hh + (mm/60) + (ss/3600)
    
    # Julian Day (JD)
    A = int(year/100)
    B = 2 - A + int(A/4)
    JD = int(365.25*(year + 4716) + int(30.6001 * (month + 1)) + day + (UT/24) + B - 1524.5)
    
    # Earth-Sun distance (dES)
    D = JD - 2451545
    g  = 357.529 + 0.98560028 * D # (in degrees)
    dES = 1.00014 - 0.01671 * math.cos(g) - 0.00014 * np.cos(np.deg2rad(2*g))
    
    return dES



# Generate top of atmosphere reflectance and surface reflectance raster from DN numbers - needs dES from previous function
def refl_wv2(in_raster, absFactors, effBandWidths, ESUN, dES, theta, out_raster_TOA):
    drv = gdal.GetDriverByName('GTiff')
    raster = gdal.Open(in_raster)
    
    referenceProj = raster.GetProjection()
    referenceTrans = raster.GetGeoTransform()
    x = raster.RasterXSize
    y = raster.RasterYSize
    n = raster.RasterCount
    
    print(n)
    
    raster_array = raster.ReadAsArray()
    data = np.single(raster_array)
    data = data.T
    data[data == 0] = np.nan
    
    # Plot 4 bans of multi raster
    fig = plt.subplots(figsize=(13,7))

    plt.subplot(221)
    plt.imshow(data[:,:,1]/255)
    plt.title('Blue')
    plt.colorbar

    plt.subplot(222)
    plt.imshow(data[:,:,3]/255)
    plt.title('Green')
    plt.colorbar

    plt.subplot(223)
    plt.imshow(data[:,:,4]/255)
    plt.title('Red')
    plt.colorbar

    plt.subplot(224)
    plt.imshow(data[:,:,6]/255)
    plt.title('NIR')
    plt.colorbar

    plt.show()

    # Emtpy matrices that will store the TOA reflectance data
    refl_TOA = np.zeros(data.shape)   

    for i in range(n):
        im = data[:,:,i]
        # Calculate DN to radiance
        L = (im * absFactors[i])/effBandWidths[i]
        L = np.squeeze(L)
        # Calculates the theoretical radiance of a dark object as 1% of the max possible radiance
        L1percent = (0.01 * ESUN[i] * np.cos(np.deg2rad(theta))) / (dES**2 * math.pi)
        # Find darkest pixel in image
        Lmini = np.nanmin(L)
        # The difference between the theoretical 1% radiance of a dark object and the radiance of the darkest image pixel is due to the atm (empirical)
        Lhaze = Lmini - L1percent
        # TOA reflectance
        refl_TOA[:, :, i] = (L * math.pi * dES**2) / (ESUN[i] * np.cos(np.deg2rad(theta)))

    
    # Save to rasters
    refl_TOA = np.float32(refl_TOA)
    output_TOA = drv.Create(out_raster_TOA, x, y, n, gdal.GDT_Float32)
    if output_TOA is None:
        print('The output raster could not be created')
        sys.exit(-1)
    output_TOA.SetGeoTransform(referenceTrans)
    output_TOA.SetProjection(referenceProj)
    
    refl_TOA = refl_TOA.T        
    for i, image in enumerate(refl_TOA, 1):
        output_TOA.GetRasterBand(i).WriteArray(image)
        output_TOA.GetRasterBand(i).SetNoDataValue(0)
    
    output_TOA = None    
            
    return output_TOA



# This creates interger id attributes in shapefile based on class names, splits the training data into trainig and test subsets and then saves them to a raster file
def rasterize_train_data(ref_ras, shp_file, csv_class_name_id, train_file_with_ids,  train_out_ras):
    # Create id in shp corresponding to class name in shapefile
    gdf = gpd.read_file(shp_file)
    class_names = gdf['Class name'].unique()
    print('class names', class_names)
    class_ids = np.arange(class_names.size) + 1
    print('class ids', class_ids)
    df = pd.DataFrame({'label':class_names, 'id':class_ids})
    df.to_csv(csv_class_name_id)
    print('gdf without ids', gdf.head())
    gdf['id'] = gdf['Class name'].map(dict(zip(class_names, class_ids)))
    print('gdf with ids', gdf.head())

    gdf.to_file(train_file_with_ids)
   
    # Open reference raster
    in_ras = gdal.Open(ref_ras)
    
    train_ds = ogr.Open(train_file_with_ids)
    train_lyr = train_ds.GetLayer()

  
    # Create new training raster layer
    driver = gdal.GetDriverByName('GTiff')
    target_ds = driver.Create(train_out_ras, in_ras.RasterXSize, in_ras.RasterYSize, 1, gdal.GDT_UInt16)
    target_ds.SetGeoTransform(in_ras.GetGeoTransform())
    target_ds.SetProjection(in_ras.GetProjection())
       
    # Rasterise training points
    options = ['ATTRIBUTE=id']
    gdal.RasterizeLayer(target_ds, [1], train_lyr, options=options)
  
    # Check generated training data raster and display basic stats, set no data value to 0
    data = target_ds.GetRasterBand(1).ReadAsArray()
    print('Train', 'min', data.min(), 'max', data.max(), 'mean', data.mean())
    target_ds.GetRasterBand(1).SetNoDataValue(0)
    
    # Save raster file
    target_ds = None



def fit_RF(img, trai, out_ras_):
    gdal.UseExceptions()
    gdal.AllRegister()

    # Read in image to be classified
    img_ds = gdal.Open(img)
    img = np.zeros((img_ds.RasterYSize, img_ds.RasterXSize, img_ds.RasterCount), \
               gdal_array.GDALTypeCodeToNumericTypeCode(img_ds.GetRasterBand(1).DataType))   
        
    for i in range(1, img_ds.RasterCount + 1):
    # set the nodata value of the band
        img_ds.GetRasterBand(i).SetNoDataValue(0)
    
    # Read each band into an array
    for b in range(img.shape[2]):
        img[:, :, b] = img_ds.GetRasterBand(b + 1).ReadAsArray()
    print('The sat image data matrix is sized: {sz}'.format(sz=img.shape))    
        
    # Read in training data raster
    trai_ds = gdal.Open(trai)        
    trai = trai_ds.GetRasterBand(1).ReadAsArray().astype(np.uint8)
    # Check NoData value = 0
    print(trai_ds.GetRasterBand(1).GetNoDataValue())
    print('The training data matrix is sized: {sz}'.format(sz=trai.shape))
    
    
    # Display satellite image and training data
    plt.subplot(121)
    plt.imshow(img[:, :, 0], cmap=plt.cm.Greys_r)
    plt.title('RS image - first band')

    plt.subplot(122)
    plt.imshow(trai, cmap=plt.cm.Spectral)
    plt.title('Training data')
    
    # Find how many non-zero entries we have -- i.e how many training data samples?
    n_samples = (trai > 0).sum()
    print('We have {n} training samples'.format(n=n_samples))


    # What are the classficiation labels?
    labels = np.unique(trai[trai>0])
    print('The training data includes {n} classes: {classes}'.format(n=labels.size, classes=labels))
        
    # We will need a 'X' matrix containing our features - data in training matrix rows and a 'y' array containing our labels - cols, these will have n sample rows
    y = trai[trai>0]
    x = img[trai>0, :]

    print('Our x (training input samples) matrix is sized: {sz}'.format(sz=x.shape))
    print('Our y (target values) matrix is sized: {sz}'.format(sz=y.shape))
        
    # Initialise Random Forest Classifier, n_jobs = -1 utilises all available cores
    rf = RandomForestClassifier(oob_score=True, n_jobs=-1, random_state=42)

    # Fit our model to training data
    rf = rf.fit(x,y)
    
    print(rf.get_params())
    
    # Check the "Out-of-Bag" (OOB) prediction score
    print('OOB prediction of accuracy from first RF setup is: {oob}%'.format(oob=rf.oob_score_ * 100))
        
    # Show the band importance:
    bands = range(1,img_ds.RasterCount+1)

    for b, imp in zip(bands, rf.feature_importances_):
        print('Band {b} importance: {imp}'.format(b=b, imp=imp))
    
    # Take our full image and reshape it into a long 2d array (nrow * ncol, nband) for classification
    # Change to img.shape[2]-1 when using 4 bands (stacked NDWI + RGB)
    new_shape = (img.shape[0] * img.shape[1], img.shape[2])
    img_as_array = img[:, :, :img.shape[2]].reshape(new_shape)
    print('Reshaped from {o} to {n}'.format(o=img.shape, n=img_as_array.shape))

    # Now predict for each pixel
    prediction_img = rf.predict(img_as_array)
    
    # Reshape back into original size
    prediction_img = prediction_img.reshape(img[:, :, 0].shape)
    print('Reshaped back to {}'.format(prediction_img.shape))

    
    # Get info to save output raster from input raster
    geo = img_ds.GetGeoTransform()
    proj = img_ds.GetProjectionRef()

    ncol = img.shape[1]
    nrow = img.shape[0]

    drv = gdal.GetDriverByName('GTiff')
    
    prediction_img.astype(np.float16)
    
    out_ras_prediction = drv.Create(out_ras_, ncol, nrow, 1, gdal.GDT_UInt16)
    out_ras_prediction.SetProjection(proj)
    out_ras_prediction.SetGeoTransform(geo)       
    out_ras_prediction.GetRasterBand(1).WriteArray(prediction_img)
    out_ras_prediction = None

    return 


 
# Function for salt and pepper cleanup
# Mark lakes by using shapefile
def cleanup(in_ras, sieved_ras, fn_clip, fn_in, fn_poly):
    """
    fn_clip: needs to be a string, not a path
    fn_in: needs to be a string, not a path
    fn_poly: needs to be a string, not a path
    
    """
    raster_original = gdal.Open(in_ras)  
    raster_band_original = raster_original.GetRasterBand(1)
    raster_array_original = raster_original.GetRasterBand(1).ReadAsArray() 
    
    geo = raster_original.GetGeoTransform()
    proj = raster_original.GetProjectionRef()  
    
    ncol = raster_array_original.shape[1]
    nrow = raster_array_original.shape[0]
    
    drv = gdal.GetDriverByName('GTiff')
    
    sieved_raster = drv.Create(sieved_ras, ncol, nrow, 1, gdal.GDT_UInt16)
    sieved_raster.SetProjection(proj)
    sieved_raster.SetGeoTransform(geo)       
    a = sieved_raster.GetRasterBand(1)    
    
    gdal.SieveFilter(srcBand=raster_band_original, maskBand=None, dstBand=a, threshold=2)
   
    a = None
    sieved_raster = None
         
    # Mark lake areas with lake shapefile, lake class will be 0
    gdal.Warp(fn_clip, fn_in, format="GTiff", cutlineDSName=fn_poly, cropToCutline=True)  
    
    fn_clip = None
    
    return    

 

# Calculate area covered by each water body based on class pixels
def calc_wb_areas(in_raster, wbs_dict, txt_file):
    # Open reprojected raster
    raster = gdal.Open(in_raster)  
    raster_array = raster.GetRasterBand(1).ReadAsArray()  

    # Create dictionary with class name : amount of pixels
    dicts = {}
    for key, value in wbs_dict.items():
        dicts[key] = len(raster_array[raster_array==value])
           
    # Get pixel size and calculate area
    gt = raster.GetGeoTransform()
    pixel_area = gt[1] * (-gt[5])
    
    total_area = sum(dicts.values()) * pixel_area

    areas = []
    # Multiply each total in list with pixel_area to get area covered by wb in km2
    # Create dictionary of classes and corresponding area in km2     
    wb_areas = {}
    for key, value in dicts.items():
       wb_areas[key] = round(value*pixel_area*10**-6,2)
    
    areas.append(wb_areas)
            
    # Create dictionary of classes and corresponding area in % of total study area   
    wb_areas_perc = {}
    for key, value in dicts.items():
       wb_areas_perc[key] = round((value*pixel_area)/total_area, 2) * 100    
    
    areas.append(wb_areas_perc)
    
    with open(txt_file, 'w') as f:
        print(areas, file=f)    
       
    
    return wb_areas, wb_areas_perc




