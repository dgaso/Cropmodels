# -*- coding: utf-8 -*-
# Copyright (c) 2020 Wageningen-UR
# Deborah Gaso Melgar and Allard de Wit, June 2020
import sys, os
import time
from itertools import product

from multiprocessing import Pool
import rasterio

from run_optimization_S2 import optimize_one_pixel
import config

import numpy as np


def gather_all_pixels():
    """Gathers the cols/rows for all pixels with data
    :return: a list of [(col, row), (col, row),...]
    """
    im=rasterio.open(os.path.join(config.s2_dir,config.farmer))    
    imarray=im.read()  
    ci = np.array(imarray) 
    meta = im.meta.copy()

    relevant_pixels = []; index=[] ; index1=[]; index2=[]; index3=[]; index4=[];index5=[];index6=[];index7=[] 
    i=0
    for (row, col) in product(range(imarray.shape[1]), range(imarray.shape[2])):
        i+=1
        ar_nan=np.isnan(imarray[:,row,col])        
        if all(flag == True for flag in ar_nan[:]) == True:
            continue
        relevant_pixels.append((col, row))
        index.append(i-1)
        index1.append(i-1)
        index2.append(i-1)
        index3.append(i-1)
        index4.append(i-1)
        index5.append(i-1)
        index6.append(i-1)
        index7.append(i-1)
        
    return relevant_pixels,index,index1,index2,index3,index4, index5,index6,index7, meta, im

def run_one_pixel(inputs):
    """Runs for one pixel with given inputs.

    Function is mainly here to handle exceptions and to unpack the inputs tuple.

    :param inputs: a tuple of type (year, (col, row))
    :return: the optimized parameters.
    """
    year, (col, row) = inputs
    try:
        r = optimize_one_pixel(year, col, row, silent=True)
    except Exception as e:
        print(f"Optimization failed on year: {year}, row: {row}, col: {col}")
        r = None
    return r


def optimizer_with_mp():
    """Runs the optimization using multiprocessing.Pool
    """
    years = [2020, ]
    pixels_for_DA, index,index1,index2, index3, index4, index5,index6,index7, meta, im = gather_all_pixels()
    relevant_years_pixels = list(product(years, pixels_for_DA))
    p = Pool(30)
    start_time = time.time()
    results = p.map(run_one_pixel, relevant_years_pixels)
    
    end_time = time.time()-start_time
    print("\n")
    print(f"Processing{len(relevant_years_pixels)} numbers took {end_time} time using multiprocessing.")
#    for result, (year, (col, row)) in zip(results, relevant_years_pixels):
#        print(f"optimized parameters for year: {year}, row: {row}, col: {col}")
#        for parname, value in zip(config.selected_parameters, result):
#            print(f" - {parname}: {value}")
    
    l=results
    
    l_len = len(l)
    print(results)
    l_item_len = len(l[1])
    
    l1=[];l2=[];l3=[];l4=[];l5=[];l6=[];l7=[];l8=[]
    for i in range(l_len):
        for j in range(l_item_len):
            v1=l[i][0]
            v2=l[i][1]
            v3=l[i][2]
            v4=l[i][3]
            v5=l[i][4]
            v6=l[i][5]
            v7=l[i][6]
            v8=l[i][7]
    
        l1.append(v1)        
        l2.append(v2)
        l3.append(v3)
        l4.append(v4)
        l5.append(v5)
        l6.append(v6)
        l7.append(v7)
        l8.append(v8)
        

    leng=im.shape[0]*im.shape[1]
    y_pred = np.full(leng, np.NaN)
    p1_pred = np.full(leng, np.NaN)
    p2_pred = np.full(leng, np.NaN)
    p3_pred = np.full(leng, np.NaN)
    p4_pred = np.full(leng, np.NaN)
    max_lai = np.full(leng, np.NaN)
    def_veg = np.full(leng, np.NaN)
    def_rep = np.full(leng, np.NaN)
    
    results_y = np.array(l1)
    results_p1 = np.array(l2)
    results_p2 = np.array(l3)
    results_p3 = np.array(l4)
    results_p4 = np.array(l5)
    results_max_lai = np.array(l6)
    results_def_veg = np.array(l7)
    results_def_rep = np.array(l8)
    
    for index, value in zip(index, results_y):
        y_pred[index] = value
    result_yield=np.reshape(y_pred[:], (im.shape[0],im.shape[1]))    

    for index1, value in zip(index1, results_p1):
        p1_pred[index1] = value
    result_p1 =np.reshape(p1_pred[:], (im.shape[0],im.shape[1]))

    for index2, value in zip(index2, results_p2):
        p2_pred[index2] = value
    result_p2 =np.reshape(p2_pred[:], (im.shape[0],im.shape[1]))

    for index3, value in zip(index3, results_p3):
        p3_pred[index3] = value
    result_p3 =np.reshape(p3_pred[:], (im.shape[0],im.shape[1]))

    for index4, value in zip(index4, results_p4):
        p4_pred[index4] = value
    result_p4 =np.reshape(p3_pred[:], (im.shape[0],im.shape[1]))
    
    for index5, value in zip(index5, results_max_lai):
        max_lai[index5] = value
    results_max_lai =np.reshape(max_lai[:], (im.shape[0],im.shape[1]))

    for index6, value in zip(index6, results_def_veg):
        def_veg[index6] = value
    results_def_veg =np.reshape(def_veg[:], (im.shape[0],im.shape[1]))
    
    for index7, value in zip(index7, results_def_rep):
        def_rep[index7] = value
    results_def_rep =np.reshape(def_rep[:], (im.shape[0],im.shape[1]))

    #Update meta from stack tif (index) to write a single tif
    meta.update({'dtype': 'float64',
                 'count': 1})
    
    with rasterio.open('Field1_yield.tif', 'w', **meta) as outds:
        outds.write(np.expand_dims(result_yield,0))
    with rasterio.open('Field1_WUE.tif', 'w', **meta) as outds:
        outds.write(np.expand_dims(result_p1,0))
    with rasterio.open('Field1_RDMAX.tif', 'w', **meta) as outds:
        outds.write(np.expand_dims(result_p2,0))
    with rasterio.open('Field1_FNTR.tif', 'w', **meta) as outds:
        outds.write(np.expand_dims(result_p3,0))
    with rasterio.open('Field1_initLAI.tif', 'w', **meta) as outds:
        outds.write(np.expand_dims(result_p4,0))
    with rasterio.open('Field1_LAImax.tif', 'w', **meta) as outds:
        outds.write(np.expand_dims(results_max_lai,0))
    with rasterio.open('Field1_defveg.tif', 'w', **meta) as outds:
        outds.write(np.expand_dims(results_def_veg,0))
    with rasterio.open('Field1_defrep.tif', 'w', **meta) as outds:
        outds.write(np.expand_dims(results_def_rep,0))
        
    
    return 

if __name__ == '__main__':
    optimizer_with_mp()
    

