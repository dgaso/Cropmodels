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

    relevant_pixels = []; index=[] 
    i=0
    for (row, col) in product(range(imarray.shape[1]), range(imarray.shape[2])):
        i+=1
        ar_nan=np.isnan(imarray[:,row,col])        
        if all(flag == True for flag in ar_nan[:]) == True:
            continue
        relevant_pixels.append((col, row))
        index.append(i-1)
              
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
    pixels_for_DA, index,meta, im = gather_all_pixels()
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
    
    l1=[]
    for i in range(l_len):
        for j in range(l_item_len):
            v1=l[i][0]         
    
        l1.append(v1)        

    leng=im.shape[0]*im.shape[1]
    y_pred = np.full(leng, np.NaN)
    
    results_y = np.array(l1)
       
    for index, value in zip(index, results_y):
        y_pred[index] = value
    result_yield=np.reshape(y_pred[:], (im.shape[0],im.shape[1]))    

    #Update meta from stack tif (index) to write a single tif
    meta.update({'dtype': 'float64',
                 'count': 1})  
    with rasterio.open('Field1_yield.tif', 'w', **meta) as outds:
        outds.write(np.expand_dims(result_yield,0))
    
    return 

if __name__ == '__main__':
    optimizer_with_mp()
    

