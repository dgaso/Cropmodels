
# Copyright (c) 2020 Wageningen-UR
# Deborah Gaso Melgar and Allard de Wit, June 2020

import sys, os
this_dir = os.getcwd()
up_dir = os.path.dirname(this_dir)
if not up_dir in sys.path:
    sys.path.append(up_dir)

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yaml
import pandas as pd
import numpy as np
import nlopt

from pcse.base import ParameterProvider
from pcse.fileinput import ExcelWeatherDataProvider, PCSEFileReader
import rasterio

from optimizer import ObjectiveFunctionCalculator
import config


def make_agromanagement(year):
    """Returns the proper crop agromanagement dates for given campaign year

    :param year: The harvest year of the cropping season
    :return: the PCSE agromanagement structure
    """

    campaign_dates = config.agromanagement[year]
    agromanagement = """
    - {campaign_start_date}:
        CropCalendar:
            crop_name: Soybean 
            variety_name: Soybean 
            crop_start_date: {crop_start_date}
            crop_start_type: sowing
            crop_end_date: {crop_end_date}
            crop_end_type: harvest
            max_duration: 300
        TimedEvents:
        StateEvents:
        """
    agro = yaml.safe_load(agromanagement.format(**campaign_dates))
    return agro

def read_Sentinel2_timeseries(col, row):
    """Reads a time series of sentinel2 LAI estimates for given pixel,

    :param col: the column number of S2 timeseries to read
    :param row: the row number of the S2 timeseries to read
    :return: a dateframe with dates and LAI values
    """
    # Read images and create a df
    
    im=rasterio.open(os.path.join(config.s2_dir,config.farmer))    
    imarray=im.read()  
    ci = np.array(imarray) 

     
    #get the date as a tuple
    date=im.descriptions
    #convert tuple into list
    date_l=list(date)
  
    dates=[]
    for fname in date_l:
        name, d, *rest = fname.split("_")
        dates.append(d.split("T")[0])
    days = pd.to_datetime(dates)

    CItimeseries = np.isnan(imarray[:,row,col])
    if all(flag == True for flag in CItimeseries[:]) == True:
        raise RuntimeError("No data available in S2 timeseries!")
        
    LAIpixel = pd.DataFrame({"day": days,
                             "LAI": (ci[:, row, col] ** config.CI_coefficient) / config.CI_offset})
    LAIpixel.index = LAIpixel.day
        
    return LAIpixel


def start_optimizer(objfunc_calculator):
    """Starts and and returns the optimizer object
    """
    # Start optimizer with the SUBPLEX algorithm for one parameter
    opt = nlopt.opt(nlopt.LN_SBPLX, len(config.selected_parameters))
    # Assign the objective function calculator
    opt.set_min_objective(objfunc_calculator)
    # lower bounds of parameters values
    lbounds = [config.parameter_settings[pname]["minimum"] for pname in config.selected_parameters]
    opt.set_lower_bounds(lbounds)
    # upper bounds of parameters values
    ubounds = [config.parameter_settings[pname]["maximum"] for pname in config.selected_parameters]
    opt.set_upper_bounds(ubounds)
    # the initial step size to compute numerical gradients
    stepsize = [config.parameter_settings[pname]["stepsize"] for pname in config.selected_parameters]
    opt.set_initial_step(stepsize)
    # Maximum number of evaluations allowed
    opt.set_maxeval(config.nlopt_maxeval)
    # Relative tolerance for convergence
    opt.set_ftol_rel(config.nlopt_ftol_rel)
    # Start the optimization with the first guess

    return opt


def optimize_one_pixel(year, col, row, silent=False):
    """Runs the optimization for given year and pixel.
    """

    # Weather data for Uruguay
    wdp = ExcelWeatherDataProvider(config.weather_fname)

    # Model parameters
    cropd = PCSEFileReader(config.crop_fname)
    soild = PCSEFileReader(config.soil_fname)          
    params = ParameterProvider(cropdata=cropd, soildata=soild,sitedata={})
 
    # Here we define the agromanagement for soybean
    agro = make_agromanagement(year)

    # Read images and create a df
    LAI_pixel = read_Sentinel2_timeseries(col, row)

    # Start the objective function calculator
    objfunc_calculator = ObjectiveFunctionCalculator(params, wdp, agro, LAI_pixel)
    defaults = [ cropd["WUE"],soild["RDMAX"],cropd["FNTR"],cropd["initLAI"] ]# ,
    error = objfunc_calculator(defaults)
    
    if not silent:
        print("Objective function value with default parameters (%s: %s): %s" %
              (config.selected_parameters, defaults, error))

    # Start and run the optimizer
    opt = start_optimizer(objfunc_calculator)
    firstguess = [config.parameter_settings[pname]["default"] for pname in config.selected_parameters]
    x = opt.optimize(firstguess)
    harvest_yield=objfunc_calculator.df_simulations['YIELD'][objfunc_calculator.df_simulations.index[-1]]
    deficit_veg=objfunc_calculator.df_simulations['CWDv'][objfunc_calculator.df_simulations.index[-1]]
    deficit_rep=objfunc_calculator.df_simulations['CWDr'][objfunc_calculator.df_simulations.index[-1]]

    #save errors from opt
    errlai=objfunc_calculator.err
    rrmse_lai = objfunc_calculator.rrmse
    LAI_max = objfunc_calculator.LAImx

   
    if not silent:
        print("\noptimum for parameters %s at %s" % (config.selected_parameters, x))
        print("minimum value = ",  opt.last_optimum_value())
        print("result code = ", opt.last_optimize_result())
        print("With %i function calls" % objfunc_calculator.n_calls)

    # Generate output figures
#    fig, axes = plt.subplots(figsize=(12, 20), nrows=3, ncols=1)
#    objfunc_calculator.df_simulations.YIELD.plot(ax=axes[0])
#    objfunc_calculator.df_simulations.TDM.plot(ax=axes[1])
#    objfunc_calculator.df_simulations.LAI.plot(ax=axes[2])
#    objfunc_calculator.df_observations.LAI.plot(ax=axes[2], color='r', marker='o')
#    fig.autofmt_xdate()
#    fname_figure = os.path.join(config.this_dir, "output", f"optimization_results_{year}_{col}_{row}.png")
#    fig.savefig(fname_figure)
#    plt.close("all")


    return harvest_yield ,x[0],x[1],x[2], x[3], LAI_max,deficit_veg,deficit_rep

if __name__ == "__main__":
    optimize_one_pixel(year=2020, col=4, row=4)
