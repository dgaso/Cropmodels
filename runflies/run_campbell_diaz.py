# -*- coding: utf-8 -*-
# Copyright (c) 2020 Wageningen-UR
# Deborah Gaso Melgar and Allard de Wit, June 2020

import sys, os
this_dir = os.getcwd()
up_dir = os.path.dirname(this_dir)
if not up_dir in sys.path:
    sys.path.append(up_dir)
    
import datetime as dt
import yaml
import pandas as pd
import numpy as np
import matplotlib
#matplotlib.style.use("ggplot")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

from pcse.fileinput import ExcelWeatherDataProvider, PCSEFileReader
from pcse.base import ParameterProvider
from campbell_diaz.model import CampbellDiazModel
import config


from matplotlib import rc
rc('mathtext', default='regular')

def make_agromanagement(year):
    
    start_dates = dict()

    start_dates[2014] = dict(campaign_start_date=dt.date(year-1, 11,10),
                             crop_start_date=dt.date(year-1,11,10),
                             crop_end_date=dt.date(2014, 4, 20))
    start_dates[2015] = dict(campaign_start_date=dt.date(year-1, 11,13),
                             crop_start_date=dt.date(year-1,11,13),
                             crop_end_date=dt.date(2015, 4, 25))
    start_dates[2016] = dict(campaign_start_date=dt.date(year-1, 11,17),
                         crop_start_date=dt.date(year-1,11,17),
                         crop_end_date=dt.date(2016, 4, 20))
    start_dates[2017] = dict(campaign_start_date=dt.date(year-1, 11,11),
                     crop_start_date=dt.date(year-1,11,11),
                     crop_end_date=dt.date(2017, 4, 20))
    start_dates[2018] = dict(campaign_start_date=dt.date(year-1, 11,23),
                         crop_start_date=dt.date(year-1,11,23),
                         crop_end_date=dt.date(2018, 4, 20))
    start_dates[2019] = dict(campaign_start_date=dt.date(year-1, 12,23),
                         crop_start_date=dt.date(year-1,12,23),
                         crop_end_date=dt.date(2019, 5, 20))
    start_dates[2020] = dict(campaign_start_date=dt.date(year-1, 11,14),
                         crop_start_date=dt.date(year-1,11,14),
                         crop_end_date=dt.date(2020, 4, 19))

    campaign_dates = start_dates[year]
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

def main():
    
    agro=make_agromanagement(2014)

    weather_data = ExcelWeatherDataProvider(config.weather_fname)

    cropd = PCSEFileReader(config.crop_fname)
    soild = PCSEFileReader(config.soil_fname)
    params = ParameterProvider(cropdata=cropd, soildata=soild,sitedata={})

    model = CampbellDiazModel(params, weather_data, agro)
    model.run_till_terminate()
    output=model.get_output()
    
    df = pd.DataFrame(model.get_output()).set_index("day")

#    # Plot results
    fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(10,5))
    for key, axis in zip(df.columns, axes.flatten()):
        df[key].plot(ax=axis, title=key)
    fig.autofmt_xdate()
    fig.savefig(os.path.join(this_dir, "output", "Campbell_soybean.png"))

    csv_fname = os.path.join(this_dir, "output", "Campbell_soybean.csv")
    df.to_csv(csv_fname, header=True)
    
    
    return model,output


if __name__ == "__main__":
    model,output=main()
    
 