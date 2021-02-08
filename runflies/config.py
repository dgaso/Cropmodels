# ---------------------------------------------------------------------#
# Configuration file for running the Campbell-Diaz model and optimizer #
#                                                                      #
# Allard de Wit and Deborah Gaso Melgar, Wageningen 2020               #
#----------------------------------------------------------------------#
import sys, os
import yaml

this_dir = os.path.dirname(__file__)
top_dir = os.path.dirname(this_dir)
data_dir = os.path.join(top_dir, "data")


# ------- SETTINGS FOR CAMPBELL-DIAZ MODEL --------
# Weather data
#weather_fname = os.path.join(data_dir, "Palmar.xlsx")
weather_fname = os.path.join(data_dir, "uy_le2.xlsx")
#weather_fname = os.path.join(data_dir, "Dolores.xlsx")
#weather_fname = os.path.join(data_dir, "VillaDarwin.xlsx") 
#weather_fname = os.path.join(data_dir, "Mendoza.xlsx")
#weather_fname = os.path.join(data_dir, "Trinidad.xlsx")
#weather_fname = os.path.join(data_dir, "25Agosto.xlsx")
#weather_fname = os.path.join(data_dir, "ElAguila.xlsx")
#weather_fname = os.path.join(data_dir, "Cardona.xlsx")
# model parameters
crop_fname = os.path.join(data_dir, "wofost_soybean_parameters.dat")
soil_fname = os.path.join(data_dir, "soil_parameters.dat")

# agromanagement
agromanagement = """
2014:
    campaign_start_date: 2013-10-30
    crop_start_date: 2013-11-15
    crop_end_date: 2014-05-10
2015:
    campaign_start_date: 2014-10-30
    crop_start_date: 2014-11-15
    crop_end_date: 2015-05-10
2016:
    campaign_start_date: 2015-10-30
    crop_start_date: 2015-11-15
    crop_end_date: 2016-05-10
2017:
    campaign_start_date: 2016-10-30
    crop_start_date: 2016-11-15
    crop_end_date: 2017-05-10
2018:
    campaign_start_date: 2017-12-20
    crop_start_date: 2017-12-26
    crop_end_date: 2018-05-10
2019:
    campaign_start_date: 2018-10-30
    crop_start_date: 2018-11-15
    crop_end_date: 2019-05-20
2020:
    campaign_start_date: 2019-10-30
    crop_start_date: 2019-11-15
    crop_end_date: 2020-04-19
"""
agromanagement = yaml.safe_load(agromanagement)


# --- SETTINGS FOR READING SENTINEL2 DATA ---
s2_dir = os.path.join(data_dir, "sentinel2")
# parameters for converting CI to LAI
CI_coefficient = 0.898
CI_offset = 0.904
# No data value
S2_nodata = -9999


# --- SETTINGS FOR THE NLOPT OPTIMIZER ----
parameter_settings = """

initLAI:
    default: 0.15
    minimum: 0.08
    maximum: 0.3
    stepsize: 0.01

RDMAX:
    default: 0.9
    minimum: 0.6
    maximum: 1.2
    stepsize: 0.02

p_fc:
    default: 0.30
    minimum: 0.28
    maximum: 0.40
    stepsize: 0.01

fr_ntr:
    default: 2
    minimum: 1
    maximum: 5
    stepsize: 0.2



"""
parameter_settings = yaml.safe_load(parameter_settings)
selected_parameters = ["initLAI","RDMAX","p_fc","fr_ntr"]# 
#selected_parameters = ["rdmax","p_pwp","p_fc"]# 
#selected_parameters = ["p_fc","rdmax","p_pwp"] 
#selected_parameters = ["p_fc","p_pwp","rdmax"] 
#selected_parameters = ["p_pwp","p_fc","rdmax"] 
#selected_parameters = ["p_pwp","rdmax","p_fc"] 


nlopt_maxeval = 250
nlopt_ftol_rel = 0.001

#p_pwp:
#    default: 0.20
#    minimum: 0.1
#    maximum: 0.25
#    stepsize: 0.05
#ntr:
#    default: 0.065
#    minimum: 0.020
#    maximum: 0.099
#    stepsize: 0.01
#lntr:
#    default: 0.04
#    minimum: 0.01
#    maximum: 0.08
#    stepsize: 0.01
#dwr:
#    default: 40
#    minimum: 30
#    maximum: 50
#    stepsize: 1
