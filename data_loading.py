# -*- coding: utf-8 -*-
"""
The following code is part of the manuscript: 
Written by M. Martyn Rosco
"""

import os
import numpy as np
import pandas as pd
import openpyxl
import datetime as dt
import math

# Create a new column with data from other columns based on condition, condit1 == condit2
def newcol(dataframe_, new, col, condit1, condit2, yes, no):
    for i in range(len(new)):
       dataframe_[new[i]] = np.where((condit1 ==  dataframe_[col]) | (dataframe_[col] == condit2) , dataframe_[yes[i]] , dataframe_[no[i]])
    return dataframe_

# Needs work to be more general. Loop that generates conditions based on how many given and choices string list as input param
def wb_type(dataframe):
    conditions = [dataframe.index.str.contains('L1', na=False), dataframe.index.str.contains('L2|L3', na=False), dataframe.index.str.contains('P', na=False), dataframe.index.str.contains('R', na=False),\
        dataframe.index.str.contains('S', na=False), dataframe.index.str.contains('T', na=False)]
    choices = ['Lake', 'Lake', 'Pond', 'Fluvial', 'Snow', 'Flood']
    dataframe['waterbody'] = np.select(conditions, choices)
    return dataframe

def year(dataframe):
    conditions = [dataframe.index.str.contains('B', na=False), dataframe.index.str.contains('C', na=False)]
    choices = ['2016', '2017']
    df['year'] = np.select(conditions, choices)
    return dataframe

def lake_type(dataframe):
    conditions = [dataframe.index.str.contains('L1', na=False), dataframe.index.str.contains('L2', na=False), dataframe.index.str.contains('L3', na=False)]
    choices = ['LakeT', 'Lake1', 'Lake2']
    df['Lake'] = np.select(conditions, choices)
    return dataframe

def fill_air_gas(dataframe_main, dataframe_fill, columns):
    df_means = dataframe_fill.groupby('Date')[columns].mean()
    #dataframe_main = dataframe_main.reset_index().set_index('Date')
    dataframe_filled = dataframe_main.reset_index().set_index('Date').combine_first(df_means)
    dataframe_main = dataframe_filled.reset_index().set_index('Site')
    return dataframe_main

def wind(dataframe_main, dataframe_fill, columns, new_col_names):
    dataframe_main = dataframe_main.reset_index()
    df_means = dataframe_fill.groupby('Date')[columns].mean()
    df_means.index = pd.to_datetime(df_means.index, format='%Y-%m-%d')    
    dataframe_main = pd.merge(dataframe_main, df_means, how='left', on='Date')
    dataframe_main = dataframe_main.set_index('Site')
    dataframe_main.rename(columns=new_col_names, inplace=True)
    
    # Convert to wind at 10m from 5m (??)
    H = 5        
    a = 1/7
    dataframe_main['wind speed'] = dataframe_main['wind speed'] * math.pow((10/H), a)
    
    return dataframe_main    
    


# Set directory
os.chdir(r'')
# Read in main spreadhseet with dissolved conc data
df = pd.read_excel('Data_Readin_2023_Tiski.xlsx', index_col=0, engine='openpyxl', skiprows=1)

# Read in spreadsheet with 2017 air gases
df_air_17 = pd.read_excel('2017_air_gases.xlsx', index_col=0, engine='openpyxl', skiprows=1)

cols = ['CH4 air', 'N2O air']
df = fill_air_gas(df, df_air_17, cols)

# Convert CH4 concentration from nmol L-1 to micromol L-1
df['diss CH4'] = df['diss CH4'] * 10**-3

# Make a new col with T in K
df['T - K'] = df['T']+273.15

# Make a new column with lake type
lake_type(df)

#  Make new columns with inland water system names
wb_type(df)


# Make year column
year(df)


# Get wind data -------------------------------------------------------------------------------------------
fname = r'...wind.csv'
df_wind = pd.read_csv(fname, delimiter=(';'))
df_wind.drop(labels='X', axis=1, inplace=True)


df_wind.dtypes
# Get only 2016 and 2017 data
# Format TIMESTAMP1 column into a date and time column. Delte the text UTC from TIMESTASMP1
df_wind[['Date', 'time', 'UTC text']] = df_wind['TIMESTAMP1'].str.split(' ', 3, expand=True)
df_wind.drop(labels='UTC text', axis=1, inplace=True)
df_wind['time'] =  pd.to_datetime(df_wind['time'], format='%H:%M:%S')
df_wind['time'] = df_wind['time'].dt.time

df_wind['Date'] =  pd.to_datetime(df_wind['Date'], format='%Y-%m-%d')
df_wind['Date'] = df_wind['Date'].dt.date
df_wind.dtypes
df_wind.drop(labels='TIMESTAMP1', axis=1, inplace=True)


# Make a mask for time, only select wind during day time
starttime = pd.to_datetime('08:30:00').time()
endtime = pd.to_datetime('17:00:00').time()
mask_time = (df_wind['time'] > starttime) & (df_wind['time'] <= endtime)

# Make a mask for dates 2016
startdate = pd.to_datetime('2016-07-30').date()
enddate = pd.to_datetime('2016-08-10').date()
mask_16 = (df_wind['Date'] > startdate) & (df_wind['Date'] <= enddate)
df_wind_16 = df_wind.loc[mask_16]
df_wind_16 = df_wind_16.loc[mask_time]


# Make a mask for dates 2017
startdate = pd.to_datetime('2017-07-07').date()
enddate = pd.to_datetime('2017-07-18').date()
mask_17 = (df_wind['Date'] > startdate) & (df_wind['Date'] <= enddate)
df_wind_17 = df_wind.loc[mask_17]
df_wind_17 = df_wind_17.loc[mask_time]

# Concatenate 2016 and 2017 wind dataframes
df_wind_measured = pd.concat([df_wind_16, df_wind_17], axis=0)

# Assigned wind speed and direction as day averages from wind dataframe to main dataframe on corresponding dates
rename_keys = {'WS_1_1_1':'wind speed', 'WD_1_1_1':'wind direction'}
df = wind(df, df_wind_measured, ['WS_1_1_1', 'WD_1_1_1'], rename_keys)









