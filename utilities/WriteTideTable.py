#!/usr/bin/env python 3.5.2
# -*- coding: utf-8 -*-
"""
@author: slawler@dewberry.com
Created on Sun Feb 19 09:35:33 2017





!!!!!------UNDER CONSTRUCTION---Script Not Ready for use-----!!!!!









###############################################################################

WriteTideTable.py                 USAGE NOTES



In User Inputs, make the following changes:
    (Assumes Repo cloned from github with directory structure in tact)
    
# Required Variables    
    1. UNDER CONSTRUCTION


###############################################################################
"""
import pandas as pd
import os



#----Required Variables

root_dir = r'C:\Users\sml\Desktop\tz' #--Root Directory
prediction_file = 'feb_cbbt_95.out'   # output file from ntp4 (station data by month)

# Optional Changes  
prediction_file_2 = 'feb_cbbt_95_2.out' # output file from ntp4 (differnt adcirc run for same station)


#----------------------------------RUN SCRIPT--------------------------------#
tides  = os.path.join(root_dir, 'predictions')
tides  = os.path.join(tides, 'outputs')
tides  = os.path.join(tides, prediction_file)


months = {'feb':2, 'mar':3, 'apr':4, 'may':5, 'jun':6,'jul':7,
          'aug':8, 'sep':9,'oct':10,'nov':11}
year=17

master=None #Initialize outside of loop

for month in months:
    
    station_file = 'ANAD2.{}{}'.format(month, year)
    start = pd.datetime(2000+year, months[month],1, 0)
    
    adcirc_file = os.path.join(tide_dir, '{}'.format(station_file))
    #noaa_file = os.path.join(tide_dir, 'noaa_prediction_feb_dc.txt')
    startline = 6
    heights = []
    
    j=0
    with open(adcirc_file) as f:
        for i in range(68):  
            try:
                if i < startline:
                    line = f.readline()
                else:
                    line = f.readline().split()
                    day = line[0]
                    h = line[1:]
                    heights = heights + h
                    
            except:
                break
    
    dtm = pd.date_range(start = start, freq = 'H' , periods= len(heights))
    df = pd.DataFrame(heights, dtype = float, columns = ['Water Surface'], index = dtm)          
    
    
   
    if isinstance(master, pd.DataFrame):
        master = pd.concat([master, df])
    else:
        master = df.copy()
    print(len(master))
    print(master.head())    