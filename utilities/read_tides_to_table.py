# -*- coding: utf-8 -*-
"""
python 3
slawler@dewberry.com

"""
import pandas as pd
import os


tide_dir = r'C:\Users\sml\Desktop\nwm\station_data\harmonics'

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