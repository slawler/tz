#!/usr/bin/env python 3.5.2
# -*- coding: utf-8 -*-
"""
@author: slawler@dewberry.com
Created on Sun Feb 19 09:35:33 2017


###############################################################################

PlotTides.py                   USAGE NOTES


In User Inputs, make the following changes:
    (Assumes Repo cloned from github with directory structure in tact)
    
# Required Variables    
    1. Assign directory path, filenames & Start date
    2. assign plot start, stop or set as entire month from tide table
    
# Optional Changes   
    2. Not Implemented: import addtional files to evaulate sensitivity

###############################################################################
"""

import pandas as pd
import matplotlib.pyplot as plt
import os

#----Required Variables

#root_dir = r'C:\Users\sml\Desktop\tz' #--Root Directory
root_dir = '/Users/slawler/Desktop/tz'
p1 = '90d.out'   # output file from ntp4 (station data by month)
p2 = '410d.out' # ntp4 output 2 (dif run for same station)

# Commented out below
#p3 = '8573364_feb17.out' #ntp4 output 3 --or---NOAA official tides


#--Set pstop = 0 to plot complete monthly time series
start = pd.datetime(2017, 2,1, 0)
save_fig_name = p1[:7] + '_' + p2[:7]
pstart, pstop = 24, 0

#----------------------------------RUN SCRIPT--------------------------------#
tides  = os.path.join(root_dir, 'predictions')
tides  = os.path.join(tides, 'outputs')

tides_p1  = os.path.join(tides, p1)
tides_p2  = os.path.join(tides, p2)
#tides_p3  = os.path.join(tides, p3)

def tides2df(start, tide_table):
    startline = 6
    heights = []
    j=0
    with open(tide_table) as f:
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
    df = pd.DataFrame(heights, dtype = float, 
                      columns = ['Water Surface'], index = dtm) 
    return df         

    
#--Tidal Prediction from p1    
df = tides2df(start, tides_p1)
x = df.index
y1 = df['Water Surface']


#--Tidal Prediction from p2   
df = tides2df(start, tides_p2)
y2 = df['Water Surface']


#--Tidal Prediction from p3   
#df = tides2df(start, tides_p3)
#y3 = df['Water Surface']


#--Initialize Plot
fig, ax = plt.subplots()

#--If plot interval not specified, plot cpomplete time series
if pstop == 0: 
    pstart = 0
    pstop = len(df)
    

#---To avoid confusion, clear df
df=None
    
#ax.plot(x, y, color = 'b')
ax.plot(x[pstart:pstop], y1[pstart:pstop], color = 'b')
ax.plot(x[pstart:pstop], y2[pstart:pstop], color = 'r')
#ax.plot(x[pstart:pstop], y3[pstart:pstop], color = 'g')

ax.legend([p1[:7],p2[:7]])
ax.grid()

'''
#--Initialize Difference Plot
fig_dif, ax_dif = plt.subplots()
dif = y2[pstart:pstop] - y1[pstart:pstop]
ax_dif.plot(x[pstart:pstop],dif , color = 'black')
ax_dif.legend(['Residuals'])
ax_dif.grid()


#fig.savefig('{}.png'.format(save_fig_name), figsize=(3, 11), dpi=600, facecolor='w', edgecolor='k')
#fig.show()
'''