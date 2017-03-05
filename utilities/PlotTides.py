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
    
# Optional Changes   
    2. Not Implemented: import addtional files to evaulate sensitivity

###############################################################################
"""

import pandas as pd
import matplotlib.pyplot as plt
import os

#----Required Variables

root_dir = r'C:\Users\sml\Desktop\tz' #--Root Directory
prediction_file = 'feb_cbbt_95.out'   # output file from ntp4 (station data by month)
start = pd.datetime(2017, 2,1, 0)


# Optional Changes  
prediction_file_2 = 'feb_cbbt_95_2.out' # output file from ntp4 (differnt adcirc run for same station)


#----------------------------------RUN SCRIPT--------------------------------#
tides  = os.path.join(root_dir, 'predictions')
tides  = os.path.join(tides, 'outputs')
tides  = os.path.join(tides, prediction_file)

startline = 6
heights = []

j=0
with open(tides) as f:
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
x = df.index
y = df['Water Surface']

fig, ax = plt.subplots()
ax.plot(x, y, color = 'b')


'''

#USE THIS SECTION FOR COMPARING AGAINST NOAA CONTROL FILE or ADCIRC RESULTS
x = df.index
y = df['Water Surface']

#--Copied from above
heights = []
j=0
with open(noaa_file) as f:
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
x2 = df.index
y2 = df['Water Surface']
ax.plot(x2, y2, color = 'r')
ax.grid()


#fig.savefig('feb.png', figsize=(3, 11), dpi=600, facecolor='w', edgecolor='k')
#fig.show()
'''