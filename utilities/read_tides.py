# -*- coding: utf-8 -*-
"""
python 3
slawler@dewberry.com

"""
import pandas as pd
import matplotlib.pyplot as plt
import os

tide_dir = r'C:\Users\sml\Desktop\tz\validation'

stationname = '8594900'

adcirc_file = os.path.join(tide_dir, '{}.txt'.format(stationname))
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

start = pd.datetime(2017, 2,1, 0)
dtm = pd.date_range(start = start, freq = 'H' , periods= len(heights))
df = pd.DataFrame(heights, dtype = float, columns = ['Water Surface'], index = dtm)          

fig, ax = plt.subplots()
ax.plot(x, y, color = 'b')


'''

#USE THIS SECTION FOR COMPARING AGAINST NOAA CONTROL FILE
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