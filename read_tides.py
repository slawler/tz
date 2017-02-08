# -*- coding: utf-8 -*-
"""
python 3
slawler@dewberry.com

"""
import pandas as pd
import matplotlib.pyplot as plt

tidefile = 'feb17.txt'
startline = 6
heights = []

j=0
with open(tidefile) as f:
    
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

x = df.index
y = df['Water Surface']

fig, ax = plt.subplots()
ax.plot(x, y)
fig.autofmt_xdate()

fig.savefig('feb.png', figsize=(3, 11), dpi=600, facecolor='w', edgecolor='k')
fig.show()





'''
j=0
with open(tidefile) as f:
    
    for i in range(68):
        print(i)
        
        if i < startline:
            line = f.readline()
            
        else:
            line = f.readline().split()
            day = line[0]
            h = line[1:]
            
            if day not in mh and j !=1:
                j = 1
                heights = h
                print('adding line 1 ', day)
                
            else:
                heights = heights + h
                mh[str(day)] = heights
                print('adding line 2', day)
                j = 0
'''