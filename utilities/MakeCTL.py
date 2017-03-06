#!/usr/bin/env python 3.5.2
# -*- coding: utf-8 -*-
"""
@author: slawler@dewberry.com
Created on Sun Feb 19 09:35:33 2017


###############################################################################

MakeControl.py                   USAGE NOTES


In User Inputs, make the following changes:
    (Assumes Repo cloned from github with directory structure in tact)
    
# Required Variables    
    1. Assign directory path & filenames 
    2. Input the datetime (month & days) for desired prediction 
    3. Assign Values for datum shift at each station, no value will result in 0
    
# Optional Changes   
    4. If year not 2017, update lines 9-13 from YEARDATA file
    5. If desired, change line 1

###############################################################################
"""
#------------Load Python Modules--------------------#
import pandas as pd
import numpy as np
import os 


#----Required Variables

#root_dir = r'C:\Users\sml\Desktop\tz' #--Root Directory
root_dir = '/Users/slawler/Desktop/tz'
poi_file = 'poi.in'                   # Stations File (poi.in)
h_constits = 'elev_hc.out'            # dbex output file (elev_hc.out)
dtm = '020128'                        # Month, start day, endday ('020128')
datum_dict = {'8594900': '001548'}    # Datum Dictionary for ctl file


#-----Optional Changes 

#CTL Header Data
line_1 = '1      000-0000 ----------------    T.M. 75 W.     Used 2001->           0000000'

# Year Data
line_9 = '2017 0 1103228621000   01032 241 906  5310642124 847283910981386 9352915   1   1'
line_10 = '2017 0 21000   010643103103211561000   01032213610321221 572260610322768   1   1'
line_11 = '2017 0 310001800 9661023 86726431109262010002017100028081032 738 6962583   1   1'
line_12 = '2017 0 4 8471133 847 2191000  2410001776 8471198100034921032 7381048 693   1   1'
line_13 = '2017 0 51165  20 9652071 78619131132 64810322862 718218710322970 8741917   1   1'

#----------------------------------RUN SCRIPT--------------------------------#
h_in_dir =  os.path.join(root_dir, 'harmonics')      
h_in_dir =  os.path.join(h_in_dir, 'outputs')      #--Input Directory

h_out_dir = os.path.join(root_dir, 'predictions')  #--Output Directory
h_out_dir = os.path.join(h_out_dir , 'ctl_files')  #--Output Directory

pois  = os.path.join(h_in_dir, poi_file)           #--Stations File Path
constits = os.path.join(h_in_dir,h_constits)       #--Constituents File Path


#---Function to convert meters to feet
def m2ft(val_m):
    val_ft = 3.28084*val_m
    return val_ft

#---Read in Station ID's & Order from poi.in
station_dict = dict()
with open(pois, 'r') as f:
    n_stations = f.readline().strip()
    n_stations = int(n_stations)
    for i in range(n_stations):
        line = f.readline().split()
        station = line[3]
        station_dict[i] = station

cols = list(station_dict.values())

#def MakeControl(s, dtm, constits, cols, datum_dict):  
idx=['M(2)','S(2)','N(2)','K(1)','M(4)','O(1)','M(6)','MK(3)',
     'S(4)','MN(4)','Nu(2)','S(6)','Mu(2)','2N(2)','OO(1)',
     'Lambda(2)','S(1)','M(1)','J(1)','Mm','Ssa','Sa','Msf',
     'Mf','Rho(1)','Q(1)','T(2)','R(2)','2Q(1)','P(1)','2SM(2)',
     'M(3)','L(2)','2MK(3)','K(2)','M(8)','MS(4)']
   
dfK = pd.DataFrame(data = 0,index =idx, columns = cols)
dfH = pd.DataFrame(data = 0,index =idx, columns = cols)

#---Read in Computed Harmonic Constituents 
with open(constits, 'r') as f:
    try:
        for i in range((n_stations + 1)*37):
            line = f.readline().split()
            constit = line[0]
            if constit in dfK.index:
                for j in range(n_stations):
                    line = f.readline().split()
                    H = m2ft(float(line[0]))*1000
                    K = float(line[1])*10
                    
                    H = int(H)
                    K = int(K)
                    
                    station = station_dict[j]
                    #print('\t',constit,station,  H,K)   
                    dfH[station].ix[constit] = H
                    dfK[station].ix[constit] = K
    except:
        print('end of file')


#---Loop throught columns (stations from poi) & Create ctl files for each        
for s in  dfH:       
    #---Format Harmonics for ctl file
    cs = []
    for c in idx:
        h = str(dfH[s].ix[c]) 
        if len(h) == 1:
            h =  '   ' + h 
        elif len(h) == 2:
            h = '  ' + h 
        elif len(h) == 3:
            h = ' ' + h
            
        k = str(dfK[s].ix[c]) 
        if len(k) == 1:
            k =  '   '+ h 
        elif len(k) == 2:
            k = '  ' + k 
        elif len(k) == 3:
            k = ' ' + k
            
        output = h+k
        cs.append(output) 
        print(c, output)        
    
    try:
        datum = datum_dict[s]
        line_2 = '{} 3 0                                                                      '.format(datum)
    
    except:
        line_2 = '{} 3 0                                                                      '.format('000000')
        
    line_3 = '2949   1' 
    for i in range(7):
        line_3 = line_3 + ' '+ cs[i]
    #print(line_3)
        
    line_4 = '2949   2' 
    for i in range(7,14):
        line_4 = line_4 + ' '+ cs[i]    
    #print(line_4)
        
    line_5 = '2949   3' 
    for i in range(14,21):
        line_5 = line_5 + ' '+ cs[i]   
    #print(line_5)
        
    line_6 = '2949   4' 
    for i in range(21,28):
        line_6 = line_6 + ' '+ cs[i]   
    #print(line_6)
        
    line_7 = '2949   5' 
    for i in range(28,35):
        line_7 = line_7 + ' '+ cs[i]   
    #print(line_7)
        
    line_8 = '2949   6' 
    for i in range(35,37):
        line_8 = line_8 + ' '+ cs[i]      
    #print(line_8)    
    
    line_14 = dtm
    line_15 = '   0   0   0'
       
    output_file = os.path.join(h_out_dir,'{}.ctl'.format(s))   
    with open(output_file, 'w') as f:
        f.write(line_1 + '\n')
        f.write(line_2 + '\n')
        f.write(line_3 + '\n')
        f.write(line_4 + '\n')
        f.write(line_5 + '\n')
        f.write(line_6 + '\n')
        f.write(line_7 + '\n')
        f.write(line_8 + '\n')
        f.write(line_9 + '\n')
        f.write(line_10 + '\n')
        f.write(line_11 + '\n')
        f.write(line_12 + '\n')
        f.write(line_13 + '\n')
        f.write(line_14 + '\n')
        f.write(line_15 + '\n')
