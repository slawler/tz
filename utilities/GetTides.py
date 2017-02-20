# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 16:16:29 2017

@author: slawler
"""

import pandas as pd
import requests
import json
from datetime import datetime
from collections import OrderedDict
import matplotlib.pyplot as plt

         
# Start date (year, month, day)
y0, m0 ,d0 = 2017, 2, 10         
y1, m1 ,d1 = 2017, 2, 27    
  
noaa_gage = '8594900'
usgs_gage = '01646500'

#USGS parameter
flow = "00060"  

# Create Datetime Objects
start     = datetime(y0, m0, d0,0)    
stop      = datetime(y1, m1 ,d1,0) 
  

def GetTidePrediction(gage, start, stop):   
    #--NOAA API https://tidesandcurrents.noaa.gov/api/
    datum     = "msl"   #"NAVD"                  #Datum
    units     = "metric"                         #Units
    time_zone = "gmt"                         #Time Zone
    fmt       = "json"                            #Format
    url       = 'http://tidesandcurrents.noaa.gov/api/datagetter'
    product   = 'predictions'                     #Product
    
    noaa_time_step = '6T'
    noaa = pd.DataFrame()
    gages = dict()
    
    t0     = start.strftime('%Y%m%d %H:%M')
    t1     = stop.strftime('%Y%m%d %H:%M')
    api_params = {'begin_date': t0, 'end_date': t1,
                'station': gage,'product':product,'datum':datum,
                'units':units,'time_zone':time_zone,'format':fmt,
                'application':'web_services' }
        
    pred=[];t=[]
    
    
    r = requests.get(url, params = api_params)
    jdata =r.json()
    
    for j in jdata['predictions']:
        t.append(str(j['t']))
        pred.append(str(j['v']))

        
    colname = str(gage)    
    noaa[colname]= pred
    noaa[colname] = noaa[colname].astype(float)
      
    idx = pd.date_range(start,periods = len(noaa.index), freq=noaa_time_step)   
    noaa = noaa.set_index(idx)  
    print('We have Predictions ' , gage) 
    return noaa

    

def GetTideObservation(gage, start, stop):   
    #--NOAA API https://tidesandcurrents.noaa.gov/api/
    datum     = "msl"   #"NAVD"                  #Datum
    units     = "metric"                         #Units
    time_zone = "gmt"                         #Time Zone
    fmt       = "json"                            #Format
    url       = 'http://tidesandcurrents.noaa.gov/api/datagetter'
    product   = 'water_level'                     #Product
    
    noaa_time_step = '6T'
    noaa = pd.DataFrame()
    gages = dict()
    
    t0     = start.strftime('%Y%m%d %H:%M')
    t1     = stop.strftime('%Y%m%d %H:%M')
    api_params = {'begin_date': t0, 'end_date': t1,
                'station': gage,'product':product,'datum':datum,
                'units':units,'time_zone':time_zone,'format':fmt,
                'application':'web_services' }
        
    pred=[];t=[]
    
    
    r = requests.get(url, params = api_params)
    jdata =r.json()
    
    for j in jdata['data']:
        t.append(str(j['t']))
        pred.append(str(j['v']))

        
    colname = str(gage)    
    noaa[colname]= pred
    noaa[colname] = noaa[colname].astype(float)
         
    idx = pd.date_range(start,periods = len(noaa.index), freq=noaa_time_step)   
    noaa = noaa.set_index(idx)  
    print('We have Observations ' , gage) 
    return noaa

    
def Get_USGS(gage, parameter, start, stop):  
               
    dformat    = "json"                                  # Data Format  
    url        = 'http://waterservices.usgs.gov/nwis/iv' # USGS API
    
    # Format Datetime Objects for USGS API
    first    =  datetime.date(start).strftime('%Y-%m-%d')
    last     =  datetime.date(stop).strftime('%Y-%m-%d') 
    

    # Ping the USGS API for data
    params = OrderedDict([('format',dformat),('sites',gage),('startDT',first), 
                ('endDT',last), ('parameterCD',parameter)])  
    
    r = requests.get(url, params = params) 
    print("\nRetrieved Data for USGS Gage: ", gage)
    data = r.content.decode()
    d = json.loads(data)
    mydict = dict(d['value']['timeSeries'][0])
        
    
    if params['parameterCD'] == '00060':
        obser = "StreamFlow"
    else:
        obser = "Stage"
        
    
    # Great, We can pull the station name, and assign to a variable for use later:
    SiteName = mydict['sourceInfo']['siteName']
    print('\n', SiteName)
    
    
    # After reveiwing the JSON Data structure, select only data we need: 
    tseries = d['value']['timeSeries'][0]['values'][0]['value'][:]
    
    # Create a Dataframe, format Datetime data,and assign numeric type to observations
    df = pd.DataFrame.from_dict(tseries)
    df.index = pd.to_datetime(df['dateTime'],format='%Y-%m-%d{}%H:%M:%S'.format('T'))

    df.value = pd.to_numeric(df.value)
    
    # Get Rid of unwanted data, rename observed data
    df = df.drop('dateTime', 1)
    df = df.rename(columns = {'value':obser})
    return df    
     
    
df = GetTidePrediction(noaa_gage, start, stop)        
#df1 = GetTideObservation(noaa_gage, start, stop)    
#df2 = Get_USGS(usgs_gage, flow, start, stop)    
#df2['StreamFlow'] = df2['StreamFlow'].astype(float) 

#--NOAA Data
x = df.index        #X-axis is time from predicted df
y1 = df[noaa_gage]  #Y-axis 1 = predicted
#y2 = df1[noaa_gage] #Y-axis 2 = observed
#y3 = df1 - df       #Y-axis 2 = difference

#--Merge Data
#merge = pd.merge(y3, df2, left_index=True, right_index=True)
#merge['coeff'] = merge['8594900']/merge['StreamFlow']
#merge = pd.merge(merge,df, left_index=True, right_index=True )
#merge.rename(columns={'8594900_x': 'observed', '8594900_y':'prediced'}, inplace=True)
#merge.drop(axis = 1,labels = ["qualifiers"], inplace=True )


plt.plot(x, y1)  
#plt.plot(x, y2)  
#plt.plot(x, y3)




    