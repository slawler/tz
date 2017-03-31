# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 11:56:16 2017

@author: slawler
"""

def AlexTidal(start, stop):
    '''
	Pulls data for USGS 0165258890 POTOMAC RIVER AT CAMERON ST DOCK AT ALEXANDRIA, VA
	'''
    url = 'https://nwis.waterdata.usgs.gov/usa/nwis/uv/?cb_62620=on&format=rdb&site_no=0165258890&period=&begin_date={}&end_date={}'.format(start, stop)
    #--Read in Daily data, convert flow to stage, Skiprows Value will vary by gage
    df = pd.read_csv(url, skiprows=26, sep='\t') 
    df.drop(0, axis=0, inplace=True)
    df.rename(columns = {'146633_62620':'Stage', '146633_62620_cd':'Qual'}, inplace=True)
    df['Stage'] = pd.to_numeric(df['Stage'])
    df['datetime'] = pd.to_datetime(df['datetime'], format= '%Y-%m-%d')
    df = df[['datetime', 'Stage']]
    df.set_index('datetime', inplace=True)
    return df

    
def AnacTidal(start, stop):
    '''
	Pulls data for USGS 0165258890 POTOMAC RIVER AT CAMERON ST DOCK AT ALEXANDRIA, VA
	'''
    url = 'https://nwis.waterdata.usgs.gov/dc/nwis/uv?cb_00065=on&format=rdb&site_no=01651750&period=&begin_date={}&end_date={}'.format(start, stop)
    #--Read in Daily data, convert flow to stage, Skiprows Value will vary by gage
    df = pd.read_csv(url, skiprows=26, sep='\t') 
    df.drop(0, axis=0, inplace=True)
    df.rename(columns = {'70044_00065':'Stage', '146633_62620_cd':'Qual'}, inplace=True)
    df['Stage'] = pd.to_numeric(df['Stage'])
    df['datetime'] = pd.to_datetime(df['datetime'], format= '%Y-%m-%d')
    df = df[['datetime', 'Stage']]
    df.set_index('datetime', inplace=True)
    return df    
    

def noaa_time(t):
    time  = datetime.strptime(t, '%Y-%m-%d')
    return time    
    
def GetTidePrediction(gage, start, stop):   
    #--NOAA API https://tidesandcurrents.noaa.gov/api/
    datum     = "msl"   #"NAVD"                  #Datum
    units     = "english"                         #Units
    time_zone = "lst_ldt"                         #Time Zone
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
    units     = "english"                         #Units
    time_zone = "lst_ldt"                         #Time Zone
    fmt       = "json"                            #Format
    url       = 'http://tidesandcurrents.noaa.gov/api/datagetter'
    product   = 'hourly_height'                     #Product
    
    noaa_time_step = '60T'
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
    
    
    
    
    
    
    
def plotmultiple(df1, df2, df3, df4, df5, df1_name, df2_name, df3_name, Tides_Only,order, start_plot, stop_plot):
    f, (ax_1, ax_2, ax_3) = plt.subplots(3, sharex=True, sharey=True)
    df = df1[['datetime', 'Stage']].copy()
    df.set_index('datetime', inplace=True)
    df = df.ix[start_plot:stop_plot]
    x1 = df.index
    y1 = df['Stage'].resample('6T').mean()
    interpolated = y1.interpolate(method='cubic', order=2)


    ax_1.set_ylabel('{}'.format(Tides_Only), color='black')     
    ax_1.plot(df4.index, df4[Tides_Only], 'r')
    ax_1.plot(df5.index, df5[Tides_Only], 'b')
    ax_1a = ax_1.twinx()
    ax_1a.plot(interpolated, 'black')
    ax_1a.set_ylabel('{}'.format(df1_name), color='black') 


    df = df2[['datetime', 'Stage']].copy()
    df.set_index('datetime', inplace=True)
    df = df.ix[start_plot:stop_plot]
    x1 = df.index
    y1 = df['Stage'].resample('6T').mean()
    interpolated = y1.interpolate(method='cubic', order=2)


    ax_2.set_ylabel('{}'.format(Tides_Only), color='black')     
    ax_2.plot(df4.index, df4[Tides_Only], 'r')
    ax_2.plot(df5.index, df5[Tides_Only], 'b')
    ax_2a = ax_2.twinx()
    ax_2a.plot(interpolated, 'black')
    ax_2a.set_ylabel('{}'.format(df2_name), color='black')  


    #--df3
    df = df2[['datetime', 'Stage']].copy()
    df.set_index('datetime', inplace=True)
    df = df.ix[start_plot:stop_plot]
    x1 = df.index
    y1 = df['Stage'].resample('6T').mean()
    interpolated = y1.interpolate(method='cubic', order=2)


    ax_3.set_ylabel('{}'.format(Tides_Only), color='black')    
    ax_3.plot(df4.index, df4[Tides_Only], 'r')
    ax_3.plot(df5.index, df5[Tides_Only], 'b')
    ax_3a = ax_3.twinx()
    ax_3a.plot(interpolated, 'black')
    ax_3a.set_ylabel('{}'.format(df3_name), color='black')  
    
    
    f.autofmt_xdate()
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)    