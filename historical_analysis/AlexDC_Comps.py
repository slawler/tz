# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 14:54:18 2017

@author: admin
"""

#--Import Python Libraries
import pandas as pd
import numpy as np
import requests
import json
from datetime import datetime
from collections import OrderedDict
import matplotlib.pyplot as plt
gage = '8594900'

# FEY
s0 = {1: '2008-08-25', 2:'2008-09-01'}

# IRENE
s1 = {1: '2011-08-25', 2:'2011-09-01'}

# LEE
s2 = {1: '2011-09-06', 2:'2011-09-15'}

# SANDY
s3 = {1: '2012-10-08', 2:'2012-10-15'}

# HERMINE
s4 = {1: '2016-09-01', 2:'2016-09-10'}

         
storm = s4

t0, t1 = storm[1], storm[2]
alex = AlexTidal(t0, t1)
gtown = GetTideObservation(gage, noaa_time(t0), noaa_time(t1))
anac = AnacTidal(t0, t1)
plt.plot(alex)
plt.plot(gtown)
plt.plot(anac)

plt.legend(labels = ['Alexandria','Georgetown','Anacostia'])
plt.grid()

plt.title('Hermine')


