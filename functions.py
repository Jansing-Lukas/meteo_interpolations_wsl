# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 13:28:31 2019

@author: lukas jansing

functions file

"""
import datetime
import math

#--------------------------------------------------------- 
# Function to create a datelist for given start, end and hstep
#---------------------------------------------------------
def make_datelist(startdate,enddate,hstep):
    nowdate = startdate
    datelist = []
    while nowdate <= enddate:
        datelist.append(nowdate)
        nowdate = nowdate + datetime.timedelta(hours=hstep)

    return datelist

#--------------------------------------------------------- 
# Function to calculate distance between two points
#---------------------------------------------------------
def distance(origin, destination):
    """
    Calculate the Haversine distance.

    Parameters
    ----------
    origin : tuple of float
        (lat, long)
    destination : tuple of float
        (lat, long)

    Returns
    -------
    distance_in_km : float

    Examples
    --------
    >>> origin = (48.1372, 11.5756)  # Munich
    >>> destination = (52.5186, 13.4083)  # Berlin
    >>> round(distance(origin, destination), 1)
    504.2
    
    credits for this function to: https://stackoverflow.com/questions/19412462/
    getting-distance-between-two-points-based-on-latitude-longitude/38187562#38187562
    
    """
    lat1, lon1 = origin
    lat2, lon2 = destination
    radius = 6371  # km

    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) * math.sin(dlat / 2) +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
         math.sin(dlon / 2) * math.sin(dlon / 2))
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    d = radius * c

    return d