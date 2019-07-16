# -*- coding: utf-8 -*-
"""
Created on Wed May 22 11:32:57 2019

@author: lukas jansing

script to import precipitation data
data with 10min resolution from meteoswiss
time span: 01.01.11 00 UTC - 21.05.19 23:50 UTC

afterwards use precip data to interpolate
interpolation in a manual way
save data for further use

"""
#--------------------------------------------------------- 
# Import modules
#---------------------------------------------------------
import os
os.environ['PROJ_LIB'] = 'set environmental path variable to directory proj4'

import numpy as np
import pandas
import datetime

#--------------------------------------------------------- 
# Import functions
#---------------------------------------------------------
from functions import make_datelist,distance
from import_data import import_precipitation_data

#--------------------------------------------------------- 
# Define paths
#---------------------------------------------------------
precippath = 'add path to MeteoSwiss 10min precipitation data'
metapath   = 'add path to Metadata'
figpath    = 'add path to your figures folder'

#--------------------------------------------------------- 
# Define stationnames and coordinates
# Note: Not all stations of meteoswiss are considered!!!
# Only the ones identified as three nearest stations already
# For the full procedure: See interpolate_precipitation_hourly.py
#---------------------------------------------------------
stationnames = ['sorniotlacinferieur', 'derborence', 'iserables', 'visperterminen',
       'ergisch', 'sierre', 'jeizinen', 'baltschiedertal', 'vercorin',
       'anzere', 'ilanz', 'thusis', 'savognin', 'weissfluhjoch', 'davos',
       'latsch', 'tschiertschen', 'valbella', 'chur', 'klostersaeuja',
       'vättis', 'lohnsh', 'schaffhausen', 'hallau', 'möhlin',
       'rünenberg', 'mervelier', 'delemont', 'baselbinningen',
       'ebnatkappel', 'hörnli', 'zürichaffoltern', 'opfikon', 'weesen',
       'doggenbenken', 'jona', 'wädenswil', 'zürichfluntern', 'sihlbrugg',
       'luzern', 'cham', 'zwillikon', 'stetten', 'interlaken', 'frutigen',
       'villarstiercelin', 'nesselboden', 'wynau', 'huttwil', 'egolzwil',
       'langenbruck', 'gösgen', 'attelwil', 'mosen', 'muriag',
       'unterbözberg', 'psiwürenlingen', 'beznau', 'oberehrendingen',
       'visp', 'leukerbad', 'montana', 'sion', 'martignyravoire', 'pully',
       'lausanne', 'ladole', 'nyonchangins', 'genevecointrin',
       'magadinocadenazzo', 'montegeneroso', 'coldrerio', 'lugano',
       'cranatorricella', 'stabio', 'buffalora', 'susch', 'scuol',
       'martina', 'kiental', 'dietikon', 'rothenbrunnen']

stationcoords_lat = np.array((46.166594 , 46.2885   , 46.1609252, 46.262883 , 46.293008 ,
                              46.298775 , 46.328389 , 46.342317 , 46.237964 , 46.305244 ,
                              46.775042 , 46.709414 , 46.594758 , 46.833328 , 46.812972 ,
                              46.627278 , 46.819667 , 46.755042 , 46.870369 , 46.865153 ,
                              46.909081 , 47.751947 , 47.689844 , 47.697278 , 47.572194 ,
                              47.434569 , 47.346617 , 47.351703 , 47.541136 , 47.273389 ,
                              47.370864 , 47.427694 , 47.437636 , 47.131622 , 47.188736 ,
                              47.223908 , 47.220961 , 47.377925 , 47.220564 , 47.036442 ,
                              47.188278 , 47.290094 , 47.396503 , 46.672236 , 46.599006 ,
                              46.621775 , 47.245167 , 47.255022 , 47.114    , 47.179428 ,
                              47.352883 , 47.363144 , 47.265233 , 47.243842 , 47.265381 ,
                              47.481697 , 47.536475 , 47.557256 , 47.507639 , 46.3029   ,
                              46.367036 , 46.298806 , 46.218647 , 46.109786 , 46.512278 ,
                              46.527817 , 46.424789 , 46.401047 , 46.247517 , 46.160031 ,
                              45.927606 , 45.853928 , 46.004228 , 46.075819 , 45.843458 ,
                              46.648139 , 46.754883 , 46.793275 , 46.884011 , 46.582578 ,
                              47.413558 , 46.765708))

stationcoords_lon = np.array((7.099811,  7.232511,  7.24875 ,  7.904211,  7.716297,  7.556392,
                              7.722478,  7.879328,  7.556286,  7.407575,  9.215344,  9.441703,
                              9.597131,  9.806383,  9.843547,  9.753694,  9.608131,  9.554422,
                              9.5305  ,  9.889689,  9.443592,  8.677622,  8.620125,  8.470458,
                              7.8779  ,  7.879406,  7.490036,  7.349561,  7.583519,  9.108486,
                              8.941633,  8.517942,  8.560325,  9.086081,  8.9612  ,  8.848278,
                              8.677694,  8.565731,  8.578272,  8.301014,  8.464631,  8.431947,
                              8.299978,  7.870192,  7.657542,  6.710075,  7.508575,  7.787469,
                              7.837358,  8.00475 ,  7.759586,  7.973725,  8.050511,  8.232819,
                              8.327597,  8.131836,  8.226931,  8.233311,  8.339211,  7.842964,
                              7.621644,  7.460822,  7.330211,  7.062525,  6.667522,  6.643272,
                              6.099461,  6.227728,  6.12775 ,  8.933667,  9.017872,  8.997431,
                              8.960319,  8.895361,  8.932364, 10.267183, 10.083   , 10.283253,
                              10.463017,  7.727097,  8.404461,  9.418456))

#--------------------------------------------------------- 
# Import coordinates of WSL network (FE Boden Netzwerk and additional TreeNet stations)
#---------------------------------------------------------
# stationcoords_2 contains a combined coordinate list of FE Boden network and TreeNet
treenetmeta = pandas.read_csv(metapath+'\stationcoords_2.csv',sep=';')
treenetmeta = np.array(treenetmeta)
        
treenetcoords_lat = treenetmeta[:,0]
treenetcoords_lon = treenetmeta[:,1]

treenetstationnames = ['Saillon_Versuch','Sihlwald','Bärschwil_tief','Neunkirch_Nord'
                       ,'Versam_Bu','Felsberg_Bu','Chamoson_Bu','Vetroz_Bu','Saillon_Bu'
                       ,'LWF-Neunkirch1','Neunkirch_SW_Bu','Bärschwil_flach','Tamins'
                       ,'Neunkirch_SW_Ei','Chamoson_Ei_Bu','Remigen','Bueren','Chamoson_Ei'
                       ,'Saillon_Ei','Chippis_Ei','Vetroz_Ei','Saillon_extrem','Chippis_top'
                       ,'Surava_8300','Surava_8306','Tarasp_9107_Fi','Ransun_Fi','Sent_9152_Fi'
                       ,'Bhutan','Pulligen','Hohtenn','Geissberg','Sent_9152_Foe','Scuol_9107_Foe'
                       ,'Ransun_Foe','Felsberg_Foe','Versam_Foe','Surava_8106','Alvaneu_8101'
                       ,'Alvaneu_8115','Alvaneu_8134','LWF-Lens3','Pfynwald','Chippis_Foe'
                       ,'Bachtel','Beatenberg','Birmensdorf','Davos','Grosswangen','Jussy'
                       ,'Laegeren_FF','Laegeren_Hut','Lausanne','Muri_Beech','Muri_Spruce'
                       ,'Muri_Meteo','Neunkirch_SE','Neunkirch_N','Neunkirch_SW','Novaggio'
                       ,'Pfynwald-Illgraben_NW','Pfynwald-Illgraben_N','Riehen_Forest'
                       ,'Riehen_Meteo','Sagno_SW','Sagno_SE','Sagno_Meteo','Saillon_1'
                       ,'Saillon_2','Saillon_3','Salgesch','Schaenis','Schmitten'
                       ,'Sempach','Surava_S','Surava_N','Visp','Vordemwald','Zürich','Wangen_SW']

#--------------------------------------------------------- 
# Import data from all meteoswiss stations
#---------------------------------------------------------
data = import_precipitation_data(stationnames=stationnames,resolution='10min',path=precippath)
   
#--------------------------------------------------------- 
# Create datelist for time period considered
#---------------------------------------------------------
station = stationnames[0]   # in order to have non-nan list
startdate = datetime.datetime(int(data[station][1,0]),int(data[station][2,0]),int(data[station][3,0]),int(data[station][4,0]),int(data[station][5,0]))   
enddate   = datetime.datetime(int(data[station][1,-1]),int(data[station][2,-1]),int(data[station][3,-1]),int(data[station][4,-1]),int(data[station][5,-1]))    
hstep     = 1/6

datelist  = make_datelist(startdate,enddate,hstep)

#for i in range(0,len(stationnames)):
#    print(data[stationnames[i]].shape)

#---------------------------------------------------------
# shift precip data into separate container
#---------------------------------------------------------
precip = np.zeros((len(stationnames),len(datelist)))
for i in range(0,len(stationnames)):
    station     = stationnames[i]
    precip[i,:] = data[station][6,:]
    
    # Filtering: Convert unrealistic high values into nans
    for j in range(0,len(datelist)):
        if precip[i,j] > 50:
            precip[i,j] = np.nan
            
#---------------------------------------------------------
# For each WSL site, identify the tree closest Meteoswiss stations
#---------------------------------------------------------
# Loop over all WSL sites
nearest          = []
distances        = []
closest_stations = []
for i in range(0,len(treenetcoords_lat)):
    
    # Calculate distance to all meteoswiss stations for each WSL site
    distance_to_meteo = []
    for j in range(0,len(stationnames)):
        distance_to_meteo.append(distance((treenetcoords_lat[i],treenetcoords_lon[i]),(stationcoords_lat[j],stationcoords_lon[j])))
    
    # Find indices of the three nearest stations
    # method using np.argpartition    
    distance_to_meteo = np.array(distance_to_meteo)
    k = 3
    idx = np.argpartition(distance_to_meteo,k)
    nearest.append(idx[:k])
    
    # Save distance of three nearest stations
    distances.append(distance_to_meteo[idx[:k]])
    
    # Print out result for manual check
    print('The 3 nearest stations for '+treenetstationnames[i]+' are '
          +stationnames[idx[0]]+', '+stationnames[idx[1]]+' and '+stationnames[idx[2]])
    
    # Save the nearest stations to a list
    closest_stations.append([stationnames[idx[0]],stationnames[idx[1]],stationnames[idx[2]]])
    
nearest = np.array(nearest)
distances = np.array(distances)
    
#--------------------------------------------------------------------------
# Perform horizontal interpolation for all WSL sites using the three nearest Meteoswiss sites
# Weighting according to the distance to the stations
# Create exceptions for NaN gaps in data
#--------------------------------------------------------------------------

#---------------------------------------------------------
# Loop over all WSL sites to create weights
# Option 2 should preferentially be chosen
#---------------------------------------------------------
weight = []
#weight_calculation = 'option1'
weight_calculation = 'option2'
for i in range(0,len(treenetcoords_lat)):
    if weight_calculation == 'option1':
        # Calculate weighting factor according to distance (normalised)
        # formula: weight = (total_distance - distance)/sum(total_distance - distance)
        # weights_summed = 2*sum(distances[i,:])
        weight_summed = (sum(distances[i,:]) - distances[i,0]) + (sum(distances[i,:]) - distances[i,1]) + (sum(distances[i,:]) - distances[i,2])
        weight.append([(sum(distances[i,:]) - distances[i,0])/weight_summed,(sum(distances[i,:]) - distances[i,1])/weight_summed,(sum(distances[i,:]) - distances[i,2])/weight_summed])
    if weight_calculation == 'option2':
        # formula: weight w1 = (b*c)/(b*c+a*c+a*b)
        weight_summed = distances[i,1]*distances[i,2]+distances[i,0]*distances[i,2]+distances[i,0]*distances[i,1]
        weight.append([(distances[i,1]*distances[i,2])/weight_summed,(distances[i,0]*distances[i,2])/weight_summed,(distances[i,0]*distances[i,1])/weight_summed])
weight = np.array(weight)

#---------------------------------------------------------
# Loop over all stations and all timesteps to interpolate precip data
#---------------------------------------------------------
precip_interpolated = np.zeros((len(treenetcoords_lat),len(datelist)))
for i in range(0,len(treenetcoords_lat)):
    print('interpolation of '+treenetstationnames[i])
    for j in range(0,len(datelist)):
        
        #---------------------------------------------------------
        # Separate treatment for cases where stations get excluded
        #---------------------------------------------------------
        if treenetstationnames[i] in ['Bärschwil_tief','Bärschwil_flach','Bueren','Sent_9152_Fi','Sent_9152_Foe','Grosswangen','Jussy',
                                      'Novaggio','Riehen_Forest','Riehen_Meteo','Salgesch','Schaenis']:
            
           # Find index of station with max distance to exclude it --> set weight to zero
           exclude_index = np.argmax(distances[i,:])
           weight[i,exclude_index] = 0
           # Find indices of remaining two stations
           include_index_1 = np.where(weight[i,:] != 0)[0][0]
           include_index_2 = np.where(weight[i,:] != 0)[0][1]
           
           # Option if none of the two station has NaN values
           if (np.isnan(precip[nearest[i,include_index_1],j]) == False and np.isnan(precip[nearest[i,include_index_2],j]) == False):
               precip_interpolated[i,j] = precip[nearest[i,include_index_1],j]*(distances[i,include_index_2]/(distances[i,include_index_1]+distances[i,include_index_2])) \
               + precip[nearest[i,include_index_2],j]*(distances[i,include_index_1]/(distances[i,include_index_1]+distances[i,include_index_2]))
               
           # Options if one of the stations has NaN values
           elif (np.isnan(precip[nearest[i,include_index_1],j]) == False):
               precip_interpolated[i,j] = precip[nearest[i,include_index_1],j]
           elif (np.isnan(precip[nearest[i,include_index_2],j]) == False):
               precip_interpolated[i,j] = precip[nearest[i,include_index_2],j]
               
           # Option if no station actually has values
           else:
               precip_interpolated[i,j] = np.nan
        #---------------------------------------------------------   
        # Separate treatment for cases where one Meteoswiss station fits best
        #---------------------------------------------------------
        elif treenetstationnames[i] in ['Felsberg_Bu','Chippis_Ei','Chippis_top','Tarasp_9107_Fi','Scuol_9107_Foe','Felsberg_Foe',
                                'Chippis_Foe','Davos','Muri_Beech','Muri_Spruce','Muri_Meteo','Sagno_SW','Sagno_SE','Sagno_Meteo','Visp']:
            
            # Find index of closest station
            closest_index = np.argmin(distances[i,:])
        
            # Option if it is not NaN
            if (np.isnan(precip[nearest[i,closest_index],j]) == False):
                precip_interpolated[i,j] = precip[nearest[i,closest_index],j]
                
            # Option if it is NaN
            else:
                precip_interpolated[i,j] = np.nan
        
        #---------------------------------------------------------
        # Separate treatment for cases where stations have to be switched
        #---------------------------------------------------------

        #---------------------------------------------------------
        # Normal treatment in case all stations are used
        #---------------------------------------------------------
        else:
            # Option if no station has NaN values
            if (np.isnan(precip[nearest[i,0],j]) == False and np.isnan(precip[nearest[i,1],j]) == False and np.isnan(precip[nearest[i,2],j]) == False): 
                precip_interpolated[i,j] = precip[nearest[i,0],j]*weight[i,0] + precip[nearest[i,1],j]*weight[i,1] + precip[nearest[i,2],j]*weight[i,2]

            # In case of one NaN station, but others not --> change weights: (sum(new_distances)-distance)/sum(distances)
        
            # Option if station 1 has NaN values but 2 and 3 not
            elif (np.isnan(precip[nearest[i,0],j]) == True and np.isnan(precip[nearest[i,1],j]) == False and np.isnan(precip[nearest[i,2],j]) == False):
                precip_interpolated[i,j] = precip[nearest[i,1],j]*(distances[i,2]/(distances[i,1]+distances[i,2])) + precip[nearest[i,2],j]*(distances[i,1]/(distances[i,1]+distances[i,2]))
            
            # Option if station 2 has NaN values but 1 and 3 not
            elif (np.isnan(precip[nearest[i,1],j]) == True and np.isnan(precip[nearest[i,0],j]) == False and np.isnan(precip[nearest[i,2],j]) == False):
                precip_interpolated[i,j] = precip[nearest[i,0],j]*(distances[i,2]/(distances[i,0]+distances[i,2])) + precip[nearest[i,2],j]*(distances[i,0]/(distances[i,0]+distances[i,2]))
        
            # Option if station 3 has NaN values but 1 and 2 not
            elif (np.isnan(precip[nearest[i,2],j]) == True and np.isnan(precip[nearest[i,1],j]) == False and np.isnan(precip[nearest[i,0],j]) == False):
                precip_interpolated[i,j] = precip[nearest[i,1],j]*(distances[i,0]/(distances[i,1]+distances[i,0])) + precip[nearest[i,0],j]*(distances[i,1]/(distances[i,1]+distances[i,0]))
        
            # In case of only one station with values --> just add those (nearest principle)
        
            elif (np.isnan(precip[nearest[i,0],j]) == False):
                precip_interpolated[i,j] = precip[nearest[i,0],j]
            
            elif (np.isnan(precip[nearest[i,1],j]) == False):
                precip_interpolated[i,j] = precip[nearest[i,1],j]
            
            elif (np.isnan(precip[nearest[i,2],j]) == False):
                precip_interpolated[i,j] = precip[nearest[i,2],j]

            # Option if all stations have NaN values
            else:
                precip_interpolated[i,j] = np.nan

#---------------------------------------------------------
# Quality check of interpolation
#---------------------------------------------------------
# Check for all-NaN timeseries
# Not the case
np.nanmax(precip_interpolated,axis=1)
    
# Check for amount of NaN values in timeseries
for i in range(0,len(treenetcoords_lat)):
    print(np.where(precip_interpolated[i,:] >= 0)[0].shape[0]/len(datelist))

#--------------------------------------------------------- 
# Save data for further use
#---------------------------------------------------------
interpolated_precipitation_data = np.concatenate((np.array([datelist]),precip_interpolated),axis=0)
np.save(precippath+'\interpolated_precipitation_10minres_newmethod.npy',interpolated_precipitation_data)

