# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 11:02:38 2019

@author: lukas jansing

script to import relative humidity data
data with 10min resolution from meteoswiss
time span: 01.01.11 00 UTC - 21.05.19 23:50 UTC

afterwards use relative humidity data to interpolate
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
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from scipy import stats

#--------------------------------------------------------- 
# Import functions
#---------------------------------------------------------
from functions import make_datelist,distance
from import_data import import_meteoswiss_data

#--------------------------------------------------------- 
# Define paths
#---------------------------------------------------------
relhumpath = 'add path to MeteoSwiss relative humidity data'
metapath   = 'add path to Metadata'
figpath    = 'add path to your figures folder'

#--------------------------------------------------------- 
# Define stationnames and coordinates
#---------------------------------------------------------
meteoswissmeta = pandas.read_csv(metapath+'\\relhum_stations.csv',sep=';')
meteoswissmeta = np.array(meteoswissmeta)
stationcoords_lat = meteoswissmeta[:,0]
stationcoords_lon = meteoswissmeta[:,1]
station_heights   = meteoswissmeta[:,2]

stationnames = ['aadorftänikon','acquarossacomprovasco','adelboden','aigle',
                'altdorf','altenrhein','andeer','andermatt','arosa','badragaz',
                'baselbinningen','bergünlatsch','bernzollikofen','beznau','biasca',
                'biere','binn','bischofszell','bivio','blattenlötschental',
                'boltigen','bouveret','buchsaarau','buffalora','bulletlafretaz',
                'cevio','cham1','cham2','chasseral','chateaudoex','chaumont',
                'chur','cimetta','coldesmosses','coldugrandstbernard','courtelary',
                'crapmasegn','cressier','davos','delemont','disentis',
                'ebnatkappel','eggishorn1','eggishorn2','egolzwil','einsiedeln',
                'elm','engelberg','evionnaz1','evionnaz2','evolenevilla','fahy',
                'flühlilu','fribourgposieux','frutigen','genevecointrin','gersau',
                'giswil','glarus','gornergrat','göschenen','gösgen','grächen',
                'grenchen','grimselhospiz','grono','gütschandermatt','güttingen',
                'hallau','hörnli','ilanz','interlaken','jungfraujoch',
                'koppigen','labrevine','lachauxdefonds','lachengalgenen','ladole',
                'lägern','langnauie','leibstadt','lemoleson','lesattelas1',
                'lesattelas2','lescharbonnieres','lesdiablerets1','lesdiablerets2',
                'lesmarecottes','locarnomonti','lugano','luzern','magadinocadenazzo',
                'marsens','mathod','matro','meiringen','möhlin','montagnierbagnes',
                'montana','montegeneroso','monterosaplattje','mosen','mottec',
                'mühleberg','nalunsschlivera','napf','neuchatel','nyonchangins',
                'oberrietkriessern','oron','passodelbernina','payerne','pilatus',
                'piotta','pizcorvatsch','pizmartegnas1','pizmartegnas2','plaffeien',
                'poschiavorobbia','pully','robiei','rünenberg','salenreutenen',
                'samedan','sanbernadino','säntis','sattelsz','schaffhausen','schiers',
                'schüpfheim','scuol','seglmaria','simplondorf','sion','stabio',
                'stamariavalmustair','stgallen','thun','titlis1','titlis2','ulrichen',
                'vaduz','valbella','vals','veveycorseaux','vicosoprano','villarstiercelin',
                'visp','wädenswil','würenlingenpsi','zürichaffoltern','zürichfluntern',
                'zürichkloten']

#--------------------------------------------------------- 
# Import coordinates of WSL network (FE Boden Netzwerk and additional TreeNet stations)
#---------------------------------------------------------
# stationcoords_2 contains a combined coordinate list of FE Boden network and TreeNet
treenetmeta = pandas.read_csv(metapath+'\stationcoords_2.csv',sep=';')
treenetmeta = np.array(treenetmeta)
        
treenetcoords_lat = treenetmeta[:,0]
treenetcoords_lon = treenetmeta[:,1]
treenet_heights   = treenetmeta[:,2]

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
# Create control datelist for time period considered
#---------------------------------------------------------
startdate = datetime.datetime(2011,1,1,0)   
enddate   = datetime.datetime(2019,5,21,23,50)
hstep     = 1/6

datelist_control  = make_datelist(startdate,enddate,hstep)

#--------------------------------------------------------- 
# Choose options:
# Path 1: Load and process station data
# Path 2: Load th processed data
#---------------------------------------------------------
data_new,stationnames_new = import_meteoswiss_data(stationnames=stationnames,variable='relhum',path=relhumpath,process='no',\
                                                   datelist_control=datelist_control)

# Check for shape and amount of values in the arrays
for i in range(0,len(stationnames_new)):
        print(data_new[()][stationnames_new[i]].shape)
        print(np.where(~np.isnan(data_new[()][stationnames_new[i]][1,:].astype(float)))[0].shape)
        
#---------------------------------------------------------
# shift relhum data into separate container
#---------------------------------------------------------
relhum = np.zeros((len(stationnames_new),len(datelist_control)))

for i in range(0,len(stationnames_new)):
    relhum[i,:] = data_new[()][stationnames_new[i]][1,:].astype(float)
    
    # Filtering: Convert unrealistic high values into nans
    for j in range(0,len(datelist_control)):
        if relhum[i,j] > 110 or relhum[i,j] < 0:
            relhum[i,j] = np.nan
            print('corrected a value!')

# Free the memory...
del data_new
    
#--------------------------------------------------------- 
# Visualize stations on map
#---------------------------------------------------------
#plotstations = 'yes'
plotstations = 'no'

if plotstations == 'yes':
    
    # create figure
    plt.close('all')
    fig = plt.figure(figsize=(22,30))
    
    # define map projection and add basic geographical features
    #map = Basemap(width=0.4e6,height=0.3e6,projection='stere',lat_ts=46.8,lat_0=46.8,lon_0=8.21,resolution='h')
    map = Basemap(projection='merc',llcrnrlat=45.7,urcrnrlat=47.95,\
            llcrnrlon=5.7,urcrnrlon=10.7,lat_ts=46.8,resolution='h',epsg=3395)
    map.drawcountries(linewidth=1.5)
    map.drawcoastlines()
    #map.drawrivers(color='blue')
    map.arcgisimage(service='World_Shaded_Relief',xpixels=2000)
    
    # plot meteoswiss automatic measurement network
    x,y = map(stationcoords_lon,stationcoords_lat)
    map.plot(x,y,'vk',markersize=10)
    
    # plot TreeNet station network
    i,j = map(treenetcoords_lon,treenetcoords_lat)
    map.plot(i,j,'ro',markersize=10)
    
    # export figure
    saveas = '\stationmap_relhum.png'
    #plt.show()
    plt.savefig(figpath+saveas,bbox_inches='tight')

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
    for j in range(0,len(stationnames_new)):
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
          +stationnames_new[idx[0]]+', '+stationnames_new[idx[1]]+' and '+stationnames_new[idx[2]])
    
    # Save the nearest stations to a list
    closest_stations.append([stationnames_new[idx[0]],stationnames_new[idx[1]],stationnames_new[idx[2]]])
    
nearest = np.array(nearest)
distances = np.array(distances)
closest_stations = np.array(closest_stations)

#--------------------------------------------------------------------------
# Perform horizontal interpolation for all WSL sites using the three nearest Meteoswiss sites
# Weighting according to the distance to the stations
# Create exceptions for NaN gaps in data
# Choose options for the vertical interpolation
#--------------------------------------------------------------------------
#vertical_interpolation = 'yes'
vertical_interpolation = 'no'
#method = 'standard_gradient'
method = 'empirical_gradient'

#---------------------------------------------------------
# Loop over all WSL sites to create weights
#---------------------------------------------------------
weight = []
for i in range(0,len(treenetcoords_lat)):
        # formula: weight w1 = (b*c)/(b*c+a*c+a*b)
        weight_summed = distances[i,1]*distances[i,2]+distances[i,0]*distances[i,2]+distances[i,0]*distances[i,1]
        weight.append([(distances[i,1]*distances[i,2])/weight_summed,(distances[i,0]*distances[i,2])/weight_summed,(distances[i,0]*distances[i,1])/weight_summed])
weight = np.array(weight)

#---------------------------------------------------------
# Loop over all stations and all timesteps to interpolate temp data
#---------------------------------------------------------
relhum_interpolated = np.zeros((len(treenetcoords_lat),len(datelist_control)))
relhum_unscaled = np.copy(relhum) # array to rescale stationdata

#visp_index = treenetstationnames.index('Visp')
for i in range(0,len(treenetcoords_lat)):
#for i in range(visp_index,visp_index+1):
    print('interpolation of '+treenetstationnames[i])
    
    # Vertical interpolation
    if vertical_interpolation == 'yes':
        
        # Use standard gradient
        if method == 'standard_gradient':
        
            print('not implemented!')

        # Use empirical gradient
        if method == 'empirical_gradient':
            
            # Prepare data arrays
            x_values    = np.concatenate((np.full((len(datelist_control)),station_heights[nearest[i,0]]),\
                                          np.full((len(datelist_control)),station_heights[nearest[i,1]]),\
                                          np.full((len(datelist_control)),station_heights[nearest[i,2]])))
            temp_linreg = np.concatenate((relhum[nearest[i,0],:],relhum[nearest[i,1],:],relhum[nearest[i,2],:]))
            
            # Calculate mean ambient temperature gradient
            slope,intercept,r_value,p_value,std_err = stats.linregress(x_values[~np.isnan(temp_linreg)],\
                                                                       temp_linreg[~np.isnan(temp_linreg)])
            # Loop over 3 closest stations
            for j in range(0,3):
                
                # Calculate height difference for each station
                height_diff = treenet_heights[i] - station_heights[nearest[i,j]]
                
                # Apply gradient
                relhum[nearest[i,j],:] = relhum[nearest[i,j],:]+height_diff*slope
            
    # Horizontal interpolation
    for j in range(0,len(datelist_control)):
        
        #---------------------------------------------------------
        # Separate treatment for cases where stations get excluded
        #---------------------------------------------------------
        if treenetstationnames[i] in ['Neunkirch_Nord','LWF-Neunkirch1','Neunkirch_SW_Bu','Neunkirch_SW_Ei','Sent_9152_Fi',
                                      'Sent_9152_Foe','LWF-Lens3','Jussy','Muri_Beech','Muri_Spruce',
                                      'Muri_Meteo','Neunkirch_SE','Neunkirch_N','Neunkirch_SW',
                                      'Riehen_Forest','Riehen_Meteo']:
            
           # Find index of station with max distance to exclude it --> set weight to zero
           exclude_index = np.argmax(distances[i,:])
           weight[i,exclude_index] = 0
           # Find indices of remaining two stations
           include_index_1 = np.where(weight[i,:] != 0)[0][0]
           include_index_2 = np.where(weight[i,:] != 0)[0][1]
           
           # Option if none of the two station has NaN values
           if (np.isnan(relhum[nearest[i,include_index_1],j]) == False and np.isnan(relhum[nearest[i,include_index_2],j]) == False):
               relhum_interpolated[i,j] = relhum[nearest[i,include_index_1],j]*(distances[i,include_index_2]/(distances[i,include_index_1]+distances[i,include_index_2])) \
               + relhum[nearest[i,include_index_2],j]*(distances[i,include_index_1]/(distances[i,include_index_1]+distances[i,include_index_2]))
               
           # Options if one of the stations has NaN values
           elif (np.isnan(relhum[nearest[i,include_index_1],j]) == False):
               relhum_interpolated[i,j] = relhum[nearest[i,include_index_1],j]
           elif (np.isnan(relhum[nearest[i,include_index_2],j]) == False):
               relhum_interpolated[i,j] = relhum[nearest[i,include_index_2],j]
               
           # Option if no station actually has values
           else:
               relhum_interpolated[i,j] = np.nan
        #---------------------------------------------------------   
        # Separate treatment for cases where one Meteoswiss station fits best
        #---------------------------------------------------------
        elif treenetstationnames[i] in ['Felsberg_Bu','Tarasp_9107_Fi','Scuol_9107_Foe','Felsberg_Foe',
                                        'Davos','Laegeren_FF','Laegeren_Hut','Visp']:
            
            # Find index of closest station
            closest_index = np.argmin(distances[i,:])
        
            # Option if it is not NaN
            if (np.isnan(relhum[nearest[i,closest_index],j]) == False):
                relhum_interpolated[i,j] = relhum[nearest[i,closest_index],j]
                
            # Option if it is NaN
            else:
                relhum_interpolated[i,j] = np.nan
        
        #---------------------------------------------------------
        # Separate treatment for cases where stations have to be switched
        #---------------------------------------------------------

        #---------------------------------------------------------
        # Normal treatment in case all stations are used
        #---------------------------------------------------------
        else:
            # Option if no station has NaN values
            if (np.isnan(relhum[nearest[i,0],j]) == False and np.isnan(relhum[nearest[i,1],j]) == False and np.isnan(relhum[nearest[i,2],j]) == False): 
                relhum_interpolated[i,j] = relhum[nearest[i,0],j]*weight[i,0] + relhum[nearest[i,1],j]*weight[i,1] + relhum[nearest[i,2],j]*weight[i,2]

            # In case of one NaN station, but others not --> change weights: (sum(new_distances)-distance)/sum(distances)
        
            # Option if station 1 has NaN values but 2 and 3 not
            elif (np.isnan(relhum[nearest[i,0],j]) == True and np.isnan(relhum[nearest[i,1],j]) == False and np.isnan(relhum[nearest[i,2],j]) == False):
                relhum_interpolated[i,j] = relhum[nearest[i,1],j]*(distances[i,2]/(distances[i,1]+distances[i,2])) + relhum[nearest[i,2],j]*(distances[i,1]/(distances[i,1]+distances[i,2]))
            
            # Option if station 2 has NaN values but 1 and 3 not
            elif (np.isnan(relhum[nearest[i,1],j]) == True and np.isnan(relhum[nearest[i,0],j]) == False and np.isnan(relhum[nearest[i,2],j]) == False):
                relhum_interpolated[i,j] = relhum[nearest[i,0],j]*(distances[i,2]/(distances[i,0]+distances[i,2])) + relhum[nearest[i,2],j]*(distances[i,0]/(distances[i,0]+distances[i,2]))
        
            # Option if station 3 has NaN values but 1 and 2 not
            elif (np.isnan(relhum[nearest[i,2],j]) == True and np.isnan(relhum[nearest[i,1],j]) == False and np.isnan(relhum[nearest[i,0],j]) == False):
                relhum_interpolated[i,j] = relhum[nearest[i,1],j]*(distances[i,0]/(distances[i,1]+distances[i,0])) + relhum[nearest[i,0],j]*(distances[i,1]/(distances[i,1]+distances[i,0]))
        
            # In case of only one station with values --> just add those (nearest principle)
        
            elif (np.isnan(relhum[nearest[i,0],j]) == False):
                relhum_interpolated[i,j] = relhum[nearest[i,0],j]
            
            elif (np.isnan(relhum[nearest[i,1],j]) == False):
                relhum_interpolated[i,j] = relhum[nearest[i,1],j]
            
            elif (np.isnan(relhum[nearest[i,2],j]) == False):
                relhum_interpolated[i,j] = relhum[nearest[i,2],j]

            # Option if all stations have NaN values
            else:
                relhum_interpolated[i,j] = np.nan
                
    #---------------------------------------------------------
    # Rescale the data
    # Reason: They might possibly be used for other station
    #---------------------------------------------------------            
    if vertical_interpolation == 'yes':
        for j in range(0,3):
            relhum[nearest[i,j],:] = relhum_unscaled[nearest[i,j],:]

#---------------------------------------------------------
# Short quality check of interpolation
#---------------------------------------------------------
# Check for all-NaN timeseries
# Not the case
np.nanmax(relhum_interpolated,axis=1)
np.nanmin(relhum_interpolated,axis=1)
    
# Check for amount of values in timeseries
for i in range(0,len(treenetcoords_lat)):
    print(np.where(~np.isnan(relhum_interpolated[i,:]))[0].shape[0]/len(datelist_control))
    
# Compare two timeseries which should be identical: Visp --> are identical
visp_index_meteo = stationnames_new.index('visp')
visp_index_wsl   = treenetstationnames.index('Visp')

np.nanmax(relhum_interpolated[visp_index_wsl,:] - relhum[visp_index_meteo,:])

#--------------------------------------------------------- 
# Save data for further use
#---------------------------------------------------------
interpolated_relhum_data = np.concatenate((np.array([datelist_control]),relhum_interpolated),axis=0)
np.save(relhumpath+'\interpolated_relhum_10minres_novertical.npy',interpolated_relhum_data)