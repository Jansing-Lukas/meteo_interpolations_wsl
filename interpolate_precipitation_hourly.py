# -*- coding: utf-8 -*-
"""
Created on Wed May  8 09:01:42 2019

@author: jansing

script to import precipitation data
data with hourly resolution from meteoswiss
time span: 01.01.11 00 UTC - 31.12.18 23 UTC

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
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
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
# Import MeteoSwiss coordinates and elevations
#---------------------------------------------------------
meteoswissmeta = pandas.read_csv(metapath+'\precip_stations.csv',sep=';')
meteoswissmeta = np.array(meteoswissmeta)
stationcoords_lat = meteoswissmeta[:,0]
stationcoords_lon = meteoswissmeta[:,1]
station_heights   = meteoswissmeta[:,2]

stationnames = ['arolla','barragegrandedixence','bricola','findelen','stafel',
                'lafouly','emosson','turtmann','moiry','salanfe','clusanfe',
                'rossberg','durnand','trient','sorniotlacinferieur',
                'derborence','tsanfleuron','nendaz','iserables','bruchji',
                'visperterminen','ergisch','blinnen','lescollons','sierre',
                'choex','bourgstpierre','saasbalen','jeizinen','baltschiedertal',
                'mattsand','vercorin','anzere','saleina','champery','disentissedrun',
                'trun','vrin','vals','ilanz','andeer','thusis','bivio','savognin',
                'pizmartegnas','weissfluhjoch','davos','latsch','arosa','tschiertschen',
                'valbella','chur','klostersaeuja','stantönien','schiers','vättis',
                'badragaz','vaduz','schaanuntereau','malbun','oberrietsg','altstättensg',
                'altenrhein','stgallen','amriswil','güttingen','salenreutenen','eschenz',
                'lohnsh','schaffhausen','hallau','leibstadt','wittnau','möhlin',
                'rünenberg','mervelier','bellelay','delemont','baselbinningen',
                'starkenbach','ebnatkappel','flawil','säntis','urnäsch','appenzell',
                'bischofszell','aadorftänikon','hörnli','winterthurseen','zürichaffoltern',
                'opfikon','zürichkloten','braunwald','elm','glarus','weesen','doggenbenken',
                'innerthal','rempen','siebnen','lachengalgenen','jona','wartau','wädenswil',
                'zürichfluntern','oberiberg','einsiedeln','sihlbrugg','gütschobandermatt',
                'andermatt','göscheneralp','göschenen','bristen','altdorf','sattelaegeri',
                'gersau','titlis','engelberg','giswil','stöckalp','luzern','pilatus','flühlilu',
                'schüpfheim','entlebuch','cham','zwillikon','stetten','grimselhospiz',
                'guttannen','meiringen','brienz','lauterbrunnen','interlaken','kandersteg',
                'adelboden','frutigen','boltigen','thun','bernzollikofen','mühleberg',
                'gsteiggstaad','chateaudoex','lemoleson','marsens','romont','fribourgpoiseux',
                'plaffeien','oron','payerne','lescharbonnieres','mathod','villarstiercelin',
                'bulletlafretaz','lauberson','couvet','labrevine','neuchatel','chaumont',
                'cressier','chasseral','magglingen','courtelary','grenchen','nesselboden',
                'napf','langnauie','affolternimemmental','riedholzwallierhof','koppigen',
                'wynau','huttwil','egolzwil','langenbruck','gösgen','buchsaarau','attelwil',
                'unterkulm','mosen','muriag','unterbözberg','psiwürenlingen','beznau',
                'oberehrendingen','ulrichen','fieschertal','binn','brig','zermatt','grächen',
                'visp','blattenlötschental','leukerbad','montana','mottec','evolenevilla',
                'sion','fionnay','coldugrandstbernard','orsieres','martignyravoire',
                'lesmarecottes','evionnaz','bex','coldesmosses','aigle','lesavants',
                'veveycorseaux','pully','lausanne','cossonay','biere','longirod','ladole',
                'nyonchangins','genevecointrin','lachauxdefonds','saignelegier','fahy',
                'mormont','airolo','piotta','faido','acquarossacomprovasco','biasca',
                'sbernardino','grono','bellinzona','magadinocadenazzo','robiei','cevio',
                'boscogurin','mosogno','cimetta','locarnomonti','montegeneroso',
                'coldrerio','lugano','cranatorricella','stabio','passodelbernina',
                'poschiavorobbia','brusiopiazzo','vicosoprano','soglio','seglmaria',
                'pizcorvatsch','berninacurtinatsch','samedan','buffalora','susch','scuol',
                'martina','stamariavalmustair','simplondorf','montagnierbagnes',
                'lavalsainte','belp','kiesen','zweisimmen','kiental','dietikon',
                'rothenbrunnen','safienplatz','zervreila']

#--------------------------------------------------------- 
# Import coordinates of WSL network (FE Boden Netzwerk and additional TreeNet stations)
#---------------------------------------------------------
# stationcoords_2 contains a combined coordinate list of FE Boden network and TreeNet
treenetmeta = pandas.read_csv(metapath+'\stationcoords_2.csv',sep=';')
treenetmeta = np.array(treenetmeta)
        
treenetcoords_lat = treenetmeta[:,0]
treenetcoords_lon = treenetmeta[:,1]
treenet_heights   = treenetmeta[:,2]

# Attribute region flag to treenet sites
# 0 = alpine_north, 1 = inneralpine, 2 = ticino
treenet_regionflag = treenetmeta[:,3]

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
data = import_precipitation_data(stationnames=stationnames,resolution='hourly',path=precippath)
   
#--------------------------------------------------------- 
# Create datelist for time period considered
#---------------------------------------------------------
station   = stationnames[0]
startdate = datetime.datetime(int(data[station][1,0]),int(data[station][2,0]),int(data[station][3,0]),int(data[station][4,0]))   
enddate   = datetime.datetime(int(data[station][1,-1]),int(data[station][2,-1]),int(data[station][3,-1]),int(data[station][4,-1]))    
hstep     = 1

datelist  = make_datelist(startdate,enddate,hstep)

#---------------------------------------------------------
# shift precip data into separate container
#---------------------------------------------------------
precip = np.zeros((len(stationnames),len(datelist)))
for i in range(0,len(stationnames)):
    station     = stationnames[i]
    precip[i,:] = data[station][6,:]
    
    # Filtering: Convert unrealistic high values into nans
    # Note: Only for one value in Kiental data (2015,10,3)
    for j in range(0,len(datelist)):
        if precip[i,j] > 200:
            precip[i,j] = np.nan

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
    saveas = '\stationmap.png'
    #plt.show()
    plt.savefig(figpath+saveas,bbox_inches='tight')
    
#--------------------------------------------------------- 
# Visualize LWF locations on map
#---------------------------------------------------------
#plotlwf  = 'yes'
plotlwf  = 'no'

if plotlwf == 'yes':
    lwfstations = ['Beatenberg','Jussy','Lausanne','LWF-Neunkirch1',
                   'Novaggio','Schaenis','Visp','Vordemwald']
    lwflons = []
    lwflats = []
    for i in range(0,len(lwfstations)):
        lwf_ind = treenetstationnames.index(lwfstations[i])
        lwflons.append(treenetcoords_lon[lwf_ind])
        lwflats.append(treenetcoords_lat[lwf_ind])
        
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
    
    # plot TreeNet station network
    i,j = map(lwflons,lwflats)
    map.plot(i,j,'ro',markersize=10)
    
    # export figure
    saveas = '\lwf_stations.png'
    #plt.show()
    plt.savefig(figpath+saveas,bbox_inches='tight')
   
#---------------------------------------------------------
# For each WSL site, identify the tree closest Meteoswiss stations
#---------------------------------------------------------
# Loop over all WSL sites
nearest          = []
distances        = []
closest_stations = []
height_diff_abs  = []
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
    
    # Calculate height differences
    height_diff_abs.append([abs(station_heights[idx[0]] - treenet_heights[i]),\
                       abs(station_heights[idx[1]] - treenet_heights[i]),\
                       abs(station_heights[idx[2]] - treenet_heights[i])])
    
nearest = np.array(nearest)
distances = np.array(distances)
height_diff_abs = np.array(height_diff_abs)
    
#--------------------------------------------------------------------------
# Perform horizontal interpolation for all WSL sites using the three nearest Meteoswiss sites
# Weighting according to the distance to the stations
# Create exceptions for NaN gaps in data
#--------------------------------------------------------------------------
vertical_interpolation = 'yes'
#vertical_interpolation = 'no'
method = 'standard_gradient'
#method = 'empirical_gradient'

# Define regional precipitation gradients
# From Uttinger (1951)
# Convert from yearly corrections to hourly
# In mm/ m
# Lower elevation gradient in first position
# Higher elevation gradient in second position
# Separation: 1700 m altitude
alpinenorth = np.array((85/(100*24*365),57/(100*24*365)))
inneralpine = np.array((27/(100*24*365),99/(100*24*365)))
ticino      = np.array((24/(100*24*365),24/(100*24*365)))

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
precip_unscaled     = np.copy(precip) # array to rescale stationdata
elevation_exclusion = 'yes'

for i in range(0,len(treenetcoords_lat)):
    print('interpolation of '+treenetstationnames[i])
    
    # Vertical interpolation
    if vertical_interpolation == 'yes':
        
        # Use standard gradient
        if method == 'standard_gradient':
 
            # Identify which gradient to choose
            if treenet_regionflag[i] == 0:
                gradient = alpinenorth
            if treenet_regionflag[i] == 1:
                gradient = inneralpine
            if treenet_regionflag[i] == 2:
                gradient = ticino
            
            # Loop over 3 closest stations
            # Note: This coding options only work for WSL sites lower than 1700 m
            for j in range(0,3):
                
                # Case if Meteoswiss station is above 1700 m
                if station_heights[nearest[i,j]] > 1700:
                  
                    # Calculate height difference for each station
                    # Separated by the 1700 m mark
                    height_diff_top   = station_heights[nearest[i,j]] - 1700
                    height_diff_below = 1700 - treenet_heights[i]
                    
                    # Apply the gradient
                    # Important: Only for non-zero values!
                    precip[nearest[i,j],:][precip[nearest[i,j],:] > 0] = precip[nearest[i,j],:][precip[nearest[i,j],:] > 0] - \
                                                                         height_diff_top*gradient[1] - height_diff_below*gradient[0]
                                             
                # Case if Meteoswiss station is below 1700 m    
                else:
                    
                    # Calculate height difference for each station
                    height_diff = station_heights[nearest[i,j]] - treenet_heights[i]
                    
                    # Apply gradient
                    precip[nearest[i,j],:][precip[nearest[i,j],:] > 0] = precip[nearest[i,j],:][precip[nearest[i,j],:] > 0]-height_diff*gradient[0]
                
                # Correct for negative precipitation values which can occur
                precip[nearest[i,j],:][precip[nearest[i,j],:] < 0] = 0

        # Use empirical gradient
        if method == 'empirical_gradient':
            
            print('not implemented!')
                
    # Horizontal interpolation
    for j in range(0,len(datelist)):
        
        #---------------------------------------------------------
        # Separate treatment for cases where stations get excluded
        #---------------------------------------------------------
        if treenetstationnames[i] in ['Bärschwil_tief','Bärschwil_flach','Bueren','Sent_9152_Fi','Sent_9152_Foe','Grosswangen','Jussy',
                                      'Novaggio','Riehen_Forest','Riehen_Meteo','Salgesch']:
            
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
        # Separate treatment for cases where elevation difference is too large
        #---------------------------------------------------------
        elif elevation_exclusion == 'yes' and \
             treenetstationnames[i] in ['Saillon_Versuch','Saillon_Bu','Saillon_Ei','Saillon_extrem','Bhutan','Pulligen','Pfynwald',
                                        'Bachtel','Beatenberg','Pfynwald-Illgraben_NW','Pfynwald-Illgraben_N','Saillon_1','Saillon_2',
                                        'Saillon_3']:
                 
             # Find index of station with largest elevation difference --> set weight to zero
             exclude_index = np.argmax(height_diff_abs[i,:])
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
    # Rescale the data
    #---------------------------------------------------------            
    if vertical_interpolation == 'yes':
        for j in range(0,3):
            precip[nearest[i,j],:] = precip_unscaled[nearest[i,j],:]

#---------------------------------------------------------
# Quality check of interpolation
#---------------------------------------------------------
# Check for all-NaN timeseries
# Not the case
np.nanmax(precip_interpolated,axis=1)
    
# Check for amount of NaN values in timeseries
# Huge improvement when including NaN handling!!!
data_availability = []
for i in range(0,len(treenetcoords_lat)):
    print(np.where(~np.isnan(precip_interpolated[i,:]))[0].shape[0]/len(datelist))
    data_availability.append(np.where(~np.isnan(precip_interpolated[i,:]))[0].shape[0]/len(datelist))
data_availability = np.array(data_availability)
print(np.mean(data_availability))
print(np.min(data_availability))
    
# Compare two timeseries which should be identical: Visp --> Are identical
visp_index_meteo = stationnames.index('visp')
visp_index_wsl   = treenetstationnames.index('Visp')

np.nanmax(precip_interpolated[visp_index_wsl,:] - precip[visp_index_meteo,:])

#--------------------------------------------------------- 
# Save data for further use
#---------------------------------------------------------
# new method: without exclusion of stations due to elevation diffs and without vertical interpolation
# stanard gradient: without exclusion of stations but with vertical interpolation
# version 3: with exclusion of stations but without vertical interpolation
# version 4: with exclusion of stations and with vertical interpolation
interpolated_precipitation_data = np.concatenate((np.array([datelist]),precip_interpolated),axis=0)
if elevation_exclusion == 'no' and vertical_interpolation == 'no':
    np.save(precippath+'\interpolated_precipitation_newmethod.npy',interpolated_precipitation_data)
elif elevation_exclusion == 'no' and vertical_interpolation == 'yes':
    np.save(precippath+'\interpolated_precipitation_standardgradient.npy',interpolated_precipitation_data)
elif elevation_exclusion == 'yes' and vertical_interpolation == 'no':
    np.save(metapath+'\interpolated_precipitation_version3.npy',interpolated_precipitation_data)
elif elevation_exclusion == 'yes' and vertical_interpolation == 'yes':
    np.save(metapath+'\interpolated_precipitation_version4.npy',interpolated_precipitation_data)


    