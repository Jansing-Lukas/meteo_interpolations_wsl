# -*- coding: utf-8 -*-
"""
Created on Thu May 16 14:29:22 2019

@author: lukas jansing

script to cluster the precipitation events
clusters can then be compared in terms of timing and magnitude

"""

#--------------------------------------------------------- 
# Import modules
#---------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import datetime
import csv
import matplotlib
import matplotlib.dates as dt
from dateutil.relativedelta import *

from import_data import import_lwf_precipitation_data,import_combiprecip


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
                       ,'Sempach','Surava_S','Surava_N','Visp','Vordemwald','Zürich']

treenetstation_id   = [1,2,3,4,5,6,7,8,1,10,11,12,13,14,7,16,17,7,1,20,21,1
                       ,20,24,25,26,27,28,29,30,31,32,28,34,27,36,37,25,39,40,41
                       ,42,43,20,45,46,47,48,49,50,51,52,53,54,54,56,14,14,14,60
                       ,43,43,63,63,65,65,67,1,1,1,71,72,73,74,24,24,77,78,79]

save_combiprecip = 'yes'

# Define which data to import from TreeNet
# 10minres data go further back in time
#treenetprecip_res = 'hourly'
treenetprecip_res = '10minres'

#---------------------------------------------------------
# !!!!!!  EDIT  !!!!!!!!! 
# Choose station here
# !!!!!!  EDIT  !!!!!!!!!
#---------------------------------------------------------
# Available stations:
# TreeNet sites of LWF: Jussy, Beatenberg, Lausanne, Lens, Neunkirch,
#                       Novaggio, Visp, Vordemwald, Schänis
# TreeNet site of Roman: Pfynwald
# TreeNet sites of IAP: Muri, Riehen, Grosswangen

treenetstation = 'Jussy'

station_id = treenetstation_id[treenetstationnames.index(treenetstation)]

#--------------------------------------------------------- 
# Define paths
#---------------------------------------------------------
stationprecippath = 'add path to MeteoSwiss hourly precipitation data'
combiprecippath   = 'add path to CombiPrecip'
treenetprecippath = 'add path to LWF precipitation data'
figpath           = 'add path to your figures folder'

#--------------------------------------------------------- 
# Load TreeNet/LWF data
#---------------------------------------------------------
precip = import_lwf_precipitation_data(treenetstation=treenetstation,treenetprecip_res=treenetprecip_res,\
                         path=treenetprecippath,process_schänis='no')

#--------------------------------------------------------- 
# Load interpolated precipitation data
#---------------------------------------------------------
print('import interpolated precipitation data')
interpolated_precipitation_data = np.load(stationprecippath+'\interpolated_precipitation_newmethod.npy')
precip['date_meteoswiss'] = interpolated_precipitation_data[0,:]

interp_ind = treenetstationnames.index(treenetstation)
precip['precip_interpolated'] = interpolated_precipitation_data[interp_ind+1,:].astype(float)

#--------------------------------------------------------- 
# Import processed combiprecip data
# Acess data with combiprecip_data[()][key]
#---------------------------------------------------------
combiprecip_data = import_combiprecip(combiprecippath=combiprecippath,processing_combiprecip='no',\
                                      save_combiprecip='yes')

#--------------------------------------------------------- 
# Extract location from combiprecip dataset
#---------------------------------------------------------
station_index = np.where(combiprecip_data[()]['combiprecip'][0,:] == station_id)[0][0]
precip['combiprecip'] = combiprecip_data[()]['combiprecip'][1:,station_index]
precip['date_combiprecip'] = np.array(combiprecip_data[()]['date_combiprecip'])

#--------------------------------------------------------- 
# Find out which timespans are covered by the 3 datasets
#---------------------------------------------------------
# Find latest starting date
if (precip['date_combiprecip'][0] > precip['date_meteoswiss'][0]) and \
   (precip['date_combiprecip'][0] > precip['date_hourly_UTC'][0]):
       earliest_date = precip['date_combiprecip'][0]
elif (precip['date_meteoswiss'][0] > precip['date_combiprecip'][0]) and \
     (precip['date_meteoswiss'][0] > precip['date_hourly_UTC'][0]):
        earliest_date = precip['date_meteoswiss'][0]
else:
    earliest_date = precip['date_hourly_UTC'][0]

# Find earliest ending date    
if (precip['date_combiprecip'][-1] < precip['date_meteoswiss'][-1]) and \
   (precip['date_combiprecip'][-1] < precip['date_hourly_UTC'][-1]):
       latest_date = precip['date_combiprecip'][-1]
elif (precip['date_meteoswiss'][-1] < precip['date_combiprecip'][-1]) and \
     (precip['date_meteoswiss'][-1] < precip['date_hourly_UTC'][-1]):
        latest_date = precip['date_meteoswiss'][-1]
else:
    latest_date = precip['date_hourly_UTC'][-1]
    
# Find indices of the season in the timeseries
wsl_ind_0    = np.where(precip['date_hourly_UTC'] == earliest_date)[0][0]
interp_ind_0 = np.where(precip['date_meteoswiss'] == earliest_date)[0][0]
combi_ind_0  = np.where(precip['date_combiprecip'] == earliest_date)[0][0]
wsl_ind_1    = np.where(precip['date_hourly_UTC'] == latest_date)[0][0]
interp_ind_1 = np.where(precip['date_meteoswiss'] == latest_date)[0][0]
combi_ind_1  = np.where(precip['date_combiprecip'] == latest_date)[0][0]

# Extract data
precip['precip_hourly']       = precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1]
precip['precip_interpolated'] = precip['precip_interpolated'][interp_ind_0:interp_ind_1+1]
precip['combiprecip']         = precip['combiprecip'][combi_ind_0:combi_ind_1+1]
precip['date_hourly_UTC']     = precip['date_hourly_UTC'][wsl_ind_0:wsl_ind_1+1]
precip['date_meteoswiss']     = precip['date_meteoswiss'][interp_ind_0:interp_ind_1+1]
precip['date_combiprecip']    = precip['date_combiprecip'][combi_ind_0:combi_ind_1+1]

#--------------------------------------------------------- 
# Cluster the major precipitation events in all datasets
#---------------------------------------------------------
# Initialise with zeros (all non-events should remain with index 0)
print('cluster precipitation events')
precip['clusterindex_wsl']         = np.zeros(np.shape(precip['precip_hourly']))
precip['clusterindex_meteoswiss']  = np.zeros(np.shape(precip['precip_interpolated']))
precip['clusterindex_combiprecip'] = np.zeros(np.shape(precip['combiprecip']))

# Index starting with 1
clusterindex_wsl         = 1
clusterindex_meteoswiss  = 1
clusterindex_combiprecip = 1

# Setting parameters
precip_threshold = 0.2  # minimal precip in order to index it
timegap = 10            # allowed time difference between two exceedences
count = 0               # counting events where there is a gap between events

# Identify first possible event
if np.any(precip['precip_hourly'][0:timegap+1] > precip_threshold) == True:
    precip['clusterindex_wsl'][0:timegap+1] = clusterindex_wsl
if np.any(precip['precip_interpolated'][0:timegap+1] > precip_threshold) == True:
    precip['clusterindex_meteoswiss'][0:timegap+1] = clusterindex_meteoswiss
if np.any(precip['combiprecip'][0:timegap+1] > precip_threshold) == True:
    precip['clusterindex_combiprecip'][0:timegap+1] = clusterindex_combiprecip

# Loop over whole timespan to cluster
for i in range (timegap,precip['precip_hourly'].shape[0]-timegap):
#for i in range (timegap,300):    
    # ---------------- Handling of WSL/LWF precipitation ----------------
    
    # Cases where precipitation events takes place
    if precip['precip_hourly'][i] > 0:
        
        # Check if it has rained in the x hours before
        if np.any(precip['clusterindex_wsl'][i-timegap:i] > 0) == True:
            # yes --> assign the same index (for the current event)
            precip['clusterindex_wsl'][i] = np.max(precip['clusterindex_wsl'][i-timegap:i])
            
        else:
            if np.any(precip['precip_hourly'][i:i+timegap+1] > precip_threshold) == True:
            # no --> assign new index (increase first by one)
                clusterindex_wsl += 1
                precip['clusterindex_wsl'][i] = clusterindex_wsl
         
    # Cases where no precipitation takes place
    else:
        
        # Check if it has rained in the x hours before
        if np.any(precip['precip_hourly'][i-timegap:i] > 0) == True:
            # yes --> find out if distance to next event is smaller than 5 hours
            last_rain_diff = i - np.argwhere(precip['precip_hourly'] > 0)
            last_rain_index = int(i - min([item for item in last_rain_diff if item > 0]))
            if i - last_rain_index <= timegap:
            #print(i)
            #print(last_rain)
            
                # Check out if it rains again in the time window considered
                if np.any(precip['precip_hourly'][last_rain_index+1:last_rain_index+timegap+2] > 0) == True:
                    # yes --> assign the same index as previous timestep
                    precip['clusterindex_wsl'][i] = precip['clusterindex_wsl'][i-1]
                    count+=1
                
    # ---------------- Handling of interpolated precipitation ----------------
    
    # Cases where precipitation events takes place
    if precip['precip_interpolated'][i] > 0:
        
        # Check if it has rained in the x hours before
        if np.any(precip['clusterindex_meteoswiss'][i-timegap:i] > 0) == True:
            # yes --> assign the same index (for the current event)
            precip['clusterindex_meteoswiss'][i] = np.max(precip['clusterindex_meteoswiss'][i-timegap:i])
            
        else:
            if np.any(precip['precip_interpolated'][i:i+timegap+1] > precip_threshold) == True:
            # no --> assign new index (increase first by one)
                clusterindex_meteoswiss += 1
                precip['clusterindex_meteoswiss'][i] = clusterindex_meteoswiss
         
    # Cases where no precipitation takes place
    else:
        
        # Check if it has rained in the x hours before
        if np.any(precip['precip_interpolated'][i-timegap:i] > 0) == True:
            # yes --> find out if distance to next event is smaller than 5 hours
            last_rain_diff = i - np.argwhere(precip['precip_interpolated'] > 0)
            last_rain_index = int(i - min([item for item in last_rain_diff if item > 0]))
            if i - last_rain_index <= timegap:
            #print(i)
            #print(last_rain)
            
                # Check out if it rains again in the time window considered
                if np.any(precip['precip_interpolated'][last_rain_index+1:last_rain_index+timegap+2] > 0) == True:
                    # yes --> assign the same index as previous timestep
                    precip['clusterindex_meteoswiss'][i] = precip['clusterindex_meteoswiss'][i-1]
                    count+=1
                
    # ---------------- Handling of CombiPrecip precipitation ----------------
    
    # Cases where precipitation events takes place
    if precip['combiprecip'][i] > 0:
        
        # Check if it has rained in the x hours before
        if np.any(precip['clusterindex_combiprecip'][i-timegap:i] > 0) == True:
            # yes --> assign the same index (for the current event)
            precip['clusterindex_combiprecip'][i] = np.max(precip['clusterindex_combiprecip'][i-timegap:i])
            
        else:
            if np.any(precip['combiprecip'][i:i+timegap+1] > precip_threshold) == True:
            # no --> assign new index (increase first by one)
                clusterindex_combiprecip += 1
                precip['clusterindex_combiprecip'][i] = clusterindex_combiprecip
         
    # Cases where no precipitation takes place
    else:
        
        # Check if it has rained in the x hours before
        if np.any(precip['combiprecip'][i-timegap:i] > 0) == True:
            # yes --> find out if distance to next event is smaller than 5 hours
            last_rain_diff = i - np.argwhere(precip['combiprecip'] > 0)
            last_rain_index = int(i - min([item for item in last_rain_diff if item > 0]))
            if i - last_rain_index <= timegap:
            #print(i)
            #print(last_rain)
            
                # Check out if it rains again in the time window considered
                if np.any(precip['combiprecip'][last_rain_index+1:last_rain_index+timegap+2] > 0) == True:
                    # yes --> assign the same index as previous timestep
                    precip['clusterindex_combiprecip'][i] = precip['clusterindex_combiprecip'][i-1]
                    count+=1

#--------------------------------------------------------- 
# Filter the cluster accordint to some criteria
#---------------------------------------------------------
clusterfiltering = 'yes'
#clusterfiltering = 'no'
                    
if clusterfiltering == 'yes':
    # Only keep clusters with total sum > 20 mm and one event over 2 mm/h
    event_threshold  = 20
    hourly_threshold = 2

    for i in range(int(np.min(precip['clusterindex_wsl'])),int(np.max(precip['clusterindex_wsl']))+1):
        # identify indices of the respective event
        event_index = np.argwhere(precip['clusterindex_wsl'] == i)[:,0]

        # check for the given condition --> set to zero if not fulfilled    
        if (np.nansum(precip['precip_hourly'][event_index]) < event_threshold) or \
           (np.any(precip['precip_hourly'][event_index] > hourly_threshold) == False):
            precip['clusterindex_wsl'][event_index] = 0
           
    for i in range(int(np.min(precip['clusterindex_meteoswiss'])),int(np.max(precip['clusterindex_meteoswiss']))+1):
        # identify indices of the respective event
        event_index = np.argwhere(precip['clusterindex_meteoswiss'] == i)[:,0]

        # check for the given condition --> set to zero if not fulfilled    
        if (np.nansum(precip['precip_interpolated'][event_index]) < event_threshold) or \
           (np.any(precip['precip_interpolated'][event_index] > hourly_threshold) == False):
               precip['clusterindex_meteoswiss'][event_index] = 0
           
    for i in range(int(np.min(precip['clusterindex_combiprecip'])),int(np.max(precip['clusterindex_combiprecip']))+1):
        # identify indices of the respective event
        event_index = np.argwhere(precip['clusterindex_combiprecip'] == i)[:,0]

        # check for the given condition --> set to zero if not fulfilled    
        if (np.nansum(precip['combiprecip'][event_index]) < event_threshold) or \
           (np.any(precip['combiprecip'][event_index] > hourly_threshold) == False):
               precip['clusterindex_combiprecip'][event_index] = 0

                
#--------------------------------------------------------- 
# Plot cluster boundaries into monthly plots
#---------------------------------------------------------             
# Add up to next month
year = earliest_date.year
if earliest_date.day != 1 or earliest_date.hour != 0:
    month = earliest_date.month+1
    earliest_date = datetime.datetime(year,month,1,0)
    
# Subtract to previous month
year = latest_date.year
if (latest_date.day != 31 or 30) or latest_date.hour != 23:
    latest_date = latest_date - relativedelta(months=1,day=31,hour=23)
    
# Create the datelist to loop over
datelist_months = [earliest_date]
nowdate = earliest_date
while nowdate <= latest_date:
    nowdate = nowdate+relativedelta(months=+1)
    datelist_months.append(nowdate)
    
# Loop over the whole timespan in order to create plots for each month
matplotlib.rcParams.update({'font.size': 20})
plt.close('all')

for nowdate in datelist_months:
    
    print(nowdate.strftime('%b %Y'))
    print('plotting')
    plt.close('all')
    f,axarr=plt.subplots(3,1)
    f.set_size_inches(25, 20)
    
    #Get latest instant of the month
    month_end = nowdate+relativedelta(day=31,hour=23)
 
    # Find indices of the month in the timeseries
    wsl_ind_0    = np.where(precip['date_hourly_UTC'] == nowdate)[0][0]
    interp_ind_0 = np.where(precip['date_meteoswiss'] == nowdate)[0][0]
    combi_ind_0  = np.where(precip['date_combiprecip'] == nowdate)[0][0]
    wsl_ind_1    = np.where(precip['date_hourly_UTC'] == month_end)[0][0]
    interp_ind_1 = np.where(precip['date_meteoswiss'] == month_end)[0][0]
    combi_ind_1  = np.where(precip['date_combiprecip'] == month_end)[0][0]

    # Get indices of wsl nans and meteoswiss nans
    wsl_nans         = np.argwhere(np.isnan(precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1]))[:,0]
    meteoswiss_nans  = np.argwhere(np.isnan(precip['precip_interpolated'][interp_ind_0:interp_ind_1+1]))[:,0]
    combiprecip_nans = np.argwhere(np.isnan(precip['combiprecip'][combi_ind_0:combi_ind_1+1]))[:,0]

    # Plot
    axarr[0].bar(precip['date_hourly_UTC'][wsl_ind_0:wsl_ind_1+1],\
             precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1],color='blue',width=0.02,align='center')
    axarr[0].scatter(precip['date_hourly_UTC'][wsl_ind_0:wsl_ind_1+1][wsl_nans],\
             np.zeros((len(wsl_nans))),color='red',s=20)
    ymin1,ymax1 = axarr[0].get_ylim()
    axarr[0].set_ylabel('Precipitation [mm/h]',fontsize=20)
    axarr[0].set_title(treenetstation+' WSL',fontsize=25)

    axarr[1].bar(precip['date_meteoswiss'][interp_ind_0:interp_ind_1+1],\
             precip['precip_interpolated'][interp_ind_0:interp_ind_1+1],color='blue',width=0.02,align='center')
    axarr[1].scatter(precip['date_meteoswiss'][interp_ind_0:interp_ind_1+1][meteoswiss_nans],\
             np.zeros((len(meteoswiss_nans))),color='red',s=20)
    ymin2,ymax2 = axarr[1].get_ylim()
    axarr[1].set_ylabel('Precipitation [mm/h]',fontsize=20)
    axarr[1].set_title(treenetstation+' Interpolation',fontsize=25)

    axarr[2].bar(precip['date_combiprecip'][combi_ind_0:combi_ind_1+1],\
             precip['combiprecip'][combi_ind_0:combi_ind_1+1],color='blue',width=0.02,align='center')
    axarr[2].scatter(precip['date_combiprecip'][combi_ind_0:combi_ind_1+1][combiprecip_nans],\
             np.zeros((len(combiprecip_nans))),color='red',s=20)
    ymin3,ymax3 = axarr[2].get_ylim()
    axarr[2].set_ylabel('Precipitation [mm/h]',fontsize=20)
    axarr[2].set_title('Combi Precip',fontsize=25)

    overall_min = 0
    overall_max = np.max([ymax1,ymax2,ymax3])
    axarr[0].set_ylim([overall_min,overall_max])
    axarr[1].set_ylim([overall_min,overall_max])
    axarr[2].set_ylim([overall_min,overall_max])
    
    for i in range(wsl_ind_0,wsl_ind_1):
        if precip['clusterindex_wsl'][i+1] > precip['clusterindex_wsl'][i]:
            axarr[0].axvline(x=precip['date_hourly_UTC'][i],ymin=overall_min,ymax=overall_max,color='r',lw=3)
        if precip['clusterindex_wsl'][i+1] < precip['clusterindex_wsl'][i]:
            axarr[0].axvline(x=precip['date_hourly_UTC'][i],ymin=overall_min,ymax=overall_max,color='g',lw=3)
    for i in range(interp_ind_0,interp_ind_1):
        if precip['clusterindex_meteoswiss'][i+1] > precip['clusterindex_meteoswiss'][i]:
            axarr[1].axvline(x=precip['date_meteoswiss'][i],ymin=overall_min,ymax=overall_max,color='r',lw=3)
        if precip['clusterindex_meteoswiss'][i+1] < precip['clusterindex_meteoswiss'][i]:
            axarr[1].axvline(x=precip['date_meteoswiss'][i],ymin=overall_min,ymax=overall_max,color='g',lw=3)
    for i in range(combi_ind_0,combi_ind_1):
        if precip['clusterindex_combiprecip'][i+1] > precip['clusterindex_combiprecip'][i]:
            axarr[2].axvline(x=precip['date_combiprecip'][i],ymin=overall_min,ymax=overall_max,color='r',lw=3)
        if precip['clusterindex_combiprecip'][i+1] < precip['clusterindex_combiprecip'][i]:
            axarr[2].axvline(x=precip['date_combiprecip'][i],ymin=overall_min,ymax=overall_max,color='g',lw=3)

    plt.suptitle(precip['date_hourly_UTC'][wsl_ind_0].strftime('%b %Y'),fontsize=40)

    # Format x-axis and grid
    datelist = make_datelist(precip['date_hourly_UTC'][wsl_ind_0],precip['date_hourly_UTC'][wsl_ind_1],1)
    for axx in axarr:
        axx.xaxis_date()
        axx.xaxis.set_ticks(datelist[0::24*7],minor=False)
        axx.xaxis.set_major_formatter(dt.DateFormatter('%d-%m-%Y %H'))
        axx.xaxis.grid(True, which='major', color='0.5',alpha=0.7, linestyle='--',lw=1.5)

    saveas='\precip_clustering_'+treenetstation+'_'
    plt.savefig(figpath+saveas+nowdate.strftime('%Y_%m')+'.png',bbox_inches='tight')
    
#--------------------------------------------------------- 
# Characterise/compare the different clusters
#---------------------------------------------------------
# Plot the major events in order to compare them to each other
matplotlib.rcParams.update({'font.size': 20})

# Loop over indexlist of CombiPrecip
for i in np.unique(precip['clusterindex_combiprecip'])[1:]:
    
    plt.close('all')
    f,axarr=plt.subplots(3,1)
    f.set_size_inches(10, 20)
    
    event_index = np.argwhere(precip['clusterindex_combiprecip'] == i)[:,0]
    event_number = i
    print('plotting event '+str(i))
    
    if int(np.min(event_index)) > 24:
        event_min = int(np.min(event_index)) - 24
    else:
        event_min = 0    
    event_max= int(np.max(event_index)) + 24
    
    # Get indices of wsl nans and meteoswiss nans
    wsl_nans         = np.argwhere(np.isnan(precip['precip_hourly'][event_min:event_max+1]))[:,0]
    meteoswiss_nans  = np.argwhere(np.isnan(precip['precip_interpolated'][event_min:event_max+1]))[:,0]
    combiprecip_nans = np.argwhere(np.isnan(precip['combiprecip'][event_min:event_max+1]))[:,0]
    
    axarr[0].bar(precip['date_hourly_UTC'][event_min:event_max+1],\
             precip['precip_hourly'][event_min:event_max+1],color='blue',width=0.02,align='center')
    axarr[0].scatter(precip['date_hourly_UTC'][event_min:event_max+1][wsl_nans],\
             np.zeros((len(wsl_nans))),color='red',s=20)
    ymin1,ymax1 = axarr[0].get_ylim()
    axarr[0].set_ylabel('Precipitation [mm/h]',fontsize=20)
    axarr[0].set_title(treenetstation+' WSL',fontsize=25)

    axarr[1].bar(precip['date_meteoswiss'][event_min:event_max+1],\
             precip['precip_interpolated'][event_min:event_max+1],color='blue',width=0.02,align='center')
    axarr[1].scatter(precip['date_meteoswiss'][event_min:event_max+1][meteoswiss_nans],\
             np.zeros((len(meteoswiss_nans))),color='red',s=20)
    ymin2,ymax2 = axarr[1].get_ylim()
    axarr[1].set_ylabel('Precipitation [mm/h]',fontsize=20)
    axarr[1].set_title(treenetstation+' Interpolation',fontsize=25)

    axarr[2].bar(precip['date_combiprecip'][event_min:event_max+1],\
             precip['combiprecip'][event_min:event_max+1],color='blue',width=0.02,align='center')
    axarr[2].scatter(precip['date_combiprecip'][event_min:event_max+1][combiprecip_nans],\
             np.zeros((len(combiprecip_nans))),color='red',s=20)
    ymin3,ymax3 = axarr[2].get_ylim()
    axarr[2].set_ylabel('Precipitation [mm/h]',fontsize=20)
    axarr[2].set_title('Combi Precip',fontsize=25)

    overall_min = 0
    overall_max = np.max([ymax1,ymax2,ymax3])
    axarr[0].set_ylim([overall_min,overall_max])
    axarr[1].set_ylim([overall_min,overall_max])
    axarr[2].set_ylim([overall_min,overall_max])
    
    axarr[0].text(precip['date_hourly_UTC'][event_min],overall_max-0.1*overall_max,\
                  str('%.1f' % np.nansum(precip['precip_hourly'][event_min:event_max+1]))+' mm')
    axarr[1].text(precip['date_meteoswiss'][event_min],overall_max-0.1*overall_max,\
                  str('%.1f' % np.nansum(precip['precip_interpolated'][event_min:event_max+1]))+' mm')
    axarr[2].text(precip['date_combiprecip'][event_min],overall_max-0.1*overall_max,\
                  str('%.1f' % np.nansum(precip['combiprecip'][event_min:event_max+1]))+' mm')
    
    for i in range(event_min,event_max):
        if precip['clusterindex_wsl'][i+1] > precip['clusterindex_wsl'][i]:
            axarr[0].axvline(x=precip['date_hourly_UTC'][i],ymin=overall_min,ymax=overall_max,color='r',lw=3)
        if precip['clusterindex_wsl'][i+1] < precip['clusterindex_wsl'][i]:
            axarr[0].axvline(x=precip['date_hourly_UTC'][i],ymin=overall_min,ymax=overall_max,color='g',lw=3)
    for i in range(event_min,event_max):
        if precip['clusterindex_meteoswiss'][i+1] > precip['clusterindex_meteoswiss'][i]:
            axarr[1].axvline(x=precip['date_meteoswiss'][i],ymin=overall_min,ymax=overall_max,color='r',lw=3)
        if precip['clusterindex_meteoswiss'][i+1] < precip['clusterindex_meteoswiss'][i]:
            axarr[1].axvline(x=precip['date_meteoswiss'][i],ymin=overall_min,ymax=overall_max,color='g',lw=3)
    for i in range(event_min,event_max):
        if precip['clusterindex_combiprecip'][i+1] > precip['clusterindex_combiprecip'][i]:
            axarr[2].axvline(x=precip['date_combiprecip'][i],ymin=overall_min,ymax=overall_max,color='r',lw=3)
        if precip['clusterindex_combiprecip'][i+1] < precip['clusterindex_combiprecip'][i]:
            axarr[2].axvline(x=precip['date_combiprecip'][i],ymin=overall_min,ymax=overall_max,color='g',lw=3)
            
    # Format x-axis and grid
    datelist = make_datelist(precip['date_hourly_UTC'][event_min],precip['date_hourly_UTC'][event_max],1)
    for axx in axarr:
        axx.xaxis_date()
        axx.xaxis.set_ticks(datelist[0::24*2],minor=False)
        axx.xaxis.set_major_formatter(dt.DateFormatter('%d-%m-%Y %H'))
        axx.xaxis.grid(True, which='major', color='0.5',alpha=0.7, linestyle='--',lw=1.5)
        
    saveas='\precip_major_events_'+treenetstation+'_'
    plt.savefig(figpath+saveas+'event_combi_'+str(int(event_number))+'.png',bbox_inches='tight')
            
# Loop over indexlist of WSL/LWF data
for i in np.unique(precip['clusterindex_wsl'])[1:]:
    
    event_index = np.argwhere(precip['clusterindex_wsl'] == i)[:,0]
    event_number = i
    
    if int(np.min(event_index)) > 24:
        event_min = int(np.min(event_index)) - 24
    else:
        event_min = 0
    event_max= int(np.max(event_index)) + 24
    
    if np.any(precip['clusterindex_combiprecip'][event_min:event_max+1] != 0) == False:
        
        print('plotting event '+str(i))
        plt.close('all')
        f,axarr=plt.subplots(3,1)
        f.set_size_inches(10, 20)

        # Get indices of wsl nans and meteoswiss nans
        wsl_nans         = np.argwhere(np.isnan(precip['precip_hourly'][event_min:event_max+1]))[:,0]
        meteoswiss_nans  = np.argwhere(np.isnan(precip['precip_interpolated'][event_min:event_max+1]))[:,0]
        combiprecip_nans = np.argwhere(np.isnan(precip['combiprecip'][event_min:event_max+1]))[:,0]
    
        axarr[0].bar(precip['date_hourly_UTC'][event_min:event_max+1],\
                 precip['precip_hourly'][event_min:event_max+1],color='blue',width=0.02,align='center')
        axarr[0].scatter(precip['date_hourly_UTC'][event_min:event_max+1][wsl_nans],\
                 np.zeros((len(wsl_nans))),color='red',s=20)
        ymin1,ymax1 = axarr[0].get_ylim()
        axarr[0].set_ylabel('Precipitation [mm/h]',fontsize=20)
        axarr[0].set_title(treenetstation+' WSL',fontsize=25)

        axarr[1].bar(precip['date_meteoswiss'][event_min:event_max+1],\
                 precip['precip_interpolated'][event_min:event_max+1],color='blue',width=0.02,align='center')
        axarr[1].scatter(precip['date_meteoswiss'][event_min:event_max+1][meteoswiss_nans],\
                 np.zeros((len(meteoswiss_nans))),color='red',s=20)
        ymin2,ymax2 = axarr[1].get_ylim()
        axarr[1].set_ylabel('Precipitation [mm/h]',fontsize=20)
        axarr[1].set_title(treenetstation+' Interpolation',fontsize=25)
        
        axarr[2].bar(precip['date_combiprecip'][event_min:event_max+1],\
                 precip['combiprecip'][event_min:event_max+1],color='blue',width=0.02,align='center')
        axarr[2].scatter(precip['date_combiprecip'][event_min:event_max+1][combiprecip_nans],\
                 np.zeros((len(combiprecip_nans))),color='red',s=20)
        ymin3,ymax3 = axarr[2].get_ylim()
        axarr[2].set_ylabel('Precipitation [mm/h]',fontsize=20)
        axarr[2].set_title('Combi Precip',fontsize=25)

        overall_min = 0
        overall_max = np.max([ymax1,ymax2,ymax3])
        axarr[0].set_ylim([overall_min,overall_max])
        axarr[1].set_ylim([overall_min,overall_max])
        axarr[2].set_ylim([overall_min,overall_max])
            
        axarr[0].text(precip['date_hourly_UTC'][event_min],overall_max-0.1*overall_max,\
                      str('%.1f' % np.nansum(precip['precip_hourly'][event_min:event_max+1]))+' mm')
        axarr[1].text(precip['date_meteoswiss'][event_min],overall_max-0.1*overall_max,\
                      str('%.1f' % np.nansum(precip['precip_interpolated'][event_min:event_max+1]))+' mm')
        axarr[2].text(precip['date_combiprecip'][event_min],overall_max-0.1*overall_max,\
                      str('%.1f' % np.nansum(precip['combiprecip'][event_min:event_max+1]))+' mm')
    
        for i in range(event_min,event_max):
            if precip['clusterindex_wsl'][i+1] > precip['clusterindex_wsl'][i]:
                axarr[0].axvline(x=precip['date_hourly_UTC'][i],ymin=overall_min,ymax=overall_max,color='r',lw=3)
            if precip['clusterindex_wsl'][i+1] < precip['clusterindex_wsl'][i]:
                axarr[0].axvline(x=precip['date_hourly_UTC'][i],ymin=overall_min,ymax=overall_max,color='g',lw=3)
        for i in range(event_min,event_max):
            if precip['clusterindex_meteoswiss'][i+1] > precip['clusterindex_meteoswiss'][i]:
                axarr[1].axvline(x=precip['date_meteoswiss'][i],ymin=overall_min,ymax=overall_max,color='r',lw=3)
            if precip['clusterindex_meteoswiss'][i+1] < precip['clusterindex_meteoswiss'][i]:
                axarr[1].axvline(x=precip['date_meteoswiss'][i],ymin=overall_min,ymax=overall_max,color='g',lw=3)
        for i in range(event_min,event_max):
            if precip['clusterindex_combiprecip'][i+1] > precip['clusterindex_combiprecip'][i]:
                axarr[2].axvline(x=precip['date_combiprecip'][i],ymin=overall_min,ymax=overall_max,color='r',lw=3)
            if precip['clusterindex_combiprecip'][i+1] < precip['clusterindex_combiprecip'][i]:
                axarr[2].axvline(x=precip['date_combiprecip'][i],ymin=overall_min,ymax=overall_max,color='g',lw=3)
            
        # Format x-axis and grid
        datelist = make_datelist(precip['date_hourly_UTC'][event_min],precip['date_hourly_UTC'][event_max],1)
        for axx in axarr:
            axx.xaxis_date()
            axx.xaxis.set_ticks(datelist[0::24*2],minor=False)
            axx.xaxis.set_major_formatter(dt.DateFormatter('%d-%m-%Y %H'))
            axx.xaxis.grid(True, which='major', color='0.5',alpha=0.7, linestyle='--',lw=1.5)
        
        saveas='\precip_major_events_'+treenetstation+'_'
        plt.savefig(figpath+saveas+'event_wsl_'+str(int(event_number))+'.png',bbox_inches='tight')

# Loop over indexlist of Interpolated data
for i in np.unique(precip['clusterindex_meteoswiss'])[1:]:
    
    event_index = np.argwhere(precip['clusterindex_meteoswiss'] == i)[:,0]
    event_number = i
    
    if int(np.min(event_index)) > 24:
        event_min = int(np.min(event_index)) - 24
    else:
        event_min = 0
    event_max= int(np.max(event_index)) + 24
    
    if np.any(precip['clusterindex_combiprecip'][event_min:event_max+1] != 0) == False and \
       np.any(precip['clusterindex_wsl'][event_min:event_max+1] != 0) == False:
           
        print('plotting event '+str(i))
        plt.close('all')
        f,axarr=plt.subplots(3,1)
        f.set_size_inches(10, 20)

        # Get indices of wsl nans and meteoswiss nans
        wsl_nans         = np.argwhere(np.isnan(precip['precip_hourly'][event_min:event_max+1]))[:,0]
        meteoswiss_nans  = np.argwhere(np.isnan(precip['precip_interpolated'][event_min:event_max+1]))[:,0]
        combiprecip_nans = np.argwhere(np.isnan(precip['combiprecip'][event_min:event_max+1]))[:,0]
    
        axarr[0].bar(precip['date_hourly_UTC'][event_min:event_max+1],\
                 precip['precip_hourly'][event_min:event_max+1],color='blue',width=0.02,align='center')
        axarr[0].scatter(precip['date_hourly_UTC'][event_min:event_max+1][wsl_nans],\
                 np.zeros((len(wsl_nans))),color='red',s=20)
        ymin1,ymax1 = axarr[0].get_ylim()
        axarr[0].set_ylabel('Precipitation [mm/h]',fontsize=20)
        axarr[0].set_title(treenetstation+' WSL',fontsize=25)

        axarr[1].bar(precip['date_meteoswiss'][event_min:event_max+1],\
                 precip['precip_interpolated'][event_min:event_max+1],color='blue',width=0.02,align='center')
        axarr[1].scatter(precip['date_meteoswiss'][event_min:event_max+1][meteoswiss_nans],\
                 np.zeros((len(meteoswiss_nans))),color='red',s=20)
        ymin2,ymax2 = axarr[1].get_ylim()
        axarr[1].set_ylabel('Precipitation [mm/h]',fontsize=20)
        axarr[1].set_title(treenetstation+' Interpolation',fontsize=25)
        
        axarr[2].bar(precip['date_combiprecip'][event_min:event_max+1],\
                 precip['combiprecip'][event_min:event_max+1],color='blue',width=0.02,align='center')
        axarr[2].scatter(precip['date_combiprecip'][event_min:event_max+1][combiprecip_nans],\
                 np.zeros((len(combiprecip_nans))),color='red',s=20)
        ymin3,ymax3 = axarr[2].get_ylim()
        axarr[2].set_ylabel('Precipitation [mm/h]',fontsize=20)
        axarr[2].set_title('Combi Precip',fontsize=25)

        overall_min = 0
        overall_max = np.max([ymax1,ymax2,ymax3])
        axarr[0].set_ylim([overall_min,overall_max])
        axarr[1].set_ylim([overall_min,overall_max])
        axarr[2].set_ylim([overall_min,overall_max])
            
        axarr[0].text(precip['date_hourly_UTC'][event_min],overall_max-0.1*overall_max,\
                      str('%.1f' % np.nansum(precip['precip_hourly'][event_min:event_max+1]))+' mm')
        axarr[1].text(precip['date_meteoswiss'][event_min],overall_max-0.1*overall_max,\
                      str('%.1f' % np.nansum(precip['precip_interpolated'][event_min:event_max+1]))+' mm')
        axarr[2].text(precip['date_combiprecip'][event_min],overall_max-0.1*overall_max,\
                      str('%.1f' % np.nansum(precip['combiprecip'][event_min:event_max+1]))+' mm')
    
        for i in range(event_min,event_max):
            if precip['clusterindex_wsl'][i+1] > precip['clusterindex_wsl'][i]:
                axarr[0].axvline(x=precip['date_hourly_UTC'][i],ymin=overall_min,ymax=overall_max,color='r',lw=3)
            if precip['clusterindex_wsl'][i+1] < precip['clusterindex_wsl'][i]:
                axarr[0].axvline(x=precip['date_hourly_UTC'][i],ymin=overall_min,ymax=overall_max,color='g',lw=3)
        for i in range(event_min,event_max):
            if precip['clusterindex_meteoswiss'][i+1] > precip['clusterindex_meteoswiss'][i]:
                axarr[1].axvline(x=precip['date_meteoswiss'][i],ymin=overall_min,ymax=overall_max,color='r',lw=3)
            if precip['clusterindex_meteoswiss'][i+1] < precip['clusterindex_meteoswiss'][i]:
                axarr[1].axvline(x=precip['date_meteoswiss'][i],ymin=overall_min,ymax=overall_max,color='g',lw=3)
        for i in range(event_min,event_max):
            if precip['clusterindex_combiprecip'][i+1] > precip['clusterindex_combiprecip'][i]:
                axarr[2].axvline(x=precip['date_combiprecip'][i],ymin=overall_min,ymax=overall_max,color='r',lw=3)
            if precip['clusterindex_combiprecip'][i+1] < precip['clusterindex_combiprecip'][i]:
                axarr[2].axvline(x=precip['date_combiprecip'][i],ymin=overall_min,ymax=overall_max,color='g',lw=3)
            
        # Format x-axis and grid
        datelist = make_datelist(precip['date_hourly_UTC'][event_min],precip['date_hourly_UTC'][event_max],1)
        for axx in axarr:
            axx.xaxis_date()
            axx.xaxis.set_ticks(datelist[0::24*2],minor=False)
            axx.xaxis.set_major_formatter(dt.DateFormatter('%d-%m-%Y %H'))
            axx.xaxis.grid(True, which='major', color='0.5',alpha=0.7, linestyle='--',lw=1.5)
        
        saveas='\precip_major_events_'+treenetstation+'_'
        plt.savefig(figpath+saveas+'event_interpolated_'+str(int(event_number))+'.png',bbox_inches='tight')            

    
#--------------------------------------------------------- 
# Determine mean time gap
#---------------------------------------------------------
index_before = 0
time_diff = []
for i in range(0,precip['clusterindex_wsl'].shape[0]):
    if precip['clusterindex_wsl'][i] > 0 and \
       precip['clusterindex_wsl'][i] != index_before:
           index_before = precip['clusterindex_wsl'][i]
           if np.any(precip['clusterindex_combiprecip'][i-12:i+13] > 0) == True:
               index_combiprecip = np.min(np.nonzero(precip['clusterindex_combiprecip'][i-24:i+25])[0])
               time_diff.append(index_combiprecip - 12)

time_diff = np.array(time_diff)        
    
index_before = 0
time_diff = []
for i in range(0,precip['clusterindex_wsl'].shape[0]):
    if precip['clusterindex_wsl'][i] > 0 and \
       precip['clusterindex_wsl'][i] != index_before:
           index_before = precip['clusterindex_wsl'][i]
           if np.any(precip['clusterindex_meteoswiss'][i-12:i+13] > 0) == True:
               index_combiprecip = np.min(np.nonzero(precip['clusterindex_meteoswiss'][i-24:i+25])[0])
               time_diff.append(index_combiprecip - 12)

time_diff = np.array(time_diff)        


