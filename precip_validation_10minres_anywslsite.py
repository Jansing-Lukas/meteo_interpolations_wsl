# -*- coding: utf-8 -*-
"""
Created on Wed May 22 2019

@author: lukas jansing

Script in order to validate the interpolated precipitation data

Compare 10minres data to 10minres data of LWF

Compare hourly sum of LWF and interpolation to CombiPrecip

"""

#--------------------------------------------------------- 
# Import modules
#---------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import datetime
from datetime import timedelta
import csv
import matplotlib
import matplotlib.dates as dt
from dateutil.relativedelta import *
import gc

from import_data import import_lwf_precipitation_data,import_combiprecip
from functions import make_datelist

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

treenetstation = 'Vordemwald'

# needed for selecting proper CombiPrecip timeseries
station_id = treenetstation_id[treenetstationnames.index(treenetstation)]

#--------------------------------------------------------- 
# Define paths
#---------------------------------------------------------
stationprecippath = 'add path to MeteoSwiss 10min precipitation data'
combiprecippath   = 'add path to CombiPrecip'
treenetprecippath = 'add path to LWF precipitation data'
figpath           = 'add path to your figures folder'

#--------------------------------------------------------- 
# Load TreeNet/LWF data
#---------------------------------------------------------
precip = import_lwf_precipitation_data(treenetstation=treenetstation,treenetprecip_res='10minres',\
                                       path=treenetprecippath,process_schänis='no')

#--------------------------------------------------------- 
# Load interpolated precipitation data
#---------------------------------------------------------
print('import interpolated precipitation data in 10min-resolution')
#interpolated_precipitation_data = np.load(stationprecippath+'\interpolated_precipitation_10minres.npy')
interpolated_precipitation_data = np.load(stationprecippath+'\interpolated_precipitation_10minres_newmethod.npy')
precip['date_meteoswiss'] = interpolated_precipitation_data[0,:]

interp_ind = treenetstationnames.index(treenetstation)
precip['precip_interpolated_10minres'] = interpolated_precipitation_data[interp_ind+1,:].astype(float)

#--------------------------------------------------------- 
# Calculate hourly sums of interpolated data
#---------------------------------------------------------
# First check if first and last entry are start and finish of an hour
# identify the indices to read out correct data
for i in range(0,6):
    if precip['date_meteoswiss'][i].minute == 10:
        starting_ind = i
        break
    else:
        continue

for i in range(-1,-7,-1):
    if precip['date_meteoswiss'][i].minute == 0:
        ending_ind = i
        break
    else:
        continue
    
if ending_ind != -1:   
    interpolated_precip_shape = precip['precip_interpolated_10minres'][starting_ind:ending_ind+1].shape[0]
    precip['precip_interpolated_10minres_reshaped']  = np.reshape(precip['precip_interpolated_10minres'][starting_ind:ending_ind+1],\
                                           (int(interpolated_precip_shape/6),6))
    precip['date_interpolated_10minres_reshaped'] = np.reshape(precip['date_meteoswiss'][starting_ind:ending_ind+1],\
                                           (int(interpolated_precip_shape/6),6))
else:
    interpolated_precip_shape = precip['precip_interpolated_10minres'][starting_ind:].shape[0]
    precip['precip_interpolated_10minres_reshaped']  = np.reshape(precip['precip_interpolated_10minres'][starting_ind:],\
                                           (int(interpolated_precip_shape/6),6))
    precip['date_interpolated_10minres_reshaped'] = np.reshape(precip['date_meteoswiss'][starting_ind:],\
                                           (int(interpolated_precip_shape/6),6))


# Calculate hourly sums
precip['precip_interpolated_hourly'] = []
for i in range(0,int(interpolated_precip_shape/6)):
    precip['precip_interpolated_hourly'].append(sum(precip['precip_interpolated_10minres_reshaped'][i,:]))

precip['precip_interpolated_hourly'] = np.array(precip['precip_interpolated_hourly'])

# Create hourly datelist
# Note: Date always refers to precip from the last hour!!!
precip['date_hourly_meteoswiss'] = make_datelist(precip['date_interpolated_10minres_reshaped'][0,-1],\
                  precip['date_interpolated_10minres_reshaped'][-1,-1],1)
precip['date_hourly_meteoswiss'] = np.array(precip['date_hourly_meteoswiss'])

#--------------------------------------------------------- 
# Import processed combiprecip data
#---------------------------------------------------------
print('import processed combiprecip data')
combiprecip_data = np.load(combiprecippath+'\combiprecip_processed.npy')

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
   (precip['date_combiprecip'][0] > precip['date_10minres_UTC'][0]):
       earliest_date = precip['date_combiprecip'][0]
elif (precip['date_meteoswiss'][0] > precip['date_combiprecip'][0]) and \
     (precip['date_meteoswiss'][0] > precip['date_10minres_UTC'][0]):
        earliest_date = precip['date_meteoswiss'][0]
else:
    earliest_date = precip['date_10minres_UTC'][0]

# Add up to next month
year = earliest_date.year
if earliest_date.day != 1 or earliest_date.hour != 0:
    month = earliest_date.month+1
    earliest_date = datetime.datetime(year,month,1,0)

# Find earliest ending date    
if (precip['date_combiprecip'][-1] < precip['date_meteoswiss'][-1]) and \
   (precip['date_combiprecip'][-1] < precip['date_10minres_UTC'][-1]):
       latest_date = precip['date_combiprecip'][-1]
elif (precip['date_meteoswiss'][-1] < precip['date_combiprecip'][-1]) and \
     (precip['date_meteoswiss'][-1] < precip['date_10minres_UTC'][-1]):
        latest_date = precip['date_meteoswiss'][-1]
else:
    latest_date = precip['date_10minres_UTC'][-1]

# Subtract to previous month
year = latest_date.year
if (latest_date.day != 31 or 30) or latest_date.hour != 23:
    latest_date = latest_date - relativedelta(months=1,day=31,hour=23)
    
# Create the datelist to loop over
datelist_months = []
nowdate = earliest_date
while nowdate <= latest_date:
    datelist_months.append(nowdate)
    nowdate = nowdate+relativedelta(months=+1)
    
#--------------------------------------------------------- 
# Loop over the whole timespan in order to create plots for each month
#---------------------------------------------------------
matplotlib.rcParams.update({'font.size': 20})

for nowdate in datelist_months:
    
    print(nowdate.strftime('%b %Y'))
    print('plotting')
    f,axarr=plt.subplots(3,1)
    f.set_size_inches(30, 24)
    
    #Get latest instant of the month
    month_end = nowdate+relativedelta(day=31,hour=23)
 
    # Find indices of the month in the timeseries
    wsl_ind_0    = np.where(precip['date_10minres_UTC'] == nowdate)[0][0]
    interp_ind_0 = np.where(precip['date_meteoswiss'] == nowdate)[0][0]
    combi_ind_0  = np.where(precip['date_combiprecip'] == nowdate)[0][0] 
    wsl_ind_1    = np.where(precip['date_10minres_UTC'] == month_end)[0][0]
    interp_ind_1 = np.where(precip['date_meteoswiss'] == month_end)[0][0]
    combi_ind_1  = np.where(precip['date_combiprecip'] == month_end)[0][0]
    
    if nowdate != earliest_date:
        wsl_ind_2    = np.where(precip['date_hourly_UTC'] == nowdate)[0][0]
        interp_ind_2 = np.where(precip['date_hourly_meteoswiss'] == nowdate)[0][0]
        wsl_ind_3    = np.where(precip['date_hourly_UTC'] == month_end)[0][0]
        interp_ind_3 = np.where(precip['date_hourly_meteoswiss'] == month_end)[0][0]
    
    if nowdate == earliest_date:    # special treatment of first month
        if earliest_date not in precip['date_hourly_UTC']:
            wsl_ind_2    = 0
        else:
            wsl_ind_2 = np.where(precip['date_hourly_UTC'] == nowdate)[0][0]
            
        if earliest_date not in precip['date_hourly_meteoswiss']:
            interp_ind_2    = 0
        else:
            interp_ind_2 = np.where(precip['date_hourly_meteoswiss'] == nowdate)[0][0]
            
        wsl_ind_3    = np.where(precip['date_hourly_UTC'] == month_end)[0][0]
        interp_ind_3 = np.where(precip['date_hourly_meteoswiss'] == month_end)[0][0]
        
    # Get indices of wsl nans and meteoswiss nans
    wsl_nans         = np.argwhere(np.isnan(precip['precip_10minres'][wsl_ind_0:wsl_ind_1+1]))[:,0]
    meteoswiss_nans  = np.argwhere(np.isnan(precip['precip_interpolated_10minres'][interp_ind_0:interp_ind_1+1]))[:,0]
    combiprecip_nans = np.argwhere(np.isnan(precip['combiprecip'][combi_ind_0:combi_ind_1+1]))[:,0]

    # Plot
    axarr[0].bar(precip['date_10minres_UTC'][wsl_ind_0:wsl_ind_1+1],\
             precip['precip_10minres'][wsl_ind_0:wsl_ind_1+1],color='blue',width=0.005,align='center',label='10min')
    axarr[0].scatter(precip['date_10minres_UTC'][wsl_ind_0:wsl_ind_1+1][wsl_nans],\
             np.zeros((len(wsl_nans))),color='red',s=20)
    axarr[0].plot(precip['date_hourly_UTC'][wsl_ind_2:wsl_ind_3+1],\
             precip['precip_hourly'][wsl_ind_2:wsl_ind_3+1],ls='--',color='orange',label='hourly')
    axarr[0].set_ylabel('Precipitation [mm/10min]',fontsize=20)
    axarr[0].set_title(treenetstation+' WSL',fontsize=25)
    axarr[0].legend(loc=2)
    ymin1,ymax1 = axarr[0].get_ylim()

    axarr[1].bar(precip['date_meteoswiss'][interp_ind_0:interp_ind_1+1],\
             precip['precip_interpolated_10minres'][interp_ind_0:interp_ind_1+1],color='blue',width=0.005,align='center',label='10min')
    axarr[1].scatter(precip['date_meteoswiss'][interp_ind_0:interp_ind_1+1][meteoswiss_nans],\
             np.zeros((len(meteoswiss_nans))),color='red',s=20)
    axarr[1].plot(precip['date_hourly_meteoswiss'][interp_ind_2:interp_ind_3+1],\
             precip['precip_interpolated_hourly'][interp_ind_2:interp_ind_3+1],ls='--',color='green',label='hourly')
    axarr[1].set_ylabel('Precipitation [mm/10min]',fontsize=20)
    axarr[1].set_title(treenetstation+' Interpolation',fontsize=25)
    axarr[1].legend(loc=2)
    ymin2,ymax2 = axarr[1].get_ylim()

    axarr[2].bar(precip['date_combiprecip'][combi_ind_0:combi_ind_1+1],\
             precip['combiprecip'][combi_ind_0:combi_ind_1+1],color='blue',width=0.02,align='center',label='CombiPrecip',zorder=100)
    axarr[2].scatter(precip['date_combiprecip'][combi_ind_0:combi_ind_1+1][combiprecip_nans],\
             np.zeros((len(combiprecip_nans))),color='red',s=20)
    axarr[2].set_ylabel('Precipitation [mm/h]',fontsize=20)
    axarr[2].set_title('Combi Precip',fontsize=25)
    ymin3,ymax3 = axarr[2].get_ylim()

    axarr[0].set_ylim([0,np.max([ymax1,ymax2,ymax3])])
    axarr[1].set_ylim([0,np.max([ymax1,ymax2,ymax3])])
    axarr[2].set_ylim([0,np.max([ymax1,ymax2,ymax3])])

    plt.suptitle(precip['date_hourly_UTC'][wsl_ind_2].strftime('%b %Y'),fontsize=40)

    # Format x-axis and grid
    datelist = make_datelist(precip['date_combiprecip'][combi_ind_0],precip['date_combiprecip'][combi_ind_1],1)
    for axx in axarr:
        axx.xaxis_date()
        axx.xaxis.set_ticks(datelist[0::24*7],minor=False)
        axx.xaxis.set_major_formatter(dt.DateFormatter('%d-%m-%Y %H'))
        axx.xaxis.grid(True, which='major', color='0.5',alpha=0.7, linestyle='--',lw=1.5)

    saveas='\precip_comparison_'+treenetstation+'_'
    #plt.savefig(figpath+saveas+nowdate.strftime('%Y_%m')+'.png',bbox_inches='tight',dpi=400)
    plt.savefig(figpath+saveas+nowdate.strftime('%Y_%m')+'.pdf',bbox_inches='tight',format='pdf')
    #plt.savefig(figpath+saveas+nowdate.strftime('%Y_%m')+'.eps',bbox_inches='tight',format='eps')
    
    # Optimize performance
    #del wsl_ind_0,wsl_ind_1,wsl_ind_2,interp_ind_0,interp_ind_1,interp_ind_2,combi_ind_0,combi_ind_1,\
    #    wsl_nans,meteoswiss_nans,combiprecip_nans,axarr,f,ymin1,ymin2,ymin3,ymax1,ymax2,ymax3,saveas,datelist
    plt.close('all')
    gc.collect()





