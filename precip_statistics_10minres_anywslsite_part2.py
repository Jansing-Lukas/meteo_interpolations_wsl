# -*- coding: utf-8 -*-
"""
Created on Mon May 27 15:13:57 2019

@author: lukas jansing

script to validate interpolation of stationdata
10min resolution of stationdata
Quantitative validation by comparing interpolation, LWF data and CombiPrecip:
    - Daily and weekly sums
    - Seasonal sums
    - Yearly sums
    - Histograms
    - Correlation between interpolation and LWF, interpolation sum and CombiPrecip sum, CombiPrecip sum and LWF sum

"""

#--------------------------------------------------------- 
# Import modules
#---------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import datetime
from datetime import timedelta
import matplotlib
import math
from dateutil.relativedelta import *

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
#                       Novaggio, Visp, Vordemwald, Schaenis
# TreeNet site of Roman: Pfynwald
# TreeNet sites of IAP: Muri, Riehen, Grosswangen

treenetstation = 'LWF-Neunkirch1'

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
# Identify seasons to start and end
#---------------------------------------------------------  
# Identify next start of a season
# Part of winter belongs to last year
year_start = earliest_date.year
winter_start = datetime.datetime(year_start-1,12,1,0)
winter_end   = datetime.datetime(year_start,2,28,23)+relativedelta(day=31,hour=23)
spring_start = datetime.datetime(year_start,3,1,0)
spring_end   = datetime.datetime(year_start,5,31,23)
summer_start = datetime.datetime(year_start,6,1,0)
summer_end   = datetime.datetime(year_start,8,31,23)
autumn_start = datetime.datetime(year_start,9,1,0)
autumn_end   = datetime.datetime(year_start,11,30,23)

if earliest_date <= winter_start:
    earliest_date = winter_start
elif earliest_date <= spring_start:
    earliest_date = spring_start
elif earliest_date <= summer_start:
    earliest_date = summer_start
elif earliest_date <= autumn_start:
    earliest_date = autumn_start
elif earliest_date <= autumn_end:
    earliest_date = autumn_end+timedelta(hours=1) # in this case winter of next year is start    
else:
    print('strange')
    
# Identify last season
year_end = latest_date.year
winter_start = datetime.datetime(year_end-1,12,1,0)
winter_end   = datetime.datetime(year_end,2,28,23)+relativedelta(day=31,hour=23)
spring_start = datetime.datetime(year_end,3,1,0)
spring_end   = datetime.datetime(year_end,5,31,23)
summer_start = datetime.datetime(year_end,6,1,0)
summer_end   = datetime.datetime(year_end,8,31,23)
autumn_start = datetime.datetime(year_end,9,1,0)
autumn_end   = datetime.datetime(year_end,11,30,23)

if latest_date >= autumn_end:
    latest_date = autumn_end
elif latest_date >= summer_end:
    latest_date = summer_end
elif latest_date >= spring_end:
    latest_date = spring_end
elif latest_date >= winter_end:
    latest_date = winter_end
elif latest_date >= winter_start:
    latest_date = winter_start-timedelta(hours=1) # in this case winter of year before is end    
else:
    print('strange')

# Create seasonal datelist to loop over
datelist_seasons = []
nowdate = earliest_date
while nowdate <= latest_date:
    datelist_seasons.append(nowdate)
    if nowdate.month+3 <= 12:
        nowdate = nowdate+relativedelta(month=nowdate.month+3,day=1,hour=0)
    else:
        nowdate = nowdate+relativedelta(year=nowdate.year+1,month=3,day=1,hour=0)

#--------------------------------------------------------- 
# Loop over datelist to create seasonal summed plots
#---------------------------------------------------------   
matplotlib.rcParams.update({'font.size': 20})
plt.close('all')

seasonlabels = ['winter','spring','summer','autumn']
    
for nowdate in datelist_seasons:
    
    print(nowdate.strftime('%b %Y'))
    print('plotting')
    
    # start new figure in beginning or in winter
    if nowdate.month == 12 or nowdate == earliest_date:
        fig = plt.figure()
        fig.set_size_inches(20, 20)
        i = 1   # reset subplot index when starting new figure
    
    #Get latest instant of the season
    if nowdate.month != 12:
        season_end = nowdate+relativedelta(month=nowdate.month+2,day=31,hour=23)
    else:
        season_end = nowdate+relativedelta(year=nowdate.year+1,month=3,day=1,hour=0)
 
    # Find indices of the season in the timeseries
    wsl_ind_0    = np.where(precip['date_10minres_UTC'] == nowdate)[0][0]
    interp_ind_0 = np.where(precip['date_meteoswiss'] == nowdate)[0][0]
    combi_ind_0  = np.where(precip['date_combiprecip'] == nowdate)[0][0]
    wsl_ind_1    = np.where(precip['date_10minres_UTC'] == season_end)[0][0]
    interp_ind_1 = np.where(precip['date_meteoswiss'] == season_end)[0][0]
    combi_ind_1  = np.where(precip['date_combiprecip'] == season_end)[0][0]

    # Get indices of wsl nans and meteoswiss nans
    wsl_nans         = np.argwhere(np.isnan(precip['precip_10minres'][wsl_ind_0:wsl_ind_1+1]))[:,0]
    meteoswiss_nans  = np.argwhere(np.isnan(precip['precip_interpolated_10minres'][interp_ind_0:interp_ind_1+1]))[:,0]
    combiprecip_nans = np.argwhere(np.isnan(precip['combiprecip'][combi_ind_0:combi_ind_1+1]))[:,0]
        
    # Calculate seasonal sums
    # Only valid if more than 5% are no-nan values
    if (wsl_nans.shape[0]/precip['precip_10minres'][wsl_ind_0:wsl_ind_1+1].shape[0]) < 0.05:
        precip['wsl_sum'] = np.nansum(precip['precip_10minres'][wsl_ind_0:wsl_ind_1+1])
    else:
        precip['wsl_sum'] = np.nan
    if (meteoswiss_nans.shape[0]/precip['precip_interpolated_10minres'][interp_ind_0:interp_ind_1+1].shape[0]) < 0.05:
        precip['interpolated_sum'] = np.nansum(precip['precip_interpolated_10minres'][interp_ind_0:interp_ind_1+1])
    else:
        precip['interpolated_sum'] = np.nan
    if (combiprecip_nans.shape[0]/precip['combiprecip'][combi_ind_0:combi_ind_1+1].shape[0]) < 0.05:
        precip['combiprecip_sum'] = np.nansum(precip['combiprecip'][combi_ind_0:combi_ind_1+1])
    else:
        precip['combiprecip_sum'] = np.nan

    # Set subplot index for first subplot arrangement
    # Depends upon in which season the start is
    if nowdate == earliest_date:
        if  earliest_date.month == 12:
            i = 1
        elif earliest_date.month == 3:
            i = 2
        elif earliest_date.month == 6:
            i = 3
        elif earliest_date.month == 9:
            i = 4
    
    # Plot bar charts
    fig.add_subplot(2,2,i)
    x = np.arange(3)+0.4
    plt.bar(x,[precip['wsl_sum'],precip['interpolated_sum'],precip['combiprecip_sum']],align='center',width=0.4)
    plt.xlim([0,2.4+0.4])
    labels = ['LWF','Interpolated','Combiprecip']
    plt.xticks(x,labels)
    plt.title(seasonlabels[i-1])
    plt.ylabel('summed precipitation [mm]')
    
    # Add table with values below
    plt.table(cellText=[['%.1f' % precip['wsl_sum'],'%.1f' % precip['interpolated_sum'],'%.1f' % precip['combiprecip_sum']]],\
              bbox = [0.0,-0.12, 1.0, 0.05],cellLoc='center',rowLoc='center',fontsize=20)
    
    # Save figure if season is full   
    if i == 4:
        plt.suptitle(nowdate.strftime('%Y'),fontsize=40)
        saveas='\precip_statistics_'+treenetstation+'_seasonalsum_'
        plt.savefig(figpath+saveas+nowdate.strftime('%Y')+'.png',bbox_inches='tight')
    
    # Increase subplot index    
    i += 1
    
#--------------------------------------------------------- 
# Create plots of yearly sums
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

# Create yearly datelist to loop over
earliest_date = earliest_date - relativedelta(year=earliest_date.year+1,month=1,day=1,hour=0)
# Only use last year if it is fully available until the end
if latest_date.month == 12 and latest_date.day == 31 and latest_date.hour == 23:
    latest_date   = latest_date - relativedelta(year=latest_date.year,month=1,day=1,hour=0)
else:
    latest_date   = latest_date - relativedelta(year=latest_date.year-1,month=1,day=1,hour=0)
datelist_years = []
nowdate = earliest_date
while nowdate <= latest_date:
    datelist_years.append(nowdate)
    nowdate = nowdate+relativedelta(year=nowdate.year+1,month=1,day=1,hour=0)
    
# Loop over datelist to create yearly summed plots
matplotlib.rcParams.update({'font.size': 20})
plt.close('all')
fig = plt.figure()
fig.set_size_inches(26, 26)
i = 1   # subplot index
        
for nowdate in datelist_years:
    
    print(nowdate.strftime('%Y'))
    print('plotting')
    
    year_end = nowdate+relativedelta(year=nowdate.year,month=12,day=31,hour=23)
    
    # Find indices of the year in the timeseries
    wsl_ind_0    = np.where(precip['date_10minres_UTC'] == nowdate)[0][0]
    interp_ind_0 = np.where(precip['date_meteoswiss'] == nowdate)[0][0]
    combi_ind_0  = np.where(precip['date_combiprecip'] == nowdate)[0][0]
    wsl_ind_1    = np.where(precip['date_10minres_UTC'] == year_end)[0][0]
    interp_ind_1 = np.where(precip['date_meteoswiss'] == year_end)[0][0]
    combi_ind_1  = np.where(precip['date_combiprecip'] == year_end)[0][0]

    # Get indices of wsl nans and meteoswiss nans
    wsl_nans         = np.argwhere(np.isnan(precip['precip_10minres'][wsl_ind_0:wsl_ind_1+1]))[:,0]
    meteoswiss_nans  = np.argwhere(np.isnan(precip['precip_interpolated_10minres'][interp_ind_0:interp_ind_1+1]))[:,0]
    combiprecip_nans = np.argwhere(np.isnan(precip['combiprecip'][combi_ind_0:combi_ind_1+1]))[:,0]
    
    # Calculate yearly sums
    # Only valid if more than 5% are no-nan values
    if (wsl_nans.shape[0]/precip['precip_10minres'][wsl_ind_0:wsl_ind_1+1].shape[0]) < 0.05:
        precip['wsl_sum'] = np.nansum(precip['precip_10minres'][wsl_ind_0:wsl_ind_1+1])
    else:
        precip['wsl_sum'] = np.nan
    if (meteoswiss_nans.shape[0]/precip['precip_interpolated_10minres'][interp_ind_0:interp_ind_1+1].shape[0]) < 0.05:
        precip['interpolated_sum'] = np.nansum(precip['precip_interpolated_10minres'][interp_ind_0:interp_ind_1+1])
    else:
        precip['interpolated_sum'] = np.nan
    if (combiprecip_nans.shape[0]/precip['combiprecip'][combi_ind_0:combi_ind_1+1].shape[0]) < 0.05:
        precip['combiprecip_sum'] = np.nansum(precip['combiprecip'][combi_ind_0:combi_ind_1+1])
    else:
        precip['combiprecip_sum'] = np.nan

    # Plot bar charts
    fig.add_subplot(math.ceil(np.sqrt(len(datelist_years))),math.ceil(np.sqrt(len(datelist_years))),i)
    x = np.arange(3)+0.4
    plt.bar(x,[precip['wsl_sum'],precip['interpolated_sum'],precip['combiprecip_sum']],align='center',width=0.4)
    plt.xlim([0,2.4+0.4])
    labels = ['LWF','Interpolated','Combiprecip']
    plt.xticks(x,labels)
    plt.title(nowdate.strftime('%Y'))
    plt.ylabel('summed precipitation [mm]')
    
    # Add table with values below
    plt.table(cellText=[['%.1f' % precip['wsl_sum'],'%.1f' % precip['interpolated_sum'],'%.1f' % precip['combiprecip_sum']]],\
              bbox = [0.0,-0.12, 1.0, 0.05],cellLoc='center',rowLoc='center',fontsize=20)
    
    # Increase subplot index    
    i += 1
    
plt.suptitle('Yearly precipitation sums at '+treenetstation,fontsize=40)
saveas='\precip_statistics_'+treenetstation+'_yearlysum_'
plt.savefig(figpath+saveas+earliest_date.strftime('%Y')+'-'+latest_date.strftime('%Y')+'.png',\
            bbox_inches='tight')

#--------------------------------------------------------- 
# Create histograms of the data
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
if (latest_date.day != 31 or 30) or latest_date.hour != 23 or latest_date.minute != 0:
    latest_date = latest_date - relativedelta(months=1,day=31,hour=23,minute=0)    

print('plotting')
plt.close('all')
fig = plt.figure()
fig.set_size_inches(27, 12)

# Find indices of the timespan in the timeseries
wsl_ind_0    = np.where(precip['date_10minres_UTC'] == earliest_date)[0][0]
interp_ind_0 = np.where(precip['date_meteoswiss'] == earliest_date)[0][0]
combi_ind_0  = np.where(precip['date_combiprecip'] == earliest_date)[0][0]
wsl_ind_1    = np.where(precip['date_10minres_UTC'] == latest_date)[0][0]
interp_ind_1 = np.where(precip['date_meteoswiss'] == latest_date)[0][0]
combi_ind_1  = np.where(precip['date_combiprecip'] == latest_date)[0][0]

# Get indices of wsl nans and meteoswiss nans
wsl_nans         = np.argwhere(np.isnan(precip['precip_10minres'][wsl_ind_0:wsl_ind_1+1]))[:,0]
meteoswiss_nans  = np.argwhere(np.isnan(precip['precip_interpolated_10minres'][interp_ind_0:interp_ind_1+1]))[:,0]
combiprecip_nans = np.argwhere(np.isnan(precip['combiprecip'][combi_ind_0:combi_ind_1+1]))[:,0]

# Get maxima of timeseries to create bin array
wsl_max         = np.nanmax(precip['precip_10minres'][wsl_ind_0:wsl_ind_1+1])
interp_max      = np.nanmax(precip['precip_interpolated_10minres'][interp_ind_0:interp_ind_1+1])
combiprecip_max = np.nanmax(precip['combiprecip'][combi_ind_0:combi_ind_1+1])
bin_max         = math.ceil(np.max([wsl_max,interp_max,combiprecip_max]))
bin_array       = np.arange(0,bin_max+1,1)

# Plot histograms
# Only if at least 50% no-nan values
# Note: Only take the non-zero = actual precipitation events!
fig.add_subplot(1,3,1)
plt.hist(precip['precip_10minres'][wsl_ind_0:wsl_ind_1+1][precip['precip_10minres'][wsl_ind_0:wsl_ind_1+1] > 0],\
         bins=bin_array,alpha=0.7,edgecolor='black',density=True)
plt.yscale('log')
plt.xlabel('precipitation [mm/10min]')
plt.ylabel('normalized frequency')
plt.title('Histogram of WSL data')
axes1 = plt.gca()   # define axis anyway
ymin1,ymax1 = axes1.get_ylim()  

    
fig.add_subplot(1,3,2)
plt.hist(precip['precip_interpolated_10minres'][interp_ind_0:interp_ind_1+1][precip['precip_interpolated_10minres'][interp_ind_0:interp_ind_1+1] > 0],\
         bins=bin_array,alpha=0.7,edgecolor='black',density=True)
plt.yscale('log')
plt.xlabel('precipitation [mm/10min]')
plt.ylabel('normalized frequency')
plt.title('Histogram of interpolated data')
axes2 = plt.gca()
ymin2,ymax2 = axes2.get_ylim()


fig.add_subplot(1,3,3)
plt.hist(precip['combiprecip'][combi_ind_0:combi_ind_1+1][precip['combiprecip'][combi_ind_0:combi_ind_1+1] > 0],\
         bins=bin_array,alpha=0.7,edgecolor='black',density=True)
plt.yscale('log')
plt.xlabel('precipitation [mm/h]')
plt.ylabel('normalized frequency')
plt.title('Histogram of CombiPrecip data')
axes3 = plt.gca()
ymin3,ymax3 = axes3.get_ylim()

    
axes1.set_ylim([np.min([ymin1,ymin2,ymin3]),np.max([ymax1,ymax2,ymax3])])
axes2.set_ylim([np.min([ymin1,ymin2,ymin3]),np.max([ymax1,ymax2,ymax3])])
axes3.set_ylim([np.min([ymin1,ymin2,ymin3]),np.max([ymax1,ymax2,ymax3])])
for ax in [axes1,axes2,axes3]:
    ax.set_xticks(bin_array[::3])

    
saveas='\precip_statistics_'+treenetstation+'_histograms_'
plt.savefig(figpath+saveas+earliest_date.strftime('%Y')+'-'+latest_date.strftime('%Y')+'.png',\
            bbox_inches='tight')

#--------------------------------------------------------- 
# Write out temporal correlation of the timeseries
#---------------------------------------------------------
# Redefine earliest and latest date based on hourly timeseries
# Find latest starting date
if (precip['date_combiprecip'][0] > precip['date_meteoswiss'][0]) and \
   (precip['date_combiprecip'][0] > precip['date_10minres_UTC'][0]):
       earliest_date = precip['date_combiprecip'][0]
elif (precip['date_hourly_meteoswiss'][0] > precip['date_combiprecip'][0]) and \
     (precip['date_hourly_meteoswiss'][0] > precip['date_hourly_UTC'][0]):
        earliest_date = precip['date_hourly_meteoswiss'][0]
else:
    earliest_date = precip['date_hourly_UTC'][0]

# Add up to next month
year = earliest_date.year
if earliest_date.day != 1 or earliest_date.hour != 0:
    month = earliest_date.month+1
    earliest_date = datetime.datetime(year,month,1,0)

# Find earliest ending date    
if (precip['date_combiprecip'][-1] < precip['date_meteoswiss'][-1]) and \
   (precip['date_combiprecip'][-1] < precip['date_hourly_UTC'][-1]):
       latest_date = precip['date_combiprecip'][-1]
elif (precip['date_hourly_meteoswiss'][-1] < precip['date_combiprecip'][-1]) and \
     (precip['date_hourly_meteoswiss'][-1] < precip['date_hourly_UTC'][-1]):
        latest_date = precip['date_hourly_meteoswiss'][-1]
else:
    latest_date = precip['date_hourly_UTC'][-1]

lwf_precip_10minres          = precip['precip_10minres'][np.where(precip['date_10minres_UTC'] == earliest_date)[0][0]\
                                      :np.where(precip['date_10minres_UTC'] == latest_date)[0][0]+1]
interpolated_precip_10minres = precip['precip_interpolated_10minres'][np.where(precip['date_meteoswiss'] == earliest_date)[0][0]\
                                      :np.where(precip['date_meteoswiss'] == latest_date)[0][0]+1]
lwf_precip_hourly            = precip['precip_hourly'][np.where(precip['date_hourly_UTC'] == earliest_date)[0][0]\
                                      :np.where(precip['date_hourly_UTC'] == latest_date)[0][0]+1]
interpolated_precip_hourly   = precip['precip_interpolated_hourly'][np.where(precip['date_hourly_meteoswiss'] == earliest_date)[0][0]\
                                      :np.where(precip['date_hourly_meteoswiss'] == latest_date)[0][0]+1]
combi_precip                 = precip['combiprecip'][np.where(precip['date_combiprecip'] == earliest_date)[0][0]\
                                      :np.where(precip['date_combiprecip'] == latest_date)[0][0]+1]

bad1 = ~np.logical_or(np.isnan(lwf_precip_10minres), np.isnan(interpolated_precip_10minres))
bad2 = ~np.logical_or(np.isnan(lwf_precip_hourly), np.isnan(interpolated_precip_hourly))
bad3 = ~np.logical_or(np.isnan(lwf_precip_hourly), np.isnan(combi_precip))
bad4 = ~np.logical_or(np.isnan(interpolated_precip_hourly), np.isnan(combi_precip))

lwf_compressed_10minres          = np.compress(bad1,lwf_precip_10minres)
interpolated_compressed_10minres = np.compress(bad1,interpolated_precip_10minres)

lwf_compressed_1 = np.compress(bad2,lwf_precip_hourly)
lwf_compressed_2 = np.compress(bad3,lwf_precip_hourly)
interpolated_compressed_1 = np.compress(bad2,interpolated_precip_hourly)
interpolated_compressed_3 = np.compress(bad4,interpolated_precip_hourly)
combi_compressed_2 = np.compress(bad3,combi_precip)
combi_compressed_3 = np.compress(bad4,combi_precip)

ccoef1 = np.corrcoef(lwf_compressed_10minres,interpolated_compressed_10minres)[0,1]
ccoef2 = np.corrcoef(lwf_compressed_1,interpolated_compressed_1)[0,1]
ccoef3 = np.corrcoef(lwf_compressed_2,combi_compressed_2)[0,1]
ccoef4 = np.corrcoef(interpolated_compressed_3,combi_compressed_3)[0,1]

print('--------------------------------------------------------------')
print('Correlation coefficients:')
print('LWF and interpolation (10minres) '+str(ccoef1))
print('LWF and interpolation (hourly)   '+str(ccoef2))
print('LWF and CombiPrecip              '+str(ccoef3))
print('Interpolation and CombiPrecip    '+str(ccoef4))
print('--------------------------------------------------------------')    
