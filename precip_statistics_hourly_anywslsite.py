# -*- coding: utf-8 -*-
"""
Created on Wed May 15 12:56:12 2019

@author: lukas jansing

script in order to create statistical plots of the different precipitation data
plotting of: winter sums, summer sums, yearly sums, maximas, histograms

"""

#--------------------------------------------------------- 
# Import modules and functions
#---------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import datetime
from datetime import timedelta
import math
import matplotlib
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

treenetstation = 'Schaenis'

station_id = treenetstation_id[treenetstationnames.index(treenetstation)]

#process_schänis = 'yes'
process_schänis = 'no'

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
                                       path=treenetprecippath,process_schänis=process_schänis)

#--------------------------------------------------------- 
# Load interpolated precipitation data
#---------------------------------------------------------
print('import interpolated precipitation data')
interpolated_precipitation_data = np.load(stationprecippath+'\interpolated_precipitation_version3.npy')
#interpolated_precipitation_data = np.load(stationprecippath+'\interpolated_precipitation_newmethod.npy')

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

# Add up to next month
year = earliest_date.year
if earliest_date.day != 1 or earliest_date.hour != 0:
    month = earliest_date.month+1
    earliest_date = datetime.datetime(year,month,1,0)

# Find earliest ending date    
if (precip['date_combiprecip'][-1] < precip['date_meteoswiss'][-1]) and \
   (precip['date_combiprecip'][-1] < precip['date_hourly_UTC'][-1]):
       latest_date = precip['date_combiprecip'][-1]
elif (precip['date_meteoswiss'][-1] < precip['date_combiprecip'][-1]) and \
     (precip['date_meteoswiss'][-1] < precip['date_hourly_UTC'][-1]):
        latest_date = precip['date_meteoswiss'][-1]
else:
    latest_date = precip['date_hourly_UTC'][-1]

# Subtract to previous month
year = latest_date.year
if (latest_date.day != 31 or 30) or latest_date.hour != 23:
    latest_date = latest_date - relativedelta(months=1,day=31,hour=23)
    
#--------------------------------------------------------- 
# Identify seasons to start and end
#---------------------------------------------------------  
# Identify next start of a season
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

# Define lists for the values
precip['wsl_winter_sum']          = []
precip['interpolated_winter_sum'] = []
precip['combiprecip_winter_sum']  = []
precip['wsl_spring_sum']          = []
precip['interpolated_spring_sum'] = []
precip['combiprecip_spring_sum']  = []
precip['wsl_summer_sum']          = []
precip['interpolated_summer_sum'] = []
precip['combiprecip_summer_sum']  = []
precip['wsl_autumn_sum']          = []
precip['interpolated_autumn_sum'] = []
precip['combiprecip_autumn_sum']  = []
    
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
    wsl_ind_0    = np.where(precip['date_hourly_UTC'] == nowdate)[0][0]
    interp_ind_0 = np.where(precip['date_meteoswiss'] == nowdate)[0][0]
    combi_ind_0  = np.where(precip['date_combiprecip'] == nowdate)[0][0]
    wsl_ind_1    = np.where(precip['date_hourly_UTC'] == season_end)[0][0]
    interp_ind_1 = np.where(precip['date_meteoswiss'] == season_end)[0][0]
    combi_ind_1  = np.where(precip['date_combiprecip'] == season_end)[0][0]

    # Get indices of wsl nans and meteoswiss nans
    wsl_nans         = np.argwhere(np.isnan(precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1]))[:,0]
    meteoswiss_nans  = np.argwhere(np.isnan(precip['precip_interpolated'][interp_ind_0:interp_ind_1+1]))[:,0]
    combiprecip_nans = np.argwhere(np.isnan(precip['combiprecip'][combi_ind_0:combi_ind_1+1]))[:,0]
        
    # Calculate seasonal sums
    # Only valid if more than 5% are no-nan values
    if (wsl_nans.shape[0]/precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1].shape[0]) < 0.05:
        precip['wsl_sum'] = np.nansum(precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1])
    else:
        precip['wsl_sum'] = np.nan
    if (meteoswiss_nans.shape[0]/precip['precip_interpolated'][interp_ind_0:interp_ind_1+1].shape[0]) < 0.05:
        precip['interpolated_sum'] = np.nansum(precip['precip_interpolated'][interp_ind_0:interp_ind_1+1])
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
    
    # Append values
    if nowdate.year >= 2013 and nowdate.month == 12:
        precip['wsl_winter_sum'].append(precip['wsl_sum'])
        precip['interpolated_winter_sum'].append(precip['interpolated_sum'])
        precip['combiprecip_winter_sum'].append(precip['combiprecip_sum'])
    elif nowdate.year >= 2014 and nowdate.month == 3:
        precip['wsl_spring_sum'].append(precip['wsl_sum'])
        precip['interpolated_spring_sum'].append(precip['interpolated_sum'])
        precip['combiprecip_spring_sum'].append(precip['combiprecip_sum'])
    elif nowdate.year >= 2014 and nowdate.month == 6:
        precip['wsl_summer_sum'].append(precip['wsl_sum'])
        precip['interpolated_summer_sum'].append(precip['interpolated_sum'])
        precip['combiprecip_summer_sum'].append(precip['combiprecip_sum'])
    elif nowdate.year >= 2014 and nowdate.month == 9:
        precip['wsl_autumn_sum'].append(precip['wsl_sum'])
        precip['interpolated_autumn_sum'].append(precip['interpolated_sum'])
        precip['combiprecip_autumn_sum'].append(precip['combiprecip_sum'])
    
    # Save figure if season is full   
    if i == 4:
        plt.suptitle(nowdate.strftime('%Y'),fontsize=40)
        saveas='\precip_statistics_'+treenetstation+'_seasonalsum_'
        plt.savefig(figpath+saveas+nowdate.strftime('%Y')+'.png',bbox_inches='tight')
    
    # Increase subplot index    
    i += 1
    
# Print out mean values
print('--------------------------------------------------------------')
print('Mean seasonal sums:')
print('LWF winter:   '+str(np.nanmean(precip['wsl_winter_sum'])))
print('LWF spring:   '+str(np.nanmean(precip['wsl_spring_sum'])))
print('LWF summer:   '+str(np.nanmean(precip['wsl_summer_sum'])))
print('LWF autumn:   '+str(np.nanmean(precip['wsl_autumn_sum'])))
print('----------------------')
print('Interpolated winter:   '+str(np.nanmean(precip['interpolated_winter_sum'])))
print('Interpolated spring:   '+str(np.nanmean(precip['interpolated_spring_sum'])))
print('Interpolated summer:   '+str(np.nanmean(precip['interpolated_summer_sum'])))
print('Interpolated autumn:   '+str(np.nanmean(precip['interpolated_autumn_sum'])))
print('----------------------')
print('CombiPrecip winter:   '+str(np.nanmean(precip['combiprecip_winter_sum'])))
print('CombiPrecip spring:   '+str(np.nanmean(precip['combiprecip_spring_sum'])))
print('CombiPrecip summer:   '+str(np.nanmean(precip['combiprecip_summer_sum'])))
print('CombiPrecip autumn:   '+str(np.nanmean(precip['combiprecip_autumn_sum'])))
print('--------------------------------------------------------------')  
    
#--------------------------------------------------------- 
# Create plots of yearly sums
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

# Add up to next month
year = earliest_date.year
if earliest_date.day != 1 or earliest_date.hour != 0:
    month = earliest_date.month+1
    earliest_date = datetime.datetime(year,month,1,0)

# Find earliest ending date    
if (precip['date_combiprecip'][-1] < precip['date_meteoswiss'][-1]) and \
   (precip['date_combiprecip'][-1] < precip['date_hourly_UTC'][-1]):
       latest_date = precip['date_combiprecip'][-1]
elif (precip['date_meteoswiss'][-1] < precip['date_combiprecip'][-1]) and \
     (precip['date_meteoswiss'][-1] < precip['date_hourly_UTC'][-1]):
        latest_date = precip['date_meteoswiss'][-1]
else:
    latest_date = precip['date_hourly_UTC'][-1]

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

# Define lists
precip['wsl_yearly_sum']          = []
precip['interpolated_yearly_sum'] = []
precip['combiprecip_yearly_sum']  = []
        
for nowdate in datelist_years:
    
    print(nowdate.strftime('%Y'))
    print('plotting')
    
    year_end = nowdate+relativedelta(year=nowdate.year,month=12,day=31,hour=23)
    
    # Find indices of the year in the timeseries
    wsl_ind_0    = np.where(precip['date_hourly_UTC'] == nowdate)[0][0]
    interp_ind_0 = np.where(precip['date_meteoswiss'] == nowdate)[0][0]
    combi_ind_0  = np.where(precip['date_combiprecip'] == nowdate)[0][0]
    wsl_ind_1    = np.where(precip['date_hourly_UTC'] == year_end)[0][0]
    interp_ind_1 = np.where(precip['date_meteoswiss'] == year_end)[0][0]
    combi_ind_1  = np.where(precip['date_combiprecip'] == year_end)[0][0]

    # Get indices of wsl nans and meteoswiss nans
    wsl_nans         = np.argwhere(np.isnan(precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1]))[:,0]
    meteoswiss_nans  = np.argwhere(np.isnan(precip['precip_interpolated'][interp_ind_0:interp_ind_1+1]))[:,0]
    combiprecip_nans = np.argwhere(np.isnan(precip['combiprecip'][combi_ind_0:combi_ind_1+1]))[:,0]
    
    # Calculate yearly sums
    # Only valid if more than 5% are no-nan values
    if (wsl_nans.shape[0]/precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1].shape[0]) < 0.05:
        precip['wsl_sum'] = np.nansum(precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1])
    else:
        precip['wsl_sum'] = np.nan
    if (meteoswiss_nans.shape[0]/precip['precip_interpolated'][interp_ind_0:interp_ind_1+1].shape[0]) < 0.05:
        precip['interpolated_sum'] = np.nansum(precip['precip_interpolated'][interp_ind_0:interp_ind_1+1])
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
    
    # Append data
    if nowdate.year >= 2014:
        precip['wsl_yearly_sum'].append(precip['wsl_sum'])
        precip['interpolated_yearly_sum'].append(precip['interpolated_sum'])
        precip['combiprecip_yearly_sum'].append(precip['combiprecip_sum'])
    
    # Increase subplot index    
    i += 1
    
plt.suptitle('Yearly precipitation sums at '+treenetstation,fontsize=40)
saveas='\precip_statistics_'+treenetstation+'_yearlysum_'
plt.savefig(figpath+saveas+earliest_date.strftime('%Y')+'-'+latest_date.strftime('%Y')+'.png',\
            bbox_inches='tight')

#--------------------------------------------------------- 
# Create plots of yearly maxima
#---------------------------------------------------------
# Loop over datelist to create yearly maxima plots
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
    wsl_ind_0    = np.where(precip['date_hourly_UTC'] == nowdate)[0][0]
    interp_ind_0 = np.where(precip['date_meteoswiss'] == nowdate)[0][0]
    combi_ind_0  = np.where(precip['date_combiprecip'] == nowdate)[0][0]
    wsl_ind_1    = np.where(precip['date_hourly_UTC'] == year_end)[0][0]
    interp_ind_1 = np.where(precip['date_meteoswiss'] == year_end)[0][0]
    combi_ind_1  = np.where(precip['date_combiprecip'] == year_end)[0][0]

    # Get indices of wsl nans and meteoswiss nans
    wsl_nans         = np.argwhere(np.isnan(precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1]))[:,0]
    meteoswiss_nans  = np.argwhere(np.isnan(precip['precip_interpolated'][interp_ind_0:interp_ind_1+1]))[:,0]
    combiprecip_nans = np.argwhere(np.isnan(precip['combiprecip'][combi_ind_0:combi_ind_1+1]))[:,0]
    
    # Calculate yearly maxima
    precip['wsl_max'] = np.nanmax(precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1])
    precip['interpolated_max'] = np.nanmax(precip['precip_interpolated'][interp_ind_0:interp_ind_1+1])
    precip['combiprecip_max'] = np.nanmax(precip['combiprecip'][combi_ind_0:combi_ind_1+1])
        
    # Plot bar charts
    fig.add_subplot(math.ceil(np.sqrt(len(datelist_years))),math.ceil(np.sqrt(len(datelist_years))),i)
    x = np.arange(3)+0.4
    plt.bar(x,[precip['wsl_max'],precip['interpolated_max'],precip['combiprecip_max']],align='center',width=0.4)
    plt.xlim([0,2.4+0.4])
    labels = ['LWF','Interpolated','Combiprecip']
    plt.xticks(x,labels)
    plt.title(nowdate.strftime('%Y'))
    plt.ylabel('maximal hourly precipitation [mm]')
    
    # Add table with values below
    plt.table(cellText=[['%.1f' % precip['wsl_max'],'%.1f' % precip['interpolated_max'],'%.1f' % precip['combiprecip_max']]],\
              bbox = [0.0,-0.12, 1.0, 0.05],cellLoc='center',rowLoc='center',fontsize=20)
    
    # Increase subplot index    
    i += 1
    
plt.suptitle('Yearly precipitation maxima at '+treenetstation,fontsize=40)
saveas='\precip_statistics_'+treenetstation+'_yearlymax_'
plt.savefig(figpath+saveas+earliest_date.strftime('%Y')+'-'+latest_date.strftime('%Y')+'.png',\
            bbox_inches='tight')

#--------------------------------------------------------- 
# Create histograms of the data
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

# Add up to next month
year = earliest_date.year
if earliest_date.day != 1 or earliest_date.hour != 0:
    month = earliest_date.month+1
    earliest_date = datetime.datetime(year,month,1,0)

# Find earliest ending date    
if (precip['date_combiprecip'][-1] < precip['date_meteoswiss'][-1]) and \
   (precip['date_combiprecip'][-1] < precip['date_hourly_UTC'][-1]):
       latest_date = precip['date_combiprecip'][-1]
elif (precip['date_meteoswiss'][-1] < precip['date_combiprecip'][-1]) and \
     (precip['date_meteoswiss'][-1] < precip['date_hourly_UTC'][-1]):
        latest_date = precip['date_meteoswiss'][-1]
else:
    latest_date = precip['date_hourly_UTC'][-1]
    
matplotlib.rcParams.update({'font.size': 20})
plt.close('all')

print('plotting')
fig = plt.figure()
fig.set_size_inches(27, 12)

# Find indices of the timespan in the timeseries
wsl_ind_0    = np.where(precip['date_hourly_UTC'] == earliest_date)[0][0]
interp_ind_0 = np.where(precip['date_meteoswiss'] == earliest_date)[0][0]
combi_ind_0  = np.where(precip['date_combiprecip'] == earliest_date)[0][0]
wsl_ind_1    = np.where(precip['date_hourly_UTC'] == latest_date)[0][0]
interp_ind_1 = np.where(precip['date_meteoswiss'] == latest_date)[0][0]
combi_ind_1  = np.where(precip['date_combiprecip'] == latest_date)[0][0]

# Get indices of wsl nans and meteoswiss nans
wsl_nans         = np.argwhere(np.isnan(precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1]))[:,0]
meteoswiss_nans  = np.argwhere(np.isnan(precip['precip_interpolated'][interp_ind_0:interp_ind_1+1]))[:,0]
combiprecip_nans = np.argwhere(np.isnan(precip['combiprecip'][combi_ind_0:combi_ind_1+1]))[:,0]

# Get maxima of timeseries to create bin array
wsl_max         = np.nanmax(precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1])
interp_max      = np.nanmax(precip['precip_interpolated'][interp_ind_0:interp_ind_1+1])
combiprecip_max = np.nanmax(precip['combiprecip'][combi_ind_0:combi_ind_1+1])
bin_max         = math.ceil(np.max([wsl_max,interp_max,combiprecip_max]))
bin_array       = np.arange(0,bin_max+1,1)

# Plot histograms
# Only if at least 50% no-nan values
# Note: Only take the non-zero = actual precipitation events!
fig.add_subplot(1,3,1)
plt.hist(precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1][precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1] > 0],\
         bins=bin_array,alpha=0.7,edgecolor='black',density=True)
plt.yscale('log')
plt.xlabel('precipitation [mm/h]')
plt.ylabel('normalized frequency')
plt.title('Histogram of WSL data')
axes1 = plt.gca()   # define axis anyway
ymin1,ymax1 = axes1.get_ylim()  

    
fig.add_subplot(1,3,2)
plt.hist(precip['precip_interpolated'][interp_ind_0:interp_ind_1+1][precip['precip_interpolated'][interp_ind_0:interp_ind_1+1] > 0],\
         bins=bin_array,alpha=0.7,edgecolor='black',density=True)
plt.yscale('log')
plt.xlabel('precipitation [mm/h]')
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
lwf_precip          = precip['precip_hourly'][np.where(precip['date_hourly_UTC'] == earliest_date)[0][0]:np.where(precip['date_hourly_UTC'] == latest_date)[0][0]+1]
interpolated_precip = precip['precip_interpolated'][np.where(precip['date_meteoswiss'] == earliest_date)[0][0]:np.where(precip['date_meteoswiss'] == latest_date)[0][0]+1]
combi_precip        = precip['combiprecip'][np.where(precip['date_combiprecip'] == earliest_date)[0][0]:np.where(precip['date_combiprecip'] == latest_date)[0][0]+1]

bad1 = ~np.logical_or(np.isnan(lwf_precip), np.isnan(interpolated_precip))
bad2 = ~np.logical_or(np.isnan(lwf_precip), np.isnan(combi_precip))
bad3 = ~np.logical_or(np.isnan(interpolated_precip), np.isnan(combi_precip))

lwf_compressed_1 = np.compress(bad1,lwf_precip)
lwf_compressed_2 = np.compress(bad2,lwf_precip)
interpolated_compressed_1 = np.compress(bad1,interpolated_precip)
interpolated_compressed_3 = np.compress(bad3,interpolated_precip)
combi_compressed_2 = np.compress(bad2,combi_precip)
combi_compressed_3 = np.compress(bad3,combi_precip)

ccoef1 = np.corrcoef(lwf_compressed_1,interpolated_compressed_1)[0,1]
ccoef2 = np.corrcoef(lwf_compressed_2,combi_compressed_2)[0,1]
ccoef3 = np.corrcoef(interpolated_compressed_3,combi_compressed_3)[0,1]

print('--------------------------------------------------------------')
print('Correlation coefficients:')
print('LWF and interpolation           '+str(ccoef1))
print('LWF and CombiPrecip             '+str(ccoef2))
print('Interpolation and CombiPrecip   '+str(ccoef3))
print('--------------------------------------------------------------')

#--------------------------------------------------------- 
# Calculate growing period (April-October) sums
#---------------------------------------------------------
# First: Reset earliest and latest date
# Find latest starting date
if (precip['date_combiprecip'][0] > precip['date_meteoswiss'][0]) and \
   (precip['date_combiprecip'][0] > precip['date_hourly_UTC'][0]):
       earliest_date = precip['date_combiprecip'][0]
elif (precip['date_meteoswiss'][0] > precip['date_combiprecip'][0]) and \
     (precip['date_meteoswiss'][0] > precip['date_hourly_UTC'][0]):
        earliest_date = precip['date_meteoswiss'][0]
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
elif (precip['date_meteoswiss'][-1] < precip['date_combiprecip'][-1]) and \
     (precip['date_meteoswiss'][-1] < precip['date_hourly_UTC'][-1]):
        latest_date = precip['date_meteoswiss'][-1]
else:
    latest_date = precip['date_hourly_UTC'][-1]

# Subtract to previous month
year = latest_date.year
if (latest_date.day != 31 or 30) or latest_date.hour != 23:
    latest_date = latest_date - relativedelta(months=1,day=31,hour=23)
    
# Create growing season datelist
if earliest_date.month > 4:
    earliest_date = earliest_date + relativedelta(month=4,day=1,hour=0)
else:
    earliest_date = earliest_date + relativedelta(year=earliest_date.year+1,month=4,day=1,hour=0)
datelist_growingperiods = []
nowdate = earliest_date
while nowdate <= latest_date:
    datelist_growingperiods.append(nowdate)
    nowdate = nowdate+relativedelta(year=nowdate.year+1,month=4,day=1,hour=0)

# Define lists
precip['wsl_growingseason_sum']          = []
precip['interpolated_growingseason_sum'] = []
precip['combiprecip_growingseason_sum']  = []

matplotlib.rcParams.update({'font.size': 20})
plt.close('all')
fig = plt.figure()
fig.set_size_inches(26, 26)
i = 1   # subplot index
   
# Loop over growing season datelist in order to create plots
for nowdate in datelist_growingperiods:
    
    print(nowdate.strftime('%Y'))
    print('plotting')
    
    year_end = nowdate+relativedelta(year=nowdate.year,month=10,day=31,hour=23)
    
    # Find indices of the year in the timeseries
    wsl_ind_0    = np.where(precip['date_hourly_UTC'] == nowdate)[0][0]
    interp_ind_0 = np.where(precip['date_meteoswiss'] == nowdate)[0][0]
    combi_ind_0  = np.where(precip['date_combiprecip'] == nowdate)[0][0]
    wsl_ind_1    = np.where(precip['date_hourly_UTC'] == year_end)[0][0]
    interp_ind_1 = np.where(precip['date_meteoswiss'] == year_end)[0][0]
    combi_ind_1  = np.where(precip['date_combiprecip'] == year_end)[0][0]

    # Get indices of wsl nans and meteoswiss nans
    wsl_nans         = np.argwhere(np.isnan(precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1]))[:,0]
    meteoswiss_nans  = np.argwhere(np.isnan(precip['precip_interpolated'][interp_ind_0:interp_ind_1+1]))[:,0]
    combiprecip_nans = np.argwhere(np.isnan(precip['combiprecip'][combi_ind_0:combi_ind_1+1]))[:,0]
    
    # Calculate yearly sums
    # Only valid if more than 5% are no-nan values
    if (wsl_nans.shape[0]/precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1].shape[0]) < 0.05:
        precip['wsl_sum'] = np.nansum(precip['precip_hourly'][wsl_ind_0:wsl_ind_1+1])
    else:
        precip['wsl_sum'] = np.nan
    if (meteoswiss_nans.shape[0]/precip['precip_interpolated'][interp_ind_0:interp_ind_1+1].shape[0]) < 0.05:
        precip['interpolated_sum'] = np.nansum(precip['precip_interpolated'][interp_ind_0:interp_ind_1+1])
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
    
    # Append data
    if nowdate.year >= 2014:
        precip['wsl_growingseason_sum'].append(precip['wsl_sum'])
        precip['interpolated_growingseason_sum'].append(precip['interpolated_sum'])
        precip['combiprecip_growingseason_sum'].append(precip['combiprecip_sum'])
    
    # Increase subplot index    
    i += 1
    
plt.suptitle('Growing season precipitation sums at '+treenetstation,fontsize=40)
saveas='\precip_statistics_'+treenetstation+'_growingseasonsum_'
plt.savefig(figpath+saveas+earliest_date.strftime('%Y')+'-'+latest_date.strftime('%Y')+'.png',\
            bbox_inches='tight')

# Print out mean values
print('--------------------------------------------------------------')
print('Mean growing season sums:')
print('LWF:   '+str(np.nanmean(precip['wsl_growingseason_sum'])))
print('----------------------')
print('Interpolated:   '+str(np.nanmean(precip['interpolated_growingseason_sum'])))
print('----------------------')
print('CombiPrecip:   '+str(np.nanmean(precip['combiprecip_growingseason_sum'])))
print('--------------------------------------------------------------')  


