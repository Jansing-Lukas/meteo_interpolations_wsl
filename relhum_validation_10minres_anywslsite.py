# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 12:40:03 2019

@author: lukas jansing

Script in order to validate the interpolated relative humidity data

Compare 10minres data to 10minres data of LWF

Compare the three different interpolations:
    - horizontally only
    - empirical local gradient

"""

#--------------------------------------------------------- 
# Import modules
#---------------------------------------------------------
import numpy as np
import datetime
import matplotlib.pyplot as plt
from dateutil.relativedelta import *
import matplotlib
import matplotlib.dates as dt

from functions import make_datelist
from import_data import import_lwf_data

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
# Available stations (so far):
# TreeNet sites of LWF: Jussy, Beatenberg, Lausanne, Lens, Neunkirch,
#                       Novaggio, Visp, Vordemwald, Schänis

treenetstation = 'Novaggio'

# needed for selecting proper CombiPrecip timeseries
station_id = treenetstation_id[treenetstationnames.index(treenetstation)]

#--------------------------------------------------------- 
# Define paths
#---------------------------------------------------------
stationrelhumpath = 'add path to MeteoSwiss relative humidity data'
treenetrelhumpath = 'add path to LWF relative humidity data'
figpath           = 'add path to your figures folder'

#--------------------------------------------------------- 
# Import precipitation data of TreeNet/LWF station
#---------------------------------------------------------
# Choose option: Load and process or load the processed data
relhum = import_lwf_data(treenetstation=treenetstation,path=treenetrelhumpath,\
                         variable='relhum',process_treenet_data='no')

#--------------------------------------------------------- 
# Load interpolated relative humidity data
#---------------------------------------------------------
print('import interpolated relative humidity data (version without vertical interpolation)')
interpolated_relhum_data = np.load(stationrelhumpath+'\interpolated_relhum_10minres_novertical.npy')
relhum['date_interpolation'] = interpolated_relhum_data[0,:]
interp_ind = treenetstationnames.index(treenetstation)
relhum['relhum_interpolated_novertical'] = interpolated_relhum_data[interp_ind+1,:].astype(float)

print('import interpolated relative humidity data (version with empirical gradient)')
interpolated_relhum_data = np.load(stationrelhumpath+'\interpolated_relhum_10minres_empiricalgradient.npy')
relhum['relhum_interpolated_empiricalgradient'] = interpolated_relhum_data[interp_ind+1,:].astype(float)

#--------------------------------------------------------- 
# Create monthly datelist
#---------------------------------------------------------
# Find latest starting date
if (relhum['treenetdate'][0] > relhum['date_interpolation'][0]):
       earliest_date = relhum['treenetdate'][0]
else:
       earliest_date = relhum['date_interpolation'][0]

# Add up to next month
year = earliest_date.year
if earliest_date.day != 1 or earliest_date.hour != 0 or earliest_date.minute!=0:
    month = earliest_date.month+1
    earliest_date = datetime.datetime(year,month,1,0,0)

# Find earliest ending date    
if (relhum['treenetdate'][-1] < relhum['date_interpolation'][-1]):
       latest_date = relhum['treenetdate'][-1]
else:
    latest_date = relhum['date_interpolation'][-1]
    
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
    fig = plt.figure()
    fig.set_size_inches(30, 14)
    
    #Get latest instant of the month
    month_end = nowdate+relativedelta(day=31,hour=23)
 
    # Find indices of the month in the timeseries
    wsl_ind_0    = np.where(relhum['treenetdate'] == nowdate)[0][0]
    interp_ind_0 = np.where(relhum['date_interpolation'] == nowdate)[0][0]
    wsl_ind_1    = np.where(relhum['treenetdate'] == month_end)[0][0]
    interp_ind_1 = np.where(relhum['date_interpolation'] == month_end)[0][0]
        
    # Get indices of wsl nans and meteoswiss nans
    wsl_nans                   = np.argwhere(np.isnan(relhum['treenetrelhum'][wsl_ind_0:wsl_ind_1+1]))[:,0]
    meteoswiss_nans_novertical = np.argwhere(np.isnan(relhum['relhum_interpolated_novertical'][interp_ind_0:interp_ind_1+1]))[:,0]
    meteoswiss_nans_empirical  = np.argwhere(np.isnan(relhum['relhum_interpolated_empiricalgradient'][interp_ind_0:interp_ind_1+1]))[:,0]

    # Plot    
    plt.plot(relhum['treenetdate'][wsl_ind_0:wsl_ind_1+1],\
             relhum['treenetrelhum'][wsl_ind_0:wsl_ind_1+1],ls='-',color='orange',label='treenet')
    plt.plot(relhum['date_interpolation'][interp_ind_0:interp_ind_1+1],\
             relhum['relhum_interpolated_novertical'][interp_ind_0:interp_ind_1+1],ls='-',color='blue',label='interpolation horizontal only')
    plt.plot(relhum['date_interpolation'][interp_ind_0:interp_ind_1+1],\
             relhum['relhum_interpolated_empiricalgradient'][interp_ind_0:interp_ind_1+1],ls='-',color='red',label='interpolation vertical empirical')
    plt.ylabel('Relative humidity [%]',fontsize=20)
    plt.title(treenetstation+' WSL',fontsize=25,loc='left')
    plt.title(nowdate.strftime('%b %Y'),fontsize=25,loc='right')
    plt.legend(loc=2)
    
    # Calculate RMSEs
    RMSE_novertical        = np.sqrt(np.nanmean((relhum['treenetrelhum'][wsl_ind_0:wsl_ind_1+1] - \
                                     relhum['relhum_interpolated_novertical'][interp_ind_0:interp_ind_1+1])**2))
    RMSE_empiricalgradient = np.sqrt(np.nanmean((relhum['treenetrelhum'][wsl_ind_0:wsl_ind_1+1] - \
                                     relhum['relhum_interpolated_empiricalgradient'][interp_ind_0:interp_ind_1+1])**2))
    
    # Add table with values below
    labels = [r'$RMSE_{novertical}$',r'$RMSE_{empirical}$']
    plt.table(cellText=[['%.2f' % RMSE_novertical+' %','%.2f' % RMSE_empiricalgradient+' %']],\
              bbox = [0.0,-0.12, 1.0, 0.07],cellLoc='center',rowLoc='center',colLabels=labels,fontsize=20)
    
    # Format x-axis and grid
    ax = plt.axes()
    datelist = make_datelist(relhum['treenetdate'][wsl_ind_0],relhum['treenetdate'][wsl_ind_1],1/6)
    ax.xaxis_date()
    plt.xticks(datelist[0::6*24*7])
    ax.xaxis.set_major_formatter(dt.DateFormatter('%d-%m-%Y %H'))
    ax.xaxis.grid(True, which='major', color='0.5',alpha=0.7, linestyle='--',lw=1.5)

    saveas='\\relhum_comparison_'+treenetstation+'_'
    plt.savefig(figpath+saveas+nowdate.strftime('%Y_%m')+'.png',bbox_inches='tight',dpi=400)
    plt.close('all')