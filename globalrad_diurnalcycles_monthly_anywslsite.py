# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 12:35:00 2019

@author: lukas jansing

Script in order to validate the interpolated radiation data

Compare 10minres data to 10minres data of LWF

Comparison in the form of mean diurnal cycles for each month

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
import gc

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
# TreeNet sites of LWF: Jussy, Beatenberg, Lausanne, Neunkirch,
#                       Novaggio, Visp, Vordemwald, Schaenis

treenetstation = 'choose station'

# needed for selecting proper CombiPrecip timeseries
station_id = treenetstation_id[treenetstationnames.index(treenetstation)]

#--------------------------------------------------------- 
# Define paths
#---------------------------------------------------------
stationradpath = 'add path to MeteoSwiss radiation data'
treenetradpath = 'add path to LWF radiation data'
figpath        = 'add path to your figures folder'

#--------------------------------------------------------- 
# Import precipitation data of TreeNet/LWF station
#---------------------------------------------------------
# Choose option: Load and process or load the processed data
globalrad = import_lwf_data(treenetstation=treenetstation,path=treenetradpath,\
                            variable='globalrad',process_treenet_data='yes')

#--------------------------------------------------------- 
# Load interpolated radiation data
#---------------------------------------------------------
print('import interpolated radiation data (version without vertical interpolation)')
interpolated_globalrad_data = np.load(stationradpath+'\interpolated_globalrad_10minres_novertical.npy')
globalrad['date_interpolation'] = interpolated_globalrad_data[0,:]
interp_ind = treenetstationnames.index(treenetstation)
globalrad['globalrad_interpolated_novertical'] = interpolated_globalrad_data[interp_ind+1,:].astype(float)

#--------------------------------------------------------- 
# Create monthly datelist
#---------------------------------------------------------
# Find latest starting date
if (globalrad['treenetdate'][0] > globalrad['date_interpolation'][0]):
       earliest_date = globalrad['treenetdate'][0]
else:
       earliest_date = globalrad['date_interpolation'][0]

# Add up to next month
year = earliest_date.year
if earliest_date.day != 1 or earliest_date.hour != 0 or earliest_date.minute!=0:
    month = earliest_date.month+1
    earliest_date = datetime.datetime(year,month,1,0,0)

# Find earliest ending date    
if (globalrad['treenetdate'][-1] < globalrad['date_interpolation'][-1]):
       latest_date = globalrad['treenetdate'][-1]
else:
    latest_date = globalrad['date_interpolation'][-1]
    
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
    fig.set_size_inches(15, 15)
    
    #Get latest instant of the month
    month_end = nowdate+relativedelta(day=31,hour=23,minute=50)
 
    # Find indices of the month in the timeseries
    wsl_ind_0    = np.where(globalrad['treenetdate'] == nowdate)[0][0]
    interp_ind_0 = np.where(globalrad['date_interpolation'] == nowdate)[0][0]
    wsl_ind_1    = np.where(globalrad['treenetdate'] == month_end)[0][0]
    interp_ind_1 = np.where(globalrad['date_interpolation'] == month_end)[0][0]
    
    # Reshape data to daily blocks
    globalrad['treenetdate_diurnal'] = np.reshape(globalrad['treenetdate'][wsl_ind_0:wsl_ind_1+1],\
             (int(globalrad['treenetdate'][wsl_ind_0:wsl_ind_1+1].shape[0]/(6*24)),6*24))
    globalrad['treenetglobalrad_diurnal'] = np.reshape(globalrad['treenetglobalrad'][wsl_ind_0:wsl_ind_1+1],\
             (int(globalrad['treenetglobalrad'][wsl_ind_0:wsl_ind_1+1].shape[0]/(6*24)),6*24))
    globalrad['date_interpolation_diurnal'] = np.reshape(globalrad['date_interpolation'][interp_ind_0:interp_ind_1+1],\
             (int(globalrad['date_interpolation'][interp_ind_0:interp_ind_1+1].shape[0]/(6*24)),6*24))
    globalrad['globalrad_interpolated_diurnal'] = np.reshape(globalrad['globalrad_interpolated_novertical'][interp_ind_0:interp_ind_1+1],\
             (int(globalrad['globalrad_interpolated_novertical'][interp_ind_0:interp_ind_1+1].shape[0]/(6*24)),6*24))
        
    # Calculate mean diurnal cycles
    globalrad_treenet_mean       = []
    globalrad_interpolation_mean = []
    for i in range(0,6*24):
        globalrad_treenet_mean.append(np.nanmean(globalrad['treenetglobalrad_diurnal'][:,i]))
        globalrad_interpolation_mean.append(np.nanmean(globalrad['globalrad_interpolated_diurnal'][:,i]))
    globalrad['treenetglobalrad_mean_diurnal'] = np.array(globalrad_treenet_mean)
    globalrad['globalrad_interpolated_mean_diurnal'] = np.array(globalrad_interpolation_mean)

    # Plot    
    plt.plot(globalrad['treenetdate_diurnal'][0,:],\
             globalrad['treenetglobalrad_mean_diurnal'],ls='-',color='orange',label='treenet',lw=2)
    plt.plot(globalrad['date_interpolation_diurnal'][0,:],globalrad['globalrad_interpolated_mean_diurnal'],\
             ls='-',color='blue',label='interpolation horizontal only',lw=2)
    plt.ylabel(r'Global radiation [$W/m^2$]',fontsize=20)
    plt.title(treenetstation+' WSL',fontsize=25,loc='left')
    plt.title(nowdate.strftime('%b %Y'),fontsize=25,loc='right')
    plt.legend(loc=2)
    
    # Calculate RMSEs
    RMSE_novertical = np.sqrt(np.nanmean((globalrad['treenetglobalrad_mean_diurnal'] - \
                              globalrad['globalrad_interpolated_mean_diurnal'])**2))
    
    # Add table with values below
    labels = [r'$RMSE_{novertical}$']
    plt.table(cellText=[['%.2f' % RMSE_novertical+r' $W/m^2$']],\
              bbox = [0.0,-0.12, 1.0, 0.07],cellLoc='center',rowLoc='center',colLabels=labels,fontsize=20)
    
    # Format x-axis and grid
    ax = plt.axes()
    ax.xaxis_date()
    plt.xticks(globalrad['treenetdate_diurnal'][0,:][::36])
    ax.xaxis.set_major_formatter(dt.DateFormatter('%m-%Y %H'))
    ax.xaxis.grid(True, which='major', color='0.5',alpha=0.7, linestyle='--',lw=1.5)

    saveas='\\globalrad_diurnal_'+treenetstation+'_'
    plt.savefig(figpath+saveas+nowdate.strftime('%Y_%m')+'.png',bbox_inches='tight',dpi=400)
    plt.close('all')
    gc.collect()