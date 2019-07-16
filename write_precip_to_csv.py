# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 09:14:50 2019

@author: lukas jansing

file to write hourly precipitation data of:
    - CombiPrecip
    - Interpolated Precip
into csv files

"""

#--------------------------------------------------------- 
# Import modules
#---------------------------------------------------------
import numpy as np
import pandas as pd

#--------------------------------------------------------- 
# Define paths and variables
#---------------------------------------------------------
precippath      = 'add path to interpolated precipitation data'
combiprecippath = 'add path to combiprecip data'

variables = ['precip_interpolated1','precip_interpolated2','combiprecip']

#--------------------------------------------------------- 
# Define ID
#---------------------------------------------------------
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

treenetstation_id   = [1,2,3,4,5,6,7,8,1,10,11,12,13,14,7,16,17,7,1,20,21,1
                       ,20,24,25,26,27,28,29,30,31,32,28,34,27,36,37,25,39,40,41
                       ,42,43,20,45,46,47,48,49,50,51,52,53,54,54,56,14,14,14,60
                       ,43,43,63,63,65,65,67,1,1,1,71,72,73,74,24,24,77,78,79,80]

#--------------------------------------------------------- 
# Loop over the variables
#---------------------------------------------------------
for variable in variables:
    
    #--------------------------------------------------------- 
    # Load interpolated data from numpy file
    #---------------------------------------------------------
    print('Import data of '+variable)
    if variable == 'precip_interpolated1':
        interpolated_data = np.load(precippath+'\\interpolated_precipitation_version3.npy')
    if variable == 'precip_interpolated2':
        interpolated_data = np.load(precippath+'\\interpolated_precipitation_version4.npy')
    if variable == 'combiprecip':
        combiprecip = np.load(combiprecippath+'\\combiprecip_processed.npy')
        
    #--------------------------------------------------------- 
    # Write csv of interpolated or CombiPrecip data
    #---------------------------------------------------------
    if variable == 'precip_interpolated1' or variable == 'precip_interpolated2':

        # Define data array to write out
        interpolated_data_new = np.zeros((interpolated_data.shape[1]+2,\
                                            interpolated_data.shape[0])).astype(object)
    
        # Add site header to object array
        interpolated_data_new[0,:] = ['Sites'] + treenetstationnames
        
        # Add date header to header line and 2nd header line to object array
        interpolated_data_new[1,:] = ['Date'] + treenetstation_id
        
        # Add date column to object array
        # Start with 3rd row (first row is site, 2nd row is ID)
        for i in range(0,interpolated_data.shape[1]):
            interpolated_data_new[i+2,0] = interpolated_data[0,i].strftime('%Y-%m-%d-%H-%M')
            
        # Add precipitation to object array
        # Start with 3rd row and 2nd column (1st column is date)
        # Don't use first two rows (site and ID row)
        interpolated_data_new[2:,1:] = interpolated_data[1:,:].transpose()
        
        # Write out data using pandas
        print('Save interpolated data to csv')
        if variable == 'precip_interpolated1':
            pd.DataFrame(interpolated_data_new).to_csv(precippath+'\\interpolation_horizontal.csv',\
                        header=None,index=None,na_rep='NaN',sep=';')
        if variable == 'precip_interpolated2':
            pd.DataFrame(interpolated_data_new).to_csv(precippath+'\\interpolation_vertical.csv',\
                        header=None,index=None,na_rep='NaN',sep=';')
        
    if variable == 'combiprecip':
        
        # Create new data array with double columns (one for each station)
        combiprecip_processed = np.zeros((combiprecip[()]['combiprecip'].shape[0]+1,\
                                          len(treenetstationnames))).astype(object)
        
        # Add site header to object array
        combiprecip_processed[0,:] = ['Sites'] + treenetstationnames[:-1]
        
        # Add date header to header line and 2nd header line to object array
        combiprecip_processed[1,:] = ['Date'] + treenetstation_id[:-1]
        
        # Add date to first column (starting from 3rd row)
        for i in range(0,len(combiprecip[()]['date_combiprecip'])):
            combiprecip_processed[i+2,0] = combiprecip[()]['date_combiprecip'][i].strftime('%d.%m.%Y %H:%M')
        
        # Add precipitation to object array
        # Start with 3rd row and 2nd column
        # Add it column-wise
        for column in range(0,len(treenetstationnames)-1):
            id_combiprecip = combiprecip_processed[1,column+1]
            column_of_id = np.where(combiprecip[()]['combiprecip'][0,:] == id_combiprecip)[0][0] 
            combiprecip_processed[2:,column+1] = combiprecip[()]['combiprecip'][1:,column_of_id]

        pd.DataFrame(combiprecip_processed).to_csv(combiprecippath+'\\combiprecip_processed.csv',\
            header=None,index=None,na_rep='-9999999',sep=';')
    
