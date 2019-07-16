# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 13:20:57 2019

@author: lukas jansing

routines in order to import data

"""

import numpy as np
import os
import datetime

#--------------------------------------------------------- 
# Import precipitation data of MeteoSwiss stations
#---------------------------------------------------------
def import_precipitation_data(stationnames,resolution,path):
    
    data = {}

    # Loop over all stations
    for i in range(0,len(stationnames)):

        # Load data of respective station
        station       = stationnames[i]
        
        if resolution == 'hourly':
            print('loading data of '+station)
            data[station] = np.loadtxt(path+'\precip_data_2011010100_2018123123_'+station+'.dat',unpack=True,skiprows=9)
        
        if resolution == '10min':
            # Note: Data is distributed over 3 files
            # Have to be appended to the same data structure
            # Some data is not available for the first time period (2011-2014)
            # Therefore append nans there
            
            data[station] = []
            print('loading data of '+station+' (1/3)')
            if (os.stat(path+'\precip_data_2011010100_2013123123_'+station+'.dat').st_size == 0) == False:
                data[station] = np.loadtxt(path+'\precip_data_2011010100_2013123123_'+station+'.dat',unpack=True,skiprows=9)
                if i == 0:
                    shape_of_file = data[station].shape # save shape of first file to fill missings with nans
            else:
                data[station] = np.full(shape_of_file,np.nan)
            print('loading data of '+station+' (2/3)')
            data[station] = np.append(data[station],np.loadtxt(path+'\precip_data_2014010100_2016123123_'+station+'.dat',unpack=True,skiprows=9),axis=1)
            print('loading data of '+station+' (3/3)')
            data[station] = np.append(data[station],np.loadtxt(path+'\precip_data_2017010100_2019052123_'+station+'.dat',unpack=True,skiprows=9),axis=1)
            
        # NaN handling: Convert the fill numbers (32767) into nans
        data[station][data[station] == 32767] = np.nan
            
    return data

#--------------------------------------------------------- 
# Import other data of MeteoSwiss (temp, rad, relhum)
#---------------------------------------------------------
def import_meteoswiss_data(stationnames,variable,path,process,datelist_control):

   if process == 'yes':

        #--------------------------------------------------------- 
        # Import data from all meteoswiss stations
        #---------------------------------------------------------
        data              = {}
        
        # Loop over all stations
        for i in range(0,len(stationnames)):
        
            # Load data of respective station
            station       = stationnames[i]
            print('loading data of '+station)
            data[station] = np.genfromtxt(path+'\\'+variable+'_data_2011010100_2019052123_'+station+'.txt',unpack=True,\
                                          missing_values='-',skip_header=3,usecols=[1,2],delimiter=';')
            
        #--------------------------------------------------------- 
        # Combine data from the stations where there are two files
        #---------------------------------------------------------
        for i in range(0,len(stationnames)):
            station = stationnames[i]
            # If part 1 of stationdata --> save to new overall key
            if '1' in station:
                data[station[:-1]] = data[station]
                del data[station]
            # If part 2 of stationdata --> add to overall key
            if '2' in station:
                data[station[:-1]] = np.hstack((data[station[:-1]],data[station]))
                del data[station]
                
        # Create new stationnames array
        stationnames_new = list(data.keys())
        
        #--------------------------------------------------------- 
        # Fill up data gaps with NaNs
        #---------------------------------------------------------
        data_new = {}
        # Loop over all stations and apply filling if not full data are available
        for station in data.keys():
        #for station in ['bergünlatsch']:
            
            data_new[station] = [[],[]]
            
            # Case where there are missing data points
            if data[station].shape[1] != len(datelist_control):
                
                print('process data of '+station)
                
                # Read out datelist of the station
                datelist_test = []
                for i in range(0,data[station].shape[1]):
                    datelist_test.append(datetime.datetime(int(str(data[station][0,i])[0:4]),\
                                                           int(str(data[station][0,i])[4:6]),\
                                                           int(str(data[station][0,i])[6:8]),\
                                                           int(str(data[station][0,i])[8:10]),\
                                                           int(str(data[station][0,i])[10:12])))
                
                # Convert the test datelist to a set
                datelist_test_unordered = set(datelist_test)
                    
                # Compare test datelist to control datelist
                ten_minutes = datetime.timedelta(minutes=10)
                test_date = datelist_control[0]
                while test_date <= datelist_control[-1]:
                
                    # In case the respective date is missing in the station data
                    if test_date not in datelist_test_unordered:
                        data_new[station][0].append(test_date)
                        data_new[station][1].append(np.nan)
                        
                    # In case the respecitve date is in the station data
                    else:
                        now_ind = datelist_test.index(test_date)
                        data_new[station][0].append(datelist_test[now_ind])
                        data_new[station][1].append(data[station][1,now_ind])
                        
                    test_date += ten_minutes
                    
            # Case where there are no missing data points
            else:
                data_new[station][0] = datelist_control
                data_new[station][1] = data[station][1,:]
        
            # Final call to convert list structure to numpy array
            data_new[station] = np.array(data_new[station])

        #--------------------------------------------------------- 
        # Save processed station data for further use
        #---------------------------------------------------------
        save_stationdata = 'yes'
        #save_stationdata = 'no'

        if save_stationdata == 'yes':
            np.save(path+'\stationdata_processed.npy',data_new)
                
   else:
        print('load processed station data')
        
        if variable == 'rad':
            data_part1 = np.load(path+'\stationdata_processed_part1.npy')
            data_part2 = np.load(path+'\stationdata_processed_part2.npy')
            data_new   = data_part1
            data_new[()].update(data_part2[()])

        else:
            data_new = np.load(path+'\stationdata_processed.npy')
    
        # Create new stationnames array
        stationnames_new = list(data_new[()].keys())
                
                
   return data_new,stationnames_new

#--------------------------------------------------------- 
# Import LWF precipitation data
#---------------------------------------------------------
def import_lwf_precipitation_data(treenetstation,treenetprecip_res,path,process_schänis):
    
    #--------------------------------------------------------- 
    # Import modules and functions
    #---------------------------------------------------------
    import csv
    from functions import make_datelist
    
    #--------------------------------------------------------- 
    # Import precipitation data of TreeNet station
    #---------------------------------------------------------
    # Load data from csv or txt file
    print('import wsl data of '+treenetstation)
    
    # Different handling for the csv files than for txt files
    if treenetstation in ['Pfynwald','Bachtel','Muri_Beech','Muri_Spruce','Muri_Meteo',
                          'Riehen_Forest','Riehen_Meteo','Grosswangen']:
        treenet_data = []
        f = open(path+'\\'+treenetstation+'_precipitation.csv','r')
        reader = csv.reader(f)
        for row in reader:
            treenet_data.append(row)
        
        treenet_data = np.array(treenet_data)
    
        # Extract precipitation and date
        precip = {}
        precip['date_UTC+1'] = []
        precip['treenetprecip']     = []
    
        for i in range(1,treenet_data.shape[0]):
            
            # Extract date from strings and convert to datetime format
            precip['date_UTC+1'].append(datetime.datetime(int(treenet_data[i,2][0:4]),int(treenet_data[i,2][5:7]),\
                                      int(treenet_data[i,2][8:10]),int(treenet_data[i,2][11:13]),int(treenet_data[i,2][14:16])))
            # Extract precipitation
            if treenet_data[i,3] != 'NA':
                precip['treenetprecip'].append(treenet_data[i,3].astype(np.float))
            elif treenet_data[i,3] == 'NA':
                precip['treenetprecip'].append(np.nan)
            else: # Only to check for strange cases
                print('strange')
    
        precip['treenetprecip'] = np.array(precip['treenetprecip'])
        
        #--------------------------------------------------------- 
        # Calculate hourly sums and convert time if data sourc is csv
        #---------------------------------------------------------
        # First check if first and last entry are start and finish of an hour
        # identify the indices to read out correct data
        for i in range(0,6):
            if precip['date_UTC+1'][i].minute == 10:
                starting_ind = i
                break
            else:
                continue
        
        for i in range(-1,-7,-1):
            if precip['date_UTC+1'][i].minute == 0:
                ending_ind = i
                break
            else:
                continue
        
        # Reshape the list into hourly blocks in order to calculate sums
        # Save shape beforehand
        if ending_ind != -1:
            treenetprecip_shape = len(precip['treenetprecip'][starting_ind:ending_ind+1])
            for key in precip.keys():
                precip[key] = np.reshape(precip[key][starting_ind:ending_ind+1],\
                      (int(treenetprecip_shape/6),6))
        else:
            for key in precip.keys():
                treenetprecip_shape = len(precip['treenetprecip'][starting_ind:])
                precip[key] = np.reshape(precip[key][starting_ind:],\
                              (int(treenetprecip_shape/6),6))
    
        # Calculate hourly sums
        precip['precip_hourly'] = []
        for i in range(0,int(treenetprecip_shape/6)):
            precip['precip_hourly'].append(sum(precip['treenetprecip'][i,:]))
    
        precip['precip_hourly'] = np.array(precip['precip_hourly'])
    
        # Create hourly datelist
        # Note: Date always refers to precip from the last hour!!!
        precip['date_hourly_UTC+1'] = make_datelist(precip['date_UTC+1'][0,-1],\
                      precip['date_UTC+1'][-1,-1],1)
    
        # Convert to UTC
        precip['date_hourly_UTC'] = []
        for i in range(0,len(precip['date_hourly_UTC+1'])):
            precip['date_hourly_UTC'].append(precip['date_hourly_UTC+1'][i]-datetime.timedelta(hours=1)) 
            
        precip['date_hourly_UTC'] = np.array(precip['date_hourly_UTC'])
        
    #--------------------------------------------------------- 
    # Import WSL LWF data
    # Handling for the LWF data (in textfiles)
    #---------------------------------------------------------     
    else:
        
        # Convert the string of the stationname for special cases
        if treenetstation == 'LWF-Lens3':
            treenetstationname_new = 'lens'
        elif 'Neunkirch' in treenetstation:
            treenetstationname_new = 'neunkirch'
        elif treenetstation == 'Schaenis':
            treenetstationname_new = 'schänis'
        else:
            treenetstationname_new = treenetstation.lower()
            
        #---------------------------------------------------------
        # Possibility to import 10minres data and convert them to hourly
        #---------------------------------------------------------
        # case where hourly data is used
        if treenetprecip_res == 'hourly':
            # Import the data
            treenet_data = {}
            treenet_data[treenetstation] = np.loadtxt(path+'\precip_data_2011010100_2018123123_'+treenetstationname_new+'.dat',\
                                                        unpack=True,skiprows=9)
            
            # NaN handling: Convert the fill numbers (32767) into nans
            treenet_data[treenetstation][treenet_data[treenetstation] == 32767] = np.nan
            
            # Create datelist
            startdate = datetime.datetime(int(treenet_data[treenetstation][1,0]),int(treenet_data[treenetstation][2,0]),\
                                          int(treenet_data[treenetstation][3,0]),int(treenet_data[treenetstation][4,0]))   
            enddate   = datetime.datetime(int(treenet_data[treenetstation][1,-1]),int(treenet_data[treenetstation][2,-1]),\
                                          int(treenet_data[treenetstation][3,-1]),int(treenet_data[treenetstation][4,-1]))    
            hstep     = 1
        
            # Save precipitation to numpy array    
            precip = {}
            precip['date_hourly_UTC'] = make_datelist(startdate,enddate,hstep)
            precip['date_hourly_UTC'] = np.array(precip['date_hourly_UTC'])
            precip['precip_hourly'] = np.array(treenet_data[treenetstation][6,:])
                
        # case where 10minres data is used
        elif treenetprecip_res == '10minres':
            
            # create exception for Schaenis (different data format)
            if treenetstation == 'Schaenis':
                
               if process_schänis == 'yes':
                   
                    # Import the data
                    treenet_data = {} 
                    treenet_data[treenetstation] = np.genfromtxt(path+'\\precip_data_2011010100_2019052123_'+treenetstationname_new+'_freiland.txt',\
                                                     unpack=True,missing_values='-',skip_header=3,usecols=[1,2],delimiter=';')
                    # Create control datelist
                    startdate = datetime.datetime(2013,1,1,0,0)   
                    enddate   = datetime.datetime(2019,5,21,23,50)
                    hstep     = 1/6
        
                    datelist_control  = make_datelist(startdate,enddate,hstep)
                    
                    print('process data of '+treenetstation)
            
                    # Read out datelist of the station
                    datelist_test = []
                    for i in range(0,treenet_data[treenetstation].shape[1]):
                        datelist_test.append(datetime.datetime(int(str(treenet_data[treenetstation][0,i])[0:4]),\
                                                       int(str(treenet_data[treenetstation][0,i])[4:6]),\
                                                       int(str(treenet_data[treenetstation][0,i])[6:8]),\
                                                       int(str(treenet_data[treenetstation][0,i])[8:10]),\
                                                       int(str(treenet_data[treenetstation][0,i])[10:12])))
                    
                    # Identify and fill gaps
                    treenet_data_new = {}
                    treenet_data_new[treenetstation] = [[],[]]
            
                    # Convert the test datelist to a set
                    datelist_test_unordered = set(datelist_test)
                
                    # Compare test datelist to control datelist
                    ten_minutes = datetime.timedelta(minutes=10)
                    test_date = datelist_control[0]
                    while test_date <= datelist_control[-1]:
            
                        # In case the respective date is missing in the station data
                        if test_date not in datelist_test_unordered:
                            treenet_data_new[treenetstation][0].append(test_date)
                            treenet_data_new[treenetstation][1].append(np.nan)
                    
                        # In case the respecitve date is in the station data
                        else:
                            now_ind = datelist_test.index(test_date)
                            treenet_data_new[treenetstation][0].append(datelist_test[now_ind])
                            treenet_data_new[treenetstation][1].append(treenet_data[treenetstation][1,now_ind])
                    
                        test_date += ten_minutes
                    
                    # Save precipitation to numpy array
                    precip = {}
                    precip['date_10minres_UTC'] = np.array(datelist_control)
                    precip['precip_10minres']   = np.array(treenet_data_new[treenetstation])[1,:].astype(float)
                    
                    # Save processed Schänis data
                    np.save(path+'\Schänis_processed.npy',precip)
                    
               if process_schänis == 'no':
                    precip = np.load(path+'\Schänis_processed.npy')
                    precip = precip[()]
                    
               # Correct for the unrealistic values
               count = 0
               for i in range(0,precip['precip_10minres'].shape[0]):
                   if precip['precip_10minres'][i] > 15:
                       if precip['precip_10minres'][i-1] > 15 and \
                          precip['precip_10minres'][i+1] > 15:
                              precip['precip_10minres'][i] = np.nan
                              count += 1
                   if precip['date_10minres_UTC'][i].year > 2017:
                     precip['precip_10minres'][i] = np.nan
                   if precip['date_10minres_UTC'][i].year < 2014:
                     precip['precip_10minres'][i] = np.nan
                     
            # Here starts the "normal" procedure
            else:
                
                # Import the data
                treenet_data = {}
                treenet_data[treenetstation] = np.loadtxt(path+'\precip_data_2011010100_2013123123_'+treenetstationname_new+'.dat',\
                                                          unpack=True,skiprows=9)
                treenet_data[treenetstation] = np.append(treenet_data[treenetstation],np.loadtxt(path+'\precip_data_2014010100_2016123123_'+treenetstationname_new+'.dat',\
                                                         unpack=True,skiprows=9),axis=1)
                treenet_data[treenetstation] = np.append(treenet_data[treenetstation],np.loadtxt(path+'\precip_data_2017010100_2019052123_'+treenetstationname_new+'.dat',\
                                                         unpack=True,skiprows=9),axis=1)
                
                # NaN handling: Convert the fill numbers (32767) into nans
                treenet_data[treenetstation][treenet_data[treenetstation] == 32767] = np.nan
                
                # Create datelist
                startdate = datetime.datetime(int(treenet_data[treenetstation][1,0]),int(treenet_data[treenetstation][2,0]),\
                                              int(treenet_data[treenetstation][3,0]),int(treenet_data[treenetstation][4,0]),int(treenet_data[treenetstation][5,0]))   
                enddate   = datetime.datetime(int(treenet_data[treenetstation][1,-1]),int(treenet_data[treenetstation][2,-1]),\
                                              int(treenet_data[treenetstation][3,-1]),int(treenet_data[treenetstation][4,-1]),int(treenet_data[treenetstation][5,-1]))    
                hstep     = 1/6
                
                # Save precipitation to numpy array
                precip = {}
                precip['date_10minres_UTC'] = make_datelist(startdate,enddate,hstep)
                precip['date_10minres_UTC'] = np.array(precip['date_10minres_UTC'])
                precip['precip_10minres'] = np.array(treenet_data[treenetstation][6,:])
    
            # Calculate hourly values
            for i in range(0,6):
                if precip['date_10minres_UTC'][i].minute == 10:
                    starting_ind = i
                    break
                else:
                    continue
    
            for i in range(-1,-7,-1):
                if precip['date_10minres_UTC'][i].minute == 0:
                    ending_ind = i
                    break
                else:
                    continue
        
            if ending_ind != -1:   
                precip_shape = precip['precip_10minres'][starting_ind:ending_ind+1].shape[0]
                precip['precip_10minres_reshaped']   = np.reshape(precip['precip_10minres'][starting_ind:ending_ind+1],\
                                                       (int(precip_shape/6),6))
                precip['date_10minres_UTC_reshaped'] = np.reshape(precip['date_10minres_UTC'][starting_ind:ending_ind+1],\
                                                       (int(precip_shape/6),6))
            else:
                precip_shape = precip['precip_10minres'][starting_ind:].shape[0]
                precip['precip_10minres_reshaped']   = np.reshape(precip['precip_10minres'][starting_ind:],\
                                                       (int(precip_shape/6),6))
                precip['date_10minres_UTC_reshaped'] = np.reshape(precip['date_10minres_UTC'][starting_ind:],\
                                                       (int(precip_shape/6),6))
    
    
            # Calculate hourly sums
            precip['precip_hourly'] = []
            for i in range(0,int(precip_shape/6)):
                precip['precip_hourly'].append(sum(precip['precip_10minres_reshaped'][i,:]))
            
            precip['precip_hourly'] = np.array(precip['precip_hourly'])
            
            # Create hourly datelist
            # Note: Date always refers to precip from the last hour!!!
            precip['date_hourly_UTC'] = make_datelist(precip['date_10minres_UTC_reshaped'][0,-1],\
                                                      precip['date_10minres_UTC_reshaped'][-1,-1],1)
            precip['date_hourly_UTC'] = np.array(precip['date_hourly_UTC'])
        
    return precip

#--------------------------------------------------------- 
# Import other LWF data (temp,rad,relhum)
#---------------------------------------------------------
def import_lwf_data(treenetstation,path,variable,process_treenet_data):
    
    from functions import make_datelist
    
    # Load data from csv file
    print('import wsl data of '+treenetstation)
    
    # Convert the string of the stationname
    if treenetstation == 'LWF-Lens3':
        treenetstationname_new = 'lens'
    elif 'Neunkirch' in treenetstation:
        treenetstationname_new = 'neunkirch'
    elif treenetstation == 'Schaenis':
        treenetstationname_new = 'schänis'
    else:
        treenetstationname_new = treenetstation.lower()
    
    if process_treenet_data == 'yes':
    
        # Import the data
        treenet_data = {}
        treenet_data[treenetstation] = np.genfromtxt(path+'\\'+variable+'_data_2011010100_2019052123_'+treenetstationname_new+'.txt',\
                                                     unpack=True,missing_values='-',skip_header=3,usecols=[1,2],delimiter=';')
        
        # Create control datelist
        startdate = datetime.datetime(2011,1,1,0,0)      # earlier none of the files have data (checked manually)   
        enddate   = datetime.datetime(2019,5,21,23,50)
        hstep     = 1/6
        
        datelist_control  = make_datelist(startdate,enddate,hstep)
        
        # Identify and fill gaps if there are any
        treenet_data_new = {}
        treenet_data_new[treenetstation] = [[],[]]
        
        if treenet_data[treenetstation].shape[1] != len(datelist_control):
                    
            print('process data of '+treenetstation)
            
            # Read out datelist of the station
            datelist_test = []
            for i in range(0,treenet_data[treenetstation].shape[1]):
                datelist_test.append(datetime.datetime(int(str(treenet_data[treenetstation][0,i])[0:4]),\
                                                       int(str(treenet_data[treenetstation][0,i])[4:6]),\
                                                       int(str(treenet_data[treenetstation][0,i])[6:8]),\
                                                       int(str(treenet_data[treenetstation][0,i])[8:10]),\
                                                       int(str(treenet_data[treenetstation][0,i])[10:12])))
            
            # Convert the test datelist to a set
            datelist_test_unordered = set(datelist_test)
                
            # Compare test datelist to control datelist
            ten_minutes = datetime.timedelta(minutes=10)
            test_date = datelist_control[0]
            
            while test_date <= datelist_control[-1]:
            
                # In case the respective date is missing in the station data
                if test_date not in datelist_test_unordered:
                    treenet_data_new[treenetstation][0].append(test_date)
                    treenet_data_new[treenetstation][1].append(np.nan)
                    
                # In case the respecitve date is in the station data
                else:
                    now_ind = datelist_test.index(test_date)
                    treenet_data_new[treenetstation][0].append(datelist_test[now_ind])
                    treenet_data_new[treenetstation][1].append(treenet_data[treenetstation][1,now_ind])
                    
                test_date += ten_minutes
                
        # Case where there are no missing data points
        else:
            treenet_data_new[treenetstation][0] = datelist_control
            treenet_data_new[treenetstation][1] = treenet_data[treenetstation][1,:]
        
        # Final call to convert list structure to numpy arrays
        treenet_final = {}
        treenet_final['treenet'+variable] = np.array(treenet_data_new[treenetstation])[1,:].astype(float)
        treenet_final['treenetdate']      = np.array(treenet_data_new[treenetstation])[0,:]
        
        # Save data to numpy file
        np.save(path+'\\'+treenetstation+'_processed.npy',treenet_final)
        
    if process_treenet_data == 'no':
        treenet_final = np.load(path+'\\'+treenetstation+'_processed.npy')[()]
        
    return treenet_final

#--------------------------------------------------------- 
# Import CombiPrecip data
#---------------------------------------------------------       
def import_combiprecip(combiprecippath,processing_combiprecip,save_combiprecip):
    
    if processing_combiprecip == 'yes':

        # Load the data
        print('import combiprecip data')
        combiprecip_data = []
        for i in range(2005,2020):
        #for i in range(2012,2013):
            print(i)
            f = open(combiprecippath+'\T_lor500_'+str(i)+'_o\prad_lor500.s'+str(i),'r')
            for line in f:
                line = line.strip()
                line = line.split()
                combiprecip_data.append(line)
            
        # Identify the unnecessary lines (except for id line: keep it once)
        combiprecip_list =  []
        trash_ind = []
        trash_ind = [[i for i,val in enumerate(combiprecip_data) if val==combiprecip_data[0]]]
        trash_ind.append([i for i,val in enumerate(combiprecip_data) if val==combiprecip_data[1]])
        trash_ind.append([i for i,val in enumerate(combiprecip_data) if val==combiprecip_data[2]])
        trash_ind = [item for sublist in trash_ind for item in sublist]
        
        # Get rid of them by copying everything into combiprecip_intermediate
        # except for the unnecessary lines
        for i in range(0,len(combiprecip_data)):
            if i in trash_ind and i != 1:
                continue
            else:
                combiprecip_list.append(combiprecip_data[i])
        combiprecip_intermediate = np.array(combiprecip_list)
        
        # Convert the strings into int and float
        combiprecip = np.zeros(np.shape(combiprecip_intermediate))
        for i in range(0,combiprecip.shape[0]):
            if i == 0:
                for j in range(4,combiprecip_intermediate.shape[1]):
                    combiprecip[i,j] = np.int(combiprecip_intermediate[i,j])
            else:
                for j in range(0,combiprecip_intermediate.shape[1]):
                    combiprecip[i,j] = np.float(combiprecip_intermediate[i,j])
        
        #--------------------------------------------------------- 
        # Identify gaps and fill them with nans
        #---------------------------------------------------------
        print('process combiprecip data')
        
        # create control datelist containing all dates from start to end
        datelist_control = []
        startdate = datetime.datetime(int(combiprecip[1,0]),int(combiprecip[1,1]),\
                                      int(combiprecip[1,2]),int(combiprecip[1,3]))
        enddate = datetime.datetime(int(combiprecip[-1,0]),int(combiprecip[-1,1]),\
                                      int(combiprecip[-1,2]),int(combiprecip[-1,3]))
        datelist_control = make_datelist(startdate,enddate,1)
        
        # extract datelist from files
        # note: change 24 UTC to 00 UTC of next day
        # start with 1 because of id line!!!
        datelist_test = []
        count = 0
        for i in range(1,combiprecip.shape[0]):
            if int(combiprecip[i,3]) == 24:
                datelist_test.append(datelist_test[i-2]+datetime.timedelta(hours=1))
            else:    
                datelist_test.append(datetime.datetime(int(combiprecip[i,0]),int(combiprecip[i,1]),\
                                                       int(combiprecip[i,2]),int(combiprecip[i,3])))
        
        # define new combiprecip container
        combiprecip_new = {}
        combiprecip_new['date_combiprecip'] = []
        combiprecip_new['combiprecip']      = []
        # Fill in id's in first line
        combiprecip_new['combiprecip']      = [combiprecip[0,4:]]
        
        # loop over datelist and check for missing values
        # fill up missing values in datelist
        # append nans at the respective dates
        count = 0
        one_hour = timedelta(hours=1)
        test_date = datelist_control[0]
        while test_date <= datelist_control[-1]:
            
            # case if the respective date is missing in combiprecip data
            if test_date not in datelist_test:
                print(test_date)
                combiprecip_new['date_combiprecip'].append(test_date)
                combiprecip_new['combiprecip'].append(np.full(combiprecip.shape[1]-4,np.nan))
                count += 1
                
            # case if the respecitve date is in the combiprecip data
            else:
                now_ind = datelist_test.index(test_date)
                combiprecip_new['date_combiprecip'].append(datelist_test[now_ind])
                combiprecip_new['combiprecip'].append(combiprecip[now_ind+1,4:])
                
            test_date += one_hour
        
        combiprecip_new['combiprecip'] = np.array(combiprecip_new['combiprecip'])
        
        if save_combiprecip == 'yes':
            np.save(combiprecippath+'\combiprecip_processed.npy',combiprecip_new)
        
        return combiprecip_new

    #--------------------------------------------------------- 
    # Import processed combiprecip data
    #---------------------------------------------------------
    if processing_combiprecip == 'no':
        print('import processed combiprecip data')
        combiprecip_data = np.load(combiprecippath+'\combiprecip_processed.npy')
        
        return combiprecip_data
    
       
