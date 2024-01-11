
import h5py 
from NuRadioReco.utilities import units
from NuRadioReco.detector import detector
from NuRadioMC.utilities import fluxes
from NuRadioMC.utilities.Veff import get_Veff_Aeff, get_Veff_Aeff_array, get_index, get_Veff_water_equivalent
import numpy as np
from dash.dependencies import Input, Output, State
from dash import dcc
from dash import html
import dash
import json
import uuid
import glob
import os
import math
import sys

temp_stdout = None
# Disable
def blockPrint():
    global temp_stdout
    temp_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    global temp_stdout
    sys.stdout = temp_stdout

#################################################################
#Print out all the attributes and keys from the hdf5 file
def get_attr_data(PATH):
    f = h5py.File(PATH, mode='r')
    attributes = f.attrs
    attrs = []
    Keys = []
    for i in attributes:
        attrs.append(i)
    for i in f.keys():
        Keys.append(i)
    print("Attributes of hdf5:",attrs)
    print("Keys of hdf5:", Keys)

#####################################################################
#Calculating the effective volume from the output file:
def get_Veff_from_file(PATH):
    try:
        blockPrint()
        data = get_Veff_Aeff(PATH)
        Veffs, energies, energies_low, energies_up, zenith_bins, utrigger_names = get_Veff_Aeff_array(data)
            # calculate the average over all zenith angle bins (in this case only one bin that contains the full sky)
        Veff = np.average(Veffs[:, :, get_index("all_triggers", utrigger_names), 0], axis=1)
            # we also want the water equivalent effective volume times 4pi
        Veff = get_Veff_water_equivalent(Veff) * 4 * np.pi
            # calculate the uncertainty for the average over all zenith angle bins. The error relative error is just 1./sqrt(N)
        Veff_error = Veff / np.sum(Veffs[:, :, get_index("all_triggers", utrigger_names), 2], axis=1) ** 0.5
        
            #energies = energies / units.eV
            #Veff / units.km ** 3 / units.sr
            #Veff_error = Veff_error / units.km ** 3 / units.sr
        enablePrint()
        print(Veff, Veff_error, energies)
        return Veff, Veff_error, energies
    except:
        return  np.array([]),  np.array([]),  np.array([])

def get_Veff_data(PATH):
    if os.path.isfile(PATH):
        Veff, Veff_error, energies = get_Veff_from_file(PATH)
        return Veff, Veff_error, energies
    
    elif os.path.isdir(PATH):
        Veff = np.array([])
        Veff_error = np.array([])
        energies = np.array([])
        
        for filename in os.listdir(PATH):
            file = os.path.join(PATH, filename)
            V_i, Ve_i, e_i = get_Veff_from_file(file)
            Veff = np.append(Veff, V_i)
            Veff_error = np.append(Veff_error, Ve_i)
            energies = np.append(energies, e_i)
        return np.mean(Veff), np.mean(Veff_error), np.mean(energies)
    
    else:
        print("No file or folder found, Please try a different path")
    

#######################################################################
#Get shower data from the output file:
def get_shower_data_from_file(PATH):
    f = h5py.File(PATH, mode='r')
    #The shower ID
    shower_id = np.array(f["shower_ids"])
    #The shower position
    xx = np.array(f["xx"])
    yy = np.array(f["yy"])
    zz = np.array(f["zz"])
    #Energy of the shower which is used to determine the radio emission
    shower_energies = np.array(f["shower_energies"])
    #Type of the shower (so far we only have “em” and “had”)
    shower_type = np.array(f["shower_type"])
    triggered = np.array(f["triggered"])
    return shower_id, shower_energies, xx, yy, zz, shower_type, triggered

def get_shower_data(PATH):
    if os.path.isfile(PATH):
        shower_id, shower_energies, xx, yy, zz, shower_type, triggered = get_shower_data_from_file(PATH)
        return shower_id, shower_energies, xx, yy, zz, shower_type, triggered
    
    elif os.path.isdir(PATH):
        shower_id, shower_energies, xx, yy, zz, shower_type, triggered = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
        for filename in os.listdir(PATH):
            file = os.path.join(PATH, filename)
            a, b, c, d, f, g, h = get_shower_data_from_file(file)
            shower_id = np.append(shower_id, a)
            shower_energies = np.append(shower_energies, b)
            xx = np.append(xx, c)
            yy = np.append(yy, d)
            zz = np.append(zz, f)
            shower_type = np.append(shower_type, g)
            triggered = np.append(triggered, h)
        return shower_id, shower_energies, xx, yy, zz, shower_type, triggered
    
    else:
        print("No file or folder found, Please try a different path")



###################################################################################
#Get station data from the out put file aka hdf5 file:
def get_station_data_from_file(PATH, station_ID, channel_ID):
    f = h5py.File(PATH, mode='r')
    shower_id = np.array([])
    time_shower_and_ray = np.array([])
    max_amp_shower_and_ray = np.array([])
    travel_distances = np.array([])
    for s_id in station_ID:
        try:
            total_shower = len(f[f"station_{s_id}"]['time_shower_and_ray'])
            for c_id in channel_ID:
                #The time travelled by each ray tracing solution to a specific channel
                #Maximum amplitude per shower, channel and ray tracing solution.
                    for i in range(total_shower):
                        shower = f[f'station_{s_id}']['shower_id'][i]
                        time_per_shower = f[f'station_{s_id}']['time_shower_and_ray'][i][c_id]
                        max_amp_per_shower = f[f'station_{s_id}']['max_amp_shower_and_ray'][i][c_id]
                        travel_distance_per_shower = f[f'station_{s_id}']['travel_distances'][i][c_id]
                    #############################################################################
                        max_amp_shower_and_ray = np.append(max_amp_shower_and_ray, max_amp_per_shower)
                        time_shower_and_ray = np.append(time_shower_and_ray, time_per_shower)
                        travel_distances = np.append(travel_distances, travel_distance_per_shower)
                        shower_id = np.append(shower_id, shower)

                    #Shorting the data:
                    time_shower_and_ray, max_amp_shower_and_ray, travel_distances, shower_id = list(zip(*sorted(zip(time_shower_and_ray, max_amp_shower_and_ray, travel_distances, shower_id))))
                    time_shower_and_ray = np.array(time_shower_and_ray)
                    max_amp_shower_and_ray = np.array(max_amp_shower_and_ray)
                    travel_distances = np.array(travel_distances)
                    shower_id = np.array(shower_id)

    
        except:
            pass
    return time_shower_and_ray, max_amp_shower_and_ray, travel_distances, shower_id 

def get_station_data(PATH, station_ID, channel_ID):
    if os.path.isfile(PATH):
        time_shower_and_ray, max_amp_shower_and_ray, travel_distances, shower_id = get_station_data_from_file(PATH, station_ID, channel_ID)
        print("Finish!")
        return time_shower_and_ray, max_amp_shower_and_ray, travel_distances, shower_id 

    elif os.path.isdir(PATH):
        time_shower_and_ray, max_amp_shower_and_ray, travel_distances, shower_id = np.array([]), np.array([]), np.array([]), np.array([])

        for filename in os.listdir(PATH):
            file = os.path.join(PATH, filename)
            a, b, c, d = get_station_data_from_file(file, station_ID, channel_ID)
            time_shower_and_ray = np.append(time_shower_and_ray, a)
            max_amp_shower_and_ray = np.append(max_amp_shower_and_ray, b)
            travel_distances = np.append(travel_distances, c)
            shower_id = np.append(shower_id, d)
        print("Finish!")
        return time_shower_and_ray, max_amp_shower_and_ray, travel_distances, shower_id
    else:
        print("No file or folder found, Please try a different path")

#########################################################################################
#Get data from event file for visualization:
def get_event_data_from_file(PATH):
    f = h5py.File(PATH, mode = "r")
    flavors = np.array(f["flavors"])
    xx = np.array(f["xx"])
    yy = np.array(f["yy"])
    zz = np.array(f["zz"])
    return xx, yy, zz, flavors

def get_event_data(PATH):
    if os.path.isfile(PATH):
        xx, yy, zz, flavors = get_event_data_from_file(PATH)
        return xx, yy, zz, flavors
    
    elif os.path.isdir(PATH):
        flavors = np.array([])
        xx = np.array([])
        yy = np.array([])
        zz = np.array([])
        for filename in os.listdir(PATH):
            file = os.path.join(PATH, filename)
            f = h5py.File(file, mode = "r")
            flavors = np.append(flavors, np.array(f["flavors"]))
            xx = np.append(xx, np.array(f["xx"]))
            yy = np.append(yy, np.array(f["yy"]))
            zz = np.append(zz, np.array(f["zz"]))
        return xx, yy, zz, flavors
    
    else:
        print("No file or folder found, Please try a different path")

########################################################################################
import NuRadioReco.detector.detector as detector
import NuRadioReco.modules.io.eventReader
import astropy
from NuRadioReco.framework.parameters import stationParameters as stnp

#Get the station and shower data from the Nur file
def get_nur_data_from_file(inputfilename, detectordescription, station_ID, channel_ID, e_id):
    # read in detector positions 
    print("Read file:", inputfilename)
    det = detector.Detector(json_filename=detectordescription)

        # initialize modules
    eventReader = NuRadioReco.modules.io.eventReader.eventReader()
    eventReader.begin(inputfilename)
    detector_time = astropy.time.Time('2040-01-01 20:00:00')
    det.update(detector_time)
    trace = np.array([])
    times = np.array([])
    S_ID = np.array([])
    C_ID = np.array([])
    event_id = np.array([])
    for event in eventReader.run():
        print("Get data for event: ",event)
        e_id += 1
        for station in event.get_stations():
            station_id = station.get_id()
            if station_id in station_ID:
                for channel in station.iter_channels():
                    channel_id = channel.get_id()
                
                    if channel_id in channel_ID:
                        # get time trace and times of bins
                        trace_i = channel.get_trace()
                        times_i = channel.get_times()
                        trace_i[np.isnan(trace_i)] = 0 #Remove nan value
                        #Save the data:
                        trace = np.append(trace, trace_i)
                        times = np.append(times, times_i)
                        num_data = len(trace_i)
                        S_ID = np.append(S_ID, np.ones(num_data)*station_id)
                        C_ID = np.append(C_ID, np.ones(num_data)*channel_id)
                        event_id = np.append(event_id, np.ones(num_data)*e_id)
                    else:
                        pass
            else:
                pass   
    return trace, times, S_ID, C_ID, event_id

def get_nur_data(PATH, detectordescription, station_ID, channel_ID): 
    if os.path.isfile(PATH):
        e_id = 0
        trace, times, S_ID, C_ID, event_id = get_nur_data_from_file(PATH, detectordescription, station_ID, channel_ID, e_id)
        print("Number of events:", max(event_id))
        return trace, times, S_ID, C_ID, event_id
    
    elif os.path.isdir(PATH):
        e_id = 0
        trace = np.array([])
        times = np.array([])
        S_ID = np.array([])
        C_ID = np.array([])
        event_id = np.array([])
        for filename in os.listdir(PATH):
            file = os.path.join(PATH, filename)
            v, t, s, c, e = get_nur_data_from_file(file, detectordescription, station_ID, channel_ID, e_id)
            trace = np.append(trace, np.array(v))	
            times = np.append(times, np.array(t))
            S_ID = np.append(S_ID, np.array(s))
            C_ID = np.append(C_ID, np.array(c))
            event_id = np.append(event_id, np.array(e))
            e_id += 1 
        print("Number of events:", e_id)
        return trace, times, S_ID, C_ID, event_id
    
    else:
        print("No file or folder found, Please try a different path")

