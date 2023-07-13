from Get_hdf5_data import get_attr_data, get_shower_data, get_station_data, get_Veff_data, get_event_data, get_nur_data
from Get_hdf5_data import get_nur_data_from_file
from scipy.signal import find_peaks
import numpy as np
import time as TIME
import datetime
import math
import os 
import astropy
import NuRadioReco.detector.detector as detector
import NuRadioReco.modules.io.eventReader
import csv
import pandas as pd
from math import isclose
from matplotlib import pyplot as plt


###############################################
#Calculate root mean square of an array
def rms(v): 
    return np.sqrt(np.mean(v**2))
##########################################################################################
#Save data:
def if_else_write(i, row, column):
    if i < len(column):
        row.append(column[i])
    else:
        row.append('')

#Save data into csv file:
def save_data_2_csv(file_path ,inputfilename, detectordescription, station_ID, channel_ID):
    # Data to be written into columns
    if os.path.isfile(inputfilename):
        trace, times, S_ID, C_ID, event_id = get_nur_data_from_file(inputfilename, detectordescription, station_ID, channel_ID, 0)
        # Determine the maximum column length
        max_length = max(len(trace), len(times), len(S_ID), len(C_ID), len(event_id))
        print("Max length of data:", max_length)
        print("Start writing data into csv file")
        # Write data into columns
        with open(file_path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Traces', 'Times', 'station', 'channel', 'event'])  # Write header row
            # Write data rows
            for i in range(max_length):
                row = []
                if_else_write(i, row, trace) 
                if_else_write(i, row, times)
                if_else_write(i, row, S_ID)
                if_else_write(i, row, C_ID)
                if_else_write(i, row, event_id)

                writer.writerow(row)
        print("Finish writing data into csv file")
    elif os.path.isdir(inputfilename):
        isExist = os.path.exists(file_path)
        if not isExist:
            with open(file_path, 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(['Traces', 'Times', 'station', 'channel', 'event'])
        with open(file_path, 'a', newline='') as file:
            writer = csv.writer(file)
            e_id = 0
            #########################################
            for filename in os.listdir(inputfilename):
                #reset data:
                event_id = np.array([])
                ###########################################
                file = os.path.join(inputfilename, filename)
                trace, times, S_ID, C_ID, e = get_nur_data_from_file(file, detectordescription, station_ID, channel_ID, e_id)
                event_id = np.append(event_id, np.array(e))
                e_id += 1 
                # Determine the maximum column length
                max_length = max(len(trace), len(times), len(S_ID), len(C_ID), len(event_id))
                print(f"Max length of data for event {e_id}:", max_length)
                print("Start writing data into csv file")
                # Write data into columns
                # Write data rows
                for i in range(max_length):
                    row = []
                    if_else_write(i, row, trace) 
                    if_else_write(i, row, times)
                    if_else_write(i, row, S_ID)
                    if_else_write(i, row, C_ID)
                    if_else_write(i, row, event_id)

                    writer.writerow(row)
        print("Finish writing data into csv file")
##########################################################################################
#Read data:
#Read column by name
def read_column_csv(file_path, column_name):
    # Columns to read (by column name)
    columns_to_read = [column_name]  # Example: read columns 'Column1' and 'Column3'

    # Read columns by name
    C = []
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            selected_columns = [row[col_name] for col_name in columns_to_read]
            try:
                C.append(float(str(selected_columns).replace("[", "").replace("]", "").replace("'", "")))
            except:
                pass
    C = np.array(C)
    return C 


#Read row by index
def read_row_csv(file_path, row_index):
    # Rows to read (by row index)
    rows_to_read = [row_index]  # Example: read rows 0, 1, 2, 3, 4, 5
    # Read rows by index
    R = []
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file)
        for row_index, row in enumerate(reader):
            if row_index in rows_to_read:
                R.append(row)
    with open('samplecsv.csv') as file_obj:
        
        # Skips the heading
        # Using next() method
        heading = next(file_obj)
        
        # Create reader object by passing the file 
        # object to reader method
        reader_obj = csv.reader(file_obj)
        
        # Iterate over each row in the csv file 
        # using reader object
        for row in reader_obj:
            print(row)
    return R
#get number of row in csv file
def get_num_row(data):
    with open(data) as f:
        num_row = sum(1 for line in f)
    return num_row

#Get data by row index:
def get_data_by_row_index(data, i): #i is row index
    data_frame = pd.read_csv(data)
    trace = float(data_frame.loc[i]['Traces'])
    times = float(data_frame.loc[i]['Times'])
    s_id = float(data_frame.loc[i]['station'])
    c_id = float(data_frame.loc[i]['channel'])
    e_id = float(data_frame.loc[i]['event'])
    #print(trace, times, s_id, c_id, e_id)
    return trace, times, s_id, c_id, e_id

#Get data by sample size:
def get_data_by_sample_size(data, row_index , sample_size):
    trace, times, s_id, c_id, e_id = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
    for i in range(sample_size):
        v, t, s, c, e = get_data_by_row_index(data, row_index+i)
        trace = np.append(trace, v)
        times = np.append(times, t)
        s_id = np.append(s_id, s)
        c_id = np.append(c_id, c)
        e_id = np.append(e_id, e)
    return trace, times, s_id, c_id, e_id
##########################################################################################
#Analyzing data:
#make array of max value of each bin
def make_max_array(data, num_subarrays): 
    subarrays = np.array_split(data, num_subarrays) #split array into smaller chunk
    max_values = np.array([np.max(subarray) for subarray in subarrays]) #max value of each bin
    return max_values

#Binning data:
def find_index(trace1, V_peak1):
    index = []
    for i in range(len(trace1)):
        if trace1[i] in V_peak1:
            index.append(i)
    return index

def binning_data(v, t, bin_size):
    v = abs(v) #absolute value
    num_bin = int(len(v)/bin_size) #number of bin
    V = make_max_array(v, num_bin) #max value of each bin
    index_max = find_index(v, V) #index of max value of each bin
    T = t[index_max] #time of each bin
    return V, T
##########################################################################################
def get_bin_heights_threshold(v):
    V = np.histogram(v, bins=len(v))
    sigma = np.std(v)
    index = np.where(abs(V[1][:-1]) > 2*sigma)
    H = V[0][index]
    return max(H)

#Get noise from waveform recieceved by each signal:
def get_noise(v, t, bin_size):
    v, t = binning_data(v, t, bin_size)
    v[np.isnan(v)] = 0 #remove nan value
    #v, t = binning_data(v, t, bin_size)
    try:
        num_bins = len(v)
        #bin_heights_threshold = int(0.001*num_bins) #Need to improve the noise calculation********
        bin_heights_threshold = get_bin_heights_threshold(v)
        print("bin_heights_threshold:", bin_heights_threshold)
        h = np.histogram(v, num_bins)
        index = np.where(h[0] > bin_heights_threshold)
        #noise = np.median(h[1][index])
        noise = rms(h[1][index])
        #noise = h[1][2]
    except:
        noise = rms(v)
    t, v = list(zip(*sorted(zip(t, v))))
    t = np.array(t)
    v = np.array(v)
    return noise , v, t

#############################################
#get the number of peak
def get_peak(v, t, noise, threshold):
    #signal to noise ratio:
    snr = v/noise
    #get the index of peak
    peaks, _ = find_peaks(snr, height=threshold)
    num_peak = len(peaks)
    #if num_peak == 0 and max(snr) <= 0.1:
    #    plt.plot(t, snr, 'b')
    #    plt.plot(t, np.ones(len(t))*threshold, 'r')
    #    plt.show()
    if max(snr) > threshold and num_peak == 0:
        num_peak = 1
        peaks = np.where(snr == max(snr))
    v_peak = v[peaks]
    t_peak = t[peaks]
    noise = np.ones(len(v_peak))*noise
    #if max(snr) < 1:
    #    v_peak, t_peak, num_peak, noise = np.array([]), np.array([]), np.array([]), np.array([])
    return v_peak, t_peak, num_peak, noise

#Get noise and peak per event:
def get_noise_and_peak_from_bin(trace, times, bin_size, threshold, NOISE):
    #noise, V, T = get_noise(trace, times, bin_size) #first method
    #v_peak, t_peak, num_peak, noise = get_peak(V, T, noise, threshold)
    #v_peak, t_peak, num_peak, noise = np.array([]), np.array([]), np.array([]), np.array([])
    #NOISE = get_noise(trace) 
    V, T = binning_data(trace, times, bin_size)
    v_peak, t_peak, num_peak, noise = get_peak(V, T, NOISE, threshold)

    return v_peak, t_peak, num_peak, noise
    #Break the array into smaller chunk
""" num_subarrays = int(len(trace)/sample_size)
    subtrace = np.array_split(trace, num_subarrays)
    subtime = np.array_split(times, num_subarrays)
    for i in range(num_subarrays): 
        V, T = binning_data(subtrace[i], subtime[i], bin_size)
        a, b, c, d = get_peak(V, T, NOISE, threshold)
        v_peak = np.append(v_peak, a)
        t_peak = np.append(t_peak, b)
        num_peak = np.append(num_peak, c)
        noise = np.append(noise, d)
"""


#Counting peak from waveform event recieceved by each station and channel:
def count_peak(data, station_ID, channel_ID, bin_size, threshold, sample_size):
    start_time = TIME.time()
    V_peak, T_peak, PEAK, NOISE = np.array([]), np.array([]), np.array([]), np.array([])
    data_frame = pd.read_csv(data)
    num_row = len(data_frame)
    E_id = 0
    #noise, _, _ = get_noise(data_frame['Traces'], data_frame['Times'], bin_size)
    noise = 1.0661514487245143e-05 #2.8e-05
    for i in np.arange(0, num_row, sample_size):
        try: 
            trace = np.array(data_frame.loc[i:i+sample_size]['Traces'])
            times = np.array(data_frame.loc[i:i+sample_size]['Times'])
            s_id = np.array(data_frame.loc[i:i+sample_size]['station'])
            c_id = np.array(data_frame.loc[i:i+sample_size]['channel'])
            e_id = np.array(data_frame.loc[i:i+sample_size]['event'])
            #print(s_id, c_id)
            if any(item in s_id for item in station_ID) and any(item in c_id for item in channel_ID):
                v, t, p, n = get_noise_and_peak_from_bin(trace, times, bin_size, threshold, noise)
                #print("Peak", p,"Channel ID:", min(c_id), "Station ID:", min(s_id), max(v)/noise)
                V_peak = np.append(V_peak, v)
                T_peak = np.append(T_peak, t)
                PEAK = np.append(PEAK, p)
                NOISE = np.append(NOISE, n)
                #print(f"Number of peak per {sample_size*0.625}ns:", p, ",Channel ID:", min(c_id), ",Station ID:", min(s_id))
            if max(e_id) > E_id:
                E_id = max(e_id)
                print("Event ID:", E_id)


        except:
            pass
    print("Finish!")
    end_time = TIME.time()
    time_second = end_time - start_time
    run_time = str(datetime.timedelta(seconds=time_second))
    print(f"Time taken to get data: ",run_time)
    print("Number of peak:", sum(PEAK), "Mean noise:", np.mean(NOISE))
    return V_peak, T_peak, PEAK, NOISE













