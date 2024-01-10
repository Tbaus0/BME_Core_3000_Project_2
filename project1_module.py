# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 09:59:48 2023

File Name: project1_module.py

project1_module.py has the functions gather_data(), plot_raw_data(), plot_events(), extract_trials(), plot_mean_and_std_trials(), 
and save_means() to manipulate ecg data.

@author: Nathan and Jackson
"""
import numpy as np
from matplotlib import pyplot as plt

#%% load_data Function 
def load_data(filename):
    """
    load_data takes the name of a .npz file and extracts the data field names. The feild names are used to 
    index the data file and load values into seperate arrays which are returned.

    Parameters
    ----------
    filename : string
        used to load external data within function.

    Returns
    -------
    ecg_voltage : 1D array of floats, size (n, 1) ; where n is every voltage value collected 
        allows for access to ecg voltage entries outside of function.
        
    fs : 1D array of int (one entry)
        used for calculation of ecg voltage sampling interval outside of function.
        
    label_samples : 1D array of ints, size (a, 1); where a is total number of events.
        voltage index where event occurs.
        
    label_symbols : 1D array of strings, size (a, 1); where a total number of events.
        "N" or "V" entries that is the same length as label_samples.
        
    subject_id : 1D array of string (one entry)
        identification label coresponding to patient ECG data.
        
    electrode : 1D array of string (one entry)
        identification label coresponding to electrode used to collect data.
        
    units_array : 1D array of string (one entry)
        units of y-axis for voltage plots.

    """
    
    # load file
    data = np.load(filename)
    
    # extract fields (specific to .npz)
    data_fields = data.files
    
    # use fields to extract variables
    ecg_voltage = data[data_fields[0]]
    fs = data[data_fields[1]]
    label_samples = data[data_fields[2]]
    label_symbols = data[data_fields[3]]
    subject_id = data[data_fields[4]]
    electrode = data[data_fields[5]]
    units = data[data_fields[6]]
    
    # print field names
    print(f'Fields in {filename}:')
    for value in data_fields:
        print(value)
    
    # Returns all the data arrays
    return ecg_voltage, fs, label_samples, label_symbols, subject_id, electrode, units


#%% plot_raw_data Function 
def plot_raw_data(signal_time, signal_voltage, units = 'V', title = ''):
    """
   plot_raw_data plots signal_voltage values on a time interval containing two normal beats and one 
   arrhythmia as a line plot with appropriately labeled axes and a descriptive title. 
   
   Parameters
   ----------
   signal_voltage : 1D array of floats, size (n , 1) ; where n is every collected voltage
       plotted on y-axis of figure 1.
       
   signal_time : 1D array of floats, size (n , 1); where n is time when voltage is sampled
       plotted on x-axis of figure 1.
       
   units : string, optional
       units of y-axis on plot. The default is 'V' for volts.
       
   title : string, optional
       descriptive title of line plot. The default is ''.

   Returns
   -------
   None.

    """
    # creating figure and plot raw ecg signal against time
    plt.figure(1, clear = True, dpi=200)
    plt.plot(signal_time, signal_voltage, label = 'signal')

    # annotate plot with correct units
    plt.xlabel('time (s)')
    plt.ylabel(f'voltage ({units})')
    plt.grid()
    plt.title(title)
    
#%% plot_events Function 
def plot_events(label_samples, label_symbols, signal_time, signal_voltage):
    """
    plot_events labels where each event ocurrs on the figure made in plot_raw_data. Different 
    types of events are distinguished by different colored dots and labeled in the legend. 
    
    Parameters
    ----------
    label_samples : 1D array of ints, size (a, 1); where a is total number of events.
        voltage index where event occurs.
        
    label_symbols : 1D array of strings, size (a, 1); where a total number of events.
        "N" or "V" entries corresponding to normal (N) or arrhythmic (V) beats.
        
    signal_time : 1D array of floats, size (n , 1); where n is time when voltage is sampled
        plotted on x-axis of figure 1..
        
    signal_voltage : 1D array of floats, size (n , 1) ; where n is every collected voltage
        plotted on y-axis of figure 1.

    Returns
    -------
    None.

    """
    
    # creates array of unique values in label_symbols
    unique_symbols = np.unique(label_symbols)
   
    # creates array of voltage idicies where event of unique label "N" occurs 
    for sample_index in label_samples:
        #reads tuple into array
        n_label_index = label_samples[label_symbols == unique_symbols[0]] 
        
    # creates array of voltage idicies where event of unique label "N" occurs 
    for sample_index in label_samples:
        #reads tuple into array
        v_label_index = label_samples[label_symbols == unique_symbols[1]] 
        
    #find times and voltages of N and V
    
    # fidning n
    n_times = signal_time[n_label_index]
    n_voltages = signal_voltage[n_label_index]
   
    # finding v
    v_times = signal_time[v_label_index]
    v_voltages = signal_voltage[v_label_index]
    
    # plotting data 
    plt.scatter(n_times, n_voltages, c = "green", marker = 'o', label = 'N')
    plt.scatter(v_times, v_voltages, c = "orange", marker = 'o', label = 'A')
    
    # creating legend for the figure 1 plot
    plt.legend()
   
    # constants used to find 
    start_xlim = 1913
    stop_xlim = 1916
    plt.xlim(start_xlim, stop_xlim)
    plt.tight_layout()

#%% extract_trials Function 
def extract_trials(signal_voltage, trial_start_samples, trial_sample_count):
    """
    extract_trials takes in the ecg voltage, the time start of every trial, and the number of samples in each trial
    and finds all the voltage data assocaited with one trial and stores it all in a 2D array that is returned 

    Parameters
    ----------
    signal_voltage : 1D array of floats, size (n, 1); n being the number of voltage data points collected 
        plotted on y-axis of figure 1

    trial_start_samples : 1D array of size (X, 1); X being the number of ecg event for a unique type of event
        this 1D array holds the start index for every individual ecg event of that unique event type

    trial_sample_count : int 
        represnts the number of samples for each trial

    Returns
    -------
    trials : 2D array of size (X by trial_sample_count)
        2D array holds all the ecg voltage values for each indiviudal trial of a unique event type


    """
    # get number of rows in trial from length of signal_voltage and trial_sample_count; convert to integers
    trial_count = int(np.size(trial_start_samples)) # number of total trials
    # makes any array of zeros from calculated values above
    trials = np.zeros((trial_count, int(trial_sample_count)))
    # fill the array with nan's 
    trials.fill(np.nan)
    # fills array of zeros with trial voltages 
    for event_index in range(0, trial_count):
        # try to add to trials, because the first and last index would not work
        try:
            # adds the 250 indeces from the trial start and on
            trials[event_index, :] = signal_voltage[int(trial_start_samples[event_index]): int(trial_start_samples[event_index]) + trial_sample_count]
        except:
            # if they dont work which are the first index and last index of trial_start_samples
            # passes to keep them as nan in trials 
            pass
    # returns filled 2D array 
    return trials

#%% plot_mean_and_std_trials Function 
def plot_mean_and_std_trials(signal_voltage, label_samples, label_symbols, trial_duration_seconds, fs, units = "", title = ""):
    """
    plot_mean_and_std_trials takes in voltages, the sampels and symbols of events, the time of one sample, repesentation of a sample size, 
    loops through all the data to find the mean and standard deviation for each unique event over trial_druation_seconds, plots them in one figure, 
    and then saves and returns the unique symbols, trial_times, and means for each unique event
   
    Parameters
    ----------
    signal_voltage : 1D array of floats, size (n , 1); n being the number of voltage data points collected 
        plotted on y-axis of figure 1

    label_samples : 1D array of ints, size (a, 1); a being the number of ecg events that occured 
        represents the index values of each event 

    label_symbols : 1D array of characters, size (a, 1); a being the number of ecg events that occured 
        character represnetation of the type of event for each event

    trial_duration_seconds : int 
        represnets the duration of one trial in seconds 

    fs : int 
        the sampeling frequency thrughout all the collected data

    units : string, optional
        Represents the units of voltage on the graph. The default is "".

    title : string, optional
        Represents the title of the graph. The default is "".

 
    Returns
    -------
    symbols : 1D array (w, 1) w being the number of unique symbols in label_symbols
        array that represents the number of unique ecg events 

    trial_time : 1D array of size (fs, 1) fs being the sample frequency
        1D array represents the time values for a trial

    mean_trial_signal : 2D array of size (symbols, fs) symbols holding the number of unique events, and fs being the sample frequency
        2D array hold the mean voltages for each time index for each unique event type

    """
    
    
    # creating the signal step, which is the size of each index 
    signal_time_step = trial_duration_seconds/fs 
    
    # creating a count which is a constant variable that represnts the number of index for each event
    trial_sample_count = int(trial_duration_seconds / signal_time_step)
    
    # time array that is the size of an event duration with the times
    trial_time = np.arange(0, trial_duration_seconds, signal_time_step)
    
    # creating an array of all the unique symbols in label_symbols
    symbols = np.unique(label_symbols)
    
    # creating an empty array of zeros that will hold the means for both Normal and Abnormal events
    mean_trial_signal = np.zeros((len(symbols), trial_sample_count))
    
    # creating a figure to plot mean and standard deviation on
    plt.figure(3, clear = True, dpi=200)
    
    # for loop to loop through the array of unique symbols
    for symbol_index in range(len(symbols)):
        
        # symbol is a variable representing 
        symbol = symbols[symbol_index]
        trials = label_samples[label_symbols == symbol]
        
        # change V to A for "Arrhythmia"
        if symbol == "V":
            symbol = "A"
        
        # create trial start index array 
        trial_start_time = trials - (trial_sample_count/2)
        
        # call function that extracts trials for unique symbol, store as variable, calculate mean voltages in trials 
        trial_event_array = extract_trials(signal_voltage, trial_start_time, trial_sample_count)
        trial_mean = np.nanmean(trial_event_array, axis=0)
        
        # Yo Nate IDK what this one does lowkey 
        mean_trial_signal[symbol_index, :] = trial_mean
        trial_standard_deviation = np.nanstd(trial_event_array, axis=0)
        
        # plot voltages in trials for specific labels centered event where event time = 0 
        plt.plot(trial_time - 0.5, trial_mean, label = symbol, lw = 2.5)
        plt.fill_between(trial_time - 0.5, (trial_mean + trial_standard_deviation), (trial_mean - trial_standard_deviation), alpha = 0.3, label = f"standard deviation {symbol}")

    # add title to plot, label x and y axis, create legend 
    plt.title(title)
    plt.xlabel('trial_duration_seconds centered at event (s)')
    plt.ylabel(f'voltage ({units})')
    plt.tight_layout()
    plt.legend()
    

    return symbols, trial_time, mean_trial_signal

#%% save_means Function
def save_means(symbols, trial_time, mean_trial_signal, out_filename = 'ecg_means.npz'):
    """
    save_means takes in necessary data and saves it to file for latter use
    
    Parameters
    ----------
    symbols : 1D array of strings the size of the number of unique event symbols
        symbols array is filled by all the unique strings that represent a type of event in the ecg signal

    trial_time : 1D array of floats size sample size count for one event trial
        trial_time is a 1D array that represents all the time indexes for one ecg event, being the size of one the sample size for one event

    mean_trial_signal : 2D aray of floats size of len(symbols) by trial_sample_count 
        mean_trial_signal is a 2D array storing the mean ecg voltage value for each unique event for the entirety of one trial 

    out_filename : string 
        out_filename stores the stirng representation of the name of the file where data is being stored

 
    Returns
    -------
    None.

    """
    np.savez(out_filename, symbols = symbols, trial_time = trial_time, mean_trial_signal = mean_trial_signal)
    
    
    
    
    
    
    
    
    