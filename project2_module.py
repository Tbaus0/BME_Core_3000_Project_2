# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 10:25:18 2023

@author: Jackson Swalley
"""
# import modules 
import numpy as np

#%% function 1
# create load means function that takes npz file 
def load_means(filename):
    '''
    load_means takes a string that is the name of a .npz file. The function loads the file, extracts and prints 
    the field names of each array of data, extracts the variable values using the field names and returns the 
    arrays. 
    
    Parameters
    ----------
    filename : string
        file name of .npz file that has fields symbols, trial_time, and mean_trial_signal.

    Returns
    -------
    symbols : 1D array of strings size (a,1) where a is the number of unique cardiac events
        contains all unique catagoriztions of cardiac events.
    trial_time : 1D array of floats size (m,1) where m is the length of time of the average beat 
        time vector that matches period and length of mean_trial_signal.
    mean_trial_signal : 2D array of floats size (a,m) where a is the number of unique event types
        contains mean normal beat and mean arrythmic beat; used to make cross corelation templates.

    '''
    # load data
    data = np.load(filename)
    
    # extract fields
    data_fields = data.files
    
   # get variables from fields
    symbols = data[data_fields[0]]
    trial_time = data[data_fields[1]]
    mean_trial_signal = data[data_fields[2]]
    
    # print field names
    print(f'Fields in {filename}:')
    for value in data_fields:
        print(value)
    
    print(' ')
    
    # return values 
    return symbols, trial_time, mean_trial_signal



#%% function 2
def decimate_data(original_signal, decimation_factor):
    '''
    decimate_data takes a signal and resamples it less frequently given an integer decimation factor.

    Parameters
    ----------
    original_signal : 1D array of float size (n,1) where n is the number of voltage values collected 
        used to create new decimated signal.
    decimation_factor : integer 
        used to resample data less frequently without losing information .

    Returns
    -------
    decimated_signal : 1D array of float size (n/decimation_factor,1) where n is the number of voltage values collected in the original signal 
        resampled data now matches the sampling frequency of the mean_trial_signal and cross correlation can be obtained.

    '''
    # define start and stop index of signal using decimation_factor to make index for decimated signal
    start = 0
    stop = len(original_signal)
    original_signal_decimated_index = np.arange(start, stop, decimation_factor).astype(np.int)
    
    # create decimated signal 
    decimated_signal = original_signal[original_signal_decimated_index]
    
    # return decimated signal
    return decimated_signal
    
#%% function 3
def normalize_template(trial_mean):
    '''
    normalize_template takes the trial mean of a heart beat, subtracts the average value in the signal 
    and divides the signal values by the energy of the signal. The template is later used to find the
    cross correlation between original signal and template

    Parameters
    ----------
    trial_mean : 1D array of float size (v,1) where v is the number of samples in an average beat given frequency
        used to create a normalized template.

    Returns
    -------
    template : 1D array of float size (v,1) where v is the number of samples in an average beat
        later used to find the cross correlation of original signal with template.

    '''
    # calculate average value in array
    average = np.mean(trial_mean, 0)
    
    # initialize enegery variable
    energy = 0
    for value in trial_mean:
        # sum the values in trial_mean according to given energy formula 
        energy += (np.absolute(value))**2
        
    # create template by subtracting average and normalizing by energy
    template = (trial_mean - average) / energy

    # return template
    return template

#%% function 4
def get_template_match(signal_voltage, template):
    '''
    This function is used to create an array with the length of the signal_voltage 
    that is the cross-correlation of the signal voltage and template inputs. The array
    will give you how similar the signal_voltage is to the template. Smaller numbers resemble
    less correlation and higher numbers resemble more similarities. 

    Parameters
    ----------
    signal_voltage : 1D array of floats any length (n,1) where n is the number of voltage data points collected 
        This array should be data from an ECG that you want to test for similarities with a template
    template : 1D Array of float any length (n,1)
        This array should be the reference you are comparing signal_voltage to in order to find cross correlation 

    Returns
    -------
    template_match : 1D array of floats with length signal_voltage
        The template match is the cross-correlation of signal_voltage and template
    '''
    # flip template according to cross corelation formula
    template = np.flip(template)
    
    # compute cross correlation (convolution)
    template_match = np.convolve(signal_voltage, template, mode = 'same')
    
    return template_match

#%% function 5
def predict_beat_time(template_match, threshold):
    """
    This function runs through the template match and picks out the first sample
    of every time the threshold is passed. It returns the sample index where this occurs. 
    

    Parameters
    ----------
    template_match : 1D array of floats length (n,1) where n is the number of voltage data points collected 
        The cross correlation between signal and template in a 1D array
    threshold : float
        a single float that represents the threshold that template_match will pass when a beat is recognized 

    Returns
    -------
    beat_samples : 1D array of integers length (m,1) where m is the number of samples in template_match just above the set threshold 
        each number in beat_samples represents the sample where the threshold given is surpassed
        for the first time of each beat.

    """

    # get array of values above thresholds
    if_threshold = template_match >= threshold
    
    # define stop time
    stop = len(template_match)
    
    # multiply array of ones by values above threshold to make all other values 0
    ones = np.ones(stop) * if_threshold
    
    # manipulate array to have a 1 at the first value when array crosses threshold
    diff = np.diff(ones)
    
    # make index of first events above threshold 
    beat_samples = np.array(np.where(diff==1))
    
    # return the index of events just above threshold 
    return beat_samples

#%% function 6
def run_beat_detection(trial_mean, signal_voltage, threshold):
    """
    This function acts as a parent function to run three consecutive function within the module. 
    With the inputs it will get a template from the trial_mean input and cross correlate it
    to the signal voltage at a given threshold. This will return values where the signal_voltage
    matches the trial_mean template.

    Parameters
    ----------
    trial_mean : 1D array of floats
        The array contains value to make a template out of
    signal_voltage : 1D array if floats
        The signal that is going to be compared to the trial_mean template
    threshold : float
        a single float that represents the threshold that template_match will pass when a beat is recognized 


    Returns
    -------
    beat_samples : 1D array of int (m,1) where m is the number of samples in template_match just above the set threshold 
        all integers represent samples that are the first to go above the threshold for each beat
    template_match : 1D array of floats
        this array shows the cross-correlation between the signals and trial_mean template
    """

    # get template from normalize_template with input trial_mean
    template = normalize_template(trial_mean)
    
    # get template_match from get_template_match with inputs signal_voltage and template
    template_match = get_template_match(signal_voltage, template)
    
    # get beat_samples from predict_beat_time with inputs template_match and threshold
    beat_samples = predict_beat_time(template_match, threshold)

    # return beat_samples and template_match 
    return beat_samples, template_match


    
    