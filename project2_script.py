# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 10:27:34 2023

@author: Jackson Swalley, Thomas Bausman 

project_2 script uses the average normal heart beat and average arrhythmic beat 
to create normalized templates and find the cross correlation between the templates and 
a raw ECG signal. The cross correlation is then used to predict if a beat is normal or arrhythmic given 
a threshold for the values. project2_script also investigates the difference between a raw signal and the 
decimated signal via visualization. 

"""

# import modules 
import numpy as np
from matplotlib import pyplot as plt
import project2_module as p2m
import project1_module as p1m

#%% Cell 1

# define file name for average beats
input_file = 'ecg_means_e0103.npz'

# call function to return variables in file and print fields 
symbols, trial_time, mean_trial_signal = p2m.load_means(input_file)

# define filename for new raw data 
new_file = 'ecg_e0103_half2.npz'

# call function to return variables, print fields of filename 
ecg_voltage, fs, label_samples, label_symbols, subject_id, electrode, units  = p1m.load_data(new_file) 

#%% Cell 2

# calculate sample rate from ecg_means_e0103.npz
sample_rate = len(mean_trial_signal[0]) # sample rate is the same for 1st and second row 

# determine which file was sampled faster 
if sample_rate > fs:
    # assign original_signal to file with faster sample rate, calculate respective decimation factor
    original_signal = mean_trial_signal
    decimation_factor = sample_rate / fs
    # assign frequency and calculate step size 
    frequency = sample_rate
    step = 1/frequency
else:
    # assign original_signal to file with faster sample rate, calculate respective decimation factor
    original_signal = ecg_voltage
    decimation_factor = fs / sample_rate
    # assign frequncy and calculate step size
    frequency = fs
    step = 1/frequency

# call function to decimate signal 
decimated_signal = p2m.decimate_data(original_signal, decimation_factor)

# define start time for time vectors 
start = 0

# time vector for original_signal where step size is variable "step"
stop_1 = len(original_signal)*step
time_vector_1 = np.arange(start, stop_1, step)

# time vector for decimated_signal where step size is "step * decimation_factor"
stop_2 = len(decimated_signal)*step*decimation_factor
time_vector_2 = np.arange(start, stop_2, step*decimation_factor)

# plot original and decimated signal, label, title and legend
plt.figure(1, clear=True, dpi=200)
plt.plot(time_vector_1, original_signal, label = 'original_signal')
plt.plot(time_vector_2, decimated_signal, label = 'decimated_signal', linestyle='dashed' )
plt.xlabel('time (seconds)')
plt.ylabel(f'{units}')
plt.legend()
plt.title('Decimated Signal versus Original Signal')

# modify x axis limits to show difference bwtween signals (specific to this data)
plt.xlim(0.5, 1)
plt.ylim(-0.12,0.07)

# save figure
plt.savefig('Decimated Signal versus Original Signal.pdf' )

# comment on assumptions of signals and implications if assumptions are not valid
print('1.) The decimation factor must change the sampling frequncy so the nyquist theorm still holds: (sampling frequency > 2* highest frequncy in signal)') 
print('2.) This is true for standard clinical ECG data, as the nyquist theorem must be true for all signals to prevent aliasing.')
print('3.) If this assumption is not true, the ECG signal would suffer from aliasing and a patients ECG data could be mis-interpreted which may have life threatening consequences.')
print(' ')
#%% Cell 3

# normalize signal to make normal beat template 
trial_mean = mean_trial_signal[0,:]
template = p2m.normalize_template(trial_mean)

# print units of template based on manipulations done to original signal to get normalized signal 
print(f'The units of the template are 1/{units}')

#%% Cell 4
# calculate template match
template_match = p2m.get_template_match(decimated_signal, template)

# plot template match on seperate figure
plt.figure(2, clear=True, figsize = (10,8), dpi=200)
plt.plot(time_vector_2, template_match, label = 'N template match (A.U.)')
plt.xlabel('time (s)')
plt.ylabel('See legend for details')

#%% Cell 5

# define threshold, call function to get beats above threshold 
threshold = 0.9
beat_samples = p2m.predict_beat_time(template_match, threshold)

# find time and voltage values where an event occurs 
beat_samples_voltage = template_match[beat_samples]
beat_samples_time = time_vector_2[beat_samples]

# plot event predictions as triangles on template 
plt.scatter(beat_samples_time, beat_samples_voltage, marker ="v", color = 'purple', label = 'N prediction')

# set x limits 
plt.xlim(2904,2910)

#%% Cell 6

# run arrhythmia detection by passing function arithmetic trial mean 
trial_mean = mean_trial_signal[1,:]
beat_samples, template_match = p2m.run_beat_detection(trial_mean, decimated_signal, threshold)

# get arrhythmic beat times and voltages 
a_beat_voltage = template_match[beat_samples]
a_beat_time = time_vector_2[beat_samples]

# plot arrhythmic predictions 
plt.scatter(a_beat_time, a_beat_voltage, marker ='v', color ='red', label = 'A prediction')

# plot labels of confirmed types of 
p1m.plot_events(label_samples, label_symbols, time_vector_1, ecg_voltage)

# plot arrhythmic template match on top of normal template match 
plt.plot(time_vector_2, template_match, label = 'A template match (A.U.)')

# plot original signal 
plt.plot(time_vector_1, ecg_voltage, linestyle='dashed', linewidth = 1, label = f'original signal ({units})')

# make legend and title
plt.title('Decimated Signal with Normal and Arrhythmic Templates labeled with Predictions and Confirmation of Beat Types.')
plt.legend() 
plt.tight_layout()

# save figure 
plt.savefig('Decimated Signal with Normal and Arrhythmic Templates labeled with Predictions and Confirmation of Beat Types.pdf')
