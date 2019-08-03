"""
Analyze one-stim data by Rafiq
"""

import numpy as np
import matplotlib.pyplot as plt
import utils
import mat4py
import data_classes
import os
import neuron_group_operations

if os.environ['COMPUTERNAME'] == 'DESKTOP-FN1P6HD':
    filedir = 'C:/Users/Sur lab/Documents/RafiqAnalysis'
elif os.environ['COMPUTERNAME'] == 'homecomputer':
    filedir = 'C:/Users/Le/Dropbox (MIT)/Sur/For Nhat'
raw_data = mat4py.loadmat(filedir + '/error_trial_data_permuted_20190715.mat')
c_corr_all = raw_data['c_corr_all']
c_incorr_all = raw_data['c_incorr_all']
i_corr_all = raw_data['i_corr_all']
i_incorr_all = raw_data['i_incorr_all']

# Make neuron groups
exp = data_classes.Experiment(dict(l_trials=[0], r_trials=[1]))
neur_group_c_corr = neuron_group_operations.make_neuron_group(c_corr_all, exp)
neur_group_i_corr = neuron_group_operations.make_neuron_group(i_corr_all, exp)
neur_group_c_incorr = neuron_group_operations.make_neuron_group(c_incorr_all, exp)
neur_group_i_incorr = neuron_group_operations.make_neuron_group(i_incorr_all, exp)

# Plot
neur_group_c_corr.plot_all_means(sort=True, style='heatmap')
neur_group_c_incorr.plot_all_means(sort=True, style='heatmap')
neur_group_i_corr.plot_all_means(sort=True, style='heatmap')
neur_group_i_incorr.plot_all_means(sort=True, style='heatmap')

# Combine groups
neur_c = neuron_group_operations.combine_groups_by_trials(neur_group_c_corr, neur_group_c_incorr)

#Plotting for validation
utils.plot_example_neurons(neur_group_c_corr, np.arange(100))
utils.plot_example_neurons(neur_c, np.arange(100))
utils.plot_example_neurons(neur_group_c_incorr, np.arange(100))

