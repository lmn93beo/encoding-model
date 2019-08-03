"""
Analyze one-stim data by Rafiq
"""

import numpy as np
import matplotlib.pyplot as plt
import mat4py
import data_classes
import os
import neuron_group_operations

if os.environ['COMPUTERNAME'] == 'DESKTOP-FN1P6HD':
    filedir = 'C:/Users/Sur lab/Documents/RafiqAnalysis'
elif os.environ['COMPUTERNAME'] == 'homecomputer':
    filedir = 'C:/Users/Le/Dropbox (MIT)/Sur/For Nhat'
raw_i = mat4py.loadmat(filedir + '/raw_i.mat')['raw_i']
raw_c = mat4py.loadmat(filedir + '/raw_c.mat')['raw_c']

exp = data_classes.Experiment(dict(l_trials=[0], r_trials=[1]))
neuron_group_c = neuron_group_operations.make_neuron_group(raw_c, exp)
neuron_group_i = neuron_group_operations.make_neuron_group(raw_i, exp)

session_lst = neuron_group_c.get_session_list()
subgroups_c = []
subgroups_i = []
for i in range(np.max(session_lst)):
    neuron_lst = np.where(session_lst == i)[0]
    subgroup_c = neuron_group_c.make_subgroup_by_neurons(neuron_lst)
    subgroup_i = neuron_group_i.make_subgroup_by_neurons(neuron_lst)
    subgroups_c.append(subgroup_c)
    subgroups_i.append(subgroup_i)

# Plot only for individual sessions
def plot_sessions():
    sessions = neuron_group.get_session_list()
    for i in range(np.max(sessions)):
        to_plot_list = np.where(sessions == i)
        mean_act = neuron_group.plot_all_means(sort=True, plotid=to_plot_list, style='heatmap')
        #plt.close('all')

# Compare neurons
def plot_neurons():
    session_id = 0
    for i in range(subgroups_i[session_id].n_neurons):
        plt.figure()
        plt.subplot('121')
        subgroups_c[0].neurons[i].plot_all_trials()
        plt.subplot('122')
        subgroups_i[0].neurons[i].plot_all_trials()

# Decide if a neuron is c or i-responsive
session_id = 3
for i in range(100):
    plt.subplot(10, 10, i + 1)
    i_trials = subgroups_i[session_id].neurons[i]
    c_trials = subgroups_c[session_id].neurons[i]
    plt.errorbar(np.arange(len(c_trials.mean_activity)), c_trials.mean_activity, c_trials.stderr_activity)
    plt.errorbar(np.arange(len(i_trials.mean_activity)), i_trials.mean_activity, i_trials.stderr_activity)
    plt.title(str(i + 1))




