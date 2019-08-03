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
raw_data = mat4py.loadmat(filedir + '/error_trial_data_20190715.mat')
c_corr_all = raw_data['c_corr_all']
c_incorr_all = raw_data['c_incorr_all']
i_corr_all = raw_data['i_corr_all']
i_incorr_all = raw_data['i_incorr_all']

# Create neurons
cell_lst_c_corr = []
cell_lst_c_incorr = []
cell_lst_i_corr = []
cell_lst_i_incorr = []

for i in range(len(c_corr_all)):
    cell_arr_c_corr = np.array(c_corr_all[i])
    cell_arr_c_incorr = np.array(c_incorr_all[i])
    cell_arr_i_corr = np.array(i_corr_all[i])
    cell_arr_i_incorr = np.array(i_incorr_all[i])

    ntrials_c_corr = cell_arr_c_corr.shape[1]
    ntrials_c_incorr = cell_arr_c_incorr.shape[1]
    ntrials_i_corr = cell_arr_i_corr.shape[1]
    ntrials_i_incorr = cell_arr_c_incorr.shape[1]

    exp = data_classes.Experiment(dict(l_trials=[0], r_trials=[1]))
    neuron_c_corr = data_classes.OneStimNeuron(i, cell_arr_c_corr, exp)
    neuron_c_incorr = data_classes.OneStimNeuron(i, cell_arr_c_incorr, exp)
    neuron_i_corr = data_classes.OneStimNeuron(i, cell_arr_i_corr, exp)
    neuron_i_incorr = data_classes.OneStimNeuron(i, cell_arr_i_incorr, exp)

    cell_lst_c_corr.append(neuron_c_corr)
    cell_lst_c_incorr.append(neuron_c_incorr)
    cell_lst_i_corr.append(neuron_i_corr)
    cell_lst_i_incorr.append(neuron_i_incorr)

# Make neuron groups
neur_group_c_corr = neuron_group_operations.NeuronGroup(cell_lst_c_corr)
neur_group_c_incorr = neuron_group_operations.NeuronGroup(cell_lst_c_incorr)
neur_group_i_corr = neuron_group_operations.NeuronGroup(cell_lst_i_corr)
neur_group_i_incorr = neuron_group_operations.NeuronGroup(cell_lst_i_incorr)

# Plot
neur_group_c_corr.plot_all_means(sort=True, style='heatmap')
neur_group_c_incorr.plot_all_means(sort=True, style='heatmap')
neur_group_i_corr.plot_all_means(sort=True, style='heatmap')
neur_group_i_incorr.plot_all_means(sort=True, style='heatmap')

'''

# Create neurons
cell_lst = []
for i in range(len(raw_i)):
    cell_arr_i = np.array(raw_i[i])
    cell_arr_c = np.array(raw_c[i])

    ntrials_i = cell_arr_i.shape[1]
    ntrials_c = cell_arr_c.shape[1]

    # Make corresponding experiment object
    exp = data_classes.Experiment(dict(l_trials=np.arange(ntrials_i),
                                       r_trials=np.arange(ntrials_i, ntrials_i + ntrials_c)))

    cell_arr = np.hstack((cell_arr_i, cell_arr_c)).T
    neuron = data_classes.OneStimNeuron(i, cell_arr, exp)
    cell_lst.append(neuron)

# Assign session to neurons
curr_session = -1
curr_ntrials = -1
for neuron in cell_lst:
    if neuron.ntrials != curr_ntrials:
        curr_session += 1
        print('Changed, current ntrials = ', neuron.ntrials)
        print('Current session = ', curr_session)
        curr_ntrials = neuron.ntrials
    neuron.session = curr_session

neuron_group = neuron_group_operations.NeuronGroup(cell_lst)

# Make subgroups based on sessions
session_lst = neuron_group.get_session_list()
subgroups_by_session = []
for i in range(np.max(session_lst)):
    neuron_lst = np.where(session_lst == i)[0]
    subgroup = neuron_group.make_subgroup_by_neurons(neuron_lst)
    subgroups_by_session.append(subgroup)

# Split each group based on left or right
subgroups_left = []
subgroups_right = []
for subgroup in subgroups_by_session:
    left_trials = subgroup.exp.l_trials
    right_trials = subgroup.exp.r_trials
    subgroup_left = subgroup.make_subgroup_by_trials(left_trials)
    subgroup_right = subgroup.make_subgroup_by_trials(right_trials)
    subgroups_left.append(subgroup_left)
    subgroups_right.append(subgroup_right)
    #subgroup.plot_all_means(sort=True, style='heatmap')

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
    for i in range(3): #(subgroups_right[session_id].n_neurons):
        plt.figure()
        plt.subplot('121')
        subgroups_left[0].neurons[i].plot_all_trials()
        plt.subplot('122')
        subgroups_right[0].neurons[i].plot_all_trials()

# Decide if a neuron is left or right-responsive
session_id = 3
for i in range(100):
    plt.subplot(10, 10, i + 1)
    left_trials = subgroups_left[session_id].neurons[i]
    right_trials = subgroups_right[session_id].neurons[i]
    plt.errorbar(np.arange(len(left_trials.mean_activity)), left_trials.mean_activity, left_trials.stderr_activity)
    plt.errorbar(np.arange(len(right_trials.mean_activity)), right_trials.mean_activity, right_trials.stderr_activity)
    plt.title(str(i + 1))



'''
