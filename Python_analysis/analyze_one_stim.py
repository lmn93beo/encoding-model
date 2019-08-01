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

# Parse into numpy arrays
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
neuron_group.classify_all()
neuron_group.get_all_means()

# Plot only for individual sessions
sessions = neuron_group.get_session_list()
for i in range(np.max(sessions)):
    to_plot_list = np.where(sessions == i)
    plt.figure()
    mean_act = neuron_group.plot_all_means(sort=True, plotid=to_plot_list, style='heatmap')
    plt.show()

plt.close('all')

neuronA = neuron_group.neurons[0]
ltrials = neuronA.exp.l_trials
neuronB = neuronA.make_subneuron(ltrials)