"""
Analyze one-stim data by Rafiq
"""

import numpy as np
import matplotlib.pyplot as plt
import mat4py
import data_classes
import neuron_group_operations

filedir = 'C:/Users/Sur lab/Documents/RafiqAnalysis'
raw_i = mat4py.loadmat(filedir + '/raw_i.mat')['raw_i']
raw_c = mat4py.loadmat(filedir + '/raw_c.mat')['raw_c']

# Parse into numpy arrays
cell_lst = []
for i in range(len(raw_i)):
    cell_arr = np.array(raw_i[i])
    neuron = data_classes.OneStimNeuron(i, cell_arr)
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
    #print(neuron.session)

neuron_group = neuron_group_operations.NeuronGroup(cell_lst)
neuron_group.classify_all()
neuron_group.get_all_means()
neuron_group.plot_all_means(sort=True, style='heatmap', normalize=True)
plt.show()

