"""
Analyze one-stim data by Rafiq
"""
import utils
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


plt.figure()
utils.plot_example_neurons(subgroups_c[0], np.arange(100))
utils.plot_example_neurons(subgroups_i[0], np.arange(100))


