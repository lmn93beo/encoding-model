import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import pandas as pd
import utils
import mat4py
import neuron_group_utils
import neuron_utils

# Load data
raw_behavior_summary = mat4py.loadmat('TB41_behavior_summary_090519.mat')
raw_encoding_struct = mat4py.loadmat('TB41_encoding_structs.mat')

# Extract structures from the raw data
one, two, three, four, correct, incorrect, left, right, prev_corr, prev_incorr, \
    stim_onset_per_trial, reward_onset_per_trial = utils.get_multiple_struct_fields_mat4py(raw_behavior_summary,
                        ['one', 'two', 'three', 'four', 'correct', 'incorrect',
                         'left', 'right', 'prev_right', 'prev_wrong', 'stim_onset_per_trial',
                         'reward_onset_per_trial'])

# Being careful about one-indexing
for i in [one, two, three, four, correct, incorrect, prev_corr, prev_incorr, left, right]:
    i -= 1

neural_matmat, predall_matmat = \
    utils.get_multiple_struct_fields_mat4py(raw_encoding_struct, ['neural_matmat', 'predall_matmat'])

# Make an experiment object to capture the experiment
exp_dict = dict(l_trials=left, r_trials=right, one_diff=one, two_diff=two, three_diff=three,
                four_diff=four, correct=correct, incorrect=incorrect, prev_corr=prev_corr,
                prev_incorr=prev_incorr)
exp = neuron_utils.Experiment(exp_dict)

# Make a neuron group
neurons = neuron_group_utils.make_neuron_list(raw_encoding_struct, 'neural_act_mat', np.arange(180), exp)

for neuron in neurons:
    neuron.align_activity(stim_onset_per_trial, [-5, 18])
neuron_group = neuron_group_utils.NeuronGroup(neurons)
neuron_group.plot_all_means()



# Make subgroup based on class
class0 = np.where(neuron_group.classes == 0)[0]
class1 = np.where(neuron_group.classes == 1)[0]
class2 = np.where(neuron_group.classes == 2)[0]
class0_group = neuron_group.make_subgroup_by_neurons(class0)
class1_group = neuron_group.make_subgroup_by_neurons(class1)
class2_group = neuron_group.make_subgroup_by_neurons(class2)

# Plot activity for each group
plt.figure()
class0_group.plot_all_means(color='b')
class1_group.plot_all_means(color='r')
class2_group.plot_all_means(color='g')

# Compare reward and non-reward
class0_group_corr = class0_group.make_subgroup_by_trials(correct)
class0_group_incorr = class0_group.make_subgroup_by_trials(incorrect)
plt.figure(figsize=(15, 15))
class0_group_corr.plot_example_neurons()
class0_group_incorr.plot_example_neurons()
plt.title('Class 0')

plt.figure(figsize=(15, 15))
class1_group_corr = class1_group.make_subgroup_by_trials(correct)
class1_group_incorr = class1_group.make_subgroup_by_trials(incorrect)
class1_group_corr.plot_example_neurons()
class1_group_incorr.plot_example_neurons()
plt.title('Class 1')

plt.figure(figsize=(15, 15))
class2_group_corr = class2_group.make_subgroup_by_trials(correct)
class2_group_incorr = class2_group.make_subgroup_by_trials(incorrect)
class2_group_corr.plot_example_neurons()
class2_group_incorr.plot_example_neurons()
plt.title('Class 2')

