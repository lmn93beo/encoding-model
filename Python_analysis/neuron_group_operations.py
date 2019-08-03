import utils
import collections
import data_classes
import numpy as np
import matplotlib.pyplot as plt

class NeuronGroup(object):
    """
    A class for a group of Neuron objects
    """
    def __init__(self, neuron_lst):
        self.neurons = neuron_lst
        self.n_neurons = len(neuron_lst)
        self.exp = neuron_lst[0].exp

        # Check that no two cells have the same id
        id_lst = [neuron.id for neuron in neuron_lst]
        assert len(np.unique(id_lst)) == len(id_lst)

        # Classify all neurons in group
        classes = []
        for neuron in self.neurons:
            if len(neuron.mean_activity) == 0:
                neuron.classify()
            classes.append(neuron.neuron_class)
        self.classes = np.array(classes)
        self.class0cells = np.where(self.classes == 0)[0]
        self.class1cells = np.where(self.classes == 1)[0]
        self.class2cells = np.where(self.classes == 2)[0]

        # Find mean activities for all neurons in group
        mean_activities = []
        for neuron in self.neurons:
            assert neuron.mean_activity != []
            mean_activities.append(neuron.mean_activity)
        self.mean_activities = np.array(mean_activities)

    def align_all(self, tpoints, window):
        """
        Align all neurons in group
        :param tpoints: time points to align to (one-indexed)
        :param window: a tuple (start, end), range of time points for cropping
        :return: nothing
        """
        for neuron in self.neurons:
            neuron.align_activity(tpoints, window)

    def get_session_list(self):
        """
        Get the list of sessions of each neuron
        :return: the np array of sessions
        """
        sessions = []
        for neuron in self.neurons:
            assert neuron.session >= 0
            sessions.append(neuron.session)
        return np.array(sessions)

    def collect_tmax(self, plot_hist=False):
        """
        Collect all the tmax and plot
        :return: array of tmax
        """
        tmax_lst = []
        for neuron in self.neurons:
            assert neuron.t_max_activity > -1
            tmax_lst.append(neuron.t_max_activity)

        if plot_hist:
            plt.hist(tmax_lst)
        return np.array(tmax_lst)

    def make_subgroup_by_neurons(self, neurons):
        """
        Make a subgroup containing only neurons specified in neurons
        :param neurons: a list of neurons to include
        :return: a NeuronGroup object
        """
        subneurons = [self.neurons[i] for i in neurons]
        return NeuronGroup(subneurons)

    def make_subgroup_by_trials(self, trials):
        """
        Make a subneuron for each neuron, and create a new group with specified trials
        :param trials: a list of trials
        :return: a NeuronGroup object
        """
        subneurons = []
        for neuron in self.neurons:
            subneuron = neuron.make_subneuron(trials)
            subneurons.append(subneuron)
        subgroup = NeuronGroup(subneurons)
        return subgroup

    def plot_all_means(self, plotid=None, normalize=False, sort=False, style='lines', side='both'):
        """
        Plot the mean activity of all neurons
        :param plotid array of cells to plot
        :param sort: if True, will sort by the peak time of the neurons
        :param style: 'lines' or 'heatmap'
        :return: an array ncells x Tpts, mean activity of the plotted cells
        """
        assert self.mean_activities != []

        mean_activities = self.mean_activities

        if plotid is not None:
            mean_activities = mean_activities[plotid]
            print('Number to plot:', len(mean_activities))

        if sort:
            tmax_lst = self.collect_tmax()
            if plotid is not None:
                tmax_lst = tmax_lst[plotid]
            sort_id = np.argsort(tmax_lst)
            mean_activities = mean_activities[sort_id]
        else:
            mean_activities = mean_activities


        # Subtract the baseline of first 5 frames
        baseline = np.mean(mean_activities[:, 0:5], axis=1)
        mean_activities = (mean_activities.T - baseline).T


        # Normalize
        if normalize:
            print('# of cells before pruning:', mean_activities.shape)
            max_activity = mean_activities.max(axis=1)
            mean_activities = mean_activities.T / max_activity
            mean_activities = mean_activities.T

            # Eliminate activities below 0
            #mean_activities = mean_activities[max_activity > 0]
            print('Shape of mean_act:', mean_activities.shape)

        print('Number to plot:', mean_activities.shape)
        # Subplot 1: plot the mean activities
        plt.figure()
        plt.subplot('211')
        if style == 'lines':
            plt.plot(mean_activities.T)
        elif style == 'heatmap':
            if normalize:
                plt.imshow(mean_activities, cmap='bwr', aspect='auto', vmin=-1, vmax=1)
            else:
                plt.imshow(mean_activities, cmap='bwr', aspect='auto', vmin=-2, vmax=2)
        else:
            raise ValueError('Invalid style')

        # Subplot 2: plot histograms
        plt.subplot('212')
        self.collect_tmax(True)

        return mean_activities


"""
Some useful functions for working with groups
"""

# Functions for creating objects from raw data
def make_neuron_obj(rawdata, field, cellid):
    """
    Create a list of neurons based on raw data
    :param rawdata: raw data returned by mat4py.loadmat
    :param field: field containing neural_act_mat
    :return: a neuron object
    """
    neuron_activity = []
    ntrials = len(rawdata[field])
    for i in range(ntrials):
        trial_activity = utils.get_struct_field_mat4py(rawdata, 'neural_act_mat', True, i)
        neuron_activity.append(trial_activity[:, cellid])
        neuron = data_classes.Neuron(cellid, neuron_activity)
    return neuron

def make_neuron_list(rawdata, field, cell_lst):
    """
    Create a list of neurons
    :param rawdata: raw data returned by mat4py.loadmat
    :param field: field containing neural_act_mat
    :param cell_lst: lst of ints, cell ids
    :return: a list of neuron objects
    """
    neurons = []
    for cellid in cell_lst:
        print('Making neuron # ', cellid, '...')
        neuron = make_neuron_obj(rawdata, field, cellid)
        neurons.append(neuron)
    return neurons

def compare_conditions(rawdata, field, tpoints, window, cellid, condA, condB):
    """
    Plot graphs to compare between condition A and B
    :param cellid: cell number
    :param condA: np array to indicate trials belonging to condition A
    :param condB: np array to indicate trials belonging to condition B
    :return: nothing
    """
    neuron = make_neuron_obj(rawdata, field, cellid)
    neuron.align_activity(tpoints, window)
    plt.subplot(121)
    neuron.plot_all_trials(condA)
    plt.ylim((-2, 2))
    plt.subplot(122)
    neuron.plot_all_trials(condB)
    plt.ylim((-2, 2))


def make_neuron_group(cell_arr, exp):
    """
    Make a group of neuron
    :param cell_arr: a list of ncells elements, each is a list of ntrials, each trial is a list of Tpts time
    :param exp: an experiment
    :return: a NeuronGroup object
    """
    cell_lst = []
    for i in range(len(cell_arr)):
        neuron_arr = np.array(cell_arr[i]).T
        neuron = data_classes.OneStimNeuron(i, neuron_arr, exp)
        neuron.classify()
        cell_lst.append(neuron)

    curr_session = -1
    curr_ntrials = -1
    for neuron in cell_lst:
        if neuron.ntrials != curr_ntrials:
            curr_session += 1
            print('Changed, current ntrials = ', neuron.ntrials)
            print('Current session = ', curr_session)
            curr_ntrials = neuron.ntrials
        neuron.session = curr_session

    return NeuronGroup(cell_lst)

# Functions for combining neuron groups
def combine_groups_by_trials(group1, group2):
    """
    For combining neuron groups from the same experiment with different trials
    The two groups need to have the same neurons (based on id)
    :param group1: a NeuronGroup object
    :param group2: a NeuronGroup object
    :return: the combined NeuronGroup object
    """
    # Check that the id's of the two groups are the same
    id_group1 = [neuron.id for neuron in group1.neurons]
    id_group2 = [neuron.id for neuron in group2.neurons]
    assert collections.Counter(id_group1) == collections.Counter(id_group2)

    # Combine
    combined_neurons = []
    for i in range(len(id_group1)):
        # Find the corresponding id in group2
        group2_id = np.where(np.array(id_group2) == i)[0][0]
        combined = data_classes.combine_neurons(group1.neurons[i], group2.neurons[group2_id])
        combined_neurons.append(combined)

    return NeuronGroup(combined_neurons)




