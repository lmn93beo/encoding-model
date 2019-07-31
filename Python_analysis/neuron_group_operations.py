import utils
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
        self.classes = []
        self.class0cells = []
        self.class1cells = []
        self.class2cells = []
        self.mean_activities = []

    def align_all(self, tpoints, window):
        """
        Align all neurons in group
        :param tpoints: time points to align to (one-indexed)
        :param window: a tuple (start, end), range of time points for cropping
        :return: nothing
        """
        for neuron in self.neurons:
            neuron.align_activity(tpoints, window)

    def classify_all(self):
        """
        Classify all neurons in group
        Keeping a list of neuron classes
        :return: an np array giving the classes for all neurons
        """
        classes = []
        for neuron in self.neurons:
            neuron.classify()
            classes.append(neuron.neuron_class)
        self.classes = np.array(classes)
        self.class0cells = np.where(self.classes == 0)[0]
        self.class1cells = np.where(self.classes == 1)[0]
        self.class2cells = np.where(self.classes == 2)[0]
        return np.array(classes)

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
        return tmax_lst

    def get_all_means(self, style=None):
        """
        :return: the mean activity of neurons in the group in one plot
        """
        # Find the mean activities and sort them
        mean_activities = []
        for neuron in self.neurons:
            assert neuron.mean_activity != []
            mean_activities.append(neuron.mean_activity)
        self.mean_activities = np.array(mean_activities)
        return mean_activities

    def plot_all_means(self, plotid=None, normalize=False, sort=False, style='lines'):
        """
        Plot the mean activity of all neurons
        :param plotid array of cells to plot
        :param sort: if True, will sort by the peak time of the neurons
        :param style: 'lines' or 'heatmap'
        :return: nothing
        """
        assert self.mean_activities != []
        if plotid is None:
            plotid = np.arange(self.n_neurons)

        if sort:
            tmax_lst = self.collect_tmax()
            sort_id = np.argsort(tmax_lst)
            mean_activities = self.mean_activities[sort_id]
        else:
            mean_activities = self.mean_activities

        # Normalize
        if normalize:
            max_activity = np.abs(mean_activities.max(axis=1))
            mean_activities = mean_activities.T / max_activity
            mean_activities = mean_activities.T
        mean_activities = mean_activities[plotid]
        print(mean_activities.shape)

        if style == 'lines':
            plt.plot(mean_activities.T)
        elif style == 'heatmap':
            if normalize:
                plt.imshow(mean_activities, cmap='bwr', aspect='auto', vmin=-1, vmax=1)
            else:
                plt.imshow(mean_activities, cmap='bwr', aspect='auto', vmin=-100, vmax=100)
        else:
            raise ValueError('Invalid style')

# Package a single neuron into an object
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


