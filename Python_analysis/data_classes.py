import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class Neuron(object):
    """
    A class for a neuron in ACC
    """
    def __init__(self, cellid, activity):
        """
        :param cellid: neuron number
        :param activity: A list of ntrials np arrays, each array is 1-d activity of that trial
        """
        self.id = cellid
        self.activity = activity
        self.ntrials = len(activity)
        self.aligned_activity = []
        self.neuron_class = -1
        self.mean_activity = []
        self.t_max_activity = -1

    def plot_all_trials(self, trials=None, style='lines'):
        """
        For plotting of all trials
        :param trials: list of trials, if none, will plot all the trials
        :return: nothing
        """
        assert self.aligned_activity != []

        if trials is None:
            trials = np.arange(self.ntrials)

        extracted_trials = self.aligned_activity[trials, :]
        print(extracted_trials.shape)
        if style == 'lines':
            plt.plot(extracted_trials.T, 'b', alpha=0.1)
            plt.plot(np.mean(extracted_trials, axis=0), 'r')
            plt.show()
        elif style == 'heatmap':
            plt.imshow(extracted_trials, cmap='bwr')


    def align_activity(self, tpoints, window):
        """
        Align neural activity to time points specified in tpoints
        :param tpoints: time points to align to (one-indexed)
        :param window: a tuple (start, end), range of time points for cropping
        :return: an np array of size ntrials x T corresponding to aligned activity
        """
        arr = []
        for i in range(self.ntrials):
            startT = tpoints[i]
            single_trial = self.activity[i][startT + window[0] - 1 : startT + window[1]]
            arr.append(single_trial)
        self.aligned_activity = np.array(arr)
        return np.array(arr)

    def classify(self):
        """
        Classify the neuron based on the peak of trial-averaged activity
        Peak time: class 0 (0 - 9 frames), class 1 (10-15 frames), class 2 (>= 16 frames)
        :return: class of the neuron
        """
        # Find peak of mean activity
        assert(self.aligned_activity != [])
        mean_activity = np.mean(self.aligned_activity, axis=0)
        self.mean_activity = mean_activity
        t_max_activity = np.argmax(mean_activity)
        self.t_max_activity = t_max_activity
        if t_max_activity < 10:
            self.neuron_class = 0
        elif t_max_activity < 16:
            self.neuron_class = 1
        else:
            self.neuron_class = 2
        return self.neuron_class

    def get_mean_activity(self):
        """
        Get the mean activity of the neuron across all trials
        :return: an np array with the mean activity
        """
        assert self.aligned_activity != []
        if self.mean_activity == []:
            self.mean_activity = np.mean(self.aligned_activity, axis=0)
        return self.mean_activity

class OneStimNeuron(Neuron):
    """
    A class for neurons in 1-stim task
    """
    def __init__(self, cellid, activity):
        """
        :param cellid: cell number
        :param activity: a list of ntrials np array, each array containing cell activity
        """
        self.id = cellid
        self.activity = activity.T
        self.ntrials = activity.shape[1]
        self.aligned_activity = activity.T
        self.neuron_class = -1
        self.session = -1


class Experiment(object):
    """
    A class of trackball experiment
    """
    def __init__(self, rawdata):
        self.rawdata = rawdata
