import numpy as np
import matplotlib.pyplot as plt
import neuron_group_utils
import seaborn as sns

class Neuron(object):
    """
    A class for a neuron in ACC
    """
    def __init__(self, cellid, activity, exp):
        """
        :param cellid: neuron number
        :param activity: A list of ntrials np arrays, each array is 1-d activity of that trial
        :param exp: An experiment object
        """
        self.id = cellid
        self.activity = activity
        self.ntrials = len(activity)
        self.aligned_activity = []
        self.neuron_class = -1
        self.mean_activity = []
        self.stderr_activity = []
        self.t_max_activity = -1
        self.exp = exp

    def plot_all_trials(self, trials=None, style='heatmap'):
        """
        For plotting of all trials
        :param trials: list of trials, if none, will plot all the trials
        :return: nothing
        """
        assert len(self.aligned_activity) > 0

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


    def align_activity(self, tpoints):
        """
        Align neural activity to time points specified in tpoints
        :param tpoints: time points to align to (one-indexed)
        :param window: a tuple (start, end), range of time points for cropping
        :return: an np array of size ntrials x T corresponding to aligned activity
        """
        arr = []
        window = self.exp.window
        for i in range(self.ntrials):
            startT = tpoints[i]

            # Handle case where window falls out of trial
            if startT + window[0] - 1 < 0 or startT + window[1] > len(self.activity[i]):
                print('Warning: trial %d alignment window falls outside of trial\n' % i)
                single_trial = np.full(window[1] - window[0] + 1, np.nan)
            else:
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
        assert len(self.aligned_activity) > 0
        self.get_mean_activity()
        #self.get_mean_activity_sides()
        t_max_activity = np.argmax(self.mean_activity)
        self.t_max_activity = t_max_activity

        # Class 0: before 1.6 s, class 1: 1.6-3s, class 2: after 3s
        BORDER1 = 0.6
        BORDER2 = 2.0
        if t_max_activity + self.exp.window[0] < BORDER1 * self.exp.rate:
            self.neuron_class = 0
        elif t_max_activity + self.exp.window[0] < BORDER2 * self.exp.rate:
            self.neuron_class = 1
        else:
            self.neuron_class = 2
        return self.neuron_class

    def get_mean_activity(self):
        """
        Get the mean activity of the neuron across all trials
        :return: an np array with the mean activity
        """
        assert len(self.aligned_activity) > 0
        if self.mean_activity == []:
            self.mean_activity = np.nanmean(self.aligned_activity, axis=0)
            self.stderr_activity = np.nanstd(self.aligned_activity, axis=0) / np.sqrt(self.ntrials)
        else:
            print('Warning: mean_activity already computed, recomputing...')
        return self.mean_activity

    def get_mean_activity_sides(self):
        """
        Get the mean activity of left and right trials, separately
        :return: nothing
        """
        assert len(self.aligned_activity) > 0

        # Get mean of left trials
        trials_left = self.exp.l_trials
        trials_right = self.exp.r_trials
        left_activity = self.aligned_activity[trials_left]
        right_activity = self.aligned_activity[trials_right]

        self.mean_activity_left = np.mean(left_activity, axis=0)
        self.mean_activity_right = np.mean(right_activity, axis=0)
        print('mean activity left dimensions : ', self.mean_activity_left.shape)


    def plot_side(self, side='l', style='heatmap'):
        """
        Plot all trial on either left or right side
        :param side: 'l' or 'r'
        :return: nothing
        """
        if side == 'l':
            trials_to_plot = self.exp.l_trials
        elif side == 'r':
            trials_to_plot = self.exp.r_trials
        else:
            raise ValueError("Side must be 'l' or 'r'")

        self.plot_all_trials(trials=trials_to_plot, style=style)

    def make_subneuron(self, trials):
        """
        Make a neuron that is identical with self, which only contains the trials specified in trials
        :param trials: a list of trials to subsample
        :return: a Neuron object
        """
        assert len(self.aligned_activity) > 0
        selected_trials = [self.activity[i] for i in trials]
        neuron = Neuron(self.id, selected_trials, self.exp)
        neuron.aligned_activity = self.aligned_activity[trials]
        neuron.classify()

        return neuron

class OneStimNeuron(Neuron):
    """
    A class for neurons in 1-stim task
    """

    def __init__(self, cellid, activity, exp):
        """
        :param cellid: cell number
        :param activity: a list of ntrials np array, each array containing cell activity
        or a ntrials x Tpts np array
        """
        self.id = cellid
        self.activity = activity
        self.exp = exp
        self.ntrials = len(activity)
        self.aligned_activity = activity
        self.stderr_activity = np.std(self.aligned_activity, axis=0) / np.sqrt(self.ntrials)
        self.mean_activity = np.mean(self.aligned_activity, axis=0)
        self.neuron_class = -1
        self.t_max_activity = -1
        self.session = -1

    def copy(self):
        copy_neuron = OneStimNeuron(self.id, self.activity[:], self.exp)
        copy_neuron.neuron_class = self.neuron_class
        copy_neuron.session = self.session
        return copy_neuron


class Experiment(object):
    """
    A class of trackball experiment
    """
    def __init__(self, exp_dict):
        #self.l_trials = exp_dict['l_trials']
        #self.r_trials = exp_dict['r_trials']
        #self.ntrials = len(self.l_trials) + len(self.r_trials)
        #self.exp_dict = exp_dict
        self.rate = exp_dict['rate']
        self.window = exp_dict['window']
        self.animal = exp_dict['animal']

def combine_neurons(neuron1, neuron2):
    """
    For combining two neurons with different trials
    :param neuron1: a Neuron object
    :param neuron2: a Neuron object
    :return: a Neuron object - with combined trials
    """
    assert neuron1.id == neuron2.id
    assert neuron1.activity.shape[1] == neuron2.activity.shape[1]
    assert neuron1.session == neuron2.session

    combined_activity = np.vstack((neuron1.activity, neuron2.activity))
    combined = OneStimNeuron(neuron1.id, combined_activity, neuron1.exp)
    combined.session = neuron1.session
    combined.classify()

    return combined

