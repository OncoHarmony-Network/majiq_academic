from __future__ import division
import json
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy
from event import Event


try:
    import cPickle as pkl
except ImportError:
    try:
        import pickle as pkl
    except ImportError:
        print "[Error] Neither pickle nor cPickle are installed. Please, check python dependencies."
        sys.exit(1)

try:
    import numpy as np
except ImportError:
    print "[Error] Numpy not installed. Please, check python dependencies."
    sys.exit(1)


def find_excl_incl_percentages(bins, threshold):
    """
    Calculate the percentage of inclusion/exclusion given the differential bins set

    @param bins: array of bins where the sum of all elements is equal to 1
    @param threshold: the absolute value of the minimum differential delta PSI (e.g. 0.2)
    @return array of exclusion and inclusion percentages.
    """
    edges = np.linspace(-1, 1, num=len(bins))
    edges_bins = edges * bins
    bins_per_threshold = int(len(bins) * (threshold / 2))
    return [-sum(edges_bins[:int(len(bins) / 2) - bins_per_threshold]),
            sum(edges_bins[int(len(bins) / 2) + bins_per_threshold:])]


# TODO: This method should be part of Juan's library within Majiq
def padding(accumulated, converted_confidence, bins, position, direction):
    """
    When there is only one way to look for the coverage.

    @param accumulated:
    @param converted_confidence:
    @param bins:
    @param right:
    @param left:
    @return: position of the bin where the confidence is reached
    """
    bins_len = bins.size
    while accumulated < converted_confidence:
        accumulated += bins.item(position)
        position += direction
        if position < 0:
            return 0
        if position >= bins_len:
            return bins_len - 1
    return position


# TODO: This method should be part of Juan's library within Majiq
def find_confidence_interval(bins, meanX, confidence=0.95):
    """
    This method returns a confidence interval around the mean
    @param bins:
    @param mean:
    """

    pos_mean = max(int(meanX) - 1, 0)
    accumulated = bins.item(pos_mean)
    if accumulated >= confidence:  # TODO: what if the mean has more than the conf. interval required?
        return [pos_mean, pos_mean + 1]

    left = pos_mean - 1
    right = pos_mean + 1

    while True:
        if left < 0:  # No more variance to catch on the left
            return [0, padding(accumulated, confidence, bins, right, 1)]

        if right >= bins.size:  # No more variance to catch on the right
            return [padding(accumulated, confidence, bins, left, -1), bins.size - 1]

        # Choose the side with bigger variance
        if bins.item(left) > bins.item(right):
            accumulated += bins.item(left)
            offset = (-1, 0)
        else:
            accumulated += bins.item(right)
            offset = (0, 1)

        if accumulated >= confidence:
            return [left, right]

        left, right = left + offset[0], right + offset[1]


def find_quartiles(bins):
    """
    Iterate the bins finding how the PSI values are distributed.
    @param bins:
    @return:
    """
    accumulated = 0
    bin_index = 0
    quartiles_set = (.10, .25, .50, .75, .90)
    quartiles_values = []
    quartiles_index = 0
    while quartiles_index < len(quartiles_set) and accumulated < 1:
        accumulated += bins[bin_index]
        while quartiles_index < len(quartiles_set) and accumulated >= quartiles_set[quartiles_index]:
            quartiles_values.append(bin_index)
            quartiles_index += 1
        bin_index += 1

    return quartiles_values


def get_mean_step(bins):
    step = 1 / bins.size
    projection_prod = bins * np.arange(step / 2, 1, step)
    return step, np.sum(projection_prod)


def create_array_bins(bins, confidence):
    """
    Recaps bins info from data previously generated and stored in a Pickle file
    @param event_id: to access the bins associated with the event
    @param num_bins: ONLY in DEBUG (it should be retrieved from file)
    @return: a tuple with:
     *.- the mean,
     *.- the coordinates of the confidence interval [coord1, coord2].

    """

    # find where the confidence interval falls
    #   *IMPORTANT*! there are 39 bins

    step, mean = get_mean_step(bins)
    conf_interval = find_confidence_interval(bins, mean / step, confidence)
    quartiles_set = find_quartiles(bins)

    return mean, conf_interval, quartiles_set


def generate_event(i, events_bins, confidence, **post_metadata):
    """Collect all event information from data files"""
    PREFIX = "../templates/static/matrix_plots/"

    # type_set = ('Exon skipping', '5-prime', '3-prime')
    random_num = random()  # Random number between 0 and 1
    bins_info = create_array_bins(events_bins, confidence)
    events_bins.tolist()

    if 'names' in post_metadata and post_metadata['names']:
        id = post_metadata['names'][i]
        try:
            coords = str(id).split(":")[1]
        except IndexError:
            # Name in unknown format, skip coordinates parsing ...
            pass

        matrix_link = "#"
        if 'keys_plots' in post_metadata:
            # Correct key names for plots
            for key_plot, name_plot in post_metadata['keys_plots'].items():
                if str(id).startswith(str(key_plot).split("_")[0]):
                    matrix_link = PREFIX + "mas" + name_plot + ".png"
    else:
        id = "chr" + str(random_num)[-1] + ":" + "10000" + str(random_num)[-2:] + "-" + "1000" + str(random_num)[-3:]
        coords = "chr" + str(random_num)[-1] + ":" + "10000" + str(random_num)[-2:] + "-" + "1000" + str(random_num)
        matrix_link = "#"

    return {
        'number': i,
        'id': id,
        # TODO: information missing
        'type': 'Exon skipping',  # TODO: information missing
        'bins': events_bins.tolist(),  # bins array
        'mean_psi': bins_info[0],
        'conf_interval': bins_info[1],
        'quartiles': bins_info[2],
        'coords': coords,
        'matrix_link': matrix_link
    }


def get_single_exp_data(majiq_bins_file=None, metadata_pre=None, metadata_post=None, confidence=.95):
    """
    Create a dictionary to summarize the information from majiq output file.
    """
    if metadata_post['names']:
        try:
            event_names = pkl.load(open(metadata_post['names'], 'rb'))
        except pkl.PickleError:
            print "[Error] :: Pickle could not load the file. Please, check that the file %s is in Pickle format." % metadata_post
            sys.exit(1)
        except IOError:
            print "[Error] :: %s doesn't exists." % metadata_post
            sys.exit(1)
    else:
        event_names = None

    try:
        bins_matrix = pkl.load(open(majiq_bins_file, 'rb'))
    except pkl.PickleError:
        print "[Error] :: Pickle could not load the file. Please, check that the file %s is in Pickle format." % majiq_bins_file
        sys.exit(1)
    except IOError:
        print "[Error] :: %s doesn't exists." % majiq_bins_file
        sys.exit(1)
    event_list = get_event_list_from_bins(bins_matrix, confidence, names=event_names)

    # Load metadata
    meta_pre = None
    if metadata_pre:
        try:
            meta_pre = pkl.load(open(metadata_pre, 'rb'))
        except pkl.PickleError:
            print "[Error] :: Pickle could not load the file. Please, check that the file %s is in Pickle format." % metadata_pre
            sys.exit(1)
        except IOError:
            print "[Error] :: %s doesn't exists." % metadata_pre
            sys.exit(1)

    return {'event_list': event_list, 'metadata_pre': meta_pre}


def get_event_list_from_bins(bins_matrix, confidence=.95, **kwargs):
    """Process bins representing a sampling from Majiq PSI distributions."""
    event_counter = 0
    event_list = []

    for bins_array in bins_matrix:
        event_list.append(Event(generate_event(event_counter, bins_array, confidence, **kwargs)))
        event_counter += 1

    return event_list


def sample_event_list(event_list, size=30):
    sample = []
    # sampled = 0
    # sample_intervals = [2 for i in range(size)]
    # for event in event_list:
    #     if sampled == 2*size:
    #         break
    #     for j in xrange(size-1):
    #         if sample_intervals[j] and 2*j/size-1 <= event.mean_psi < 2*(j+1)/size-1:
    #             sample.append(event)
    #             sampled += 1
    #             sample_intervals[j] -= 1
    #
    # import random
    # sample.extend([event_list[random.randrange(0, len(event_list))] for i in range(10)])
    selected_list = [
        1847,
        103,
        2275,
        1984,
        918,
        2813,
        3040,
        4013,
        4,
        1811,
        41,
        711,
        513,
        1029,
        311,
        1033,
        895,
        3738,
        824,
        872,
        1,
        439,
        3244,
        4087,
        9,
        3893,
        2484,
        4214,
        2038,
        286,
        286,
        4130,
        4173
    ]

    selected_list = [i-1 for i in selected_list]

    for event in event_list:
        if event.number in selected_list:
            sample.append(event)

    # for event in event_list:
    #     excl_incl = event.get_excl_incl()
    #     min_bins = 1
    #     for bin in event.bins:
    #         if bin > 0.2:
    #             min_bins -= 1
    #             if min_bins == 0: break
    #
    #     if min_bins == 0:
    #         sample.append(event)
    #         # if excl_incl[0]>0.01 and excl_incl[1]>0.01:
    #         #     sample.append(event)

    return sample


def get_delta_exp_data(majiq_out_file, metadata_post=None, confidence=.95, threshold=.2):
    """
    Create a dictionary to summarize the delta information.
    """
    # Collapse matrix in diagonal
    try:
        matrix_paired = np.array(pkl.load(open(majiq_out_file, 'rb')))
    except pkl.PickleError, e:
        print "[Error] :: Loading the file %s: %s." % (majiq_out_file, e.message)

    if metadata_post['names']:
        try:
            event_names = pkl.load(open(metadata_post['names'], 'rb'))
        except pkl.PickleError:
            print "[Error] :: Pickle could not load the file. Please, check that the file %s is in Pickle format." % metadata_post
            sys.exit(1)
        except IOError:
            print "[Error] :: %s doesn't exists." % metadata_post['names']
            sys.exit(1)

        keys_plots = None
        if metadata_post['keys_plots']:
            try:
                keys_plots = pkl.load(open(metadata_post['keys_plots'], 'rb'))
            except pkl.PickleError:
                print "[Error] :: Pickle could not load the file. Please, check that the file %s is in Pickle format." % metadata_post
                sys.exit(1)
            except IOError:
                print "[Error] :: %s doesn't exists." % metadata_post['keys_plots']
                sys.exit(1)

    else:
        event_names = None
        keys_plots = None


    bins_list = []
    events_list = []

    for event in matrix_paired:
        bins_list.append(collapse_matrix(np.array(event)))

    # Add information exclusive from Delta Psi files: experiments info, percentages of incl. excl., etc.
    for event in get_event_list_from_bins(bins_list, confidence, names=event_names, keys_plots=keys_plots):
        event.set_excl_incl(find_excl_incl_percentages(event.get_bins(), threshold))
        event.mean_psi = event.mean_psi * 2 - 1
        events_list.append(event)

    events_list = sample_event_list(events_list)
    # events_indexes = []
    # for event in events_list:
    #     events_indexes.append(int(event.number)-1)

    # TODO: Extract experiments info from Majiq output file
    experiments_info = [{'name': 'experiment1', 'link': '#', 'color': '#e41a1c'},
                        {'name': 'experiment2', 'link': '#', 'color': '#377e80'}]

    return {'event_list': events_list, 'experiments_info': experiments_info}


def get_lsv_single_exp_data(majiq_bins_file, meta_preprocess, confidence):
    """
    Create a dictionary to summarize the information from majiq output file.
    """

    try:
        bins_matrix = pkl.load(open(majiq_bins_file, 'rb'))
    except pkl.PickleError:
        print "[Error] :: Pickle could not load the file. Please, check that the file %s is in Pickle format." % majiq_bins_file
        sys.exit(1)
    except IOError:
        print "[Error] :: %s doesn't exists." % + majiq_bins_file
        sys.exit(1)

    event_counter = 0
    event_list = []

    for bins_array in bins_matrix[range(0, 50)]:
        event_counter += 1
        event_list.append(Event(generate_event(event_counter, bins_array, confidence)))

    metadata = None
    # Load metadata
    if meta_preprocess:
        try:
            metadata = pkl.load(open(meta_preprocess, 'rb'))
        except pkl.PickleError:
            print "[Error] :: Pickle could not load the file. Please, check that the file %s is in Pickle format." % meta_preprocess
            sys.exit(1)
        except IOError:
            print "[Error] :: %s doesn't exists." % + meta_preprocess
            sys.exit(1)

    return {'event_list': event_list,
            'metadata': metadata}


def collapse_matrix(matrix):
    """Collapse the diagonals probabilities in 1-D and return them"""
    collapse = []

    matrix_corner = matrix.shape[0]
    for i in xrange(-matrix_corner, matrix_corner):
        collapse.append(diagonal(matrix, offset=i).sum())

    return np.array(collapse)


# So far, this is not called anywhere cos the data should be coming in python format already. This is an ad-hoc solution
# to read Matlab data. Using ipython, load this function, use it with a Matlab matrix and dump it using json.
def load_matlab_mat(matlab_mat):
    import scipy.io
    from scipy.stats.mstats import mquantiles
    from collections import defaultdict

    mat = scipy.io.loadmat(matlab_mat)
    to_voila = defaultdict(lambda: defaultdict())

    for key in mat:
        if not key.startswith('__'):
            for way in range(mat[key].shape[1]):
                # print key, way, np.mean(mat[key][:,way]), mquantiles(mat[key][:,way], prob=(.1,.25, .5, .75, .9))
                to_voila[key]["PSI"+str(way+1)]['mean'] = np.mean(mat[key][:, way])
                to_voila[key]["PSI"+str(way+1)]['quantiles'] = mquantiles(mat[key][:, way], prob=(.1, .25, .5, .75, .9))
                to_voila[key]["PSI"+str(way+1)]['var'] = np.var(mat[key][:, way])


class PickleEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        if isinstance(obj, numpy.ndarray):
            return list(obj)
        if isinstance(obj, tuple):
            return list(obj)
        if isinstance(obj, numpy.int64):
            return int(obj)

        return json.JSONEncoder.default(self, obj)
