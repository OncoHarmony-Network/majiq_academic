from __future__ import division
import fnmatch
import logging
import os
from collections import defaultdict
import json

from voila.lsv import Lsv
from voila.splice_graphics.exonGraphic import ExonGraphic
from voila.splice_graphics.junctionGraphic import JunctionGraphic
from voila.splice_graphics.geneGraphic import GeneGraphic

import shutil
import errno


try:
    import cPickle as pkl
except ImportError:
    try:
        import pickle as pkl
    except ImportError:
        print "[Error] :: Neither pickle nor cPickle are installed. Please, check python dependencies."
        import sys
        sys.exit(1)

try:
    import numpy as np
except ImportError:
    print "[Error] :: Numpy not installed. Please, check python dependencies."
    sys.exit(1)


class PickleEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        if isinstance(obj, np.ndarray):
            return list(obj)
        if isinstance(obj, tuple):
            return list(obj)
        if isinstance(obj, np.int64):
            return int(obj)
        if isinstance(obj, Lsv):
            return obj.to_JSON(PickleEncoder)

        return json.JSONEncoder.default(self, obj)


class LsvGraphicEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        if isinstance(obj, np.ndarray):
            return list(obj)
        if isinstance(obj, tuple):
            return list(obj)
        if isinstance(obj, np.int64):
            return int(obj)
        if isinstance(obj, ExonGraphic):
            return obj.to_JSON(PickleEncoder)
        if isinstance(obj, JunctionGraphic):
            return obj.to_JSON(PickleEncoder)
        if isinstance(obj, GeneGraphic):
            return obj.to_JSON(PickleEncoder)

        return json.JSONEncoder.default(self, obj)


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
    @param position:
    @param direction:
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


def get_mean_step(bins): # TODO: rename
    bins = np.array(bins)
    step = 1 / bins.size
    projection_prod = bins * np.arange(step / 2, 1, step)
    return step, np.sum(projection_prod)


def get_variance(bins, mean):
    """Compute the variance = E[X^2] - (E[X])^2"""
    bins = np.array(bins)
    step_bins = 1 / bins.size
    projection_prod = bins * np.arange(step_bins / 2, 1, step_bins)**2
    return np.sum(projection_prod) - mean**2


def get_prob_delta_psi_greater_v(bins, expected, V=.2):
    """Calculate probability of delta psi outside the acceptable area"""
    bins = np.array(bins)
    step = 2.0 / bins.size
    left = 0
    right = bins.size*2 - 1
    for i, w in enumerate(np.arange(-1+step/2, 1, step)):
        if not left and w > (expected - abs(expected*V)):
            left = i-1
        if right == bins.size*2 and w > (expected + abs(expected*V)):
            right = i
    # print np.sum(bins[:left]), np.sum(bins[right:])
    return (np.sum(bins[:left]) + np.sum(bins[right:]))


def create_array_bins(bins, confidence):
    """
    Recaps bins info from data previously generated and stored in a Pickle file
    @param event_id: to access the bins associated with the event
    @param num_bins: ONLY in DEBUG (it should be retrieved from file)
    @return: a tuple with:
     *.- the mean,
     *.- the coordinates of the confidence interval [coord1, coord2].

    """
    bins = np.array(bins)
    step, mean = get_mean_step(bins)
    conf_interval = find_confidence_interval(bins, mean / step, confidence)
    quartiles_set = find_quartiles(bins)
    variance = get_variance(bins, mean)
    return mean, conf_interval, quartiles_set, variance


def generate_lsv(i, lsvs_bins, confidence, **post_metadata):
    """Parse the info generated by Majiq into a LSV representable object"""

    prefix_matrix_plot = "../templates/static/matrix_plots/"

    # type_set = ('Exon skipping', '5-prime', '3-prime')
    random_num = np.random.random()  # Random number between 0 and 1
    means_psi_list = []
    conf_interval_list = []
    quartile_list = []
    variance_list = []
    coords = None
    for lsv_bins in lsvs_bins:
        m, c, q, v = create_array_bins(lsv_bins, confidence)
        means_psi_list.append(m)
        conf_interval_list.append(c)
        quartile_list.append(q)
        variance_list.append(v)

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
                    matrix_link = prefix_matrix_plot + "mas" + name_plot + ".png"
    else:
        id = "chr" + str(random_num)[-1] + ":" + "10000" + str(random_num)[-2:] + "-" + "1000" + str(random_num)[-3:]
        coords = "chr" + str(random_num)[-1] + ":" + "10000" + str(random_num)[-2:] + "-" + "1000" + str(random_num)
        matrix_link = "#"

    lsv_type = ''
    if 'lsv_type' in post_metadata:
        lsv_type = post_metadata['lsv_type']

    return {
        'name': 'PS'+str(i),
        'id': id,
        'type': lsv_type,
        'bins_list': lsvs_bins,  # bins array
        'mean_psi': means_psi_list,
        'conf_interval': conf_interval_list,
        'quartiles': quartile_list,
        'coords': coords,
        'matrix_link': matrix_link,
        'variance_list': variance_list
    }


def generate_event(i, events_bins, confidence, **post_metadata):
    """Collect all event information from data files"""
    prefix_matrix_plot = "../templates/static/matrix_plots/"

    # type_set = ('Exon skipping', '5-prime', '3-prime')
    random_num = np.random.random()  # Random number between 0 and 1
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
                    matrix_link = prefix_matrix_plot + "mas" + name_plot + ".png"
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
        'matrix_link': matrix_link,
    }


def get_lsv_single_exp_data(majiq_bins_file, confidence, gene_name_list=None, lsv_types=None, logger=None):
    """
    Create a dictionary to summarize the information from majiq output file.
    """
    majiq_data = None
    try:
        majiq_data = pkl.load(open(majiq_bins_file, 'rb'))
    except pkl.PickleError:
        logger.error("Pickle could not load the file. Please, check that the file %s is in Pickle format." % majiq_bins_file, exc_info=1)

    except IOError:
        logger.error("%s doesn't exists." % + majiq_bins_file, exc_info=1)
        sys.exit(1)

    # Load metadata
    metadata_pre = majiq_data[1]
    # metadata = []

    lsv_counter = 0
    lsv_list = []

    genes_dict = defaultdict(list)

    nofilter_genes = not gene_name_list and not lsv_types
    if gene_name_list is None:
        gene_name_list = []

    for i, lsv_meta in enumerate(metadata_pre):
        if nofilter_genes or str(lsv_meta[1]).split(':')[0] in gene_name_list or lsv_meta[2] in lsv_types:
            # print lsv_meta[0], lsv_meta[1], lsv_meta[2]

            bins_array_list = majiq_data[0][i]

            # In 1-way LSVs, create the additional bins set for commodity
            if len(bins_array_list) == 1:
                bins_array_list.append(bins_array_list[-1][::-1])

            lsv_counter += 1
            try:
                # metadata.append([lsv_meta[0], lsv_meta[1], lsv_meta[2]]) #collapse_lsv(lsv_meta[2])])
                lsv_list.append(Lsv(generate_lsv(lsv_counter, bins_array_list, confidence, lsv_type=lsv_meta[2])))
                lsv_list[-1].set_coords(lsv_meta[0])
                lsv_list[-1].set_id(lsv_meta[1])
                lsv_list[-1].set_type(lsv_meta[2])

            except ValueError, e:
                logger.warning("%s produced an error:\n%s. Skipped." % (bins_array_list, e))
                continue

            genes_dict[str(lsv_meta[1]).split(':')[0]].append([lsv_list[-1], lsv_meta])

    return {'event_list':   lsv_list,
            'metadata':     metadata_pre,
            'genes_dict':   genes_dict,
            'meta_exps':    majiq_data[2]}


def extract_bins_info(lsv, threshold, include_lsv):
    expected_psis_bins = []
    excl_inc_perc_list = []
    collapsed_matrices = []

    for junc_matrix in lsv:
        collapsed_matrices.append(collapse_matrix(np.array(junc_matrix)))

    if len(collapsed_matrices)<2:
        collapsed_matrices.append(collapsed_matrices[-1][::-1])

    for bins in collapsed_matrices:
        expected_psis_bins.append(list(bins))
        excl_inc_tuple = find_excl_incl_percentages(bins, threshold)
        excl_inc_perc_list.append(excl_inc_tuple)

        # If the delta is significant (over the threshold) or 'show-all' option, include LSV
        include_lsv = include_lsv or np.any(np.array(excl_inc_tuple)[np.array(excl_inc_tuple)>threshold])
    return expected_psis_bins, excl_inc_perc_list, include_lsv


def get_lsv_delta_exp_data(majiq_out_file, confidence=.95, threshold=.2, show_all=False, gene_name_list=None, logger=None):
    """
    Load lsv delta psi pickle file. It contains a list with 2 elements:
        [0] List with LSV bins matrices
        [1] List with info per LSV

    :param majiq_out_file:
    :param metadata_post:
    :param confidence:
    :param threshold:
    :param logger:
    @return: dictionary
    """
    # Collapse matrix in diagonal
    lsv_matrix_list_info_list = None
    try:
        lsv_matrix_list_info_list = np.array(pkl.load(open(majiq_out_file, 'rb')))
    except pkl.PickleError, e:
        logger.error("Loading the file %s: %s." % (majiq_out_file, e.message), exc_info=1)

    genes_dict = defaultdict(list)

    lsv_list = lsv_matrix_list_info_list[0]
    lsv_info = lsv_matrix_list_info_list[1]

    for i, lsv in enumerate(lsv_list):
        include_lsv = show_all
        gene_name = str(lsv_info[i][1]).split(':')[0]
        if not gene_name_list or gene_name in gene_name_list:
            collapsed_bins, excl_inc_perc_list, include_lsv = extract_bins_info(lsv, threshold, include_lsv)
            if not include_lsv: continue
            try:
                lsv_o = Lsv(generate_lsv(i, collapsed_bins, confidence, lsv_type=lsv_info[i][2]))
                lsv_o.set_excl_incl(excl_inc_perc_list)
                # lsv_list.append(lsv)
                genes_dict[gene_name].append([lsv_o, lsv_info[i]])

            except ValueError, e:
                print e.message
                logger.warning("%s produced an error:\n%s. Skipped." % (lsv_info[i], e))

    logger.info("Number of genes added: %d" % len(genes_dict.keys()))
    # TODO: Extract experiments info from Majiq output file
    experiments_info = [{'name': 'exp1', 'link': '#', 'color': '#e41a1c'},
                        {'name': 'exp2', 'link': '#', 'color': '#377e80'}]

    return {'genes_dict': genes_dict, 'experiments_info': experiments_info}


def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc:  # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        if exc.errno == errno.EEXIST:  # Static folder exists
            shutil.rmtree(dst)
            copyanything(src, dst)
        else:
            raise


def collapse_lsv(lsv_type):
    tab = lsv_type.split('|')
    if len(tab) < 3:
        return lsv_type + '|1'
    min_sss = 20
    min_sst = [20]*20
    res = tab[0]

    ss_list = set()
    ex_list = set()
    for tt in tab[1:]:
        pp= tt.split('e')
        ss2 = int(pp[0])
        min_sss = min(min_sss, ss2)
        try:
            ss3 = int(pp[1].split('.')[1])
        except IndexError, e:
            ss3 = 1
        ext = int(pp[1].split('.')[0])

        min_sst[ext] = min(min_sst[ext], ss3)
    for tt in tab[1:]:
        tab2 = tt.split('e')
        new_ss = int(tab2[0])-min_sss +1
        tab3 = tab2[1].split('.')
        if len(tab3) == 1:
            tab3.append(1)
        new_sst = int(tab3[1])-min_sst[int(tab3[0])] + 1
        ss_list.add(new_ss)
        ex_list.add( '%s.%s'%(tab3[0],new_sst))
        #ss_list += '%s,'%new_ss
        #ex_list += '%s.%s,'%(tab3[0],tab3[1])
    #res += '|%se%s.%s'%(new_ss,tab3[0],tab3[1])
    ss = ','.join([str(x) for x in sorted(ss_list)])
    exs = ','.join([str(x) for x in sorted(ex_list)])
    res += '|%se%s|%s'%(ss,exs, len(tab[1:]))

    return res


def list_files_or_dir(file_or_dir_list, suffix='*', containing='*'):

    if type(file_or_dir_list) != list: return [file_or_dir_list]
    files = []
    for file_or_dir in file_or_dir_list:
        if os.path.isdir(file_or_dir):
            for root, dirnames, filenames in os.walk(file_or_dir):
                for filename in fnmatch.filter(filenames, '*%s*%s' % (containing, suffix)):
                    files.append(os.path.join(root, filename))
            # for file in os.listdir(file_or_dir):
            #     if not suffix or file.endswith(suffix):
            #         files.append(file_or_dir+'/'+file)
        else:
            files.append(file_or_dir)
    return files


def get_logger(logger_name, silent=False, debug=False):
    """
    Returns a logger instance. verbose = False will silence the logger, debug will give
    more information intended for debugging purposes.
    """
    logging_format = "%(asctime)s (PID:%(process)s) - %(levelname)s - %(message)s"
    logging.basicConfig(filename=logger_name, format=logging_format)
    logger = logging.getLogger(logger_name)
    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    ch = logging.StreamHandler()
    if debug:
        ch.setLevel(logging.DEBUG)
    elif not silent:
        ch.setLevel(logging.INFO)
    else:
        ch.setLevel(logging.WARNING)

    formatter = logging.Formatter("%(asctime)s (PID:%(process)s) - %(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger


def create_if_not_exists(my_dir, logger=False):
    """Create a directory path if it does not exist"""
    try:
        if logger:
            logger.info("\nCreating directory %s..." % my_dir)
        os.makedirs(my_dir)
    except OSError:
        if logger:
            logger.info("\nDirectory %s already exists..." % my_dir)


def collapse_matrix(matrix):
    "Collapse the diagonals probabilities in 1-D and return them"
    collapse = []
    #FOR TEST matrix = array([[0, 1, 2, 3, 4, 500], [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5], [100, 1, 2, 3, 4, 5], ])

    matrix_corner = matrix.shape[0]+1
    for i in xrange(-matrix_corner, matrix_corner):
        collapse.append(np.diagonal(matrix, offset=i).sum())

    return np.array(collapse)
