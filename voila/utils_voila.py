from __future__ import division
import matplotlib
matplotlib.use('Agg')

from collections import defaultdict
import json

import sys
from analysis.matrix import collapse_matrix

from lsv import Lsv
from splice_graphics.exonGraphic import ExonGraphic
from splice_graphics.junctionGraphic import JunctionGraphic
from splice_graphics.geneGraphic import GeneGraphic

import shutil
import errno

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


def get_mean_step(bins):
    bins = numpy.array(bins)
    step = 1 / bins.size
    projection_prod = bins * np.arange(step / 2, 1, step)
    return step, np.sum(projection_prod)


def get_variance(bins, mean):
    """Compute the variance = E[X^2] - (E[X])^2"""
    return 0  # TODO: for now, we are going to skip the variance display
    bins = numpy.array(bins)
    step_bins = 1 / bins.size
    projection_prod = bins * np.arange(step_bins / 2, 1, step_bins)**2
    return np.sum(projection_prod) - mean**2


def create_array_bins(bins, confidence):
    """
    Recaps bins info from data previously generated and stored in a Pickle file
    @param event_id: to access the bins associated with the event
    @param num_bins: ONLY in DEBUG (it should be retrieved from file)
    @return: a tuple with:
     *.- the mean,
     *.- the coordinates of the confidence interval [coord1, coord2].

    """
    bins = numpy.array(bins)
    step, mean = get_mean_step(bins)
    conf_interval = find_confidence_interval(bins, mean / step, confidence)
    quartiles_set = find_quartiles(bins)
    variance = get_variance(bins, mean)
    return mean, conf_interval, quartiles_set, variance


def generate_lsv(i, lsvs_bins, confidence, **post_metadata):
    """Parse the info generated by Majiq into a LSV representable object"""

    PREFIX = "../templates/static/matrix_plots/"

    # type_set = ('Exon skipping', '5-prime', '3-prime')
    random_num = np.random.random()  # Random number between 0 and 1
    means_psi_list = []
    conf_interval_list = []
    quartile_list = []
    variance_list = []
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
                    matrix_link = PREFIX + "mas" + name_plot + ".png"
    else:
        id = "chr" + str(random_num)[-1] + ":" + "10000" + str(random_num)[-2:] + "-" + "1000" + str(random_num)[-3:]
        coords = "chr" + str(random_num)[-1] + ":" + "10000" + str(random_num)[-2:] + "-" + "1000" + str(random_num)
        matrix_link = "#"

    return {
        'name': 'PS'+str(i),
        'id': id,
        # TODO: information missing
        'type': 'Exon skipping',  # TODO: information missing
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
    PREFIX = "../templates/static/matrix_plots/"

    # type_set = ('Exon skipping', '5-prime', '3-prime')
    random_num = numpy.random.random()  # Random number between 0 and 1
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
        'matrix_link': matrix_link,
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

def get_delta_exp_data(majiq_out_file, metadata_post=None, confidence=.95, threshold=.2):
    """
    Create a dictionary to summarize the delta information.
    """
    # Collapse matrix in diagonal
    try:
        matrix_paired = np.array(pkl.load(open(majiq_out_file, 'rb')))
    except pkl.PickleError, e:
        print "[Error] :: Loading the file %s: %s." % (majiq_out_file, e.message)
        sys.exit(1)

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

    # events_list = sample_event_list(events_list)
    # events_indexes = []
    # for event in events_list:
    #     events_indexes.append(int(event.number)-1)

    # TODO: Extract experiments info from Majiq output file
    experiments_info = [{'name': 'experiment1', 'link': '#', 'color': '#e41a1c'},
                        {'name': 'experiment2', 'link': '#', 'color': '#377e80'}]

    return {'event_list': events_list, 'experiments_info': experiments_info}


def get_lsv_single_exp_data(majiq_bins_file, confidence, gene_name_list=None, lsv_types=None, bed_file=None):
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

    # Load metadata
    metadata_pre = bins_matrix[1]
    metadata = []

    lsv_counter = 0
    lsv_list = []

    genes_dict = defaultdict(list)

    bed_dict = defaultdict(str)
    if bed_file:
        try:
            #'/Users/abarrera/workspace/majiq/data/builder_output/ensambl.mm10.sorted.bed.tailored'
            with open(bed_file) as bedfile:
                for line in bedfile:
                    fields = line.rstrip().split()
                    bed_dict[fields[3]] = fields[0]
        except Exception, e:
            print e.message
            pass

    nofilter_genes = not gene_name_list and not lsv_types
    if gene_name_list is None:
        gene_name_list = []

    for i, lsv_meta in enumerate(metadata_pre):
        if nofilter_genes or str(lsv_meta[1]).split(':')[0] in gene_name_list or lsv_meta[2] in lsv_types:
            print lsv_meta[0], lsv_meta[1], lsv_meta[2]

            metadata.append([lsv_meta[0], lsv_meta[1], lsv_meta[2]]) #collapse_lsv(lsv_meta[2])])
            bins_array_list = bins_matrix[0][i]

            # In 1-way LSVs, create the additional bins set for commodity
            if len(bins_array_list) == 1:
                bins_array_list.append(bins_array_list[-1][::-1])

            lsv_counter += 1
            lsv_list.append(Lsv(generate_lsv(lsv_counter, bins_array_list, confidence)))
            try:
                lsv_list.append(Lsv(generate_lsv(lsv_counter, bins_array_list, confidence)))
            except ValueError, e:
                print "[WARNING] :: %s produced an error:\n%s (Skipped)" % (bins_array_list, e)
                continue

            genes_dict[str(lsv_meta[1]).split(':')[0]].append([lsv_list[-1], [lsv_meta[0], lsv_meta[1], lsv_meta[2], bed_dict[lsv_meta[1].split(':')[0]]]])

    return {'event_list':   lsv_list,
            'metadata':     metadata,
            'genes_dict':    genes_dict }


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


def get_lsv_delta_exp_data(majiq_out_file, confidence=.95, threshold=.2, show_all=False, gene_name_list=None):
    """
    Load lsv delta psi pickle file. It contains a list with 2 elements:
        [0] List with LSV bins matrices
        [1] List with info per LSV

    :param majiq_out_file:
    :param metadata_post:
    :param confidence:
    :param threshold:
    @return: dictionary
    """
    # Collapse matrix in diagonal
    try:
        lsv_matrix_list_info_list = np.array(pkl.load(open(majiq_out_file, 'rb')))
    except pkl.PickleError, e:
        print "[Error] :: Loading the file %s: %s." % (majiq_out_file, e.message)
        sys.exit(1)

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
                lsv_o = Lsv(generate_lsv(i, collapsed_bins, confidence))
                lsv_o.set_excl_incl(excl_inc_perc_list)
                # lsv_list.append(lsv)
                genes_dict[gene_name].append([lsv_o, lsv_info[i]])

            except ValueError, e:
                print "[WARNING] :: %s produced an error:\n%s (Skipped)" % (repr(lsv_o), e)


    print "Number of genes added: %d" % len(genes_dict.keys())
    # TODO: Extract experiments info from Majiq output file
    experiments_info = [{'name': 'experiment1', 'link': '#', 'color': '#e41a1c'},
                        {'name': 'experiment2', 'link': '#', 'color': '#377e80'}]

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