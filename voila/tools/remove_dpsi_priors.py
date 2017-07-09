from voila.tools import Tool
from voila.io_voila import Voila
from voila.tools.io_voila_caleb import import_dpsi_pandas
from voila.tools.utils import find_files
import pickle as pkl
import numpy as np
import pdb
import os

from voila.tools.utils.merge_dicts import merge_dicts
from voila.tools.utils.percent_through_list import percent_through_list
from voila.utils.voila_log import voila_log

LOG = voila_log()


class ThisisRemoveDpsiPriors(Tool):
    help = 'Intelligently find the voila_deltapsi file, priormatrix.pkl file, ' \
           'and the tab file given a directory and a comparison name, then remove' \
           'the prior from the expected dPSI values, and write the results to file' \
           'as a pickle dictionary.'

    def arguments(self):
        parser = self.get_parser()
        parser.add_argument('directory',
                            type=str,
                            help='Directory where voila texts are.')
        mutually_excl_grp = parser.add_mutually_exclusive_group(required=True)
        mutually_excl_grp.add_argument('--comparison_name',
                                       type=str,
                                       help='Comparison name you want to remove priors for.')
        mutually_excl_grp.add_argument('--deltapsi_run_lines',
                                       type=str,
                                       help='The file path to the scriptthat was used to run majiq and '
                                            'voila deltapsi. All of the\'deltapsi comparisons in the script '
                                            'will be automatically processed.')
        return parser

    def run(self, args):
        second_arg = args.comparison_name if args.comparison_name else args.delta_psi_run_lines
        wrapper(directory=os.path.abspath(args.directory), second_arg=second_arg)


def wrapper(directory, second_arg):
    """

    :param directory: directory where majiq was run
    :param second_arg: either a single deltapsi comparison name
                        OR a file path to the script used to run majiq/voila deltapsi
    :return: nothing
    """
    if os.path.isfile(second_arg):
        comparisons = get_comparisons_from_runlines(runline_fp=second_arg)
    else:  # else a single comparison provided
        comparisons = [second_arg]
    for comparison in comparisons:
        dpsi_voila, dpsi_prior, dpsi_tsv = get_matched_voila_files(directory=directory,
                                                                   comparison_name=comparison)
        LOG.info("Removing priors using:\n%s\n%s\n%s" % (
            "\n".join(dpsi_tsv),
            dpsi_prior,
            dpsi_voila))
        unique_outpaths = set()
        LOG.info("Note: there are %s voila tab files to be processed." % len(dpsi_tsv))
        pdb.set_trace()
        priors_removed = remove_dpsi_priors(deltapsi_voila=dpsi_voila[0],
                                            deltapsi_prior=dpsi_prior[0],
                                            deltapsi_tabfile=dpsi_tsv)
        for tab_file, prior_rem in zip(dpsi_tsv, priors_removed):
            outname = os.path.basename(tab_file)
            outname = outname.replace("deltapsi_deltapsi", "deltapsi_no_prior")
            outname = outname.replace("tsv", "pickle")
            outdir = os.path.dirname(tab_file)
            outpath = os.path.join(outdir, outname)
            if outpath in unique_outpaths:
                LOG.error("Oops, the same exact outpath is being used more than once.. files"
                          "from this analysis will be overwritten...")
                exit(1)
            # make sure you aren't going to overwrite a result you just created...
            unique_outpaths.add(outpath)
            pkl.dump(prior_rem, open(outpath, "w"))
            LOG.info("Wrote results to file: %s" % outpath)


def get_comparisons_from_runlines(runline_fp):
    """

    :param runline_fp: path to majiq deltapsi and voila deltapsi run lines script
    :return: [deltapsi comparison names]
    """
    comparisons = list()
    with open(runline_fp, "r") as handle:
        for line in handle:
            if line.startswith("majiq deltapsi"):
                line = line.rstrip("\n\r")
                line_split = line[line.find("--name"):].split(" ")
                comp1 = line_split[1]
                comp2 = line_split[2]
                comparison = comp1+"_"+comp2
                comparisons.append(comparison)
    return comparisons


def get_matched_voila_files(directory, comparison_name):
    """
    Given a comparison name and a directory, return file paths for the
            voila_deltapsi file,
            priormatrix.pkl file,
            and the voila tab-separated file(s) <- there could be more than one...

    :rtype: str
    :param directory: where majiq was run
    :param comparison_name: which deltapsi data to get
    :return: (deltapsi_voila_fp, deltapsi_prior_fp, deltapsi_tabfile_fp)
    """
    deltapsi_voila_fp = find_files.find_files(path=directory,
                                              pattern="%s.deltapsi.voila" % comparison_name,
                                              recursive=True)
    deltapsi_prior_fp = find_files.find_files(path=directory,
                                              pattern="%s.priormatrix.pkl" % comparison_name,
                                              recursive=True)
    deltapsi_tabfile_fp = find_files.find_files(path=directory,
                                                pattern="%s.deltapsi_deltapsi.tsv" % comparison_name,
                                                recursive=True)
    if len(deltapsi_voila_fp) != 1:
        raise RuntimeError("Didn't find 1 deltapsi_voila file, instead found %s" % len(deltapsi_voila_fp))
    if len(deltapsi_prior_fp) != 1:
        raise RuntimeError("Didn't find 1 deltapsi_prior file, instead found %s" % len(deltapsi_prior_fp))
    if len(deltapsi_tabfile_fp) == 0:
        raise RuntimeError("Didn't find any deltapsi_tabfile files")
    return deltapsi_voila_fp, deltapsi_prior_fp, deltapsi_tabfile_fp


def remove_dpsi_priors(deltapsi_voila, deltapsi_prior, deltapsi_tabfile):
    """
    :param deltapsi_voila: *.deltapsi.voila
    :param deltapsi_prior: *.priomatrix.pkl
    :param deltapsi_tabfile: list of [*tsv tab output text file of voila]
     Returns a list of results; each elem in same order as elem in [deltapsi_tabfile]
    :return: [ {LSV_IDs:[junction_dPSIs_minus_prior]} ]
    """
    res = list()
    with Voila(deltapsi_voila, 'r') as v:
        LOG.info('loading hdf5 file: %s' % deltapsi_voila)
        tissue_lsvs = v.get_voila_lsvs()
        LOG.info("done")
        LOG.info("Loading priormatrix file %s" % deltapsi_prior)
        tissue_priors = np.array(get_file(deltapsi_prior))
        LOG.info("done")
        tsv_dict = dict()
        for tsv in deltapsi_tabfile:
            this_tsvs_juncs = get_tsv_junctions(tsv)
            tsv_dict[tsv] = this_tsvs_juncs
        # As long as they came from the same build, we know any given LSV will have
        # all the possible junctions, so we just need a master list of all LSVs here
        master_junction_dict = merge_dicts(*tsv_dict.values())
        lsv_ids = list(master_junction_dict.keys())
        n_lsvs = float(len(lsv_ids))  # floating for some reason
        dist_no_priors_edpsi_all = dict()
        i = 0.0
        indeces_at_x_percent = percent_through_list(len(lsv_ids), 0.01)
        for tissue_lsv in tissue_lsvs:
            if i in indeces_at_x_percent:
                perc = indeces_at_x_percent[i]
                LOG.info("Processed %s%% of the LSVs... " % perc)
            dist_no_priors_edpsi = list()
            # This script will work even if you truncate the tsv file
            # It will simply skip trying to process lsvs not in the provided tsv file..
            if tissue_lsv.lsv_id not in lsv_ids:
                continue
            elif i == n_lsvs:
                break
            for junction_number in range(len(master_junction_dict[tissue_lsv.lsv_id])):
                prior_info = np.log(collapse_matrix(tissue_priors[0, :, :]))
                mle_bins = np.exp(np.log(np.array(tissue_lsv.bins[int(junction_number)])) - prior_info)
                mle_bins /= 1. * np.sum(mle_bins)
                dist_no_priors_edpsi.append(expected_dpsi(mle_bins))
            dist_no_priors_edpsi_all[tissue_lsv.lsv_id] = dist_no_priors_edpsi
            i += 1.0
        for tsv in deltapsi_tabfile:
            # get all the prior-removed dPSI data for only the LSV IDs
            # in this given tsv
            this_res = {k: dist_no_priors_edpsi_all.get(k) for k in tsv_dict[tsv].keys()}
            res.append(this_res)
    return res


def get_tsv_junctions(tsv_file):
    """
    :param tsv_file: path to tsv
    :return: {LSV_ID : junctions}
    """
    LOG.info("Getting junction coords from voila text file %s" % tsv_file)
    # Just read the lsv ids and junction coordinate column into memory
    deltapsi_txt = import_dpsi_pandas(tsv_file, columns=[2, 16])
    # column 0 is LSV IDs
    lsv_ids = deltapsi_txt[deltapsi_txt.columns[0]]
    # convert to list
    lsv_ids = [x for x in lsv_ids]
    # junction_list = deltapsi_txt[:, 16]
    # column 1 is junction coordinates
    junction_list = deltapsi_txt[deltapsi_txt.columns[1]]
    junctions = dict()
    for ii in range(len(lsv_ids)):
        junctions[lsv_ids[ii]] = junction_list[ii].split(';')
    LOG.info("done")
    return junctions


def collapse_matrix(matrix):
    """Collapse the diagonals probabilities in 1-D and return them"""
    collapse = []
    matrix_corner = matrix.shape[0]
    for i in range(-matrix_corner + 1, matrix_corner):
        collapse.append(np.diagonal(matrix, offset=i).sum())
    return np.array(collapse)


def find_delta_border(v, numbins):
    """Finds the border index to which a V corresponds in its delta_space given
    the number of bins the matrix will have"""
    delta_space = list(np.linspace(-1, 1, num=numbins + 1))
    delta_space.pop(0)  # first border to the left is -1, and we are not interested in it
    # get the index position that corresponds to the V threshold
    for i, value in enumerate(delta_space):
        if value > v:
            return i
    # if nothing hit, V = 1
    return numbins


def matrix_area(matrix, thresh=0.2, absolute=True, collapsed_mat=True):
    """Returns the probability of an event to be above a certain threshold.
    The absolute flag describes if the value is absolute"""
    collapse = matrix
    if not collapsed_mat:
        collapse = collapse_matrix(matrix)
    # get the delta psi histogram borders based on the size of 'collapse'
    border = find_delta_border(thresh, collapse.shape[0])
    # grab the values inside the area of interest
    area = []
    if thresh < 0:
        area.append(collapse[0:border + 1])
        if absolute:  # if absolute V, pick the other side of the array
            area.append(collapse[-border - 1:])
    else:
        area.append(collapse[border:])
        if absolute:  # if absolute V, pick the other side of the array
            area.append(collapse[0:len(collapse) - border])
    return np.sum(area)


def expected_dpsi(bins):
    return sum(np.array(bins) * np.arange(-1 + 1. / len(bins), 1., 2. / len(bins)))


def get_file_voila(file_name):
    print('loading hdf5 file: ', file_name)
    with Voila(file_name, 'r') as v:
        print('done')
        return v.get_voila_lsvs()


def get_file(file_name):
    print('loading pickle file: ', file_name)
    pkl_file = pkl.load(open(file_name, 'rb'), encoding='bytes')
    print('done')
    return pkl_file
