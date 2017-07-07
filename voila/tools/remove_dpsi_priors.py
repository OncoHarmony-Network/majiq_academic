from voila.tools import Tool
from voila.io_voila import Voila
from voila.tools.utils import find_files
import pickle as pkl
import numpy as np
import pdb
import os

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
        parser.add_argument('comparison_name',
                            type=str,
                            help='Comparison name you want to remove priors for.')
        # TODO : optionally instead of giving directory + comp name, directly give file paths
        # parser.add_argument('--voila_deltapsi_object',
        #                     type=str,
        #                     help='Optional *.deltapsi.voila file path if)
        # parser.add_argument('--prior_matrix',
        #                     type=str,
        #                     help='*.priormatrix.pkl file')
        # parser.add_argument('--voila_txt',
        #                     type=str,
        #                     help='Voila tab output text file')
        parser.add_argument('out_dir',
                            type=str,
                            help='Directory save pickle file results. Name will be same as tsv, but '
                                 'with deltapsi_no_prior instead of deltapsi_deltapsi')
        # help_mes = 'Optional flag: return comparison names ?'
        # parser.add_argument('--return-names',
        #                     action='store_true',
        #                     help=help_mes)

        return parser

    def run(self, args):
        deltapsi_voila_fp = find_files.find_files(path=args.directory,
                                                  pattern="%s.deltapsi.voila" % args.comparison_name,
                                                  recursive=True)
        deltapsi_prior_fp = find_files.find_files(path=args.directory,
                                                  pattern="%s.priormatrix.pkl" % args.comparison_name,
                                                  recursive=True)
        deltapsi_tabfile_fp = find_files.find_files(path=args.directory,
                                                    pattern="%s.deltapsi_deltapsi.tsv" % args.comparison_name,
                                                    recursive=True)
        if len(deltapsi_voila_fp) != 1:
            raise RuntimeError("Didn't find 1 deltapsi_voila file, instead found %s" % len(deltapsi_voila_fp))
        if len(deltapsi_prior_fp) != 1:
            raise RuntimeError("Didn't find 1 deltapsi_prior file, instead found %s" % len(deltapsi_prior_fp))
        if len(deltapsi_tabfile_fp) != 1:
            raise RuntimeError("Didn't find 1 deltapsi_tabfile file, instead found %s" % len(deltapsi_tabfile_fp))
        outname = os.path.basename(deltapsi_tabfile_fp[0])
        LOG.info("removing prior from %s " % outname)
        res = remove_dpsi_priors(deltapsi_voila=deltapsi_voila_fp[0],
                                 deltapsi_prior=deltapsi_prior_fp[0],
                                 deltapsi_tabfile=deltapsi_tabfile_fp[0])
        LOG.info("Saving results to %s " % args.out_dir)
        outname.replace("deltapsi_deltapsi", "deltapsi_no_prior")
        outname.replace("tsv", "pickle")
        outpath = os.path.join(args.out_dir, outname)
        pkl.dump(res, open(outpath, "wb"))
        LOG.info("Finished removing removing prior from %s" % os.path.basename(deltapsi_tabfile_fp[0]))


def remove_dpsi_priors(deltapsi_voila, deltapsi_prior, deltapsi_tabfile):
    """
    :param deltapsi_voila: *.deltapsi.voila
    :param deltapsi_prior: *.priomatrix.pkl
    :param deltapsi_tabfile: *tsv tab output text file of voila
    :return:
    """
    with Voila(deltapsi_voila, 'r') as v:
        print('loading hdf5 file: ', deltapsi_voila)
        tissue_lsvs = v.get_voila_lsvs()
        print("done")
        tissue_priors = np.array(get_file(deltapsi_prior))
        print("Loading voila text file %s" % deltapsi_tabfile)
        deltapsi_txt = np.loadtxt(deltapsi_tabfile, delimiter='\t', dtype=bytes).astype(str)
        print("done")
        lsv_ids = deltapsi_txt[:, 2]
        junction_list = deltapsi_txt[:, 16]
        junctions = dict()
        for ii in range(len(lsv_ids)):
            junctions[lsv_ids[ii]] = junction_list[ii].split(';')

        dist_no_priors_edpsi_all = dict()
        i = 1
        for tissue_lsv in tissue_lsvs:
            print(i)
            dist_no_priors_edpsi = list()
            if tissue_lsv.lsv_id not in lsv_ids:
                continue
            elif i == len(lsv_ids):
                break
            for junction_number in range(len(junctions[tissue_lsv.lsv_id])):
                prior_info = np.log(collapse_matrix(tissue_priors[0, :, :]))
                mle_bins = np.exp(np.log(np.array(tissue_lsv.bins[int(junction_number)])) - prior_info)
                mle_bins /= 1. * np.sum(mle_bins)
                dist_no_priors_edpsi.append(expected_dpsi(mle_bins))
            dist_no_priors_edpsi_all[tissue_lsv.lsv_id] = dist_no_priors_edpsi
            i += 1
    return dist_no_priors_edpsi_all


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
