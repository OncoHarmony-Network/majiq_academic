from collections import defaultdict
import os
import numpy as np
from voila import constants
from voila.utils import utils_voila

__author__ = 'abarrera'
import cPickle as pkl


class VoilaInput(object):
    """Standard input interface by experiment used by Voila"""
    def __init__(self, lsvs=(), metainfo=None):
        self.lsvs = lsvs
        self.metainfo = metainfo

    def get_lsvs(self):
        return self.lsvs

    def add_lsv(self, l):
        self.lsvs.append(l)

    def get_metainfo(self):
        """Retrieve experiment metainfo (primarily samples information, expandable)"""
        return self.metainfo

    def samples_metainfo(self):
        """Sample information generator, one by one."""
        for sample_info in self.metainfo:
            yield sample_info


def dump_voila_input(voila_input, target, logger=None, protocol=-1):
    import os
    base_path, name = os.path.split(target)
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    try:
        pkl.dump(voila_input, target, protocol)
    except pkl.PickleError, e:
        if logger:
            logger.error("Dumping file %s in %s:\n\t%s." % (voila_input, target, e.message), exc_info=1)


def load_voila_input(voila_input_file, logger=None):
    try:
        with open(voila_input_file, 'rb') as vif:
            voila_input = pkl.load(vif)
            return voila_input
    except pkl.PickleError, e:
        if logger:
            logger.error("Loading the file %s:\n\t%s." % (voila_input_file, e.message), exc_info=1)


def load_dpairs(pairwise_dir, majiq_output, logger):
    """
    Load pairwise files from MAJIQ analysis.

    :param str pairwise_dir: directory containing pairwise comparisons produced by MAJIQ.
    :param majiq_output: parsed data from majiq.
    :param logger: logger instance.
    :return: list of deltapsi lsvs
    :return: name of condition 1
    :return: name of condition 2
    """
    meta_exps = majiq_output['meta_exps']
    lmajiq_pairs = [[None for i in range(len(meta_exps[1])) ] for j in range(len(meta_exps[0]))]

    lsv_names = majiq_output['genes_dict'].keys()

    group1_name = meta_exps[0][0]['group']
    group2_name = meta_exps[1][0]['group']

    for idx1 in range(len(meta_exps[0])):
        for idx2 in range(len(meta_exps[1])):
            pairwise_file = "%s/%s_%d_%s_%d.deltapsi.pickle" % (pairwise_dir, group1_name, idx1+1, group2_name, idx2+1)
            try:
                lmajiq_pairs[idx1][idx2] = utils_voila.get_lsv_delta_exp_data(pairwise_file,
                                                                              show_all=True,
                                                                              gene_name_list=lsv_names,
                                                                              logger=logger)
            except IOError:
                pass
    return lmajiq_pairs, group1_name, group2_name


def write_tab_output(output_dir, output_html, majiq_output, type_summary, logger=None, pairwise_dir=False, threshold=0.2):
    """
    Create tab-delimited output file summarizing all the LSVs detected and quantified with MAJIQ.

    :param output_dir: output directory for the file.
    :param output_html: name for the output html file used to create a *.txt version.
    :param majiq_output: parsed data from majiq.
    :param type_summary: type of analysis performed.
    :param logger: logger instance.
    :param pairwise_dir: whether pairwise comparisons are included or not.
    :param threshold: minimum change considered as significant (in deltapsi analysis).
    :return: nothing.
    """
    ofile_str = "%s%s.%s" % (output_dir, output_html.rsplit('.html', 1)[0], constants.EXTENSION)
    tlb_categx = {'A5SS': 'prime5', 'A3SS': 'prime3', 'Num. Junctions': 'njuncs', 'Num. Exons': 'nexons', 'ES': 'ES'}

    logger.info("Creating Tab-delimited output file in %s..." % ofile_str)

    if pairwise_dir:
        # In deltapsi, add columns with pairwise comparisons between group members
        logger.info("Load pairwise comparison files from %s..." % pairwise_dir)
        lmajiq_pairs, group1_name, group2_name = load_dpairs(pairwise_dir, majiq_output, logger=logger)

    with open(ofile_str, 'w+') as ofile:
        headers = ['#Gene Name',
                   'Gene ID',
                   'LSV ID',
                   'E(PSI) per LSV junction',
                   'Var(E(PSI)) per LSV junction',
                   'LSV Type',
                   'A5SS',
                   'A3SS',
                   'ES',
                   'Num. Junctions',
                   'Num. Exons',
                   'De Novo Junctions?',
                   'chr',
                   'strand',
                   'Junctions coords',
                   'Exons coords',
                   'Exons Alternative Start',
                   'Exons Alternative End',
                   'IR coords']
        if 'voila_links' in majiq_output.keys():
            headers.append('Voila link')

        if 'delta' in type_summary:
            headers[3] = 'E(Delta(PSI)) per LSV junction'
            headers[4] = 'P(Delta(PSI)>%.2f) per LSV junction' % threshold

            if pairwise_dir:
                for idx1 in range(len(lmajiq_pairs)):
                    for idx2 in range(len(lmajiq_pairs[0])):
                        headers.append("%s_%d_%s_%d" % (group1_name, idx1+1, group2_name, idx2+1))

                exp_names_map = ['#Group names and file names mapping']
                for iexp in range(len(lmajiq_pairs)):
                    exp_names_map.append("#%s_%d=%s" % (group1_name, iexp+1, lmajiq_pairs[0][0]['meta_exps'][0][iexp]['experiment']))
                for iexp in range(len(lmajiq_pairs[0])):
                    exp_names_map.append("#%s_%d=%s" % (group2_name, iexp+1, lmajiq_pairs[0][0]['meta_exps'][1][iexp]['experiment']))
                ofile.write('\n'.join(exp_names_map))
                ofile.write('\n')
                ofile.write('#\n')

        # ofile.write("#Tab-delimited file\n#\n")
        ofile.write(constants.DELIMITER.join(headers))
        ofile.write('\n')

        for gene in majiq_output['genes_dict']:
            for llsv_dict in majiq_output['genes_dict'][gene]:
                llsv = llsv_dict
                if type(llsv_dict) == dict:
                    llsv = llsv_dict['lsv']
                lline = []
                lline.extend([llsv.lsv_graphic.get_name(), gene, llsv.get_id()])
                lexpected = []
                lconfidence = []
                for i, bins in enumerate(llsv.get_bins()):
                    if 'delta' in type_summary:
                        lexpected.append(str(-llsv.get_excl_incl()[i][0] + llsv.get_excl_incl()[i][1]))
                        lconfidence.append(str(utils_voila.get_prob_delta_psi_greater_v(bins, float(lexpected[-1]), threshold)))
                    else:
                        lexpected.append(repr(llsv.get_means()[i]))
                        lconfidence.append(repr(llsv.get_variances()[i]))

                lline.append(';'.join(lexpected))
                lline.append(';'.join(lconfidence))

                lline.append(llsv.get_type())
                lline.append(repr(llsv.get_categories()[tlb_categx['A5SS']]))
                lline.append(repr(llsv.get_categories()[tlb_categx['A3SS']]))
                lline.append(repr(llsv.get_categories()[tlb_categx['ES']]))
                lline.append(repr(llsv.get_categories()[tlb_categx['Num. Junctions']]))
                lline.append(repr(llsv.get_categories()[tlb_categx['Num. Exons']]))
                lline.append(str(int(np.any([junc.get_type() == 1 for junc in llsv.lsv_graphic.get_junctions()]))))

                lline.append(llsv.lsv_graphic.get_chrom())
                lline.append(llsv.lsv_graphic.get_strand())

                lline.append(';'.join(['-'.join(str(c) for c in junc.get_coords()) for junc in llsv.lsv_graphic.get_junctions()]))
                lline.append(';'.join(['-'.join(str(c) for c in exon.get_coords()) for exon in llsv.lsv_graphic.get_exons()]))

                try:
                    lline.append(';'.join(['|'.join([str(c) for c in exon.get_alt_starts()]) for exon in llsv.lsv_graphic.get_exons()]))
                    lline.append(';'.join(['|'.join([str(c) for c in exon.get_alt_ends()]) for exon in llsv.lsv_graphic.get_exons()]))
                except TypeError:
                    pass

                lline.append(';'.join([repr(exon.coords) for exon in llsv.lsv_graphic.get_exons() if exon.intron_retention]))

                if pairwise_dir:
                    llpairwise = []
                    for idx1 in range(len(lmajiq_pairs)):
                        for idx2 in range(len(lmajiq_pairs[0])):
                            lpairwise = []
                            if gene in lmajiq_pairs[idx1][idx2]['genes_dict']:
                                for llsv_tmp in lmajiq_pairs[idx1][idx2]['genes_dict'][gene]:
                                    if llsv_tmp[0].get_id() == llsv.get_id():
                                        lsv_pair = llsv_tmp[0]
                                        break
                                else:
                                    logger.warning("LSV %s present in deltagroup but missing in %s." %
                                                   (llsv.get_id(), "%s_%d_%s_%d" % (group1_name, idx1+1,
                                                                                        group2_name, idx2+1)))
                                    lpairwise.append('N/A')
                                    continue
                                for iway in range(len(llsv.get_bins())):
                                    lpairwise.append(str(sum(lsv_pair.get_excl_incl()[iway])))
                            else:
                                lpairwise.append('N/A')
                            llpairwise.append(';'.join(lpairwise))
                    lline.extend(llpairwise)
                if 'voila_links' in majiq_output.keys():
                    summary_path = majiq_output['voila_links'][llsv.get_gene_name()]
                    if not os.path.isabs(summary_path):
                        summary_path = "%s/%s/%s" % (os.getcwd(), output_dir, summary_path)
                    lline.append(constants.URL_COMPOSITE % (summary_path, llsv.get_gene_name()))
                ofile.write(constants.DELIMITER.join(lline))
                ofile.write('\n')

    logger.info("Delimited output file successfully created in: %s" % ofile_str)


def load_dpsi_tab(tab_files_list, sample_names, thres_change=None):
    """Load LSV delta psi information from tab-delimited file."""
    # #Tab-delimited file
    # #
    # #Gene Name      Gene ID LSV ID  E(Delta(PSI)) per LSV junction  P(Delta(PSI)>0.20) per LSV junction     LSV Type        A5SS    A3SS    ES      Num. Junctions  Num. Exons      De Novo Junctions?      chr     strand  Junctions coords        Exons coords    Exons Alternative Start Exons Alternative End   Voila link
    # SOCS4   ENSG00000180008 ENSG00000180008:55493948-55494189:source        -6.31549134148e-06;-6.21425294136e-06;-0.0014750236202;1.11410403456e-05;0.00251893273428       6.19368325476e-06;6.17099789797e-06;0.0152398891153;8.65809173803e-06;0.0113650379196   s|1e1.1o4|1e1.2o4|1e1.3o4|1e1.4o4|1e2.1o1       False   True    True    5       3       1       chr14   +  55494189-55498522;55494189-55498559;55494189-55498581;55494189-55498607;55494189-55509670        55493948-55494189;55498522-55498712;55498522-55498712;55498522-55498712;55498522-55498712;55509670-55516206     ;;;;;   ;;;;;   file:///data/JennieLin/voila/1298.M1_M2//summaries/783_1298_M1_1298_M2.deltapsi_deltapsi.html#SOCS4
    # DDX20   ENSG00000064703 ENSG00000064703:112299268-112299362:source      1.68483233714e-05;-1.68483233698e-05    7.62215489154e-06;4.49754345183e-06     s|1e1.1o1|1e2.2o2       False   False   True    2       3       0       chr1    +       112299362-112302022;112299362-112303096 112299268-112299362;112302022-112302190;112303019-112303210     ;;      ;;      file:///data/JennieLin/voila/1298.M1_M2//summaries/57_1298_M1_1298_M2.deltapsi_deltapsi.html#DDX20
    lsvs_dict = defaultdict(lambda: defaultdict(lambda: None))
    for idx, tab_file in enumerate(tab_files_list):
        with open(tab_file, 'r') as tabf:
            for line in tabf:
                if line.startswith("#"): continue
                fields = line.split()
                expecs = [float(aa) for aa in fields[3].split(";")]
                if np.max([abs(ee) for ee in expecs]) < thres_change: continue
                if lsvs_dict[fields[2]]['expecs'] is None:
                    lsvs_dict[fields[2]]['expecs'] = [None]*len(sample_names)
                    lsvs_dict[fields[2]]['links'] = [None]*len(sample_names)
                lsvs_dict[fields[2]]['expecs'][idx] = expecs[np.argmax([abs(ee) for ee in expecs])]  #TODO: What happens when the most changing junction is not the same?!
                lsvs_dict[fields[2]]['links'][idx] = fields[-1].split('voila/')[1]
                lsvs_dict[fields[2]]['gene'] = fields[0]
    return lsvs_dict