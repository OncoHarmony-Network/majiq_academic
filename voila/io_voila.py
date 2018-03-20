from voila.utils import utils_voila

__author__ = 'abarrera'


def load_dpairs(pairwise_dir, majiq_output):
    """
    Load pairwise files from MAJIQ analysis.

    :param str pairwise_dir: directory containing pairwise comparisons produced by MAJIQ.
    :param majiq_output: parsed data from old_majiq.
    :return: list of deltapsi lsvs
    :return: name of condition 1
    :return: name of condition 2
    """
    meta_exps = majiq_output['meta_exps']
    lmajiq_pairs = [[None for i in range(len(meta_exps[1]))] for j in range(len(meta_exps[0]))]

    lsv_names = majiq_output['genes_dict'].keys()

    group1_name = meta_exps[0][0]['group']
    group2_name = meta_exps[1][0]['group']

    for idx1 in range(len(meta_exps[0])):
        for idx2 in range(len(meta_exps[1])):
            pairwise_file = "%s/%s_%d_%s_%d.deltapsi.pickle" % (
                pairwise_dir, group1_name, idx1 + 1, group2_name, idx2 + 1)
            try:
                lmajiq_pairs[idx1][idx2] = utils_voila.get_lsv_delta_exp_data(
                    pairwise_file,
                    show_all=True,
                    gene_name_list=lsv_names
                )
            except IOError:
                pass
    return lmajiq_pairs, group1_name, group2_name

# def generic_feature_format_txt_files(args, out_gff3=False):
#     """
#     Create GFF3 files for each LSV.
#     :param majiq_output: majiq data
#     :param args: parsed input data
#     :param out_gff3: output as a GFF3 file
#     :return: None
#     """
#
#     log = voila_log()
#     output_dir = args.output
#
#     if out_gff3:
#         log.info("Create GFF files for LSVs")
#     else:
#         log.info("Create GTF files for LSVs")
#
#     header = "##gff-version 3"
#
#     odir = join(output_dir, "static/doc/lsvs")
#     utils_voila.create_if_not_exists(odir)
#
#     with Voila(args.voila_file, 'r') as v:
#         for lsv in v.get_voila_lsvs(args):
#             lsv_file_basename = "%s/%s" % (odir, lsv.lsv_id)
#             try:
#                 lsv_gff3_str = lsv.to_gff3()
#                 utils_voila.gff2gtf(lsv_gff3_str.split('\n'), "%s.gtf" % lsv_file_basename)
#
#                 # not accessible from command line
#                 if out_gff3:
#                     gff_file = "%s.gff3" % (lsv_file_basename)
#                     with open(gff_file, 'w') as ofile:
#                         ofile.write(header + "\n")
#                         ofile.write(lsv_gff3_str + "\n")
#
#             except UnboundLocalError as e:
#                 log.warning("problem generating GTF file for %s" % lsv.id)
#                 log.error(e)
