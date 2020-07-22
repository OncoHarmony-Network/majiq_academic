from os.path import join

from rna_voila.utils import utils_voila
from rna_voila.voila_log import voila_log

__author__ = 'abarrera'


def generic_feature_format_txt_files(args, out_gff3=False):
    """
    This is legecy code that might be needed in the future to write GTFs. There is a issue to implement this feature.

    Create GFF3 files for each LSV.

    :param majiq_output: majiq data
    :param args: parsed input data
    :param out_gff3: output as a GFF3 file
    :return: None
    """

    log = voila_log()
    output_dir = args.output

    if out_gff3:
        log.info("Create GFF files for LSVs")
    else:
        log.info("Create GTF files for LSVs")

    header = "##gff-version 3"

    odir = join(output_dir, "static/doc/lsvs")
    utils_voila.create_if_not_exists(odir)

    with Voila(args.voila_file, 'r') as v:
        for lsv in v.get_voila_lsvs(args):
            lsv_file_basename = "%s/%s" % (odir, lsv.lsv_id)
            try:
                lsv_gff3_str = lsv.to_gff3()
                utils_voila.gff2gtf(lsv_gff3_str.split('\n'), "%s.gtf" % lsv_file_basename)

                # not accessible from command line
                if out_gff3:
                    gff_file = "%s.gff3" % (lsv_file_basename)
                    with open(gff_file, 'w') as ofile:
                        ofile.write(header + "\n")
                        ofile.write(lsv_gff3_str + "\n")

            except UnboundLocalError as e:
                log.warning("problem generating GTF file for %s" % lsv.id)
                log.error(e)
