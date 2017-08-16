import pdb
import os

from voila.tools import Tool
from voila.tools.utils import io_caleb
from voila.utils.voila_log import voila_log

# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'

LOG = voila_log()


class ThisisLookup(Tool):
    help = 'Given appropriate directories and a Gene Name or ID, write bash script to generate' \
           ' voila htmls. Automatically generates htmls for each voila threshold you ran.'

    def arguments(self):
        parser = self.get_parser()
        parser.add_argument('indir',
                            type=str,
                            help='Directory where majiq was run.')
        parser.add_argument('gene',
                            type=str,
                            help='Gene Name or Gene IDs to lookup. Maybe be more than 1 separated by commas'
                                 ' without any spaces. Case-insensitive. (e.g. gapdh,RUNX1,Ptprc)')
        parser.add_argument('outdir',
                            type=str,
                            help='Directory to write bash script to.')
        help_mes = 'Optional pattern matching to identify the voila text files'
        parser.add_argument('--pattern',
                            default="*tsv",
                            type=str,
                            help=help_mes)
        help_mes = 'Which comparisons or samples to lookup ID in? Single space or comma separated please.'
        parser.add_argument('--names',
                            '--comparisons',
                            type=str,
                            help=help_mes)
        help_mes = 'Set this if default location is wrong. Location relative to indir.'
        parser.add_argument('--majiq_dpsi_dir',
                            default="majiq/dpsi/",
                            type=str,
                            help=help_mes)
        help_mes = 'Set this if default location is wrong. Location relative to indir.'
        parser.add_argument('--splicegraph_dir',
                            default="build/",
                            type=str,
                            help=help_mes)
        help_mes = 'Set this if default joiner is wrong.'
        parser.add_argument('--comparison_joiner',
                            default="_",
                            type=str,
                            help=help_mes)
        return parser

    def run(self, args):
        # parse the comparisons argument
        if args.names:
            if "," in args.names or " " in args.names:
                args.names.replace(" ", ",")
                to_lookup = args.names.split(",")
            else:
                to_lookup = [args.names]
        else:
            to_lookup = None
        imported = io_caleb.quick_import(input=args.indir,
                                         cutoff_d_psi=0,
                                         cutoff_prob=0,
                                         pattern=args.pattern,
                                         keep_ir=True,
                                         just_one=False,
                                         stop_at=args.gene,
                                         comparisons=to_lookup)
        bash_lines, filepath = write_voila_bash(imported,
                                                args.gene,
                                                args.indir,
                                                args.outdir,
                                                deltapsi_voila_loc=args.majiq_dpsi_dir,
                                                splicegraph_loc=args.splicegraph_dir,
                                                comparisons=to_lookup,
                                                comp_joiner=args.comparison_joiner)
        LOG.info("Printed to this file: %s run lines:\n:%s" % (filepath, bash_lines))


def write_voila_bash(data,
                     gene,
                     indir,
                     outdir,
                     deltapsi_voila_loc="majiq/dpsi/",
                     splicegraph_loc="build/",
                     comparisons=False,
                     overwrite=True,
                     comp_joiner="_"):
    """
    :param data: Quick import
    :param gene: Gene name or ID
        string or list of strings
    :param outdir: Path to write file to
    :param deltapsi_voila_loc: location of deltapsi voila file relative to outdir
    :param splicegraph_loc: location of splicegraph file relative to outdir
    :param comparisons: Boolean. If provided as string or
        list of strings only generate voila lines for these comparisons
    :param overwrite: Boolean. Overwrite outfile if it already exists?
    :param comp_joiner: What are the comparison names joined by in directory?
    :return: String that is written to file..
    """
    io_caleb.check_is_quick_import(data)
    from voila.tools.lookup import lookup
    if comparisons:
        if isinstance(comparisons, str):
            comparison_list = [comparisons]
        elif isinstance(comparisons, list) and isinstance(comparisons[0], str):
            comparison_list = comparisons
        else:
            raise ValueError("Comparison must be a string or list of strings")
    else:
        comparison_list = data.keys()
    if isinstance(gene, str):
        gene_list = [gene]
    elif isinstance(gene, list) and isinstance(gene[0], str):
        gene_list = gene
    else:
        raise ValueError("gene must be a string or list of strings")
    for thisgene in gene_list:
        if "ENSG" in thisgene or "ENSMUSG" in thisgene:
            gene_list.remove(thisgene)
            gene_list.append(io_caleb.genename_from_id(data, thisgene))
    gene_list_join = "_".join(gene_list)
    filename = "html_gen_" + gene_list_join
    if not os.path.exists(outdir):
        raise ValueError("%s path doesn't exist.." % (outdir))
    out_file = os.path.join(outdir, filename + ".sh")
    if overwrite and os.path.exists(out_file):
        LOG.info("Warning! Over-writing existing file...")
    elif not overwrite and os.path.exists(out_file):
        raise RuntimeError("Filename already exists and you don't want to overwrite...")
    voila_outdir = os.path.join(outdir, gene_list_join + "_htmls/")
    deltapsi_voila_loc_abs = os.path.join(indir, deltapsi_voila_loc)
    splicegraph_loc_abs = os.path.join(indir, splicegraph_loc, "splicegraph.hdf5")
    run_lines = []
    with open(out_file, "w") as handle:
        # Don't need to do this anymore..
        #handle.write('source /opt/venv/majiq_hdf5/bin/activate')
        good_to_go = False
        good_comps = list()
        good_p_threshs = list()
        for lsv_dict_name in data.keys():
            if comparisons:
                if io_caleb.comp_without_dup(lsv_dict_name) not in \
                        [io_caleb.comp_without_dup(x) for x in comparison_list]:
                    continue
            found_data = lookup(data[lsv_dict_name],
                                name=gene,
                                printable=True,
                                not_found_error=False)
            if found_data != "gene_not_found" and found_data != "lsv_id_not_found":
                good_to_go = True
                good_comps.append(lsv_dict_name)
                good_p_threshs = io_caleb.get_prob_threshold(data[lsv_dict_name])
            else:
                LOG.info("%s not found in %s" % (gene, lsv_dict_name))
        if not good_to_go:
            raise RuntimeError("None of your genes found in the data ...")
        for comp, prob in zip(good_comps, good_p_threshs):
            cname = io_caleb.gen_comparison_name(data[comp], comp_joiner)
            runline = 'voila deltapsi '
            pdb.set_trace()
            this_deltapsi_voila_loc_abs = os.path.join(deltapsi_voila_loc_abs,
                                                       cname,
                                                       io_caleb.comp_without_dup(comp) + ".deltapsi.voila")
            if not os.path.exists(this_deltapsi_voila_loc_abs):
                raise RuntimeError("Couldn't find %s deltapsi voila file..." % (this_deltapsi_voila_loc_abs))
            runline += this_deltapsi_voila_loc_abs
            runline += " --show-all --no-tsv --gene-names %s --splice-graph " % (" ".join(gene_list))
            runline += splicegraph_loc_abs
            runline += " --threshold %s" % prob
            voila_outfile = os.path.join(voila_outdir, io_caleb.comp_without_dup(comp))
            prob = str(int(float(prob)*100.0))
            runline += " -o " + voila_outfile + "_prob%s" % prob
            handle.writelines(runline)
            run_lines.append(runline)
    return run_lines, out_file
