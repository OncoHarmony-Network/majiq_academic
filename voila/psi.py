from collections import OrderedDict
from math import ceil
from os import path

from voila import constants
from voila import io_voila
from voila.io_voila import Voila
from voila.utils import utils_voila
from voila.utils.run_voila_utils import parse_gene_graphics, VoilaNoLSVsException, table_marks_set, get_env, \
    get_output_html, \
    grouper, copy_static, get_summary_template
from voila.utils.utils_voila import create_if_not_exists
from voila.utils.voila_log import voila_log


def psi(args):
    with Voila(args.majiq_quantifier, 'r') as v:
        voila_lsvs = v.get_voila_lsvs(lsv_types=args.lsv_types,
                                      lsv_ids=args.lsv_ids,
                                      gene_names=args.gene_names)
        metainfo = v.get_metainfo()

    majiq_output = utils_voila.lsvs_to_gene_dict(voila_lsvs, metainfo)

    gene_ids_list = majiq_output['genes_dict'].keys()

    if not len(voila_lsvs):
        raise VoilaNoLSVsException()

    # Get gene info
    majiq_output['genes_exp'] = parse_gene_graphics(args.splice_graph, metainfo, gene_ids_list)

    majiq_output['lsv_list'] = [ll for g in majiq_output['genes_dict'].viewvalues() for ll in g]

    if not args.no_html:
        render_html(args, majiq_output)

    if not args.no_tsv:
        io_voila.tab_output(args, majiq_output)

    if args.gtf:
        io_voila.generic_feature_format_txt_files(args)

    if args.gff:
        io_voila.generic_feature_format_txt_files(args, out_gff3=True)


def render_html(args, majiq_output):
    log = voila_log()
    gene_keys = sorted(majiq_output['genes_dict'].keys())
    log.info("Number of genes detected in Voila: %d." % len(gene_keys))
    log.info("Number of LSVs detected in Voila: %d." % sum(
        [len(majiq_output['genes_dict'][g]) for g in majiq_output['genes_dict']]))
    log.info("Creating HTML5 with splice graphs summaries ...")

    output_dir = args.output
    output_html = get_output_html(args, args.majiq_quantifier)
    summaries_subfolder = path.join(output_dir, constants.SUMMARIES_SUBFOLDER)
    create_if_not_exists(summaries_subfolder)

    env = get_env()
    sum_template = get_summary_template(args, env)

    last_page = ceil(len(majiq_output['genes_dict']) / constants.MAX_GENES)
    comb_spliceg_cond1 = majiq_output['genes_exp'][0].keys()[0]

    links_dict = {}

    for page_number, subset_keys in enumerate(grouper(gene_keys, constants.MAX_GENES)):

        log.info("Processing %d out of %d genes ..." % (
            min((page_number + 1) * constants.MAX_GENES, len(gene_keys)), len(majiq_output['genes_dict'])))

        subset_keys = [x for x in subset_keys if x]
        genes_dict = OrderedDict((k, majiq_output['genes_dict'][k]) for k in subset_keys)
        page_name = '{0}_{1}'.format(page_number, output_html)

        prev_page = None
        next_page = None

        full_path = path.join(summaries_subfolder, page_name)

        if page_number != last_page:
            next_page = '{0}_{1}'.format(page_number + 1, output_html)

        if page_number != 0:
            prev_page = '{0}_{1}'.format(page_number - 1, output_html)

        with open(full_path, 'w') as voila_output:
            voila_output.write(
                sum_template.render(
                    tableMarks=[table_marks_set(len(gene_set)) for gene_set in genes_dict],
                    genes_dict=genes_dict,
                    prevPage=prev_page,
                    nextPage=next_page,
                    namePage=page_name,
                    lexps=majiq_output['meta_exps'],
                    genes_exps_list=majiq_output['genes_exp'],
                    comb_spliceg_cond1=comb_spliceg_cond1
                )
            )

        for key in genes_dict:
            links_dict[genes_dict[key][0].name] = path.join(constants.SUMMARIES_SUBFOLDER,
                                                            page_name)

    majiq_output['voila_links'] = links_dict

    # Generate index
    log.info("Creating HTML5 index summary ...")
    sum_template = env.get_template("index_single_summary_template.html")
    with open(path.join(output_dir, 'index.html'), 'w') as voila_output:
        voila_output.write(
            sum_template.render(
                lsvList=majiq_output['lsv_list'],
                tableMarks=table_marks_set(len(majiq_output['lsv_list'])),
                lexps=majiq_output['meta_exps'],
                links_dict=links_dict,
                maxLsvs=constants.MAX_LSVS_PSI_INDEX
            )
        )

    copy_static(args)
