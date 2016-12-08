from collections import OrderedDict
from os import path

from voila import constants
from voila.splice_graphics import SpliceGraph
from voila.utils.run_voila_utils import get_env, get_output_html, copy_static
from voila.utils.voila_log import voila_log


def splice_graphs(args):
    voila_log().info("Loading %s." % args.splice_graph)
    majiq_output = {
        'gene_dicts': parse_gene_graphics_obj(args.splice_graph, args.gene_names, args.limit)[0]
    }
    render_summary(args, majiq_output)


def render_summary(args, majiq_output):
    count_pages = 0
    log = voila_log()
    output_html = get_output_html(args, args.splice_graph)
    output_dir = args.output

    env = get_env()
    template_file_name = args.type_analysis.replace("-", "_") + "_summary_template.html"
    sum_template = env.get_template(template_file_name)

    genes = majiq_output['gene_dicts'][majiq_output['gene_dicts'].keys()[0]].values()
    log.info("Number of genes detected in Voila: %d." % len(genes))

    while count_pages * constants.MAX_GENES < len(genes):
        prev_page = None
        next_page = None

        log.info("Processing %d out of %d genes ..." % (
            min((count_pages + 1) * constants.MAX_GENES, len(genes)), len(genes)))

        if (count_pages + 1) * constants.MAX_GENES < len(genes):
            next_page = str(count_pages + 1) + "_" + output_html

        if not count_pages == 0:
            prev_page = str(count_pages - 1) + "_" + output_html

        name_page = '{0}_{1}'.format(count_pages, output_html)

        with open(path.join(output_dir, name_page), 'w') as voila_output:
            voila_output.write(sum_template.render(
                genes=genes[count_pages * constants.MAX_GENES:(count_pages + 1) * constants.MAX_GENES],
                gene_dicts=majiq_output['gene_dicts'],
                first_sample=[sam for sam in majiq_output['gene_dicts'].keys()][0],
                prevPage=prev_page,
                nextPage=next_page,
                namePage=name_page
            ))

        count_pages += 1

    copy_static(args)


def parse_gene_graphics_obj(splice_graph, gene_names_list=None, limit=None):
    genes_exp1_exp2 = []
    log = voila_log()
    genes_exp = {}

    log.info("Parsing splice graph information files ...")

    with SpliceGraph(splice_graph, 'r') as sg:
        genes = sg.get_genes_list(limit=limit)
        experiments = sg.get_experiments_list()

    gene_names_list = [g.lower() for g in gene_names_list]

    for gene in genes:
        for experiment_number, experiment_name in enumerate(experiments):

            if not gene_names_list or gene.name.lower() in gene_names_list:
                experiment_dict = gene.get_experiment(experiment_number)

                try:
                    genes_exp[experiment_name][gene.gene_id] = experiment_dict
                except KeyError:
                    genes_exp[experiment_name] = {gene.gene_id: experiment_dict}

        # todo: there will only be one "exp"... this needs to be refactored.
        genes_exp1_exp2.append(OrderedDict(sorted(genes_exp.items(), key=lambda t: t[0])))

    log.info("Splice graph information files correctly loaded.")
    return genes_exp1_exp2
