from collections import OrderedDict
from os import path

from voila import constants
from voila.splice_graphics import SpliceGraph
from voila.utils.run_voila_utils import get_env, get_output_html, copy_static, grouper, get_prev_next_pages
from voila.utils.voila_log import voila_log
from voila.view.psi import get_prev_next_pages


def splice_graphs(args):
    """
    Render splice graph htmls.
    :param args: arguments parsed from the
    :return:
    """

    log = voila_log()
    log.info("Loading %s." % args.splice_graph)
    gene_dicts = get_gene_dicts(args.splice_graph, args.gene_names, args.limit)

    output_html = get_output_html(args, args.splice_graph)
    output_dir = args.output

    env = get_env()
    template_file_name = args.type_analysis.replace("-", "_") + "_summary_template.html"
    sum_template = env.get_template(template_file_name)

    genes = gene_dicts[gene_dicts.keys()[0]].values()
    genes_length = len(genes)
    log.info("Number of genes detected in Voila: %d." % genes_length)

    for page_number, subset_genes in enumerate(grouper(genes, constants.MAX_GENES)):
        log.info("Processing %d out of %d genes ..." % (
            min((page_number + 1) * constants.MAX_GENES, genes_length), genes_length))

        subset_genes = [x for x in subset_genes if x]

        prev_page, next_page = get_prev_next_pages(page_number, genes_length, output_html, args.limit)

        name_page = '{0}_{1}'.format(page_number, output_html)

        with open(path.join(output_dir, name_page), 'w') as voila_output:
            voila_output.write(sum_template.render(
                genes=subset_genes,
                gene_dicts=gene_dicts,
                first_sample=[sam for sam in gene_dicts.keys()][0],
                prevPage=prev_page,
                nextPage=next_page,
                namePage=name_page
            ))

    copy_static(args)


def get_gene_dicts(splice_graph, gene_names_list, limit):
    log = voila_log()
    genes_exp = {}

    log.info("Parsing splice graph information files ...")

    with SpliceGraph(splice_graph, 'r') as sg:
        gene_ids = sg.get_gene_ids()
        genes = sg.get_genes_list(gene_ids, limit=limit)
        experiments = sg.get_experiments_list()

    gene_names_list = [g.lower() for g in gene_names_list]

    combined_genes_exp = {}

    log.info('Getting gene experiment data...')

    for experiment_index, experiment_name in enumerate(experiments):
        for gene in genes:

            if not gene_names_list or gene.name.lower() in gene_names_list:
                experiment_dict = gene.get_experiment(experiment_index)

                try:
                    genes_exp[experiment_name][gene.gene_id] = experiment_dict
                except KeyError:
                    genes_exp[experiment_name] = {gene.gene_id: experiment_dict}

                combined_genes_exp[gene.gene_id] = gene.combine(experiment_index,
                                                                combined_genes_exp.get(gene.gene_id, None))

    if len(experiments) > 1:
        genes_exp['Combined'] = {gene_id: combined_genes_exp[gene_id] for gene_id in combined_genes_exp}

    log.info("Splice graph information files correctly loaded.")

    return OrderedDict(sorted(genes_exp.items(), key=lambda t: t[0]))
