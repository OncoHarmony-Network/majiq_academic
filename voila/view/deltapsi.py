from collections import OrderedDict
from os import path

from voila import constants
from voila import io_voila
from voila.constants import LSV_TEXT_VERSION
from voila.io_voila import Voila
from voila.producer_consumer import ProducerConsumer
from voila.utils import utils_voila
from voila.utils.run_voila_utils import VoilaNoLSVsException, parse_gene_graphics, table_marks_set, get_output_html, \
    grouper, copy_static
from voila.utils.voila_log import voila_log
from voila.view.psi import get_prev_next_pages
from voila.view.splice_graphs import get_env


class Deltapsi(ProducerConsumer):
    def __init__(self, args):
        super(Deltapsi, self).__init__()
        self.args = args
        with Voila(args.majiq_quantifier, 'r') as v:
            metainfo = v.get_metainfo()
            voila_lsvs = v.get_voila_lsvs(lsv_ids=args.lsv_ids,
                                          gene_names=args.gene_names,
                                          lsv_types=args.lsv_types)

        self.majiq_output = utils_voila.lsvs_to_gene_dict(voila_lsvs,
                                                          metainfo,
                                                          threshold=args.threshold,
                                                          show_all=args.show_all)

        majiq_output = self.majiq_output

        gene_ids_list = majiq_output['genes_dict'].keys()

        if not len(voila_lsvs):
            raise VoilaNoLSVsException()

        majiq_output['genes_exp'] = parse_gene_graphics(args.splice_graph, metainfo, gene_ids_list)

        majiq_output['lsv_list'] = [ll['lsv'] for g in majiq_output['genes_dict'].viewvalues() for ll in g]

        if not args.no_html:
            self.render_html()

        if not args.no_tsv:
            io_voila.tab_output(args, majiq_output)

        if args.gtf:
            io_voila.generic_feature_format_txt_files(args, majiq_output)

        if args.gff:
            io_voila.generic_feature_format_txt_files(args, majiq_output, out_gff3=True)

    def __enter__(self):
        pass

    def _producer(self):
        gene_keys = sorted(self.majiq_output['genes_dict'].keys())

        for page_number, subset_keys in enumerate(grouper(gene_keys, constants.MAX_GENES)):
            self.queue.put((page_number, subset_keys))

    def _worker(self):
        pass

    def render_html(self):
        log = voila_log()
        majiq_output = self.majiq_output
        args = self.args
        gene_keys = sorted(majiq_output['genes_dict'].keys())
        output_dir = args.output
        output_html = get_output_html(args, args.majiq_quantifier)
        summaries_subfolder = path.join(output_dir, constants.SUMMARIES_SUBFOLDER)
        utils_voila.create_if_not_exists(summaries_subfolder)
        env = get_env()
        template_file_name = args.type_analysis.replace("-", "_") + "_summary_template.html"
        sum_template = env.get_template(template_file_name)
        comb_spliceg_cond1 = majiq_output['genes_exp'][0].keys()[0]
        comb_spliceg_cond2 = majiq_output['genes_exp'][1].keys()[0]
        links_dict = {}

        log.info("Number of genes detected in Voila: %d." % len(gene_keys))
        log.info("Number of LSVs detected in Voila: %d." % sum(
            [len(majiq_output['genes_dict'][g]) for g in majiq_output['genes_dict']]))
        log.info("Creating HTML5 with splice graphs summaries ...")

        # Generate summary
        for page_number, subset_keys in enumerate(grouper(gene_keys, constants.MAX_GENES)):

            log.info("Processing %d out of %d genes ..." % (
                min((page_number + 1) * constants.MAX_GENES, len(gene_keys)), len(majiq_output['genes_dict'])))

            subset_keys = [x for x in subset_keys if x]

            genes_dict = OrderedDict((k, majiq_output['genes_dict'][k]) for k in subset_keys)

            name_page = '{0}_{1}'.format(page_number, output_html)

            prev_page, next_page = get_prev_next_pages(page_number, len(gene_keys), output_html)

            full_path = path.join(summaries_subfolder, name_page)

            with open(full_path, 'w') as voila_output:
                voila_output.write(
                    sum_template.render(
                        tableMarks=[table_marks_set(len(gene_set)) for gene_set in genes_dict],
                        genes_dict=genes_dict,
                        genes_exps_list=majiq_output['genes_exp'],
                        prevPage=prev_page,
                        nextPage=next_page,
                        namePage=name_page,
                        threshold=args.threshold,
                        lexps=majiq_output['meta_exps'],
                        comb_spliceg_cond1=comb_spliceg_cond1,
                        comb_spliceg_cond2=comb_spliceg_cond2,
                        lsv_text_version=LSV_TEXT_VERSION
                    )
                )

            for glsv_list in genes_dict.values():
                links_dict[glsv_list[0]['lsv'].name] = "%s/%s" % (constants.SUMMARIES_SUBFOLDER, name_page)

        majiq_output['voila_links'] = links_dict

        # Generate index
        log.info("Creating HTML5 index summary ...")
        sum_template = env.get_template("index_delta_summary_template.html")

        with open(path.join(output_dir, 'index.html'), 'w') as voila_output:
            voila_output.write(
                sum_template.render(
                    lsvList=majiq_output['lsv_list'],
                    tableMarks=table_marks_set(len(majiq_output['lsv_list'])),
                    threshold=args.threshold,
                    lexps=majiq_output['meta_exps'],
                    links_dict=links_dict,
                    maxLsvs=constants.MAX_LSVS_DELTAPSI_INDEX
                )
            )

        copy_static(args)

    def close(self):
        pass
