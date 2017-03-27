import os
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
        self.comb_spliceg_cond1 = None
        self.comb_spliceg_cond2 = None
        self.template_file_name = None
        self.output_html = None
        self.gene_keys_length = None
        self.sum_template = None

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
        majiq_output = self.majiq_output
        gene_keys = sorted(majiq_output['genes_dict'].keys())

        for page_number, subset_keys in enumerate(grouper(gene_keys, constants.MAX_GENES)):
            subset_keys = (x for x in subset_keys if x)
            genes_dict = OrderedDict((k, majiq_output['genes_dict'][k]) for k in subset_keys)
            self.queue.put((page_number, genes_dict))

    def _worker(self):
        majiq_output = self.majiq_output

        while True:
            page_number, genes_dict = self.queue.get()
            page_name = '{0}_{1}'.format(page_number, self.output_html)
            prev_page, next_page = get_prev_next_pages(page_number, self.gene_keys_length, self.output_html)

            with open(path.join(self.args.output, constants.SUMMARIES_SUBFOLDER, page_name), 'w') as voila_output:
                voila_output.write(
                    self.sum_template.render(
                        tableMarks=[table_marks_set(len(gene_set)) for gene_set in genes_dict],
                        genes_dict=genes_dict,
                        genes_exps_list=majiq_output['genes_exp'],
                        prevPage=prev_page,
                        nextPage=next_page,
                        namePage=page_name,
                        threshold=self.args.threshold,
                        lexps=majiq_output['meta_exps'],
                        comb_spliceg_cond1=self.comb_spliceg_cond1,
                        comb_spliceg_cond2=self.comb_spliceg_cond2,
                        lsv_text_version=LSV_TEXT_VERSION
                    )
                )

            for _, value in genes_dict.iteritems():
                self.dict(value[0]['lsv'].name, os.path.join(constants.SUMMARIES_SUBFOLDER, page_name))

            self.queue.task_done()

    def render_html(self):
        args = self.args
        majiq_output = self.majiq_output

        self.comb_spliceg_cond1 = majiq_output['genes_exp'][0].keys()[0]
        self.comb_spliceg_cond2 = majiq_output['genes_exp'][1].keys()[0]
        self.template_file_name = args.type_analysis.replace("-", "_") + "_summary_template.html"
        self.output_html = get_output_html(args, args.majiq_quantifier)
        self.gene_keys_length = len(majiq_output['genes_dict'].keys())
        self.sum_template = get_env().get_template(self.template_file_name)

        utils_voila.create_if_not_exists(os.path.join(args.output, constants.SUMMARIES_SUBFOLDER))

        self.run()

        majiq_output['voila_links'] = self.get_dict()

        # Generate index
        voila_log().info("Creating HTML5 index summary ...")
        sum_template = get_env().get_template("index_delta_summary_template.html")

        with open(path.join(args.output, 'index.html'), 'w') as voila_output:
            voila_output.write(
                sum_template.render(
                    lsvList=majiq_output['lsv_list'],
                    tableMarks=table_marks_set(len(majiq_output['lsv_list'])),
                    threshold=args.threshold,
                    lexps=majiq_output['meta_exps'],
                    links_dict=majiq_output['voila_links'],
                    maxLsvs=constants.MAX_LSVS_DELTAPSI_INDEX
                )
            )

        copy_static(args)

    def close(self):
        pass
