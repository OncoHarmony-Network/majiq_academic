from collections import OrderedDict
from os import path

from voila import constants
from voila import io_voila
from voila.constants import LSV_TEXT_VERSION
from voila.io_voila import Voila
from voila.producer_consumer import ProducerConsumer
from voila.utils import utils_voila
from voila.utils.run_voila_utils import parse_gene_graphics, VoilaNoLSVsException, table_marks_set, get_env, \
    get_output_html, grouper, copy_static, get_summary_template, get_prev_next_pages
from voila.utils.utils_voila import create_if_not_exists
from voila.utils.voila_log import voila_log


class Psi(ProducerConsumer):
    def __init__(self, args):
        """
        Render psi output.
        :param args: command line arguments
        :return: None
        """
        super(Psi, self).__init__()

        self.args = args
        self.majiq_output = None

        with Voila(args.voila_file, 'r') as v:
            metainfo = v.get_metainfo()
            voila_lsvs = v.get_voila_lsvs(lsv_types=args.lsv_types,
                                          lsv_ids=args.lsv_ids,
                                          gene_names=args.gene_names)

        self.majiq_output = utils_voila.lsvs_to_gene_dict(voila_lsvs, metainfo)
        majiq_output = self.majiq_output

        gene_ids_list = majiq_output['genes_dict'].keys()

        if not len(voila_lsvs):
            raise VoilaNoLSVsException()

        majiq_output['genes_exp'] = parse_gene_graphics(args.splice_graph, metainfo, gene_ids_list)

        majiq_output['lsv_list'] = [ll for g in majiq_output['genes_dict'].viewvalues() for ll in g]

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
        majiq_output = self.majiq_output
        args = self.args
        output_dir = args.output
        gene_keys = sorted(self.majiq_output['genes_dict'].keys())
        output_html = get_output_html(args, args.voila_file)
        summaries_subfolder = path.join(output_dir, constants.SUMMARIES_SUBFOLDER)
        sum_template = get_summary_template(args, get_env())
        comb_spliceg_cond1 = majiq_output['genes_exp'][0].keys()[0]

        while True:
            page_number, subset_keys = self.queue.get()

            subset_keys = [x for x in subset_keys if x]
            genes_dict = OrderedDict((k, majiq_output['genes_dict'][k]) for k in subset_keys)
            page_name = '{0}_{1}'.format(page_number, output_html)

            full_path = path.join(summaries_subfolder, page_name)

            prev_page, next_page = get_prev_next_pages(page_number, len(gene_keys), output_html)

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
                        comb_spliceg_cond1=comb_spliceg_cond1,
                        lsv_text_version=LSV_TEXT_VERSION
                    )
                )

            for _, value in genes_dict.iteritems():
                self.dict(value[0].name, path.join(constants.SUMMARIES_SUBFOLDER, page_name))

            self.queue.task_done()

    def close(self):
        pass

    def render_html(self):
        """
        Render psi html output.
        :param args: command line arguments
        :param majiq_output: dictionary containing majiq output values
        :return: None
        """

        log = voila_log()
        args = self.args
        majiq_output = self.majiq_output
        gene_keys = sorted(majiq_output['genes_dict'].keys())
        output_dir = args.output
        summaries_subfolder = path.join(output_dir, constants.SUMMARIES_SUBFOLDER)
        create_if_not_exists(summaries_subfolder)
        env = get_env()

        log.info("Number of genes detected in Voila: %d." % len(gene_keys))
        log.info("Number of LSVs detected in Voila: %d." % sum(
            [len(majiq_output['genes_dict'][g]) for g in majiq_output['genes_dict']]))
        log.info("Creating HTML5 with splice graphs summaries ...")

        self.run()
        links_dict = self.get_dict()
        majiq_output['voila_links'] = links_dict
        self.manager_shutdown()

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
