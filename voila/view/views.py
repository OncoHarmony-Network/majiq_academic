import gunicorn.app.base
from flask import jsonify
from gunicorn.six import iteritems

from voila import constants
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.config import ViewConfig
from voila.exceptions import UnknownAnalysisType
from voila.index import Index
from voila.view import deltapsi, heterogen, psi, splicegraph


def run_service():
    Index()
    config = ViewConfig()
    analysis_type = config.analysis_type

    options = {
        'bind': '127.0.0.1:' + str(config.port),
        'workers': number_of_workers(),
        'threads': number_of_threads(),
        'worker_class': 'gthread'
    }

    if not analysis_type:
        StandaloneApplication(splicegraph.app, options).run()

    elif analysis_type == constants.ANALYSIS_PSI:
        StandaloneApplication(psi.app, options).run()

    elif analysis_type == constants.ANALYSIS_DELTAPSI:
        StandaloneApplication(deltapsi.app, options).run()

    elif analysis_type == constants.ANALYSIS_HETEROGEN:
        StandaloneApplication(heterogen.app, options).run()

    else:
        raise UnknownAnalysisType(analysis_type)


def number_of_workers():
    return (ViewConfig().nproc * 2) + 1


def number_of_threads():
    return ViewConfig().nproc * 2


class StandaloneApplication(gunicorn.app.base.BaseApplication):

    def init(self, parser, opts, args):
        raise NotImplementedError()

    def __init__(self, application, options=None):
        self.options = options or {}
        self.application = application
        super(StandaloneApplication, self).__init__()

    def load_config(self):
        config = dict([(key, value) for key, value in iteritems(self.options)
                       if key in self.cfg.settings and value is not None])
        for key, value in iteritems(config):
            self.cfg.set(key.lower(), value)

    def load(self):
        return self.application


def copy_lsv(lsv_id, view_matrix):
    with ViewSpliceGraph() as sg, view_matrix() as m:
        dpsi = m.lsv(lsv_id)
        gene = sg.gene(dpsi.gene_id)
        lsv_junctions = dpsi.junctions.tolist()
        lsv_exons = sg.lsv_exons(gene, lsv_junctions)

        juncs = list(j for j in sg.junctions(gene) if [j.start, j.end] in lsv_junctions)
        exons = list(e._asdict() for e in sg.exons(gene) if (e.start, e.end) in lsv_exons)

        lsv_type = dpsi.lsv_type
        strand = gene.strand
        chromosome = gene.chromosome
        name = gene.name
        genome = sg.genome

        group_bins = dict(dpsi.group_bins)
        sample_names = m.group_names
        bins = [group_bins[sample_name] for sample_name in sample_names]

        splice_graphs = []

        for exp_names in m.experiment_names:
            junc_dict = {
                'junctions': [],
                'exons': exons
            }
            for junc in juncs:
                junction_reads = list(sg.junction_reads_exp(junc, exp_names))
                reads = sum(r.reads for r in junction_reads)
                junc = junc._asdict()
                junc['reads'] = reads
                junc_dict['junctions'].append(junc)

            splice_graphs.append(junc_dict)

        juncs = [j._asdict() for j in juncs]

    return jsonify({
        'sample_names': sample_names,
        'genome': genome,
        'lsv_text_version': constants.LSV_TEXT_VERSION,
        'splice_graphs': splice_graphs,
        'lsv': {
            'exons': exons,
            'junctions': juncs,
            'lsv_type': lsv_type,
            'strand': strand,
            'chromosome': chromosome,
            'lsv_id': lsv_id,
            'name': name,
            'bins': bins
        }
    })
