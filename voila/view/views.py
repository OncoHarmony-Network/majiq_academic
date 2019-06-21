from operator import itemgetter
from urllib.parse import urlencode

from flask import jsonify, redirect, url_for, session, render_template, request
from waitress import serve

from voila import constants
from voila.api.view_splice_graph import ViewSpliceGraph
from voila.config import ViewConfig
from voila.exceptions import UnknownAnalysisType
from voila.index import Index
from voila.view import deltapsi, heterogen, psi, splicegraph
import os



if os.name != 'nt':
    import gunicorn.app.base
    from gunicorn.six import iteritems

    class GunicornStandaloneApplication(gunicorn.app.base.BaseApplication):

        def __init__(self, app, options=None):
            self.options = options or {}
            self.application = app
            super(GunicornStandaloneApplication, self).__init__()

        def load_config(self):
            config = dict([(key, value) for key, value in iteritems(self.options)
                           if key in self.cfg.settings and value is not None])
            for key, value in iteritems(config):
                self.cfg.set(key.lower(), value)

        def load(self):
            return self.application










def run_service():
    port = ViewConfig().port
    host = ViewConfig().host
    run_app = get_app()
    web_server = ViewConfig().web_server

    if web_server == 'waitress':
        serve(run_app, port=port, host=host)
    elif web_server == 'gunicorn':
        if os.name == 'nt':
            raise Exception("Gunicorn is unsupported on windows")
        options = {
            'bind': '%s:%s' % (host, port),
            'workers': ViewConfig().num_web_workers,
            'timeout': 9999999
        }
        GunicornStandaloneApplication(run_app, options).run()
    elif web_server == 'flask':
        run_app.run(host=host, port=port, debug=True)
    else:
        raise Exception("Unsupported web server %s specified" % web_server)



def get_app():
    Index()
    analysis_type = ViewConfig().analysis_type

    if not analysis_type:
        run_app = splicegraph.app

    elif analysis_type == constants.ANALYSIS_PSI:
        run_app = psi.app

    elif analysis_type == constants.ANALYSIS_DELTAPSI:
        run_app = deltapsi.app

    elif analysis_type == constants.ANALYSIS_HETEROGEN:
        run_app = heterogen.app

    else:
        raise UnknownAnalysisType(analysis_type)

    return run_app


def copy_lsv(lsv_id, view_matrix, voila_file=None):
    with ViewSpliceGraph() as sg, view_matrix(voila_file=voila_file) as m:
        lsv = m.lsv(lsv_id)
        gene_id = lsv.gene_id
        gene = sg.gene(gene_id)
        lsv_junctions = lsv.junctions.tolist()
        lsv_exons = sg.lsv_exons(gene_id, lsv_junctions)

        juncs = list(j for j in sg.junctions(gene_id) if [j['start'], j['end']] in lsv_junctions)
        exons = list(e for e in sg.exons(gene_id) if (e['start'], e['end']) in lsv_exons)

        lsv_type = lsv.lsv_type
        strand, chromosome, name = itemgetter('strand', 'chromosome', 'name')(gene)
        genome = sg.genome

        if request.get_json():
            sample_name = request.get_json().get('group_name', None)
            sample_names = [sample_name]
        else:
            sample_names = m.group_names

        group_bins = dict(lsv.group_bins)
        bins = [group_bins[sample_name] for sample_name in sample_names]

        splice_graphs = []

        for exp_names in m.experiment_names:
            junc_dict = {
                'junctions': [],
                'exons': exons
            }
            for junc in juncs:
                junction_reads = list(sg.junction_reads_exp(junc, exp_names))
                reads = sum(r['reads'] for r in junction_reads)
                junc['reads'] = reads
                junc_dict['junctions'].append(junc)

            splice_graphs.append(junc_dict)

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


def ucsc_href(genome, chromosome, start, end):
    query_string = {
        'db': genome,
        'position': chromosome + ':' + str(start) + '-' + str(end)
    }

    return 'http://genome.ucsc.edu/cgi-bin/hgTracks?' + urlencode(query_string)


def lsv_boundries(lsv_exons):
    lsv_exons = list(e if e[1] != -1 else (e[0], e[0] + 10) for e in lsv_exons)
    lsv_exons = list(e if e[0] != -1 else (e[1] - 10, e[1]) for e in lsv_exons)
    start = max(e for es in lsv_exons for e in es if e != -1)
    end = min(e for es in lsv_exons for e in es if e != -1)
    return start, end


def gene_view(summary_template, gene_id, view_matrix, **kwargs):

    with ViewSpliceGraph() as sg:
        gene = sg.gene(gene_id)

        # For this gene, remove any already selected highlight/weighted lsvs from session.
        highlight = session.get('highlight', {})
        lsv_ids = [h for h in highlight if h.startswith(gene_id)]
        for lsv_id in lsv_ids:
            del highlight[lsv_id]
        session['highlight'] = highlight

        exons = list(sg.exons(gene_id))
        if not exons:
            return redirect(url_for('index'))
        start = min(e['start'] for e in exons if e['start'] != -1)
        end = max(e['end'] for e in exons if e['end'] != -1)
        href = ucsc_href(sg.genome, gene['chromosome'], start, end)

        gene.update({
            'start': start,
            'end': end,
            'href': href,
            'overlap': sg.gene_overlap(gene_id)
        })

        kwargs['gene'] = gene

        return render_template(summary_template, **kwargs)


def find_exon_number(exons, ref_exon, strand):
    exons = filter(lambda e: -1 not in [e['start'], e['end']], exons)
    exons = list(exons)

    for idx, exon in enumerate(exons):

        if (exon['start'], exon['end']) == ref_exon:
            if strand == '-':
                return len(exons) - idx
            else:
                return idx + 1
    return 'unk'


if __name__ == '__main__':
    app = get_app()
    app.config.update(
        DEBUG=True,
        TEMPLATES_AUTO_RELOAD=True,
        ENV='development'
    )
    app.run()
