from operator import itemgetter
from urllib.parse import urlencode

from flask import jsonify, redirect, url_for, session, render_template
from waitress import serve

from voila import constants
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.config import ViewConfig
from voila.exceptions import UnknownAnalysisType
from voila.index import Index
from voila.view import deltapsi, heterogen, psi, splicegraph


def run_service():
    port = ViewConfig().port
    run_app = get_app()
    serve(run_app, port=port)


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


def copy_lsv(lsv_id, view_matrix):
    with ViewSpliceGraph() as sg, view_matrix() as m:
        dpsi = m.lsv(lsv_id)
        gene = sg.gene(dpsi.gene_id)
        lsv_junctions = dpsi.junctions.tolist()
        lsv_exons = sg.lsv_exons(gene, lsv_junctions)

        juncs = list(j for j in sg.junctions(gene) if [j['start'], j['end']] in lsv_junctions)
        exons = list(e for e in sg.exons(gene) if (e['start'], e['end']) in lsv_exons)

        lsv_type = dpsi.lsv_type
        strand, chromosome, name = itemgetter('strand', 'chromosome', 'name')(gene)
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
    start = max(e for es in lsv_exons for e in es if e != -1)
    end = min(e for es in lsv_exons for e in es if e != -1)
    return start, end


def gene_view(summary_template, gene_id, view_matrix):
    with view_matrix() as m:
        if gene_id not in m.gene_ids:
            return redirect(url_for('index'))

    with ViewSpliceGraph() as sg:
        gene = sg.gene(gene_id)

        # For this gene, remove any already selected highlight/weighted lsvs from session.
        highlight = session.get('highlight', {})
        lsv_ids = [h for h in highlight if h.startswith(gene_id)]
        for lsv_id in lsv_ids:
            del highlight[lsv_id]
        session['highlight'] = highlight

        exons = list(sg.exons(gene))
        start = min(e['start'] for e in exons if e['start'] != -1)
        end = max(e['end'] for e in exons if e['end'] != -1)
        href = ucsc_href(sg.genome, gene['chromosome'], start, end)

        gene.update({
            'start': start,
            'end': end,
            'href': href
        })

        return render_template(summary_template, gene=gene)


if __name__ == '__main__':
    import git

    repo = git.Repo(search_parent_directories=True)
    sha = repo.head.object.hexsha
    print(sha[:7])

    app = get_app()
    app.config.update(
        DEBUG=True,
        TEMPLATES_AUTO_RELOAD=True,
        ENV='development'
    )
    app.run()
