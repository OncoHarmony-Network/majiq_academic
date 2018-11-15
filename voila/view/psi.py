import os
from bisect import bisect

from flask import render_template, url_for, jsonify, request, session, Flask, Response

from voila.api.view_matrix import ViewPsi
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.view.datatables import DataTables
from voila.view.forms import LsvFiltersForm
from voila.index import Index

app = Flask(__name__)
app.secret_key = os.urandom(16)


@app.route('/')
def index():
    form = LsvFiltersForm()
    return render_template('psi_index.html', form=form)


@app.route('/gene/<gene_id>/')
def gene(gene_id):
    # For this gene, remove any already selected highlight/weighted lsvs from session.
    highlight = session.get('highlight', {})
    lsv_ids = [h for h in highlight if h.startswith(gene_id)]
    for lsv_id in lsv_ids:
        del highlight[lsv_id]
    session['highlight'] = highlight

    return render_template('psi_summary.html', gene_id=gene_id)


@app.route('/index-table', methods=('POST',))
def index_table():
    with ViewPsi() as v:
        grp_name = v.group_names[0]

        dt = DataTables(Index.psi(), ('gene_name', 'lsv_id'))

        for idx, index_row, records in dt.callback():
            gene_name = index_row['gene_name'].decode('utf-8')
            gene_id = index_row['gene_id'].decode('utf-8')
            lsv_id = index_row['lsv_id'].decode('utf-8')

            records[idx] = [
                [url_for('gene', gene_id=gene_id), gene_name],
                lsv_id,
                v.lsv(lsv_id).lsv_type,
                grp_name,
                ''
            ]

        return jsonify(dict(dt))


@app.route('/nav/<gene_id>', methods=('POST',))
def nav(gene_id):
    with ViewPsi() as h:
        gene_ids = list(sorted(h.gene_ids))
        idx = bisect(gene_ids, gene_id)

        return jsonify({
            'next': url_for('gene', gene_id=gene_ids[idx % len(gene_ids)]),
            'prev': url_for('gene', gene_id=gene_ids[(idx % len(gene_ids)) - 2])
        })


@app.route('/metadata', methods=('POST',))
def metadata():
    with ViewPsi() as v:
        return jsonify({
            'group_names': v.group_names,
            'experiment_names': v.experiment_names
        })


@app.route('/splice-graph/<gene_id>', methods=('POST', 'GET'))
def splice_graph(gene_id):
    with ViewSpliceGraph() as sg, ViewPsi() as v:
        g = sg.gene(gene_id)
        exp_names = v.splice_graph_experiment_names
        gd = sg.gene_experiment(g, exp_names)
        gd['experiment_names'] = exp_names
        gd['group_names'] = v.group_names
        return jsonify(gd)


@app.route('/summary-table/<gene_id>', methods=('POST',))
def summary_table(gene_id):
    with ViewPsi() as v:
        grp_name = v.group_names[0]
        index_data = Index.psi(gene_id)

        dt = DataTables(index_data, ('highlight', 'lsv_id'), sort=False, slice=False)

        dt.add_sort('highlight', DataTables.highlight)
        dt.add_sort('lsv_id', DataTables.lsv_id)

        dt.sort()
        dt.slice()

        for idx, record, records in dt.callback():
            lsv_id = record['lsv_id'].decode('utf-8')
            psi = v.lsv(lsv_id)
            lsv_type = psi.lsv_type
            try:
                highlight = session['highlight'][lsv_id]
            except KeyError:
                highlight = [False, False]

            records[idx] = [
                highlight,
                lsv_id,
                lsv_type,
                grp_name,
                ''
            ]

        return jsonify(dict(dt))


@app.route('/psi-splice-graphs', methods=('POST',))
def psi_splice_graphs():
    with ViewPsi() as v:
        try:
            sg_init = session['psi_init_splice_graphs']
        except KeyError:
            sg_init = [[v.group_names[0], v.splice_graph_experiment_names[0][0]]]

        json_data = request.get_json()

        if json_data:
            if 'add' in json_data:
                if all(s != json_data['add'] for s in sg_init):
                    sg_init.append(json_data['add'])

            if 'remove' in json_data:
                sg_init = filter(lambda s: s != json_data['remove'], sg_init)
                sg_init = list(sg_init)

        session['psi_init_splice_graphs'] = sg_init

        return jsonify(sg_init)


@app.route('/lsv-data', methods=('POST',))
@app.route('/lsv-data/<lsv_id>', methods=('POST',))
def lsv_data(lsv_id):
    gene_id = ':'.join(lsv_id.split(':')[:-2])
    ref_exon = list(map(int, lsv_id.split(':')[-1].split('-')))

    def find_exon_number(exons):
        exons = filter(lambda e: -1 not in [e.start, e.end], exons)
        exons = list(exons)

        for idx, exon in enumerate(exons):
            if [exon.start, exon.end] == ref_exon:
                if strand == '-':
                    return len(exons) - idx
                else:
                    return idx + 1

    with ViewSpliceGraph() as sg, ViewPsi() as m:
        gene = sg.gene(gene_id)
        strand = gene.strand
        exons = sg.exons(gene)
        exon_number = find_exon_number(exons)

        lsv = m.lsv(lsv_id)

        return jsonify({
            'lsv': {
                'name': m.group_names[0],
                'junctions': lsv.junctions.tolist(),
                'group_means': dict(lsv.group_means),
                'group_bins': dict(lsv.group_bins)
            },
            'exon_number': exon_number
        })


@app.route('/lsv-highlight', methods=('POST',))
def lsv_highlight():
    json_data = request.get_json()

    with ViewPsi() as m:

        lsvs = []
        highlight_dict = session.get('highlight', {})

        for lsv_id, highlight, weighted in json_data:
            highlight_dict[lsv_id] = [highlight, weighted]

        session['highlight'] = highlight_dict

        for lsv_id, (highlight, weighted) in highlight_dict.items():
            if highlight:
                psi = m.lsv(lsv_id)
                junctions = psi.junctions.tolist()

                if psi.lsv_type[-1] == 'i':
                    intron_retention = junctions[-1]
                    junctions = junctions[:-1]
                else:
                    intron_retention = []

                lsvs.append({
                    'junctions': junctions,
                    'intron_retention': intron_retention,
                    'reference_exon': list(psi.reference_exon),
                    'weighted': weighted,
                    'group_means': dict(psi.group_means)

                })

        return jsonify(lsvs)


@app.route('/download-lsvs', methods=('POST',))
def download_lsvs():
    dt = DataTables(Index.psi(), ('gene_name', 'lsv_id'), slice=False)

    data = (d['lsv_id'].decode('utf-8') for d in dict(dt)['data'])
    data = '\n'.join(data)

    return Response(data, mimetype='text/plain')


@app.route('/download-genes', methods=('POST',))
def download_genes():
    dt = DataTables(Index.psi(), ('gene_name', 'lsv_id'), slice=False)

    data = set(d['gene_id'].decode('utf-8') for d in dict(dt)['data'])
    data = '\n'.join(data)

    return Response(data, mimetype='text/plain')
