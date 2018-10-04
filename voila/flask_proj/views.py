import os
from bisect import bisect

from flask import Flask, render_template, jsonify, url_for, session, request, redirect
from waitress import serve

from voila.api.view_matrix import ViewPsi
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.config import Config
from voila.flask_proj.datatables import DataTables
from voila.flask_proj.html import Html

app = Flask(__name__)
app.secret_key = os.urandom(16)


def run_service(args):
    serve(app, listen='127.0.0.1:55555')


def get_meta():
    with ViewPsi() as m:
        return m.view_metadata


@app.route('/gene/<gene_id>/')
def gene(gene_id):
    return render_template('psi_summary.html', gene_id=gene_id)


@app.route('/')
def index():
    return render_template('psi_index.html')


@app.route('/index-table', methods=('POST',))
def index_table():
    config = Config()
    meta = get_meta()
    with ViewSpliceGraph() as sg, ViewPsi() as p:
        records = []

        lsv_ids = list(p.lsv_ids())

        for lsv_id in lsv_ids:
            psi = p.psi(lsv_id)
            gene_name = sg.gene(psi.gene_id).name
            records.append([gene_name, psi.lsv_id, psi.lsv_type, 'psi per junc', 'links'])

        def add_html(rs):
            html = Html()

            for r in rs:
                gene_id = ':'.join(r[1].split(':')[0:-2])

                html.tag('a', href=url_for('gene', gene_id=gene_id))
                html.text(r[0])
                r[0] = str(html)
                html.reset()

                html.tag('canvas', classes=['lsv-cartoon'], data_lsv_id=r[1], data_lsv_type=r[2])
                r[2] = str(html)
                html.reset()

                html.tag('canvas',
                         classes=['lsv-single-compact-percentiles'],
                         data_group=meta['group_names'][0],
                         data_lsv_id=r[1]
                         )
                r[3] = str(html)
                html.reset()
                html.tag('svg', classes=['psi-violin-plot'], data_lsv_id=r[1], data_group=meta['group_names'][0])
                r[3] += str(html)
                html.reset()

    dt = DataTables(records, add_html)

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
    return jsonify(get_meta())


@app.route('/splice-graph/<gene_id>', methods=('POST', 'GET'))
def splice_graph(gene_id):
    meta = get_meta()
    if request.method == 'GET':
        return redirect(url_for('index'))

    with ViewSpliceGraph() as sg:
        g = sg.gene(gene_id)
        gd = sg.gene_experiment(g, meta['experiment_names'])
        return jsonify(gd)


@app.route('/psi/<gene_id>', methods=('POST',))
def psi(gene_id):
    records = []
    meta = get_meta()
    with ViewPsi() as m:
        lsv_ids = list(m.lsv_ids(gene_ids=[gene_id]))
        for lsv_id in lsv_ids:
            lsv = m.psi(lsv_id)
            records.append(['', lsv_id, lsv.lsv_type, '', ''])

        def add_html(rs):
            html = Html()
            for r in rs:
                html.tag('canvas', classes=['lsv-cartoon'], data_lsv_id=r[1], data_lsv_type=r[2])
                r[2] = str(html)
                html.reset()

                html.tag('canvas',
                         classes=['lsv-single-compact-percentiles'],
                         data_group=meta['group_names'][0],
                         data_lsv_id=r[1]
                         )
                r[3] = str(html)
                html.reset()
                html.tag('svg', classes=['psi-violin-plot'], data_lsv_id=r[1], data_group=meta['group_names'][0])
                r[3] += str(html)
                html.reset()

        dt = DataTables(records, add_html)

        return jsonify(dict(dt))


@app.route('/psi-splice-graphs', methods=('POST',))
def psi_splice_graphs():
    meta = get_meta()

    try:
        sg_init = session['psi_init_splice_graphs']
    except KeyError:
        sg_init = [[meta['group_names'][0], meta['experiment_names'][0][0]]]

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
@app.route('/lsv-data/<gene_id>', methods=('POST',))
def lsv_data(gene_id=None):
    meta = get_meta()
    with ViewSpliceGraph() as sg, ViewPsi() as m:
        exon_numbers = {}

        if gene_id:
            genes = [sg.gene(gene_id)]
        else:
            genes = sg.genes()

        for gene in genes:
            strand = gene.strand
            exons = sg.exons(gene)
            exons = filter(lambda e: -1 not in [e.start, e.end], exons)
            exons = list(exons)

            for idx, exon in enumerate(exons):
                if strand == '+':
                    value = idx + 1
                else:
                    value = len(exons) - idx

                key = str(exon.start) + '-' + str(exon.end)
                assert key not in exon_numbers
                exon_numbers[key] = value

        lsvs = {}
        for lsv_id in m.lsv_ids():
            lsv = m.lsv(lsv_id)
            lsvs[lsv_id] = {
                'name': meta['group_names'][0],
                'junctions': lsv.junctions.tolist(),
                'group_means': dict(lsv.group_means),
                'group_bins': dict(lsv.group_bins)
            }

    return jsonify({'lsvs': lsvs, 'exon_numbers': exon_numbers})
