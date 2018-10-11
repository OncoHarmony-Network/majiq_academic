import os
from bisect import bisect

from flask import Flask, render_template, jsonify, url_for, session, request, redirect
from waitress import serve

from voila import constants
from voila.api.view_matrix import ViewPsi
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.config import Config
from voila.flask_proj.datatables import DataTables
from voila.flask_proj.html import Html

app = Flask(__name__)
app.secret_key = os.urandom(16)


def run_service():
    serve(app, listen='127.0.0.1:55555')


def get_meta():
    with ViewPsi() as m:
        return m.view_metadata


def lsv_cartoon(lsv_type):
    html = Html()
    html.tag('canvas', classes=['lsv-cartoon'], width=200, height=80, data_lsv_type=lsv_type)
    return str(html)


def single_compact(group_name):
    html = Html()
    html.tag('canvas',
             classes=['lsv-single-compact-percentiles'],
             data_group=group_name)
    return str(html)


def psi_violin(group_name):
    html = Html()
    html.tag('svg', classes=['psi-violin-plot'], width=0, height=0, data_group=group_name)
    return str(html)


@app.route('/gene/<gene_id>/')
def gene(gene_id):
    config = Config()

    if config.voila_files == constants.ANALYSIS_PSI:
        return render_template('psi_summary.html', gene_id=gene_id)

    if config.voila_files == constants.ANALYSIS_DELTAPSI:
        return render_template('deltapsi_summary.html', gene_id=gene_id)


@app.route('/')
def index():
    config = Config()
    print(config.analysis_type)
    if config.voila_files == constants.ANALYSIS_PSI:
        return render_template('psi_index.html')

    if config.voila_files == constants.ANALYSIS_DELTAPSI:
        return render_template('deltapsi_index.html')


@app.route('/index-table', methods=('POST',))
def index_table():
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
            group_name = meta['group_names'][0]
            for r in rs:
                gene_id = ':'.join(r[1].split(':')[0:-2])

                html.tag('a', href=url_for('gene', gene_id=gene_id))
                html.text(r[0])
                r[0] = str(html)
                html.reset()

                r[2] = lsv_cartoon(r[2])
                r[3] = single_compact(group_name) + psi_violin(group_name)

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
            group_name = meta['group_names'][0]
            for r in rs:
                """
                <div class="highlight-form"><div><label>Highlight<input class="highlight" type="checkbox"></label></div><div><label>PSI Weighted<input class="psi-weighted" type="checkbox"></label></div></div>
                """

                labels = []
                for c in ['highlight', 'psi-weighted']:
                    label = Html().tag('div')
                    label = label.tag('label')
                    label.children(
                        Html().text(c.replace('-', ' ').title()),
                        Html().tag('input', classes=[c], type='checkbox')
                    )
                    labels.append(label)
                hl_form = Html().tag('div', classes=['highlight-form']).children(*labels)

                r[0] = str(hl_form)
                r[2] = lsv_cartoon(r[2])
                r[3] = single_compact(group_name) + psi_violin(group_name)

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
@app.route('/lsv-data/<lsv_id>', methods=('POST',))
def lsv_data(lsv_id):
    meta = get_meta()
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
                'name': meta['group_names'][0],
                'junctions': lsv.junctions.tolist(),
                'group_means': dict(lsv.group_means),
                'group_bins': dict(lsv.group_bins)
            },
            'exon_number': exon_number
        })


@app.route('/lsv-highlight', methods=('POST',))
def lsv_highlight():
    json_data = request.get_json()

    assert 'weighted' in json_data
    assert 'highlight' in json_data

    with ViewPsi() as m:
        lsvs = []

        for lsv_id in set(j for js in json_data.values() for j in js):
            lsv = m.lsv(lsv_id)
            lsvs.append({
                'junctions': lsv.junctions.tolist(),
                'reference_exon': list(lsv.reference_exon),
                'target': lsv.target,
                'weighted': lsv_id in json_data['weighted']
            })

        return jsonify(lsvs)
