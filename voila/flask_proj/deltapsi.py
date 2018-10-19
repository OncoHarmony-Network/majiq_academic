import os
from bisect import bisect

from flask import Flask, render_template, jsonify, url_for, request, redirect, session

from voila.api.view_matrix import ViewDeltaPsi
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.flask_proj.datatables import DataTables

app = Flask(__name__)
app.secret_key = os.urandom(16)


@app.route('/')
def index():
    return render_template('dpsi_index.html')


@app.route('/gene/<gene_id>/')
def gene(gene_id):
    # For this gene, remove any already selected highlight/weighted lsvs from session.
    highlight = session.get('highlight', {})
    lsv_ids = [h for h in highlight if h.startswith(gene_id)]
    for lsv_id in lsv_ids:
        del highlight[lsv_id]
    session['highlight'] = highlight

    return render_template('dpsi_summary.html', gene_id=gene_id)


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

    with ViewSpliceGraph() as sg, ViewDeltaPsi() as m:
        gene = sg.gene(gene_id)
        strand = gene.strand
        exons = sg.exons(gene)
        exon_number = find_exon_number(exons)

        dpsi = m.lsv(lsv_id)

        return jsonify({
            'lsv': {
                'excl_incl': list(dpsi.excl_incl),
                'junctions': dpsi.junctions.tolist(),
                'means': list(dpsi.means),
                'bins': dpsi.bins,
                'group_bins': dict(dpsi.group_bins),
                'group_means': dict(dpsi.group_means),
            },
            'exon_number': exon_number
        })


@app.route('/index-table', methods=('POST',))
def index_table():
    with ViewSpliceGraph() as sg, ViewDeltaPsi() as p:
        def create_record(lsv_ids):
            for lsv_id in lsv_ids:
                dpsi = p.lsv(lsv_id)
                gene_id = dpsi.gene_id
                gene_name = sg.gene(gene_id).name
                gene_name_col = {'sort': gene_name, 'display': [url_for('gene', gene_id=gene_id), gene_name]}

                excl_incl = dpsi.excl_incl
                excl_incl = max(abs(a - b) for a, b in excl_incl)

                yield [gene_name_col, dpsi.lsv_id, dpsi.lsv_type, excl_incl, 'links']

        records = list(create_record(p.lsv_ids()))

        dt = DataTables(records)

        return jsonify(dict(dt))


@app.route('/nav/<gene_id>', methods=('POST',))
def nav(gene_id):
    with ViewDeltaPsi() as h:
        gene_ids = list(sorted(h.gene_ids))
        idx = bisect(gene_ids, gene_id)

        return jsonify({
            'next': url_for('gene', gene_id=gene_ids[idx % len(gene_ids)]),
            'prev': url_for('gene', gene_id=gene_ids[(idx % len(gene_ids)) - 2])
        })


@app.route('/splice-graph/<gene_id>', methods=('POST', 'GET'))
def splice_graph(gene_id):
    if request.method == 'GET':
        return redirect(url_for('index'))

    with ViewSpliceGraph() as sg, ViewDeltaPsi() as v:
        meta = v.metadata
        g = sg.gene(gene_id)
        gd = sg.gene_experiment(g, meta['experiment_names'])
        return jsonify(gd)


@app.route('/psi-splice-graphs', methods=('POST',))
def psi_splice_graphs():
    with ViewDeltaPsi() as v:
        meta = v.metadata
        try:
            sg_init = session['psi_init_splice_graphs']
        except KeyError:
            sg_init = [[meta['group_names'][0], meta['experiment_names'][0][0]],
                       [meta['group_names'][1], meta['experiment_names'][1][0]]]

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


@app.route('/lsv-highlight', methods=('POST',))
def lsv_highlight():
    json_data = request.get_json()

    with ViewDeltaPsi() as m:

        lsvs = []
        highlight_dict = session.get('highlight', {})

        for lsv_id, highlight, weighted in json_data:
            highlight_dict[lsv_id] = [highlight, weighted]

        session['highlight'] = highlight_dict

        for lsv_id, (highlight, weighted) in highlight_dict.items():
            if highlight:
                lsv = m.lsv(lsv_id)
                lsvs.append({
                    'junctions': lsv.junctions.tolist(),
                    'reference_exon': list(lsv.reference_exon),
                    'target': lsv.target,
                    'weighted': weighted
                })

        return jsonify(lsvs)


@app.route('/summary-table/<gene_id>', methods=('POST',))
def summary_table(gene_id):
    with ViewDeltaPsi() as v:
        def create_records(lsv_ids):
            for lsv_id in lsv_ids:
                dpsi = v.lsv(lsv_id)
                excl_incl = dpsi.excl_incl
                excl_incl = max(abs(a - b) for a, b in excl_incl)
                ref_exon = dpsi.reference_exon
                lsv_id_col = {'sort': list(ref_exon), 'display': lsv_id}
                lsv_type = dpsi.lsv_type

                try:
                    highlight = session['highlight'][lsv_id]
                except KeyError:
                    highlight = [False, False]

                yield [highlight, lsv_id_col, lsv_type, grp_names[0], excl_incl, grp_names[1], 'links']

        grp_names = v.metadata['group_names']
        lsv_ids = v.lsv_ids(gene_ids=[gene_id])
        records = list(create_records(lsv_ids))

        dt = DataTables(records)
        dt = dict(dt)

        return jsonify(dt)


@app.route('/metadata', methods=('POST',))
def metadata():
    with ViewDeltaPsi() as v:
        return jsonify(v.metadata)
