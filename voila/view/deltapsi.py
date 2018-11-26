import os
from bisect import bisect

from flask import Flask, render_template, jsonify, url_for, request, session, Response, redirect

from voila.api.view_matrix import ViewDeltaPsi
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.index import Index
from voila.view import views
from voila.view.datatables import DataTables
from voila.view.forms import LsvFiltersForm, DeltaPsiFiltersForm

app = Flask(__name__)
app.secret_key = os.urandom(16)


@app.route('/')
def index():
    form = LsvFiltersForm()
    dpsi_form = DeltaPsiFiltersForm()
    return render_template('dpsi_index.html', form=form, dpsi_form=dpsi_form)


@app.route('/gene/<gene_id>/')
def gene(gene_id):
    with ViewDeltaPsi() as m:
        if gene_id not in m.gene_ids:
            return redirect(url_for('index'))

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
    with ViewSpliceGraph() as sg, ViewDeltaPsi() as m:
        def find_exon_number(exons):
            exons = filter(lambda e: -1 not in [e.start, e.end], exons)
            exons = list(exons)

            for idx, exon in enumerate(exons):
                if [exon.start, exon.end] == ref_exon:
                    if strand == '-':
                        return len(exons) - idx
                    else:
                        return idx + 1

        dpsi = m.lsv(lsv_id)
        ref_exon = dpsi.reference_exon
        gene_id = dpsi.gene_id
        gene = sg.gene(gene_id)
        strand = gene.strand
        exons = list(sg.exons(gene))
        exon_number = find_exon_number(exons)

        excl_incl = list(dpsi.excl_incl)
        lsv_junctions = dpsi.junctions.tolist()
        means = list(dpsi.means)
        bins = dpsi.bins
        group_bins = dict(dpsi.group_bins)
        group_means = dict(dpsi.group_means)

        return jsonify({
            'lsv': {
                'excl_incl': excl_incl,
                'junctions': lsv_junctions,
                'means': means,
                'bins': bins,
                'group_bins': group_bins,
                'group_means': group_means,
            },
            'exon_number': exon_number
        })


@app.route('/index-table', methods=('POST',))
def index_table():
    with ViewDeltaPsi() as p:
        dt = DataTables(Index.delta_psi(), ('gene_name', 'lsv_id', '', 'excl_incl'), slice=False)
        dt.delta_psi_filters()
        dt.slice()

        for idx, index_row, records in dt.callback():
            lsv_id = index_row['lsv_id'].decode('utf-8')
            excl_incl = index_row['excl_incl'].item()
            gene_id = index_row['gene_id'].decode('utf-8')
            gene_name = index_row['gene_name'].decode('utf-8')

            records[idx] = [
                [url_for('gene', gene_id=gene_id), gene_name],
                lsv_id,
                p.lsv(lsv_id).lsv_type,
                excl_incl,
                ''
            ]

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
    with ViewSpliceGraph() as sg, ViewDeltaPsi() as v:
        g = sg.gene(gene_id)
        exp_names = v.splice_graph_experiment_names
        gd = sg.gene_experiment(g, exp_names)
        gd['group_names'] = v.group_names
        gd['experiment_names'] = exp_names
        return jsonify(gd)


@app.route('/psi-splice-graphs', methods=('POST',))
def psi_splice_graphs():
    with ViewDeltaPsi() as v:
        grp_names = v.group_names
        exp_names = v.splice_graph_experiment_names

        try:
            sg_init = session['psi_init_splice_graphs']
        except KeyError:
            sg_init = [[grp_names[0], exp_names[0][0]],
                       [grp_names[1], exp_names[1][0]]]

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
                dpsi = m.lsv(lsv_id)

                junctions = dpsi.junctions.tolist()

                if dpsi.lsv_type[-1] == 'i':
                    intron_retention = junctions[-1]
                    junctions = junctions[:-1]
                else:
                    intron_retention = []

                lsvs.append({
                    'junctions': junctions,
                    'intron_retention': intron_retention,
                    'reference_exon': list(dpsi.reference_exon),
                    'weighted': weighted,
                    'group_means': dict(dpsi.group_means)

                })

        return jsonify(lsvs)


@app.route('/summary-table/<gene_id>', methods=('POST',))
def summary_table(gene_id):
    with ViewDeltaPsi() as v:

        grp_names = v.group_names
        index_data = Index.delta_psi(gene_id)

        dt = DataTables(index_data, ('highlight', 'lsv_id', '', '', 'excl_incl'), sort=False, slice=False)

        dt.add_sort('highlight', DataTables.highlight)
        dt.add_sort('lsv_id', DataTables.lsv_id)

        dt.sort()
        dt.slice()

        for idx, record, records in dt.callback():
            lsv_id = record['lsv_id'].decode('utf-8')
            excl_incl = record['excl_incl'].item()
            dpsi = v.lsv(lsv_id)
            lsv_type = dpsi.lsv_type

            try:
                highlight = session['highlight'][lsv_id]
            except KeyError:
                highlight = [False, False]

            records[idx] = [
                highlight,
                lsv_id,
                lsv_type,
                grp_names[0],
                excl_incl,
                grp_names[1],
                ''
            ]

        return jsonify(dict(dt))


@app.route('/download-lsvs', methods=('POST',))
def download_lsvs():
    dt = DataTables(Index.delta_psi(), ('gene_name', 'lsv_id', '', 'excl_incl'), slice=False)
    dt.delta_psi_filters()

    data = (d['lsv_id'].decode('utf-8') for d in dict(dt)['data'])
    data = '\n'.join(data)

    return Response(data, mimetype='text/plain')


@app.route('/download-genes', methods=('POST',))
def download_genes():
    dt = DataTables(Index.delta_psi(), ('gene_name', 'lsv_id', '', 'excl_incl'), slice=False)
    dt.delta_psi_filters()

    data = set(d['gene_id'].decode('utf-8') for d in dict(dt)['data'])
    data = '\n'.join(data)

    return Response(data, mimetype='text/plain')


@app.route('/copy-lsv', methods=('POST',))
@app.route('/copy-lsv/<lsv_id>', methods=('POST',))
def copy_lsv(lsv_id):
    return views.copy_lsv(lsv_id, ViewDeltaPsi)
