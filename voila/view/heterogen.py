import os
from bisect import bisect
from operator import itemgetter

from flask import Flask, render_template, jsonify, url_for, request, session

from voila.api.view_matrix import ViewDeltaPsi, ViewHeterogens
from voila.api.view_splice_graph import ViewSpliceGraph
from voila.index import Index
from voila.view.datatables import DataTables

app = Flask(__name__)
app.secret_key = os.urandom(16)


@app.route('/')
def index():
    return render_template('het_index.html')


@app.route('/gene/<gene_id>/')
def gene(gene_id):
    # For this gene, remove any already selected highlight/weighted lsvs from session.

    highlight = session.get('highlight', {})
    lsv_ids = [h for h in highlight if h.startswith(gene_id)]
    for lsv_id in lsv_ids:
        del highlight[lsv_id]
    session['highlight'] = highlight

    with ViewHeterogens() as m, ViewSpliceGraph() as sg:
        gene = sg.gene(gene_id)
        lsv_data = list((lsv_id, m.lsv(lsv_id).lsv_type) for lsv_id in m.lsv_ids(gene_ids=[gene_id]))
        lsv_data.sort(key=lambda x: len(x[1].split('|')))
        return render_template('het_summary.html', gene=gene, lsv_data=lsv_data)


@app.route('/lsv-data', methods=('POST',))
@app.route('/lsv-data/<lsv_id>', methods=('POST',))
def lsv_data(lsv_id):
    def find_exon_number(exons):
        ref_exon = list(map(int, lsv_id.split(':')[-1].split('-')))
        exons = filter(lambda e: -1 not in [e['start'], e['end']], exons)
        exons = list(exons)

        for idx, exon in enumerate(exons):
            if [exon['start'], exon['end']] == ref_exon:
                if strand == '-':
                    return len(exons) - idx
                else:
                    return idx + 1

    with ViewSpliceGraph() as sg, ViewHeterogens() as m:
        gene_id = ':'.join(lsv_id.split(':')[:-2])
        gene = sg.gene(gene_id)
        strand = gene['strand']
        exons = sg.exons(gene_id)

        # return empty string when reference exon is a half exon
        try:
            exon_number = find_exon_number(exons)
        except ValueError:
            exon_number = ''

        dpsi = m.lsv(lsv_id)

        return jsonify({
            'lsv': {
                # 'excl_incl': list(dpsi.excl_incl),
                'junctions': dpsi.junctions.tolist(),
                # 'means': list(dpsi.means),
                # 'bins': dpsi.bins,
                # 'group_bins': dict(dpsi.group_bins),
                # 'group_means': dict(dpsi.group_means),
            },
            'exon_number': exon_number
        })


@app.route('/index-table', methods=('POST',))
def index_table():
    with ViewHeterogens() as p:
        dt = DataTables(Index.heterogen(), ('gene_name', 'lsv_id'))

        for idx, index_row, records in dt.callback():
            values = itemgetter('lsv_id', 'gene_id', 'gene_name')(index_row)
            values = [v.decode('utf-8') for v in values]
            lsv_id, gene_id, gene_name = values
            het = p.lsv(lsv_id)

            records[idx] = [(url_for('gene', gene_id=gene_id), gene_name), lsv_id, het.lsv_type, '']

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


@app.route('/splice-graph/<gene_id>', methods=('POST',))
def splice_graph(gene_id):
    with ViewSpliceGraph() as sg, ViewHeterogens() as v:
        exp_names = v.splice_graph_experiment_names
        gd = sg.gene_experiment(gene_id, exp_names)
        gd['group_names'] = v.group_names
        gd['experiment_names'] = exp_names
        return jsonify(gd)


@app.route('/psi-splice-graphs', methods=('POST',))
def psi_splice_graphs():
    try:
        sg_init = session['psi_init_splice_graphs']
    except KeyError:
        sg_init = []

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
                het = m.lsv(lsv_id)
                junctions = het.junctions.tolist()

                if het.lsv_type[-1] == 'i':
                    intron_retention = junctions[-1]
                    junctions = junctions[:-1]
                else:
                    intron_retention = []

                lsvs.append({
                    'junctions': junctions,
                    'intron_retention': intron_retention,
                    'reference_exon': list(het.reference_exon),
                    'weighted': weighted,
                    'group_means': dict(het.group_means)

                })

        return jsonify(lsvs)


@app.route('/summary-table', methods=('POST',))
@app.route('/summary-table/<lsv_id>', methods=('POST',))
def summary_table(lsv_id):
    with ViewHeterogens() as v:
        exp_names = v.experiment_names
        grp_names = v.group_names
        stat_names = v.stat_names
        stat_name = stat_names[0]

        het = v.lsv(lsv_id)
        juncs = het.junctions
        mu_psis = het.mu_psi
        mean_psis = het.mean_psi

        table_data = []

        for idx, (junc, mean_psi, mu_psi) in enumerate(zip(juncs, mean_psis, mu_psis)):
            junc = map(str, junc)
            junc = '-'.join(junc)
            heatmap = het.junction_heat_map(stat_name, idx)

            table_data.append({
                'junc': junc,
                'mean_psi': mean_psi,
                'mu_psi': mu_psi,
                'heatmap': heatmap
            })

        dt = DataTables(table_data, ('junc', '', ''))

        for idx, row_data, records in dt.callback():
            junc, mean_psi, mu_psi, heatmap = itemgetter('junc', 'mean_psi', 'mu_psi', 'heatmap')(row_data)

            records[idx] = [
                junc,
                {
                    'group_names': grp_names,
                    'experiment_names': exp_names,
                    'junction_idx': idx,
                    'mean_psi': mean_psi,
                    'mu_psi': mu_psi,
                },
                {
                    'heatmap': heatmap,
                    'group_names': grp_names,
                    'stat_name': stat_name
                }
            ]

        return jsonify(dict(dt))
