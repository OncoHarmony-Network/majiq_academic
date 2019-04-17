import os
from bisect import bisect
from operator import itemgetter
from statistics import median

import numpy as np
from flask import Flask, render_template, jsonify, url_for, request, session, Response

from voila.api.view_matrix import ViewDeltaPsi, ViewHeterogens
from voila.api.view_splice_graph import ViewSpliceGraph
from voila.index import Index
from voila.view import views
from voila.view.datatables import DataTables
from voila.view.forms import LsvFiltersForm

app = Flask(__name__)
app.secret_key = os.urandom(16)


@app.route('/')
def index():
    form = LsvFiltersForm()
    return render_template('het_index.html', form=form)

@app.route('/toggle-simplified', methods=('POST',))
def toggle_simplified():
    if not 'omit_simplified' in session:
        session['omit_simplified'] = True
    else:
        session['omit_simplified'] = not session['omit_simplified']
    return jsonify({'ok':1})

@app.route('/gene/<gene_id>/')
def gene(gene_id):
    with ViewHeterogens() as m, ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:


        ucsc = {}
        exon_numbers = {}

        for het in m.lsvs(gene_id):
            lsv_junctions = het.junctions
            lsv_exons = sg.lsv_exons(gene_id, lsv_junctions)
            start, end = views.lsv_boundries(lsv_exons)
            gene = sg.gene(gene_id)
            ucsc[het.lsv_id] = views.ucsc_href(sg.genome, gene['chromosome'], start, end)
            exon_numbers[het.lsv_id] = views.find_exon_number(sg.exons(gene_id), het.reference_exon, gene['strand'])


        lsv_data = []
        lsv_is_source = {}
        for lsv_id in m.lsv_ids(gene_ids=[gene_id]):
            lsv = m.lsv(lsv_id)

            lsv_data.append( [lsv_id, lsv.lsv_type] )
            lsv_is_source[lsv_id] = 1 if lsv.source else 0


        # this is the default sort, so modify the list, and add the indexes
        lsv_data.sort(key=lambda x: (exon_numbers[x[0]], lsv_is_source[x[0]]))

        type_length_idx = [i[0] for i in sorted(enumerate(lsv_data), key=lambda x: len(x[1][1].split('|')))]

        for i, lsv in enumerate(lsv_data):
            # appending exon number
            lsv.append(exon_numbers[lsv[0]])
            # appending default sort index
            lsv.append(i)
            # appending other sort indexes
            lsv.append(type_length_idx[i])

        return views.gene_view('het_summary.html', gene_id, ViewDeltaPsi,
                               lsv_data=lsv_data,
                               group_names=m.group_names,
                               ucsc=ucsc,
                               stat_names=m.stat_names)


@app.route('/lsv-data', methods=('POST',))
@app.route('/lsv-data/<lsv_id>', methods=('POST',))
def lsv_data(lsv_id):
    with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg, ViewHeterogens() as m:
        het = m.lsv(lsv_id)
        gene_id = het.gene_id
        gene = sg.gene(gene_id)
        strand = gene['strand']
        exons = sg.exons(gene_id)
        ref_exon = het.reference_exon
        exon_number = views.find_exon_number(exons, ref_exon, strand)

        return jsonify({
            'lsv': {
                'junctions': het.junctions.tolist(),
            },
            'exon_number': exon_number
        })


@app.route('/index-table', methods=('POST',))
def index_table():
    with ViewHeterogens() as p, ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:
        dt = DataTables(Index.heterogen(), ('gene_name', 'lsv_id'))

        for idx, index_row, records in dt.callback():
            values = itemgetter('lsv_id', 'gene_id', 'gene_name')(index_row)
            values = [v.decode('utf-8') for v in values]
            lsv_id, gene_id, gene_name = values
            het = p.lsv(lsv_id)
            lsv_junctions = het.junctions
            lsv_exons = sg.lsv_exons(gene_id, lsv_junctions)
            start, end = views.lsv_boundries(lsv_exons)
            gene = sg.gene(gene_id)
            ucsc = views.ucsc_href(sg.genome, gene['chromosome'], start, end)
            records[idx] = [
                (url_for('gene', gene_id=gene_id),
                 gene_name),
                lsv_id,
                het.lsv_type,
                {
                    'ucsc': ucsc,
                    'group_names': p.group_names
                }
            ]

        return jsonify(dict(dt))


@app.route('/nav/<gene_id>', methods=('POST',))
def nav(gene_id):
    with ViewDeltaPsi() as h:
        gene_ids = list(sorted(h.gene_ids))
        if len(gene_ids) == 1:
            return jsonify({
                'next': url_for('gene', gene_id=gene_ids[0]),
                'prev': url_for('gene', gene_id=gene_ids[0])
            })
        idx = bisect(gene_ids, gene_id)

        return jsonify({
            'next': url_for('gene', gene_id=gene_ids[idx % len(gene_ids)]),
            'prev': url_for('gene', gene_id=gene_ids[(idx % len(gene_ids)) - 2])
        })


@app.route('/splice-graph/<gene_id>', methods=('POST',))
def splice_graph(gene_id):
    with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg, ViewHeterogens() as v:
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

    with ViewHeterogens() as m:

        lsvs = []
        highlight_dict = session.get('highlight', {})

        for lsv_id, highlight, weighted in json_data:
            highlight_dict[lsv_id] = [highlight, weighted]

        session['highlight'] = highlight_dict
        splice_graphs = session.get('psi_init_splice_graphs', {})

        if splice_graphs:
            for lsv_id, (highlight, weighted) in highlight_dict.items():
                if highlight:
                    group_means = {}

                    het = m.lsv(lsv_id)
                    junctions = het.junctions.tolist()

                    if het.lsv_type[-1] == 'i':
                        intron_retention = junctions[-1]
                        junctions = junctions[:-1]
                    else:
                        intron_retention = []

                    means = np.array(het.mu_psi).transpose((1, 2, 0))

                    for gn, ens, xs in zip(m.group_names, m.experiment_names, means):

                        if any(sg[0] == gn for sg in splice_graphs):

                            for en, x in zip(ens, xs):

                                if any(sg[1] == en for sg in splice_graphs):

                                    x = x.tolist()

                                    try:
                                        group_means[gn][en] = x
                                    except KeyError:
                                        group_means[gn] = {en: x}

                    for junc in het.mu_psi:
                        for grp_name, exp in zip(m.group_names, junc):
                            if any(sg[0] == grp_name and sg[1].endswith(' Combined') for sg in splice_graphs):
                                comb_name = grp_name + ' Combined'

                                if grp_name not in group_means:
                                    group_means[grp_name] = {}

                                if comb_name not in group_means[grp_name]:
                                    group_means[grp_name][comb_name] = []

                                group_means[grp_name][comb_name].append(median(exp))

                    lsvs.append({
                        'junctions': junctions,
                        'intron_retention': intron_retention,
                        'reference_exon': het.reference_exon,
                        'weighted': weighted,
                        'group_means': group_means
                    })

        return jsonify(lsvs)


@app.route('/summary-table', methods=('POST',))
def summary_table():
    lsv_id, stat_name = itemgetter('lsv_id', 'stat_name')(request.form)
    if 'hidden_idx' in request.form:
        # this is reversed because we are removing these indexes from lists later, and that only works
        # consistently if we do it backwards
        hidden_idx = sorted([int(x) for x in request.form['hidden_idx'].split(',')], reverse=True)
    else:
        hidden_idx = []

    with ViewHeterogens() as v:
        exp_names = v.experiment_names
        grp_names = v.group_names
        for _idx in hidden_idx:
            del grp_names[_idx]
            del exp_names[_idx]

        het = v.lsv(lsv_id)
        juncs = het.junctions
        mu_psis = het.mu_psi
        mean_psis = het.mean_psi

        table_data = []

        skipped_idx = 0
        for idx, (junc, mean_psi, mu_psi) in enumerate(zip(juncs, mean_psis, mu_psis)):
            if idx in hidden_idx:
                skipped_idx += 1
                continue
            junc = map(str, junc)
            junc = '-'.join(junc)
            heatmap = het.junction_heat_map(stat_name, idx)

            table_data.append({
                'junc': junc,
                'junc_idx': idx - skipped_idx,
                'mean_psi': mean_psi,
                'mu_psi': mu_psi,
                'heatmap': heatmap,
            })

        dt = DataTables(table_data, ('junc', '', ''))

        for idx, row_data, records in dt.callback():
            junc, junc_idx, mean_psi = itemgetter('junc', 'junc_idx', 'mean_psi')(row_data)
            mu_psi, heatmap = itemgetter('mu_psi', 'heatmap')(row_data)
            for _idx in hidden_idx:
                heatmap = np.delete(heatmap, _idx, axis=0)
                heatmap = np.delete(heatmap, _idx, axis=1).tolist()
                del mu_psi[_idx]
                del mean_psi[_idx]

            records[idx] = [
                junc,
                {
                    'group_names': grp_names,
                    'experiment_names': exp_names,
                    'junction_idx': junc_idx,
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


@app.route('/download-lsvs', methods=('POST',))
def download_lsvs():
    dt = DataTables(Index.heterogen(), ('gene_name', 'lsv_id'), slice=False)

    data = (d['lsv_id'].decode('utf-8') for d in dict(dt)['data'])
    data = '\n'.join(data)

    return Response(data, mimetype='text/plain')


@app.route('/download-genes', methods=('POST',))
def download_genes():
    dt = DataTables(Index.heterogen(), ('gene_name', 'lsv_id'), slice=False)

    data = set(d['gene_id'].decode('utf-8') for d in dict(dt)['data'])
    data = '\n'.join(data)

    return Response(data, mimetype='text/plain')


@app.route('/copy-lsv', methods=('POST',))
@app.route('/copy-lsv/<lsv_id>', methods=('POST',))
def copy_lsv(lsv_id):
    return views.copy_lsv(lsv_id, ViewHeterogens)
