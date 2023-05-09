import os
from bisect import bisect

from flask import render_template, jsonify, url_for, request, session, Response

from rna_voila.api.view_matrix import ViewDeltaPsi
from rna_voila.api.view_splice_graph import ViewSpliceGraph
from rna_voila.index import Index
from rna_voila.view import views
from rna_voila.view.datatables import DataTables
from rna_voila.view.forms import LsvFiltersForm, DeltaPsiFiltersForm
from rna_voila.config import ViewConfig

app, bp = views.get_bp(__name__)

@bp.before_request
def init_session():
    if not 'omit_simplified' in session:
        session['omit_simplified'] = True

@bp.route('/')
def index():
    form = LsvFiltersForm()
    dpsi_form = DeltaPsiFiltersForm()
    return render_template('dpsi_index.html', form=form, dpsi_form=dpsi_form)

@bp.route('/dismiss-warnings', methods=('POST',))
def dismiss_warnings():
    session['warnings'] = []
    return jsonify({'ok':1})

@bp.route('/toggle-simplified', methods=('POST',))
def toggle_simplified():
    session['omit_simplified'] = not session['omit_simplified']
    return jsonify({'ok':1})

@bp.route('/gene/<gene_id>/')
def gene(gene_id):
    with ViewDeltaPsi() as m, ViewSpliceGraph() as sg:
        filter_exon_numbers = {}
        for lsv in m.lsvs(gene_id):
            gene = sg.gene(gene_id)
            exon_num = views.find_exon_number(sg.exons(gene_id), lsv.reference_exon, gene['strand'])
            if type(exon_num) is int:
                # accounting for 'unk' exon numbers
                if not exon_num in filter_exon_numbers:
                    filter_exon_numbers[exon_num] = [lsv.lsv_id]
                else:
                    filter_exon_numbers[exon_num].append(lsv.lsv_id)

    return views.gene_view('dpsi_summary.html', gene_id, ViewDeltaPsi, filter_exon_numbers=filter_exon_numbers)


@bp.route('/lsv-data', methods=('POST',))
@bp.route('/lsv-data/<lsv_id>', methods=('POST',))
def lsv_data(lsv_id):
    with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg, ViewDeltaPsi() as m:
        dpsi = m.lsv(lsv_id)
        ref_exon = dpsi.reference_exon
        gene_id = dpsi.gene_id
        gene = sg.gene(gene_id)
        strand = gene['strand']
        exons = list(sg.exons(gene_id))
        exon_number = views.find_exon_number(exons, ref_exon, strand)

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


@bp.route('/index-table', methods=('POST',))
def index_table():
    with ViewDeltaPsi() as p, ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:
        dt = DataTables(Index.delta_psi(), ('gene_name', 'lsv_id', '', 'excl_incl'), slice=False)
        dt.delta_psi_filters()
        dt.slice()

        for idx, index_row, records in dt.callback():
            lsv_id = index_row['lsv_id'].decode('utf-8')
            excl_incl = index_row['excl_incl'].item()
            gene_id = index_row['gene_id'].decode('utf-8')
            gene_name = index_row['gene_name'].decode('utf-8')
            dpsi = p.lsv(lsv_id)

            gene = sg.gene(gene_id)
            lsv_junctions = dpsi.junctions
            lsv_exons = sg.lsv_exons(gene_id, lsv_junctions)

            start, end = views.lsv_boundries(lsv_exons)
            ucsc = url_for('main.generate_ucsc_link', lsv_id=lsv_id)

            records[idx] = [
                [url_for('main.gene', gene_id=gene_id), gene_name],
                lsv_id,
                dpsi.lsv_type,
                excl_incl,
                ucsc
            ]

        return jsonify(dict(dt))


@bp.route('/nav/<gene_id>', methods=('POST',))
def nav(gene_id):
    with ViewDeltaPsi() as h:
        gene_ids = list(sorted(h.gene_ids))
        if len(gene_ids) == 1:
            return jsonify({
                'next': url_for('main.gene', gene_id=gene_ids[0]),
                'prev': url_for('main.gene', gene_id=gene_ids[0])
            })
        idx = bisect(gene_ids, gene_id)

        return jsonify({
            'next': url_for('main.gene', gene_id=gene_ids[idx % len(gene_ids)]),
            'prev': url_for('main.gene', gene_id=gene_ids[(idx % len(gene_ids)) - 2])
        })


@bp.route('/splice-graph/<gene_id>', methods=('POST', 'GET'))
def splice_graph(gene_id):
    with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg, ViewDeltaPsi() as v:
        exp_names = v.splice_graph_experiment_names
        gd = sg.gene_experiment(gene_id, exp_names)
        gd['group_names'] = v.group_names
        gd['experiment_names'] = exp_names
        return jsonify(gd)


@bp.route('/psi-splice-graphs', methods=('POST',))
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


@bp.route('/lsv-highlight', methods=('POST',))
def lsv_highlight():
    json_data = request.get_json()

    with ViewDeltaPsi() as m:

        lsvs = []
        highlight_dict = session.get('highlight', {})

        for lsv_id, highlight, weighted in json_data:
            highlight_dict[lsv_id] = [highlight, weighted]

        session['highlight'] = highlight_dict

        splice_graphs = session.get('psi_init_splice_graphs', {})

        if splice_graphs:
            for lsv_id, (highlight, weighted) in highlight_dict.items():
                if highlight:
                    dpsi = m.lsv(lsv_id)

                    junctions = dpsi.junctions.tolist()

                    if dpsi.lsv_type[-1] == 'i':
                        intron_retention = junctions[-1]
                        junctions = junctions[:-1]
                    else:
                        intron_retention = []

                    means = dict(dpsi.group_means)
                    group_means = {}

                    for sg in splice_graphs:
                        grp_name, exp_name = sg
                        if grp_name not in group_means:
                            group_means[grp_name] = {}
                        group_means[grp_name][exp_name] = means[grp_name]

                    lsvs.append({
                        'junctions': junctions,
                        'intron_retention': intron_retention,
                        'reference_exon': list(dpsi.reference_exon),
                        'weighted': weighted,
                        'group_means': group_means

                    })

        return jsonify(lsvs)


@bp.route('/summary-table/<gene_id>', methods=('POST',))
def summary_table(gene_id):
    with ViewDeltaPsi() as v, ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:

        grp_names = v.group_names
        index_data = Index.delta_psi(gene_id)

        dt = DataTables(index_data, ('highlight', 'lsv_id', '', '', 'excl_incl'), sort=False, slice=False)

        dt.add_sort('highlight', DataTables.Sort.highlight)
        dt.add_sort('lsv_id', DataTables.Sort.lsv_id)

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

            gene = sg.gene(gene_id)
            lsv_junctions = dpsi.junctions
            lsv_exons = sg.lsv_exons(gene_id, lsv_junctions)
            start, end = views.lsv_boundries(lsv_exons)
            ucsc = url_for('main.generate_ucsc_link', lsv_id=lsv_id)

            records[idx] = [
                highlight,
                lsv_id,
                lsv_type,
                grp_names[0],
                excl_incl,
                grp_names[1],
                ucsc
            ]

        return jsonify(dict(dt))


@bp.route('/download-lsvs', methods=('POST',))
def download_lsvs():
    dt = DataTables(Index.delta_psi(), ('gene_name', 'lsv_id', '', 'excl_incl'), slice=False)
    dt.delta_psi_filters()

    data = (d['lsv_id'].decode('utf-8') for d in dict(dt)['data'])
    data = '\n'.join(data)

    return Response(data, mimetype='text/plain')


@bp.route('/download-genes', methods=('POST',))
def download_genes():
    dt = DataTables(Index.delta_psi(), ('gene_name', 'lsv_id', '', 'excl_incl'), slice=False)
    dt.delta_psi_filters()

    data = set(d['gene_id'].decode('utf-8') for d in dict(dt)['data'])
    data = '\n'.join(data)

    return Response(data, mimetype='text/plain')


@bp.route('/copy-lsv', methods=('POST',))
@bp.route('/copy-lsv/<lsv_id>', methods=('POST',))
def copy_lsv(lsv_id):
    return views.copy_lsv(lsv_id, ViewDeltaPsi)

@bp.route('/generate_ucsc_link', methods=('GET',))
def generate_ucsc_link():
    return views._generate_ucsc_link(request.args, ViewDeltaPsi)


app.register_blueprint(bp)