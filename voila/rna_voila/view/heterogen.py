from bisect import bisect
from operator import itemgetter
from statistics import median
import numpy as np
from flask import render_template, jsonify, url_for, request, session, Response

from rna_voila.api.view_matrix import ViewHeterogens, ViewHeterogen
from rna_voila.api.view_splice_graph import ViewSpliceGraph
from rna_voila.index import Index
from rna_voila.view import views
from rna_voila.view.datatables import DataTables
from rna_voila.view.forms import LsvFiltersForm, HeterogenFiltersForm
from rna_voila.config import ViewConfig
import json

app, bp = views.get_bp(__name__)

def get_group_name_voila_file_map():
    group_map = {}
    config = ViewConfig()
    for voila_file in config.voila_files:
        with ViewHeterogen(voila_file) as m:
            group_map[frozenset(m.group_names)] = [voila_file, m.num_lsv_ids, True]

    return group_map

@bp.before_request
def init_session():
    if not 'omit_simplified' in session:
        session['omit_simplified'] = True

@bp.route('/')
def index():
    form = LsvFiltersForm()
    with ViewHeterogens() as m:
        het_form = HeterogenFiltersForm(m.stat_names)
        group_names = m.group_names

    if not 'group_name_voila_file_map' in session:
        session['group_map'] = get_group_name_voila_file_map()
    return render_template('het_index.html', form=form, het_form=het_form,
                           group_names=group_names, frozenset=frozenset,
                           enable_het_comparison_chooser=ViewConfig().enable_het_comparison_chooser)

@bp.route('/reindex', methods=('POST',))
def reindex():

    enabled_voila_files = []
    if 'comparisons' in request.json:
        for comp_files in request.json['comparisons']:
            # format goes file1, file2, enable or disable
            key = frozenset((comp_files[0], comp_files[1]))
            if comp_files[2]:
                enabled_voila_files.append(session['group_map'][key][0])
            session['group_map'][key][2] = comp_files[2]

    Index(force_create=True, voila_files=enabled_voila_files)

    return jsonify({'ok':1})

@bp.route('/dismiss-warnings', methods=('POST',))
def dismiss_warnings():
    session['warnings'] = []
    return jsonify({'ok':1})

@bp.route('/toggle-simplified', methods=('POST',))
def toggle_simplified():
    session['omit_simplified'] = not session['omit_simplified']
    return jsonify({'ok':1})

@bp.route('/reset-group-settings', methods=('POST',))
def reset_group_settings():
    del session['group_order_override']
    del session['group_display_name_override']
    del session['group_visibility']
    return jsonify({'ok':1})

@bp.route('/update-group-order', methods=('POST',))
def update_group_list():
    session['group_order_override'] = request.json
    return jsonify({'ok':1})

@bp.route('/update-group-display-names', methods=('POST',))
def update_group_display_names():
    session['group_display_name_override'] = request.json
    return jsonify({'ok':1})

@bp.route('/update-group-visibility', methods=('POST',))
def update_group_visibility():
    session['group_visibility'] = request.json
    return jsonify({'ok':1})

@bp.route('/gene/<gene_id>/')
def gene(gene_id):

    with ViewHeterogens() as m, ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:


        ucsc = {}
        exon_numbers = {}

        for het in m.lsvs(gene_id):
            lsv_junctions = het.junctions
            lsv_exons = sg.lsv_exons(gene_id, lsv_junctions)
            start, end = views.lsv_boundries(lsv_exons)
            gene = sg.gene(gene_id)
            ucsc[het.lsv_id] = url_for('main.generate_ucsc_link', lsv_id=het.lsv_id)
            exon_num = views.find_exon_number(sg.exons(gene_id), het.reference_exon, gene['strand'])
            if type(exon_num) is int:
                # accounting for 'unk' exon numbers
                exon_numbers[het.lsv_id] = exon_num
            else:
                exon_numbers[het.lsv_id] = 0

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

        return views.gene_view('het_summary.html', gene_id, ViewHeterogens,
                               lsv_data=lsv_data,
                               group_names=m.group_names,
                               ucsc=ucsc,
                               stat_names=m.stat_names,
                               analysis_type='heterogen')


@bp.route('/lsv-data', methods=('POST',))
@bp.route('/lsv-data/<lsv_id>', methods=('POST',))
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


@bp.route('/index-table', methods=('POST',))
def index_table():
    with ViewHeterogens() as p, ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:
        dt = DataTables(Index.heterogen(), ('gene_name', 'lsv_id', 'lsv_type', 'links', 'dpsi_threshold', 'stat_threshold'), sort=False, slice=False)
        stat_i = dt.heterogen_filters()

        dt.sort()
        dt.slice()

        for idx, index_row, records in dt.callback():
            _values = itemgetter('lsv_id', 'gene_id', 'gene_name', 'dpsi_threshold', 'stat_threshold')(index_row)
            values = [v.decode('utf-8') for v in _values[:-2]]
            values.append(round(float(_values[-2].decode('utf-8')), 3))  # dpsi_threshold
            values.append(round(json.loads(_values[-1].decode('utf-8'))[stat_i], 3))  # stat_threshold
            lsv_id, gene_id, gene_name, max_dpsi, min_stats = values

            het = p.lsv(lsv_id)
            lsv_junctions = het.junctions
            lsv_exons = sg.lsv_exons(gene_id, lsv_junctions)
            start, end = views.lsv_boundries(lsv_exons)
            gene = sg.gene(gene_id)
            ucsc = views.ucsc_href(sg.genome, gene['chromosome'], start, end)
            records[idx] = [
                (url_for('main.gene', gene_id=gene_id),
                 gene_name),
                lsv_id,
                het.lsv_type,
                {
                    'ucsc': ucsc,
                    'group_names': p.group_names
                },
                max_dpsi,
                min_stats
            ]

        return jsonify(dict(dt))


@bp.route('/nav/<gene_id>', methods=('POST',))
def nav(gene_id):
    with ViewHeterogens() as h:
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


@bp.route('/splice-graph/<gene_id>', methods=('POST',))
def splice_graph(gene_id):
    with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg, ViewHeterogens() as v:
        exp_names = v.splice_graph_experiment_names
        if ViewConfig().disable_reads:
            gd = sg.gene_experiment(gene_id, [])
            exp_names = [['splice graph']]
        else:
            gd = sg.gene_experiment(gene_id, exp_names)
        gd['group_names'] = v.group_names
        gd['experiment_names'] = exp_names
        return jsonify(gd)


@bp.route('/psi-splice-graphs', methods=('POST',))
def psi_splice_graphs():

    sg_init = session.get('psi_init_splice_graphs', [])
    json_data = request.get_json()

    if json_data:
        if 'add' in json_data:
            if all(s != json_data['add'] for s in sg_init):
                sg_init.append(json_data['add'])

        if 'remove' in json_data:
            sg_init = filter(lambda s: s != json_data['remove'], sg_init)
            sg_init = list(sg_init)


    if sg_init:
        session['psi_init_splice_graphs'] = sg_init

    return jsonify(sg_init)


@bp.route('/lsv-highlight', methods=('POST',))
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

                    output_lsv_data = {
                        'junctions': junctions,
                        'intron_retention': intron_retention,
                        'reference_exon': het.reference_exon,
                        'weighted': weighted
                    }

                    if weighted:
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

                                    group_means[grp_name][comb_name].append(median((v for v in exp if v != -1)))

                        output_lsv_data['group_means'] = group_means

                    lsvs.append(output_lsv_data)

        return jsonify(lsvs)

def rename_groups(group_names):

    if not session.get('group_display_name_override', None):
        return group_names

    group_names_new = [session['group_display_name_override'][n] for n in group_names]
    return group_names_new


@bp.route('/summary-table', methods=('POST',))
def summary_table():
    lsv_id, stat_name = itemgetter('lsv_id', 'stat_name')(request.form)

    with ViewHeterogens(group_order_override=session.get('group_order_override', None)) as v:
        exp_names = v.experiment_names
        grp_names = v.group_names

        if 'group_visibility' in session:
            # this is reversed because we are removing these indexes from lists later, and that only works
            # consistently if we do it backwards
            hidden_idx_unsorted = [grp_names.index(name) for name in session['group_visibility'] if session['group_visibility'][name] is False]
            hidden_idx = sorted([int(x) for x in hidden_idx_unsorted], reverse=True)
        else:
            hidden_idx = []

        het = v.lsv(lsv_id)
        juncs = het.junctions
        mu_psis = het.mu_psi
        mean_psis = het.mean_psi
        median_psis = het.median_psi.T

        table_data = []

        skipped_idx = 0
        for idx, (junc, mean_psi, mu_psi, median_psi) in enumerate(zip(juncs, mean_psis, mu_psis, median_psis)):

            junc = map(str, junc)
            junc = '-'.join(junc)
            heatmap = het.junction_heat_map(stat_name, idx)

            table_data.append({
                'junc': junc,
                'junc_idx': idx - skipped_idx,
                'mean_psi': mean_psi,
                'mu_psi': mu_psi,
                'median_psi': list(median_psi),
                'heatmap': heatmap,
            })

        dt = DataTables(table_data, ('junc', '', ''))

        o_grp_names = grp_names.copy()
        o_exp_names = exp_names.copy()
        for _idx in hidden_idx:
            del o_grp_names[_idx]
            del o_exp_names[_idx]

        for idx, row_data, records in dt.callback():
            junc, junc_idx, mean_psi = itemgetter('junc', 'junc_idx', 'mean_psi')(row_data)
            mu_psi, heatmap, median_psi = itemgetter('mu_psi', 'heatmap', 'median_psi')(row_data)
            for _idx in hidden_idx:
                heatmap = np.delete(heatmap, _idx, axis=0)
                heatmap = np.delete(heatmap, _idx, axis=1).tolist()
                del mu_psi[_idx]
                del mean_psi[_idx]

            records[idx] = [
                junc,
                {
                    'group_names': rename_groups(o_grp_names),
                    'experiment_names': o_exp_names,
                    'junction_idx': junc_idx,
                    'mean_psi': mean_psi,
                    'mu_psi': mu_psi,
                    'median_psi': median_psi
                },
                {
                    'heatmap': heatmap,
                    'group_names': rename_groups(o_grp_names),
                    'stat_name': stat_name
                }
            ]

        return jsonify(dict(dt))


@bp.route('/download-lsvs', methods=('POST',))
def download_lsvs():
    dt = DataTables(Index.heterogen(), ('gene_name', 'lsv_id'), slice=False)

    data = (d['lsv_id'].decode('utf-8') for d in dict(dt)['data'])
    data = '\n'.join(data)

    return Response(data, mimetype='text/plain')


@bp.route('/download-genes', methods=('POST',))
def download_genes():
    dt = DataTables(Index.heterogen(), ('gene_name', 'lsv_id'), slice=False)

    data = set(d['gene_id'].decode('utf-8') for d in dict(dt)['data'])
    data = '\n'.join(data)

    return Response(data, mimetype='text/plain')


@bp.route('/copy-lsv', methods=('POST',))
@bp.route('/copy-lsv/<lsv_id>', methods=('POST',))
def copy_lsv(lsv_id):
    return views.copy_lsv(lsv_id, ViewHeterogens)

@bp.route('/generate_ucsc_link', methods=('GET',))
def generate_ucsc_link():
    return views._generate_ucsc_link(request.args, ViewHeterogens)

app.register_blueprint(bp)