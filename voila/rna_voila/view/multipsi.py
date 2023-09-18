import os
from bisect import bisect
from operator import itemgetter

from flask import render_template, url_for, jsonify, request, session, Response, redirect, abort


from rna_voila.api.view_matrix import ViewPsi, ViewPsis
from rna_voila.api.view_splice_graph import ViewSpliceGraph
from rna_voila.index import Index
from rna_voila.view import views
from rna_voila.view.datatables import DataTables
from rna_voila.view.forms import LsvFiltersForm
from rna_voila.config import ViewConfig
from rna_voila.exceptions import LsvIdNotFoundInVoilaFile, LsvIdNotFoundInAnyVoilaFile, GeneIdNotFoundInVoilaFile


app, bp = views.get_bp(__name__)

@bp.before_request
def init_session():
    if not 'omit_simplified' in session:
        session['omit_simplified'] = True

@bp.route('/')
def index():
    form = LsvFiltersForm()
    return render_template('multipsi_index.html', form=form)

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

    with ViewPsis() as m, ViewSpliceGraph() as sg:
        ucsc = {}
        exon_numbers = {}
        filter_exon_numbers = {}
        for lsv in m.lsvs(gene_id):
            lsv_junctions = lsv.junctions
            lsv_exons = sg.lsv_exons(gene_id, lsv_junctions)
            start, end = views.lsv_boundries(lsv_exons)
            gene = sg.gene(gene_id)
            ucsc[lsv.lsv_id] = url_for('main.generate_ucsc_link', lsv_id=lsv.lsv_id)
            exon_num = views.find_exon_number(sg.exons(gene_id), lsv.reference_exon, gene['strand'])
            if type(exon_num) is int:
                # accounting for 'unk' exon numbers
                exon_numbers[lsv.lsv_id] = exon_num
                if not exon_num in filter_exon_numbers:
                    filter_exon_numbers[exon_num] = [lsv.lsv_id]
                else:
                    filter_exon_numbers[exon_num].append(lsv.lsv_id)
            else:
                exon_numbers[lsv.lsv_id] = 0

        lsv_data = []
        lsv_is_source = {}
        for lsv_id in m.lsv_ids(gene_ids=[gene_id]):
            lsv = m.lsv(lsv_id)

            lsv_data.append([lsv_id, lsv.lsv_type])
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



        return views.gene_view('multipsi_summary.html', gene_id, ViewPsis,
                               lsv_data=lsv_data,
                               group_names=m.group_names,
                               ucsc=ucsc,
                               filter_exon_numbers=filter_exon_numbers
                               )


@bp.route('/index-table', methods=('POST',))
def index_table():

    with ViewPsis() as v, ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:
        grp_name = v.group_names[0]

        dt = DataTables(Index.psi(), ('gene_name', 'lsv_id'))

        for idx, index_row, records in dt.callback():
            values = itemgetter('gene_name', 'gene_id', 'lsv_id')(index_row)
            values = [x.decode('utf-8') for x in values]
            gene_name, gene_id, lsv_id = values

            psi = v.lsv(lsv_id)
            gene = sg.gene(gene_id)
            lsv_exons = sg.lsv_exons(gene_id, psi.junctions)
            start, end = views.lsv_boundries(lsv_exons)
            ucsc = views.ucsc_href(sg.genome, gene['chromosome'], start, end)

            records[idx] = [
                {'href': url_for('main.gene', gene_id=gene_id), 'gene_name': gene_name},
                lsv_id,
                psi.lsv_type
            ]
            if len(ViewConfig().voila_files) == 1:
                records[idx].append(grp_name)
            records[idx].append(ucsc)

        return jsonify(dict(dt))


@bp.route('/nav/<gene_id>', methods=('POST',))
def nav(gene_id):
    with ViewPsis() as h:
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
    with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg, ViewPsis() as v:
        exp_names = v.splice_graph_experiment_names
        gd = sg.gene_experiment(gene_id, exp_names)
        gd['experiment_names'] = exp_names
        gd['group_names'] = v.group_names
        return jsonify(gd)


@bp.route('/summary-table/<gene_id>', methods=('POST',))
def summary_table(gene_id):

    with ViewPsis() as v, ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:
        grp_name = v.group_names[0]
        index_data = Index.psi(gene_id)

        dt = DataTables(index_data, ('highlight', 'lsv_id'), sort=False, slice=False)

        dt.add_sort('highlight', DataTables.Sort.highlight)
        dt.add_sort('lsv_id', DataTables.Sort.lsv_id)

        dt.sort()
        dt.slice()

        for idx, record, records in dt.callback():
            lsv_id = record['lsv_id'].decode('utf-8')
            psi = v.lsv(lsv_id)
            lsv_type = psi.lsv_type

            gene = sg.gene(gene_id)
            lsv_exons = sg.lsv_exons(gene_id, psi.junctions)
            start, end = views.lsv_boundries(lsv_exons)
            ucsc = views.ucsc_href(sg.genome, gene['chromosome'], start, end)

            try:
                highlight = session['highlight'][lsv_id]
            except KeyError:
                highlight = [False, False]

            records[idx] = [
                highlight,
                lsv_id,
                lsv_type,
                grp_name,
                ucsc
            ]

        return jsonify(dict(dt))


@bp.route('/psi-splice-graphs', methods=('POST',))
def psi_splice_graphs():
    with ViewPsis() as v:
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


@bp.route('/lsv-data', methods=('POST',))
@bp.route('/lsv-data/<lsv_id>', methods=('POST',))
def lsv_data(lsv_id):

    with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg, ViewPsis() as m:
        psi = m.lsv(lsv_id)
        ref_exon = psi.reference_exon
        gene_id = psi.gene_id
        gene = sg.gene(gene_id)
        strand = gene['strand']
        exons = sg.exons(gene_id)
        exon_number = views.find_exon_number(exons, ref_exon, strand)

        return jsonify({
            'lsv': {
                'name': m.group_names[0],
                'junctions': psi.junctions.tolist(),
                'group_means': dict(psi.group_means),
                'group_bins': dict(psi.group_bins)
            },
            'exon_number': exon_number
        })

@bp.route('/violin-data', methods=('POST',))
@bp.route('/violin-data/<lsv_id>', methods=('POST',))
def violin_data(lsv_id):
    config = ViewConfig()
    if 'hidden_idx' in request.form:
        hidden_idx = sorted([int(x) for x in request.form['hidden_idx'].split(',')], reverse=True)
    else:
        hidden_idx = []

    """
    Expected workflow:
    For each lsv, we first get the union of all possible junctions from all voila files
    For each found junction, we represent a table row
    Then, For each voila file, we look for that LSV and check if the junction is available in it.
    If so, we  add group bins / means to that table row
    """
    with ViewPsis() as v:
        exp_names = v.experiment_names
        grp_names = v.group_names
        files = config.voila_files[:]
        for idx in hidden_idx:
            grp_names.pop(idx)
            exp_names.pop(idx)
            files.pop(idx)

        all = v.lsv(lsv_id)

        table_data = []

        for i, _junc in enumerate(all.junctions.tolist()):


            """
            In one table row, "group_bins" is a 2d array. Outer array index refers to the test group index
            (voila file index). 
            
            For each DataTable box (one cell in the table), violin plots are made for all horizontal elements 
            at once. So we need to provide the array of group_bins in terms of junction rather than groups
            For example, the first element of group_bins should represent all data from the first junction, 
            and the value of this first element should be an array of data for each group (from the first junction)
            
            'group_means': [ <junc1> , <junc2> ]
            'group_means': [ [ <group1>, <group2> ] , <junc2> ]
            """
            table_data.append([
                _junc,
                {
                    'junction_idx': i,
                    'junction_name': _junc,
                    "group_names": grp_names,
                    "experiment_names": exp_names,
                    'group_means': [[] for _ in range(len(all.junctions.tolist()))],
                    'group_bins': [[] for _ in range(len(all.junctions.tolist()))],
                }
            ])

            for j, grp in enumerate(grp_names):

                with ViewPsi(files[j]) as m:

                    try:
                        psi = m.lsv(lsv_id)
                        means = dict(psi.group_means)[grp][i]
                        bins = dict(psi.group_bins)[grp][i]
                        juncs = psi.junctions.tolist()
                    except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
                        means = []
                        bins = []
                        juncs = []

                    if _junc in juncs:

                        table_data[-1][1]['group_means'][i].append(means)
                        table_data[-1][1]['group_bins'][i].append(bins)

                    else:
                        table_data[-1][1]['group_means'][i].append([])
                        table_data[-1][1]['group_bins'][i].append([])

        dt = DataTables(table_data, (), sort=False)


        return jsonify(dict(dt))


@bp.route('/lsv-highlight', methods=('POST',))
def lsv_highlight():
    json_data = request.get_json()

    with ViewPsis() as m:

        lsvs = []
        highlight_dict = session.get('highlight', {})

        for lsv_id, highlight, weighted in json_data:
            highlight_dict[lsv_id] = [highlight, weighted]

        session['highlight'] = highlight_dict
        splice_graphs = session.get('psi_init_splice_graphs', {})

        if splice_graphs:
            for lsv_id, (highlight, weighted) in highlight_dict.items():
                if highlight:
                    psi = m.lsv(lsv_id)
                    junctions = psi.junctions.tolist()

                    if psi.lsv_type[-1] == 'i':
                        intron_retention = junctions[-1]
                        junctions = junctions[:-1]
                    else:
                        intron_retention = []

                    means = dict(psi.all_group_means)
                    group_means = {}

                    for sg in splice_graphs:
                        grp_name, exp_name = sg
                        if grp_name not in group_means:
                            group_means[grp_name] = {}
                        if grp_name in means:
                            group_means[grp_name][exp_name] = means[grp_name]

                    lsvs.append({
                        'junctions': junctions,
                        'intron_retention': intron_retention,
                        'reference_exon': list(psi.reference_exon),
                        'weighted': weighted,
                        'group_means': group_means

                    })

        return jsonify(lsvs)


@bp.route('/download-lsvs', methods=('POST',))
def download_lsvs():
    dt = DataTables(Index.psi(), ('gene_name', 'lsv_id'), slice=False)

    data = (d['lsv_id'].decode('utf-8') for d in dict(dt)['data'])
    data = '\n'.join(data)

    return Response(data, mimetype='text/plain')


@bp.route('/download-genes', methods=('POST',))
def download_genes():
    dt = DataTables(Index.psi(), ('gene_name', 'lsv_id'), slice=False)

    data = set(d['gene_id'].decode('utf-8') for d in dict(dt)['data'])
    data = '\n'.join(data)

    return Response(data, mimetype='text/plain')


@bp.route('/copy-lsv', methods=('POST',))
@bp.route('/copy-lsv/<lsv_id>', methods=('POST',))
def copy_lsv(lsv_id):
    return views.copy_lsv(lsv_id, ViewPsi, voila_file=ViewConfig().voila_files[0])


@bp.route('/generate_ucsc_link', methods=('GET',))
def generate_ucsc_link():
    return views._generate_ucsc_link(request.args, ViewPsis)

@bp.route('/transcripts/<gene_id>', methods=('POST', 'GET'))
def transcripts(gene_id):
    with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:
        return jsonify(sg.gene_transcript_exons(gene_id))


app.register_blueprint(bp)