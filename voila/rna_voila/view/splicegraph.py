import os
from bisect import bisect

from flask import url_for, jsonify, request, session, redirect

from rna_voila.api.view_splice_graph import ViewSpliceGraph
from rna_voila.api.splice_graph_lr import SpliceGraphLR
from rna_voila.view import views
from rna_voila.config import ViewConfig

app, bp = views.get_bp(__name__)

@bp.before_request
def init_session():
    if not 'omit_simplified' in session:
        session['omit_simplified'] = True

@bp.route('/')
def index():
    with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:
        first_gene_id = sorted(sg.gene_ids)[0]
        return redirect(url_for('main.gene', gene_id=first_gene_id))

@bp.route('/toggle-simplified', methods=('POST',))
def toggle_simplified():
    session['omit_simplified'] = not session['omit_simplified']
    return jsonify({'ok':1})

@bp.route('/gene/<gene_id>/')
@bp.route('/gene')
def gene(gene_id=None):
    with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:
        if gene_id not in sg.gene_ids:
            return redirect(url_for('main.gene_not_found', gene_id=gene_id))
    return views.gene_view('sg_summary.html', gene_id, ViewSpliceGraph)


@bp.route('/gene-not-found/<gene_id>/')
def gene_not_found(gene_id):
    with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:
        if gene_id in sg.gene_ids:
            return redirect(url_for('main.gene', gene_id=gene_id))

    return '<h1>' + gene_id + '</h1>' + '<h3>Gene ID was not found in splice graph.</h3>'


@bp.route('/nav/<gene_id>', methods=('POST',))
def nav(gene_id):
    with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:
        gene_ids = sorted(sg.gene_ids)
        idx = bisect(gene_ids, gene_id)

        try:
            return jsonify({
                'next': url_for('main.gene', gene_id=gene_ids[idx % len(gene_ids)]),
                'prev': url_for('main.gene', gene_id=gene_ids[(idx % len(gene_ids)) - 2])
            })
        except:
            return jsonify({'next': '#', 'prev': '#'})


@bp.route('/splice-graph/<gene_id>', methods=('POST', 'GET'))
def splice_graph(gene_id):
    with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:
        exp_names = [sg.experiment_names]
        gd = sg.gene_experiment(gene_id, exp_names)
        gd['experiment_names'] = exp_names
        gd['group_names'] = ['splice graph']
        return jsonify(gd)


@bp.route('/splice-graph/lr/<gene_id>', methods=('POST', 'GET'))
def splice_graph_lr(gene_id):
    with SpliceGraphLR(ViewConfig().long_read_file) as sgl:
        with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:
            annot_exons = [(x['annotated_start'], x['annotated_end'],) for x in sg.exons(gene_id) if x['annotated']]
            gd = sgl.gene(gene_id, annot_exons)

            #print(gd)
            return jsonify(gd)

@bp.route('/splice-graph/combined/<gene_id>', methods=('POST', 'GET'))
def splice_graph_combined(gene_id):

    with SpliceGraphLR(ViewConfig().long_read_file) as sgl:
        with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:

            sr = sg.gene_experiment(gene_id, [])
            gd = sgl.combined_gene(gene_id, sr)

            # if not sg.experiment_names:
            #     for exon in gd['exons']:
            #         exon['color'] =

            return jsonify(gd)

@bp.route('/psi-splice-graphs', methods=('POST',))
def psi_splice_graphs():
    with ViewSpliceGraph(omit_simplified=session.get('omit_simplified', False)) as sg:
        try:
            sg_init = session['psi_init_splice_graphs']
        except KeyError:
            try:
                sg_init = [['splice graph', sg.experiment_names[0]]]
            except IndexError:
                sg_init = [['splice graph', 'no experiment']]

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


app.register_blueprint(bp)