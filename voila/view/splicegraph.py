import os
from bisect import bisect

from flask import render_template, url_for, jsonify, request, session, Flask, redirect

from voila.api.view_splice_graph_sqlite import ViewSpliceGraph

app = Flask(__name__)
app.secret_key = os.urandom(16)


@app.route('/')
def index():
    with ViewSpliceGraph() as sg:
        first_gene_id = sorted(sg.gene_ids())[0]
        return redirect(url_for('gene', gene_id=first_gene_id))


@app.route('/gene/<gene_id>/')
@app.route('/gene')
def gene(gene_id=None):
    # For this gene, remove any already selected highlight/weighted lsvs from session.
    highlight = session.get('highlight', {})
    lsv_ids = [h for h in highlight if h.startswith(gene_id)]
    for lsv_id in lsv_ids:
        del highlight[lsv_id]
    session['highlight'] = highlight

    return render_template('sg_summary.html', gene_id=gene_id)


@app.route('/nav/<gene_id>', methods=('POST',))
def nav(gene_id):
    with ViewSpliceGraph() as sg:
        gene_ids = list(sorted(sg.gene_ids()))
        idx = bisect(gene_ids, gene_id)

        return jsonify({
            'next': url_for('gene', gene_id=gene_ids[idx % len(gene_ids)]),
            'prev': url_for('gene', gene_id=gene_ids[(idx % len(gene_ids)) - 2])
        })


@app.route('/splice-graph/<gene_id>', methods=('POST', 'GET'))
def splice_graph(gene_id):
    with ViewSpliceGraph() as sg:
        g = sg.gene(gene_id)
        exp_names = [sg.experiment_names]
        gd = sg.gene_experiment(g, exp_names)
        gd['experiment_names'] = exp_names
        gd['group_names'] = ['splice graph']
        return jsonify(gd)


@app.route('/psi-splice-graphs', methods=('POST',))
def psi_splice_graphs():
    with ViewSpliceGraph() as sg:
        try:
            sg_init = session['psi_init_splice_graphs']
        except KeyError:
            sg_init = [['splice graph', sg.experiment_names[0]]]

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
