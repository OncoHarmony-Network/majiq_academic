from flask import Flask, render_template, jsonify, request, url_for

from voila.api import Matrix

app = Flask(__name__)

voila = '/Users/cjgreen/Development/het_test/majiq_psi/Adr.psi.voila'
splice_graph = '/Users/cjgreen/Development/het_test/majiq_build/splicegraph.sql'


@app.route('/')
def hello_world():
    return render_template('psi_index.html')


@app.route('/<gene_id>/')
def gene(gene_id):
    return render_template('psi_summary.html')


@app.route('/results', methods=('POST',))
def results():
    with Matrix(voila) as p:
        records = []
        form = request.form
        draw = form['draw']
        start = int(form['start'])
        length = int(form['length'])
        column_sort = int(form['order[0][column]'])
        sort_direction = form['order[0][dir]']
        lsv_ids = list(p.lsv_ids())
        search_value = form['search[value]']

        for lsv_id in lsv_ids:
            psi = p.psi(lsv_id)
            records.append([psi.gene_id, psi.lsv_id, psi.lsv_type, 'psi per junc', 'links'])

        records.sort(key=lambda d: d[column_sort], reverse=sort_direction == 'desc')

        records = list(filter(lambda x: any(search_value in y for y in x), records))

        filtered_len = len(records)

        records = records[start:start + length]

        for r in records:
            r[0] = '<a href="' + url_for('gene', gene_id=r[0]) + '">' + r[0] + '</a>'

        return jsonify({
            'data': records,
            'draw': draw,
            'recordsTotal': len(lsv_ids),
            'recordsFiltered': filtered_len
        })


if __name__ == "__main__":
    app.run(host='127.0.0.1')
