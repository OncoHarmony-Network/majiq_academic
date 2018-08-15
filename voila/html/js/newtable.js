class NewTable {
    constructor(table_selector, db_gene, db_lsv, opts) {
        this.db_gene = db_gene;
        this.db_lsv = db_lsv;
        this.table_selector = table_selector;
        this.show_data = opts.show_data;
        this.options = {limit: 5, include_docs: true};
        this.page_number = 0;
        this.urlParams = new URLSearchParams(document.location.search);

        this.next();
        this.previous();
        this.update();

        this.table.onclick = (event) => {
            const col = event.target;
            if (col.tagName === 'TH') {
                const col_idx = Array.from(this.table.querySelectorAll('th')).indexOf(col);
                col.asc = !col.asc;
                this.sort(col_idx, col.asc);
            }
        }
    }

    get table() {
        return document.querySelector(this.table_selector)
    }

    get body() {
        return this.table.querySelector('tbody');
    }

    get rows() {
        return this.body.querySelectorAll('tr')
    }

    static col_value(tr, col_idx) {
        const cell = tr.children[col_idx];
        return 'sortValue' in cell.dataset ? cell.dataset.sortValue : cell.textContent
    }

    load_db(gene_id) {
        return new Promise(resolve => {
            const scriptTag = document.createElement('script');
            scriptTag.src = `${gene_id}.js`;
            scriptTag.onload = () => resolve();
            scriptTag.onreadystatechange = () => resolve();
            document.body.appendChild(scriptTag);
        })
    };

    retrieve_data() {
        return new Promise(resolve => {
            const filter = {};
            const lsv_filter = document.querySelector('.lsv-filters');
            lsv_filter.querySelector('#prime-5').checked ? filter.A5SS = true : null;
            lsv_filter.querySelector('#prime-3').checked ? filter.A3SS = true : null;
            lsv_filter.querySelector('#exon-skipping').checked ? filter.exon_skipping = true : null;
            lsv_filter.querySelector('#target').checked ? filter.target = true : null;
            lsv_filter.querySelector('#source').checked ? filter.target = false : null;
            lsv_filter.querySelector('#binary').checked ? filter.binary = true : null;
            lsv_filter.querySelector('#complex').checked ? filter.binary = false : null;
            const gene_id = this.urlParams.get('gene_id');

            const lsv_ids = lsvs_arr
                .filter(l => {
                    if (gene_id)
                        return l.gene_id === gene_id;
                    else
                        return true
                })
                .filter(l => {
                    const ks = Object.keys(filter);
                    if (ks.length)
                        return ks.every(k => l[k] === filter[k]);
                    else
                        return true
                })
                .sort((a, b) => a._id.localeCompare(b._id))
                .slice(this.page_number * 5, (this.page_number + 1) * 5);
            resolve(lsv_ids)
        })
    }


    curate_data(results) {
        const gene_ids = Array.from(new Set(results.map(l => l.gene_id)));

        return Promise.all(gene_ids.map(gene_id => this.load_db(gene_id))).then(() => {
            return this.db_lsv.allDocs({
                keys: results.map(r => r._id),
                include_docs: true
            })
                .then(lsvs => {
                    if (lsvs.rows.some(r => r.error))
                        return this.curate_data(results);
                    else
                        return lsvs.rows.map(r => r.doc)
                })
        })
    }

    update() {
        console.time('update');
        this.retrieve_data()
            .then(data => this.curate_data(data))
            .then(data => this.show_data(data, this, this.body))
            .then(() => console.timeEnd('update'));
    }

    comparator(col_idx, asc) {
        return (a, b) => {
            const a1 = NewTable.col_value(a, col_idx);
            const b1 = NewTable.col_value(b, col_idx);
            return asc ? a1.toString().localeCompare(b1) : b1.toString().localeCompare(a1);
        };
    }

    sort(col_idx, asc) {
        Array.from(this.table.querySelectorAll('tbody tr'))
            .sort(this.comparator(col_idx, asc))
            .forEach(tr => this.table.querySelector('tbody').appendChild(tr));
    }


    next() {
        document.querySelector('.next').onclick = () => {
            this.page_number++;
            this.update();
        }
    }

    previous(table) {
        document.querySelector('.previous').onclick = () => {
            this.page_number--;
            this.update();
        }
    }

    highlight_form(cell, sgs) {
        const fieldset = cell
            .append('form')
            .on('change', () => {
                const hl = Array.from(document.querySelectorAll('input#highlight:checked'))
                    .map(el => el.closest('tr').dataset.lsvId);
                const w = Array.from(document.querySelectorAll('input#psi-weighted:checked'))
                    .map(el => el.closest('tr').dataset.lsvId);
                sgs.highlight(hl, w);
            })
            .attr('class', 'lsv-form pure-form pure-form-aligned')
            .append('fieldset');


        const highlight = fieldset.append('div')
            .attr('class', 'pure-control-group');

        highlight.append('label')
            .attr('for', 'highlight')
            .text('Highlight');

        highlight.append('input')
            .attr('id', 'highlight')
            .attr('type', 'checkbox');

        const psi_weighted = fieldset.append('div')
            .attr('class', 'pure-control-group');

        psi_weighted.append('label')
            .attr('for', 'psi-weighted')
            .text('PSI Weighted');

        psi_weighted.append('input')
            .attr('id', 'psi-weighted')
            .attr('type', 'checkbox')
            .on('change', (d, i, a) => {
                if (a[i].checked)
                    a[i].closest('fieldset').querySelector('input#highlight').checked = true;
            });
    }

    ucsc_link(els) {
        return els
            .each((d, i, a) => {
                db_gene.get(d.gene_id).then(gene => {
                    // The order of the junctions array is super important, so we'll copy the array before sorting.
                    const juncs = Array.from(d.junctions);
                    juncs.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
                    const first_junc = juncs[0];
                    const last_junc = juncs[juncs.length - 1];
                    const first_exon = gene.exons.find(e => first_junc[0] >= e.start & first_junc[0] <= e.end);
                    const last_exon = gene.exons.find(e => last_junc[1] >= e.start & last_junc[1] <= e.end);
                    a[i].setAttribute('href', `http://genome.ucsc.edu/cgi-bin/hgTracks?db=${gene.genome}&position=${gene.chromosome}:${first_exon.start}-${last_exon.end}`);
                })
            })
    };

}