class Table {
    constructor(table, opts) {

        if (typeof table === 'string')
            this.table = document.querySelector(table);
        else
            this.table = table;


        if (opts !== undefined) {
            this.db_gene = opts.db_gene;
            this.db_lsv = opts.db_lsv;
            this.show_data = opts.show_data;

            this.options = {limit: 5, include_docs: true};
            this.page_number = 1;
            this.total_pages = 1;
            this.curr_page_rows = 10;
            this.urlParams = new URLSearchParams(document.location.search);
            this.lsv = new Lsv(db_lsv, db_gene);
            this.violin = new Violin(db_lsv);
            this.next();
            this.previous();
            this.page_rows();
            this.update();
        }


        this.table.onclick = (event) => {
            const col = event.target;
            if (col.tagName === 'TH') {
                const col_idx = Array.from(this.table.querySelectorAll('th')).indexOf(col);
                col.asc = !col.asc;
                this.sort(col_idx, col.asc);
            }
        }
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

    static load_db(db_gene, gene_id) {
        const load = new Promise(resolve => {
            const scriptTag = document.createElement('script');
            scriptTag.src = `${gene_id}.js`;
            scriptTag.onload = () => resolve();
            scriptTag.onreadystatechange = () => resolve();
            document.body.appendChild(scriptTag);
        });

        const check = new Promise(resolve => {
            const re_check = () => {
                db_gene.get(gene_id).then(() => {
                    resolve()
                }).catch(err => {
                    if (err.reason === 'missing')
                        re_check();
                });
            };
            re_check();
        });

        return Promise.all([load, check])
    };

    filter_data() {
        return new Promise(resolve => {

            const filter = {};
            const gene_id = this.urlParams.get('gene_id');
            const lsv_filter = document.querySelector('.lsv-filters');

            lsv_filter.querySelector('#prime-5').checked ? filter.A5SS = true : null;
            lsv_filter.querySelector('#prime-3').checked ? filter.A3SS = true : null;
            lsv_filter.querySelector('#exon-skipping').checked ? filter.exon_skipping = true : null;
            lsv_filter.querySelector('#target').checked ? filter.target = true : null;
            lsv_filter.querySelector('#source').checked ? filter.target = false : null;
            lsv_filter.querySelector('#binary').checked ? filter.binary = true : null;
            lsv_filter.querySelector('#complex').checked ? filter.binary = false : null;


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
                .sort((a, b) => a._id.localeCompare(b._id));

            this.total_pages = Math.ceil(lsv_ids.length / this.curr_page_rows);

            if (this.total_pages < this.page_number || this.page_number === 0)
                this.page_number = this.total_pages;

            const slice_start = (this.page_number - 1) * this.curr_page_rows;
            const slice_end = this.page_number * this.curr_page_rows;
            const lsv_ids_page = lsv_ids.slice(slice_start, slice_end);

            resolve(lsv_ids_page)
        })
    }


    retrieve_data(results) {
        const gene_ids = Array.from(new Set(results.map(l => l.gene_id)));

        return Promise.all(gene_ids.map(gene_id => Table.load_db(this.db_gene, gene_id))).then(() => {
            return this.db_lsv.allDocs({
                keys: results.map(r => r._id),
                include_docs: true
            }).then(lsvs => lsvs.rows.map(r => r.doc))
        })
    }

    update() {
        console.time('update');
        this.filter_data()
            .then(data => this.retrieve_data(data))
            .then(data => this.show_data(data, this, this.body))
            .then(() => this.update_toolbar())
            .then(() => console.timeEnd('update'));
    }

    comparator(col_idx, asc) {
        return (a, b) => {
            const a1 = Table.col_value(a, col_idx);
            const b1 = Table.col_value(b, col_idx);
            return asc ? a1.toString().localeCompare(b1) : b1.toString().localeCompare(a1);
        };
    }

    sort(col_idx, asc) {
        Array.from(this.table.querySelectorAll('tbody tr'))
            .sort(this.comparator(col_idx, asc))
            .forEach(tr => this.table.querySelector('tbody').appendChild(tr));
    }

    update_toolbar() {
        const prev_btn = document.querySelector('.previous');
        const next_btn = document.querySelector('.next');
        const page = document.querySelector('.current-page');
        next_btn.disabled = this.page_number === this.total_pages;
        prev_btn.disabled = this.page_number === 1;
        page.innerHTML = `Page ${this.page_number} of ${this.total_pages}`
    }

    next() {
        const next_btn = document.querySelector('.next');
        next_btn.onclick = (event) => {
            event.preventDefault();
            this.page_number++;
            this.update();
        }
    }

    previous() {
        const prev_btn = document.querySelector('.previous');
        prev_btn.onclick = (event) => {
            event.preventDefault();
            this.page_number--;
            this.update();
        }
    }

    page_rows() {
        document.querySelector('.rows').onchange = (event) => {
            const s = event.target;
            this.curr_page_rows = parseInt(s.options[s.selectedIndex].value);
            this.update();
        }
    }

    highlight_form(cell, sgs) {
        const div = cell.append('div')
            .attr('class', 'highlight-form');

        div.append('div')
            .append('label')
            .text('Highlight')
            .append('input')
            .attr('class', 'highlight')
            .attr('type', 'checkbox')
            .on('change', () => {
                const hl = Array.from(document.querySelectorAll('input.highlight:checked'))
                    .map(el => el.closest('tr').dataset.lsvId);
                const w = Array.from(document.querySelectorAll('input.psi-weighted:checked'))
                    .map(el => el.closest('tr').dataset.lsvId);
                sgs.highlight(hl, w);
            });

        div.append('div')
            .append('label')

            .text('PSI Weighted')
            .append('input')
            .attr('class', 'psi-weighted')
            .attr('type', 'checkbox')
            .on('change', (d, i, a) => {
                if (a[i].checked) {
                    const hl = a[i].closest('.highlight-form').querySelector('.highlight');
                    hl.checked = true;
                    hl.dispatchEvent(new Event('change'))
                }
            })
    }


    ucsc_link(a_tag, d) {
        db_gene.get(d.gene_id).then(gene => {
            // The order of the junctions array is super important, so we'll copy the array before sorting.
            const juncs = Array.from(d.junctions);
            juncs.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
            const first_junc = juncs[0];
            const last_junc = juncs[juncs.length - 1];
            const first_exon = gene.exons.find(e => first_junc[0] >= e.start & first_junc[0] <= e.end);
            const last_exon = gene.exons.find(e => last_junc[1] >= e.start & last_junc[1] <= e.end);
            a_tag.setAttribute('href', `http://genome.ucsc.edu/cgi-bin/hgTracks?db=${gene.genome}&position=${gene.chromosome}:${first_exon.start}-${last_exon.end}`);
        })
    };

    majiq_spel_button(btn, d) {
        btn.onclick = () => {
            const textArea = document.createElement("textarea");
            textArea.value = `test value`;
            document.body.appendChild(textArea);
            textArea.focus();
            textArea.select();
            document.execCommand('copy');
            document.body.removeChild(textArea);
        }
    };

    psi_summary_compact(canvas, group_name) {
        if (group_name)
            canvas.dataset.group = group_name;
        canvas.style.display = 'block';
        this.lsv.draw_lsv_compact_stack_bars(canvas, 1);
        canvas.onclick = () => {
            const violin = canvas.parentNode.querySelector('.psi-violin-plot');
            violin.style.display = 'block';
            canvas.style.display = 'none';
        }
    }

    psi_summary_violin(svg, group_name) {
        if (group_name)
            svg.dataset.group = group_name;
        svg.style.display = 'None';
        this.violin.psi(svg);
        svg.onclick = () => {
            const comp = svg.parentNode.querySelector('.lsv-single-compact-percentiles');
            svg.style.display = 'none';
            comp.style.display = 'block';
        }
    }

    dpsi_summary_excl_incl(div, datum) {
        div.style.display = 'block';
        this.lsv.draw_delta_lsv_compact_svg(div, datum);
        div.onclick = () => {
            const svg = div.parentNode.querySelector('.deltapsi-violin-plot');
            div.style.display = 'none';
            svg.style.display = 'block'
        }
    }

    dpsi_summary_violin(svg) {
        this.violin.deltapsi(svg);
        svg.style.display = 'none';
        svg.onclick = () => {
            const div = svg.parentNode.querySelector('.excl-incl-rect');
            svg.style.display = 'none';
            div.style.display = 'block'
        }
    }
}