class LsvTools {
    constructor(sgs) {
        this.sgs = sgs;
        this.db_lsv = sgs.db_lsv;
        this.db_gene = sgs.db_gene;
        this.lsv = new Lsv(this.db_lsv, this.db_gene);
        this.heatmap = new HeatMap(this.db_lsv)
    }

    load_lsvs() {
        return db_lsv.createIndex({
            index: {fields: ['gene_id']}
        }).then(() => {
            return db_lsv.find({
                selector: {
                    gene_id: urlParams.get('gene_id')
                },
                fields: ['_id', 'junctions', 'A3SS', 'A5SS', 'exon_skipping', 'is_target']
            })
        }).then(results => {
            return LsvTools.lsv_filter(results.docs)
        });
    }


    static lsv_filter(lsvs) {
        const prime5 = document.querySelector('#prime-5').checked;
        const prime3 = document.querySelector('#prime-3').checked;
        const exon_skipping = document.querySelector('#exon-skipping').checked;
        const target = document.querySelector('#target').checked;
        const source = document.querySelector('#source').checked;
        return lsvs
            .filter(e => !prime3 ? true : e['A3SS'])
            .filter(e => !prime5 ? true : e['A5SS'])
            .filter(e => !exon_skipping ? true : e.exon_skipping)
            .filter(e => !target ? true : e.target)
            .filter(e => !source ? true : !e.target);
    }

    het_enter_lsv(lsv_divs) {

        const l = lsv_divs
            .enter()
            .append('div')
            .attr('class', 'lsv')
            .attr('data-lsv-id', d => d._id);

        const header = l
            .append('div')
            .attr('class', 'lsv-header');

        const fieldset = header
            .append('form')
            .on('change', () => {
                const hl = Array.from(document.querySelectorAll('input#highlight:checked'))
                    .map(el => el.closest('.lsv').dataset.lsvId);
                const w = Array.from(document.querySelectorAll('input#psi-weighted:checked'))
                    .map(el => el.closest('.lsv').dataset.lsvId);
                this.sgs.highlight(hl, w);
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

        const copy_lsv = fieldset.append('div')
            .attr('class', 'pure-control-group');

        copy_lsv.append('button')
            .attr('id', 'copy_lsv')
            .attr('class', 'pure-button')
            .text('Copy LSV');

        header
            .append('div')
            .attr('class', 'lsv-cartoon')
            .each((d, i, a) => {
                this.lsv.cartoon(a[i], d._id);
            });


        header
            .append('div')
            .text(d => {
                return d._id
            });


        const table = l
            .append('table')
            .classed('pure-table', true)
            .attr('data-lsv-id', d => {
                return d._id
            }).each((d, i, a) => new Table(a[i]));

        table.append('thead')
            .append('tr')
            .selectAll('th')
            .data(['Junctions', 'Violin', 'Heat Map'])
            .enter()
            .append('th')
            .text(d => {
                return d
            });

        const row = table
            .append('tbody')
            .selectAll('tr')
            .data(d => {
                return d.junctions
            })
            .enter()
            .append('tr')
            .attr('data-junction-index', (d, i) => i);

        row
            .append('td')
            .text(d => {
                return d
            });

        row
            .append('td')
            .append('svg')
            .attr('class', 'het-violin-plot')
            .attr('data-type', 'swarm')
            .each((d, i, a) => new Violin(db_lsv, a[i]));

        row
            .append('td')
            .append('svg')
            .attr('class', 'heat-map')
            .attr('data-stat-name', 'wilcoxon')
            .each((d, i, a) => this.heatmap.plot(a[i]));
    };

    het_show_lsvs(lsvs) {
        lsvs.sort((a, b) => {
            return a.junctions.length - b.junctions.length || a._id.localeCompare(b._id)
        });

        const lsv_divs = d3.select('.lsv-container')
            .selectAll('.lsv')
            .data(lsvs);

        this.het_enter_lsv(lsv_divs);

    };

}

