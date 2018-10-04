class SpliceGraphTools {
    constructor(sgs) {
        console.log('top');
        this.sgs = sgs;
        this._init();
    }

    _init() {
        const gene = this.sgs.gene;

        // menu drop downs
        document.querySelector('.splice-graph-tools.tools-menu-btn').onclick = (event) => {
            event.preventDefault();
            d3.selectAll('.lsv-tools.tools-menu')
                .classed('hide-tools-menu', true);
            document.querySelector('.splice-graph-tools.tools-menu').classList.toggle('hide-tools-menu');
        };

        document.querySelector('.lsv-tools.tools-menu-btn').onclick = (event) => {
            event.preventDefault();
            d3.selectAll('.splice-graph-tools.tools-menu')
                .classed('hide-tools-menu', true);
            document.querySelector('.lsv-tools.tools-menu').classList.toggle('hide-tools-menu');
        };

        // header info
        document.querySelector('.gene-header .gene-name').textContent = `Gene name: ${gene.name}; ${ gene.chromosome }:${ gene.strand }:${ gene.start }-${ gene.end };`;
        document.querySelector('.gene-header .gene-id').textContent = `Gene ID: ${gene.id};`;
        document.querySelector('.ucsc-gene').setAttribute('href', `http://genome.ucsc.edu/cgi-bin/hgTracks?db=${gene.genome}&position=${ gene.chromosome }:${ gene.start }-${ gene.end }`)


        json_ajax('/metadata').then(meta => {

            // populate splice graph selector groups
            d3.select('.groups select')
                .selectAll('option')
                .data(meta.group_names)
                .enter()
                .append('option')
                .text(d => {
                    return d
                });


            // populate splice graph experiments when group is changed
            document.querySelector('#groups').onchange = (event) => {
                const group_name = meta.group_names[event.target.selectedIndex];

                const shown_exps = Array.from(document.querySelectorAll('.splice-graph'))
                    .filter(sg => sg.dataset.group === group_name)
                    .map(sg => sg.dataset.experiment);

                const exps = meta.experiment_names[event.target.selectedIndex]
                    .filter(e => !shown_exps.includes(e));

                const s = d3.select('.experiments select')
                    .selectAll('option')
                    .data(exps);

                s.text(d => d);

                s.enter()
                    .append('option')
                    .text(d => d);

                s.exit().remove();
            };

            // force change event to populate experiments on initial page load
            SpliceGraphTools._populate_sg_form();
        });


        // submit event for splice graph selector
        document.querySelector('.splice-graph-form').onsubmit = (event) => {
            event.preventDefault();
            const f = event.target;
            const g = f.querySelector('#groups').value;
            const e = f.querySelector('#experiments').value;
            if (g && e) {
                f.querySelector('button').disabled = true;
                this.sgs.create(g, e);
                SpliceGraphTools._populate_sg_form();
                f.querySelector('button').disabled = false;
            }
        };

        // toggle splice graph scale
        document.querySelector('.toggle-scale').onclick = () => {
            document.querySelector('.splice-graph-container').classList.toggle('default-view');
            this.sgs.update(250);
        };

        // zoom in on splice graph
        document.querySelector('.zoom-in').onclick = () => {
            let zoom = parseInt(document.querySelector('.splice-graph-container').dataset.zoom);
            zoom++;
            document.querySelector('.splice-graph-container').dataset.zoom = zoom;
            this.sgs.update(250);
        };

        // zoom out on splice graph
        document.querySelector('.zoom-out').onclick = () => {
            let zoom = parseInt(document.querySelector('.splice-graph-container').dataset.zoom);
            if (zoom !== 1) {
                zoom--;
                document.querySelector('.splice-graph-container').dataset.zoom = zoom;
                this.sgs.update(250);
            }
        };

        // reset zoom for splice graph
        document.querySelector('.zoom-reset').onclick = () => {
            document.querySelector('.splice-graph-container').dataset.zoom = 1;
            this.sgs.update(250);
        };

        // activate/deactivate junction reads filter
        document.querySelector('#junction-reads-filter').onchange = (event) => {
            document.querySelectorAll('#reads-greater-than, #reads-less-than').forEach(el => el.disabled = !el.disabled);
            if (event.target.checked) {
                document.querySelectorAll('#reads-greater-than, #reads-less-than').forEach(e => e.dispatchEvent(new Event('input')));
            } else {
                d3.selectAll('.junction-grp')
                    .classed('reads-filter', false)
            }
        };

        // adjust greater than and less than fields in junction filter
        const junctions_filter = () => {
            const gt = parseInt(document.querySelector('#reads-greater-than').value);
            const lt = parseInt(document.querySelector('#reads-less-than').value);
            d3.selectAll('.junction-grp')
                .classed('reads-filter', (d, i, a) => {
                    let r = parseInt(a[i].querySelector('.junction-reads').textContent);
                    r = isNaN(r) ? 0 : r;
                    return (!isNaN(gt) && !isNaN(lt) && r <= gt || r >= lt) || (!isNaN(gt) && r <= gt) || (!isNaN(lt) && r >= lt);
                })
        };
        document.querySelector('#reads-greater-than').oninput = junctions_filter;
        document.querySelector('#reads-less-than').oninput = junctions_filter;


    }

    copy_select_ucsc(el) {
        this.sgs.gene.then(gene => {
            const d = d3.select(el).datum();
            const textArea = document.createElement("textarea");
            textArea.value = `${gene.chromosome}:${d.start}-${d.end}`;
            document.body.appendChild(textArea);
            textArea.focus();
            textArea.select();
            document.execCommand('copy');
            document.body.removeChild(textArea);
        })
    }

    static _populate_sg_form() {
        document.querySelector('#groups').dispatchEvent(new Event('change'));
    }

}

