class SpliceGraphTools {
    constructor(sgs, gene, highlight_lsvs) {
        this.grp_names = gene.group_names;
        this.exp_names = gene.experiment_names;
        this.sgs = sgs;
        this.highlight_lsvs = highlight_lsvs;
        this._init();
        this.download_svg()
    }

    _init() {
        const gene = this.sgs.gene;


        function toggle_splicegraph_opts(){
            d3.selectAll('.lsv-tools.tools-menu')
                .classed('hide-tools-menu', true);
            document.querySelector('.splice-graph-tools.tools-menu').classList.toggle('hide-tools-menu');
            $('#splice-graph-menu-btn').toggleClass('pure-menu-active')
        }

        // menu drop downs
        document.querySelector('.splice-graph-tools.tools-menu-btn').onclick = (event) => {
            event.preventDefault();
            event.stopPropagation();
            toggle_splicegraph_opts();
        };

        // handling click outside of splice graph options hides the options
        $("body").click(function(event){
           if($('#splice-graph-menu-btn').hasClass('pure-menu-active')){
                toggle_splicegraph_opts();
           }
        });

        // (and click inside the options does not hide it)
        $(".splice-graph-tools").click(function () {
           event.stopPropagation();
        });

        document.querySelectorAll('.lsv-tools.tools-menu-btn')
            .forEach(l => l.onclick = (event) => {
                event.preventDefault();
                d3.selectAll('.splice-graph-tools.tools-menu')
                    .classed('hide-tools-menu', true);
                document.querySelector('.lsv-tools.tools-menu').classList.toggle('hide-tools-menu');
            });


        // populate splice graph selector groups
        d3.select('.groups select')
            .selectAll('option')
            .data(this.grp_names)
            .enter()
            .append('option')
            .text(d => {
                return d
            });


        // populate splice graph experiments when group is changed
        document.querySelector('#groups').onchange = (event) => {
            const group_name = this.grp_names[event.target.selectedIndex];

            const shown_exps = Array.from(document.querySelectorAll('.splice-graph'))
                .filter(sg => sg.dataset.group === group_name)
                .map(sg => sg.dataset.experiment);

            const exps = this.exp_names[event.target.selectedIndex]
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
                junctions_filter();
                this.highlight_lsvs();
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
                junctions_filter()
            } else {
                this.sgs.junctions_filter()
            }
        };

        // adjust greater than and less than fields in junction filter
        const junctions_filter = () => {
            if (document.querySelector('#junction-reads-filter').checked) {
                const gt = document.querySelector('#reads-greater-than').value;
                const lt = document.querySelector('#reads-less-than').value;
                this.sgs.junctions_filter(gt, lt)
            }
        };

        document.querySelector('#reads-greater-than').oninput = junctions_filter;
        document.querySelector('#reads-less-than').oninput = junctions_filter;


    }


    static _populate_sg_form() {
        document.querySelector('#groups').dispatchEvent(new Event('change'));
    }


    download_svg() {
        window.addEventListener('click', e => {
            if (e.target.classList.contains('splice-graph-download')) {
                const sg = e.target.closest('.splice-graph');
                const exp = sg.dataset.experiment;
                const grp = sg.dataset.group;
                const svg = sg.querySelector('svg').outerHTML;

                const element = document.createElement('a');
                element.setAttribute('href', 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svg));
                element.setAttribute('download', `${grp}_${exp}_sg.svg`);

                element.style.display = 'none';
                document.body.appendChild(element);

                element.click();

                document.body.removeChild(element);
            }
        })
    }

}

