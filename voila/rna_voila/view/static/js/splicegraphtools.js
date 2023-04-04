class SpliceGraphTools {
    constructor(sgs, gene, highlight_lsvs) {
        this.grp_names = gene.group_names;
        this.exp_names = gene.experiment_names;
        // for use with natural sorting
        this.collator = new Intl.Collator(undefined, {numeric: true, sensitivity: 'base'});
        this.sgs = sgs;
        this.highlight_lsvs = highlight_lsvs;
        this._init();
        this.download_svg()
        this.rescale_sg()
    }

    _init() {
        const gene = this.sgs.gene;

        // when clicked, will toggle all of the other active menus / buttons off and enable only ours
        function hide_others(){
            d3.selectAll('.tools-menu')
                .classed('hide-tools-menu', true);
            d3.selectAll('.tools-menu-btn')
                .classed('pure-menu-active', false);
        }

        function show_ours(){
            document.getElementById('splice-graph-tools-box').classList.remove('hide-tools-menu');
            $('#splice-graph-menu-btn').addClass('pure-menu-active')
        }

        // enable this menu when the button is clicked
        document.querySelector('#splice-graph-menu-btn').onclick = (event) => {
            event.preventDefault();
            event.stopPropagation();
            const already_active = $('#splice-graph-menu-btn').hasClass('pure-menu-active');
            hide_others();
            if(!already_active){
                show_ours();
            }
        };

        // handling click outside of our options hides them
        $("body").click(function(event){
            if($('#splice-graph-menu-btn').hasClass('pure-menu-active')){
                hide_others();
            }
        });

        // (and click inside the options does not hide it)
        $("#splice-graph-tools-box").click(function (event) {
            event.stopPropagation();
        });


        if($('.splice-graph-container .splice-graph').length === 0){
            setTimeout(function(){$('.splice-graph-form button').eq(0).trigger('click')}, 0);
        }



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
                .filter(e => !shown_exps.includes(e)).sort(this.collator.compare);

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
            const _type = 'short_read'
            if (g && e) {
                f.querySelector('button').disabled = true;
                this.sgs.create(g, e, _type);
                SpliceGraphTools._populate_sg_form();
                f.querySelector('button').disabled = false;
                junctions_filter();
                this.highlight_lsvs();
            }
        };

        document.querySelector('.splice-graph-form-lr').onsubmit = (event) => {
            event.preventDefault();
            const f = event.target;
            const g = 'LR';
            const e = 'LR';
            const _type = 'long_read'
            if (g && e) {
                f.querySelector('button').disabled = true;
                for(let transcript of this.sgs.gene_lr){
                    this.sgs.create(g, e, _type, transcript);
                }
                //SpliceGraphTools._populate_sg_form();
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

        // toggle locking psi filter disabling inputs
        const filter_lock_check = () => {
            const _readonly = document.querySelector('#lock-lr-psi-filter').checked;
            document.querySelector('#lr-psi-greater-than').readOnly = _readonly;
            document.querySelector('#lr-psi-less-than').readOnly = _readonly;
            if(_readonly){
                document.querySelector('#lr-psi-greater-than').value = document.querySelector('#psi-greater-than').value;
                document.querySelector('#lr-psi-less-than').value = document.querySelector('#psi-less-than').value;
            }
        }

        // adjust greater than and less than fields in junction filter
        const junctions_filter = () => {
            let gt, lt, gtl, ltl;
            let gtp, ltp, gtpl, ltpl;
            filter_lock_check();
            if (document.querySelector('#junction-reads-filter').checked) {
                gt = document.querySelector('#reads-greater-than').value;
                lt = document.querySelector('#reads-less-than').value;
            }
            if (document.querySelector('#junction-lr-reads-filter').checked) {
                gtl = document.querySelector('#lr-reads-greater-than').value;
                ltl = document.querySelector('#lr-reads-less-than').value;
            }
            if (document.querySelector('#junction-psi-filter').checked) {
                gtp = document.querySelector('#psi-greater-than').value;
                ltp = document.querySelector('#psi-less-than').value;
            }
            if (document.querySelector('#junction-lr-psi-filter').checked) {
                gtpl = document.querySelector('#lr-psi-greater-than').value;
                ltpl = document.querySelector('#lr-psi-less-than').value;
            }
            let presence = document.querySelector('input[name="sg-toggles"]:checked').value;
            if(presence){
                presence = presence.split(',');
            }
            console.log(presence)
            this.sgs.junctions_filter(gt, lt, gtl, ltl, gtp, ltp, gtpl, ltpl, presence);
        };


        document.querySelector('#lock-lr-psi-filter').onchange = filter_lock_check;
        $('input[name="sg-toggles"]').click(junctions_filter);
        document.querySelector('#junction-reads-filter').onchange = junctions_filter;
        document.querySelector('#junction-lr-reads-filter').onchange = junctions_filter;
        document.querySelector('#junction-psi-filter').onchange = junctions_filter;
        document.querySelector('#junction-lr-psi-filter').onchange = junctions_filter;
        document.querySelector('#reads-greater-than').oninput = junctions_filter;
        document.querySelector('#reads-less-than').oninput = junctions_filter;
        document.querySelector('#lr-reads-greater-than').oninput = junctions_filter;
        document.querySelector('#lr-reads-less-than').oninput = junctions_filter;
        document.querySelector('#psi-greater-than').oninput = junctions_filter;
        document.querySelector('#psi-less-than').oninput = junctions_filter;
        document.querySelector('#lr-psi-greater-than').oninput = junctions_filter;
        document.querySelector('#lr-psi-less-than').oninput = junctions_filter;

        junctions_filter();



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

                download_svg_elem(svg, `${grp}_${exp}_sg.svg`);

            }
        })
    }

    rescale_sg() {
        window.addEventListener('click', e => {
            if (e.target.classList.contains('splice-graph-rescale')) {
                const sg = e.target.closest('.splice-graph');

                this.sgs.scaling_transcript = sg.transcript;
                this.sgs.update(250);

            }
        })
    }

}

