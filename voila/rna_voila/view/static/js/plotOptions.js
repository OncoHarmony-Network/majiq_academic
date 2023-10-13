class PlotOptions {
    constructor(gene) {

        //this.exp_names = gene.experiment_names;


        if(group_order_override){
            this.grp_names = group_order_override;
        }else{
            if('sortable_group_names' in gene){
                this.grp_names = gene.sortable_group_names
            }else{
                this.grp_names = gene.group_names;
            }
        }

        if(group_display_name_override) {
            this.group_display_name_override = group_display_name_override;
        }else{
            this.group_display_name_override = {};
        }

        if(group_visibility) {
            this.group_visibility = group_visibility;
        }else{
            this.group_visibility = {};
        }

        this.violin_fixed_width = violin_fixed_width;


        this._init();

    }

    _init() {

        // when clicked, will toggle all of the other active menus / buttons off and enable only ours
        function hide_others(){
            d3.selectAll('.tools-menu')
                .classed('hide-tools-menu', true);
            d3.selectAll('.tools-menu-btn')
                .classed('pure-menu-active', false);

        }

        function show_ours(){
            document.getElementById('plot-options-box').classList.remove('hide-tools-menu');
            $('#plot-options-menu-btn').addClass('pure-menu-active')
        }

        // enable this menu when the button is clicked
        document.querySelector('#plot-options-menu-btn').onclick = (event) => {
            event.preventDefault();
            event.stopPropagation();
            const already_active = $('#plot-options-menu-btn').hasClass('pure-menu-active');
            hide_others();
            if(!already_active){
                show_ours();
            }
        };

        // handling click outside of our options hides them
        $("body").click(function(event){
            if($('#plot-options-menu-btn').hasClass('pure-menu-active')){
                hide_others();
            }
        });

        // (and click inside the options does not hide it)
        $("#plot-options-box").click(function (event) {
            event.stopPropagation();
        });

        // populate groups list with default names

        const self = this;
        $.each(this.grp_names, function(i, v){
            let dispname = v;
            let visible = 'checked';
            if(v in self.group_display_name_override){
                dispname = self.group_display_name_override[v];
            }
            if(v in self.group_visibility){
                visible = self.group_visibility[v] ? 'checked' : '';
            }
            const template = `
                <div class="drag">
                <dt><img class="handle" src="${base_url}/static/img/drag_blue.png" title="Re-order this group"> <span>${v}</span>
                <input type="checkbox" class="group-visible" ${visible} title="Display this group?"></dt>
                <dd><input class="group-name-override" value="${dispname}" title="Change the display name of this group"></dd>
                </div>`;
            $('#group-controls').append(template)

        })

        // set up group sortable
        var el = document.getElementById('group-controls');
        var sortable = Sortable.create(el, {
            draggable: ".drag",
            handle: ".handle",
            direction: 'vertical',
            animation: 150,
            onEnd: function (evt) {
                let elem_list = $('#group-controls dt span').map(function(){
                    return $.trim($(this).text());
                }).get();

                send_ajax(base_url + '/update-group-order', elem_list).then(ret => {
                    $('.lsv-table').DataTable().ajax.reload();
                })
            },
        });

        $('.group-visible').on('change', ($.debounce(700, function(){
            group_visibility = {};
            $('#group-controls .drag').each(function(i, v){
                const key = $(v).find('dt span').text();
                const value = $(v).find('dt .group-visible').prop('checked');
                group_visibility[key] = value;

            })

            send_ajax(base_url + '/update-group-visibility', group_visibility).then(ret => {
                $('.lsv-table').DataTable().ajax.reload();
            })
        })));

        $('.group-name-override').on('input', ($.debounce(700, function(){
            group_display_name_override = {};
            $('#group-controls .drag').each(function(i, v){
                const key = $(v).find('dt span').text();
                const value = $(v).find('dd input').val();
                group_display_name_override[key] = value;

            })

            send_ajax(base_url + '/update-group-display-names', group_display_name_override).then(ret => {
                $('.lsv-table').DataTable().ajax.reload();
            })
            //$('.lsv-table').DataTable().ajax.reload();
        })));

        // button for resetting these settings
        $('#resetGroupNamesButton').click(function(){
            send_ajax(base_url + '/reset-group-settings', {}).then(ret => {
                location.reload();
            })
        })

        // button for resetting visibility
        $('#showAllGroupNamesButton').click(function(){
            group_visibility = {};
            $.each(group_names, function(i, v){
                group_visibility[v] = true;
            })
            $('#group-controls .drag').each(function(i, v){
                $(v).find('dt .group-visible').prop('checked', true);
            })
            send_ajax(base_url + '/update-group-visibility', group_visibility).then(ret => {
                $('.lsv-table').DataTable().ajax.reload();
            })
        });

        $('#hideAllGroupNamesButton').click(function(){
            group_visibility = {};
            $.each(group_names, function(i, v){
                group_visibility[v] = false;
            })
            $('#group-controls .drag').each(function(i, v){
                $(v).find('dt .group-visible').prop('checked', false);
            })
            send_ajax(base_url + '/update-group-visibility', group_visibility).then(ret => {
                $('.lsv-table').DataTable().ajax.reload();
            })
        });

        $('#changeViolinFixedWidth').change(function(){
            violin_fixed_width = $('#changeViolinFixedWidth').prop('checked');

            send_ajax(base_url + '/update-violin-fixed-width', violin_fixed_width).then(ret => {
                $('.lsv-table').DataTable().ajax.reload();
            })
        });

        // detect changes

        // document.querySelectorAll('.tools-menu-btn')
        //     .forEach(l => l.onclick = (event) => {
        //         event.preventDefault();
        //         d3.selectAll('.plot-options.tools-menu')
        //             .classed('hide-tools-menu', true);
        //         document.querySelector('.lsv-tools.tools-menu').classList.toggle('hide-tools-menu');
        //     });


        // if($('.splice-graph-container .splice-graph').length === 0){
        //     setTimeout(function(){$('.splice-graph-form button').eq(0).trigger('click')}, 0);
        // }

        //
        //
        // // populate splice graph selector groups
        // d3.select('.groups select')
        //     .selectAll('option')
        //     .data(this.grp_names)
        //     .enter()
        //     .append('option')
        //     .text(d => {
        //         return d
        //     });
        //
        //
        // // populate splice graph experiments when group is changed
        // document.querySelector('#groups').onchange = (event) => {
        //     const group_name = this.grp_names[event.target.selectedIndex];
        //
        //     const shown_exps = Array.from(document.querySelectorAll('.splice-graph'))
        //         .filter(sg => sg.dataset.group === group_name)
        //         .map(sg => sg.dataset.experiment);
        //
        //     const exps = this.exp_names[event.target.selectedIndex]
        //         .filter(e => !shown_exps.includes(e));
        //
        //     const s = d3.select('.experiments select')
        //         .selectAll('option')
        //         .data(exps);
        //
        //     s.text(d => d);
        //
        //     s.enter()
        //         .append('option')
        //         .text(d => d);
        //
        //     s.exit().remove();
        // };
        //
        // // force change event to populate experiments on initial page load
        // SpliceGraphTools._populate_sg_form();
        //
        // // submit event for splice graph selector
        // document.querySelector('.splice-graph-form').onsubmit = (event) => {
        //     event.preventDefault();
        //     const f = event.target;
        //     const g = f.querySelector('#groups').value;
        //     const e = f.querySelector('#experiments').value;
        //     if (g && e) {
        //         f.querySelector('button').disabled = true;
        //         this.sgs.create(g, e);
        //         SpliceGraphTools._populate_sg_form();
        //         f.querySelector('button').disabled = false;
        //         junctions_filter();
        //         this.highlight_lsvs();
        //     }
        // };
        //
        // // toggle splice graph scale
        // document.querySelector('.toggle-scale').onclick = () => {
        //     document.querySelector('.splice-graph-container').classList.toggle('default-view');
        //     this.sgs.update(250);
        // };
        //
        // // zoom in on splice graph
        // document.querySelector('.zoom-in').onclick = () => {
        //     let zoom = parseInt(document.querySelector('.splice-graph-container').dataset.zoom);
        //     zoom++;
        //     document.querySelector('.splice-graph-container').dataset.zoom = zoom;
        //     this.sgs.update(250);
        // };
        //
        // // zoom out on splice graph
        // document.querySelector('.zoom-out').onclick = () => {
        //     let zoom = parseInt(document.querySelector('.splice-graph-container').dataset.zoom);
        //     if (zoom !== 1) {
        //         zoom--;
        //         document.querySelector('.splice-graph-container').dataset.zoom = zoom;
        //         this.sgs.update(250);
        //     }
        // };
        //
        // // reset zoom for splice graph
        // document.querySelector('.zoom-reset').onclick = () => {
        //     document.querySelector('.splice-graph-container').dataset.zoom = 1;
        //     this.sgs.update(250);
        // };
        //
        // // activate/deactivate junction reads filter
        // document.querySelector('#junction-reads-filter').onchange = (event) => {
        //     document.querySelectorAll('#reads-greater-than, #reads-less-than').forEach(el => el.disabled = !el.disabled);
        //     if (event.target.checked) {
        //         junctions_filter()
        //     } else {
        //         this.sgs.junctions_filter()
        //     }
        // };
        //
        // // adjust greater than and less than fields in junction filter
        // const junctions_filter = () => {
        //     if (document.querySelector('#junction-reads-filter').checked) {
        //         const gt = document.querySelector('#reads-greater-than').value;
        //         const lt = document.querySelector('#reads-less-than').value;
        //         this.sgs.junctions_filter(gt, lt)
        //     }
        // };
        //
        // document.querySelector('#reads-greater-than').oninput = junctions_filter;
        // document.querySelector('#reads-less-than').oninput = junctions_filter;
        // junctions_filter();


    }


    static _populate_form() {
        //document.querySelector('#groups').dispatchEvent(new Event('change'));
    }



}

