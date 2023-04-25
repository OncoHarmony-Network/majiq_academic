function perform_swap(el1, el2, transform){

    // this swaps the DOM position of the nodes only
    var aparent = el1.parentNode;
    var asibling = el1.nextSibling === el2 ? el1 : el1.nextSibling;
    el2.parentNode.insertBefore(el1, el2);
    aparent.insertBefore(el2, asibling);

    // this part does the visual swap (because everything is aligned with transform, x, etc)
    // el2 will be static, we can trust the x values from it
    // el1 could be anywhere, so we use its origin x value
    const data_x = el1.getAttribute('data-orig-x');

    el1.setAttribute('data-orig-x', el2.getAttribute('data-orig-x'));
    el1.setAttribute('data-x', '');
    if (transform === true) el1.setAttribute('transform', el2.getAttribute('transform'));
    el1.setAttribute('x', el2.getAttribute('x'));
    el2.setAttribute('data-orig-x', data_x);
    el2.setAttribute('data-x', '');
    el2.setAttribute('x', data_x);
    if (transform === true) el2.setAttribute('transform', `translate(${data_x})`);

    if (transform !== true){
        el2.setAttribute('transform', '');
        el1.setAttribute('transform', '');
    }
}

function dotp(x,y) {
    function dotp_sum(a,b) { return a + b; }
    function dotp_times(a,i) { return x[i] * y[i]; }
    if (x.length != y.length)
        throw "can't find dot product: arrays have different lengths";
    return x.map(dotp_times).reduce(dotp_sum,0);
}

function expectation_value(bins){
    var step = 1/bins.length;
    var cur = step / 2;
    var comp_arr = [];
    while(cur < 1){
        comp_arr.push(cur);
            cur += step;
    }
    return dotp(comp_arr, bins);
}

function draw_pairwise_plot(main_svg){

        function pairwise_significance(v1_idx, v2_idx){
            const i1 = parseInt($(main_svg.node()).find(`.violin`).eq(v1_idx).attr('data-group-idx'));
            const i2 = parseInt($(main_svg.node()).find(`.violin`).eq(v2_idx).attr('data-group-idx'));
            let comp_val;

            // we want the row index to always be larger than the col index to get the low half of the heatmap
            if(i1 > i2){
                comp_val = parseFloat($(main_svg.node()).closest('tr').find('.heat-map g').eq(i1).find('rect').eq(i2).attr('data-value'));
            }else{
                comp_val = parseFloat($(main_svg.node()).closest('tr').find('.heat-map g').eq(i2).find('rect').eq(i1).attr('data-value'));
            }

            if(comp_val < 0){
                return null;
            }
            //comp_val = Math.log(comp_val);
            //comp_val = 10 ** comp_val;



            if(comp_val > 0.05){
                return null
            }else{
                if(comp_val > 0.001){
                    return "*"
                }else{
                    var exp = 4;
                    var stars = '*';
                    while(stars !== '***'){
                        stars = stars + '*';
                        exp += 1;
                        if(comp_val > (1/(10 ** exp))){
                            return stars;
                        }else if(exp > 100){
                            return "Overflow!";
                        }
                    }
                    return stars;
                }
            }
        }

        // each bar needs to be drawn from the inside out, with max height being the maximum number of
        // plots which are away to the left or the right of the chosen
        // we will minimize vertical distance by using the same height for plots to the left and the right
        // when applicable

        const cur_index = $(main_svg.node()).find('[data-checked="true"]').parent().index();
        main_svg.select('.overbars').remove();
        if(cur_index === -1){
            return;
        }

        // get violin plots center lines
        var x_vals = [];
        main_svg.selectAll(".violin").each(function() {
            x_vals.push(parseInt(this.getAttribute('x')));
        });


        var overbars = main_svg.append('svg')
            .attr('class', 'overbars')
            .style('overflow', 'visible')
            .attr('x', x_vals[0])

        // for lines skipped due to no correlation, we subtract this index from the main index
        var y_offset1 = 0;

        // draw lines to the left
        $.each(x_vals.slice(0, cur_index).reverse(), function(i, x){
            const significance = pairwise_significance(cur_index, cur_index - i - 1);
            if(significance === null){
                y_offset1 ++;
            }else{
                //horizontal
                overbars.append('line')
                    .attr('stroke', 'black')
                    .attr('x1', x_vals[cur_index])
                    .attr('y1', -(((i-y_offset1) * 15) + 15))
                    .attr('x2', x)
                    .attr('y2', -(((i-y_offset1) * 15) + 15))
                //vertical
                overbars.append('line')
                    .attr('stroke', 'black')
                    .attr('x1', x)
                    .attr('y1', -(((i-y_offset1) * 15) + 15))
                    .attr('x2', x)
                    .attr('y2', 0)
                // 'data point'
                overbars.append('text')
                    .text(significance)
                    .attr('font-size', 12)
                    .attr('text-anchor', 'middle')
                    .attr('x', x + (Math.abs(x - x_vals[cur_index]) / 2))
                    .attr('y', -(((i-y_offset1) * 15) + 17))
            }

        });

        var y_offset2 = 0;
        // draw lines to the right
        $.each(x_vals.slice(cur_index + 1, x_vals.length), function(i, x){
            const significance = pairwise_significance(cur_index, cur_index + i + 1);
            if(significance === null){
                y_offset2 ++;
            }else {
                //horizontal
                overbars.append('line')
                    .attr('stroke', 'black')
                    .attr('x1', x_vals[cur_index])
                    .attr('y1', -(((i - y_offset2) * 15) + 15))
                    .attr('x2', x)
                    .attr('y2', -(((i - y_offset2) * 15) + 15))
                //vertical
                overbars.append('line')
                    .attr('stroke', 'black')
                    .attr('x1', x)
                    .attr('y1', -(((i - y_offset2) * 15) + 15))
                    .attr('x2', x)
                    .attr('y2', 0)
                // 'data point'
                overbars.append('text')
                    .text(significance)
                    .attr('font-size', 12)
                    .attr('text-anchor', 'middle')
                    .attr('x', x_vals[cur_index] + (Math.abs(x - x_vals[cur_index]) / 2))
                    .attr('y', -(((i - y_offset2) * 15) + 17))
            }
        });

        //final vertical line in the 'center'
        //height determined by either items to the left or right of the selected item having more commections
        const max_y = Math.max((cur_index-y_offset1), x_vals.length - cur_index - 1 - y_offset2) * 15;
        overbars.append('line')
                .attr('stroke', 'black')
                .attr('x1', x_vals[cur_index])
                .attr('y1', 0)
                .attr('x2', x_vals[cur_index])
                .attr('y2', -max_y)

        main_svg.attr('height', 160 + max_y);
        overbars.attr('y', max_y);
        $(main_svg.node()).children('g').attr('transform', `translate(40, ${5 + max_y})`);
    }

class Violin {
    constructor(violin_data) {
        this.data = violin_data;
        this.violin_width = 75;
        this.violin_pad = 5;
        this.violin_height = 135;
        // this.x_axis_height = 125;
        this.x_axis_height = 20;
        this.y_axis_width = 40;
        this.top_padding = 5;
    };

    get svg_height() {
        return this.violin_height + this.x_axis_height + this.top_padding
    }

    get svg_width() {
        return this.y_axis_width + this.violin_count * (this.violin_width + this.violin_pad)
    }

    psi(svg) {
        const data = this.data;

        this.violin_count = data.junctions.length;
        svg.setAttribute('height', this.svg_height);
        svg.setAttribute('width', this.svg_width);

        const group = svg.dataset.group;
        const violin_data = data.group_bins[group];
        const color = new Colors();
        const g = d3.select(svg)
            .append('g')
            .attr('transform', `translate(${this.y_axis_width}, ${this.top_padding})`);

        const hist = g
            .append('g')
            .attr('class', 'histograms');

        this.draw_histograms(hist, violin_data);

        hist
            .selectAll('.violin')
            // .attr('stroke', (d, i) => color.brewer(i))
            .attr('stroke', null)
            .attr('stroke-width', 1)
            .attr('fill', (d, i) => color.brewer(i))
            .attr('fill-opacity', 1);

        this.draw_x_axis(g, data.group_means[group]);
        this.draw_psi_y_axis(g);
        this.box_plots(g, data.group_bins[group])
    }

    multipsi(svg) {

        this.violin_count = this.data.group_names.length;
        svg.setAttribute('height', this.svg_height);
        svg.setAttribute('width', this.svg_width);

        // const junc_idx = svg.closest('tr').dataset.junctionIndex;
        const junc_idx = this.data.junction_idx;
        const color = new Colors().brewer(junc_idx);

        const bins = this.data.group_bins[junc_idx];


        if(junc_idx === 0)
            svg.setAttribute('height', this.svg_height + 10);

        const g = d3.select(svg)
                .append('g');

        if(junc_idx === 0)
            g.attr('transform', `translate(${this.y_axis_width}, ${this.top_padding+10})`);
        else
            g.attr('transform', `translate(${this.y_axis_width}, ${this.top_padding})`);

        const hist = g
            .append('g')
            .attr('class', 'histograms');

        this.draw_histograms(hist, bins);

        hist
            .attr('stroke', color)
            .attr('stroke-width', 1)
            .attr('fill', color)
            .attr('fill-opacity', .1);



        //this.swarm(g, color);

        if(junc_idx === 0)
            this.draw_names_above(g, this.data.group_names);
        this.draw_x_axis(g, this.data.group_means[junc_idx]);
        this.draw_psi_y_axis(g);
        this.box_plots(g, this.data.group_bins[junc_idx])


    }

    deltapsi(svg) {
        const data = this.data;

        this.violin_count = data.junctions.length;
        svg.setAttribute('height', this.svg_height);
        svg.setAttribute('width', this.svg_width);

        const color = new Colors();
        const g = d3.select(svg)
            .append('g')
            .attr('transform', `translate(${this.y_axis_width}, ${this.top_padding})`);

        const hist = g
            .append('g')
            .attr('class', 'histograms');

        this.draw_histograms(hist, data.bins);

        hist
            .selectAll('.violin')
            // .attr('stroke', (d, i) => color.brewer(i))
            .attr('stroke', null)
            .attr('stroke-width', 1)
            .attr('fill', (d, i) => color.brewer(i))
            .attr('fill-opacity', 1);

        // trying to read the user-selected threshold value, first from the selection box if it exists
        var absEdPSI_Threshold_elem = $('#dpsi_threshold');
        var absEdPSI_Threshold_value = null;
        if(absEdPSI_Threshold_elem.length){
            absEdPSI_Threshold_value = parseFloat(absEdPSI_Threshold_elem.val());
        }else{
            // else in the get arguments
            let searchParams = new URLSearchParams(window.location.search);
            if(searchParams.has('absEdPSI_Threshold')){
                absEdPSI_Threshold_value = parseFloat(searchParams.get('absEdPSI_Threshold'));
            }
        }

        // if we found a threshold value specified, draw the horizontal lines
        if(absEdPSI_Threshold_value !== null && absEdPSI_Threshold_value > 0.0){
            this.draw_horizontal_line(g, (this.violin_height / 2) + ((this.violin_height / 2)
                * absEdPSI_Threshold_value), '#9a9a9a');
            this.draw_horizontal_line(g, (this.violin_height / 2) - ((this.violin_height / 2)
                * absEdPSI_Threshold_value), '#9a9a9a');
        }
        this.draw_x_axis(g, data.means.map(n => n.toFixed(3)));
        this.draw_dpsi_y_axis(g);
        this.box_plots(g, data.bins)

    }

    heterogen(svg) {
        const data = this.data;


        
        this.violin_count = data.group_names.length;
        svg.setAttribute('height', this.svg_height);
        svg.setAttribute('width', this.svg_width);


        // const junc_idx = svg.closest('tr').dataset.junctionIndex;
        const junc_idx = data.junction_idx;
        const color = new Colors().brewer(junc_idx);
        const bins = data.mean_psi;
        const medians = data.median_psi;


        const g = d3.select(svg)
            .append('g')
            .attr('transform', `translate(${this.y_axis_width}, ${this.top_padding})`);


        const hist = g
            .append('g')
            .attr('class', 'histograms');

        this.draw_histograms(hist, bins, medians);

        hist
            .attr('stroke', color)
            .attr('stroke-width', 1)
            .attr('fill', color)
            .attr('fill-opacity', .1);

        this.draw_psi_y_axis(g);

        const g2 = d3.select(svg)
            .append('g')
            .attr('class', 'swarms')
            .attr('transform', `translate(${this.y_axis_width}, ${this.top_padding})`);

        this.swarm(g2, color);
        this.draw_x_axis(g, data.group_names);
        this.draw_view_icons(g, data.group_names);
        //this.pairwise_plot_triggers(g);
    }

    transform_plot(i) {
        return 'translate(' + i * (this.dim.group.width + this.dim.group.pad) + ')';
    };

    _get_parent(child, classname) {
        var parent = child.select(function() {
          var element = this;
          while (!d3.select(element).classed(classname))
            element = element.parentElement;
          return element;
        });
        return parent;
    };

    draw_histograms(g, bins, medians) {

        // removing data from bins...with no data!
        // (setting to an empty array makes it so nothing is shown on the violin plot)
        $.each(bins, function(i, bin){
            if(bin.every( (val, i, arr) => val === -1 )){
                bins[i] = [];
            }
        });

        const x = d3.scaleLinear()
            .range([0, this.violin_width / 2]);

        const y = d3.scaleLinear()
            .range([this.violin_height, 0]);

        const area = d3.area()
            .curve(d3.curveCatmullRom)
            // .defined(function (d) {
            //     if (d > (x.domain()[0]))
            //         return d;
            // })
            .x1(function (d) {
                return x(d);
            })
            .x0(function (d) {
                return -x(d);
            })
            .y(function (d, i) {
                return y(i);
            });


            var deltaX;
            const colWidth = this.violin_width + this.violin_pad;
            var self = this;

            var dragHandler = d3.drag()
                .on("start", function () {
                    var current = d3.select(this);
                    deltaX = current.attr("data-x") - d3.event.x;
                    current.style('cursor', 'grabbing');
                })
                .on("drag", function () {

                    var current = d3.select(this);
                    //for getting lsv_id of dragged element, this searches parents to find it
                    var parent = self._get_parent(current, "lsv-table");
                    const lsv_id = parent.attr('data-lsv-id');

                    d3.selectAll(`.lsv-table[data-lsv-id="${lsv_id}"] text[data-group-idx="${current.attr("data-group-idx")}"]`)
                        .attr("data-x", d3.event.x + deltaX)
                        .attr("transform", "translate(" + (d3.event.x + deltaX) + ")");
                    d3.selectAll(`.lsv-table[data-lsv-id="${lsv_id}"] path[data-group-idx="${current.attr("data-group-idx")}"]`)
                        .attr("data-x", d3.event.x + deltaX)
                        .attr("transform", (d, i, el) => `translate(${d3.event.x + deltaX + parseInt(el[0].getAttribute("x"))})`);
                    d3.selectAll(`.lsv-table[data-lsv-id="${lsv_id}"] g[data-group-idx="${current.attr("data-group-idx")}"]`)
                        .attr("data-x", d3.event.x + deltaX)
                        .attr("transform", (d, i, el) => `translate(${d3.event.x + deltaX + parseInt(el[0].getAttribute("x"))})`);
                    d3.selectAll(`.lsv-table[data-lsv-id="${lsv_id}"] svg[data-group-idx="${current.attr("data-group-idx")}"]`)
                        .attr("data-x", d3.event.x + deltaX)
                        .attr("x", (d, i, el) => d3.event.x + deltaX + parseInt(el[0].getAttribute("data-orig-x")));


                })
                .on("end", function (d, i) {

                    var current = d3.select(this);
                    var lsv_id = self._get_parent(current, "lsv-table").attr('data-lsv-id');
                    var prev_index = $(d3.selectAll(`.lsv-table[data-lsv-id="${lsv_id}"] path[data-group-idx="${current.attr("data-group-idx")}"]`).node()).index();
                    var finalIndex = (Math.floor(((d3.event.x + deltaX) + (colWidth/2)) / colWidth)) + prev_index;

                    if(finalIndex > self.violin_count - 1){
                        finalIndex = self.violin_count - 1;
                    }else if(finalIndex < 0){
                        finalIndex = 0;
                    }

                    // user released dragging, so we apply the swap
                    // even if not dragged far enough for swap, we swap "with itself" to reset the position
                    // we deal with quite a few element attributes here, mainly because we are swapping three different pieces
                    // (violin plot, swarm, x-axis text) which are each set up in their own way.
                    // attrs "x" and "transform" are used in the other plot generation functions to position the plots
                    // "x" is used by the x-axis and "transform" is used by swarm and violin
                    // "data-x" is used by the dragging logic to properly align the cursor with the element while dragging
                    // "data-orig-x" is the original position of an element before dragging, so that it can be 'snapped back
                    // to grid' after the drag is complete.


                    $.each($(`.lsv-table[data-lsv-id="${lsv_id}"] path[data-group-idx="${current.attr("data-group-idx")}"]`), function(i, el1){
                        var el2 = $(el1).parent().children()[finalIndex];
                        perform_swap(el1, el2, true)
                    })
                    $.each($(`.lsv-table[data-lsv-id="${lsv_id}"] svg[data-group-idx="${current.attr("data-group-idx")}"]`), function(i, el1){
                        var el2 = $(el1).parent().children()[finalIndex];
                        perform_swap(el1, el2)
                    })
                    $.each($(`.lsv-table[data-lsv-id="${lsv_id}"] g[data-group-idx="${current.attr("data-group-idx")}"]`), function(i, el1){
                        var el2 = $(el1).parent().children()[finalIndex];
                        perform_swap(el1, el2, true)
                    })
                    $.each($(`.lsv-table[data-lsv-id="${lsv_id}"] text[data-group-idx="${current.attr("data-group-idx")}"]`), function(i, el1){
                        var el2 = $(el1).parent().children()[finalIndex];
                        perform_swap(el1, el2)
                    })
                    current.style('cursor', 'ew-resize');

                });


        g.selectAll('.violin')

            .data(bins)
            .enter()
            .append('path')
            .style('cursor', 'ew-resize')
            .attr('class', 'violin')
            .attr('transform', (d, i) => `translate(${(this.violin_width + this.violin_pad) * (i + .5)})`)
            .attr("x", (d, i) => `${(this.violin_width + this.violin_pad) * (i + .5)}`)
            .attr("data-orig-x", (d, i) => `${(this.violin_width + this.violin_pad) * (i + .5)}`)
            .attr('d', function (d) {
                x.domain(d3.extent(d));
                y.domain([0, d.length - 1]);
                return area(d)
            })
            .attr('data-group-idx', (d, i) => i)
            .attr('data-expected', (d, i) => medians ? medians[i] : expectation_value(d));

        // console.log(g.selectAll('.violin'))
        // dragHandler(g.selectAll('.violin'))
        if($.inArray(view_type, ['multipsi', 'het']) !== -1){
            d3.selectAll(".violin").call(dragHandler);
        }


    }

    mean_psi(m_psi, junc_idx) {
        return m_psi.map(function (arr) {
            try {
                return arr[junc_idx]
            } catch (TypeError) {
                return []
            }
        });
    };

    swarm(svg, color) {
        const circle_radius = 3;
        const mu_psi = this.data.mu_psi;
        const experiment_names = this.data.experiment_names;

        var tool_tip = d3.select('.violin-tool-tip');
        if (tool_tip.empty()) {
            tool_tip = d3.select("body")
                .append("div")
                .attr('class', 'violin-tool-tip')
                .style("display", "none");
            tool_tip.append('div')
                .attr('class', 'sample');
            tool_tip.append('div')
                .attr('class', 'value')
        }

        var x = d3.scaleLinear()
            .domain([0, 1])
            .range([this.violin_height, 0]);

        var swarm_fn = d3.beeswarm()
            .distributeOn((d) => x(d))
            .radius(circle_radius)
            .orientation('vertical')
            .side('symetric');

        svg
            .selectAll('.swarm-group')
            .data(mu_psi)
            .enter()
            .append('g')
            .attr('class', 'swarm-group')
            .attr('data-group-idx', (d, i) => i)
            .attr('transform', (d, i) => `translate(${i * (this.violin_width + this.violin_pad)})`)
            .attr("x", (d, i) => `${i * (this.violin_width + this.violin_pad)}`)
            .attr("data-orig-x", (d, i) => `${i * (this.violin_width + this.violin_pad)}`)
            .selectAll('circle')
            .data(d => swarm_fn.data(d).arrange())
            .enter()
            .select(function(d) { return d.datum !== -1 ? this : null; })
            .append("circle")
            .attr('fill', color)
            .attr('stroke', null)
            .attr("cx", bee => (bee.x % ((this.violin_width - circle_radius) / 2)) + ((this.violin_width + this.violin_pad) / 2))
            .attr("cy", bee => bee.y)
            .attr("r", circle_radius)
            .attr('data-mu', d => d.datum)
            .attr('data-exp-name', (d, i, a) => {
                const grp_idx = a[i].closest('g').dataset.groupIdx;
                return experiment_names[grp_idx][i]
            })
    };

    translate_lsv_bins(lsv_bins) {
        const numSamples = 40;
        const binsSize = lsv_bins.length;
        let numCopies;
        let tmpBins = [];

        lsv_bins.forEach((b, i) => {
            numCopies = Math.round(numSamples * b);
            tmpBins = tmpBins.concat(new Array(numCopies).fill((1 / binsSize) / 2 + (i / binsSize)))
        });

        return tmpBins;
    };

    box_plots(svg_outer, data) {

        const svg = svg_outer
                    .append('g')
                    .attr('class', 'box-plots');

        data.forEach((d, i) => {
            if(d.length > 0) {
                const trans_d = this.translate_lsv_bins(d);

                const q = d3.scaleQuantile()
                    .domain([0, 1])
                    .range(trans_d);

                const y = d3.scaleLinear()
                    .domain([0, 1])
                    .range([this.violin_height, 0]);

                const x = d3.scaleLinear()
                    .domain([0, 1])
                    .range([0, this.violin_width + this.violin_pad]);

                const g = svg.append('g')
                    .attr('transform', `translate(${x(i)})`)
                    .attr("data-orig-x", x(i))
                    .attr("x", x(i))
                    .attr("data-x", x(i))
                    .attr("data-group-idx", i);

                g
                    .selectAll('.h-line')
                    .data([.05, .5, .95].map(d => q(d)))
                    .enter()
                    .append('line')
                    .attr('stroke', 'black')
                    .attr('class', 'h-line')
                    .attr('x1', x(.4))
                    .attr('x2', x(.6))
                    .attr('y1', d => y(d))
                    .attr('y2', d => y(d));

                g
                    .append('rect')
                    .attr('stroke-width', 0)
                    .attr('width', x(.55) - x(.45))
                    .attr('height', y(q(.25)) - y(q(.75)))
                    .attr('x', (d, i, a) => x(.5) - (a[i].getAttribute('width') / 2))
                    .attr('y', y(q(.75)));

                g
                    .append('line')
                    .attr('stroke', 'black')
                    .attr('x1', x(.5))
                    .attr('x2', x(.5))
                    .attr('y1', y(q(.05)))
                    .attr('y2', y(q(.95)));

                g
                    .append('circle')
                    .attr('stroke', 'black')
                    .attr('fill', 'white')
                    .attr("cx", x(.5))
                    .attr("cy", y(d3.mean(trans_d)))
                    .attr("r", 3);
            }
        })
    };

    draw_horizontal_line(svg, y, color) {
        svg
            .append("line")
            .attr("x1", 0)
            .attr("x2", this.svg_width)
            .attr("y1", y)
            .attr("y2", y)
            .style("stroke-dasharray","5,5")
            .style("stroke", color);
    }
    
    draw_view_icons(svg, x_axis_data) {

  }

    draw_names_above(svg, x_axis_data) {
        svg
            .append('g')
            .attr('class', 'x-axis')
            .attr('transform', 'translate(0, -5)')
            .selectAll('text')
            .data(x_axis_data)
            .enter()
            .append('text')
            .attr('y', this.svg_height - this.x_axis_height + 6)
            .attr('font-size', 12)
            .text(d => {
                try {
                    return parseFloat(d.toPrecision(3))
                } catch (TypeError) {
                    const max_length = 20;
                    if (d.length > max_length) {
                        d = d.slice(0, max_length - 3) + '...'
                    }
                    return d
                }
            })
            .each((d, i, a) => {
                const el = a[i];
                if (d.length > 20) {
                    el.setAttribute('x', (this.violin_width + this.violin_pad) * (i + .45));
                    el.setAttribute("data-x", (this.violin_width + this.violin_pad) * (i + .45));
                    el.setAttribute("data-orig-x", (this.violin_width + this.violin_pad) * (i + .45));
                    el.setAttribute('y', 0);
                    el.setAttribute('transform', `rotate(90,${a[i].getAttribute('x')},${a[i].getAttribute('y')})`);
                    el.setAttribute('text-anchor', 'left');

                } else {
                    el.setAttribute('x', (this.violin_width + this.violin_pad) * (i + .5));
                    el.setAttribute("data-x", (this.violin_width + this.violin_pad) * (i + .5));
                    el.setAttribute("data-orig-x", (this.violin_width + this.violin_pad) * (i + .45));
                    el.setAttribute('y', 0);
                    el.setAttribute('text-anchor', 'middle');
                }
                el.setAttribute('data-group-idx', i)
            })
    }

    draw_x_axis(svg, x_axis_data) {
        var self = this;
        var chain = svg
            .append('g')
            .attr('class', 'x-axis')
            .selectAll('text')
            .data(x_axis_data)
            .enter()
            .append('svg')
            .style('overflow', 'visible')
            .attr('class', 'x-axis-grp')

            if(isInArr(view_type, ['het'])) {
                chain = chain
                    .append("rect")
                    .attr('class', 'pairwise-check')
                    .attr('y', 7)
                    .attr('x', -4)
                    .attr('width', 8)
                    .attr('height', 8)
                    .attr('fill', 'white')
                    .attr('stroke', 'black')
                    .attr('stroke-width', 1)
                    .style("cursor", "pointer")
                    .select(d3_parent)
            }

            chain = chain
            .attr('y', this.svg_height - this.x_axis_height + 6)
            .append('text')
            .attr('font-size', 12)
            .attr('textLength', d =>{
                if(d.length > 10){
                    return "77px";
                }
            })
            .attr('lengthAdjust', "spacingAndGlyphs")
            .text(d => {
                try {
                    return parseFloat(d.toPrecision(3))
                } catch (TypeError) {
                    const max_length = 20;
                    if (d.length > max_length) {
                        d = d.slice(0, max_length - 3) + '...'
                    }
                    return d
                }
            })
            .select(function() { return this.parentNode; })
            .each((d, i, a) => {
                const el = a[i];
                // if (d.length > 7) {
                //     el.setAttribute('x', (this.violin_width + this.violin_pad) * (i + .45));
                //     el.setAttribute("data-x", (this.violin_width + this.violin_pad) * (i + .45));
                //     el.setAttribute("data-orig-x", (this.violin_width + this.violin_pad) * (i + .45));
                //     el.setAttribute('y', this.svg_height - this.x_axis_height + 6);
                //     el.setAttribute('transform', `rotate(90,${a[i].getAttribute('x')},${a[i].getAttribute('y')})`);
                //     el.setAttribute('text-anchor', 'left');
                //
                // } else {
                    el.setAttribute('x', (this.violin_width + this.violin_pad) * (i + .5));
                    el.setAttribute("data-x", (this.violin_width + this.violin_pad) * (i + .5));
                    el.setAttribute("data-orig-x", (this.violin_width + this.violin_pad) * (i + .5));
                    el.setAttribute('y', this.svg_height - this.x_axis_height + 10);
                    el.setAttribute('text-anchor', 'middle');
                //                }
                el.setAttribute('textLength', '40px')
                el.setAttribute('data-group-idx', i)
            })

        // upon clicking the hide button
        if(!isInArr(view_type, ['multipsi', 'het'])) {
            d3.selectAll("image.hide-btn").style('display', 'none');
        }

    }



    draw_psi_y_axis(svg) {
        var y = d3.scaleLinear().domain([0, 1]).range([this.violin_height, 0]);
        var height = this.violin_height / 2;
        var label_pad = -28;
        var axis = d3.axisLeft(y).ticks(3);

        const g = svg
            .append('g')
            .attr('class', 'y-axis');

        g
            .append('g')
            .call(axis);

        g
            .append('text')
            .text('E(Ψ)')
            .attr('font-size', 12)
            .attr('text-anchor', 'middle')
            .attr('transform', 'rotate(-90,' + label_pad + ',' + height + ')')
            .attr('y', height)
            .attr('x', label_pad)
    }

    draw_dpsi_y_axis(svg) {
        var y = d3.scaleLinear().domain([-1, 1]).range([this.violin_height, 0]);
        var height = this.violin_height / 2;
        var label_pad = -28;
        var axis = d3.axisLeft(y).ticks(3);

        const g = svg
            .append('g')
            .attr('class', 'y-axis');

        g
            .append('g')
            .call(axis);

        g
            .append('text')
            .text('E(ΔΨ)')
            .attr('font-size', 12)
            .attr('text-anchor', 'middle')
            .attr('transform', 'rotate(-90,' + label_pad + ',' + height + ')')
            .attr('y', height)
            .attr('x', label_pad)
    }


    x_axis(svg) {
        this.metadata()
            .then(data => {
                d3.select(svg)
                    .append('g')
                    .selectAll('text')
                    .data(data.group_names)
                    .enter()
                    .append('text')
                    .attr('transform', function (d, i) {
                        return v.transform_plot(i)
                    })
                    .attr('text-anchor', 'middle')
                    .attr('y', height)
                    .attr('x', width / 2)
                    .text(function (d) {
                        return d
                    })
            })
    }

    y_axis() {
        var y = d3.scaleLinear().domain([0, 1]).range([this.dim.group.height, 0]);
        var height = this.dim.group.height / 2 + this.dim.pad.top;
        var label_pad = this.dim.y_axis.label - 6;
        var axis = d3.axisLeft(y).ticks(3);

        var g = this.svg.append('g');

        g
            .append('g')
            .attr('transform', 'translate(' + (this.dim.y_axis.width + this.dim.y_axis.label) + ',' + this.dim.pad.top + ')')
            .call(axis);

        g
            .append('text')
            .text('E(PSI)')
            .attr('text-anchor', 'middle')
            .attr('transform', 'rotate(-90,' + label_pad + ',' + height + ')')
            .attr('y', height)
            .attr('x', label_pad)
    }
}
