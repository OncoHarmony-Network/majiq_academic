class Violin {
    constructor(violin_data) {
        this.data = violin_data;
        this.violin_width = 50;
        this.violin_pad = 5;
        this.violin_height = 135;
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

        this.draw_x_axis(g, data.means.map(n => n.toFixed(3)));
        this.draw_dpsi_y_axis(g);
        this.box_plots(g, data.bins)

    }

    async heterogen(svg) {
        this.violin_count = this.data.group_names.length;
        svg.setAttribute('height', this.svg_height);
        svg.setAttribute('width', this.svg_width);

        // const junc_idx = svg.closest('tr').dataset.junctionIndex;
        const junc_idx = this.data.junction_idx;
        const color = new Colors().brewer(junc_idx);
        const bins = this.data.mean_psi;

        const g = d3.select(svg)
            .append('g')
            .attr('transform', `translate(${this.y_axis_width}, ${this.top_padding})`);

        const hist = g
            .append('g')
            .attr('class', 'histograms');

        this.draw_histograms(hist, bins);


        hist
            .selectAll('.violin')
            .attr('stroke', color)
            .attr('stroke-width', 1)
            .attr('fill', color)
            .attr('fill-opacity', .1);

        this.draw_psi_y_axis(g);

        this.swarm(g, color);
        this.draw_x_axis(g, this.data.group_names);
    }

    transform_plot(i) {
        return 'translate(' + i * (this.dim.group.width + this.dim.group.pad) + ')';
    };

    draw_histograms(g, bins) {
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

        g.selectAll('.violin')
            .data(bins)
            .enter()
            .append('path')
            .attr('class', 'violin')
            .attr('transform', (d, i) => `translate(${(this.violin_width + this.violin_pad) * (i + .5)})`)
            .attr('d', function (d) {
                x.domain(d3.extent(d));
                y.domain([0, d.length - 1]);
                return area(d)
            });
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
            .attr('data-group-index', (d, i) => i)
            .attr('transform', (d, i) => 'translate(' + i * (this.violin_width + this.violin_pad) + ')')
            .selectAll('circle')
            .data(d => swarm_fn.data(d).arrange())
            .enter()
            .append("circle")
            .attr('fill', color)
            .attr('stroke', null)
            .attr("cx", bee => bee.x + ((this.violin_width + this.violin_pad) / 2))
            .attr("cy", bee => bee.y)
            .attr("r", circle_radius)
            .attr('data-mu', d => d.datum)
            .attr('data-exp-name', (d, i, a) => {
                const grp_idx = a[i].closest('g').dataset.groupIndex;
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

    box_plots(svg, data) {
        data.forEach((d, i) => {
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
                .attr('transform', `translate(${x(i)})`);

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
        })
    };

    draw_x_axis(svg, x_axis_data) {
        svg
            .append('g')
            .attr('class', 'x-axis')
            .selectAll('text')
            .data(x_axis_data)
            .enter()
            .append('text')
            .attr('text-anchor', 'middle')
            .attr('y', this.svg_height - 6)
            .attr('x', (d, i) => (this.violin_width + this.violin_pad) * (i + .5))
            .attr('font-size', 12)
            .text(d => {
                try {
                    return parseFloat(d.toPrecision(3))
                } catch (TypeError) {
                    return d
                }
            })
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
