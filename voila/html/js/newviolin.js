class Violin {
    constructor(db) {
        this.db = db;
        this.violin_width = 50;
        this.violin_pad = 5;
        this.violin_height = 135;
        this.x_axis_height = 20;
        this.y_axis_width = 40;
        this.top_padding = 5;
        this._metadata = undefined;
    };

    get svg_height() {
        return this.violin_height + this.x_axis_height + this.top_padding
    }

    get svg_width() {
        return this.y_axis_width + this.violin_count * (this.violin_width + this.violin_pad)
    }

    get metadata() {
        if (this._metadata === undefined) {
            return this.db.get('metadata')
                .then(metadata => {
                    this._metadata = metadata;
                    return this._metadata
                })
        } else {
            console.log('cached');
            return new Promise(this._metadata)
        }
    }

    psi(svg) {
        this.db.get(svg.closest('tr').dataset.lsvId)
            .then(results => {
                d3.select(svg).selectAll('*').remove();
                this.violin_count = results.junctions.length;
                svg.setAttribute('height', this.svg_height);
                svg.setAttribute('width', this.svg_width);

                const group = svg.dataset.group;
                const violin_data = results.group_bins[group];
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

                this.draw_x_axis(g, results.group_means[group]);
                this.draw_psi_y_axis(g);
                this.box_plots(g, results.group_bins[group])
            })
    }

    deltapsi(svg) {
        this.db.get(svg.closest('tr').dataset.lsvId)
            .then(results => {
                d3.select(svg).selectAll('*').remove();
                this.violin_count = results.junctions.length;
                svg.setAttribute('height', this.svg_height);
                svg.setAttribute('width', this.svg_width);

                const color = new Colors();
                const g = d3.select(svg)
                    .append('g')
                    .attr('transform', `translate(${this.y_axis_width}, ${this.top_padding})`);

                const hist = g
                    .append('g')
                    .attr('class', 'histograms');

                this.draw_histograms(hist, results.bins);

                hist
                    .selectAll('.violin')
                    // .attr('stroke', (d, i) => color.brewer(i))
                    .attr('stroke', null)
                    .attr('stroke-width', 1)
                    .attr('fill', (d, i) => color.brewer(i))
                    .attr('fill-opacity', 1);

                this.draw_x_axis(g, results.means);
                this.draw_dpsi_y_axis(g);
                this.box_plots(g, results.bins)
            })
    }

    heterogen(el) {
        this.junc_idx = el.closest('tr').dataset.junctionIndex;
        this.lsv_id = el.closest('table').dataset.lsvId;
        this.type = el.dataset.type;


        this.svg = d3.select(el);
        this.svg.selectAll('*').remove();

        this.size_svg();

        this.plot = this.svg
            .append('g')
            .attr('transform', 'translate(' + (this.dim.y_axis.width + this.dim.group.pad + this.dim.y_axis.label) + ',' + this.dim.pad.top + ')');

        var color = new Colors().brewer(this.junc_idx);

        this.het_histograms(color);
        this.x_axis();
        this.y_axis();
        if (this.type === 'box')
            this.box_plots();
        if (this.type === 'swarm')
            this.swarm(color)
    }

    size_svg(svg) {
        this.metadata()
            .then(results => {
                const groups_num = results.group_names.length;
                const height = this.dim.group.height + this.dim.pad.top + this.dim.pad.bottom + this.dim.x_axis.height;
                const width = (this.dim.group.width * groups_num) + (this.dim.group.pad * (groups_num + 1)) + this.dim.y_axis.width + this.dim.y_axis.label;
                svg.setAttribute('height', height);
                svg.setAttribute('width', width)
            });
    };

    transform_plot(i) {
        return 'translate(' + i * (this.dim.group.width + this.dim.group.pad) + ')';
    };

    draw_histograms(g, violin_data) {
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
            .data(violin_data)
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

    het_histograms(color) {
        var junc_idx = this.junc_idx;
        var width = this.dim.group.width;
        var height = this.dim.group.height;
        var pad = this.dim.group.pad;
        var plot = this.plot;


        this.db.get(this.lsv_id).then(data => {

            var x = d3.scaleLinear()
                .rangeRound([0, width / 2]);

            var y = d3.scaleLinear()
                .range([height, 0]);

            var area = d3.area()
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

            var mean_psi = this.mean_psi(data, junc_idx);

            plot
                .selectAll('.violin')
                .data(mean_psi)
                .enter()
                .append('path')
                .attr('class', 'violin')
                .attr('transform', function (d, i) {
                    return 'translate(' + ((width / 2) + (i * (width + pad))) + ')'
                })
                .attr('stroke', color)
                .attr('stroke-width', 1)
                .attr('fill', color)
                .attr('fill-opacity', .05)
                .attr('d', function (d) {
                    x.domain(d3.extent(d));
                    y.domain([0, d.length - 1]);
                    return area(d)
                });
        })
    };

    mean_psi(data, junc_idx) {
        return data.mean_psi.map(function (arr) {
            try {
                return arr[junc_idx]
            } catch (TypeError) {
                return []
            }
        });
    };

    swarm(color) {
        var height = this.dim.group.height;
        var width = this.dim.group.width;
        var padding = this.dim.group.pad;
        var svg = this.plot;
        var junc_idx = this.junc_idx;
        var circle_radius = 2;

        this.db.allDocs({
            keys: ['metadata', this.lsv_id],
            include_docs: true
        }, function (err, response) {
            var meta = response.rows[0].doc;
            var data = response.rows[1].doc;

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
                .range([height, 0]);

            var swarm_fn = d3.beeswarm()
                .distributeOn(function (d) {
                    return x(d);
                })
                .radius(circle_radius)
                .orientation('vertical')
                .side('symetric');

            svg
                .selectAll('.swarm-group')
                .data(data.mu_psi[junc_idx])
                .enter()
                .append('g')
                .attr('class', 'swarm-group')
                .attr('data-group-index', function (d, i) {
                    return i
                })
                .attr('transform', function (d, i) {
                    return 'translate(' + i * (width + padding) + ')'
                })
                .selectAll('circle')
                .data(function (d) {
                    return swarm_fn
                        .data(d)
                        .arrange();
                })
                .enter()
                .append("circle")
                .attr('fill', color)
                .attr('stroke', null)
                .attr("cx", function (bee) {
                    return bee.x + (width / 2);
                })
                .attr("cy", function (bee) {
                    return bee.y;
                })
                .attr("r", circle_radius)
                .attr('data-mu', function (d) {
                    return d.datum
                })
                .on("mouseover", function (d, i) {
                    d3.select(this).style('fill', 'orange');
                    tool_tip.selectAll('.value').text(this.getAttribute('data-mu'));
                    var group_idx = parseInt(this.parentNode.getAttribute('data-group-index'));
                    var exp_names = meta.experiment_names[group_idx].reduce(function (acc, curr) {
                        if (!curr.includes('Combined'))
                            acc.push(curr);
                        return acc
                    }, []);
                    tool_tip.selectAll('.sample').text(exp_names[i]);
                    tool_tip.style("display", "block");
                })
                .on("mousemove", function () {
                    tool_tip.style("top", (event.pageY - 35) + "px").style("left", (event.pageX + 10) + "px");
                })
                .on("mouseout", function () {
                    d3.select(this).style('fill', '');
                    tool_tip.style("display", "none");
                });
        });
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
            .text(d => parseFloat(d.toPrecision(3)))
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
