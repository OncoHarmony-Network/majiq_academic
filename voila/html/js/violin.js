var Violin = function (db, el) {
    this.db = db;
    this.junc_idx = el.closest('tr').dataset.junctionIndex;
    this.lsv_id = el.closest('table').dataset.lsvId;
    this.type = el.dataset.type;

    this.dim = {
        x_axis: {
            height: 18
        },
        y_axis: {
            width: 25,
            label: 18
        },
        pad: {
            top: 6,
            bottom: 6
        },
        group: {
            height: 114,
            width: 80,
            pad: 6
        }
    };

    this.svg = d3.select(el);
    this.svg.selectAll('*').remove();

    this.size_svg();

    this.plot = this.svg
        .append('g')
        .attr('transform', 'translate(' + (this.dim.y_axis.width + this.dim.group.pad + this.dim.y_axis.label) + ',' + this.dim.pad.top + ')');

    var color = new Colors().brewer(this.junc_idx);

    this.histograms(color);
    this.x_axis();
    this.y_axis();
    if (this.type === 'box')
        this.box_plots();
    if (this.type === 'swarm')
        this.swarm(color)
};

Violin.prototype.size_svg = function () {
    var v = this;
    this.db.get('metadata').then(function (data) {
        var groups_num = data.group_names.length;
        var height = v.dim.group.height + v.dim.pad.top + v.dim.pad.bottom + v.dim.x_axis.height;
        var width = (v.dim.group.width * groups_num) + (v.dim.group.pad * (groups_num + 1)) + v.dim.y_axis.width + v.dim.y_axis.label;
        v.svg.attr('height', height).attr('width', width);
    });
};

Violin.prototype.transform_plot = function (i) {
    return 'translate(' + i * (this.dim.group.width + this.dim.group.pad) + ')';
};

Violin.prototype.histograms = function (color) {
    var junc_idx = this.junc_idx;
    var width = this.dim.group.width;
    var height = this.dim.group.height;
    var pad = this.dim.group.pad;
    var plot = this.plot;
    var v = this;

    this.db.get(this.lsv_id).then(function (data) {

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

        var mean_psi = v.mean_psi(data, junc_idx);

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

Violin.prototype.mean_psi = function (data, junc_idx) {
    return data.mean_psi.map(function (arr) {
        try {
            return arr[junc_idx]
        } catch (TypeError) {
            return []
        }
    });
};


Violin.prototype.swarm = function (color) {
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


Violin.prototype.box_plots = function () {
    var mean_psi_fn = this.mean_psi;
    var junc_idx = this.junc_idx;
    var height = this.dim.group.height;
    var width = this.dim.group.width;
    var plot = this.plot;
    var r = 3;
    var v = this;

    var translateLsvBins = function (lsvBins) {
        var numSamples = 40;
        var tmpBins = [];
        var binsSize = lsvBins.length;
        var numCopies;
        lsvBins.forEach(function (b, i) {
            numCopies = Math.round(numSamples * b);
            tmpBins = tmpBins.concat(new Array(numCopies).fill((1 / binsSize) / 2 + (i / binsSize)))
        });
        return tmpBins;
    };

    this.db.get(this.lsv_id).then(data => {
        this.mean_psi(data, junc_idx).forEach((d, i) => {
            if (d.length) {
                var trans_d = translateLsvBins(d);

                var q = d3.scaleQuantile()
                    .domain([0, 1])
                    .range(trans_d);

                var y = d3.scaleLinear()
                    .domain([0, 1])
                    .range([height, 0]);

                var x = d3.scaleLinear()
                    .domain([0, 1])
                    .range([0, width]);

                var g = plot.append('g')
                    .attr('transform', v.transform_plot(i));

                g
                    .selectAll('.h-line')
                    .data([.05, .5, .95].map(function (d) {
                        return q(d)
                    }))
                    .enter()
                    .append('line')
                    .attr('stroke', 'black')
                    .attr('class', 'h-line')
                    .attr('x1', x(.4))
                    .attr('x2', x(.6))
                    .attr('y1', function (d) {
                        return y(d)
                    })
                    .attr('y2', function (d) {
                        return y(d)
                    });


                g
                    .append('rect')
                    .attr('stroke-width', 0)
                    .attr('width', x(.55) - x(.45))
                    .attr('height', y(q(.25)) - y(q(.75)))
                    .attr('x', function () {
                        return x(.5) - (this.getAttribute('width') / 2)
                    })
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
                    .attr("r", r);


            }
        })
    })
};

Violin.prototype.x_axis = function () {
    var plot = this.plot;
    var height = this.dim.group.height + this.dim.pad.top + (this.dim.x_axis.height / 2);
    var width = this.dim.group.width;
    var v = this;
    this.db.get('metadata').then(function (data) {
        plot
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
};

Violin.prototype.y_axis = function () {
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
};