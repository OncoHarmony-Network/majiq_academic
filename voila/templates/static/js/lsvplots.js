var LSVPlots = function (el) {
    this.db = db;
    this.height = 126;
    this.histo_width = 80;
    this.left_padding = 41;
    this.top_padding = 10;
    this.bottom_padding = 25;
    this.svg_height = this.height + this.top_padding + this.bottom_padding;
    this.el = el;
    this.svg = d3.select(el);
    this.lsv_id = el.getAttribute('data-lsv-id');
    this.junc_idx = parseInt(el.getAttribute('data-junction-index'));
    this.group_names = el.getAttribute('data-group-names') ? el.getAttribute('data-group-names').split(',') : [];


};

LSVPlots.prototype.heterogen = function () {
    var junc_idx = this.junc_idx;
    var lsvp = this;

    db.get(this.lsv_id).then(function (data) {
        lsvp.bins = data.mean_psi.map(function (arr) {
            try {
                return arr[junc_idx]
            } catch (TypeError) {
                return []
            }
        });

        lsvp.mu_psi = data.mu_psi[junc_idx];

        lsvp.size_svg()
            .violins()
            .box_plots()
            .axes([0, 1], [0, .5, 1], 'E(PSI)')
            // .swarm()
    })


    // this.db.get(this.lsv_id).then(function (data) {

    //     var group_width = bins.length * bp.histo_width;
    //     bp.svg.attr('group_width', group_width + bp.left_padding);
    //     bp.violins(bins, new Colors().brewer_het(junc_idx), group_names);
    //     bp.axes([0, 1], [0, .5, 1], group_width, 'E(PSI)');
    //     bp.boxPlots();
    //     bp.swarm(group_width)
    // }
    // )
};

LSVPlots.prototype.box_plots = function () {
    var lsvp = this;
    this.svg.selectAll('.violin')
        .each(function (d) {
            d = translateLsvBins(d);
            if (d.length)
                d3.select(this)
                    .append('g')
                    .attr('class', 'box-plot')
                    .attr('stroke', 'black')
                    .attr('fill', 'black')
                    .attr('transform', 'translate(0,' + lsvp.histo_width / 2 + ')')
                    .voilinBoxPlots(d, lsvp.height);
        });
    return this;
};

LSVPlots.prototype.size_svg = function () {
    var lsvp = this;
    var width = this.bins.length * lsvp.histo_width;
    lsvp.svg.attr('height', lsvp.svg_height);
    lsvp.svg.attr('width', width + lsvp.left_padding);
    return this;
};

LSVPlots.prototype.psi = function () {
    var lsv_id = this.el.getAttribute('data-lsv-id');
    var group = this.el.getAttribute('data-group');
    var bp = this;


    this.db.get(lsv_id).then(function (data) {
        var bins = data.group_bins[group];
        var width = bins.length * bp.histo_width;
        var means_rounded = data.group_means_rounded[group];

        bp.svg = d3.select(el)
            .attr('width', width + bp.left_padding)
            .attr('height', bp.svg_height);
        if (bins)
            bp.violins(bins, new Colors().brewer, means_rounded);
        bp.axes([0, 1], [0, .5, 1], width);
    });
};


LSVPlots.prototype.delta_psi = function () {
    var el = this.el;
    var lsv_id = el.getAttribute('data-lsv-id');
    var bp = this;

    this.db.get(lsv_id).then(function (data) {
        var width = data.bins.length * bp.histo_width;

        bp.svg = d3.select(el)
            .attr('width', width + bp.left_padding)
            .attr('height', bp.svg_height);
        bp.violins(data.bins, new Colors().brewer.data.means_rounded);
        bp.axes([-1, 1], [-1, 0, 1], width);
    });
};

translateLsvBins = function (lsvBins) {
    var numSamples = 100;
    var tmpBins = [];
    var binsSize = lsvBins.length;
    var numCopies;
    lsvBins.forEach(function (b, i) {
        numCopies = Math.round(numSamples * b);
        tmpBins = tmpBins.concat(new Array(numCopies).fill((1 / binsSize) / 2 + (i / binsSize)))
    });
    return tmpBins;
};


LSVPlots.prototype.axes = function (domain, tick_values, y_label) {
    var bp = this;
    var yScale = d3.scaleLinear().domain(domain).range([this.height, 0]);
    var width = this.bins.length * this.histo_width;
    var g = this.svg
        .append('g')
        .classed('axes', true)
        .attr('transform', 'translate(' + bp.left_padding + ',' + bp.top_padding + ')');

    g
        .call(d3.axisLeft(yScale).tickValues(tick_values))
        .append('text')
        .classed('y-axis-label', true)
        .text(y_label)
        .attr('text-anchor', 'middle')
        .attr('font-size', 12)
        .attr('font-family', 'sans-serif')
        .attr('fill', 'black')
        .attr('transform', 'rotate(-90) translate(-' + bp.height / 2 + ', -30)');
    g
        .append('line')
        .attr('class', 'x-axis')
        .attr('x1', 0)
        .attr('x2', width)
        .attr('y1', yScale(0) + .5)
        .attr('y2', yScale(0) + .5)
        .attr('stroke', 'black');

    return this
};

LSVPlots.prototype.violins = function () {
    var lsvp = this;
    var junc_idx = this.junc_idx;
    var colors = new Colors().brewer_het(junc_idx);
    var means_rounded = this.group_names;

    lsvp.svg
        .selectAll('.violin')
        .data(this.bins)
        .enter()
        .append('g')
        .attr('class', 'violin')
        .each(function (d, i) {
            d = translateLsvBins(d);

            var el = d3.select(this)
                .attr('transform', 'translate(' + (i * lsvp.histo_width + lsvp.left_padding) + ',' + (lsvp.height + lsvp.top_padding) + ') rotate(-90)')
                .attr('fill', 'transparent')
                .attr('stroke', colors(i));

            el
                .append('g')
                .attr('class', 'histograms')
                .violinHistograms(d, lsvp.height, lsvp.histo_width);


            el
                .append('g')
                .attr('class', 'title')
                .attr('transform', 'rotate(90) translate(' + lsvp.histo_width / 2 + ', 20)')
                .append('text')
                .attr('text-anchor', 'middle')
                .text(means_rounded[i])
        });
    return this;
};


d3.selection.prototype.violinHistograms = function (data, height, width) {

    var bins = d3.histogram()
        .domain([0, 1])
        .thresholds(d3.ticks(0, 1, 20))
        (data);

    var y = d3.scaleLinear()
        .domain(d3.extent(bins, function (d) {
            return d.length
        }))
        .range([width / 2, 0]);

    var x = d3.scaleLinear()
        .domain(d3.extent(bins, function (d) {
            return d.x0
        }))
        .range([0, height]);

    var area = d3.area()
        .curve(d3.curveBasisOpen)
        .x(function (d) {
            return x(d.x0)
        })
        .y1(function (d) {
            return y(d.length)
        })
        .y0(width / 2);

    // left violin
    this.append('path')
        .datum(bins)
        .classed('left', true)
        .attr('d', area);

    // right violin
    // this.append('path')
    //     .datum(bins)
    //     .classed('right', true)
    //     .attr('transform', 'translate(0,' + group_width + ') scale(1,-1)')
    //     .attr('d', area);

    return this;
};


d3.selection.prototype.voilinBoxPlots = function (data, height) {
    var width = 10;


    var q = d3.scaleQuantile()
        .domain([0, 100])
        .range(data);

    var bins = d3.histogram()
        .domain([0, 1])
        .thresholds(d3.ticks(0, 1, 20))
        (data);

    var y = d3.scaleLinear()
        .domain(d3.extent(bins, function (d) {
            return d.length
        }))
        .range([.5, height]);

    var x = d3.scaleLinear()
        .range([.5, height]);

    // box plot center horizontal line
    this.append('line')
        .attr('stroke', 'black')
        .classed('center-line', true)
        .attr('x1', x(q(5)))
        .attr('x2', x(q(95)))
        .attr('y1', y(0))
        .attr('y2', y(0));

    // box plot box (25% -> 75%)
    this.append('rect')
        .classed('box', true)
        .attr('width', x(q(75)) - x(q(25)))
        .attr('height', width + .5)
        .attr('x', x(q(25)))
        .attr('y', function () {
            return y(0) - (width / 2) - .5
        });

    // 5% line
    this.append('line')
        .attr('stroke', 'black')
        .classed('5-percentile', true)
        .attr('x1', x(q(5)))
        .attr('x2', x(q(5)))
        .attr('y1', y(0) - (width / 2))
        .attr('y2', y(0) + (width / 2));

    // 95% line
    this.append('line')
        .attr('stroke', 'black')
        .classed('95-percentile', true)
        .attr('x1', x(q(95)))
        .attr('x2', x(q(95)))
        .attr('y1', y(0) - (width / 2))
        .attr('y2', y(0) + (width / 2));

    // median line
    this.append('line')
        .classed('median', true)
        .attr('stroke', 'black')
        .attr('stroke-opacity', 1)
        .attr('x1', x(q(50)))
        .attr('x2', x(q(50)))
        .attr('y1', y(0) - (width / 2) - 4)
        .attr('y2', y(0) + (width / 2) + 4);

    // mean circle
    this.append('circle')
        .classed('mean', true)
        .attr('fill', 'white')
        .attr('stroke', 'None')
        .attr('cx', x(d3.mean(data)))
        .attr('cy', y(0))
        .attr('r', 2);

    return this;
};


LSVPlots.prototype.swarm = function () {
    var histo_width = this.histo_width;
    var left_padding = this.left_padding;
    var circle_radius = 2;
    var width = this.bins.length * this.histo_width;
    var svg = this.svg;
    var colors = new Colors().brewer_het(this.junc_idx);
    var x = d3.scaleLinear()
        .domain([0, 1])
        .range([this.height, this.top_padding]);

    var swarm_fn = d3.beeswarm()
        .distributeOn(function (d) {
            return x(d);
        })
        .radius(circle_radius)
        .orientation('vertical')
        .side('symetric');

    svg
        .selectAll('.swarm-group')
        .data(this.mu_psi)
        .enter()
        .append('g')
        .attr('class', 'swarm-group')
        .attr('transform', function (d, i) {
            return 'translate(' + ((histo_width * i) - left_padding) + ')'
        })
        .selectAll('circle')
        .data(function (d) {
            return swarm_fn
                .data(d)
                .arrange();
        })
        .enter()
        .append("circle")
        .attr('fill', function (d, i) {
            return colors(0)
        })
        .attr("cx", function (bee) {
            return bee.x + (width / 2) + 2;
        })
        .attr("cy", function (bee) {
            return bee.y;
        })
        .attr("r", circle_radius);
};