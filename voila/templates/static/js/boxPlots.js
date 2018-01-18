var BoxPlots = function (db) {
    this.db = db;
    this.height = 125;
    this.histo_width = 80;
    this.left_padding = 41;
    this.top_padding = 10;
    this.bottom_padding = 25;
    this.svg_height = this.height + this.top_padding + this.bottom_padding
};

BoxPlots.prototype.psi = function (el, lsv_id) {
    var bp = this;
    this.db.get(lsv_id).then(function (data) {
        bp.svg = d3.select(el)
            .attr('width', (data.bins.length * bp.histo_width) + bp.left_padding)
            .attr('height', bp.svg_height);
        bp.violins(data);
        bp.drawYAxis();
    });
};

translateLsvBins = function (lsvBins) {
    var numSamples = 30;
    var tmpBins = [];
    var binsSize = lsvBins.length;
    var numCopies;
    lsvBins.forEach(function (b, i) {
        numCopies = Math.round(numSamples * b);
        tmpBins = tmpBins.concat(new Array(numCopies).fill((1 / binsSize) / 2 + (i / binsSize)))
    });
    return tmpBins;
};


BoxPlots.prototype.drawYAxis = function () {
    var bp = this;
    var yScale = d3.scaleLinear()
        .domain([0, 1])
        .range([bp.height, 0]);

    return this.svg
        .append('g')
        .classed('y-axis', true)
        .attr('transform', 'translate(' + bp.left_padding + ',' + bp.top_padding + ')')
        .call(d3.axisLeft(yScale).tickValues([0, .5, 1]))
        .append('text')
        .classed('y-axis-label', true)
        .text('\u03A8')
        .attr('text-anchor', 'middle')
        .attr('font-size', 12)
        .attr('font-family', 'sans-serif')
        .attr('fill', 'black')
        .attr('transform', 'rotate(-90) translate(-' + bp.height / 2 + ', -30)');
};

BoxPlots.prototype.violins = function (data) {
    var colors = new Colors().toRGBArray();
    var bp = this;

    return this.svg
        .selectAll('.violin')
        .data(data.bins)
        .enter()
        .append('g')
        .attr('class', 'violin')
        .each(function (d, i) {
            var el = d3.select(this)
                .attr('transform', 'translate(' + (i * bp.histo_width + bp.left_padding) + ',' + (bp.height + bp.top_padding) + ') rotate(-90)')
                .attr('fill', colors[i]);

            el
                .append('g')
                .attr('class', 'histograms')
                .violinHistograms(d, bp.height, bp.histo_width);

            el
                .append('g')
                .attr('class', 'box-plot')
                .attr('stroke', 'black')
                .attr('fill', 'black')
                .attr('transform', 'translate(0,' + bp.histo_width / 2 + ')')
                .voilinBoxPlots(d, bp.height);

            el
                .append('g')
                .attr('class', 'title')
                .attr('transform', 'rotate(90) translate(' + bp.histo_width / 2 + ', 20)')
                .append('text')
                .attr('text-anchor', 'middle')
                .text(data.means_rounded[i].toFixed(3))
        });
};


d3.selection.prototype.violinHistograms = function (data, height, width) {
    var bins = d3.histogram()
        .domain([0, 1])
        .thresholds(d3.ticks(0, 1, 20))
        (translateLsvBins(data));

    var y = d3.scaleLinear()
        .domain(d3.extent(bins, function (d) {
            return d.length
        }))
        .rangeRound([width / 2, 0]);

    var x = d3.scaleLinear()
        .domain(d3.extent(bins, function (d) {
            return d.x0
        }))
        .rangeRound([0, height]);

    var area = d3.area()
        .curve(d3.curveBasis)
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
    this.append('path')
        .datum(bins)
        .classed('right', true)
        .attr('transform', 'translate(0,' + width + ') scale(1,-1)')
        .attr('d', area);

    return this;
};

d3.selection.prototype.voilinBoxPlots = function (data, height) {
    data = translateLsvBins(data);

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

