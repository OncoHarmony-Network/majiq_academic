var BoxPlots = function (opts) {
    var svg_height = 150;

    opts.elements.forEach(function (el) {
        d3v4.select(el).selectAll('*').remove();
        d3v4.select(el)
            .attr('height', svg_height)
            // .attr('width', '100%')
            .attr('width', '2000')
            .violins(JSON.parse(d3v4.select(el).attr('box-plot-data')))
    });
    this.svg_height = svg_height;
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

d3v4.selection.prototype.violins = function (data) {
    var colors = new Colors().toRGBArray();
    var svgHeight = this.svg_height;
    var histogramWidth = 80;
    this
        .selectAll('.group')
        .data(data)
        .enter().append('g')
        .attr('class', 'group')
        .attr('transform', function (d, i) {
            return 'translate(' + i * d.length * histogramWidth + ')'
        })
        .selectAll('.het-violin')
        .data(function (d) {
            return d
        })
        .enter().append('g')
        .attr('class', 'het-violin')
        .each(function (d, i) {
            var el = d3v4.select(this)
                .attr('transform', 'translate(' + i * histogramWidth + ',' + svgHeight + ') rotate(-90)')
                .attr('fill', colors[i]);

            el
                .append('g')
                .attr('class', 'histograms')
                .violinHistograms(d, svgHeight, histogramWidth);

            el
                .append('g')
                .attr('class', 'box-plot')
                .attr('stroke', 'black')
                .attr('fill', 'black')
                .attr('transform', 'translate(0,' + histogramWidth / 2 + ')')
                .voilinBoxPlots(d, svgHeight);
        })
};


d3v4.selection.prototype.violinHistograms = function (data, height, width) {
    var bins = d3v4.histogram()
        .domain([0, 1])
        .thresholds(d3v4.ticks(0, 1, 20))
        (translateLsvBins(data));

    var y = d3v4.scaleLinear()
        .domain(d3v4.extent(bins, function (d) {
            return d.length
        }))
        .rangeRound([width / 2, 0]);

    var x = d3v4.scaleLinear()
        .domain(d3v4.extent(bins, function (d) {
            return d.x0
        }))
        .rangeRound([0, height]);

    var area = d3v4.area()
        .curve(d3v4.curveBasis)
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

d3v4.selection.prototype.voilinBoxPlots = function (data, height) {
    data = translateLsvBins(data);

    var width = 10;

    var q = d3v4.scaleQuantile()
        .domain([0, 100])
        .range(data);

    var bins = d3v4.histogram()
        .domain([0, 1])
        .thresholds(d3v4.ticks(0, 1, 20))
        (data);

    var y = d3v4.scaleLinear()
        .domain(d3v4.extent(bins, function (d) {
            return d.length
        }))
        .range([.5, height]);

    var x = d3v4.scaleLinear()
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
        .attr('cx', x(d3v4.mean(data)))
        .attr('cy', y(0))
        .attr('r', 2);

    return this;
};

