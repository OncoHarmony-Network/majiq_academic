var Violin = function (db) {
    // height and width before the svg transformations
    this.width = 110;
    this.height = 50;
    this.paddingLeft = 40;
    this.paddingTop = 15;
    this.boxPlotHeight = 10;
    this.db = db;

    // this.data = opts.data;
    // this.color = opts.color;
    // this.sampleName = opts.sampleName;

    // this.bins = d3.histogram()
    //     .domain([0, 1])
    //     .thresholds(d3.ticks(0, 1, 20))
    //     (this.data);
    //
    // this.bx = d3.scaleLinear()
    //     .rangeRound([0, this.width]);
    //
    // this.hx = d3.scaleLinear()
    //     .domain(d3.extent(this.bins, function (d) {
    //         return d.x0
    //     }))
    //     .rangeRound([0, this.width]);
    //
    // this.y = d3.scaleLinear()
    //     .domain(d3.extent(this.bins, function (d) {
    //         return d.length
    //     }))
    //     .rangeRound([this.height, 0]);
    //
    // this.yScale = d3.scaleLinear()
    //     .domain([0, 1])
    //     .range([this.width, 0]);
    //
    // this.q = d3.scaleQuantile()
    //     .domain([0, 100])
    //     .range(this.data);
    //
    // this.area = d3.area()
    //     .curve(d3.curveBasis)
    //     .x(function (d) {
    //         return violin.hx(d.x0)
    //     })
    //     .y1(function (d) {
    //         return violin.y(d.length)
    //     })
    //     .y0(this.height);

};

Violin.prototype.translateLsvBins = function (lsvBins) {
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


Violin.prototype.psi = function (el, lsv_id) {
    var violin = this;
    this.element = el;
    this.db.get(lsv_id).then(function (data) {
        console.log(data.bins);
        this.data = violin.translateLsvBins(data.bins);

        this.bins = d3.histogram()
            .domain([0, 1])
            .thresholds(d3.ticks(0, 1, 20))
            (this.data);

        this.bx = d3.scaleLinear()
            .rangeRound([0, this.width]);

        this.hx = d3.scaleLinear()
            .domain(d3.extent(this.bins, function (d) {
                return d.x0
            }))
            .rangeRound([0, this.width]);

        this.y = d3.scaleLinear()
            .domain(d3.extent(this.bins, function (d) {
                return d.length
            }))
            .rangeRound([this.height, 0]);

        this.yScale = d3.scaleLinear()
            .domain([0, 1])
            .range([this.width, 0]);

        this.q = d3.scaleQuantile()
            .domain([0, 100])
            .range(this.data);

        this.area = d3.area()
            .curve(d3.curveBasis)
            .x(function (d) {
                return violin.hx(d.x0)
            })
            .y1(function (d) {
                return violin.y(d.length)
            })
            .y0(this.height);

        this.draw()
    });


};

Violin.prototype.draw = function () {
    this.drawYAxis();
    this.drawPsi();
    var violinPlot = this.drawViolins();
    this.drawBoxPlot(violinPlot);
    this.drawTitle();
};


Violin.prototype.drawTitle = function () {
    this.element
        .append('g')
        .classed('title', true)
        .append('text')
        .text(this.sampleName)
        .attr('text-anchor', 'middle')
        .attr('font-size', 12)
        .attr('font-family', 'sans-serif')
        .attr('fill', 'black')
        .attr('transform', 'translate(' + (this.paddingLeft + this.height) + ',10)')
};

Violin.prototype.drawYAxis = function () {
    var yAxis = this.element
        .append('g')
        .classed('y-axis', true)
        .attr('transform', 'translate(' + this.paddingLeft + ',' + this.paddingTop + ')')
        .call(d3.axisLeft(this.yScale).tickValues([0, .5, 1]))
        .append('text')
        .classed('y-axis-label', true)
        .text('Ψ')
        .attr('text-anchor', 'middle')
        .attr('font-size', 12)
        .attr('font-family', 'sans-serif')
        .attr('fill', 'black')
        .attr('transform', 'rotate(-90) translate(-55, -30)')
};

Violin.prototype.drawPsi = function () {
    this.element
        .append('g')
        .classed('psi-value', true)
        .append('text')
        .attr('transform', 'translate(' + (this.paddingLeft + this.height) + ',' + (this.width + this.paddingTop + 25) + ')')
        .attr('text-anchor', 'middle')
        .attr('font-size', 12)
        .attr('font-family', 'sans-serif')
        .text(d3.format(".3r")(d3.mean(this.data)) + ' E(Ψ)');
};


Violin.prototype.drawViolins = function () {
    var violinPlot = this.element
        .append('g')
        .classed('violin-plot', true)
        .attr('transform', 'translate(' + (this.paddingLeft + 1) + ',' + (this.width + this.paddingTop) + ') rotate(-90)');

    var histograms = violinPlot.append('g')
        .classed('histograms', true)
        .attr('fill', this.color);

    // left violin
    histograms.append('path')
        .datum(this.bins)
        .classed('left', true)
        .attr('d', this.area);

    // right violin
    histograms.append('path')
        .datum(this.bins)
        .classed('right', true)
        .attr('transform', 'translate(0,' + this.height * 2 + ') scale(1,-1)')
        .attr('d', this.area);

    return violinPlot;
};

Violin.prototype.drawBoxPlot = function (violinPlot) {
    var boxPlot = violinPlot.append('g')
        .classed('box-plot', true)
        .attr('fill', 'black')
        .attr('stroke', 'black');

    // box plot box (25% -> 75%)
    boxPlot.append('rect')
        .classed('box', true)
        .attr('stroke', 'none')
        .attr('width', this.bx(this.q(75)) - this.bx(this.q(25)))
        .attr('height', this.boxPlotHeight)
        .attr('x', this.bx(this.q(25)))
        .attr('y', this.y(0) - (this.boxPlotHeight / 2));

    // box plot center horizontal line
    boxPlot.append('line')
        .classed('center-line', true)
        .attr('x1', this.bx(this.q(5)))
        .attr('x2', this.bx(this.q(95)))
        .attr('y1', this.y(0))
        .attr('y2', this.y(0));

    // 5% line
    boxPlot.append('line')
        .classed('5-percentile', true)
        .attr('x1', this.bx(this.q(5)))
        .attr('x2', this.bx(this.q(5)))
        .attr('y1', this.y(0) - (this.boxPlotHeight / 2))
        .attr('y2', this.y(0) + (this.boxPlotHeight / 2));

    // 95% line
    boxPlot.append('line')
        .classed('95-percentile', true)
        .attr('x1', this.bx(this.q(95)))
        .attr('x2', this.bx(this.q(95)))
        .attr('y1', this.y(0) - (this.boxPlotHeight / 2))
        .attr('y2', this.y(0) + (this.boxPlotHeight / 2));

    // median line
    boxPlot.append('line')
        .classed('median', true)
        .attr('stroke', 'black')
        .attr('stroke-opacity', 1)
        .attr('x1', this.bx(this.q(50)))
        .attr('x2', this.bx(this.q(50)))
        .attr('y1', this.y(0) - (this.boxPlotHeight / 2) - 2)
        .attr('y2', this.y(0) + (this.boxPlotHeight / 2) + 2);

    // mean circle
    boxPlot.append('circle')
        .classed('mean', true)
        .attr('fill', 'white')
        .attr('stroke', 'None')
        .attr('cx', this.bx(d3.mean(this.data)))
        .attr('cy', this.y(0))
        .attr('r', 2);
};
