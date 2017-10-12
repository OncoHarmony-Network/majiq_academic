var Het = function (opts) {
    this.svg_height = 500;
    this.ajax_opts = {url: opts.url, data: opts.query};
    this.svg = d3v4.select(opts.element)
        .attr('width', '100%')
        .attr('height', this.svg_height);
    this.colors = new Colors().toArray();
    this.svg.selectAll('*').remove();
    var het = this;

    $.ajax(this.ajax_opts).done(function (json) {
        d3v4.select(opts.element)

            .selectAll('.group')
            .data(json.data)
            .enter().append('g')
            .attr('class', 'group')
            .attr('transform', function (d, i) {
                return 'translate(' + 500 * i + ')'
            })

            .selectAll('.junction')
            .data(function (d) {
                return d
            })
            .enter().append('g')
            .attr('class', 'junction')
            .attr('transform', function (d, i) {
                return 'translate(' + 50 * i + ')'
            })
            .attr('fill', function (d, i) {
                return het.colors[i]
            })

            .selectAll('.bin')
            .data(function (d) {
                return d
            })
            .enter().append('g')
            .attr('class', function (d) {
                return 'bin ' + d.bin
            })
            .attr('transform', function (d) {
                return 'translate(0, ' + ((1 - d.bin) * (het.svg_height - 3)) + ')'
            })
            .on('mouseover', function (d) {
                d3v4.select(this)
                    .append('text')
                    .attr('class', 'overlay')
                    .attr('fill', 'black')
                    .attr('x', 25)
                    .attr('y', -10)
                    .attr('font-size', 12)
                    .attr('font-family', 'sans-serif')
                    .text(d.values.map(function (x) {
                        return x.experimentName + ': ' + x.expectedPsi
                    }).join(' '));
            })
            .on('mouseout', function () {
                d3v4.select(this).select('.overlay').remove()
            })

            .selectAll('.expected-psi')
            .data(function (d) {
                return d.values
            })
            .enter().append('circle')
            .attr("cx", function (d, i) {
                return 3 + (i * 5)
            })
            .attr("cy", 0)
            .attr("r", 2)
    })
};

Het.prototype.init = function () {
    var het = this;
    $.ajax(this.ajax_opts).done(function (json) {
        het.svg.selectAll('.group')
            .data(json.data)
            .enter().append('g')
            .attr('class', 'group')
            .attr('transform', function (d, i) {
                return 'translate(' + 200 * i + ')'
            })
            .selectAll('.sample')
            .data(function (d) {
                return d.samples
            })
            .enter().append('g')
            .attr('class', 'sample')
            .selectAll('.expected-psi')
            .data(function (d) {
                return d.expected_psi
            })
            .enter().append('circle')
            .attr("r", 2)
            .attr('class', 'expected-psi')
            .style('fill', function (d, i) {
                return het.colors[i]
            })
            .attr("cx", function (d, i) {
                return 25 + (25 * i )
            })
            .attr("cy", function (d) {
                return het.svg_height * (1 - d)
            })
    })
};

Het.prototype.update = function () {
    var het = this;

    $.ajax(this.ajax_opts).done(function (json) {

        var groups = het.svg.selectAll('.group')
            .data(json.data);

        groups.exit().remove();

        groups.enter().append('g')
            .attr('class', 'group')
            .attr('transform', function (d, i) {
                return 'translate(' + 200 * i + ')'
            });

        var sample = groups.selectAll('.sample')
            .data(function (d) {
                return d.samples
            });

        sample.exit().remove();

        sample.enter().append('g')
            .attr('class', 'sample');

        var epsi = sample.selectAll('.expected-psi')
            .data(function (d) {
                return d.expected_psi
            });

        epsi.exit().transition().remove();

        epsi.enter().append('circle')
            .attr("r", 2)
            .attr('class', 'expected-psi')
            .merge(epsi)
            .style('fill', function (d, i) {
                return het.colors[i]
            })
            .transition()
            .attr("cx", function (d, i) {
                return 25 + (25 * i)
            })
            .transition()
            .attr("cy", function (d) {
                return het.svg_height * (1 - d)
            })
    })
};