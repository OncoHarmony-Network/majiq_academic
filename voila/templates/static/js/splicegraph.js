var SpliceGraph = function (db) {
    this.db = db;
    this.width = 1000;
    this.height = 300;

    this.x = {};
    this.y = null;

    this.exon_height = 20;
    this.font_size = 12;

    this.zoom = 1;

    this.max_height = this.height - 5;
};

SpliceGraph.prototype.yScale = function () {
    if (this.y)
        return this.y;
    this.y = d3.scaleLinear().domain([0, this.max_height]).range([this.max_height, 0]);
    return this.yScale();
};

SpliceGraph.prototype.xScale = function (gene, default_view, reverse_range, experiment) {
    var x_dom = [];
    var x_range = [];
    var min_width = 10;
    var max_width = (this.width * this.zoom) - 10;
    var i;
    var j;
    var length;
    var max;
    var min;
    var offset;

    // if we're not using the default view, the x-scale if very simple.
    if (!default_view) {
        x_range = [min_width, max_width];
        if (reverse_range)
            x_range.reverse();
        return d3.scaleLinear().domain([gene.start, gene.end]).range(x_range);
    }

    // general x-scale
    var x = d3.scaleLinear().domain([gene.start, gene.end]).range([min_width, max_width]);

    // get the start and end of each exon/ir for both the domain and range
    gene.exons.forEach(function (exon) {
        if (!exon.intron_retention) {
            x_dom.push(exon.start);
            x_dom.push(exon.end);
            x_range.push(x(exon.start));
            x_range.push(x(exon.end));
        }
    });

    if (reverse_range) {
        x_range.reverse()
    }

    // adjust exon sizes
    gene.exons.filter(function (d) {
        return !d.intron_retention
    }).forEach(function (e, i) {
        var start = x_range[i * 2];
        var end_idx = i * 2 + 1;
        var end = x_range[end_idx];
        length = end - start;
        offset = 0;


        if ([4, 5].includes(gene.exon_types[e.start][e.end][experiment])) {
            min = 4;
            max = 4;
        } else {
            min = 5;
            max = 100;

        }

        if (length < min)
            offset = min - length;

        if (length > max)
            offset = max - length;

        if (offset !== 0)
            for (j = end_idx; j < x_range.length; j++)
                x_range[j] += offset;
    });

    // adjust spaces between exons
    var ir_count = 0;
    gene.exons.forEach(function (e, i) {
        if (e.intron_retention) ir_count++;

        i -= ir_count;

        length = x_range[i * 2 + 2] - x_range[i * 2 + 1];
        offset = 0;


        if (e.intron_retention) {
            min = 75;
            if (length < min)
                offset = min - length;
        } else {
            max = 30;
            min = 10;
            if (length > max)
                offset = max - length;

            if (length < min)
                offset = min - length;
        }


        if (offset !== 0)
            for (j = (i * 2) + 2; j < x_range.length; j++)
                x_range[j] = x_range[j] + offset
    });


    if (reverse_range) {
        x_range.reverse()
    }

    // scale back to view width
    x = d3.scaleLinear().domain([x_range[0], x_range[x_range.length - 1]]).range([min_width, max_width]);
    x_range = x_range.reduce(function (accu, curr) {
        accu.push(x(curr));
        return accu
    }, []);

    if (reverse_range) {
        x_range.reverse();
    }
    return d3.scaleLinear().domain(x_dom).range(x_range);
};

SpliceGraph.prototype.find_min_max = function (j, x) {
    var center = x(j.start + ((j.end - j.start) / 2));
    var min = center - 20;
    var max = center + 20;
    return [min, max]
};

SpliceGraph.prototype.juncs_overlap = function (a, b, x) {
    var a_min_max = this.find_min_max(a, x);
    var a_min = a_min_max[0];
    var a_max = a_min_max[1];

    var b_min_max = this.find_min_max(b, x);
    var b_min = b_min_max[0];
    var b_max = b_min_max[1];

    return (a_min <= b_min && b_min <= a_max) || (a_min <= b_max && b_max <= a_max);
};

SpliceGraph.prototype.overlaps = function (grp, j_idx, juncs, x) {
    var j = juncs[j_idx];
    var member;
    var rtn = false;
    grp.forEach(function (member_idx) {
        member = juncs[member_idx];
        if (juncs_overlap(member, j, x)) {
            rtn = true;
        }
    });
    return rtn
};

SpliceGraph.prototype.junctions_no_ir = function (gene, x) {
    var sg = this;
    var leader;
    var unavailable = [];
    var grps = {};
    var juncs_no_ir = gene.junctions
        .filter(function (j) {
            return !j.intron_retention
        });

    juncs_no_ir.forEach(function (a, a_idx) {
        if (!unavailable.includes(a_idx)) {
            leader = a_idx;
            grps[leader] = [leader];
            unavailable.push(leader);
            juncs_no_ir.forEach(function (b, b_idx) {
                if (!unavailable.includes(b_idx)) {
                    if (sg.juncs_overlap(juncs_no_ir[leader], b, x)) {
                        grps[leader].push(b_idx);
                        unavailable.push(b_idx)
                    }
                }
            })
        }
    });

    for (const key in grps) {
        var grp = grps[key];

        grp.sort(function (a, b) {
            var a_junc = juncs_no_ir[a];
            var b_junc = juncs_no_ir[b];
            return (a_junc.end - a_junc.start) - (b_junc.end - b_junc.start)
        });

        grp.forEach(function (j_idx, i) {
            juncs_no_ir[j_idx].bin = i + 1
        })
    }

    return juncs_no_ir
};

SpliceGraph.prototype.init = function (sg_div, experiment) {
    var gene_id = sg_div.getAttribute('data-gene-id');
    var sg = this;
    sg_div.setAttribute('data-zoom', 1.0);
    this.db.get(gene_id).then(function (gene) {

        sg.sg_div = d3.select(sg_div);
        sg.gene = gene;

        var svg = sg.sg_div.append('svg')
            .attr('width', sg.width)
            .attr('height', sg.height);

        var exon_height = sg.exon_height;
        var font_size = sg.font_size;
        var default_view = sg.sg_div.classed('default-view');
        var reversed_range = gene.strand === '-';

        var y = sg.yScale();
        var x = sg.xScale(gene, default_view, reversed_range, experiment);

        var juncs_no_ir = sg.junctions_no_ir(gene, x);
        var exons = gene.exons.filter(function (d) {
            return !d.intron_retention && ![4, 5].includes(gene.exon_types[d.start][d.end][experiment])
        });

        svg.selectAll('.intron-retention')
            .data(gene.exons.filter(function (e) {
                return e.intron_retention
            }))
            .enter()
            .append('polygon')
            .intron_retention(x, y, exon_height, gene.exon_types, experiment);


        svg.selectAll('.exon')
            .data(exons)
            .enter()
            .append('polygon')
            .exons(x, y, exon_height, gene.exon_types, experiment);

        svg.selectAll('.exon-number')
            .data(exons)
            .enter()
            .append('text')
            .exon_numbers(x, y, exon_height, font_size, gene.strand);

        svg.selectAll('.junction-grp')
            .data(juncs_no_ir)
            .enter()
            .append('g')
            .attr('class', 'junction-grp')
            .each(function (d, i) {
                d3.select(this)
                    .selectAll('.junction')
                    .data([d])
                    .enter()
                    .append('path')
                    .junctions(x, y, exon_height, gene.strand, gene.reads, gene.junction_types, experiment);

                d3.select(this)
                    .selectAll('.reads')
                    .data([d])
                    .enter()
                    .append('text')
                    .reads(x, y, exon_height, font_size, gene.reads, experiment);

                d3.select(this)
                    .selectAll('.ss3p')
                    .data([d])
                    .enter()
                    .append('line')
                    .ss3p(x, y, exon_height, gene.reads, gene.junction_types, experiment);

                d3.select(this)
                    .selectAll('.ss5p')
                    .data([d])
                    .enter()
                    .append('line')
                    .ss5p(x, y, exon_height, gene.reads, gene.junction_types, experiment)
            });

        svg.selectAll('.ir-line')
            .data(gene.junctions.filter(function (j) {
                return j.intron_retention
            }))
            .enter()
            .append('polyline')
            .ir_lines(x, y, gene.strand);

        svg.selectAll('.half-exon')
            .data(gene.exons.filter(function (d) {
                return [4, 5].includes(gene.exon_types[d.start][d.end][experiment])
            }))
            .enter()
            .append('polyline')
            .half_exons(x, y, exon_height, gene.exon_types, experiment);

        svg.selectAll('.ir-reads')
            .data(gene.junctions.filter(function (j) {
                return j.intron_retention
            }))
            .enter()
            .append('text')
            .ir_reads(x, y, exon_height, font_size, gene.strand, gene.reads, experiment);

    })
};

d3.transition.prototype.ss3p =
    d3.selection.prototype.ss3p =
        function (x, y, exon_height, reads, junction_types, experiment) {
            return this
                .attr('class', 'ss3p')
                .attr('x1', function (d) {
                    return x(d.start)
                })
                .attr('x2', function (d) {
                    return x(d.start)
                })
                .splice_site(y, exon_height, reads, junction_types, experiment)
        };

d3.transition.prototype.ss5p =
    d3.selection.prototype.ss5p =
        function (x, y, exon_height, reads, junction_types, experiment) {
            return this
                .attr('class', 'ss5p')
                .attr('x1', function (d) {
                    return x(d.end)
                })
                .attr('x2', function (d) {
                    return x(d.end)
                })
                .splice_site(y, exon_height, reads, junction_types, experiment)
        };

d3.transition.prototype.splice_site =
    d3.selection.prototype.splice_site =
        function (y, exon_height, reads, junction_types, experiment) {
            return this
                .attr('y1', y(0))
                .attr('y2', y(exon_height))
                .attr('stroke', 'black')
                .style_junctions(reads, junction_types, experiment)
                .attr('stroke-dasharray', '2,2')

        };

d3.transition.prototype.ir_reads =
    d3.selection.prototype.ir_reads =
        function (x, y, exon_height, font_size, strand, reads, experiment) {
            return this
                .attr('class', 'ir-reads')
                .text(function (d) {
                    var r = reads[d.start][d.end][experiment];
                    if (r)
                        return r
                })
                .attr('x', function (d) {
                    if (strand === '-') {
                        if (d.intron_retention === 2)
                            return x(d.start) + 7;
                        if (d.intron_retention === 1)
                            return x(d.start) - 7
                    } else {
                        if (d.intron_retention === 2)
                            return x(d.start) - 7;
                        if (d.intron_retention === 1)
                            return x(d.start) + 7
                    }
                })
                .attr('y', y(exon_height - 3))
                .attr('text-anchor', function (d) {
                    if (strand === '-') {
                        if (d.intron_retention === 2)
                            return 'start';
                        else
                            return 'end'
                    } else {
                        if (d.intron_retention === 2)
                            return 'end';
                        else
                            return 'start'
                    }
                })
                .attr('font-family', 'sans-serif')
                .attr('font-size', font_size);
        };

d3.transition.prototype.half_exons =
    d3.selection.prototype.half_exons =
        function (x, y, exon_height, exon_types, experiment) {
            return this
                .attr('class', 'half-exon')
                .style_exons(exon_types, experiment)
                .attr('points', function (d) {
                    var exon_type = exon_types[d.start][d.end][experiment];

                    if (exon_type === 4)
                        return [
                            [x(d.start), y(0)].join(' '),
                            [x(d.end), y(0)].join(' '),
                            [x(d.end), y(exon_height)].join(' '),
                            [x(d.start), y(exon_height)].join(' ')
                        ].join(', ');

                    if (exon_type === 5)
                        return [
                            [x(d.end), y(0)].join(' '),
                            [x(d.start), y(0)].join(' '),
                            [x(d.start), y(exon_height)].join(' '),
                            [x(d.end), y(exon_height)].join(' ')
                        ].join(', ');
                });
        };

d3.transition.prototype.style_exons =
    d3.selection.prototype.style_exons =
        function (exon_types, experiment) {
            return this
                .attr('fill-opacity', .3)
                .attr('stroke-linejoin', 'round')
                .attr('fill', function (d) {
                    switch (exon_types[d.start][d.end][experiment]) {
                        case 0:
                            return 'grey';
                        case 1:
                        case 4:
                        case 5:
                            return 'green';
                        default:
                            return 'transparent'
                    }
                })
                .attr('stroke', function (d) {
                    switch (exon_types[d.start][d.end][experiment]) {
                        case 1:
                        case 4:
                        case 5:
                            return 'green';
                        default:
                            return 'black'
                    }
                })
                .attr('stroke-dasharray', function (d) {
                    if (exon_types[d.start][d.end][experiment] === 2) {
                        return '5, 2'
                    }
                })
        };

d3.transition.prototype.exons =
    d3.selection.prototype.exons =
        function (x, y, exon_height, exon_types, experiment) {
            return this
                .attr('class', 'exon')
                .style_exons(exon_types, experiment)
                .attr('points', function (d) {
                    return [
                        [x(d.start), y(0)].join(' '),
                        [x(d.end), y(0)].join(' '),
                        [x(d.end), y(exon_height)].join(' '),
                        [x(d.start), y(exon_height)].join(' ')
                    ].join(', ')
                });
        };


d3.transition.prototype.intron_retention =
    d3.selection.prototype.intron_retention =
        function (x, y, exon_height, exon_types, experiment) {
            return this
                .attr('class', 'intron-retention')
                .style_exons(exon_types, experiment)
                .attr('points', function (d) {
                    return [
                        [x(d.start - 1), y(exon_height / 4)].join(' '),
                        [x(d.end + 1), y(exon_height / 4)].join(' '),
                        [x(d.end + 1), y(exon_height * (3 / 4))].join(' '),
                        [x(d.start - 1), y(exon_height * (3 / 4))].join(' ')
                    ].join(', ')
                })
        };

d3.transition.prototype.exon_numbers =
    d3.selection.prototype.exon_numbers =
        function (x, y, exon_height, font_size, strand) {
            var size = this.size();
            return this
                .attr('class', 'exon-number')
                .text(function (d, i) {
                    if (strand === '+')
                        return i + 1;
                    else
                        return size - i
                })
                .attr('y', y((exon_height / 2) - (font_size / 2) + 2))
                .attr('x', function (d) {
                    return x(d.start + (d.end - d.start) / 2)
                })
                .attr('text-anchor', 'middle')
                .attr('font-family', 'sans-serif')
                .attr('font-size', font_size);
        };

d3.transition.prototype.style_junctions =
    d3.selection.prototype.style_junctions =
        function (reads, junction_types, experiment) {
            return this
                .attr('stroke-width', 1.5)
                .attr('stroke-dasharray', function (d) {
                    switch (junction_types[d.start][d.end][experiment]) {
                        case 3:
                        case 2:
                            return '5,2';
                        case 1:
                            if (reads[d.start][d.end][experiment] === 0)
                                return '5,2';
                    }
                })
                .attr('stroke', function (d) {
                    switch (junction_types[d.start][d.end][experiment]) {
                        case 3:
                        case 0:
                            return 'red';
                        case 1:
                            return 'green';
                        case 2:
                            return 'grey';
                        default:
                            return 'black'
                    }
                })
                .attr('fill', 'None')
        };

d3.transition.prototype.junctions =
    d3.selection.prototype.junctions =
        function (x, y, exon_height, strand, reads, junction_types, experiment) {
            var junc_height = 20;
            return this
                .attr('class', 'junction')
                .style_junctions(reads, junction_types, experiment)
                .attr('d', function (d) {
                    var sweep_flag = strand === '+' ? 1 : 0;
                    var junc_length = x(d.end) - x(d.start);
                    // where junctions are very long... put them one bin higher.
                    var long_junc = Math.abs(junc_length) > 200 ? 1 : 0;
                    return 'M' + [x(d.start), y(exon_height)].join(',') +
                        'A' + [junc_length / 2, junc_height * (d.bin + long_junc), 0, 0, sweep_flag, x(d.end), y(exon_height)].join(' ')
                })
        };


d3.transition.prototype.reads =
    d3.selection.prototype.reads =
        function (x, y, exon_height, font_size, reads, experiment) {
            var junc_height = 20;
            return this
                .attr('class', 'reads')
                .text(function (d) {
                    try {
                        var r = reads[d.start][d.end][experiment];
                        if (r)
                            return r
                    } catch (TypeError) {
                        return '';
                    }
                })
                .attr('x', function (d) {
                    return x(d.start) + (x(d.end) - x(d.start)) / 2
                })
                .attr('y', function (d) {
                    // where junctions are very long... put them one bin higher.
                    var junc_length = x(d.end) - x(d.start);
                    var long_junc = Math.abs(junc_length) > 200 ? 1 : 0;
                    return y(exon_height + (junc_height * (d.bin + long_junc)) + 3)
                })
                .attr('text-anchor', 'middle')
                .attr('font-family', 'sans-serif')
                .attr('font-size', font_size);
        };


d3.transition.prototype.ir_lines =
    d3.selection.prototype.ir_lines =
        function (x, y, strand) {
            return this
                .attr('class', 'ir-line')
                .attr('stroke', 'black')
                .attr('points', function (d) {
                    if (d.intron_retention === 1)
                        if (strand === '+')
                            return [
                                [x(d.start), y(18)].join(' '),
                                [x(d.start) + 5, y(18)].join(' ')
                            ].join(', ');
                        else if (strand === '-')
                            return [
                                [x(d.start), y(18)].join(' '),
                                [x(d.start) - 5, y(18)].join(' ')
                            ].join(', ');
                    if (d.intron_retention === 2)
                        if (strand === '+')
                            return [
                                [x(d.end), y(18)].join(' '),
                                [x(d.end) - 5, y(18)].join(' ')
                            ].join(', ');
                        else if (strand === '-')
                            return [
                                [x(d.end), y(18)].join(' '),
                                [x(d.end) + 5, y(18)].join(' ')
                            ].join(', ');
                })
        };


SpliceGraph.prototype.update = function (sg_div, experiment) {
    // var gene_id = sg_div.getAttribute('data-gene-id');
    var sg = this;
    this.zoom = sg_div.getAttribute('data-zoom');

    this.db.get(sg_div.getAttribute('data-gene-id')).then(function (gene) {

        var default_view = sg_div.classList.contains('default-view');
        var svg = d3.select(sg_div).select('svg');
        var y = sg.yScale();
        var exon_height = sg.exon_height;
        var reverse_range = gene.strand === '-';
        var x = sg.xScale(gene, default_view, reverse_range);
        var font_size = sg.font_size;
        var juncs_no_ir = sg.junctions_no_ir(gene, x);

        svg
            .interrupt()
            .transition()
            .attr('width', sg.width * sg.zoom);


        svg.selectAll('.exon')
            .interrupt()
            .transition()
            .exons(x, y, exon_height, gene.exon_types, experiment);

        svg.selectAll('.half-exon')
            .interrupt()
            .transition()
            .half_exons(x, y, exon_height, gene.exon_types, experiment);

        svg.selectAll('.intron-retention')
            .interrupt()
            .transition()
            .intron_retention(x, y, exon_height, gene.exon_types, experiment);

        svg.selectAll('.exon-number')
            .interrupt()
            .transition()
            .exon_numbers(x, y, exon_height, font_size, gene.strand);

        svg.selectAll('.junction')
            .interrupt()
            .data(juncs_no_ir)
            .transition()
            .junctions(x, y, exon_height, gene.strand, gene.reads, gene.junction_types, experiment);


        svg.selectAll('.reads')
            .interrupt()
            .data(juncs_no_ir)
            .transition()
            .reads(x, y, exon_height, font_size, gene.reads, experiment);

        svg.selectAll('.ir-line')
            .interrupt()
            .transition()
            .ir_lines(x, y, gene.strand);

        svg.selectAll('.ir-reads')
            .interrupt()
            .data(gene.junctions.filter(function (d) {
                return d.intron_retention
            }))
            .transition()
            .ir_reads(x, y, exon_height, font_size, gene.strand, gene.reads, experiment);

        svg.selectAll('.ss3p')
            .interrupt()
            .data(juncs_no_ir)
            .transition(sg.t)
            .ss3p(x, y, exon_height, gene.reads, gene.junction_types, experiment);

        svg.selectAll('.ss5p')
            .interrupt()
            .data(juncs_no_ir)
            .transition()
            .ss5p(x, y, exon_height, gene.reads, gene.junction_types, experiment);
    })

};