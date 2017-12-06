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

    this.t = d3.transition()
        .duration(750);
};

SpliceGraph.prototype.yScale = function () {
    if (this.y)
        return this.y;
    this.y = d3.scaleLinear().domain([0, this.max_height]).range([this.max_height, 0]);
    return this.yScale();
};

SpliceGraph.prototype.xScale = function (default_view, reverse_range) {
    // create domain and range for x scale
    var x_dom = [];
    var x_range = [];
    var gene = this.gene;
    var min_width = 10;
    var max_width = (this.width * this.zoom) - 10;

    var x = d3.scaleLinear().domain([gene.start, gene.end]).range([min_width, max_width]);
    this.gene.exons.forEach(function (exon) {
        if (!exon.intron_retention) {

            x_dom.push(exon.start);
            x_dom.push(exon.end);

            x_range.push(x(exon.start));
            x_range.push(x(exon.end));
        }
    });


    if (!default_view) {
        if (reverse_range)
            x_range.reverse();
        return d3.scaleLinear().domain(x_dom).range(x_range);
    }

    // adjust exon spacing
    x_range = x_range.reduce(function (accu, curr, index) {
        if (index && (accu.length + 1) % 2) {
            var prev = accu[accu.length - 1];
            var length = curr - prev;
            var offset = 0;
            var max = 100;
            var min = 10;

            if (length > max)
                offset += max - length;

            if (length < min)
                offset += min - length;

            if (offset !== 0) {
                curr += offset;
                for (var i = index + 1; i < x_range.length; i++) {
                    x_range[i] += offset
                }
            }
        }
        accu.push(curr);
        return accu
    }, []);


    // adjust exon width
    x_range = x_range.reduce(function (accu, curr, index) {
        if (index && accu.length % 2) {
            var prev = accu[accu.length - 1];
            var length = curr - prev;
            var offset = 0;
            var min = 20;
            var max = 200;

            if (length < min)
                offset = min - length;

            if (length > max)
                offset = max - length;

            if (offset !== 0) {
                curr += offset;
                for (var i = index + 1; i < x_range.length; i++) {
                    x_range[i] += offset
                }
            }
        }
        accu.push(curr);
        return accu
    }, []);

    // scale back to view width
    x = d3.scaleLinear().domain([x_range[0], x_range[x_range.length - 1]]).range([min_width, max_width]);
    x_range = x_range.reduce(function (accu, curr) {
        accu.push(x(curr));
        return accu
    }, []);

    if (reverse_range)
        x_range.reverse();
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
    var juncs_no_ir = gene.junctions.filter(function (j) {
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
    this.db.get(gene_id + '_' + experiment).then(function (gene) {

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
        var x = sg.xScale(default_view, reversed_range);

        var juncs_no_ir = sg.junctions_no_ir(gene, x);

        svg.selectAll('.exon')
            .data(function () {
                return gene.exons.filter(function (e) {
                    return !e.intron_retention
                })
            })
            .enter()
            .append('polygon')
            .exons(x, y, exon_height);

        svg.selectAll('.intron-retention')
            .data(function () {
                return gene.exons.filter(function (e) {
                    return e.intron_retention
                })
            })
            .enter()
            .append('polygon')
            .intron_retention(x, y, exon_height);

        var exons_no_ir = gene.exons.filter(function (e) {
            return !e.intron_retention
        });

        svg.selectAll('.exon-number')
            .data(exons_no_ir)
            .enter()
            .append('text')
            .exon_numbers(x, y, exon_height, font_size, gene.strand, exons_no_ir.length);

        svg.selectAll('.junction')
            .data(juncs_no_ir)
            .enter()
            .append('path')
            .junctions(x, y, exon_height, gene.strand);

        svg.selectAll('.reads')
            .data(juncs_no_ir)
            .enter()
            .append('text')
            .reads(x, y, exon_height, font_size);

        svg.selectAll('.ir-line')
            .data(function () {
                return gene.junctions.filter(function (j) {
                    return j.intron_retention !== 0
                })
            })
            .enter()
            .append('line')
            .ir_lines(x, y)

    })
};

d3.transition.prototype.exons =
    d3.selection.prototype.exons =
        function (x, y, exon_height) {
            return this
                .attr('class', 'exon')
                .attr('fill-opacity', .3)
                .attr('stroke-linejoin', 'round')
                .attr('fill', function (d) {
                    if (d.exon_type === 0) {
                        return 'grey'
                    } else if (d.exon_type === 1) {
                        return 'green'
                    }
                    return 'None'
                })
                .attr('stroke', function (d) {
                    if (d.exon_type === 1) {
                        return 'green'
                    }
                    // else if (d.exon_type === 2)
                    //     return
                    return 'black'
                })
                .attr('stroke-dasharray', function (d) {
                    if (d.exon_type === 2) {
                        return '5, 2'
                    }
                })
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
        function (x, y, exon_height) {
            return this
                .attr('class', 'intron-retention')
                .attr('fill', 'lightblue')
                .attr('stroke', 'black')
                .attr('points', function (d) {
                    return [
                        [x(d.start), y(exon_height / 4)].join(' '),
                        [x(d.end), y(exon_height / 4)].join(' '),
                        [x(d.end), y(exon_height * (3 / 4))].join(' '),
                        [x(d.start), y(exon_height * (3 / 4))].join(' ')
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
                    return x(d.start + d.length / 2)
                })
                .attr('text-anchor', 'middle')
                .attr('font-family', 'sans-serif')
                .attr('font-size', font_size);
        };


d3.transition.prototype.junctions =
    d3.selection.prototype.junctions = function (x, y, exon_height, strand) {
        var junc_height = 20;
        return this
            .attr('class', 'junction')
            .attr('stroke-width', 1.5)
            .attr('stroke-dasharray', function (d) {
                switch (d.junction_type) {
                    case 3:
                    case 2:
                        return '5,2';
                    case 1:
                        if (d.reads === 0)
                            return '5,2';
                }
            })
            .attr('stroke', function (d) {
                switch (d.junction_type) {
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
            .attr('d', function (d) {
                var sweep_flag = strand === '+' ? 1 : 0;
                return 'M' + [x(d.start), y(exon_height)].join(',') +
                    'A' + [(x(d.end) - x(d.start)) / 2, 10 + junc_height * (d.bin - 1), 0, 0, sweep_flag, x(d.end), y(exon_height)].join(' ')
            })
    };


d3.transition.prototype.reads =
    d3.selection.prototype.reads =
        function (x, y, exon_height, font_size) {
            var junc_height = 20;
            return this
                .filter(function (d) {
                    return !d.intron_retention
                })
                .attr('class', 'reads')
                .text(function (d) {
                    if (d.reads)
                        return d.reads
                })
                .attr('x', function (d) {
                    return x(d.start) + (x(d.end) - x(d.start)) / 2
                })
                .attr('y', function (d) {
                    return y(exon_height + (10 + junc_height * (d.bin - 1)) + 3)
                })
                .attr('text-anchor', 'middle')
                .attr('font-family', 'sans-serif')
                .attr('font-size', font_size);

        };


d3.transition.prototype.ir_lines =
    d3.selection.prototype.ir_lines =
        function (x, y) {
            return this
                .attr('class', 'ir-line')
                .attr('y1', y(18))
                .attr('y2', y(18))
                .attr('stroke', 'black')
                .attr('x1', function (d) {
                    if (d.intron_retention === 1)
                        return x(d.start);
                    else
                        return x(d.start) - 5

                })
                .attr('x2', function (d) {
                    if (d.intron_retention === 1)
                        return x(d.start) + 5;
                    else
                        return x(d.start)
                })
        };


SpliceGraph.prototype.update = function (sg_div, experiment) {
    var gene_id = sg_div.getAttribute('data-gene-id');
    var sg = this;
    this.zoom = sg_div.getAttribute('data-zoom');


    this.db.get(gene_id + '_' + experiment).then(function (gene) {
        var default_view = sg_div.classList.contains('default-view');
        var svg = d3.select(sg_div).select('svg');
        var y = sg.yScale();
        var exon_height = sg.exon_height;
        var reverse_range = sg.gene.strand === '-';
        var x = sg.xScale(default_view, reverse_range);
        var font_size = sg.font_size;
        var juncs_no_ir = sg.junctions_no_ir(gene, x);

        svg
            .interrupt()
            .transition(sg.t)
            .attr('width', sg.width * sg.zoom);

        svg.selectAll('.exon')
            .interrupt()
            .transition(sg.t)
            .exons(x, y, exon_height);

        svg.selectAll('.intron-retention')
            .interrupt()
            .transition(sg.t)
            .intron_retention(x, y, exon_height);

        svg.selectAll('.exon-number')
            .interrupt()
            .transition(sg.t)
            .exon_numbers(x, y, exon_height, font_size, gene.strand);

        svg.selectAll('.junction')
            .interrupt()
            .data(juncs_no_ir)
            .transition(sg.t)
            .junctions(x, y, exon_height, gene.strand);

        svg.selectAll('.reads')
            .interrupt()
            .data(juncs_no_ir)
            .transition(sg.t)
            .reads(x, y, exon_height, font_size);

        svg.selectAll('.ir-line')
            .interrupt()
            .transition(sg.t)
            .ir_lines(x, y)
    })
};