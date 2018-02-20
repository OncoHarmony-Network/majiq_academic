$(function () {
    $('.splice-graph').each(function () {
        sg.init(this)
    })
});

var SpliceGraph = function (db) {
    this.db = db;
    this.width = 1000;
    this.height = 300;
    this.junction_height = 25;
    this.exon_height = 20;
    this.font_size = 12;
    this.zoom = 1;
    this.max_height = this.height - 5;
    this.x = {};
    this.y = d3.scaleLinear().domain([0, this.max_height]).range([this.max_height, 0]);
};

var array_equal = function (a, b) {
    if (a.length !== b.length)
        return false;
    for (var i = 0, l = a.length; i < l; i++) {
        if (a[i] !== b[i])
            return false
    }
    return true
};

var coord_in_exon = function (exon, coord) {
    return coord >= exon.start && coord <= exon.end
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
    var exon;

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
    for (i = 0; i < gene.exons.length; i++) {
        exon = gene.exons[i];
        if (!exon.intron_retention) {
            x_dom.push(exon.start);
            x_dom.push(exon.end);
            x_range.push(x(exon.start));
            x_range.push(x(exon.end));
        }
    }

    // adjust exon sizes
    var filterd_exons = gene.exons.filter(function (d) {
        return !d.intron_retention
    });

    for (i = 0; i < filterd_exons.length; i++) {
        exon = filterd_exons[i];
        var start = x_range[i * 2];
        var end_idx = i * 2 + 1;
        var end = x_range[end_idx];
        length = end - start;
        offset = 0;

        var exon_type = gene.exon_types[experiment][exon.start][exon.end];
        if ([4, 5].includes(exon_type)) {
            min = 1;
            max = 1;
        } else {
            min = 2;
            max = 100;

        }

        if (length < min)
            offset = min - length;

        if (length > max)
            offset = max - length;

        if (offset !== 0)
            for (j = end_idx; j < x_range.length; j++)
                x_range[j] += offset;
    }

    // adjust spaces between exons
    var ir_count = 0;
    for (i = 0; i < gene.exons.length; i++) {
        exon = gene.exons[i];
        if (exon.intron_retention) ir_count++;

        var idx = i - ir_count;

        length = x_range[idx * 2 + 2] - x_range[idx * 2 + 1];
        offset = 0;


        if (exon.intron_retention) {
            min = 20;
            if (length < min)
                offset = min - length;
        } else {
            max = 10;
            min = 1;
            if (length > max)
                offset = max - length;

            if (length < min)
                offset = min - length;
        }


        if (offset !== 0)
            for (j = (idx * 2) + 2; j < x_range.length; j++)
                x_range[j] = x_range[j] + offset
    }

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

SpliceGraph.prototype.distance = function (x, j1, j2) {
    var y = this.y;
    var x1 = x(j1.start) + (x(j1.end) - x(j1.start)) / 2;
    var x2 = x(j2.start) + (x(j2.end) - x(j2.start)) / 2;
    var y1 = y(this.exon_height + (this.junction_height * j1.bin) + 3);
    var y2 = y(this.exon_height + (this.junction_height * j2.bin) + 3);
    return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
};

SpliceGraph.prototype.junction_bins = function (junctions, reads, x) {
    var i;
    var j;
    var small_junc;
    var junc;
    var changed;
    var sg = this;
    var sentinel = 0;

    for (i = 0; i < junctions.length; i++)
        junctions[i].bin = 1;

    junctions.sort(function (a, b) {
        var a_length = Math.abs(x(a.start) - x(a.end));
        var b_length = Math.abs(x(b.end) - x(b.start));
        return a_length - b_length;
    });

    do {
        changed = false;
        sentinel++;

        // Nest larger junctions around smaller ones.
        for (i = 0; i < junctions.length; i++) {
            small_junc = junctions[i];
            for (j = i + 1; j < junctions.length; j++) {
                junc = junctions[j];
                if ((junc.start <= small_junc.start) && (junc.end >= small_junc.end))
                    junc.bin = Math.max(junc.bin, small_junc.bin + 1);
            }
        }

        // Move junctions that are too close.
        for (i = 0; i < junctions.length; i++) {
            small_junc = junctions[i];
            for (j = i + 1; j < junctions.length; j++) {
                junc = junctions[j];
                var small_junc_r = reads[small_junc.start][small_junc.end];
                var junc_r = reads[junc.start][junc.end];
                var reads_length = small_junc_r.toString().length + junc_r.toString().length;
                if (junc.bin === small_junc.bin && sg.distance(x, junc, small_junc) < reads_length * 4) {
                    junc.bin += 1;
                    changed = true;
                }
            }
        }
    } while (changed && sentinel < 10);

    junctions.sort(function (a, b) {
        return a.start - b.start || a.end - b.end;
    });

    return junctions
};

SpliceGraph.prototype.junctions_no_ir = function (gene, x, experiment) {

    var juncs_no_ir = gene.junctions
        .filter(function (j) {
            return !j.intron_retention
        });

    return this.junction_bins(juncs_no_ir, gene.reads[experiment], x);
};

SpliceGraph.prototype.init = function (sg_div) {
    var gene_id = sg_div.getAttribute('data-gene-id');
    var experiment = sg_div.getAttribute('data-experiment');
    var default_view = sg_div.classList.contains('default-view');
    var sg = this;

    sg_div.setAttribute('data-zoom', 1.0);

    this.db.get(gene_id).then(function (gene) {
        var svg = d3.select(sg_div).append('svg')
            .attr('width', sg.width)
            .attr('height', sg.height);

        var reversed_range = gene.strand === '-';

        sg.x = sg.xScale(gene, default_view, reversed_range, experiment);

        var juncs_no_ir = sg.junctions_no_ir(gene, sg.x, experiment);
        var exons = gene.exons.filter(function (d) {
            return !d.intron_retention && ![4, 5].includes(gene.exon_types[experiment][d.start][d.end])
        });

        svg.selectAll('.half-exon')
            .data(gene.exons.filter(function (d) {
                return [4, 5].includes(gene.exon_types[experiment][d.start][d.end])
            }))
            .enter()
            .append('polyline')
            .attr('class', 'half-exon');


        svg.selectAll('.intron-retention')
            .data(gene.exons.filter(function (e) {
                return e.intron_retention
            }))
            .enter()
            .append('polygon')
            .attr('class', 'intron-retention');

        svg.selectAll('.exon-grp')
            .data(exons)
            .enter()
            .append('g')
            .attr('class', 'exon-grp')
            .each(function (d) {
                var exon_grp = d3.select(this);
                exon_grp
                    .datum(d)
                    .append('polygon')
                    .attr('class', 'exon');

                exon_grp
                    .datum(d)
                    .append('text')
                    .attr('class', 'exon-number');
            });


        svg.selectAll('.junction-grp')
            .data(juncs_no_ir)
            .enter()
            .append('g')
            .attr('class', 'junction-grp')
            .each(function (d) {
                var junction_grp = d3.select(this);

                junction_grp
                    .datum(d)
                    .append('path')
                    .attr('class', 'junction');

                junction_grp
                    .datum(d)
                    .append('text')
                    .attr('class', 'reads');

                junction_grp
                    .datum(d)
                    .append('line')
                    .attr('class', 'splice-site p3');

                junction_grp
                    .datum(d)
                    .append('line')
                    .attr('class', 'splice-site p5')
            });

        svg.selectAll('.ir-grp')
            .data(gene.junctions.filter(function (j) {
                return j.intron_retention
            }))
            .enter()
            .append('g')
            .attr('class', 'ir-grp')
            .each(function (d) {
                var ir_grp = d3.select(this);
                ir_grp
                    .datum(d)
                    .append('polyline')
                    .attr('class', 'ir-line');

                ir_grp
                    .datum(d)
                    .append('text')
                    .attr('class', 'ir-reads')
            });
    });

    this.update(sg_div, undefined, 0)
};

d3.transition.prototype.ss3p =
    d3.selection.prototype.ss3p =
        function (sg, gene, lsvs) {
            return this
                .attr('x1', function (d) {
                    return sg.x(d.start)
                })
                .attr('x2', function (d) {
                    return sg.x(d.start)
                })
                .splice_site(sg, gene, lsvs)
        };

d3.transition.prototype.ss5p =
    d3.selection.prototype.ss5p =
        function (sg, gene, lsvs) {
            return this
                .attr('x1', function (d) {
                    return sg.x(d.end)
                })
                .attr('x2', function (d) {
                    return sg.x(d.end)
                })
                .splice_site(sg, gene, lsvs)
        };

d3.transition.prototype.splice_site =
    d3.selection.prototype.splice_site =
        function (sg, gene, lsvs) {
            return this
                .attr('y1', sg.y(0))
                .attr('y2', sg.y(sg.exon_height))
                .attr('stroke', 'black')
                .style_junctions(sg, gene, lsvs)

        };

d3.transition.prototype.ir_reads =
    d3.selection.prototype.ir_reads =
        function (sg, gene, lsvs) {
            var reads = gene.reads[sg.experiment];
            var x = sg.x;
            var y = sg.y;
            var strand = gene.strand;
            var exon_height = sg.exon_height;
            var font_size = sg.font_size;

            return this
                .text(function (d) {
                    var r = reads[d.start][d.end];
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
                .attr('font-size', font_size)
                .highlight_junctions(lsvs)

        };

d3.transition.prototype.half_exons =
    d3.selection.prototype.half_exons =
        function (sg, gene, lsvs) {
            var exon_types = gene.exon_types;
            var experiment = sg.experiment;
            var x = sg.x;
            var y = sg.y;
            var exon_height = sg.exon_height;

            return this
                .style_exons(sg, gene, lsvs)
                .attr('points', function (d) {
                    var exon_type = exon_types[experiment][d.start][d.end];
                    if (exon_type === 4)
                        return [
                            [x(d.end - 10), y(0)].join(' '),
                            [x(d.end), y(0)].join(' '),
                            [x(d.end), y(exon_height)].join(' '),
                            [x(d.end - 10), y(exon_height)].join(' ')
                        ].join(', ');
                    else if (exon_type === 5)
                        return [
                            [x(d.start + 10), y(0)].join(' '),
                            [x(d.start), y(0)].join(' '),
                            [x(d.start), y(exon_height)].join(' '),
                            [x(d.start + 10), y(exon_height)].join(' ')
                        ].join(', ');
                })

        };

d3.transition.prototype.style_exons =
    d3.selection.prototype.style_exons =
        function (sg, gene, lsvs) {
            var exon_types = gene.exon_types;
            var experiment = sg.experiment;

            return this
                .attr('fill-opacity', .3)
                .attr('stroke-linejoin', 'round')
                .highlight_exons(lsvs)
                .each(function (d) {
                    if (lsvs.length) {
                        if (lsvs.some(function (lsv) {
                                return array_equal(lsv.doc.reference_exon, [d.start, d.end])
                            })) {
                            this.setAttribute('stroke', 'orange');
                            this.setAttribute('fill', 'orange');
                            this.setAttribute('stroke-dasharray', '');
                            return
                        }

                        var colors = new Colors();
                        var hl = lsvs.reduce(function (acc, lsv) {
                            if (gene.strand === '+') {
                                if (lsv.doc.is_target) {
                                    if (lsv.doc.reference_exon[0] - 1 === d.end)
                                        acc.push(colors.brewer(lsv.doc.junctions.length - 1));
                                } else {
                                    if (lsv.doc.reference_exon[1] + 1 === d.start)
                                        acc.push(colors.brewer(lsv.doc.junctions.length - 1));
                                }
                            } else {
                                if (lsv.doc.is_target) {
                                    if (lsv.doc.reference_exon[1] + 1 === d.start)
                                        acc.push(colors.brewer(lsv.doc.junctions.length - 1));
                                } else {
                                    if (lsv.doc.reference_exon[0] - 1 === d.end)
                                        acc.push(colors.brewer(lsv.doc.junctions.length - 1));
                                }
                            }
                            return acc;
                        }, []);

                        if (hl.length > 1) {
                            this.setAttribute('stroke', 'black');
                            this.setAttribute('fill', 'grey');
                            this.setAttribute('stroke-dasharray', '');
                            return
                        }

                        if (hl.length === 1) {
                            this.setAttribute('stroke', hl[0]);
                            this.setAttribute('fill', hl[0]);
                            this.setAttribute('stroke-dasharray', '');
                            return
                        }

                    }

                    var exon_type = exon_types[experiment][d.start][d.end];

                    switch (exon_type) {
                        case 0:
                            this.setAttribute('fill', 'grey');
                            break;
                        case 1:
                        case 4:
                        case 5:
                            this.setAttribute('fill', 'green');
                            break;
                        default:
                            this.setAttribute('fill', 'transparent');
                    }

                    switch (exon_type) {
                        case 1:
                        case 4:
                        case 5:
                            this.setAttribute('stroke', 'green');
                            break;
                        default:
                            this.setAttribute('stroke', 'black');
                    }

                    switch (exon_type) {
                        case 2:
                            this.setAttribute('stroke-dasharray', '5,2')
                    }
                })
        };

d3.transition.prototype.exons =
    d3.selection.prototype.exons =
        function (sg, gene, lsvs) {
            var x = sg.x;
            var y = sg.y;

            return this
                .style_exons(sg, gene, lsvs)
                .attr('points', function (d) {
                    return [
                        [x(d.start), y(0)].join(' '),
                        [x(d.end), y(0)].join(' '),
                        [x(d.end), y(sg.exon_height)].join(' '),
                        [x(d.start), y(sg.exon_height)].join(' ')
                    ].join(', ')
                })
        };


d3.transition.prototype.intron_retention =
    d3.selection.prototype.intron_retention =
        function (sg, gene, lsvs) {
            var x = sg.x;
            var y = sg.y;
            var exon_height = sg.exon_height;

            return this
                .attr('points', function (d) {
                    return [
                        [x(d.start - 1), y(exon_height / 4)].join(' '),
                        [x(d.end + 1), y(exon_height / 4)].join(' '),
                        [x(d.end + 1), y(exon_height * (3 / 4))].join(' '),
                        [x(d.start - 1), y(exon_height * (3 / 4))].join(' ')
                    ].join(', ')
                })
                .style_exons(sg, gene, lsvs)

        };


d3.transition.prototype.exon_numbers =
    d3.selection.prototype.exon_numbers =
        function (sg, gene, lsvs) {
            var size = this.size();
            var strand = gene.strand;
            var y = sg.y;
            var x = sg.x;
            var exon_height = sg.exon_height;
            var font_size = sg.font_size;

            return this
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
                .attr('font-size', font_size)
                .highlight_exons(lsvs)
        };

d3.transition.prototype.style_junctions =
    d3.selection.prototype.style_junctions =
        function (sg, gene, lsvs) {
            var colors = new Colors();
            var experiment = sg.experiment;

            return this
                .attr('stroke-width', function (d) {
                    if (lsvs.length) {
                        var hl = lsvs.reduce(function (acc, lsv, idx) {
                            if (sg.weighted[idx])
                                acc = acc.concat(lsv.doc.junctions.reduce(function (acc, junc, idx) {
                                    if (array_equal(junc, [d.start, d.end])) {
                                        acc.push(lsv.doc.group_means_rounded[sg.group][idx] * 3)
                                    }
                                    return acc
                                }, []));
                            return acc
                        }, []);

                        if (hl.length === 1)
                            return hl[0];
                    }
                    return 1.5
                })
                .attr('stroke-dasharray', function (d) {
                    if (this.classList.contains('splice-site'))
                        return '2,2';

                    if (gene.reads[experiment][d.start][d.end] === 0)
                        return '5,2';
                })
                .attr('stroke', function (d) {
                    if (lsvs.length) {
                        var hl = lsvs.reduce(function (acc, lsv) {
                            return acc.concat(lsv.doc.junctions.reduce(function (acc, junc, idx) {
                                if (array_equal(junc, [d.start, d.end])) {
                                    acc.push(colors.brewer(idx))
                                }
                                return acc
                            }, []))
                        }, []);

                        if (hl.length === 1)
                            return hl[0];
                        else if (hl.length > 1) {
                            return 'black'
                        }
                    }

                    switch (gene.junction_types[experiment][d.start][d.end]) {
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
                .attr('fill', 'none')
                .highlight_junctions(lsvs)
        };

d3.transition.prototype.highlight_exons =
    d3.selection.prototype.highlight_exons = function (lsvs) {
        return this.attr('opacity', function (d) {
            if (lsvs.length) {
                if (lsvs.every(function (lsv) {
                        return lsv.doc.junctions.every(function (junc) {
                            return !coord_in_exon(d, junc[0]) && !coord_in_exon(d, junc[1])
                        })
                    })) {
                    return 0.2
                }
            }
            return 1
        })
    };

d3.transition.prototype.highlight_junctions =
    d3.selection.prototype.highlight_junctions = function (lsvs) {
        return this.attr('opacity', function (d) {
            if (lsvs.length) {
                if (lsvs.every(function (lsv) {
                        return lsv.doc.junctions.every(function (junc) {
                            return !array_equal(junc, [d.start, d.end])
                        })
                    })) {
                    return 0.2
                }
            }
            return 1
        })
    };

d3.transition.prototype.junctions =
    d3.selection.prototype.junctions =
        function (sg, gene, lsvs) {

            var x = sg.x;
            var y = sg.y;
            var exon_height = sg.exon_height;
            return this
                .style_junctions(sg, gene, lsvs)
                .attr('d', function (d) {
                    var sweep_flag = gene.strand === '+' ? 1 : 0;
                    var junc_length = x(d.end) - x(d.start);
                    return 'M' + [x(d.start), y(exon_height)].join(',') +
                        'A' + [junc_length / 2, sg.junction_height * d.bin, 0, 0, sweep_flag, x(d.end), y(exon_height)].join(' ')
                })
        };


d3.transition.prototype.reads =
    d3.selection.prototype.reads =
        function (sg, gene, lsvs) {

            var reads = gene.reads[sg.experiment];
            var x = sg.x;
            var y = sg.y;
            var exon_height = sg.exon_height;
            var font_size = sg.font_size;

            return this
                .text(function (d) {
                    try {
                        var r = reads[d.start][d.end];
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
                    // var junc_length = x(d.end) - x(d.start);
                    // var long_junc = Math.abs(junc_length) > 200 ? 1 : 0;
                    var long_junc = 0;
                    return y(exon_height + (sg.junction_height * (d.bin + long_junc)) + 3)
                })
                .attr('text-anchor', 'middle')
                .attr('font-family', 'sans-serif')
                .attr('font-size', font_size)
                .highlight_junctions(lsvs)
        };


d3.transition.prototype.ir_lines =
    d3.selection.prototype.ir_lines =
        function (sg, gene, lsvs) {
            var strand = gene.strand;
            var x = sg.x;
            var y = sg.y;

            return this
                .attr('stroke', 'black')
                .attr('points', function (d) {
                    var offset = sg.exon_height - 2;
                    if (d.intron_retention === 1)
                        if (strand === '+')
                            return [
                                [x(d.start), y(offset)].join(' '),
                                [x(d.start) + 5, y(offset)].join(' ')
                            ].join(', ');
                        else if (strand === '-')
                            return [
                                [x(d.start), y(offset)].join(' '),
                                [x(d.start) - 5, y(offset)].join(' ')
                            ].join(', ');
                    if (d.intron_retention === 2)
                        if (strand === '+')
                            return [
                                [x(d.end), y(offset)].join(' '),
                                [x(d.end) - 5, y(offset)].join(' ')
                            ].join(', ');
                        else if (strand === '-')
                            return [
                                [x(d.end), y(offset)].join(' '),
                                [x(d.end) + 5, y(offset)].join(' ')
                            ].join(', ');
                })
                .highlight_junctions(lsvs)
        };


SpliceGraph.prototype.update = function (sg_div, lsv_ids, duration) {
    lsv_ids = !lsv_ids ? [] : lsv_ids;

    var sg = this;
    var keys = [sg_div.getAttribute('data-gene-id')].concat(lsv_ids.map(function (x) {
        return x[0]
    }));


    db.allDocs({
        keys: keys,
        include_docs: true
    }, function (err, data) {

        if (err)
            console.error(err);

        var svg = d3.select(sg_div).select('svg');
        var gene = data.rows[0].doc;
        var lsvs = data.rows.slice(1);

        var default_view = sg_div.classList.contains('default-view');
        var reverse_range = gene.strand === '-';

        sg.zoom = sg_div.getAttribute('data-zoom');
        sg.experiment = sg_div.getAttribute('data-experiment');
        sg.x = sg.xScale(gene, default_view, reverse_range, sg.experiment);
        sg.weighted = lsv_ids.map(function (x) {
            return x[1]
        });
        sg.group = sg_div.getAttribute('data-group');

        var juncs_no_ir = sg.junctions_no_ir(gene, sg.x, sg.experiment);

        var t = d3.transition()
            .duration(function () {
                if (duration === undefined)
                    return 125;
                return duration
            });

        svg
            .interrupt()
            .transition(t)
            .attr('width', sg.width * sg.zoom);

        svg.selectAll('.exon')
            .interrupt()
            .transition(t)
            .exons(sg, gene, lsvs);

        svg.selectAll('.half-exon')
            .interrupt()
            .transition(t)
            .half_exons(sg, gene, lsvs);

        svg.selectAll('.intron-retention')
            .interrupt()
            .transition(t)
            .intron_retention(sg, gene, lsvs);

        svg.selectAll('.exon-number')
            .interrupt()
            .transition(t)
            .exon_numbers(sg, gene, lsvs);

        svg.selectAll('.junction')
            .interrupt()
            .data(juncs_no_ir)
            .transition(t)
            .junctions(sg, gene, lsvs);

        svg.selectAll('.reads')
            .interrupt()
            .data(juncs_no_ir)
            .transition(t)
            .reads(sg, gene, lsvs);

        svg.selectAll('.ir-line')
            .interrupt()
            .transition(t)
            .ir_lines(sg, gene, lsvs);

        svg.selectAll('.ir-reads')
            .interrupt()
            .data(gene.junctions.filter(function (d) {
                return d.intron_retention
            }))
            .transition(t)
            .ir_reads(sg, gene, lsvs);

        svg.selectAll('.splice-site.p3')
            .interrupt()
            .data(juncs_no_ir)
            .transition(t)
            .ss3p(sg, gene, lsvs);

        svg.selectAll('.splice-site.p5')
            .interrupt()
            .data(juncs_no_ir)
            .transition(t)
            .ss5p(sg, gene, lsvs);
    })

};