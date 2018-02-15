$(function () {
    $('.splice-graph').each(function () {
        sg.init(this)
    })
});

var SpliceGraph = function (db) {
    this.db = db;
    this.width = 1000;
    this.height = 300;

    this.x = {};
    this.y = null;
    this.junction_height = 25;
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
    });

    // adjust spaces between exons
    var ir_count = 0;
    gene.exons.forEach(function (e, i) {
        if (e.intron_retention) ir_count++;

        i -= ir_count;

        length = x_range[i * 2 + 2] - x_range[i * 2 + 1];
        offset = 0;


        if (e.intron_retention) {
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


SpliceGraph.prototype.set_junction_height = function (junctions, x) {
    junctions.sort(function (a, b) {
        var a_length = Math.abs(x(a.start) - x(a.end));
        var b_length = Math.abs(x(b.end) - x(b.start));
        return b_length - a_length;
    });

    var binned_junctions = [];
    while (junctions.length) {
        var junc = junctions.pop();
        junctions.forEach(function (j) {
            if ((j.start <= junc.start) && (j.end >= junc.end))
                j.bin = Math.max(j.bin, junc.bin + 1);
        });
        binned_junctions.push(junc);
    }
    return binned_junctions;
};

SpliceGraph.prototype.distance = function (x, j1, j2) {
    var y = this.yScale();
    var x1 = x(j1.start) + (x(j1.end) - x(j1.start)) / 2;
    var x2 = x(j2.start) + (x(j2.end) - x(j2.start)) / 2;
    var y1 = y(this.exon_height + (this.junction_height * j1.bin) + 3);
    var y2 = y(this.exon_height + (this.junction_height * j2.bin) + 3);
    return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
};

SpliceGraph.prototype.set_junction_distance = function (junctions, x) {
    junctions.sort(function (a, b) {
        var a_length = Math.abs(x(a.start) - x(a.end));
        var b_length = Math.abs(x(b.end) - x(b.start));
        return b_length - a_length;
    });

    var distance_junctions = [];
    var sg = this;
    var changed = false;
    while (junctions.length) {
        var junc = junctions.pop();
        junctions.forEach(function (j) {
            if (sg.distance(x, j, junc) < 25) {
                j.bin += 1;
                changed = true;
            }
        });
        distance_junctions.push(junc)
    }

    return [distance_junctions, changed]
};

SpliceGraph.prototype.junction_bins = function (junctions, x) {
    junctions.forEach(function (j) {
        j.bin = 1;
    });

    junctions = this.set_junction_height(junctions, x);
    var distance_junctions = this.set_junction_distance(junctions, x);

    while(distance_junctions[1]){
        console.log('repeat');
        junctions = this.set_junction_height(distance_junctions[0], x);
        distance_junctions = this.set_junction_distance(junctions, x);
    }


    distance_junctions[0].sort(function (a, b) {
        return a.start - b.start || a.end - b.end;
    });

    return distance_junctions[0]
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

    return this.junction_bins(juncs_no_ir, x);
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

        sg.yScale();
        sg.x = sg.xScale(gene, default_view, reversed_range, experiment);

        var juncs_no_ir = sg.junctions_no_ir(gene, sg.x);
        var exons = gene.exons.filter(function (d) {
            return !d.intron_retention && ![4, 5].includes(gene.exon_types[d.start][d.end][experiment])
        });

        svg.selectAll('.half-exon')
            .data(gene.exons.filter(function (d) {
                return [4, 5].includes(gene.exon_types[d.start][d.end][experiment])
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

        svg.selectAll('.exon')
            .data(exons)
            .enter()
            .append('polygon')
            .attr('class', 'exon');

        svg.selectAll('.exon-number')
            .data(exons)
            .enter()
            .append('text')
            .attr('class', 'exon-number');


        svg.selectAll('.junction-grp')
            .data(juncs_no_ir)
            .enter()
            .append('g')
            .attr('class', 'junction-grp')
            .each(function (d) {
                var junction_grp = d3.select(this);
                junction_grp
                    .selectAll('.junction')
                    .data([d])
                    .enter()
                    .append('path')
                    .attr('class', 'junction');

                junction_grp
                    .selectAll('.reads')
                    .data([d])
                    .enter()
                    .append('text')
                    .attr('class', 'reads');

                junction_grp
                    .selectAll('.splice-site.p3')
                    .data([d])
                    .enter()
                    .append('line')
                    .attr('class', 'splice-site p3');

                junction_grp
                    .selectAll('.splice-site.p5')
                    .data([d])
                    .enter()
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
                    .selectAll('.ir-line')
                    .data([d])
                    .enter()
                    .append('polyline')
                    .attr('class', 'ir-line');

                ir_grp
                    .selectAll('.ir-reads')
                    .data([d])
                    .enter()
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
            var reads = gene.reads;
            var experiment = sg.experiment;
            var x = sg.x;
            var y = sg.y;
            var strand = gene.strand;
            var exon_height = sg.exon_height;
            var font_size = sg.font_size;

            return this
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
                })
                .highlight_exons(lsvs)
        };

d3.transition.prototype.style_exons =
    d3.selection.prototype.style_exons =
        function (sg, gene, lsvs) {
            var exon_types = gene.exon_types;
            var experiment = sg.experiment;

            return this
                .attr('fill-opacity', .3)
                .attr('stroke-linejoin', 'round')
                .each(function (d) {
                    if (lsvs.length) {
                        var colors = new Colors();
                        var hl = lsvs.reduce(function (acc, lsv) {
                            if ((lsv.doc.is_target && lsv.doc.reference_exon[0] - 1 === d.end) ||
                                (!lsv.doc.is_target && lsv.doc.reference_exon[1] + 1 === d.start)) {
                                acc.push(colors.brewer(lsv.doc.junctions.length - 1));
                            }
                            return acc;
                        }, []);

                        if (hl.length > 1)
                            return 'black';

                        if (hl.length === 1) {
                            this.setAttribute('stroke', hl[0]);
                            this.setAttribute('fill', hl[0]);
                            this.setAttribute('stroke-dasharray', '');
                            return
                        }

                    }

                    var exon_type = exon_types[d.start][d.end][experiment];

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
                .highlight_exons(lsvs)
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
                .attr('opacity', function (d) {
                    if (lsvs.length)
                        if (lsvs.every(function (lsv) {
                                if (lsv.doc.is_target)
                                    return lsv.doc.reference_exon[0] - 1 !== d.end;
                                else
                                    return lsv.doc.reference_exon[1] + 1 !== d.start;
                            }))
                            return '0.05'
                })
                .style_exons(sg, gene, lsvs)

        };

d3.transition.prototype.highlight_junctions =
    d3.selection.prototype.highlight_junctions = function (lsvs) {
        return this
            .attr('opacity', function (d) {
                if (lsvs.length)
                    if (lsvs.every(function (lsv) {
                            return lsv.doc.junctions.every(function (junc) {
                                return JSON.stringify(junc) !== JSON.stringify([d.start, d.end])
                            })
                        }))
                        return '0.05'
            })
    };

d3.transition.prototype.highlight_exons =
    d3.selection.prototype.highlight_exons = function (lsvs) {
        return this
            .attr('opacity', function (d) {
                if (lsvs.length)
                    if (lsvs.every(function (lsv) {
                            return JSON.stringify(lsv.doc.reference_exon) !== JSON.stringify([d.start, d.end])
                        }))
                        return '0.05';
            })
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
                                    if (JSON.stringify([d.start, d.end]) === JSON.stringify(junc)) {
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

                    if (lsvs.length) {
                        if (lsvs.some(function (lsv) {
                                return lsv.doc.junctions.some(function (junc) {
                                    return JSON.stringify([d.start, d.end]) === JSON.stringify(junc)
                                })
                            })) {
                            return
                        }
                    }

                    if (this.classList.contains('splice-site'))
                        return '2,2';

                    switch (gene.junction_types[d.start][d.end][experiment]) {
                        case 3:
                        case 2:
                            return '5,2';
                        case 1:
                            if (gene.reads[d.start][d.end][experiment] === 0)
                                return '5,2';
                    }
                })
                .attr('opacity', function (d) {
                    if (lsvs.length) {
                        if (lsvs.some(function (lsv) {
                                return lsv.doc.junctions.some(function (junc) {
                                    return JSON.stringify([d.start, d.end]) === JSON.stringify(junc)
                                })
                            })) {
                            return 'none'
                        } else {
                            return '0.05'
                        }
                    }
                    return 'none';

                })
                .attr('stroke', function (d) {

                    var hl = lsvs.reduce(function (acc, lsv) {
                        return acc.concat(lsv.doc.junctions.reduce(function (acc, junc, idx) {
                            if (JSON.stringify([d.start, d.end]) === JSON.stringify(junc)) {
                                acc.push(colors.brewer(idx))
                            }
                            return acc
                        }, []))
                    }, []);

                    if (hl.length === 1)
                        return hl[0];

                    if (hl.length > 1)
                        return 'black';

                    switch (gene.junction_types[d.start][d.end][experiment]) {
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
                    // where junctions are very long... put them one bin higher.
                    // var long_junc = Math.abs(junc_length) > 200 ? 1 : 0;
                    var long_junc = 0;
                    return 'M' + [x(d.start), y(exon_height)].join(',') +
                        'A' + [junc_length / 2, sg.junction_height * (d.bin + long_junc), 0, 0, sweep_flag, x(d.end), y(exon_height)].join(' ')
                })
        };


d3.transition.prototype.reads =
    d3.selection.prototype.reads =
        function (sg, gene, lsvs) {

            var reads = gene.reads;
            var experiment = sg.experiment;
            var x = sg.x;
            var y = sg.y;
            var exon_height = sg.exon_height;
            var font_size = sg.font_size;

            return this
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
    var sg = this;

    lsv_ids = !lsv_ids ? [] : lsv_ids;
    db.allDocs({
        keys: [sg_div.getAttribute('data-gene-id')].concat(lsv_ids.map(function (x) {
            return x[0]
        })),
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

        var juncs_no_ir = sg.junctions_no_ir(gene, sg.x);

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