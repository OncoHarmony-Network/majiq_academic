var SpliceGraph = function (db, div, opts) {
    var sg = this;
    opts = !opts ? {} : opts;
    var lsv_ids = !opts.lsv_ids ? [] : opts.lsv_ids;
    this.remove_fn = opts.remove_fn;
    this.callback = opts.callback;
    this.db = db;
    this.sg_div = div;
    this.width = div.parentNode.offsetWidth - 50;
    this.junction_height = 25;
    this.exon_height = 20;
    this.font_size = 12;
    this.max_height = this.height - 5;
    this.max_bin = 0;
    this.weighted = lsv_ids.map(function (x) {
        return x[1]
    });

    var sg_container = div.closest('.splice-graph-container');

    this.experiment = div.dataset.experiment;
    this.group = div.getAttribute('data-group');

    this.zoom = sg_container.getAttribute('data-zoom');
    var gene_id = sg_container.getAttribute('data-gene-id');

    var default_view = sg_container.classList.contains('default-view');

    var keys = [gene_id].concat(lsv_ids.map(function (l) {
        return l[0]
    }));

    db.allDocs({
        keys: keys,
        include_docs: true
    }).then((results) => {

        const gene = results.rows[0].doc;
        const lsvs = results.rows.slice(1).map(function (lsv) {
            return lsv.doc
        });

        sg.x = sg.xScale(gene, default_view);
        sg.juncs_no_ir = sg.junctions_no_ir(gene);

        if (opts.hasOwnProperty('init')) {
            sg.init(gene);
        } else {
            sg.update(gene, lsvs)
        }

            this.callback()
    })


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



var start_end_sort = function (a, b) {
    return a.start - b.start || a.end - b.end;
};

SpliceGraph.prototype.svg_height = function () {
    return (this.max_bin * this.junction_height) + this.exon_height + 20
};

SpliceGraph.prototype.svg_width = function () {
    return this.width * this.zoom
};

SpliceGraph.prototype.yScale = function () {
    var height = this.svg_height() - 5;
    return d3.scaleLinear()
        .domain([0, height])
        .range([height, 0]);
};

SpliceGraph.prototype.xScale = function (gene, default_view) {
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
    var reverse_range = gene.strand === '-';

    // if we're not using the default view, the x-scale if very simple.
    if (!default_view) {
        x_range = [min_width, max_width];
        if (reverse_range)
            x_range.reverse();
        return d3.scaleLinear().domain([gene.start, gene.end]).range(x_range);
    }

    // general x-scale
    var x = d3.scaleLinear().domain([gene.start, gene.end]).range([min_width, max_width]);

    gene.exons.sort(start_end_sort);

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

        if (exon.half_exon) {
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

    // scale back to view group_width
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
    var y = this.yScale();
    var x1 = x(j1.start) + (x(j1.end) - x(j1.start)) / 2;
    var x2 = x(j2.start) + (x(j2.end) - x(j2.start)) / 2;
    var y1 = y(this.exon_height + (this.junction_height * j1.bin) + 3);
    var y2 = y(this.exon_height + (this.junction_height * j2.bin) + 3);
    return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
};

SpliceGraph.prototype.junction_bins = function (junctions, reads) {
    var x = this.x;
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

    this.max_bin = 0;

    do {
        changed = false;
        sentinel++;

        // Nest larger junctions around smaller ones.
        for (i = 0; i < junctions.length; i++) {
            small_junc = junctions[i];
            for (j = i + 1; j < junctions.length; j++) {
                junc = junctions[j];
                if ((junc.start <= small_junc.start) && (junc.end >= small_junc.end)) {
                    junc.bin = Math.max(junc.bin, small_junc.bin + 1);
                    this.max_bin = Math.max(this.max_bin, junc.bin)
                }
            }
        }

        // Move junctions that are too close.
        for (i = 0; i < junctions.length; i++) {
            small_junc = junctions[i];
            for (j = i + 1; j < junctions.length; j++) {
                junc = junctions[j];
                var small_junc_r = reads[small_junc.start][small_junc.end];
                var junc_r = reads[junc.start][junc.end];
                if (small_junc_r && junc_r) {
                    var reads_length = small_junc_r.toString().length + junc_r.toString().length;
                    if (junc.bin === small_junc.bin && sg.distance(x, junc, small_junc) < reads_length * 4) {
                        junc.bin += 1;
                        changed = true;
                        this.max_bin = Math.max(this.max_bin, junc.bin)
                    }
                }
            }
        }
    } while (changed && sentinel < 10);

    junctions.sort(start_end_sort);

    return junctions
};

SpliceGraph.prototype.junctions_no_ir = function (gene) {
    var juncs_no_ir = gene.junctions
        .filter(function (j) {
            return !j.intron_retention
        });

    return this.junction_bins(juncs_no_ir, gene.junction_reads[this.experiment]);
};

SpliceGraph.prototype.intron_retention = function () {
    var x = this.x;
    var y = this.yScale();
    var exon_height = this.exon_height;

    var ir = this.svg.selectAll('.intron-retention')
        .interrupt()
        .transition(this.t)
        .attr('points', function (d) {
            return [
                [x(d.start - 1), y(exon_height / 4)].join(' '),
                [x(d.end + 1), y(exon_height / 4)].join(' '),
                [x(d.end + 1), y(exon_height * (3 / 4))].join(' '),
                [x(d.start - 1), y(exon_height * (3 / 4))].join(' ')
            ].join(', ')
        });

    this.style_exons(ir);

};

SpliceGraph.prototype.style_exons = function (exons) {
    var lsvs = this.lsvs;
    var gene = this.gene;

    exons
        .attr('fill-opacity', .3)
        .attr('stroke-linejoin', 'round')
        .each(function (d) {
            if (lsvs.length) {
                if (lsvs.some(function (lsv) {
                    return array_equal(lsv.reference_exon, [d.start, d.end])
                })) {
                    this.setAttribute('stroke', 'orange');
                    this.setAttribute('fill', 'orange');
                    this.setAttribute('stroke-dasharray', '');
                    return
                }

                var colors = new Colors();
                var hl = lsvs.reduce(function (acc, lsv) {
                    if (gene.strand === '+') {
                        if (lsv.is_target) {
                            if (lsv.reference_exon[0] - 1 === d.end)
                                acc.push(colors.brewer(lsv.junctions.length - 1));
                        } else {
                            if (lsv.reference_exon[1] + 1 === d.start)
                                acc.push(colors.brewer(lsv.junctions.length - 1));
                        }
                    } else {
                        if (lsv.is_target) {
                            if (lsv.reference_exon[1] + 1 === d.start)
                                acc.push(colors.brewer(lsv.junctions.length - 1));
                        } else {
                            if (lsv.reference_exon[0] - 1 === d.end)
                                acc.push(colors.brewer(lsv.junctions.length - 1));
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

            switch (d.color) {
                case 'green':
                    this.setAttribute('fill', 'green');
                    this.setAttribute('stroke', 'green');
                    break;
                case 'grey':
                    this.setAttribute('fill', 'grey');
                    this.setAttribute('stroke', 'black');
                    break;
                default:
                    this.setAttribute('fill', 'transparent');
                    this.setAttribute('stroke', 'black');
                    this.setAttribute('stroke-dasharray', '5,2');
                    break;
            }

        });
    this.highlight_exons(exons)
};

SpliceGraph.prototype.ss3p = function () {
    var x = this.x;
    var ss = this.svg.selectAll('.splice-site.p3')
        .interrupt()
        .data(this.juncs_no_ir)
        .transition(this.t)
        .attr('x1', function (d) {
            return x(d.start)
        })
        .attr('x2', function (d) {
            return x(d.start)
        });
    this.splice_site(ss);
};


SpliceGraph.prototype.ss5p = function () {
    var x = this.x;

    var ss = this.svg.selectAll('.splice-site.p5')
        .interrupt()
        .data(this.juncs_no_ir)
        .transition(this.t)
        .attr('x1', function (d) {
            return x(d.end)
        })
        .attr('x2', function (d) {
            return x(d.end)
        });
    this.splice_site(ss)
};


SpliceGraph.prototype.splice_site = function (ss) {
    var y = this.yScale();
    var exon_height = this.exon_height;
    ss
        .attr('y1', y(0))
        .attr('y2', y(exon_height))
        .attr('stroke', 'black');
    this.style_junctions(ss)
};


SpliceGraph.prototype.half_exons = function () {
    var x = this.x;
    var y = this.yScale();
    var exon_height = this.exon_height;

    var exons = this.svg.selectAll('.half-exon')
        .interrupt()
        .transition(this.t)
        .attr('points', function (d) {
            if (d.half_exon === 'start')
                return [
                    [x(d.end - 10), y(0)].join(' '),
                    [x(d.end), y(0)].join(' '),
                    [x(d.end), y(exon_height)].join(' '),
                    [x(d.end - 10), y(exon_height)].join(' ')
                ].join(', ');
            else if (d.half_exon === 'end')
                return [
                    [x(d.start + 10), y(0)].join(' '),
                    [x(d.start), y(0)].join(' '),
                    [x(d.start), y(exon_height)].join(' '),
                    [x(d.start + 10), y(exon_height)].join(' ')
                ].join(', ');
        });
    this.style_exons(exons)
};


SpliceGraph.prototype.exons = function () {
    var x = this.x;
    var y = this.yScale();
    var exon_height = this.exon_height;

    var exons = this.svg.selectAll('.exon')
        .interrupt()
        .transition(this.t)
        .attr('points', function (d) {
            return [
                [x(d.start), y(0)].join(' '),
                [x(d.end), y(0)].join(' '),
                [x(d.end), y(exon_height)].join(' '),
                [x(d.start), y(exon_height)].join(' ')
            ].join(', ')
        });
    this.style_exons(exons)
};


SpliceGraph.prototype.exon_numbers = function () {
    var exons_nums = this.svg.selectAll('.exon-number');

    var size = exons_nums.size();
    var strand = this.gene.strand;
    var y = this.yScale();
    var x = this.x;
    var exon_height = this.exon_height;
    var font_size = this.font_size;

    exons_nums
        .interrupt()
        .transition(this.t)
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
    this.highlight_exons(exons_nums)
};


SpliceGraph.prototype.highlight_exons = function (exons) {
    var lsvs = this.lsvs;
    exons.attr('opacity', function (d) {
        if (lsvs.length) {
            if (lsvs.every(function (lsv) {
                return lsv.junctions.every(function (junc) {
                    return !coord_in_exon(d, junc[0]) && !coord_in_exon(d, junc[1])
                })
            })) {
                return 0.2
            }
        }
        return 1
    })
};


SpliceGraph.prototype.highlight_junctions = function (junctions) {
    var lsvs = this.lsvs;

    junctions
        .attr('opacity', function (d) {
            if (lsvs.length) {
                if (lsvs.every(function (lsv) {
                    return lsv.junctions.every(function (junc) {
                        return !array_equal(junc, [d.start, d.end])
                    })
                })) {
                    return 0.2
                }
            }
            return 1
        })
};


SpliceGraph.prototype.junctions = function () {
    var x = this.x;
    var y = this.yScale();
    var exon_height = this.exon_height;
    var gene = this.gene;
    var junction_height = this.junction_height;

    juncs = this.svg.selectAll('.junction')
        .interrupt()
        .data(this.juncs_no_ir)
        .transition(this.t)
        .attr('d', function (d) {
            var sweep_flag = gene.strand === '+' ? 1 : 0;
            var junc_length = x(d.end) - x(d.start);
            return 'M' + [x(d.start), y(exon_height)].join(',') +
                'A' + [junc_length / 2, junction_height * d.bin, 0, 0, sweep_flag, x(d.end), y(exon_height)].join(' ')
        });

    this.style_junctions(juncs);
};

SpliceGraph.prototype.style_junctions = function (junctions) {
    var colors = new Colors();
    var experiment = this.experiment;
    var lsvs = this.lsvs;
    var gene = this.gene;
    var weighted = this.weighted;

    junctions
        .attr('stroke-group_width', function (d) {
            if (lsvs.length) {
                var hl = lsvs.reduce(function (acc, lsv, idx) {
                    if (weighted[idx])
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

            if (gene.junction_reads[experiment][d.start][d.end] === 0)
                return '5,2';
        })
        .attr('stroke', function (d) {
            if (lsvs.length) {
                var hl = lsvs.reduce(function (acc, lsv) {
                    return acc.concat(lsv.junctions.reduce(function (acc, junc, idx) {
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

            return d.color
        })
        .attr('fill', 'none');

    this.highlight_junctions(junctions);

};


SpliceGraph.prototype.junction_reads = function () {

    var reads = this.gene.junction_reads[this.experiment];
    var x = this.x;
    var y = this.yScale();
    var exon_height = this.exon_height;
    var font_size = this.font_size;
    var junc_height = this.junction_height;

    var juncs = this.svg.selectAll('.junction-reads')
        .interrupt()
        .data(this.juncs_no_ir)
        .transition(this.t)
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
            var long_junc = 0;
            return y(exon_height + (junc_height * (d.bin + long_junc)) + 3)
        })
        .attr('text-anchor', 'middle')
        .attr('font-family', 'sans-serif')
        .attr('font-size', font_size);
    this.highlight_junctions(juncs)
};


SpliceGraph.prototype.init = function (gene) {
    var sg_div = this.sg_div;

    var sg_header = d3.select(sg_div).append('div');
    sg_header
        .append('img')
        .style('float', 'right')
        .attr('src', 'img/remove.svg')
        .attr('height', '16px')
        .on('click', () => this.remove_fn(sg_div));

    sg_header
        .append('div')
        .text(`Group: ${this.group}; Experiment: ${this.experiment};`);

    var svg = d3.select(sg_div).append('svg')
        .attr('width', this.width)
        .attr('height', this.svg_height());

    var exons = gene.exons.filter(function (d) {
        return !d.intron_retention && !d.half_exon
    });

    svg.selectAll('.half-exon')
        .data(gene.exons.filter(function (d) {
            return Boolean(d.half_exon)
        }))
        .enter()
        .append('polyline')
        .attr('class', 'half-exon');


    svg.selectAll('.intron-retention')
        .data(gene.intron_retention)
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
        .data(this.juncs_no_ir)
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
                .attr('class', 'junction-reads');

            junction_grp
                .datum(d)
                .append('line')
                .attr('class', 'splice-site p3');

            junction_grp
                .datum(d)
                .append('line')
                .attr('class', 'splice-site p5')
        });


    this.update(gene, [], 0)
};

SpliceGraph.prototype.update = function (gene, lsvs, duration) {
    var t = d3.transition()
        .duration(function () {
            if (duration === undefined)
            // return 125;
                return undefined;
            return duration
        });

    var svg = d3.select(this.sg_div).select('svg');

    this.gene = gene;
    this.lsvs = lsvs;
    this.svg = svg;
    this.t = t;

    svg
        .interrupt()
        .transition(t)
        .attr('width', this.svg_width())
        .attr('height', this.svg_height());

    this.exons();
    this.half_exons();
    this.intron_retention();
    this.exon_numbers();
    this.junctions();
    this.junction_reads();
    this.ss3p();
    this.ss5p();
};