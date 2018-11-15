class SpliceGraphs {
    constructor(container, opts) {
        this.container_selector = container;
        this.remove_img = opts.remove_img;
        this.remove_fn = opts.remove_fn;
        this.gene = opts.gene;
        this.lsv_ids = [];
        this.zoom = 1;
        this.max_bin = 1;
        this.d = undefined;
        this.lsvs = undefined;

        //constants
        this.junction_height = 25;
        this.exon_height = 20;
        this.font_size = 12;

        //resize event listener
        window.addEventListener('resize', () => this.update());

        // splice_graph_update sg container dataset
        this.container.dataset.geneId = this.gene_id;
        this.container.dataset.zoom = '1';

        this.mutation_observer();
    }

    get container() {
        return document.querySelector(this.container_selector)
    }

    get default_view() {
        return this.container.classList.contains('default-view')
    }

    get width() {
        return this.container.clientWidth - 40
        // return this.container.parentNode.offsetWidth
    }

    get svg_height() {
        return (this.max_bin * this.junction_height) + this.exon_height + 20
    }

    get svg_width() {
        return this.width * this.zoom
    }


    static start_end_sort(a, b) {
        return a.start - b.start || a.end - b.end;
    }

    static coord_in_exon(exon, coord) {
        return coord >= exon.start && coord <= exon.end
    }

    static array_equal(a, b) {
        if (a.length !== b.length)
            return false;
        for (let i = 0, l = a.length; i < l; i++) {
            if (a[i] !== b[i])
                return false
        }
        return true
    };

    t() {
        return d3.transition().duration(this.d)
    }

    y_scale() {
        const height = this.svg_height - 5;
        return d3.scaleLinear()
            .domain([0, height])
            .range([height, 0]);
    };

    x_scale(gene) {
        const x_dom = [];
        let x_range = [];
        const min_width = 10;
        const max_width = (this.width * this.zoom) - 10;
        let i;
        let j;
        let length;
        let max;
        let min;
        let offset;
        let exon;
        // const gene = this.gene;
        const reverse_range = gene.strand === '-';

        // if we're not using the default view, the x-scale if very simple.
        if (!this.default_view) {
            x_range = [min_width, max_width];
            if (reverse_range)
                x_range.reverse();
            return d3.scaleLinear().domain([gene.start, gene.end]).range(x_range);
        }

        // general x-scale
        let x = d3.scaleLinear().domain([gene.start, gene.end]).range([min_width, max_width]);

        gene.exons.sort(SpliceGraphs.start_end_sort);

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
        const filterd_exons = gene.exons.filter(function (d) {
            return !d.intron_retention
        });

        for (i = 0; i < filterd_exons.length; i++) {
            exon = filterd_exons[i];
            const start = x_range[i * 2];
            const end_idx = i * 2 + 1;
            const end = x_range[end_idx];
            length = end - start + 1;
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
        let ir_count = 0;
        for (i = 0; i < gene.exons.length; i++) {
            exon = gene.exons[i];
            if (exon.intron_retention) ir_count++;

            const idx = i - ir_count;

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

    }

    find_reads(reads, junc) {
        if (junc.start in reads)
            if (junc.end in reads[junc.start])
                return reads[junc.start][junc.end];
        return 0;
    }

    distance(x, j1, j2) {
        const y = this.y_scale();
        const x1 = x(j1.start) + (x(j1.end) - x(j1.start)) / 2;
        const x2 = x(j2.start) + (x(j2.end) - x(j2.start)) / 2;
        const y1 = y(this.exon_height + (this.junction_height * j1.bin) + 3);
        const y2 = y(this.exon_height + (this.junction_height * j2.bin) + 3);
        return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
    }

    junction_bins(experiment, gene) {
        const reads = gene.junction_reads[experiment];
        const junctions = gene.junctions;
        const x = this.x;
        let i;
        let j;
        let small_junc;
        let junc;
        let changed;
        const sg = this;
        let sentinel = 0;
        let small_junc_r;
        let junc_r;

        for (i = 0; i < junctions.length; i++)
            junctions[i].bin = 1;

        junctions.sort(function (a, b) {
            const a_length = Math.abs(x(a.start) - x(a.end));
            const b_length = Math.abs(x(b.end) - x(b.start));
            return a_length - b_length;
        });

        this.max_bin = 1;

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

                    small_junc_r = this.find_reads(reads, small_junc);
                    junc_r = this.find_reads(reads, junc);

                    const reads_length = small_junc_r.toString().length + junc_r.toString().length;

                    if (junc.bin === small_junc.bin && sg.distance(x, junc, small_junc) < reads_length * 4) {
                        junc.bin += 1;
                        changed = true;
                        this.max_bin = Math.max(this.max_bin, junc.bin)
                    }
                }
            }
        } while (changed && sentinel < 10);

        junctions.sort(SpliceGraphs.start_end_sort);

    };

    intron_retention(sg) {
        return new Promise(resolve => {
            const x = this.x;
            const y = this.y;
            const exon_height = this.exon_height;

            d3.select(sg).selectAll('.intron-retention')
                .interrupt()
                .transition(this.t())
                .attr('points', function (d) {
                    return [
                        [x(d.start), y(exon_height / 4)].join(' '),
                        [x(d.end), y(exon_height / 4)].join(' '),
                        [x(d.end), y(exon_height * (3 / 4))].join(' '),
                        [x(d.start), y(exon_height * (3 / 4))].join(' ')
                    ].join(', ')
                });
            resolve()
        })
    };

    intron_retention_reads(sg, gene) {
        return new Promise(resolve => {
            const reads = gene.intron_retention_reads[sg.dataset.experiment];
            d3.select(sg).selectAll('.intron-retention-reads')
                .interrupt()
                .transition(this.t())
                .text(d => {
                    const x = this.find_reads(reads, d);
                    return x ? x : null;
                })
                .attr('x', d => this.x(d.start + ((d.end - d.start + 1) / 2)))
                .attr('y', this.y((this.exon_height * (3 / 4)) + 3))
                .attr('text-anchor', 'middle')
                .attr('font-family', 'sans-serif')
                .attr('font-size', this.font_size);
            resolve()
        })
    }

    ir_in_lsv(d, lsv) {
        const ir = lsv.intron_retention;
        return ir && d.start === ir[0] && d.end === ir[1]
    }

    style_intron_retention(sg, gene, lsvs) {
        const colors = new Colors();

        d3.select(sg).selectAll('.intron-retention-grp')
            .attr('opacity', d => lsvs.length && !lsvs.some(lsv => this.ir_in_lsv(d, lsv)) ? .2 : null);

        d3.select(sg).selectAll('.intron-retention')
            .attr('fill-opacity', .3)
            .attr('stroke-linejoin', 'round')
            .each((d, i, a) => {
                const el = a[i];
                const filter_lsvs = lsvs.filter(lsv => this.ir_in_lsv(d, lsv));
                if (filter_lsvs.length) {
                    if (filter_lsvs.length === 1) {
                        const lsv = filter_lsvs[0];
                        el.setAttribute('fill', colors.brewer(lsv.junctions.length));
                        el.setAttribute('stroke', colors.brewer(lsv.junctions.length));
                    } else {
                        el.setAttribute('fill', 'black');
                        el.setAttribute('stroke', 'black');
                    }
                } else {
                    el.setAttribute('fill', d.color);
                    el.setAttribute('stroke', d.color);
                }
            });
    }

    style_exons(sg, gene, lsvs) {
        return new Promise(resolve => {
            // change opacity for 'hidden' elements
            d3.select(sg).selectAll('.exon, .half-exon, .exon-number')
                .attr('opacity', d => {
                    if (lsvs.length) {
                        if (lsvs.every(lsv => {
                            return lsv.junctions.every(junc => {
                                return !SpliceGraphs.coord_in_exon(d, junc[0]) && !SpliceGraphs.coord_in_exon(d, junc[1])
                            })
                        })) {
                            return 0.2
                        }
                    }
                    return 1
                });


            d3.select(sg).selectAll('.exon, .half-exon')
                .attr('fill-opacity', .3)
                .attr('stroke-linejoin', 'round')
                .each(function (d) {

                    if (lsvs.length) {
                        if (lsvs.some(function (lsv) {
                            return SpliceGraphs.array_equal(lsv.reference_exon, [d.start, d.end])
                        })) {
                            this.setAttribute('stroke', 'orange');
                            this.setAttribute('fill', 'orange');
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
                    resolve()
                });
        })
    };

    ss3p(sg, gene) {
        return new Promise(resolve => {
            const x = this.x;
            const y = this.y;
            const exon_height = this.exon_height;
            d3.select(sg).selectAll('.splice-site.p3')
                .interrupt()
                .data(gene.junctions)
                .transition(this.t())
                .attr('x1', function (d) {
                    return x(d.start)
                })
                .attr('x2', function (d) {
                    return x(d.start)
                }).attr('y1', y(0))
                .attr('y2', y(exon_height))
                .attr('stroke', 'black');
            resolve()
        })
    };

    ss5p(sg, gene) {
        return new Promise(resolve => {
            const x = this.x;
            const y = this.y;
            const exon_height = this.exon_height;
            d3.select(sg).selectAll('.splice-site.p5')
                .interrupt()
                .data(gene.junctions)
                .transition(this.t())
                .attr('x1', function (d) {
                    return x(d.end)
                })
                .attr('x2', function (d) {
                    return x(d.end)
                }).attr('y1', y(0))
                .attr('y2', y(exon_height))
                .attr('stroke', 'black');
            resolve()
        });
    };

    half_exons(sg) {
        return new Promise(resolve => {
            const x = this.x;
            const y = this.y;
            const exon_height = this.exon_height;

            d3.select(sg).selectAll('.half-exon')
                .interrupt()
                .transition(this.t())
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
            resolve()
        })

    };


    exons(sg) {
        return new Promise(resolve => {
            const x = this.x;
            const y = this.y;
            const exon_height = this.exon_height;

            d3.select(sg).selectAll('.exon')
                .interrupt()
                .transition(this.t())
                .attr('points', function (d) {
                    return [
                        [x(d.start), y(0)].join(' '),
                        [x(d.end), y(0)].join(' '),
                        [x(d.end), y(exon_height)].join(' '),
                        [x(d.start), y(exon_height)].join(' ')
                    ].join(', ')
                });
            resolve();
        });
    };


    exon_numbers(sg, gene) {
        return new Promise(resolve => {
            const exons_nums = d3.select(sg).selectAll('.exon-number');
            const size = exons_nums.size();
            const strand = gene.strand;
            const y = this.y;
            const x = this.x;
            const exon_height = this.exon_height;
            const font_size = this.font_size;

            exons_nums
                .interrupt()
                .transition(this.t())
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
            resolve()
        })
    };

    junctions(sg, gene) {
        return new Promise(resolve => {
            const x = this.x;
            const y = this.y;
            const exon_height = this.exon_height;
            const junction_height = this.junction_height;

            d3.select(sg).selectAll('.junction')
                .interrupt()
                .data(gene.junctions)
                .transition(this.t())
                .attr('d', d => {
                    const sweep_flag = gene.strand === '+' ? 1 : 0;
                    const junc_length = x(d.end) - x(d.start);
                    return 'M' + [x(d.start), y(exon_height)].join(',') +
                        'A' + [junc_length / 2, junction_height * d.bin, 0, 0, sweep_flag, x(d.end), y(exon_height)].join(' ')
                });
            resolve()
        });
    };

    style_junctions(sg, gene, lsvs) {
        return new Promise(resolve => {
            const colors = new Colors();
            const exp = sg.dataset.experiment;
            const grp = sg.dataset.group;

            d3.select(sg).selectAll('.junction-grp')
                .attr('opacity', d => lsvs.length && lsvs
                    .every(lsv => lsv.junctions
                        .every(junc => !SpliceGraphs.array_equal(junc, [d.start, d.end]))) ? 0.2 : null);

            d3.select(sg).selectAll('.junction, .splice-site.p5, .splice-site.p3')
                .attr('stroke-width', d => {
                    const x = lsvs
                        .filter(l => l.weighted)
                        .filter(l => l.junctions.some(j => SpliceGraphs.array_equal(j, [d.start, d.end])));

                    let w = 1.5;

                    if (x.length === 1) {
                        x[0].junctions.forEach((j, i) => {
                            if (SpliceGraphs.array_equal(j, [d.start, d.end])) {
                                w = x[0].group_means[grp][i] * 3;
                            }
                        });
                    }

                    return w

                })
                .attr('stroke-dasharray', function (d) {
                    if (this.classList.contains('splice-site'))
                        return '2,2';

                    const junc_reads = gene.junction_reads;
                    if (!(exp in junc_reads && d.start in junc_reads[exp] && d.end in junc_reads[exp][d.start]))
                        return '5,2';

                })
                .attr('stroke', function (d) {
                    if (lsvs.length) {
                        const hl = lsvs.reduce(function (acc, lsv) {
                            return acc.concat(lsv.junctions.reduce(function (acc, junc, idx) {
                                if (SpliceGraphs.array_equal(junc, [d.start, d.end])) {
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
            resolve()
        })
    };


    junction_reads(sg, gene) {
        return new Promise(resolve => {
            const experiment = sg.dataset.experiment;
            const reads = gene.junction_reads[experiment];
            const x = this.x;
            const y = this.y;
            const exon_height = this.exon_height;
            const font_size = this.font_size;
            const junc_height = this.junction_height;

            d3.select(sg).selectAll('.junction-reads')
                .interrupt()
                .data(gene.junctions)
                .transition(this.t())
                .text(function (d) {
                    try {
                        const r = reads[d.start][d.end];
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
                    const long_junc = 0;
                    return y(exon_height + (junc_height * (d.bin + long_junc)) + 3)
                })
                .attr('text-anchor', 'middle')
                .attr('font-family', 'sans-serif')
                .attr('font-size', font_size);
            resolve()
        });
    };

    create(group, experiment) {
        const sg = document.createElement('div');
        this.container.appendChild(sg);

        sg.dataset.group = group;
        sg.dataset.experiment = experiment;
        sg.classList.add('splice-graph');

        this.splice_graph_init(sg);
        send_ajax('/psi-splice-graphs', {'add': [sg.dataset.group, sg.dataset.experiment]});

        // if there's a scroll bar, then run update one more time to remove it.
        if (document.querySelector('.top').scrollWidth > document.querySelector('.top').clientWidth)
            this.update();
    }

    init_create() {
        return json_ajax('/psi-splice-graphs')
            .then(json => json.forEach(x => this.create(x[0], x[1])))
            .then(() => this)
    }

    highlight(lsvs) {
        this.lsvs = lsvs;
        this.update()
    }


    remove(sg) {
        this.remove_localstorage(sg);
        sg.remove();
    }

    splice_graph_init(sg) {
        const gene = this.gene;
        const sg_header = d3.select(sg).append('div').attr('class', 'splice-graph-header');

        sg_header
            .append('img')
            .attr('src', this.remove_img)
            .attr('class', 'splice-graph-remove')
            .attr('height', '16px');

        sg_header
            .append('div')
            .text(`Group: ${sg.dataset.group}; Experiment: ${sg.dataset.experiment};`);

        this.x = this.x_scale(gene);
        this.junction_bins(sg.dataset.experiment, gene);
        this.y = this.y_scale();

        const svg = d3.select(sg).append('svg')
            .attr('width', this.svg_width)
            .attr('height', this.svg_height);

        const exons = gene.exons.filter(function (d) {
            return !d.intron_retention && !d.half_exon
        });

        svg.selectAll('.half-exon')
            .data(gene.exons.filter(function (d) {
                return Boolean(d.half_exon)
            }))
            .enter()
            .append('polyline')
            .attr('class', 'half-exon');

        const ir_grps = svg.selectAll('.intron-retention-grp')
            .data(gene.intron_retention)
            .enter()
            .append('g')
            .attr('class', 'intron-retention-grp');

        ir_grps
            .append('polygon')
            .attr('class', 'intron-retention');

        ir_grps
            .append('text')
            .attr('class', 'intron-retention-reads');

        const exon_grps = svg.selectAll('.exon-grp')
            .data(exons)
            .enter()
            .append('g')
            .attr('class', 'exon-grp');

        exon_grps
            .append('polygon')
            .attr('class', 'exon');

        exon_grps
            .append('text')
            .attr('class', 'exon-number');

        const junc_grps = svg.selectAll('.junction-grp')
            .data(gene.junctions)
            .enter()
            .append('g')
            .attr('class', 'junction-grp');

        junc_grps
            .append('path')
            .attr('class', 'junction');

        junc_grps
            .append('text')
            .attr('class', 'junction-reads');

        junc_grps
            .append('line')
            .attr('class', 'splice-site p3');

        junc_grps
            .append('line')
            .attr('class', 'splice-site p5');

        this.splice_graph_update(sg, gene, []);
    }

    svg(sg) {
        return new Promise(resolve => {
            d3.select(sg).select('svg')
                .interrupt()
                .transition(this.t())
                .attr('width', this.svg_width)
                .attr('height', this.svg_height);
            resolve()
        });
    }

    splice_graph_update(sg, gene, lsvs) {
        //update some values
        this.zoom = parseInt(sg.parentNode.dataset.zoom);
        this.x = this.x_scale(gene);
        this.junction_bins(sg.dataset.experiment, gene);
        this.y = this.y_scale();

        // update splice graph
        return Promise.all([
            this.svg(sg),
            this.exons(sg),
            this.half_exons(sg),
            this.intron_retention(sg),
            this.intron_retention_reads(sg, gene),
            this.exon_numbers(sg, gene),
            this.junctions(sg, gene),
            this.junction_reads(sg, gene),
            this.ss3p(sg, gene),
            this.ss5p(sg, gene)
        ]).then(() => {
            this.style_exons(sg, gene, lsvs);
            this.style_junctions(sg, gene, lsvs);
            this.style_intron_retention(sg, gene, lsvs);
        });
    }

    mutation_observer() {
        new MutationObserver(mutation_list => {

            const added_nodes = Array.from(mutation_list)
                .reduce((acc, curr) => acc.concat(Array.from(curr.addedNodes)), []);

            // highlight junctions and intron retentions when you mouse over them
            added_nodes
                .filter(el => el.classList && (el.classList.contains('junction-grp') || el.classList.contains('intron-retention-grp') || el.classList.contains('exon-grp')))
                .forEach(el => {
                    const datum = d3.select(el).datum();
                    el.onmouseover = () => {
                        const coords = document.querySelector('.coordinates');
                        if (!coords.classList.contains('select')) {
                            coords.innerHTML = `Coordinates: ${datum.start}-${datum.end}; Length: ${datum.end - datum.start + 1}`;

                            el.classList.add('mouseover');
                            document.querySelectorAll('.junction-grp, .intron-retention-grp').forEach(el => {
                                d3.select(el)
                                    .classed('mouseover-filter', d => d.start !== datum.start || d.end !== datum.end)
                            });
                        }
                    };

                    el.onmouseout = () => {
                        const coords = document.querySelector('.coordinates');
                        if (!coords.classList.contains('select'))
                            coords.innerHTML = null;

                        el.classList.remove('mouseover');
                        document.querySelectorAll('.junction-grp, .intron-retention-grp').forEach(el => el.classList.remove('mouseover-filter'));
                    };

                    el.onclick = () => {
                        const click_new = !el.classList.contains('select');

                        document.querySelectorAll('.select-filter, .select').forEach(x => {
                            x.classList.remove('select-filter');
                            x.classList.remove('select')
                        });

                        if (click_new) {
                            el.dispatchEvent(new Event('mouseover'));
                            document.querySelectorAll('.mouseover-filter').forEach(el => el.classList.add('select-filter'));
                            el.classList.add('select');
                            document.querySelector('.coordinates').classList.add('select');
                            this.copy_select_ucsc(el)
                        } else {
                            el.dispatchEvent(new Event('mouseout'));
                        }
                    }
                });

            // add click event to remove icon
            added_nodes
                .filter(el => el.classList && el.classList.contains('splice-graph-remove'))
                .forEach(el => el.onclick = this.remove_fn)
        })
            .observe(document.querySelector(this.container_selector), {
                childList: true,
                subtree: true
            });
    }

    update(duration) {
        const gene = this.gene;
        this.d = duration;
        if (this.lsvs) {
            this.container.querySelectorAll('.splice-graph').forEach(sg => this.splice_graph_update(sg, gene, this.lsvs));
            this.d = undefined;
        } else {
            this.container.querySelectorAll('.splice-graph').forEach(sg => this.splice_graph_update(sg, gene, []));
            this.d = undefined
        }
    }
}