class SpliceGraphs {
    constructor(db_gene, db_lsv, container, gene_id) {
        this.db_gene = db_gene;
        this.db_lsv = db_lsv;
        this.container_selector = container;
        this.lsv_ids = [];
        this.gene_id = gene_id;
        this.zoom = 1;
        this.max_bin = 1;
        this.d = undefined;
        this.highlight_lsvs = [];
        this.weighted_lsvs = [];

        //constants
        this.junction_height = 25;
        this.exon_height = 20;
        this.font_size = 12;

        //resize event listener
        window.addEventListener('resize', () => this.update());
    }

    get gene() {
        if (this.gene_data)
            return new Promise(resolve => resolve(this.gene_data));
        else
            return this.db_gene.get(this.gene_id).then(result => {
                this.gene_data = result;
                return result
            })
    }

    get container() {
        return document.querySelector(this.container_selector)
    }

    get default_view() {
        return this.container.classList.contains('default-view')
    }

    get width() {
        return this.container.parentNode.clientWidth
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
        for (var i = 0, l = a.length; i < l; i++) {
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

    distance(x, j1, j2) {
        const y = this.y_scale();
        const x1 = x(j1.start) + (x(j1.end) - x(j1.start)) / 2;
        const x2 = x(j2.start) + (x(j2.end) - x(j2.start)) / 2;
        const y1 = y(this.exon_height + (this.junction_height * j1.bin) + 3);
        const y2 = y(this.exon_height + (this.junction_height * j2.bin) + 3);
        return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
    }

    junction_bins(experiment, gene) {
        const junctions = gene.junctions;
        const reads = gene.junction_reads[experiment];
        const x = this.x;
        let i;
        let j;
        let small_junc;
        let junc;
        let changed;
        const sg = this;
        let sentinel = 0;

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
                    const small_junc_r = reads[small_junc.start][small_junc.end];
                    const junc_r = reads[junc.start][junc.end];
                    if (small_junc_r && junc_r) {
                        const reads_length = small_junc_r.toString().length + junc_r.toString().length;
                        if (junc.bin === small_junc.bin && sg.distance(x, junc, small_junc) < reads_length * 4) {
                            junc.bin += 1;
                            changed = true;
                            this.max_bin = Math.max(this.max_bin, junc.bin)
                        }
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
                        [x(d.start - 1), y(exon_height / 4)].join(' '),
                        [x(d.end + 1), y(exon_height / 4)].join(' '),
                        [x(d.end + 1), y(exon_height * (3 / 4))].join(' '),
                        [x(d.start - 1), y(exon_height * (3 / 4))].join(' ')
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
                .text(d => reads[d.start][d.end] ? reads[d.start][d.end] : null)
                .attr('x', d => this.x(d.start + ((d.end - d.start) / 2)))
                .attr('y', this.y((this.exon_height * (3 / 4)) + 3))
                .attr('text-anchor', 'middle')
                .attr('font-family', 'sans-serif')
                .attr('font-size', this.font_size);
            resolve()
        })
    }

    style_intron_retention(sg, lsvs) {
        d3.select(sg).selectAll('.intron-retention-reads')
            .attr('opacity', lsvs.length ? 0.2 : null);

        d3.select(sg).selectAll('.intron-retention')
            .attr('fill-opacity', .3)
            .attr('stroke-linejoin', 'round')
            .attr('fill', d => d.color)
            .attr('stroke', d => d.color)
            .attr('opacity', lsvs.length ? 0.2 : null)
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

                        const colors = new Colors();
                        const hl = lsvs.reduce(function (acc, lsv) {
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

    splice_site(sg) {
        const y = this.y;
        const exon_height = this.exon_height;
        d3.select(sg).selectAll('.splice-site.p5, .splice-site.p3')
            .attr('y1', y(0))
            .attr('y2', y(exon_height))
            .attr('stroke', 'black');
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


    highlight_exons(sg, gene, lsvs) {
        return new Promise(resolve => {

            d3.select(sg).selectAll('.exon, .half-exon, .intron-retention')
                .each((d, i, a) => {
                    const el = a[i];
                    if (lsvs.length) {
                        if (lsvs.some(lsv => SpliceGraphs.array_equal(lsv.reference_exon, [d.start, d.end]))) {
                            el.setAttribute('stroke', 'orange');
                            el.setAttribute('fill', 'orange');
                            el.setAttribute('stroke-dasharray', '');
                            return
                        }

                        const colors = new Colors();

                        const hl = lsvs.reduce(function (acc, lsv) {
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
                        }

                    }
                });


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


            resolve()
        });
    };


    highlight_junctions(sg, lsvs) {
        return new Promise(resolve => {
            d3.select(sg).selectAll('.junction, .junction-reads')
                .attr('opacity', d => {
                    console.log(lsvs);
                    if (lsvs.length) {
                        if (lsvs.every(lsv => lsv.junctions.every(junc => !SpliceGraphs.array_equal(junc, [d.start, d.end])))) {
                            return 0.2
                        }
                    }
                    return 1
                });
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
                .attr('d', function (d) {
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
            const experiment = sg.dataset.experiment;
            // const weighted = this.weighted;

            d3.select(sg).selectAll('.junction-grp')
                .attr('opacity', d => lsvs.length && lsvs.every(lsv => lsv.junctions.every(junc => !SpliceGraphs.array_equal(junc, [d.start, d.end]))) ? 0.2 : null);

            d3.select(sg).selectAll('.junction, .splice-site.p5, .splice-site.p3')
                .attr('stroke-width', function (d) {
                    if (lsvs.length) {
                        const hl = lsvs.reduce(function (acc, lsv, idx) {
                            // if (weighted[idx])
                            //     acc = acc.concat(lsv.doc.junctions.reduce(function (acc, junc, idx) {
                            //         if (array_equal(junc, [d.start, d.end])) {
                            //             acc.push(lsv.doc.group_means_rounded[sg.group][idx] * 3)
                            //         }
                            //         return acc
                            //     }, []));
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

    add_localstorage(sg_div) {
        const group = sg_div.dataset.group;
        const experiment = sg_div.dataset.experiment;
        const sg_key = 'splice_graphs';
        let sg = JSON.parse(localStorage.getItem(sg_key));
        if (!sg)
            sg = [];

        const duplicate = !sg.some(el => {
            return el[0] === group && el[1] === experiment;
        });

        if (duplicate) {
            sg.push([group, experiment]);
            localStorage.setItem(sg_key, JSON.stringify(sg));
        }
        return duplicate
    };

    remove_localstorage(sg_div) {
        const group = sg_div.dataset.group;
        const experiment = sg_div.dataset.experiment;
        const sg_key = 'splice_graphs';
        let sg = JSON.parse(localStorage.getItem(sg_key));
        if (sg) {
            sg = sg.filter(el => el[0] !== group || el[1] !== experiment);
            localStorage.setItem(sg_key, JSON.stringify(sg))
        }
    };

    create(group, experiment) {
        const sg = document.createElement('div');
        this.container.appendChild(sg);

        sg.dataset.group = group;
        sg.dataset.experiment = experiment;
        sg.classList.add('splice-graph');


        this.splice_graph_init(sg);
        this.add_localstorage(sg);

        // if there's a scroll bar, then run update one more time to remove it.
        if (document.querySelector('.top').scrollWidth > document.querySelector('.top').clientWidth)
            this.update();
    }

    highlight(highlight, weighted) {
        this.highlight_lsvs = highlight;
        this.weighted_lsvs = weighted;
        this.update()
    }

    create_localstorage() {
        const sg = JSON.parse(localStorage.getItem('splice_graphs'));
        if (sg)
            sg.forEach(el => this.create(el[0], el[1]));
    }

    remove(sg) {
        this.remove_localstorage(sg);
        sg.remove();
    }

    splice_graph_init(sg) {
        this.gene.then(gene => {
            // const gene = this.gene;
            const sg_header = d3.select(sg).append('div').attr('class', 'splice-graph-header');

            sg_header
                .append('img')
                // .style('float', 'right')
                .attr('src', 'img/remove.svg')
                .attr('height', '16px')
                .on('click', () => this.remove(sg));

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

            svg.selectAll('.intron-retention')
                .data(gene.intron_retention)
                .enter()
                .append('polygon')
                .attr('class', 'intron-retention');

            svg.selectAll('.intron-retention-reads')
                .data(gene.intron_retention)
                .enter()
                .append('text')
                .attr('class', 'intron-retention-reads');

            svg.selectAll('.exon-grp')
                .data(exons)
                .enter()
                .append('g')
                .attr('class', 'exon-grp')
                .each(function (d) {
                    const exon_grp = d3.select(this);
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
                .data(gene.junctions)
                .enter()
                .append('g')
                .attr('class', 'junction-grp')
                .each(function (d) {
                    const junction_grp = d3.select(this);

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

            this.splice_graph_update(sg, gene, [])
        })
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
            this.style_intron_retention(sg, lsvs);
        });
    }

    update(duration) {
        this.d = duration;
        this.gene.then(gene => {
            return this.db_lsv.allDocs({
                keys: this.highlight_lsvs,
                include_docs: true
            }).then(result => {
                const lsvs = result.rows.map(r => r.doc);
                return this.container.querySelectorAll('.splice-graph').forEach(sg => this.splice_graph_update(sg, gene, lsvs));
            }).then(() => this.d = undefined);
        })
    }
}