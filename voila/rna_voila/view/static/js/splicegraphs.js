const copy_text = str => {
    const el = document.createElement('textarea');
    el.value = str;
    document.body.appendChild(el);
    el.select();
    document.execCommand('copy');
    document.body.removeChild(el);
};

class SpliceGraphs {
    constructor(container, opts) {
        this.container_selector = container;
        this.remove_img = opts.remove_img;
        this.download_img = opts.download_img;
        this.resize_img = opts.resize_img;
        this.remove_fn = opts.remove_fn;
        this.gene = opts.gene;

        this.gene_lr = opts.gene_lr;
        this.gene_c = opts.gene_c;
        this.scaling_transcript = this.gene
        this.lr_sg_height = 40;


        // migrate reads from LR for general look-up
        if(!isEmpty(this.gene_c)){
            this.gene['junction_reads']['combined'] = this.gene_c['junction_reads']['combined']
            this.gene['intron_retention_reads']['combined'] = this.gene_c['intron_retention_reads']['combined']
        }
        if(!isEmpty(this.gene_lr)){
            for(let lr_gene of this.gene_lr){
                this.gene['junction_reads'][lr_gene.experiment] = lr_gene['junction_reads'][lr_gene.experiment]
                this.gene['intron_retention_reads'][lr_gene.experiment] = lr_gene['intron_retention_reads'][lr_gene.experiment]
            }
        }

        this.lsv_ids = [];
        this.zoom = 1;
        this.max_bin = 1;
        this.d = undefined;
        this.lsvs = [];
        this.transcripts = opts.transcripts;
        this.transcript_selector = null;
        this.transcript_sg = null;
        this.transcript_svg = null;
        this.transcript_gene = null;
        this.lr_sg_height = 40;

        if(!isEmpty(this.gene_lr)){
            this.gene_lr.start = this.gene.start;
            this.gene_lr.end = this.gene.end;
        }

        //constants
        this.junction_height = 25;
        this.exon_height = 20;
        this.font_size = 12;
        this.bottom_icons = 10;

        this.initWidth = null;

        //resize event listener
        window.addEventListener('resize', () => this.update());

        // splice_graph_update sg container dataset
        this.container.dataset.geneId = this.gene_id;
        this.container.dataset.zoom = '1';

        this.mutation_observer();
    }

    create_transcripts() {
        if (Object.keys(this.transcripts).length === 0) return;

        this.transcript_sg = document.createElement('div');
        this.container.prepend(this.transcript_sg);

        this.transcript_sg.dataset.group = '__annotation__';
        //sg.dataset.experiment = experiment;
        this.transcript_sg.classList.add('splice-graph');

        const sg_header = d3.select(this.transcript_sg).append('div').attr('class', 'splice-graph-header');

        sg_header
            .append('img')
            .attr('src', this.remove_img)
            //.attr('class', 'splice-graph-remove')
            .attr('height', '16px')
            .style('visibility', 'hidden');

        sg_header
            .append('img')
            .attr('class', 'splice-graph-download')
            .attr('src', this.download_img)
            .attr('height', '16px');

        sg_header
            .append('div')
            .text(`Annotation`);

        this.transcript_selector = sg_header.append("select")
            .attr("name", "transcript-list");

        const options = this.transcript_selector.selectAll("option")
            .data(d3.keys(this.transcripts))
            .enter()
            .append("option");

        options.text(function(d) {
            return d;
        })
        .attr("value", function(d) {
            return d;
        });

        this.transcript_selector.node().onchange = () => {
            this.render_transcript();
        }
        this.render_transcript();


    }

    exons_to_junctions(exons){
        const junctions = [];
        for(let i=0;i<exons.length-1;i++){
            junctions.push({
                'start':exons[i].end,
                'end':exons[i+1].start,
                'color': 'red'
            })
        }
        return junctions;
    }

    render_transcript(){

        if(this.transcript_svg){
            this.transcript_svg.remove();
        }

        this.x = this.x_scale(this.gene);
        //this.junction_bins(sg.dataset.experiment, gene);
        this.y = this.y_scale(this.lr_sg_height + 10);

        this.transcript_svg = d3.select(this.transcript_sg).append('svg')
            .attr('width', this.svg_width)
            .attr('height', this.lr_sg_height)
            .attr("xmlns", "http://www.w3.org/2000/svg");

        const exons = this.transcripts[this.transcript_selector.node().value];

        const g = this.transcript_svg.append('g')
            .attr('transform', `translate(0, ${-this.bottom_icons})`);


        const exon_grps = g.selectAll('.exon-grp')
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

        const junctions = this.exons_to_junctions(exons);

        this.transcript_gene = Object.create(this.gene)
        this.transcript_gene.junctions = junctions;
        this.transcript_gene.exons = exons;

        const junc_grps = g.selectAll('.junction-grp')
            .data(junctions)
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

        this.transcripts_update(this.transcript_sg, this.transcript_gene, []);

        // if there's a scroll bar, then run update one more time to remove it.
        if (document.querySelector('.top').scrollWidth > document.querySelector('.top').clientWidth)
            this.update();
    }

    get container() {
        return document.querySelector(this.container_selector)
    }

    get default_view() {
        return this.container.classList.contains('default-view')
    }

    get width() {
        if(!this.initWidth){
            this.initWidth = this.container.clientWidth - 40;
        }
        return this.initWidth;
        // return this.container.parentNode.offsetWidth
    }

    get svg_height() {
        return (this.max_bin * this.junction_height) + this.exon_height + this.bottom_icons + 20
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

    y_scale(manual_height) {
        const height = (manual_height === undefined ? this.svg_height : manual_height) - 5;
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
        if (reads === undefined) return 0;
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

    intron_retention(sg, lsvs) {
        const x = this.x;
        const y = this.y;
        const exon_height = this.exon_height;

         d3.select(sg).selectAll('.intron-retention')
            .interrupt()
            .transition(this.t())
            .attr('points', function (d) {

                const found = lsvs
                    .filter(l => l.weighted)
                    .filter(l => SpliceGraphs.array_equal(l.intron_retention, [d.start, d.end]));

                if (found.length === 1) {
                    // last group mean should always be intron retention
                    return [
                    [x(d.start), y(exon_height / 2)].join(' '),
                    [x(d.end), y(exon_height / 2)].join(' '),
                    [x(d.end), y(exon_height / 2)].join(' '),
                    [x(d.start), y(exon_height / 2)].join(' ')
                    ].join(', ')
                }else{
                    return [
                    [x(d.start), y(exon_height / 4)].join(' '),
                    [x(d.end), y(exon_height / 4)].join(' '),
                    [x(d.end), y(exon_height * (3 / 4))].join(' '),
                    [x(d.start), y(exon_height * (3 / 4))].join(' ')
                    ].join(', ')
                }
            });
    }

    intron_retention_reads(sg, gene, style) {
        const reads = gene.intron_retention_reads[sg.dataset.experiment];
        const format_reads = this.format_reads;
        d3.select(sg).selectAll('.intron-retention-reads')
            .interrupt()
            .transition(this.t())
            .text(d => {
                const x = this.find_reads(reads, d);
                return x ? format_reads(x) : null;
            })
            .attr('x', d => this.x(d.start + ((d.end - d.start + 1) / 2)))
            .attr('y', (this.y((this.exon_height * (3 / 4)) + 3) - ((style === 'flat') ? 5 : 0)))
            .attr('text-anchor', 'middle')
            .attr('font-family', 'sans-serif')
            .attr('font-size', this.font_size);
    }

    ir_in_lsv(d, lsv) {
        const ir = lsv.intron_retention;
        return ir && d.start === ir[0] && d.end === ir[1]
    }

    style_intron_retention(sg, gene, lsvs) {
        const colors = new Colors();
        const exp = sg.dataset.experiment;
        const grp = sg.dataset.group;

        d3.select(sg).selectAll('.intron-retention-grp')
            .attr('opacity', d => lsvs.length && !lsvs.some(lsv => this.ir_in_lsv(d, lsv)) ? .2 : null);

        d3.select(sg).selectAll('.intron-retention')
            .attr('stroke-width', d => {
                const x = lsvs
                    .filter(l => l.weighted)
                    .filter(l => SpliceGraphs.array_equal(l.intron_retention, [d.start, d.end]));

                let w = 1.5;

                if (x.length === 1) {
                    // last group mean should always be intron retention
                    w = x[0].group_means[grp][exp][x[0].group_means[grp][exp].length-1] * 3;
                }
                return w;

            })
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

    style_denovo_exts(sg) {
        d3.select(sg)
            .selectAll('.denovo-ext')
            .attr('fill', d => {
                return d.ext_color || 'green';
            })
            .attr('fill-opacity', 0.3)
    }

    style_exons(sg, gene, lsvs, noread_opacity_override) {
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
                    case combined_colors['ao']:
                        this.setAttribute('fill', combined_colors['ao']);
                        this.setAttribute('stroke', 'black');
                        break;
                    case combined_colors['l']:
                        this.setAttribute('fill', combined_colors['l']);
                        this.setAttribute('stroke', 'black');
                        break;
                    default:
                        this.setAttribute('fill', 'transparent');
                        this.setAttribute('stroke', 'black');
                        this.setAttribute('stroke-dasharray', '5,2');
                        if(noread_opacity_override){
                            this.setAttribute('opacity', noread_opacity_override);
                        }
                        break;
                }
            });
    }

    ss3p(sg, gene) {
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
    }

    ss5p(sg, gene) {
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
    }

    half_exons(sg) {
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

    };


    exons(sg) {
        const x = this.x;
        const y = this.y;
        const exon_height = this.exon_height;

        d3.select(sg).selectAll('.exon')
            .interrupt()
            .transition(this.t())
            .attr('points', d => {
                return [
                    [x(d.start), y(0)].join(' '),
                    [x(d.end), y(0)].join(' '),
                    [x(d.end), y((exon_height))].join(' '),
                    [x(d.start), y(exon_height)].join(' ')
                ].join(', ')
            });
    }

    denovo_ext(sg) {
        const x = this.x;
        const y = this.y;
        const exon_height = this.exon_height;

        d3.select(sg)
            .selectAll('.denovo-ext.end')
            .interrupt()
            .transition(this.t())
            .attr('points', d => {
                return [
                    [x(d.annotated_end), y(0)].join(' '),
                    [x(d.end), y(0)].join(' '),
                    [x(d.end), y((exon_height))].join(' '),
                    [x(d.annotated_end), y(exon_height)].join(' ')
                ].join(', ');
            });

        d3.select(sg)
            .selectAll('.denovo-ext.start')
            .interrupt()
            .transition(this.t())
            .attr('points', d => {
                return [
                    [x(d.start), y(0)].join(' '),
                    [x(d.annotated_start), y(0)].join(' '),
                    [x(d.annotated_start), y((exon_height))].join(' '),
                    [x(d.start), y(exon_height)].join(' ')
                ].join(', ');
            });
    }

    exon_numbers(sg, gene) {
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
                if (d.color == 'hidden'){
                    return '';
                }
                else if (d.number !== undefined){
                    return d.number;
                }
                else if (strand === '+')
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
    }

    exon_reads(sg, gene, style) {

        const experiment = sg.dataset.experiment;
        const reads = gene.exon_reads[experiment];
        const x = this.x;
        const y = this.y;
        const exon_height = this.exon_height;
        const font_size = this.font_size;
        const junc_height = this.junction_height;

        d3.select(sg).selectAll('.exon-reads')
            .interrupt()
            .data(gene.exons)
            .transition(this.t())
            .text(function (d) {
                try {
                    const r = reads[d.start][d.end];
                    if (r)
                        return r;
                } catch (TypeError) {
                    return '';
                }
            })
            .attr('x', function (d) {
                return x(d.start) + (x(d.end) - x(d.start)) / 2
            })
            .attr('y', function (d) {
                if(style == 'flat'){
                    return y(exon_height) - 5;
                }else{
                    const long_junc = 0;
                    return y(exon_height + (junc_height * (d.bin + long_junc)) + 3)
                }
            })
            .attr('text-anchor', 'middle')
            .attr('font-family', 'monospace')
            .attr('font-size', font_size)
            .attr('fill', '#5c22d0');
    }

    junctions(sg, gene, style) {
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
                if (style === "spike"){
                    return 'M' + [x(d.start), y(exon_height)].join(',') +
                        'L' + [x(d.start) + (junc_length / 2), y(junction_height*1.4 * d.bin)].join(',') +
                        'L' + [x(d.end), y(exon_height)].join(',');
                }else if(style === "flat"){
                    return 'M' + [x(d.start), y(exon_height/2)].join(',') +
                           'L' + [x(d.end), y(exon_height/2)].join(',');
                }{
                    return 'M' + [x(d.start), y(exon_height)].join(',') +
                        'A' + [junc_length / 2, junction_height * d.bin, 0, 0, sweep_flag, x(d.end), y(exon_height)].join(' ');
                }
            });
    }

    style_junctions(sg, gene, lsvs) {

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
                            if(!$.isEmptyObject(x[0].group_means[grp]))
                                w = x[0].group_means[grp][exp][i] * 3;
                        }
                    });
                }

                return w;

            })
            .attr('stroke-dasharray', function (d) {
                if (grp === "__annotation__") return;

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
                            if(grp === 'combined'){
                                if(junc[1] === d.end && lsv.dir === 't'){
                                    acc.push(colors.brewer(idx))
                                }else if(junc[0] === d.start && lsv.dir === 's'){
                                    acc.push(colors.brewer(idx))
                                }
                            }else{
                                if (SpliceGraphs.array_equal(junc, [d.start, d.end])) {
                                    acc.push(colors.brewer(idx))
                                }
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
    }

    alt_starts(sg) {
        const strand = this.gene.strand;
        d3.select(sg)
            .selectAll('.alt_start')
            .interrupt()
            .transition(this.t())
            .attr('x', d => this.x(d))
            .attr('y', this.y(this.bottom_icons - 13))
            .attr('text-anchor', strand === '+' ? 'start' : 'middle')
            .attr('font-weight', 'bold')
            .attr('font-size', 12)
            .text(() => strand === '+' ? '↳' : '^')
    }

    alt_ends(sg) {
        const strand = this.gene.strand;
        d3.select(sg)
            .selectAll('.alt_end')
            .interrupt()
            .transition(this.t())
            .attr('x', d => this.x(d))
            .attr('y', this.y(this.bottom_icons - 13))
            .attr('text-anchor', strand === '+' ? 'middle' : 'start')
            .attr('font-weight', 'bold')
            .attr('font-size', 12)
            .text(() => strand === '+' ? '^' : '↳')
    }

    format_reads(r){
        if (r){
            if(Array.isArray(r)){
                const sr_letter_count = String(r[0]).length;
                const lr_letter_count = String(r[1]).length;
                const max_lettercount = Math.max(sr_letter_count, lr_letter_count);
                const _sr_reads = ' '.repeat(max_lettercount - sr_letter_count) + String(r[0]);
                const _lr_reads = String(r[1]) + ' '.repeat(max_lettercount - lr_letter_count);
                return `${_sr_reads}╦${_lr_reads}`;
            }
            return r;
        }
    }

    junction_reads(sg, gene, style) {

        const experiment = sg.dataset.experiment;
        const reads = gene.junction_reads[experiment];
        const x = this.x;
        const y = this.y;
        const exon_height = this.exon_height;
        const font_size = this.font_size;
        const junc_height = this.junction_height;
        const format_reads = this.format_reads;
        d3.select(sg).selectAll('.junction-reads')
            .interrupt()
            .data(gene.junctions)
            .transition(this.t())
            .text(function (d) {
                try {
                    const r = reads[d.start][d.end];
                    if (r)
                        return format_reads(r);
                } catch (TypeError) {
                    return '';
                }
            })
            .attr('x', function (d) {
                return x(d.start) + (x(d.end) - x(d.start)) / 2
            })
            .attr('y', function (d) {
                if(style == 'flat'){
                    return y(exon_height) - 5;
                }else{
                    const long_junc = 0;
                    return y(exon_height + (junc_height * (d.bin + long_junc)) + 3)
                }
            })
            .attr('text-anchor', 'middle')
            .attr('font-family', 'monospace')
            .attr('font-size', font_size);
    }

    create(group, experiment, _type, transcript) {
        const sg = document.createElement('div');
        this.container.appendChild(sg);

        sg.dataset.group = group;
        sg.dataset.experiment = experiment;
        sg.transcript = (transcript === undefined) ? this.gene : transcript;
        sg.classList.add('splice-graph');

        if(_type == "short_read" || _type == "combined"){
            this.splice_graph_init.call(this, sg, transcript)
        }else if(_type == "long_read"){
            this.splice_graph_init_lr_only.call(this, sg, transcript)
        }

        if(_type == "short_read"){
            send_ajax(base_url+'/psi-splice-graphs', {'add': [sg.dataset.group, sg.dataset.experiment]});
        }

        // if there's a scroll bar, then run update one more time to remove it.
        // if (document.querySelector('.top').scrollWidth > document.querySelector('.top').clientWidth)
        //     this.update();
    }

    init_create() {
        return json_ajax(base_url+'/psi-splice-graphs')
            .then(json => json.forEach(x => this.create(x[0], x[1], "short_read")))
            //.then(() => this.create('combined', 'combined', "combined", this.gene_c))
            .then(() => isEmpty(this.gene_c) ? function(){} : this.create('combined', 'combined', "combined", this.gene_c))

            .then(() => isEmpty(this.gene_lr) ? function(){} : this.gene_lr.forEach(x => this.create('Long Reads', x.experiment, "long_read", x)))
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

    splice_graph_init_lr_only(sg, transcript) {
        const gene = transcript;
        const sg_header = d3.select(sg).append('div').attr('class', 'splice-graph-header');

        sg_header
            .append('img')
            .attr('src', this.remove_img)
            .attr('class', 'splice-graph-remove')
            .attr('height', '16px');

        sg_header
            .append('img')
            .attr('class', 'splice-graph-rescale')
            .attr('src', this.resize_img)
            .attr('height', '16px');

        sg_header
            .append('img')
            .attr('class', 'splice-graph-download')
            .attr('src', this.download_img)
            .attr('height', '16px');

        sg_header
            .append('div')
            .text(`Group: ${sg.dataset.group}; Experiment: ${sg.dataset.experiment};`);

        sg_header
            .append('div')
            .attr('class', 'splice-graph-lr-transcript-reads')
            .text(`Transcript Total Reads: ${transcript.reads}`);

        this.x = this.x_scale(this.gene);
        this.junction_bins(sg.dataset.experiment, gene);
        this.y = this.y_scale(this.lr_sg_height+10);

        const svg = d3.select(sg).append('svg')
            .attr('width', this.svg_width)
            .attr('height', this.lr_sg_height)
            .attr("xmlns", "http://www.w3.org/2000/svg");

        const exons = gene.exons.filter(function (d) {
            return !d.intron_retention && !d.half_exon
        });

        const g = svg.append('g')
            .attr('transform', `translate(0, ${-this.bottom_icons})`);

        const ir_grps = g.selectAll('.intron-retention-grp')
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

        const exon_grps = g.selectAll('.exon-grp')
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

        exon_grps
            .append('text')
            .attr('class', 'exon-reads');

        const denovo_ext_ends = gene.exons
            .filter(e => e.annotated)
            .filter(e => e.end > e.annotated_end);

        const denovo_ext_starts = gene.exons
            .filter(e => e.annotated)
            .filter(e => e.start < e.annotated_start);

        g.selectAll('.denovo-ext-end')
            .data(denovo_ext_ends)
            .enter()
            .append('polygon')
            .attr('class', 'denovo-ext end');

        g.selectAll('.denovo-ext-start')
            .data(denovo_ext_starts)
            .enter()
            .append('polygon')
            .attr('class', 'denovo-ext start');

        const junc_grps = g.selectAll('.junction-grp')
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

        this.splice_graph_update_lr(sg, gene, []);
    }

    splice_graph_init(sg, transcript) {
        const gene = transcript === undefined ? this.gene : transcript;
        const sg_header = d3.select(sg).append('div').attr('class', 'splice-graph-header');

        sg_header
            .append('img')
            .attr('src', this.remove_img)
            .attr('class', 'splice-graph-remove')
            .attr('height', '16px');

        sg_header
            .append('img')
            .attr('class', 'splice-graph-rescale')
            .attr('src', this.resize_img)
            .attr('height', '16px');

        sg_header
            .append('img')
            .attr('class', 'splice-graph-download')
            .attr('src', this.download_img)
            .attr('height', '16px');

        sg_header
            .append('div')
            .text(`Group: ${sg.dataset.group}; Experiment: ${sg.dataset.experiment};`);

        this.x = this.x_scale(gene);
        this.junction_bins(sg.dataset.experiment, gene);
        this.y = this.y_scale();

        const svg = d3.select(sg).append('svg')
            .attr('width', this.svg_width)
            .attr('height', this.svg_height)
            .attr("xmlns", "http://www.w3.org/2000/svg");

        const exons = gene.exons.filter(function (d) {
            return !d.intron_retention && !d.half_exon
        });

        const g = svg.append('g')
            .attr('transform', `translate(0, ${-this.bottom_icons})`);

        g.selectAll('.half-exon')
            .data(gene.exons.filter(function (d) {
                return Boolean(d.half_exon)
            }))
            .enter()
            .append('polyline')
            .attr('class', 'half-exon');

        const ir_grps = g.selectAll('.intron-retention-grp')
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

        const exon_grps = g.selectAll('.exon-grp')
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

        const denovo_ext_ends = gene.exons
            .filter(e => e.annotated)
            .filter(e => e.end > e.annotated_end);

        const denovo_ext_starts = gene.exons
            .filter(e => e.annotated)
            .filter(e => e.start < e.annotated_start);

        g.selectAll('.denovo-ext-end')
            .data(denovo_ext_ends)
            .enter()
            .append('polygon')
            .attr('class', 'denovo-ext end');

        g.selectAll('.denovo-ext-start')
            .data(denovo_ext_starts)
            .enter()
            .append('polygon')
            .attr('class', 'denovo-ext start');

        const junc_grps = g.selectAll('.junction-grp')
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

        svg.selectAll('.alt_start')
            .data(gene.alt_starts)
            .enter()
            .append('text')
            .attr('class', 'alt_start');

        svg.selectAll('.alt_end')
            .data(gene.alt_ends)
            .enter()
            .append('text')
            .attr('class', 'alt_end');

        this.splice_graph_update(sg, gene, []);
    }

    svg(sg, manual_height) {
        d3.select(sg).select('svg')
            .interrupt()
            .transition(this.t())
            .attr('width', this.svg_width)
            .attr('height', manual_height === undefined ? this.svg_height : manual_height);
    }

    transcripts_update(sg, gene, lsvs){
        //update some values
        this.zoom = parseInt(sg.parentNode.dataset.zoom);
        this.x = this.x_scale(this.gene);

        this.y = this.y_scale(this.lr_sg_height + 10);

        // update splice graph
        this.svg(sg, this.lr_sg_height);
        this.exons(sg);
        this.junctions(sg, gene, 'flat');
        this.exon_numbers(sg, gene);


        // add style to Splice Graph elements
        this.style_exons(sg, gene, lsvs);
        this.style_junctions(sg, gene, lsvs);

    }

    splice_graph_update_lr(sg, gene, lsvs) {
        //update some values
        this.zoom = parseInt(sg.parentNode.dataset.zoom);
        this.x = this.x_scale(this.scaling_transcript);
        this.y = this.y_scale(this.lr_sg_height+10);

        // update splice graph
        this.svg(sg, this.lr_sg_height);
        this.exons(sg);
        // this.half_exons(sg);
        this.intron_retention(sg, lsvs);
        this.intron_retention_reads(sg, gene, 'flat');
        this.exon_numbers(sg, gene);
        this.exon_reads(sg, gene, 'flat');
        this.junctions(sg, gene, 'flat');
        this.junction_reads(sg, gene, 'flat');
        // this.ss3p(sg, gene);
        // this.ss5p(sg, gene);
        // this.alt_starts(sg);
        // this.alt_ends(sg);
        this.denovo_ext(sg);

        // add style to Splice Graph elements
        this.style_exons(sg, gene, lsvs, 0.5);
        this.style_junctions(sg, gene, lsvs);
        this.style_intron_retention(sg, gene, lsvs);
        this.style_denovo_exts(sg);
    }

    splice_graph_update(sg, gene, lsvs) {
        //update some values
        this.zoom = parseInt(sg.parentNode.dataset.zoom);
        this.x = this.x_scale(this.scaling_transcript);
        this.junction_bins(sg.dataset.experiment, gene);
        this.y = this.y_scale();

        // update splice graph
        this.svg(sg);
        this.exons(sg);
        this.half_exons(sg);
        this.intron_retention(sg, lsvs);
        this.intron_retention_reads(sg, gene);
        this.exon_numbers(sg, gene);
        this.junctions(sg, gene);
        this.junction_reads(sg, gene);
        this.ss3p(sg, gene);
        this.ss5p(sg, gene);
        this.alt_starts(sg);
        this.alt_ends(sg);
        this.denovo_ext(sg);

        // add style to Splice Graph elements
        this.style_exons(sg, gene, lsvs);
        this.style_junctions(sg, gene, lsvs);
        this.style_intron_retention(sg, gene, lsvs);
        this.style_denovo_exts(sg);
    }

    mutation_observer() {
        window.addEventListener('click', e => {
            const el = e.target.parentNode;
            if (!el.classList.contains('junction-grp') && !el.classList.contains('exon-grp') && !el.classList.contains('intron-retention-grp'))
                document.querySelectorAll('.select, .select-filter')
                    .forEach(x => {
                        x.classList.remove('select');
                        x.classList.remove('select-filter')
                    })
        });

        new MutationObserver(mutation_list => {

            const added_nodes = mutation_list.map(m => m.addedNodes[0]);

            // highlight junctions and intron retentions when you mouse over them
            added_nodes
                .filter(el => el && el.classList && (el.classList.contains('junction-grp') || el.classList.contains('intron-retention-grp') || el.classList.contains('exon-grp')))
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
                            document.querySelectorAll('.exon-grp').forEach(el => {
                                d3.select(el)
                                    .classed('mouseover-filter', d => !((datum.end >= d.start-1 && datum.end <= d.end) || (datum.start >= d.start && datum.start <= d.end + 1)))
                            });
                            if(!el.classList.contains('exon-grp')){
                                document.querySelectorAll('.splice-graph').forEach(el => {
                                    if(el.dataset.group == "Long Reads"){
                                        let found_junc = false;
                                        el.querySelectorAll('.junction-grp, .intron-retention-grp').forEach(el2 => {
                                            const d = d3.select(el2).data()[0];
                                            if(!(d.start !== datum.start || d.end !== datum.end)){
                                                found_junc = true;
                                            }
                                        });
                                        if(!found_junc){
                                            el.style.display = "none";
                                        }
                                    }
                                });
                            }

                        }
                    };

                    el.onmouseout = () => {
                        const coords = document.querySelector('.coordinates');
                        if (!coords.classList.contains('select'))
                            coords.innerHTML = null;

                        el.classList.remove('mouseover');
                        document.querySelectorAll('.junction-grp, .intron-retention-grp').forEach(el => el.classList.remove('mouseover-filter'));
                        document.querySelectorAll('.exon-grp').forEach(el => el.classList.remove('mouseover-filter'));
                        document.querySelectorAll('.splice-graph').forEach(el => el.style.display = "");
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
                            const d = d3.select(el).datum();
                            copy_text(`${this.gene.chromosome}:${d.start}-${d.end}`)
                        } else {
                            el.dispatchEvent(new Event('mouseout'));
                        }
                    }
                });

            // add click event to remove icon
            added_nodes
                .filter(el => el && el.classList && el.classList.contains('splice-graph-remove'))
                .forEach(el => el.onclick = this.remove_fn);


        })
            .observe(document.querySelector(this.container_selector), {
                childList: true,
                subtree: true
            });
    }

    update(duration) {
        this.d = duration;
        this.container
            .querySelectorAll('.splice-graph')
            .forEach(sg => {
                if(sg.dataset.group === 'Long Reads'){
                    this.splice_graph_update_lr(sg, sg.transcript, this.lsvs)
                }else if (sg.dataset.group == "__annotation__"){
                    this.transcripts_update(sg, this.transcript_gene, this.lsvs);
                }else{
                    this.splice_graph_update(sg, sg.transcript, this.lsvs)
                }
            });
        this.d = undefined;
    }

    junctions_filter(gts, lts, gtl, ltl, gtp, ltp, gtpl, ltpl, presence) {
        // args:
        // for reads, greater than short, less than short, greater than long, less than long
        // then for psi, similar order.

        const gene = this.gene;
        gts = parseInt(gts);
        lts = parseInt(lts);
        gtl = parseInt(gtl);
        ltl = parseInt(ltl);
        const gtd = isNaN(gts);
        const ltd = isNaN(lts);
        const gtld = isNaN(gtl);
        const ltld = isNaN(ltl);
        gtp = parseFloat(gtp);
        ltp = parseFloat(ltp);
        gtpl = parseFloat(ltpl);
        ltpl = parseFloat(ltpl);
        const gtpd = isNaN(gtp);
        const ltpd = isNaN(ltp);
        const gtpld = isNaN(gtpl);
        const ltpld = isNaN(ltpl);

        function _determine_reads(sr, lr){

            const passed_gt = gtd || sr >= gts;
            const passed_lt = ltd || sr <= lts;
            const passed_gtl = gtld || lr >= gtl;
            const passed_ltl = ltld || lr <= lts;

            // if both short read filters undefined, only judge on long reads and vice versa
            // otherwise, either long or short may pass in order for junction to pass.
            if(sr === undefined){
                return passed_gtl && passed_ltl
            }else if(lr === undefined){
                return passed_gt && passed_lt;
            }else if(gtld && ltld){
                return passed_gt && passed_lt;
            }else if(gtd && ltd){
                return passed_gtl && passed_ltl;
            }else{
                return (passed_gt && passed_lt) || (passed_gtl && passed_ltl);
            }
        }

        function _determine_psis(sr, lr){
            const passed_gt = gtpd || sr >= gtp;
            const passed_lt = ltpd || sr <= ltp;
            const passed_gtl = gtpld || lr >= gtpl;
            const passed_ltl = ltpld || lr <= ltps;

            // if both short read filters undefined, only judge on long reads and vice versa
            // otherwise, either long or short may pass in order for junction to pass.
            if(sr === undefined){
                return passed_gtl && passed_ltl
            }else if(lr === undefined){
                return passed_gt && passed_lt;
            }else if(gtpld && ltpld){
                return passed_gt && passed_lt;
            }else if(gtpd && ltpd){
                return passed_gtl && passed_ltl;
            }else{
                return (passed_gt && passed_lt) || (passed_gtl && passed_ltl);
            }
        }

        function _determine(sr, lr, sp, lp, _j_presence){
            if (presence && !isInArr(_j_presence, presence)){
                return false;
            }
            return _determine_reads(sr, lr) && _determine_psis(sp, lp);
        }


        this.container
            .querySelectorAll('.splice-graph')
            .forEach(sg => {
                const experiment = sg.dataset.experiment;
                const junction_reads = gene.junction_reads[experiment];
                const intron_retention_reads = gene.intron_retention_reads[experiment];
                const junction_psis = gene.junction_psis;

                d3.selectAll(sg.querySelectorAll('.junction-grp'))
                    .classed('reads-filter', d => {
                        let sr, lr, sp, lp;


                        try {
                            if (Array.isArray(junction_reads[d.start][d.end])){
                                sr = junction_reads[d.start][d.end][0];
                                lr = junction_reads[d.start][d.end][1];
                                sp = junction_psis[d.start][d.end];
                                lp = junction_psis[d.start][d.end];
                            } else {
                                if(sg.dataset.group === "Long Reads"){
                                    sr = undefined;
                                    lr = parseInt(junction_reads[d.start][d.end]) || 0;
                                    sp = undefined;
                                    lp = junction_psis[d.start][d.end] || 0;
                                }else{
                                    sr = parseInt(junction_reads[d.start][d.end]) || 0;
                                    lr = undefined;
                                    sp = junction_psis[d.start][d.end] || 0;
                                    lp = undefined;
                                }
                            }
                        } catch (TypeError) {
                            sr = 0;
                            lr = 0;
                            sp = 0;
                            lp = 0;
                        }
                        return !(_determine(sr, lr, sp, lp, d.presence));
                    })
                d3.selectAll(sg.querySelectorAll('.intron-retention-grp'))
                    .classed('reads-filter', d => {
                        let sr, lr, sp, lp;
                        try {
                            if (Array.isArray(intron_retention_reads[d.start][d.end])){
                                sr = intron_retention_reads[d.start][d.end][0];
                                lr = intron_retention_reads[d.start][d.end][1];
                                sp = junction_psis[d.start][d.end];
                                lp = junction_psis[d.start][d.end];
                            } else {
                                if(sg.dataset.group === "Long Reads"){
                                    sr = undefined;
                                    lr = parseInt(intron_retention_reads[d.start][d.end]) || 0;
                                    sp = undefined;
                                    lp = junction_psis[d.start][d.end] || 0;
                                }else{
                                    sr = parseInt(intron_retention_reads[d.start][d.end]) || 0;
                                    lr = undefined;
                                    sp = junction_psis[d.start][d.end] || 0;
                                    lp = undefined;
                                }
                            }
                        } catch (TypeError) {
                            sr = 0;
                            lr = 0;
                            sp = 0;
                            lp = 0;
                        }
                        return !(_determine(sr, lr, sp, lp, d.presence));
                    })
            })

    }
}