const stat_color = d3.scaleLog()
    .domain([1e-40, 0.05, 1])
    .range(['blue', 'lightblue', 'white'])
    .interpolate(d3.interpolateCubehelixLong);

const dpsi_color = d3.scaleLinear()
    .domain([0, 1])
    .range(['white', 'brown'])
    .interpolate(d3.interpolateCubehelixLong);

const HMData = function (value, stat_name) {
    this.value = value === undefined ? -1 : value;
    this.name = stat_name;
    if (stat_name.toLowerCase() === 'dpsi')
        this.color_fn = dpsi_color;
    else
        this.color_fn = stat_color
};

HMData.prototype.color = function () {
    return this.color_fn(this.value)
};


const flip = function (arr) {
    return arr.map(function (a) {
        return a.reverse()
    })
};

const rotate = function (arr) {
    return flip(arr)[0].map(function (c, i) {
        return arr.map(function (r) {
            return r[i]
        })
    })
};


class HeatMap {
    constructor(data) {
        this.data = data;
        this.color = new Colors();
    }

    scale(scale) {
        const scale_width = 124;
        const height = 25;

        const svg = d3.select(this.el)
            .attr('width', scale_width + 2)
            .attr('height', height);

        svg.append('g')
            .attr('transform', 'translate(1)')
            .selectAll('line')
            .data(d3.range(0, 1, 1 / scale_width))
            .enter()
            .append('line')
            .attr('x1', function (d, i) {
                return i
            })
            .attr('x2', function (d, i) {
                return i
            })
            .attr('y1', 0)
            .attr('y2', height)
            .attr('stroke', function (d) {
                return scale(d)
            });

        svg
            .append('rect')
            .attr('height', height)
            .attr('width', scale_width)
            .attr('stroke', 'black')
            .attr('fill', 'transparent')
            .attr('x', 0)
            .attr('y', 0)
            .attr('stroke-width', 1);
    };

    dpsi_scale() {
        this.scale(dpsi_color)
    };

    stat_scale() {
        this.scale(stat_color)
    };

    plot(el) {
        const hm = this.data.heatmap;
        const grp_names = this.data.group_names;
        const cell_size = 20;

        // we can not use one pattern in the page outside of the SVG to render the hatch pattern, because the
        // chrome PDF renderer is not able to pick it up correctly for export.
        // Therefore, we need to define the pattern in each heatmap SVG. SVG patterns can only be referenced by the ID
        // attribute, so to comply with HTML standards we need a unique id for each of these pattern elements.
        // So we are expecting the LSV ID (which should be unique) to be set on each of these instances to use
        const uniq = this.lsv_id;

        d3.select(el)
            .attr('height', cell_size * hm.length).attr('width', cell_size * hm.length)
            .attr('class', 'heat-map')
            .attr('data-stat-name', this.data.stat_name)

            // this section generates the (repetitive) pattern element for the hatch style
            .append("defs")

            // first block here generates the background pattern on the heatmap from a small raster image
            // (commented out) , second block uses a path. Both solutions look good in the browser but
            // look strange in different ways when previewing pdf. Unsure of a good solution so far.

            // .append('pattern')
            // .attr('patternUnits', "userSpaceOnUse")
            // .attr('id', 'diagonalHatch-' + uniq)
            // .attr('width', '20')
            // .attr('height', '20')
            // .append('image')
            // .attr('xlink:href', "/static/img/hatched.png")
            // .attr('width', '20')
            // .attr('height', '20')

            .append('pattern')
            .attr('patternUnits', "userSpaceOnUse")
            .attr('id', 'diagonalHatch-' + uniq)
            .attr('width', '4')
            .attr('height', '4')
            .append('path')
            .attr('d', "M-1,1 l2,-2 M0,4 l4,-4 M3,5 l2,-2")
            .attr('stroke', 'grey')
            .attr('stroke-width', '1')

            .select(function() { return this.parentNode; })
            .select(function() { return this.parentNode; })
            .select(function() { return this.parentNode; })

            // this to draw the diagonal line
            .append('path')
            .attr('d', `M 0 0 L ${cell_size * hm.length} ${cell_size * hm.length}`)
            .attr('stroke', "black")
            .attr('stroke-width', "1")
            .select(function() { return this.parentNode; })


            .selectAll('g')
            .data(hm)
            .enter()
            .append('g')
            .attr('transform', function (d, i) {
                return 'translate(0,' + (i * cell_size) + ')'
            })
            .attr('data-row-idx', function (d, i) {
                return i;
            })
            .attr('data-row', (d, i) => grp_names[i])
            .selectAll('rect')
            .data(function (d) {
                return d
            })
            .enter()
            .append('rect')
            .attr('class', 'cell')
            .attr('data-column', (d, i) => grp_names[i])
            .attr('data-column-idx', (d, i) => i)
            .attr('data-value', d => d)
            .attr('x', (d, i) => cell_size * i)
            .attr('y', 0)
            .attr('width', cell_size)
            .attr('height', cell_size)
            .attr('fill', function (d, i, a) {
                const x = i;
                const y = a[i].closest('g').dataset.rowIdx;
                if (d >= 0) {
                    if (x - y < 0)
                        return stat_color(d);
                    else
                        return dpsi_color(d);
                }
                else
                    return `url(#diagonalHatch-${uniq})`
                    //return '';
            })
    };

    summary(el, lsv, metadata) {

        const create_matrix = (obj) => {
            const m = [];
            for (const gn1 in obj) {
                for (const gn2 in obj[gn1]) {
                    m.push(obj[gn1][gn2])
                }
            }
            return m
        };

        const row_height = 30;
        const row_width = 30;
        const row_color_height = 2;
        const row_color_width = 12;
        const stroke_width = 1;
        const padding = 10;
        const juncs_count = lsv.junctions.length;
        const x_axis = 45;
        const tests_count = metadata.stat_names.length;
        const dpsi = create_matrix(lsv.dpsi);
        const stat = {};

        metadata.stat_names.forEach(s => stat[s] = create_matrix(lsv[s]));

        const svg = d3.select(el)
            .append('svg')
            .attr('height', (juncs_count * (row_height + (stroke_width * 2))) + (padding * 2) + x_axis)
            .attr('width', ((tests_count + 2) * (row_width + (stroke_width * 2))) + (padding * 2));

        const g = svg
            .append('g')
            .attr('transform', `translate(${padding},${padding})`);

        const j = g.selectAll('.junction')
            .data(lsv.junctions);

        const jg = j.enter().append('g').attr('class', 'junction');

        jg
            .append('rect')
            .attr('width', row_color_width)
            .attr('height', row_color_height)
            .attr('y', (d, i) => (row_height * i) + (row_height / 2) - (row_color_height / 2))
            .attr('x', (row_width / 2) - (row_color_width / 2))
            .attr('fill', (d, i) => this.color.brewer(i));

        jg
            .append('rect')
            .attr('width', row_width)
            .attr('height', row_height)
            .attr('x', row_width)
            .attr('y', (d, i) => row_height * i)
            .attr('fill', (d, i) => dpsi_color(Math.max.apply(null, dpsi.map(r => r[i]).filter(a => a > 0))))
            .attr('stroke', 'lightgrey')
            .attr('stroke-width', stroke_width)
            .attr('shape-rendering', 'crispEdges')
            .attr('stroke-dasharray', (d, i) => i !== 0 ? `0,${row_width},${(row_height * 2) + row_width}` : null);


        metadata.stat_names.forEach((s, si) => {
            jg
                .append('rect')
                .attr('width', row_width)
                .attr('height', row_height)
                .attr('x', (row_width * (si + 2)))
                .attr('y', (d, i) => row_height * i)
                .attr('fill', (d, i) => stat_color(Math.min.apply(null, stat[s].map(r => r[i]).filter(a => a > 0))))
                .attr('stroke', 'lightgrey')
                .attr('stroke-width', stroke_width)
                .attr('shape-rendering', 'crispEdges')
                .attr('stroke-dasharray', (d, i) => i === 0 ? `${row_height + (row_width * 2)},${row_height}` : `0,${row_width},${row_height + row_width},${row_height}`)
        });

        g.append('g')
            .attr('class', 'view-names')
            .selectAll('text')
            .data(['dpsi'].concat(metadata.stat_names))
            .enter()
            .append('text')
            .attr('transform', (d, i) => `rotate(45,${row_width * (i + 1.33)},${row_height * (juncs_count + .33)})`)
            .attr('x', (d, i) => row_width * (i + 1.33))
            .attr('y', row_height * (juncs_count + .33))
            .attr('font-size', 12)
            .text(d => d);
    }
}


