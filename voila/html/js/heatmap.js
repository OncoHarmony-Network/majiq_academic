var stat_color = d3.scaleLog()
    .base(Math.LN10)
    .domain([-1, 1e-10, 1])
    .range(['grey', 'blue', 'white'])
    .interpolate(d3.interpolateCubehelix);

var dpsi_color = d3.scaleLinear()
    .domain([-1, 0, 1])
    .range(['grey', 'white', 'brown'])
    .interpolate(d3.interpolateCubehelix);

var HMData = function (value, stat_name) {
    this.value = value;
    this.name = stat_name;
    if (stat_name.toLowerCase() === 'dpsi')
        this.color_fn = dpsi_color;
    else
        this.color_fn = stat_color
};

HMData.prototype.color = function () {
    return this.color_fn(this.value)
};


var flip = function (arr) {
    return arr.map(function (a) {
        return a.reverse()
    })
};

var rotate = function (arr) {
    return flip(arr)[0].map(function (c, i) {
        return arr.map(function (r) {
            return r[i]
        })
    })
};


class HeatMap {
    constructor(db) {
        this.db = db;
        this.color = new Colors();
    }

    scale(scale) {
        var scale_width = 124;
        var height = 25;

        var svg = d3.select(this.el)
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
        var lsv_id = el.closest('table').dataset.lsvId;
        var junc_idx = el.closest('tr').dataset.junctionIndex;
        var stat_name = el.dataset.statName;
        var cell_size = 20;

        this.db.allDocs({
            keys: ['metadata', lsv_id],
            include_docs: true
        }).then(response => {
            var meta = response.rows[0].doc;
            var data = response.rows[1].doc;

            var group_names = meta.group_names;

            var ws = data.dpsi[junc_idx];
            var matrix = Array(group_names.length).fill(Array(group_names.length).fill(new HMData(-1, 'dpsi')));

            matrix = matrix.map(function (a, i) {
                var n = ws.slice(0, i);
                ws = ws.slice(i);
                return n.map(function (value) {
                    return new HMData(value, 'dpsi')
                }).concat(a.slice(i))
            });

            // transpose matrix
            matrix = rotate(flip(matrix));

            ws = data[stat_name][junc_idx];

            matrix = matrix.map(function (a, i) {
                var n = ws.slice(0, i);
                ws = ws.slice(i);
                return n.map(function (value) {
                    return new HMData(value, stat_name)
                }).concat(a.slice(i));
            });

            var tool_tip = d3.select('.heat-map-tool-tip');
            if (tool_tip.empty()) {
                tool_tip = d3.select("body")
                    .append("div")
                    .attr('class', 'heat-map-tool-tip')
                    .style("display", "none");
                tool_tip.append('div')
                    .attr('class', 'versus');
                tool_tip.append('div')
                    .attr('class', 'stat-name');
                tool_tip.append('div')
                    .attr('class', 'value')
            }


            var svg = d3.select(el);

            svg.selectAll("*").remove();

            svg.attr('height', cell_size * group_names.length).attr('width', cell_size * group_names.length);
            // svg
            //     .append('defs')
            //     .append('pattern')
            //     .attr('patternUnits', 'userSpaceOnUse')
            //     .attr('class', 'diagonalHatch')
            //     .attr('width', 4)
            //     .attr('height', 4)
            //     .append('path')
            //     .attr('d', "M-1,1 l2,-2 M0,4 l4,-4 M3,5 l2,-2")
            //     .attr('stroke', 'grey')
            //     .attr('stroke-width', 1);

            svg.selectAll('g')
                .data(matrix)
                .enter()
                .append('g')
                .attr('transform', function (d, i) {
                    return 'translate(0,' + (i * cell_size) + ')'
                })
                .attr('data-row', function (d, i) {
                    return group_names[i];
                })
                .selectAll('rect')
                .data(function (d) {
                    return d
                })
                .enter()
                .append('rect')

                .attr('class', 'cell')
                .attr('x', function (d, i) {
                    return cell_size * i
                })
                .attr('y', 0)
                .attr('width', cell_size)
                .attr('height', cell_size)
                .attr('fill', function (d) {
                    if (d.value !== -1)
                        return d.color();
                    else
                        return 'url(#diagonalHatch)'

                })
                .attr('data-column', function (d, i) {
                    if (d.value !== -1)
                        return group_names[i]
                })
                .attr('data-name', function (d) {
                    if (d.value !== -1)
                        return d.name
                })
                .attr('data-value', function (d) {
                    if (d.value !== -1)
                        return d.value
                })
                .on("mouseover", function (d) {
                    if (d.value !== -1) {
                        tool_tip.selectAll('.stat-name').text(this.getAttribute('data-name'));
                        tool_tip.selectAll('.versus').text(this.closest('g').getAttribute('data-row') + '/' + this.getAttribute('data-column'));
                        tool_tip.selectAll('.value').text(this.getAttribute('data-value'));
                        tool_tip.style("display", "block");
                    }
                })
                .on("mousemove", function () {
                    tool_tip.style("top", (event.pageY - 50) + "px").style("left", (event.pageX + 10) + "px");
                })
                .on("mouseout", function () {
                    tool_tip.style("display", "none");
                });
        });

    };

    summary(el, lsv, metadata) {
        const row_height = 30;
        const row_width = 30;
        const row_color_height = 2;
        const row_color_width = 12;
        const stroke_width = 1;
        const padding = 10;
        const juncs_count = lsv.junctions.length;
        const x_axis = 45;
        const tests_count = metadata.stat_names.length;
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
            .attr('fill', (d, i) => dpsi_color(lsv.dpsi[i].reduce((a, b) => Math.max(a, b))))
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
                .attr('fill', (d, i) => stat_color(lsv[s][i].filter(a => a > 0).reduce((a, b) => Math.min(a, b))))
                .attr('stroke', 'lightgrey')
                .attr('stroke-width', stroke_width)
                .attr('shape-rendering', 'crispEdges')
                .attr('stroke-dasharray', (d, i) => i === 0 ? `${row_height + (row_width * 2)},${row_height}` : `0,${row_width},${row_height + row_width},${row_height}`)
                .on('mouseover', (d, i) => console.log(lsv[s][i]));
        });

        g.append('g')
            .attr('class', 'test-names')
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


