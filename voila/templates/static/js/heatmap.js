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

var HeatMap = function (el) {
    this.el = el;
};

HeatMap.prototype.scale = function (scale) {
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

HeatMap.prototype.dpsi_scale = function () {
    this.scale(dpsi_color)
};

HeatMap.prototype.stat_scale = function () {
    this.scale(stat_color)
};

HeatMap.prototype.plot = function () {
    var el = this.el;
    var lsv_id = el.getAttribute('data-lsv-id');
    var junc_idx = el.closest('tr').getAttribute('data-junction-index');
    var stat_name = el.getAttribute('data-stat-name');
    var cell_size = 20;

    db.allDocs({
        keys: ['metadata', lsv_id],
        include_docs: true
    }, function (err, response) {
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
                .style("visibility", "hidden");
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
        svg.append('pattern')
            .attr('patternUnits', 'userSpaceOnUse')
            .attr('id', 'diagonalHatch')
            .attr('width', 4)
            .attr('height', 4)
            .append('path')
            .attr('d', "M-1,1 l2,-2 M0,4 l4,-4 M3,5 l2,-2")
            .attr('stroke', 'grey')
            .attr('stroke-width', 1);

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
                    tool_tip.style("visibility", "visible");
                }
            })
            .on("mousemove", function () {
                tool_tip.style("top", (event.pageY - 50) + "px").style("left", (event.pageX + 10) + "px");
            })
            .on("mouseout", function () {
                tool_tip.style("visibility", "hidden");
            });
    });
};

HeatMap.prototype.compact = function (el) {
    var lsv_id = el.getAttribute('data-lsv-id');
    var junc_idx = el.getAttribute('data-junction-index');
    var cell_size = 35;
    db.get(lsv_id).then(function (data) {
        var a = [
            {value: data.tnom[junc_idx], name: 'tnom'},
            {value: data.infoscore[junc_idx], name: 'infoscore'},
            {value: data.wilcoxon[junc_idx], name: 'wilcoxon'}
        ];
        var colors = d3.scaleSqrt()
            .domain([0, 1])
            .range(['blue', 'white'])
            .interpolate(d3.interpolateRgb);

        var tt = d3.select('.tool-tip');
        if (tt.empty())
            tt = d3.select('body')
                .append('div')
                .attr('class', 'tool-tip')
                .style('position', 'absolute')
                .style('z-index', '10')
                .style('cursor', 'default')
                .style('visibility', "hidden")
                .style('background-color', 'rgba(255,255,255,.80)')
                .style('padding-left', '3px')
                .style('padding-top', '1px')
                .style('padding-bottom', '1px')
                .style('padding-right', '3px');

        var svg = d3.select(el)
            .attr('height', cell_size);

        svg
            .selectAll('.cell')
            .data(a)
            .enter()
            .append('rect')
            .attr('x', function (d, i) {
                return (cell_size * i)
            })
            .attr('y', 0)
            .attr('width', cell_size)
            .attr('height', cell_size)
            .attr('fill', function (d) {
                return colors(d.value);
            })
            .on('mouseover', function (d) {
                tt.text(d.name + ": " + d.value.toFixed(3)).style("visibility", "visible");
            })
            .on('mousemove', function () {
                tt.style('top', (event.pageY - 10) + 'px').style('left', (event.pageX + 10) + 'px')
            })
            .on("mouseout", function () {
                tt.style("visibility", "hidden");
            });

        svg
            .append('rect')
            .attr('height', cell_size)
            .attr('width')
    })
};