

const stat_color = function(val){
    const stopMin = 'blue';
    const stopMax = 'white';
    const rangeMin = 1E-40;
    const rangeMax = 1;
    if(val < rangeMin){
        return stopMin;
    }else if(val > rangeMax){
        return stopMax;
    }else{
        return d3.scaleLog()
        .domain([rangeMin, 0.05, rangeMax])
        .range([stopMin, 'lightblue', stopMax])
        .interpolate(d3.interpolateCubehelixLong)(val);
    }
};

const calc_dpsi_color = function(val, rangeMin, rangeMax){
    const stopMin = 'white';
    const stopMax = 'brown';
    if(val < rangeMin){
        return stopMin;
    }else if(val > rangeMax){
        return stopMax;
    }else{
        return d3.scaleLinear()
        .domain([rangeMin, rangeMax])
        .range([stopMin, stopMax])
        .interpolate(d3.interpolateCubehelixLong)(val);
    }
};

const dpsi_color = function(val, rangeMin, rangeMax){
    if(rangeMin === undefined){
        rangeMin = -0.5;
    }
    if(rangeMax === undefined){
        rangeMax = 0.5;
    }
    return calc_dpsi_color(val, rangeMin, rangeMax);
};

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

    static change_axis_scale(dPSI, stat) {
        // static method
        // modifies all existing heat maps to use a new scale for their coloring
        // also changes the scale in the legends to reflect what is shown

        $('.heat-map-outer').each(function(){
            const dPSI_legend = $(this).find('.hm-upper-legend');
            dPSI_legend.find('text').eq(1).text(dPSI['min']);
            dPSI_legend.find('text').eq(3).text(dPSI['max']);

            $(this).find('.heat-map g').each(function(i, g){
                $(g).find('rect').each(function(j, _rect){
                    const rect = $(_rect);
                    if (rect.attr('data-value') !== '-1'){
                        const val = parseFloat(rect.attr('data-value'));
                        if (i - j < 0){
                                rect.attr('fill', dpsi_color(val, dPSI['min'], dPSI['max']));
                        }
                    }

                })
            })

        });
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

            .attr('class', 'heat-map-outer')
            .attr('height', (cell_size * hm.length) + 100).attr('width', (cell_size * hm.length) + 125)

            // this section generates the (repetitive) pattern element for the hatch style
            .append("defs")
            .append('linearGradient')
            .attr('id', 'upperGradient')
            .append('stop')
            .attr('offset', '0%')
            .attr('style', 'stop-color:white;stop-opacity:1')
            .select(d3_parent)
            .append('stop')
            .attr('offset', '100%')
            .attr('style', 'stop-color:brown;stop-opacity:1')
            .select(d3_parent)
            .select(d3_parent)
            .append('linearGradient')
            .attr('id', 'lowerGradient')
            .append('stop')
            .attr('offset', '0%')
            .attr('style', 'stop-color:white;stop-opacity:1')
            .select(d3_parent)
            .append('stop')
            .attr('offset', '50%')
            .attr('style', 'stop-color:lightblue;stop-opacity:1')
            .select(d3_parent)
            .append('stop')
            .attr('offset', '100%')
            .attr('style', 'stop-color:blue;stop-opacity:1')
            .select(d3_parent)
            .select(d3_parent)

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

            .select(d3_parent)
            .select(d3_parent)
            .select(d3_parent)

            // upper legend bar

            .append("svg")
            .attr('class', "hm-upper-legend")
            .style('overflow', 'visible')
            .attr('x',  (cell_size * hm.length) + 50 + 10).attr('y', (((cell_size * hm.length) + 100) / 2) - 5)
            .append('rect')
            .attr('x', 0).attr('y', 0).attr('width', 60).attr('height', 15)
            .attr('fill', 'url(#upperGradient)')
            .attr('stroke', 'black').attr('strokeWidth', '1px')
            .select(d3_parent)
            .append('text').style('font-size', '10px').text("dPSI").attr('text-anchor', 'middle').attr('fill', 'black').attr('x', 30).attr('y', "-2")
            .select(d3_parent)
            .append('text').style('font-size', '8px').text("-0.5").attr('text-anchor', 'middle').attr('fill', 'black').attr('x', 0).attr('y', 25)
            .select(d3_parent)
            .append('text').style('font-size', '8px').text("0.0").attr('text-anchor', 'middle').attr('fill', 'black').attr('x', 60 / 2).attr('y', 25)
            .select(d3_parent)
            .append('text').style('font-size', '8px').text("0.5").attr('text-anchor', 'middle').attr('fill', 'black').attr('x', 60).attr('y', 25)
            .select(d3_parent)
            .select(d3_parent)

            // lower legend bar

            .append("svg")
            .attr('class', "hm-lower-legend")
            .style('overflow', 'visible')
            .attr('x', (((cell_size * hm.length) + 100) / 2) - 30).attr('y', ((cell_size * hm.length) + 70))
            .append('rect')
            .attr('x', 0).attr('y', 0).attr('width', 60).attr('height', 15)
            .attr('fill', 'url(#lowerGradient)')
            .attr('stroke', 'black').attr('strokeWidth', '1px')
            .select(d3_parent)
            .append('text').style('font-size', '10px').text(this.data.stat_name).attr('text-anchor', 'middle').attr('fill', 'black').attr('x', 30).attr('y', "-2")
            .select(d3_parent)
            .append('text').style('font-size', '8px').text("1.0").attr('text-anchor', 'middle').attr('fill', 'black').attr('x', 0).attr('y', 25)
            .select(d3_parent)
            .append('text').style('font-size', '8px').text("0.05").attr('text-anchor', 'middle').attr('fill', 'black').attr('x', 60 / 2).attr('y', 25)
            .select(d3_parent)
            .append('text').style('font-size', '8px').text("1E-40").attr('text-anchor', 'middle').attr('fill', 'black').attr('x', 60).attr('y', 25)
            .select(d3_parent)
            .select(d3_parent)

            // left side titles
            .append("svg")
            .attr('class', "hm-left-titles")
            .style('overflow', 'visible')
            .attr('x', 40).attr('y', 68)
            .selectAll('text')
            .data(grp_names)
            .enter()
            .append('text').style('font-size', '10px').text(function(d){return d}).attr('text-anchor', 'middle').attr('fill', 'black').attr('x', 0).attr('y', 0)
            .attr('transform', function(d, i){
                return `translate(0, ${i*cell_size})rotate(-45)`;
            })
            .select(d3_parent)
            .select(d3_parent)

        d3.select(el)
            // top side titles
            .append("svg")
            .attr('class', "hm-top-titles")
            .style('overflow', 'visible')
            .attr('x', 68).attr('y', 40)
            .selectAll('text')
            .data(grp_names)
            .enter()
            .append('text').style('font-size', '10px').text(function(d){return d}).attr('text-anchor', 'middle').attr('fill', 'black').attr('x', 0).attr('y', 0)
            .attr('transform', function(d, i){
                return `translate(${i*cell_size}, 0)rotate(-45)`;
            })
            .select(d3_parent)
            .select(d3_parent)

        d3.select(el)
            // main heatmap svg
            .append("svg")
            .attr('class', 'heat-map')
            .attr('x', 50).attr('y', 50)
            .attr('height', cell_size * hm.length).attr('width', cell_size * hm.length)
            .attr('data-stat-name', this.data.stat_name)






            // this to draw the diagonal line
            .append('path')
            .attr('d', `M 0 0 L ${cell_size * hm.length} ${cell_size * hm.length}`)
            .attr('stroke', "black")
            .attr('stroke-width', "1")
            .select(d3_parent)



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


