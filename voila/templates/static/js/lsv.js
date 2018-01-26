$(document).on('click', '.lsv-single-compact-percentiles, .psi-violin-plot', function () {
    $(this.parentNode)
        .children('.lsv-single-compact-percentiles, .psi-violin-plot')
        .animate({height: 'toggle', width: 'toggle'});
});

$(document).on('click', '.deltapsi-violin-plot, .excl-incl-rect', function () {
    $(this.parentNode)
        .children('.deltapsi-violin-plot, .excl-incl-rect')
        .animate({height: 'toggle', width: 'toggle'});
});


var Lsv = function (db) {
    this.db = db
};

Lsv.prototype.renderLsvSpliceGraph = function (canvas) {
    var lsv = this;
    this.db.get(canvas.getAttribute('data-gene-id')).then(function (gene) {
        if (canvas.getContext) {
            // Render LSV representation from a string text representing the LSV i.e.: s|1e1.3o3|2e1.2o3|3e1.2o3|i|4e1.1o3  NEW: Intron Retention (i field)
            var ctx = canvas.getContext("2d");

            // Clear previous draw
            ctx.clearRect(0, 0, canvas.width, canvas.height);
            var margins = [1, 1, 1, 1],
                pixel_factor = 1,
                percentage_exon = .2,
                percentage_intron = .15;

            var exon_lsv_number = '';
            var exon_lsv_coords = canvas.getAttribute('data-coord-exon');
            var experiment = canvas.getAttribute('data-experiment');
            var exons_gene = gene.exons.filter(function (d) {
                var exon_type = gene.exon_types[d.start][d.end][experiment];
                return exon_type < 3 && !d.intron_retention;
            });

            if (exon_lsv_coords) {
                exon_lsv_coords = exon_lsv_coords.split(',').map(function (x) {
                    return parseInt(x)
                });


                // Find LSV exon number in the splice graph
                var lsv_exon = exons_gene.find(function (element) {
                    return element.start === exon_lsv_coords[0] && element.end === exon_lsv_coords[1]
                });


                if (lsv_exon) {
                    exon_lsv_number = exons_gene.indexOf(lsv_exon) + 1;

                    // If negative strand, complement the exon ordinal
                    if (gene.strand === '-') {
                        exon_lsv_number = exons_gene.length - exon_lsv_number + 1;
                    }
                }
            }

            var lsv_data = canvas.getAttribute('data-lsv-type');
            var lsvs = lsv_data.split('|');
            var ir_marker = 'i';

            var intron_ret_i = lsvs.indexOf(ir_marker);
            if (intron_ret_i > -1) {
                lsvs.splice(lsvs.indexOf(ir_marker), 1);  // Modifies the array in place
                var num_alt_start_end = 0;
                for (var ff = 1; ff < lsvs.length; ff++) {
                    if (lsvs[ff].indexOf('.') === -1) {
                        num_alt_start_end++;
                    }
                }
                intron_ret_i -= num_alt_start_end;
            }

            // Num exons_obj
            var num_exons = 0;
            var num_ss = 0;
            var ss_reg = {};

            for (var n_ways = 1; n_ways < lsvs.length; n_ways++) {
                var lsvs_fields = lsvs[n_ways].split('e');
                num_exons = Math.max(num_exons, parseInt(lsvs_fields[1][0]));
                num_ss = Math.max(num_ss, parseInt(lsvs_fields[0]));

                if (lsvs_fields[1] === 0) continue;
                var exonNum_ss = lsvs_fields[1].split('.');

                if (exonNum_ss[1].indexOf('o') > 0) {
                    ss_reg[exonNum_ss[0]] = parseInt(exonNum_ss[1].split('o')[1]);

                } else {
                    if (ss_reg[exonNum_ss[0]])
                        ss_reg[exonNum_ss[0]] = Math.max(ss_reg[exonNum_ss[0]], parseInt(exonNum_ss[1]));
                    else
                        ss_reg[exonNum_ss[0]] = parseInt(exonNum_ss[1]);
                }
            }
            num_exons++;  // source or target exon is implicit

            var sourceOrTarget = lsvs[0];
            var area = [canvas.width - margins[0] - margins[1], canvas.height - margins[2] - margins[3]];

            percentage_intron = Math.min(percentage_intron, .4 / (num_exons - 1));
            var exon_width = (area[0] - (num_exons - 1) * area[0] * percentage_intron - margins[0] - margins[1]) / num_exons;
            var start = margins[0];
            var direction = sourceOrTarget === 's' ? 1 : -1;

            // Render exons_obj
            var exons = [];
            for (var i = 0; i < num_exons; i++) {
                var exon = {
                    'coords': [Math.round(start), Math.round(start + exon_width)],
                    'type': (direction > 0 && i === 0 || direction < 0 && i === num_exons - 1 ? 1 : 0)
                };
                exons.push(exon);
                var number_exon = (direction > 0 ? i : i + 1);
                number_exon = (number_exon === 0 || number_exon === num_exons ? exon_lsv_number : '');
                lsv.render_exon(canvas, exon, pixel_factor, margins, percentage_exon, number_exon);
                start += exon_width + percentage_intron * area[0];
            }

            // Render IR
            if (intron_ret_i > -1) {
                var intron = {
                    'coords': (direction > 0
                        ? [margins[0] + exon_width, margins[0] + exon_width + percentage_intron * area[0]]
                        : [margins[0] + (num_exons - 1) * exon_width + (num_exons - 2) * percentage_intron * area[0],
                            margins[0] + (num_exons - 1) * exon_width + (num_exons - 1) * percentage_intron * area[0]])
                };
                lsv.render_intron_retained(canvas, intron, pixel_factor, margins, percentage_exon, intron_ret_i - 1);
            }

            // Render junctions
            var count_starts_ends = 0;
            var index_first = direction > 0 ? 0 : exons.length - 1;
            var coords = exons[index_first].coords;
            var exon_height = percentage_exon * canvas.height;

            // If source, we move the coords to the first splice site
            if (direction > 0) {
                coords[0] += exon_width / 2 + direction * (exon_width / 2) / num_ss;
            }
            coords[1] = canvas.height - exon_height - margins[3];

            // For rendering all splice sites even if they don't have a junction jumping in or out of them,
            // keep a registry of the exons whose ssites have been already rendered
            var rendered_exons = {};
            var previous_ctx = [ctx.lineWidth, ctx.strokeStyle];
            for (n_ways = 1; n_ways < lsvs.length; n_ways++) {
                lsvs_fields = lsvs[n_ways].split('e');
                var ss = parseInt(lsvs_fields[0]),
                    target_e = lsvs_fields[1].split('.')[0],
                    target_ss = parseInt(lsvs_fields[1].split('.')[1]),
                    target_num_ss = null;
                if (lsvs_fields[1].indexOf('o') > 0) {
                    target_num_ss = parseInt(lsvs_fields[1].split('.')[1].split('o')[1]);
                }

                var coords_x_start_e = coords[0] + ((exon_width / 2) / num_ss) * (ss - 1);
                var coords_x_target_e = null;

                if (lsvs_fields[1] === 0) {
                    coords_x_target_e = coords_x_start_e;
                    count_starts_ends++;
                } else {
                    var target_exon = (direction > 0 ? exons[lsvs_fields[1][0]] : exons[lsvs_fields[1][0] - 1]);
                    coords_x_target_e = target_exon.coords[direction > 0 ? 0 : 1];
                }

                var offset_ss_step = (exon_width * 2 / 3) / ss_reg[target_e];
                var offset_ss = 0;
                if (direction > 0) {
                    offset_ss = (target_ss - 1) * offset_ss_step;
                } else {
                    offset_ss = (ss_reg[target_e] - target_ss) * offset_ss_step;
                }
                var coords_x_target_ref = coords_x_target_e;
                coords_x_target_e += direction * offset_ss;


                // Now, we mark all possible splice sites, either if they have junctions jumping in/or out of them or not
                // splice sites dashed lines in LSV exon
                if (direction > 0 && ss !== num_ss || direction < 0 && ss !== 1) {

                    if (!rendered_exons[lsvs_fields[1][0]]) {  // Render ssites only if they haven't been rendered
                        ctx.strokeStyle = "rgba(0, 0, 0, 0.6)";
                        var offset_ss_aux = null,
                            coords_x_target_ss = null;
                        for (var ii = 1; ii < target_num_ss; ii++) {
                            coords_x_target_ss = coords_x_target_ref;
                            if (direction > 0) {
                                offset_ss_aux = ii * offset_ss_step;
                            } else {
                                offset_ss_aux = (ss_reg[target_e] - ii) * offset_ss_step;
                            }

                            if (offset_ss_aux) {
                                coords_x_target_ss += direction * offset_ss_aux;
                                drawDashedLine(ctx, Math.round(coords_x_target_ss), Math.round(coords[1]), Math.round(coords_x_target_ss), Math.round(coords[1] + exon_height), 2);
                            }
                        }

                        rendered_exons[lsvs_fields[1][0]] = 1;
                    }
                }

                var mid_x = (coords_x_start_e + coords_x_target_e) / 2;
                ctx.lineWidth = 2;
                ctx.strokeStyle = getColor(Math.max(0, n_ways - 1 - count_starts_ends), BREWER_PALETTE, 0.9);

                // splice sites dashed lines
                if (direction > 0 && ss !== num_ss || direction < 0 && ss !== 1) {
                    // Check if is a special exon (starter or finisher)
                    drawDashedLine(ctx, Math.round(coords_x_start_e), Math.round(coords[1]), Math.round(coords_x_start_e), Math.round(coords[1] + exon_height), 2);
                }

                // splice sites dashed lines in target exon
                if (target_ss !== 1 && direction > 0 || target_ss !== ss_reg[target_e] && direction < 0) {
                    drawDashedLine(ctx, Math.round(coords_x_target_e), Math.round(coords[1]), Math.round(coords_x_target_e), Math.round(coords[1] + exon_height), 2);
                }

                var junc_h_pos = Math.round(margins[3] * 8 * ss); // Math.round((1 - (Math.abs(coords_x_start_e - coords_x_target_e)/canvas.width)) * (canvas.height*(1-percentage_exon)));
                // junctions lines
                drawLine(ctx, Math.round(coords_x_start_e), Math.round(coords[1]), Math.round(mid_x), junc_h_pos);
                drawLine(ctx, Math.round(mid_x), junc_h_pos, Math.round(coords_x_target_e), Math.round(coords[1]));

                // render special marker for exon alternative start/end
                if (lsvs_fields[1].indexOf('.') === -1) {
                    ctx.strokeStyle = "rgba(0, 0, 0, 0.6)";
                    drawLine(ctx, Math.round(coords_x_start_e), Math.round(coords[1]), Math.round(coords_x_start_e), Math.round(coords[1] + exon_height));
                    drawArrow(ctx, Math.round(coords_x_start_e + direction * Math.max(10, percentage_exon / 2 * exon_width)), Math.round(coords[1] + exon_height / 2), Math.round(coords_x_start_e + direction * 2), Math.round(coords[1] + exon_height / 2), Math.max(5, Math.round((percentage_exon / 2 * exon_width) / 2)));

                }

            }


            ctx.lineWidth = previous_ctx[0];
            ctx.strokeStyle = previous_ctx[1];
        }
    })


};

Lsv.prototype.render_exon = function (canvas, exon, pixel_factor, margin, percen_exon, counter_exon) {

    var exon_height = percen_exon * canvas.height; // canvas.height-2*margin[2];

    var ctx = canvas.getContext("2d");
    if (exon.type === 1) {
        ctx.strokeStyle = "rgba(255, 165, 0, 1)"; //getColor(2, BREWER_PALETTE, .8);
        ctx.fillStyle = "rgba(255, 165, 0, 0.2)"; //getColor(2, BREWER_PALETTE, .2);
    } else {
        ctx.strokeStyle = "rgba(0, 0, 0, 0.8)";
        ctx.fillStyle = "rgba(0, 0, 0, 0.2)";
    }

    if (exon.type === 2) {
        drawDashedRectangle(ctx,
            margin[0] + Math.round((exon.coords[0]) * pixel_factor),
            canvas.height - margin[3],
            Math.round((exon.coords[1] - exon.coords[0]) * pixel_factor),
            exon_height, 4);
    } else {
        ctx.strokeRect(margin[0] + Math.round((exon.coords[0]) * pixel_factor),
            canvas.height - margin[3],
            Math.round((exon.coords[1] - exon.coords[0]) * pixel_factor),
            -exon_height);

        ctx.fillRect(
            margin[0] + Math.round((exon.coords[0]) * pixel_factor),
            canvas.height - margin[3],
            Math.round((exon.coords[1] - exon.coords[0]) * pixel_factor),
            -exon_height
        );
    }

    // Draw exon name (number)
    ctx.textAlign = "center";
    ctx.font = "11px Arial";
    ctx.strokeStyle = "rgba(0, 0, 0, 1)";
    ctx.fillStyle = "rgba(0, 0, 0, 1)";
    ctx.fillText(counter_exon,
        margin[0] + Math.round((exon.coords[0]) * pixel_factor) + Math.round((exon.coords[1] - exon.coords[0]) * pixel_factor / 2),
        canvas.height - margin[3] - exon_height / 4
    );
    ctx.strokeText(counter_exon,
        margin[0] + Math.round((exon.coords[0]) * pixel_factor) + Math.round((exon.coords[1] - exon.coords[0]) * pixel_factor / 2),
        canvas.height - margin[3] - exon_height / 4
    );

    function render_extra_exon(canvas, coords_extra, pixel_factor, margin, exon_height) {
        var ctx = canvas.getContext("2d");
        ctx.strokeStyle = getColor(2, BREWER_PALETTE, .8);
        ctx.fillStyle = getColor(2, BREWER_PALETTE, .2);
        drawDashedRectangle(ctx,
            margin[0] + Math.round((coords_extra[0]) * pixel_factor),
            canvas.height - margin[3],
            Math.round((coords_extra[1] - coords_extra[0]) * pixel_factor),
            exon_height, 2
        );

        ctx.fillRect(
            margin[0] + Math.round((coords_extra[0]) * pixel_factor),
            canvas.height - margin[3],
            Math.round((coords_extra[1] - coords_extra[0]) * pixel_factor),
            -exon_height
        );

    }

    if (exon.coords_extra) {
        for (var i = 0; i < exon.coords_extra.length; i++) {
            render_extra_exon(canvas, exon.coords_extra[i], pixel_factor, margin, exon_height);
        }
    }

};

Lsv.prototype.render_intron_retained = function (canvas, intron, pixel_factor, margin, percen_exon, junc_count) {
    var intron_height = Math.round(percen_exon * 2 / 4 * canvas.height);

    var ctx = canvas.getContext("2d");
    ctx.strokeStyle = getColor(junc_count, BREWER_PALETTE, 1);
    ctx.fillStyle = getColor(junc_count, BREWER_PALETTE, .8);

    ctx.fillRect(
        margin[0] + Math.round((intron.coords[0]) * pixel_factor),
        canvas.height * (1 - percen_exon * 1 / 4) - margin[3],
        Math.round((intron.coords[1] - intron.coords[0]) * pixel_factor),
        -intron_height
    );
};

Lsv.prototype.drawLSVCompactStackBars = function (canvas, fillMode) {

    function createGradientLSVGroupsCompact(coords, group, count, fillMode, hue) {
        //create a gradient object from the canvas context

        var gradient = ctx.createLinearGradient(coords.x1, coords.y1, coords.x2, coords.y2);
        // Add the colors with fixed stops.
        if (fillMode < 2 || fillMode === 4) {
            gradient.addColorStop(0, getColor(count, BREWER_PALETTE, 1));  // No gradient
            gradient.addColorStop(1, getColor(count, BREWER_PALETTE, 1));

        } else if (fillMode >= 2) {
            if (hue) {
                gradient.addColorStop(0, getColor(count, BREWER_PALETTE, hue));  // No gradient
                gradient.addColorStop(1, getColor(count, BREWER_PALETTE, hue));
            } else {
                var newHue = 1;
                if (fillMode === 2) {
                    newHue = 1 - (group.quartiles[count][4] - group.quartiles[count][0]);
                } else {
                    newHue = 1 - group.variances[count];
                }
                gradient.addColorStop(0, getColor(count, BREWER_PALETTE, newHue));  // No gradient
                gradient.addColorStop(1, getColor(count, BREWER_PALETTE, newHue));
            }
        }
        return gradient;
    }

    var fillingPSIArea = function (ctx, fillMode, group, lsv_count, x1, y1, x2, y2) {
        var area = [];
        if (fillMode === 0) {
            // Fill from center to both left and right
            // First approach, use 1 - the 25 to 75 percentile to fill the area
            area[0] = Math.round((1 - (group.quartiles[lsv_count][4] - group.quartiles[lsv_count][0])) * (x2 - x1));
            area[1] = y2 - y1;
            ctx.strokeStyle = ctx.fillStyle;
            drawRectangle(ctx, x1 + (x2 - x1 - area[0]) / 2, y1, area[0], area[1], true);
//            ctx.strokeRect(x1+(x2-x1-area[0])/2, y1, area[0], area[1]);
        } else if (fillMode === 1) {
            // Fill from center to both left and right
            // Second approach, use 1 - the variance
            //area[0] = (1 - group.variances[lsv_count]) * (x2 - x1);
            area[0] = x2 - x1;
            area[1] = y2 - y1;
            ctx.strokeStyle = ctx.fillStyle;
            drawRectangle(ctx, Math.round(x1 + (x2 - x1 - area[0]) / 2), Math.round(y1), Math.floor(area[0]), Math.floor(area[1]), true);
//            ctx.strokeRect(x1+(x2-x1-area[0])/2, y1, area[0], area[1]);
        } else if (fillMode === 2 || fillMode === 3) {
            // Fill all area, use the hue to represent variance
            area[0] = x2 - x1;
            area[1] = y2 - y1;
            ctx.strokeStyle = ctx.fillStyle;
            drawRectangle(ctx, x1, y1, area[0], area[1], true);
//            ctx.strokeRect(x1, y1, area[0], area[1]);
        } else if (fillMode === 4) {
            // Fill from left to right
            // Using 1 - the 25 to 75 percentile to fill the area
            area[0] = Math.round((x2 - x1)); //(1 - (group.quartiles[lsv_count][4] - group.quartiles[lsv_count][0])) *
            area[1] = y2 - y1;
            ctx.strokeStyle = ctx.fillStyle;
            drawRectangle(ctx, x1, y1, area[0], area[1], true);
        }
    };

    if (canvas.getContext) {  // check for support
        var ctx = canvas.getContext("2d");
        var lsv_id = canvas.getAttribute('data-lsv-id');
        var group_name = canvas.getAttribute('data-group');

        // Clear canvas
        ctx.clearRect(0, 0, canvas.width, canvas.height);


        this.db.get(lsv_id).then(function (data) {
            // var groups = JSON.parse(groups_str.replace(/'/g, '"'));
            var groups = [data];
            // Calculate origins_coords
            var header_height = 0; // canvas.height*.1;
            var num_groups = groups.length;

            var sub_canvas_w = canvas.width / num_groups;
            var sub_canvas_h = canvas.height - header_height;
            // var sub_canvas_margins = [sub_canvas_w * .00, 0, header_height, sub_canvas_h * .21];
            var sub_canvas_margins = [0, 0, 0, 0];
            var sub_canvas_pixels = [sub_canvas_w - sub_canvas_margins[0] - sub_canvas_margins[1], canvas.height - sub_canvas_margins[3] - sub_canvas_margins[2]];

            // Draw sub-boxes separators and headers
            ctx.textAlign = "center";
            ctx.font = "8pt Arial";

            var tmpStyle = ctx.strokeStyle;
            ctx.strokeStyle = "rgba(0, 0, 0, .5)";
            var i = 0;
            for (var ii = 1; ii < num_groups; ii++) {
                if (i) {
                    drawLine(ctx, sub_canvas_w * i, 0, sub_canvas_w * i, canvas.height);
                }
                i++;
                ctx.fillText(groups[ii].name, sub_canvas_w * (i - 1 / 2), sub_canvas_h - sub_canvas_margins[3] + sub_canvas_h * .2 - 1);
            }
            ctx.strokeStyle = tmpStyle;

            var origins_coords = [];
            for (var count_groups = 0; count_groups < num_groups; count_groups++) {
                origins_coords[count_groups] = [sub_canvas_margins[0] + count_groups * sub_canvas_w, canvas.height - sub_canvas_margins[3]];
            }


            // Separators
            for (var count = 0; count < num_groups; count++) {
                var offset = 0;
                var acc_height = 0;
                var group = groups[count];

                for (var lsv_count = 0; lsv_count < group.group_bins[group_name].length; lsv_count++) {
                    // Calculate the height of the accumulated mean
                    acc_height += group.group_means_rounded[group_name][lsv_count];

                    var coords_gradient = {
                        'x1': origins_coords[count][0] + offset,
                        'y1': origins_coords[count][1],
                        'x2': origins_coords[count][0] + offset + sub_canvas_pixels[0],
                        'y2': origins_coords[count][1] - sub_canvas_pixels[1]
                    };
                    ctx.fillStyle = createGradientLSVGroupsCompact(coords_gradient, group, lsv_count, fillMode, 1);
                    ctx.strokeStyle = createGradientLSVGroupsCompact(coords_gradient, group, lsv_count, fillMode, 1);

                    // Filling the PSI area
                    fillingPSIArea(ctx,
                        fillMode,
                        group,
                        lsv_count,
                        origins_coords[count][0],
                        origins_coords[count][1] - (acc_height - group.group_means_rounded[group_name][lsv_count]) * sub_canvas_pixels[1],
                        origins_coords[count][0] + sub_canvas_pixels[0],
                        origins_coords[count][1] - (acc_height) * sub_canvas_pixels[1]
                    );
                }
            }
        });
    }
};


Lsv.prototype.drawDeltaLSVCompactSVG = function (el, lsv) {
    var width = 200;
    var height = 20;
    var margin = {top: 1, bottom: 8, left: 2, right: 2};
    var border_frame = 2;
    var MIN_DELTAPSI = .05;

    var svgContainer = d3.select(el)
        .append("svg")
        .attr("class", "excl-incl-rect")
        .attr("width", width)
        .attr("height", height);

    var markerWidth = 6;
    var markerHeight = 6;
    var cRadius = 30; // play with the cRadius value
    var refX = cRadius + (markerWidth * 2);
    var refY = -Math.sqrt(cRadius);

    // Define arrow markers, right
    svgContainer.append("defs")
        .append("marker")
        .attr("id", "arrowhead-right")
        .attr("viewBox", "0 -5 5 10")
        .attr("refX", 5)
        .attr("refY", 0)
        .attr("markerWidth", 4)
        .attr("markerHeight", 3)
        .attr("orient", "auto")
        .append("path")
        .attr("d", "M0,-5L10,0L0,5");

    // Define arrow markers, left
    svgContainer.append("defs")
        .append("marker")
        .attr("id", "arrowhead-left")
        .attr("viewBox", "-5 -5 5 10")
        .attr("refX", -5)
        .attr("refY", 0)
        .attr("markerWidth", 4)
        .attr("markerHeight", 3)
        .attr("orient", "auto")
        .append("path")
        .attr("d", "M0,-5L-10,0L0,5");


    svgContainer.append("line")
        .attr("x1", margin.left)
        .attr("y1", Math.round((height - margin.bottom) / 2))
        .attr("x2", width - margin.right - margin.left)
        .attr("y2", Math.round((height - margin.bottom) / 2))
        .style('stroke', 'black')
        .style('stroke-width', border_frame)
        //        .attr("stroke-opacity", .5)
        .style('fill', 'none')
        .attr("marker-end", "url(#arrowhead-right)")
        .attr("marker-start", "url(#arrowhead-left)");


    // Draw x-axis ticks
    svgContainer.append("text")
        .attr("x", width / 2)
        .attr("y", height)
        .attr("text-anchor", "middle")
        .attr("font-size", "8px")
        .attr("fill", "black")
        .text("0");

    // Draw excl-incl bars
    var last_excl_pos = width / 2,
        last_incl_pos = width / 2;
    for (var ii = 0; ii < lsv.excl_incl.length; ii++) {
        svgContainer.append("rect")
            .attr("x", last_excl_pos - Math.round((width / 2 - margin.left) * lsv.excl_incl[ii][0]))
            .attr("y", margin.top)
            .attr("width", Math.round((width / 2 - margin.left) * lsv.excl_incl[ii][0]))
            .attr("height", height - margin.bottom - margin.top)
            .style('fill', getColor(ii, BREWER_PALETTE, 1));
        svgContainer.append("rect")
            .attr("x", last_incl_pos)
            .attr("y", margin.top)
            .attr("width", Math.round((width / 2 - margin.right) * lsv.excl_incl[ii][1]))
            .attr("height", height - margin.bottom - margin.top)
            .style('fill', getColor(ii, BREWER_PALETTE, 1));

        if (lsv.excl_incl[ii][0] < MIN_DELTAPSI && lsv.excl_incl[ii][1] < MIN_DELTAPSI)
            continue;

        // Draw percentages text
        if (Math.round((width / 2 - margin.left) * lsv.excl_incl[ii][0]) >= 1) {
            last_excl_pos -= Math.round((width / 2 - margin.left) * lsv.excl_incl[ii][0]);
            svgContainer.append("text")
                .attr("x", last_excl_pos)
                .attr("y", height)
                .attr("text-anchor", "middle")
                .attr("font-size", "9px")
                .attr("fill", getColor(ii, BREWER_PALETTE, 1))
                .text(Math.round((width / 2 - margin.left) * lsv.excl_incl[ii][0]));

        }

        if (Math.round((width / 2 - margin.right) * lsv.excl_incl[ii][1]) >= 1) {
            last_incl_pos += Math.round((width / 2 - margin.right) * lsv.excl_incl[ii][1]);
            svgContainer.append("text")
                .attr("x", last_incl_pos)
                .attr("y", height)
                .attr("text-anchor", "middle")
                .attr("font-size", "9px")
                .attr("fill", getColor(ii, BREWER_PALETTE, 1))
                .text(Math.round((width / 2 - margin.right) * lsv.excl_incl[ii][1]));
        }
    }

    // Draw separator
    svgContainer.append("line")
        .attr("x1", width / 2)
        .attr("y1", 0)
        .attr("x2", width / 2)
        .attr("y2", height - margin.bottom + 2)
        .attr("stroke-width", 2)
        .attr("stroke-opacity", .8)
        .attr("stroke", "black");


    return svgContainer;

};