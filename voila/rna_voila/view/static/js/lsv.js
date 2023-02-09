// Adding dash lines to the Canvas Rendering
CanvasRenderingContext2D.prototype.dashedLine = function (x1, y1, x2, y2, dashLen) {

    if (dashLen === undefined) dashLen = 2;
    this.moveTo(x1, y1);

    const dX = x2 - x1;
    const dY = y2 - y1;
    const dashes = Math.floor(Math.sqrt(dX * dX + dY * dY) / dashLen);
    const dashX = dX / dashes;
    const dashY = dY / dashes;

    let q = 0;
    while (q++ < dashes) {
        x1 += dashX;
        y1 += dashY;
        this[q % 2 === 0 ? 'moveTo' : 'lineTo'](x1, y1);
    }
    this[q % 2 === 0 ? 'moveTo' : 'lineTo'](x2, y2);
};

// Get Color from Brewer Palette
const BREWER_PALETTE = [
    [228, 26, 28],
    [55, 126, 184],
    [77, 175, 74],
    [152, 78, 163],
    [255, 127, 0],
//    [255,255,51],
    [166, 86, 40],
    [247, 129, 191],
    [153, 153, 153],

    [28, 126, 128],
    [155, 226, 29],
    [177, 275, 19],
    [252, 178, 8],
    [55, 227, 100],
//    [55,55,151],
    [11, 186, 140],
    [47, 229, 36],
    [253, 253, 253]
];


class Lsv {
    constructor(lsv_data) {
        if (lsv_data) {
            this.exon_number = lsv_data.exon_number;
            this.lsv = lsv_data.lsv
        }
    }

    static draw_dashed_line(contextO, startx, starty, endx, endy, dashLen) {
        contextO.beginPath();
        contextO.dashedLine(startx, starty, endx, endy, dashLen);
        contextO.stroke();
    }

    static get_color(colorNumber, palette, hue) {
        colorNumber = colorNumber % 16;
        return "rgba(" + palette[colorNumber].toString() + ", " + hue + ")";
        // return rgbToHex(palette[colorNumber][0], palette[colorNumber][1], palette[colorNumber][2]);
    }

    static draw_line(contextO, startx, starty, endx, endy) {
        contextO.beginPath();
        contextO.moveTo(startx, starty);
        contextO.lineTo(endx, endy);
        contextO.closePath();
        contextO.stroke();
    }

    async cartoon(canvas) {
        canvas.setAttribute('height', 80);
        canvas.setAttribute('width', 200);
        this.render_lsv_splice_graph(canvas)
    }

    render_lsv_splice_graph(canvas) {
        if (canvas.getContext) {
            // Render LSV representation from a string text representing the LSV i.e.: s|1e1.3o3|2e1.2o3|3e1.2o3|i|4e1.1o3  NEW: Intron Retention (i field)
            let ctx = canvas.getContext("2d");

            const margins = [1, 1, 1, 1],
                pixel_factor = 1,
                percentage_exon = .2;
            let percentage_intron = .15;

            const exon_lsv_number = this.exon_number;

            const lsv_data = canvas.getAttribute('data-lsv-type');
            const lsvs = lsv_data.split('|');
            const ir_marker = 'i';

            let intron_ret_i = lsvs.indexOf(ir_marker);
            if (intron_ret_i > -1) {
                lsvs.splice(lsvs.indexOf(ir_marker), 1);  // Modifies the array in place
                let num_alt_start_end = 0;
                for (let ff = 1; ff < lsvs.length; ff++) {
                    if (lsvs[ff].indexOf('.') === -1) {
                        num_alt_start_end++;
                    }
                }
                intron_ret_i -= num_alt_start_end;
            }

            // Num exons_obj
            let num_exons = 0;
            let num_ss = 0;
            const ss_reg = {};

            for (let n_ways = 1; n_ways < lsvs.length; n_ways++) {
                const lsvs_fields = lsvs[n_ways].split('e');
                num_exons = Math.max(num_exons, parseInt(lsvs_fields[1][0]));
                num_ss = Math.max(num_ss, parseInt(lsvs_fields[0]));

                if (lsvs_fields[1] === 0) continue;
                const exonNum_ss = lsvs_fields[1].split('.');

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

            const sourceOrTarget = lsvs[0];
            const area = [canvas.width - margins[0] - margins[1], canvas.height - margins[2] - margins[3]];

            percentage_intron = Math.min(percentage_intron, .4 / (num_exons - 1));
            const exon_width = (area[0] - (num_exons - 1) * area[0] * percentage_intron - margins[0] - margins[1]) / num_exons;
            let start = margins[0];
            const direction = sourceOrTarget === 's' ? 1 : -1;

            // Render exons_obj
            const exons = [];
            for (let i = 0; i < num_exons; i++) {
                const exon = {
                    'coords': [Math.round(start), Math.round(start + exon_width)],
                    'type': (direction > 0 && i === 0 || direction < 0 && i === num_exons - 1 ? 1 : 0)
                };
                exons.push(exon);
                let number_exon = (direction > 0 ? i : i + 1);
                number_exon = (number_exon === 0 || number_exon === num_exons ? exon_lsv_number : '');
                this.render_exon(canvas, exon, pixel_factor, margins, percentage_exon, number_exon);
                start += exon_width + percentage_intron * area[0];
            }

            // Render IR
            if (intron_ret_i > -1) {
                const intron = {
                    'coords': (direction > 0
                        ? [margins[0] + exon_width, margins[0] + exon_width + percentage_intron * area[0]]
                        : [margins[0] + (num_exons - 1) * exon_width + (num_exons - 2) * percentage_intron * area[0],
                            margins[0] + (num_exons - 1) * exon_width + (num_exons - 1) * percentage_intron * area[0]])
                };
                this.render_intron_retained(canvas, intron, pixel_factor, margins, percentage_exon, intron_ret_i - 1);
            }

            // Render junctions
            let count_starts_ends = 0;
            const index_first = direction > 0 ? 0 : exons.length - 1;
            const coords = exons[index_first].coords;
            const exon_height = percentage_exon * canvas.height;

            // If source, we move the coords to the first splice site
            if (direction > 0) {
                coords[0] += exon_width / 2 + direction * (exon_width / 2) / num_ss;
            }
            coords[1] = canvas.height - exon_height - margins[3];

            // For rendering all splice sites even if they don't have a junction jumping in or out of them,
            // keep a registry of the exons whose ssites have been already rendered
            const rendered_exons = {};
            const previous_ctx = [ctx.lineWidth, ctx.strokeStyle];
            for (let n_ways = 1; n_ways < lsvs.length; n_ways++) {
                const lsvs_fields = lsvs[n_ways].split('e');
                const ss = parseInt(lsvs_fields[0]),
                    target_e = lsvs_fields[1].split('.')[0],
                    target_ss = parseInt(lsvs_fields[1].split('.')[1]);
                let target_num_ss = null;
                if (lsvs_fields[1].indexOf('o') > 0) {
                    target_num_ss = parseInt(lsvs_fields[1].split('.')[1].split('o')[1]);
                }

                const coords_x_start_e = coords[0] + ((exon_width / 2) / num_ss) * (ss - 1);
                let coords_x_target_e = null;

                if (lsvs_fields[1] === 0) {
                    coords_x_target_e = coords_x_start_e;
                    count_starts_ends++;
                } else {
                    const target_exon = (direction > 0 ? exons[lsvs_fields[1][0]] : exons[lsvs_fields[1][0] - 1]);
                    coords_x_target_e = target_exon.coords[direction > 0 ? 0 : 1];
                }

                const offset_ss_step = (exon_width * 2 / 3) / ss_reg[target_e];
                let offset_ss = 0;
                if (direction > 0) {
                    offset_ss = (target_ss - 1) * offset_ss_step;
                } else {
                    offset_ss = (ss_reg[target_e] - target_ss) * offset_ss_step;
                }
                const coords_x_target_ref = coords_x_target_e;
                coords_x_target_e += direction * offset_ss;


                // Now, we mark all possible splice sites, either if they have junctions jumping in/or out of them or not
                // splice sites dashed lines in LSV exon
                if (direction > 0 && ss !== num_ss || direction < 0 && ss !== 1) {

                    if (!rendered_exons[lsvs_fields[1][0]]) {  // Render ssites only if they haven't been rendered
                        ctx.strokeStyle = "rgba(0, 0, 0, 0.6)";
                        let offset_ss_aux = null;
                        let coords_x_target_ss = null;
                        for (let ii = 1; ii < target_num_ss; ii++) {
                            coords_x_target_ss = coords_x_target_ref;
                            if (direction > 0) {
                                offset_ss_aux = ii * offset_ss_step;
                            } else {
                                offset_ss_aux = (ss_reg[target_e] - ii) * offset_ss_step;
                            }

                            if (offset_ss_aux) {
                                coords_x_target_ss += direction * offset_ss_aux;
                                Lsv.draw_dashed_line(ctx, Math.round(coords_x_target_ss), Math.round(coords[1]), Math.round(coords_x_target_ss), Math.round(coords[1] + exon_height), 2);
                            }
                        }

                        rendered_exons[lsvs_fields[1][0]] = 1;
                    }
                }

                const mid_x = (coords_x_start_e + coords_x_target_e) / 2;
                ctx.lineWidth = 2;
                ctx.strokeStyle = Lsv.get_color(Math.max(0, n_ways - 1 - count_starts_ends), BREWER_PALETTE, 0.9);

                // splice sites dashed lines
                if (direction > 0 && ss !== num_ss || direction < 0 && ss !== 1) {
                    // Check if is a special exon (starter or finisher)
                    Lsv.draw_dashed_line(ctx, Math.round(coords_x_start_e), Math.round(coords[1]), Math.round(coords_x_start_e), Math.round(coords[1] + exon_height), 2);
                }

                // splice sites dashed lines in target exon
                if (target_ss !== 1 && direction > 0 || target_ss !== ss_reg[target_e] && direction < 0) {
                    Lsv.draw_dashed_line(ctx, Math.round(coords_x_target_e), Math.round(coords[1]), Math.round(coords_x_target_e), Math.round(coords[1] + exon_height), 2);
                }

                const junc_h_pos = Math.round(margins[3] * 8 * ss); // Math.round((1 - (Math.abs(coords_x_start_e - coords_x_target_e)/canvas.group_width)) * (canvas.group_height*(1-percentage_exon)));
                // junctions lines
                Lsv.draw_line(ctx, Math.round(coords_x_start_e), Math.round(coords[1]), Math.round(mid_x), junc_h_pos);
                Lsv.draw_line(ctx, Math.round(mid_x), junc_h_pos, Math.round(coords_x_target_e), Math.round(coords[1]));

                // render special marker for exon alternative start/end
                if (lsvs_fields[1].indexOf('.') === -1) {
                    ctx.strokeStyle = "rgba(0, 0, 0, 0.6)";
                    Lsv.draw_line(ctx, Math.round(coords_x_start_e), Math.round(coords[1]), Math.round(coords_x_start_e), Math.round(coords[1] + exon_height));
                    drawArrow(ctx, Math.round(coords_x_start_e + direction * Math.max(10, percentage_exon / 2 * exon_width)), Math.round(coords[1] + exon_height / 2), Math.round(coords_x_start_e + direction * 2), Math.round(coords[1] + exon_height / 2), Math.max(5, Math.round((percentage_exon / 2 * exon_width) / 2)));

                }

            }
            ctx.lineWidth = previous_ctx[0];
            ctx.strokeStyle = previous_ctx[1];
        }

    }

    render_exon(canvas, exon, pixel_factor, margin, percen_exon, counter_exon) {

        const exon_height = percen_exon * canvas.height; // canvas.group_height-2*margin[2];

        const ctx = canvas.getContext("2d");
        if (exon.type === 1) {
            ctx.strokeStyle = "rgba(255, 165, 0, 1)"; //Lsv.get_color(2, BREWER_PALETTE, .8);
            ctx.fillStyle = "rgba(255, 165, 0, 0.2)"; //Lsv.get_color(2, BREWER_PALETTE, .2);
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
            const ctx = canvas.getContext("2d");
            ctx.strokeStyle = Lsv.get_color(2, BREWER_PALETTE, .8);
            ctx.fillStyle = Lsv.get_color(2, BREWER_PALETTE, .2);
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
            for (let i = 0; i < exon.coords_extra.length; i++) {
                render_extra_exon(canvas, exon.coords_extra[i], pixel_factor, margin, exon_height);
            }
        }
    }

    render_intron_retained(canvas, intron, pixel_factor, margin, percen_exon, junc_count) {
        const intron_height = Math.round(percen_exon * 2 / 4 * canvas.height);

        const ctx = canvas.getContext("2d");
        ctx.strokeStyle = Lsv.get_color(junc_count, BREWER_PALETTE, 1);
        ctx.fillStyle = Lsv.get_color(junc_count, BREWER_PALETTE, .8);

        ctx.fillRect(
            margin[0] + Math.round((intron.coords[0]) * pixel_factor),
            canvas.height * (1 - percen_exon * 1 / 4) - margin[3],
            Math.round((intron.coords[1] - intron.coords[0]) * pixel_factor),
            -intron_height
        );
    }

    draw_lsv_compact_stack_bars(canvas) {

        const fillMode = 1;
        const createGradientLSVGroupsCompact = (coords, group, count, fillMode, hue, ctx) => {
            //create a gradient object from the canvas context

            const gradient = ctx.createLinearGradient(coords.x1, coords.y1, coords.x2, coords.y2);
            // Add the colors with fixed stops.
            if (fillMode < 2 || fillMode === 4) {
                gradient.addColorStop(0, Lsv.get_color(count, BREWER_PALETTE, 1));  // No gradient
                gradient.addColorStop(1, Lsv.get_color(count, BREWER_PALETTE, 1));

            } else if (fillMode >= 2) {
                if (hue) {
                    gradient.addColorStop(0, Lsv.get_color(count, BREWER_PALETTE, hue));  // No gradient
                    gradient.addColorStop(1, Lsv.get_color(count, BREWER_PALETTE, hue));
                } else {
                    let newHue = 1;
                    if (fillMode === 2) {
                        newHue = 1 - (group.quartiles[count][4] - group.quartiles[count][0]);
                    } else {
                        newHue = 1 - group.constiances[count];
                    }
                    gradient.addColorStop(0, Lsv.get_color(count, BREWER_PALETTE, newHue));  // No gradient
                    gradient.addColorStop(1, Lsv.get_color(count, BREWER_PALETTE, newHue));
                }
            }
            return gradient;
        };

        const fillingPSIArea = (ctx, fillMode, group, lsv_count, x1, y1, x2, y2) => {
            const area = [];
            if (fillMode === 0) {
                // Fill from center to both left and right
                // First approach, use 1 - the 25 to 75 percentile to fill the area
                area[0] = Math.round((1 - (group.quartiles[lsv_count][4] - group.quartiles[lsv_count][0])) * (x2 - x1));
                area[1] = y2 - y1;
                ctx.strokeStyle = ctx.fillStyle;
                this.draw_rectangle(ctx, x1 + (x2 - x1 - area[0]) / 2, y1, area[0], area[1], true);
//            ctx.strokeRect(x1+(x2-x1-area[0])/2, y1, area[0], area[1]);
            } else if (fillMode === 1) {
                // Fill from center to both left and right
                // Second approach, use 1 - the constiance
                //area[0] = (1 - group.constiances[lsv_count]) * (x2 - x1);
                area[0] = x2 - x1;
                area[1] = y2 - y1;
                ctx.strokeStyle = ctx.fillStyle;
                this.draw_rectangle(ctx, Math.round(x1 + (x2 - x1 - area[0]) / 2), Math.round(y1), Math.floor(area[0]), Math.floor(area[1]), true);
//            ctx.strokeRect(x1+(x2-x1-area[0])/2, y1, area[0], area[1]);
            } else if (fillMode === 2 || fillMode === 3) {
                // Fill all area, use the hue to represent constiance
                area[0] = x2 - x1;
                area[1] = y2 - y1;
                ctx.strokeStyle = ctx.fillStyle;
                this.draw_rectangle(ctx, x1, y1, area[0], area[1], true);
//            ctx.strokeRect(x1, y1, area[0], area[1]);
            } else if (fillMode === 4) {
                // Fill from left to right
                // Using 1 - the 25 to 75 percentile to fill the area
                area[0] = Math.round((x2 - x1)); //(1 - (group.quartiles[lsv_count][4] - group.quartiles[lsv_count][0])) *
                area[1] = y2 - y1;
                ctx.strokeStyle = ctx.fillStyle;
                this.draw_rectangle(ctx, x1, y1, area[0], area[1], true);
            }
        };

        if (canvas.getContext) {  // check for support
            const ctx = canvas.getContext("2d");
            const group_name = canvas.dataset.group;


            const groups = [this.lsv];
            // Calculate origins_coords
            const header_height = 0; // canvas.group_height*.1;
            const num_groups = groups.length;

            const sub_canvas_w = canvas.width / num_groups;
            const sub_canvas_h = canvas.height - header_height;
            // const sub_canvas_margins = [sub_canvas_w * .00, 0, header_height, sub_canvas_h * .21];
            const sub_canvas_margins = [0, 0, 0, 0];
            const sub_canvas_pixels = [sub_canvas_w - sub_canvas_margins[0] - sub_canvas_margins[1], canvas.height - sub_canvas_margins[3] - sub_canvas_margins[2]];

            // Draw sub-boxes separators and headers
            ctx.textAlign = "center";
            ctx.font = "8pt Arial";

            const tmpStyle = ctx.strokeStyle;
            ctx.strokeStyle = "rgba(0, 0, 0, .5)";
            let i = 0;
            for (let ii = 1; ii < num_groups; ii++) {
                if (i) {
                    this.drawLine(ctx, sub_canvas_w * i, 0, sub_canvas_w * i, canvas.height);
                }
                i++;
                ctx.fillText(groups[ii].name, sub_canvas_w * (i - 1 / 2), sub_canvas_h - sub_canvas_margins[3] + sub_canvas_h * .2 - 1);
            }
            ctx.strokeStyle = tmpStyle;

            const origins_coords = [];
            for (let count_groups = 0; count_groups < num_groups; count_groups++) {
                origins_coords[count_groups] = [sub_canvas_margins[0] + count_groups * sub_canvas_w, canvas.height - sub_canvas_margins[3]];
            }


            // Separators
            for (let count = 0; count < num_groups; count++) {
                const offset = 0;
                let acc_height = 0;
                const group = groups[count];

                for (let lsv_count = 0; lsv_count < group.group_bins[group_name].length; lsv_count++) {
                    // Calculate the group_height of the accumulated mean
                    acc_height += group.group_means[group_name][lsv_count];

                    const coords_gradient = {
                        'x1': origins_coords[count][0] + offset,
                        'y1': origins_coords[count][1],
                        'x2': origins_coords[count][0] + offset + sub_canvas_pixels[0],
                        'y2': origins_coords[count][1] - sub_canvas_pixels[1]
                    };
                    ctx.fillStyle = createGradientLSVGroupsCompact(coords_gradient, group, lsv_count, fillMode, 1, ctx);
                    ctx.strokeStyle = createGradientLSVGroupsCompact(coords_gradient, group, lsv_count, fillMode, 1, ctx);

                    // Filling the PSI area
                    fillingPSIArea(ctx,
                        fillMode,
                        group,
                        lsv_count,
                        origins_coords[count][0],
                        origins_coords[count][1] - (acc_height - group.group_means[group_name][lsv_count]) * sub_canvas_pixels[1],
                        origins_coords[count][0] + sub_canvas_pixels[0],
                        origins_coords[count][1] - (acc_height) * sub_canvas_pixels[1]
                    );
                }
            }
        }

    }

    draw_rectangle(contextO, x, y, w, h, fill) {
        contextO.beginPath();
        contextO.rect(x, y, w, h);
        contextO.closePath();
        contextO.stroke();
        if (fill) contextO.fill();
    }

    draw_delta_lsv_compact_svg(svg) {
        const lsv = this.lsv;
        const width = 200;
        const height = 20;
        const margin = {top: 1, bottom: 8, left: 2, right: 2};
        const border_frame = 2;
        const MIN_DELTAPSI = .05;

        const svgContainer = d3.select(svg)
            .attr("width", width)
            .attr("height", height);

        const markerWidth = 6;
        const markerHeight = 6;
        const cRadius = 30; // play with the cRadius value
        const refX = cRadius + (markerWidth * 2);
        const refY = -Math.sqrt(cRadius);

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
            .style('stroke-group_width', border_frame)
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
        let last_excl_pos = width / 2,
            last_incl_pos = width / 2;
        for (let ii = 0; ii < lsv.excl_incl.length; ii++) {
            svgContainer.append("rect")
                .attr("x", last_excl_pos - Math.round((width / 2 - margin.left) * lsv.excl_incl[ii][0]))
                .attr("y", margin.top)
                .attr("width", Math.round((width / 2 - margin.left) * lsv.excl_incl[ii][0]))
                .attr("height", height - margin.bottom - margin.top)
                .style('fill', Lsv.get_color(ii, BREWER_PALETTE, 1));
            svgContainer.append("rect")
                .attr("x", last_incl_pos)
                .attr("y", margin.top)
                .attr("width", Math.round((width / 2 - margin.right) * lsv.excl_incl[ii][1]))
                .attr("height", height - margin.bottom - margin.top)
                .style('fill', Lsv.get_color(ii, BREWER_PALETTE, 1));

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
                    .attr("fill", Lsv.get_color(ii, BREWER_PALETTE, 1))
                    .text(Math.round((width / 2 - margin.left) * lsv.excl_incl[ii][0]));

            }

            if (Math.round((width / 2 - margin.right) * lsv.excl_incl[ii][1]) >= 1) {
                last_incl_pos += Math.round((width / 2 - margin.right) * lsv.excl_incl[ii][1]);
                svgContainer.append("text")
                    .attr("x", last_incl_pos)
                    .attr("y", height)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "9px")
                    .attr("fill", Lsv.get_color(ii, BREWER_PALETTE, 1))
                    .text(Math.round((width / 2 - margin.right) * lsv.excl_incl[ii][1]));
            }
        }

        // Draw separator
        svgContainer.append("line")
            .attr("x1", width / 2)
            .attr("y1", 0)
            .attr("x2", width / 2)
            .attr("y2", height - margin.bottom + 2)
            .attr("stroke-group_width", 2)
            .attr("stroke-opacity", .8)
            .attr("stroke", "black");


        return svgContainer;

    }

    renderFloatingLegend(canvas) {
        var ctx = canvas.getContext("2d");

        // Clear previous draw
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        var MARGINS = [10, 2, 2, 2];
        var SEP_FIG_TEXT = canvas.height * .05;
        var SEP_FIG = canvas.width * .02;
        var num_fig = 10;
        var area_figures = [
            canvas.width - MARGINS[0] - MARGINS[1] - (num_fig - 1) * SEP_FIG,
            canvas.height * .7 - MARGINS[2] - SEP_FIG_TEXT
        ];
        var area_texts = [canvas.width - MARGINS[0] - MARGINS[1], canvas.height * .3 - MARGINS[2] - SEP_FIG_TEXT];
        var legend_line_length = 20;
        var x = MARGINS[0];
        var y = MARGINS[2];
        ctx.font = "7pt Arial";
        ctx.textAlign = "center";

        /**
         * Legend exons_obj
         * */
        // DB & RNASeq
        ctx.strokeStyle = "rgba(0, 0, 0, 0.8)";
        ctx.fillStyle = "rgba(0, 0, 0, 0.2)";
        this.draw_rectangle(ctx, x, y, Math.round(area_figures[0] / num_fig - SEP_FIG), Math.round(area_figures[1]), true);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("DB & RNASeq", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // RNASeq Only
        ctx.strokeStyle = Lsv.get_color(2, BREWER_PALETTE, .8);
        ctx.fillStyle = Lsv.get_color(2, BREWER_PALETTE, .2);
        this.draw_rectangle(ctx, x, y, Math.round(area_figures[0] / num_fig - SEP_FIG), Math.round(area_figures[1]), true);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("RNASeq Only", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // DB Only
        ctx.strokeStyle = "rgba(0, 0, 0, 0.8)";
        ctx.lineWidth = 2;
        ctx.fillStyle = "rgba(255, 255, 255, .5)";
        ctx.setLineDash([5, 5]);
        this.draw_rectangle(ctx, Math.round(x), y, Math.round(area_figures[0] / num_fig - SEP_FIG), Math.round(area_figures[1]), true);
        ctx.setLineDash([]);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("DB Only", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        /**
         * Legend junctions
         * */
        // DB & RNASeq
        ctx.strokeStyle = combined_colors['sa'];
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.arc(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round(y + area_figures[1]), (area_figures[0] / num_fig - SEP_FIG) / 2, -Math.PI, 0);
        ctx.stroke();
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("DB & RNASeq", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // RNASeq Only
        ctx.strokeStyle = combined_colors['s'];
        ctx.beginPath();
        ctx.arc(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round(y + area_figures[1]), (area_figures[0] / num_fig - SEP_FIG) / 2, -Math.PI, 0);
        ctx.stroke();
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("RNASeq Only", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // DB Only
        ctx.strokeStyle = combined_colors['ao'];
        ctx.setLineDash([5, 5]);
        ctx.beginPath();
        ctx.arc(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round(y + area_figures[1]), (area_figures[0] / num_fig - SEP_FIG) / 2, -Math.PI, 0);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("DB Only", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        /**
         * Legend number of reads
         * */
        // DB & RNASeq example chosen
        ctx.strokeStyle = combined_colors['sa'];
        ctx.lineWidth = 1.2;
        ctx.font = "8pt Arial";
        var font_height = 9;
        ctx.beginPath();
        ctx.arc(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round(y + area_figures[1]), (area_figures[0] / num_fig - 2 * SEP_FIG) / 2, -Math.PI, 0);
        ctx.stroke();
        Lsv.renderNumReads(ctx, Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), MARGINS[2] + font_height, 32);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.font = "7pt Arial";
        ctx.fillText("RNASeq reads", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        /**
         * Legend Intron Retention
         * */

        ctx.lineWidth = 1.2;
        ctx.strokeStyle = combined_colors['sa'];
        ctx.fillStyle = combined_colors['sa'] + "66";
        this.draw_rectangle(ctx, Math.round(x + (area_figures[0] / num_fig) / 3 - SEP_FIG), y + area_figures[1] / 4, Math.round((area_figures[0] / num_fig - SEP_FIG) * 2 / 3) + 4, Math.round(area_figures[1] / 2), true);
        ctx.strokeStyle = "rgba(0, 0, 0, 0.8)";
        ctx.fillStyle = "rgba(0, 0, 0, 0.2)";
        this.draw_rectangle(ctx, x, y, Math.round((area_figures[0] / num_fig) / 3 - SEP_FIG), Math.round(area_figures[1]), true);
        this.draw_rectangle(ctx, Math.round(x + (area_figures[0] / num_fig) * 2 / 3), y, Math.round((area_figures[0] / num_fig) / 3 - SEP_FIG), Math.round(area_figures[1]), true);

        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("Intron Ret.", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        ctx.font = "22pt Arial";
        ctx.fillText("↳", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), y + area_figures[1]);
        ctx.font = "7pt Arial";
        ctx.fillText("DB TSS", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        ctx.font = "22pt Arial";
        ctx.fillText("^", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), y + area_figures[1] + 4);
        ctx.font = "7pt Arial";
        ctx.fillText("DB TES", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);



        ctx.lineWidth = 1;
    }

    renderFloatingLegendLr(canvas) {
        var ctx = canvas.getContext("2d");

        // Clear previous draw
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        var MARGINS = [10, 2, 2, 2];
        var SEP_FIG_TEXT = canvas.height * .05;
        var SEP_FIG = canvas.width * .01;
        console.log(SEP_FIG)
        var num_fig = 12;
        var area_figures = [
            canvas.width - MARGINS[0] - MARGINS[1] - (num_fig - 1) * SEP_FIG,
            canvas.height * .7 - MARGINS[2] - SEP_FIG_TEXT
        ];
        console.log(area_figures)
        //var area_texts = [canvas.width - MARGINS[0] - MARGINS[1], canvas.height * .3 - MARGINS[2] - SEP_FIG_TEXT];
        //var legend_line_length = 20;
        var x = MARGINS[0];
        var y = MARGINS[2];
        ctx.font = "7pt Arial";
        ctx.textAlign = "center";

        /**
         * Legend exons_obj
         * */

        // DB & RNASeq
        ctx.strokeStyle = "rgba(0, 0, 0, 0.8)";
        ctx.fillStyle = "rgba(0, 0, 0, 0.2)";
        this.draw_rectangle(ctx, x, y, Math.round(area_figures[0] / num_fig - SEP_FIG), Math.round(area_figures[1]), true);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("DB & SR/LR", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // sla
        ctx.strokeStyle = combined_colors['sla'];
        ctx.fillStyle = combined_colors['sla']+'66';
        //this.draw_rectangle(ctx, x, y + Math.round(area_figures[1]/2), Math.round(area_figures[0] / num_fig - SEP_FIG + 1), Math.round(area_figures[1]/2), true);
        this.draw_rectangle(ctx, x, y + Math.round(area_figures[1]/2), Math.round((area_figures[0] / num_fig - SEP_FIG + 1)/4), Math.round(area_figures[1]/2), true);
        this.draw_rectangle(ctx, x + Math.round((area_figures[0] / num_fig - SEP_FIG + 1)*(3/4)), y + Math.round(area_figures[1]/2), Math.round((area_figures[0] / num_fig - SEP_FIG + 1)/4), Math.round(area_figures[1]/2), true);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("DB+SR+LR", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        ctx.strokeStyle = combined_colors['sla'];
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        //ctx.ellipse(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round((y + area_figures[1])/2), (area_figures[0] / num_fig - SEP_FIG) / 2, (area_figures[0] / num_fig - SEP_FIG) / 4, 0, -Math.PI, 0);
        ctx.ellipse(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round((y + area_figures[1])/2), (area_figures[0] / num_fig - SEP_FIG) / 4, (area_figures[0] / num_fig - SEP_FIG) / 4, 0, -Math.PI, 0);
        ctx.stroke();
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // sl
        ctx.strokeStyle = combined_colors['sl'];
        ctx.fillStyle = combined_colors['sl']+'66';
        this.draw_rectangle(ctx, x, y + Math.round(area_figures[1]/2), Math.round((area_figures[0] / num_fig - SEP_FIG + 1)/4), Math.round(area_figures[1]/2), true);
        this.draw_rectangle(ctx, x + Math.round((area_figures[0] / num_fig - SEP_FIG + 1)*(3/4)), y + Math.round(area_figures[1]/2), Math.round((area_figures[0] / num_fig - SEP_FIG + 1)/4), Math.round(area_figures[1]/2), true);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("SR+LR", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        ctx.strokeStyle = combined_colors['sl'];
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.ellipse(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round((y + area_figures[1])/2), (area_figures[0] / num_fig - SEP_FIG) / 4, (area_figures[0] / num_fig - SEP_FIG) / 4, 0, -Math.PI, 0);
        ctx.stroke();
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // l
        ctx.strokeStyle = combined_colors['l'];
        ctx.fillStyle = combined_colors['l']+'66';
        this.draw_rectangle(ctx, x, y + Math.round(area_figures[1]/2), Math.round((area_figures[0] / num_fig - SEP_FIG + 1)/4), Math.round(area_figures[1]/2), true);
        this.draw_rectangle(ctx, x + Math.round((area_figures[0] / num_fig - SEP_FIG + 1)*(3/4)), y + Math.round(area_figures[1]/2), Math.round((area_figures[0] / num_fig - SEP_FIG + 1)/4), Math.round(area_figures[1]/2), true);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("LR only", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        ctx.strokeStyle = combined_colors['l'];
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.ellipse(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round((y + area_figures[1])/2), (area_figures[0] / num_fig - SEP_FIG) / 4, (area_figures[0] / num_fig - SEP_FIG) / 4, 0, -Math.PI, 0);
        ctx.stroke();
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // la
        ctx.strokeStyle = combined_colors['la'];
        ctx.fillStyle = combined_colors['la']+'66';
        this.draw_rectangle(ctx, x, y + Math.round(area_figures[1]/2), Math.round((area_figures[0] / num_fig - SEP_FIG + 1)/4), Math.round(area_figures[1]/2), true);
        this.draw_rectangle(ctx, x + Math.round((area_figures[0] / num_fig - SEP_FIG + 1)*(3/4)), y + Math.round(area_figures[1]/2), Math.round((area_figures[0] / num_fig - SEP_FIG + 1)/4), Math.round(area_figures[1]/2), true);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("DB+LR", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        ctx.strokeStyle = combined_colors['la'];
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.ellipse(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round((y + area_figures[1])/2), (area_figures[0] / num_fig - SEP_FIG) / 4, (area_figures[0] / num_fig - SEP_FIG) / 4, 0, -Math.PI, 0);
        ctx.stroke();
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // s
        ctx.strokeStyle = combined_colors['s'];
        ctx.fillStyle = combined_colors['s']+'66';
        this.draw_rectangle(ctx, x, y + Math.round(area_figures[1]/2), Math.round((area_figures[0] / num_fig - SEP_FIG + 1)/4), Math.round(area_figures[1]/2), true);
        this.draw_rectangle(ctx, x + Math.round((area_figures[0] / num_fig - SEP_FIG + 1)*(3/4)), y + Math.round(area_figures[1]/2), Math.round((area_figures[0] / num_fig - SEP_FIG + 1)/4), Math.round(area_figures[1]/2), true);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("SR only", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        ctx.strokeStyle = combined_colors['s'];
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.ellipse(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round((y + area_figures[1])/2), (area_figures[0] / num_fig - SEP_FIG) / 4, (area_figures[0] / num_fig - SEP_FIG) / 4, 0, -Math.PI, 0);
        ctx.stroke();
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // sa
        ctx.strokeStyle = combined_colors['sa'];
        ctx.fillStyle = combined_colors['sa']+'66';
        this.draw_rectangle(ctx, x, y + Math.round(area_figures[1]/2), Math.round((area_figures[0] / num_fig - SEP_FIG + 1)/4), Math.round(area_figures[1]/2), true);
        this.draw_rectangle(ctx, x + Math.round((area_figures[0] / num_fig - SEP_FIG + 1)*(3/4)), y + Math.round(area_figures[1]/2), Math.round((area_figures[0] / num_fig - SEP_FIG + 1)/4), Math.round(area_figures[1]/2), true);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("DB+SR", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        ctx.strokeStyle = combined_colors['sa'];
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.ellipse(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round((y + area_figures[1])/2), (area_figures[0] / num_fig - SEP_FIG) / 4, (area_figures[0] / num_fig - SEP_FIG) / 4, 0, -Math.PI, 0);
        ctx.stroke();
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // ao
        ctx.strokeStyle = "rgba(0, 0, 0, 0.8)";
        ctx.fillStyle = "rgba(255, 255, 255, .5)";
        ctx.lineWidth = 2;
        ctx.setLineDash([5, 5]);
        this.draw_rectangle(ctx, x, y + Math.round(area_figures[1]/2), Math.round((area_figures[0] / num_fig - SEP_FIG + 1)/4), Math.round(area_figures[1]/2), true);
        this.draw_rectangle(ctx, x + Math.round((area_figures[0] / num_fig - SEP_FIG + 1)*(3/4)), y + Math.round(area_figures[1]/2), Math.round((area_figures[0] / num_fig - SEP_FIG + 1)/4), Math.round(area_figures[1]/2), true);
        ctx.strokeStyle = "rgba(0, 0, 0, 0.8)";
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.ellipse(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round((y + area_figures[1])/2), (area_figures[0] / num_fig - SEP_FIG) / 4, (area_figures[0] / num_fig - SEP_FIG) / 4, 0, -Math.PI, 0);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("DB only", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // DB & RNASeq example chosen
        ctx.strokeStyle = combined_colors['s'];
        ctx.lineWidth = 1.2;
        ctx.font = "8pt Arial";
        var font_height = 9;
        ctx.beginPath();
        ctx.arc(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round(y + area_figures[1]), (area_figures[0] / num_fig - 2 * SEP_FIG) * (1/3), -Math.PI, 0);
        ctx.stroke();
        Lsv.renderNumReads(ctx, Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), MARGINS[2] + font_height, "69|42");
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.font = "7pt Arial";
        ctx.fillText("SR|LR reads", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        ctx.strokeStyle = "rgba(0, 0, 0, 0.8)";
        ctx.fillStyle = "rgba(0, 0, 0, 0.2)";
        ctx.lineWidth = 1.2;
        this.draw_rectangle(ctx, x, y, Math.round((area_figures[0] / num_fig) / 3 - SEP_FIG), Math.round(area_figures[1]), true);
        this.draw_rectangle(ctx, Math.round(x + (area_figures[0] / num_fig) * 2 / 3), y, Math.round((area_figures[0] / num_fig) / 3 - SEP_FIG), Math.round(area_figures[1]), true);
        this.draw_rectangle(ctx, Math.round(x + (area_figures[0] / num_fig) / 3 - SEP_FIG), y + area_figures[1] / 4, Math.round(((area_figures[0] / num_fig - SEP_FIG) * 2 / 3) - 6), Math.round(area_figures[1] / 2), true);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("Intron Ret.", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        ctx.font = "22pt Arial";
        ctx.fillText("↳", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), y + area_figures[1]);
        ctx.font = "7pt Arial";
        ctx.fillText("DB TSS", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        ctx.font = "22pt Arial";
        ctx.fillText("^", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), y + area_figures[1] + 4);
        ctx.font = "7pt Arial";
        ctx.fillText("DB TES", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);



        ctx.lineWidth = 1;
    }

    static renderNumReads(ctx, x, y, num_reads) {
        if (parseInt(num_reads) === 0) return;
        ctx.fillStyle = "rgba(0, 0, 0, .8)";
        ctx.font = "9pt Arial";
        ctx.textAlign = "center";
        ctx.fillText(num_reads, x, y - 2);
    }
}


