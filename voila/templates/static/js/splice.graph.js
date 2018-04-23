/**
 * Created by abarrera on 2/18/14.
 */

//$( document ).ready(function(){

window.splicegraph = function () {
    // Adding dash lines to the Canvas Rendering
    CanvasRenderingContext2D.prototype.dashedLine = function (x1, y1, x2, y2, dashLen) {

        if (dashLen == undefined) dashLen = 2;
        this.moveTo(x1, y1);

        var dX = x2 - x1;
        var dY = y2 - y1;
        var dashes = Math.floor(Math.sqrt(dX * dX + dY * dY) / dashLen);
        var dashX = dX / dashes;
        var dashY = dY / dashes;

        var q = 0;
        while (q++ < dashes) {
            x1 += dashX;
            y1 += dashY;
            this[q % 2 == 0 ? 'moveTo' : 'lineTo'](x1, y1);
        }
        this[q % 2 == 0 ? 'moveTo' : 'lineTo'](x2, y2);
    };


    // Array Remove - By John Resig (MIT Licensed)
    Array.prototype.remove = function (from, to) {
        var rest = this.slice((to || from) + 1 || this.length);
        this.length = from < 0 ? this.length + from : from;
        return this.push.apply(this, rest);
    };

    function drawLine(contextO, startx, starty, endx, endy) {
        contextO.beginPath();
        contextO.moveTo(startx, starty);
        contextO.lineTo(endx, endy);
        contextO.closePath();
        contextO.stroke();
    }

    // drawRectangle - draws a rectangle on a canvas context using the dimensions specified
    function drawRectangle(contextO, x, y, w, h, fill) {
        contextO.beginPath();
        contextO.rect(x, y, w, h);
        contextO.closePath();
        contextO.stroke();
        if (fill) contextO.fill();
    }

    function drawDashedLine(contextO, startx, starty, endx, endy, dashLen) {
        contextO.beginPath();
        contextO.dashedLine(startx, starty, endx, endy, dashLen);
        contextO.stroke();
    }

    function drawDashedRectangle(contextO, x, y, w, h, dashLen, fill) {
        contextO.beginPath();
        drawDashedLine(contextO, x, y, x + w, y, dashLen);
        drawDashedLine(contextO, x + w, y, x + w, y - h, dashLen);
        drawDashedLine(contextO, x + w, y - h, x, y - h, dashLen);
        drawDashedLine(contextO, x, y - h, x, y, dashLen);
        contextO.stroke();
        if (fill) contextO.fill();
    }


    function drawArrow(context, fromx, fromy, tox, toy, headlen) {
        var angle = Math.atan2(toy - fromy, tox - fromx);
        context.beginPath();
        context.moveTo(fromx, fromy);
        context.lineTo(tox, toy);
        context.lineTo(tox - headlen * Math.cos(angle - Math.PI / 6), toy - headlen * Math.sin(angle - Math.PI / 6));
        context.moveTo(tox, toy);
        context.lineTo(tox - headlen * Math.cos(angle + Math.PI / 6), toy - headlen * Math.sin(angle + Math.PI / 6));
        context.stroke();
    }

    // Utils functions
    function forEach(array, action) {
        for (var i = 0; i < array.length; i++)
            action(array[i]);
    }

    function reduce(combine, base, array) {
        forEach(array, function (element) {
            base = combine(base, element);
        });
        return base;
    }

    function add(a, b) {
        return a + b;
    }

    function sum(numbers) {
        return reduce(add, 0, numbers);
    }

    function length_junction(junction) {
        return junction.coords[1] - junction.coords[0];
    }

    function map(func, array) {
        var result = [];
        forEach(array, function (element) {
            result.push(func(element));
        });
        return result;
    }

    // Get Color from Brewer Palette
    BREWER_PALETTE = [
        [228, 26, 28],
        [55, 126, 184],
        [77, 175, 74],
        [152, 78, 163],
        [255, 127, 0],
//            [255, 255, 51],
        [166, 86, 40],
        [247, 129, 191],
        [153, 153, 153],

        [28, 126, 128],
        [155, 226, 29],
        [177, 275, 19],
        [252, 178, 8],
        [55, 227, 100],
//            [55,55,151],
        [266, 186, 140],
        [47, 229, 36],
        [253, 253, 253]
    ];

    function getColor(colorNumber, palette, hue) {
        colorNumber = colorNumber % 16;
        return "rgba(" + palette[colorNumber].toString() + ", " + hue + ")";
    }

    function render_intron(canvas, pixel_factor, margin, percen_intron, coords) {

        var intron_height = (percen_intron * canvas.height) / 2; // canvas.group_height-2*margin[2];

        var ctx = canvas.getContext("2d");
        ctx.strokeStyle = "rgba(0, 0, 0, 0.6)";
        ctx.fillStyle = "rgba(0, 0, 0, 0.1)";

        ctx.strokeRect(margin[0] + Math.round((coords[0]) * pixel_factor),
            canvas.height - margin[3] - intron_height / 2,
            Math.round((coords[1] - coords[0]) * pixel_factor),
            -intron_height
        );

        ctx.fillRect(
            margin[0] + Math.round((coords[0]) * pixel_factor),
            canvas.height - margin[3] - intron_height / 2,
            Math.round((coords[1] - coords[0]) * pixel_factor),
            -intron_height
        );

    }

    function render_intron_retained(canvas, intron, pixel_factor, margin, percen_exon, junc_count) {
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

    }

    function render_exon(canvas, exon, pixel_factor, margin, percen_exon, counter_exon) {

        var exon_height = percen_exon * canvas.height; // canvas.group_height-2*margin[2];

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


    }

    function renderFloatingLegend(canvas) {
        var ctx = canvas.getContext("2d");

        // Clear previous draw
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        var MARGINS = [10, 2, 2, 2];
        var SEP_FIG_TEXT = canvas.height * .05;
        var SEP_FIG = canvas.width * .02;
        var num_fig = 8;
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
        drawRectangle(ctx, x, y, Math.round(area_figures[0] / num_fig - SEP_FIG), Math.round(area_figures[1]), true);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("DB & RNASeq", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // RNASeq Only
        ctx.strokeStyle = getColor(2, BREWER_PALETTE, .8);
        ctx.fillStyle = getColor(2, BREWER_PALETTE, .2);
        drawRectangle(ctx, x, y, Math.round(area_figures[0] / num_fig - SEP_FIG), Math.round(area_figures[1]), true);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("RNASeq Only", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // DB Only
        ctx.strokeStyle = "rgba(0, 0, 0, 0.8)";
        ctx.lineWidth = 2;
        ctx.fillStyle = "rgba(255, 255, 255, .5)";
        ctx.setLineDash([5, 5]);
        drawRectangle(ctx, Math.round(x), y, Math.round(area_figures[0] / num_fig - SEP_FIG), Math.round(area_figures[1]), true);
        ctx.setLineDash([]);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("DB Only", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        /**
         * Legend junctions
         * */
        // DB & RNASeq
        ctx.strokeStyle = 'red';
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.arc(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round(y + area_figures[1]), (area_figures[0] / num_fig - SEP_FIG) / 2, -Math.PI, 0);
        ctx.stroke();
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("DB & RNASeq", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // RNASeq Only
        ctx.strokeStyle = getColor(2, BREWER_PALETTE, .8);
        ctx.beginPath();
        ctx.arc(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round(y + area_figures[1]), (area_figures[0] / num_fig - SEP_FIG) / 2, -Math.PI, 0);
        ctx.stroke();
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("RNASeq Only", Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        // DB Only
        ctx.strokeStyle = "rgba(0, 0, 0, 0.8)";
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
        ctx.strokeStyle = 'red';
        ctx.lineWidth = 1.2;
        ctx.font = "8pt Arial";
        var font_height = 9;
        ctx.beginPath();
        ctx.arc(Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), Math.round(y + area_figures[1]), (area_figures[0] / num_fig - 2 * SEP_FIG) / 2, -Math.PI, 0);
        ctx.stroke();
        renderNumReads(ctx, Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), MARGINS[2] + font_height, 32);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.font = "7pt Arial";
        ctx.fillText("RNASeq reads", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);
        x = x + area_figures[0] / num_fig + SEP_FIG;

        /**
         * Legend Intron Retention
         * */
        ctx.strokeStyle = "rgba(0, 0, 0, 0.8)";
        ctx.fillStyle = "rgba(0, 0, 0, 0.2)";
        ctx.lineWidth = 1.2;
        drawRectangle(ctx, x, y, Math.round((area_figures[0] / num_fig) / 3 - SEP_FIG), Math.round(area_figures[1]), true);
        drawRectangle(ctx, Math.round(x + (area_figures[0] / num_fig) * 2 / 3), y, Math.round((area_figures[0] / num_fig) / 3 - SEP_FIG), Math.round(area_figures[1]), true);
        drawRectangle(ctx, Math.round(x + (area_figures[0] / num_fig) / 3 - SEP_FIG), y + area_figures[1] / 4, Math.round((area_figures[0] / num_fig - SEP_FIG) * 2 / 3), Math.round(area_figures[1] / 2), true);
        ctx.fillStyle = "rgba(0, 0, 0, 1)";
        ctx.fillText("Intron Ret.", x + Math.round((area_figures[0] / num_fig - SEP_FIG) / 2), canvas.height - MARGINS[3]);


        ctx.lineWidth = 1;

    }

    // render number of reads
    function renderNumReads(ctx, x, y, num_reads) {
        // render number as text
        if (parseInt(num_reads) === 0) return;
        ctx.fillStyle = "rgba(0, 0, 0, .8)";
        ctx.font = "9pt Arial";
        ctx.textAlign = "center";
        ctx.fillText(num_reads, x, y - 2);
    }

    function render_junction(canvas, junction, largest_junction, pixel_factor, margin, percen_exon, length_region) {

        var min_height = percen_exon * canvas.height + 10;
        var ctx = canvas.getContext("2d");

        var height_junction = canvas.height - margin[2] - margin[3];
        height_junction -=
            Math.min(
                (Math.max(min_height, (height_junction * .9) * ((junction.coords[1] - junction.coords[0]) / largest_junction))),
                height_junction * .9
            );

        var middle_point = Math.round((junction.coords[0] + junction.coords[1]) / 2);

        // line color
        ctx.strokeStyle = 'red';
        ctx.lineWidth = 1.2;

        // render Junction
        if (junction.type_junction === 1) {
            ctx.strokeStyle = getColor(2, BREWER_PALETTE, .8);
            ctx.fillStyle = getColor(2, BREWER_PALETTE, .2);
        }

        if (junction.type_junction === 2) {
            ctx.strokeStyle = "rgba(0, 0, 0, .8)";
            drawDashedLine(ctx,
                Math.round(margin[0] + junction.coords[0] * pixel_factor),
                canvas.height - margin[3] - percen_exon * canvas.height,
                Math.round(margin[0] + middle_point * pixel_factor),
                Math.round(height_junction), 2
            );
            drawDashedLine(ctx,
                Math.round(margin[0] + middle_point * pixel_factor),
                Math.round(height_junction),
                Math.round(margin[0] + junction.coords[1] * pixel_factor),
                canvas.height - margin[3] - percen_exon * canvas.height, 2
            );
        } else {

            drawLine(ctx,
                Math.round(margin[0] + junction.coords[0] * pixel_factor),
                canvas.height - margin[3] - percen_exon * canvas.height,
                Math.round(margin[0] + middle_point * pixel_factor),
                Math.round(height_junction)
            );
            drawLine(ctx,
                Math.round(margin[0] + middle_point * pixel_factor),
                Math.round(height_junction),
                Math.round(margin[0] + junction.coords[1] * pixel_factor),
                canvas.height - margin[3] - percen_exon * canvas.height
            );
        }

        // mark 3-primes and 5-primes
        drawDashedLine(ctx,
            Math.round(margin[0] + junction.coords[1] * pixel_factor), canvas.height - margin[3] - percen_exon * canvas.height,
            Math.round(margin[0] + junction.coords[1] * pixel_factor), canvas.height - margin[3], 2);

        drawDashedLine(ctx,
            Math.round(margin[0] + junction.coords[0] * pixel_factor), canvas.height - margin[3] - percen_exon * canvas.height,
            Math.round(margin[0] + junction.coords[0] * pixel_factor), canvas.height - margin[3], 2);

        // Num Reads
        ctx.font = "8pt Arial";
        renderNumReads(ctx, Math.round(margin[0] + middle_point * pixel_factor), Math.round(height_junction), junction.num_reads);

        ctx.lineWidth = 1;
    }


    function reshape_intron(exon1, exon2, reduce_exon) {

        function constant_size(intron_size) {
            var MAX_INTRON = 300;
            if (intron_size > MAX_INTRON) {
                return intron_size - MAX_INTRON;
            } else {
                return 0;
            }
        }

        var reshaped_intron = constant_size(exon2.coords[0] - exon1.coords[1]);
        if (reduce_exon) {
            if ((exon1.coords[1] - exon1.coords[0]) > MAX_INTRON) {
                reshaped_intron += MAX_INTRON;
            }
        }
        return reshaped_intron;
    }

    function resize_exon(exon, reduce_exon) {
        var MAX_EXON = 300;
        if (reduce_exon) {
            var exon_length = exon[1] - exon[0];
            if (exon_length > MAX_EXON) {
                exon[1] -= exon_length - MAX_EXON;
                return exon
            }
        }

        return exon;
    }

    function map_exon_list(exons, junctions) {
        var i;
        var k;
        var j;
        var reshape_exons = false,
            exons_mapped_tmp = [],
            exon_tmp,
            index_exons_to_update = [];  // For exons_obj that need to be updated (shorted)

        // First offset (first exon can have  long intron
        var offset = reshape_intron({'coords': [0, 0]}, exons[0], reshape_exons);
        var acc_offset = offset;

        // Note: to account for overlapping exons_obj where exons_obj within a very large exon have long introns, we should
        // ^^^^^ store the last
        var coords_extra = [];
        for (k = 0; k < exons[0].coords_extra.length; k++) {
            coords_extra.push(map(function (x) {
                return add(x, -acc_offset);
            }, exons[0].coords_extra[k]));
        }
        exon_tmp = {
            'coords': map(function (x) {
                return add(x, -acc_offset);
            }, resize_exon(exons[0].coords, reshape_exons)),
            'type': exons[0].type_exon,
            'intron_retention': exons[0].intron_retention,
            'lsv_type': exons[0].lsv_type
        };
        exon_tmp.coords_extra = coords_extra;

        exons_mapped_tmp[0] = exon_tmp;
        last_end = exon_tmp.coords[1];


        for (i = 0; i < exons[0]['a3'].length; i++) {
            junctions[exons[0]['a3'][i]].coords[1] -= offset;
        }

        for (i = 0; i < exons[0]['a5'].length; i++) {
            junctions[exons[0]['a5'][i]].coords[0] -= offset;
        }

        for (i = 1; i < exons.length; i++) {
            offset = reshape_intron(exons[i - 1], exons[i], reshape_exons);
            acc_offset += offset;

            // Check if there are exons_obj to make shorter (intron retention)
            coords_extra = [];
            for (k = 0; k < exons[i].coords_extra.length; k++) {
                coords_extra.push(map(function (x) {
                    return add(x, -acc_offset);
                }, exons[i].coords_extra[k]));
            }
            exon_tmp = {
                'coords': map(function (x) {
                    return add(x, -acc_offset);
                }, resize_exon(exons[i].coords, reshape_exons)),
                'type': exons[i].type_exon,
                'intron_retention': exons[i].intron_retention,
                'lsv_type': exons[i].lsv_type
            };
            exon_tmp.coords_extra = coords_extra;

            exons_mapped_tmp[i] = exon_tmp;

            // Check if any previous exon needs to be updated
            if (exons[i].coords[0] - exons[i - 1].coords[1] < 0) {
                index_exons_to_update.push(i - 1);
            }

            if (offset) {
                for (k = 0; k < index_exons_to_update.length; k++) {
                    if (exons[index_exons_to_update[k]].coords[1] > exons[i].coords[0]) {
                        exons_mapped_tmp[index_exons_to_update[k]].coords[1] -= offset;
                    }
                }
            }

            for (j = 0; j < exons[i]['a3'].length; j++) {
                junctions[exons[i]['a3'][j]].coords[1] -= acc_offset;
            }

            for (j = 0; j < exons[i]['a5'].length; j++) {
                junctions[exons[i]['a5'][j]].coords[0] -= acc_offset;
            }
        }

        return exons_mapped_tmp;
    }


    var render_lsv_marker = function (canvas, exon, exon2, pixel_factor, margin) {
        var ctx = canvas.getContext("2d");
        var x1, x2,
            exon_length = exon.coords[1] - exon.coords[0];

        if (exon.lsv_type === 1) {
            ctx.strokeStyle = "rgba(255, 140, 50, 0.3)";
            ctx.fillStyle = "rgba(255, 140, 50, 0.2)";
            x1 = exon.coords[1] - exon_length * .4;
            x2 = exon.coords[1] + Math.min(exon_length * .4, (exon2.coords[0] - exon.coords[1]) / 2);
        } else {
            ctx.strokeStyle = "rgba(50, 140, 255, 0.3)";
            ctx.fillStyle = "rgba(50, 140, 255, 0.2)";
            x1 = exon.coords[0] - Math.min(exon_length * .4, (exon.coords[0] - exon2.coords[1]) / 2);
            x2 = exon.coords[0] + (exon.coords[1] - exon.coords[0]) * .4;
        }


        ctx.strokeRect(margin[0] + Math.round(x1 * pixel_factor), canvas.height - margin[3] / 2, Math.round((x2 - x1) * pixel_factor), -(canvas.height * .1 + margin[3]));
        ctx.fillRect(margin[0] + Math.round(x1 * pixel_factor), canvas.height - margin[3] / 2, Math.round((x2 - x1) * pixel_factor), -(canvas.height * .1 + margin[3]));

    };

    function renderSpliceGraph(canvas) {
        if (canvas.getContext) {

            // Clear previous draw
            var ctx = canvas.getContext("2d");
            ctx.clearRect(0, 0, canvas.width, canvas.height);

            // NOTE: The list of exons_obj is assumed to be sorted.
            // NEW 20140318: Exons have type: 0=Both reads & annotations; 1=New (no annotations); 2=Missing (annotated, no reads)
            var genomic_data = JSON.parse(canvas.getAttribute('data-exon-list').replace(/'/g, "\""));   // RESTORE ME WHEN DEBUGGING FINISH!

            var MARGIN = [10, 10, 0, 10];    // Global E, W, N, S

            var junctions = genomic_data.junctions;
            var exons = map_exon_list(genomic_data.exons, junctions, canvas, MARGIN);

            var length_region = exons[exons.length - 1].coords[1];

            // Until Jordi finish debugging his code, here is a tweak to fix orphan junctions:
            for (var junc_index = 0; junc_index < junctions.length; junc_index++) {
                junctions[junc_index].coords[0] = Math.min(length_region, Math.max(junctions[junc_index].coords[0], 0));
                junctions[junc_index].coords[1] = Math.min(length_region, Math.max(junctions[junc_index].coords[1], 0));
            }

            // Compute scaled pixel factor
            var pixel_factor = Math.min(canvas.width - (MARGIN[0] + MARGIN[1]), length_region) / length_region;
            var largest_junct = reduce(Math.max, 0, map(length_junction, junctions));

            var num_exons = exons.length;
            var percentage_exon = 0.1;

            for (var i = 0; i < junctions.length; i++) {
                render_junction(canvas, junctions[i], largest_junct, pixel_factor, MARGIN, percentage_exon, length_region);
            }

            for (var i = 0; i < num_exons; i++) {
                var exon = exons[i];
                render_exon(canvas, exon, pixel_factor, MARGIN, percentage_exon, i + 1);

                // Render intron retention
                if (i < num_exons - 1 && exon.intron_retention && exons[i + 1].intron_retention) {
                    render_intron(canvas, pixel_factor, MARGIN, percentage_exon, [exon.coords[1] + 1, exons[i + 1].coords[0] - 1]);
                }

                // Render lsv shadow rectangle
                if (exon.lsv_type) {
                    exon.lsv_type = 3;
                    if (exon.lsv_type === 1) {  // Single source
                        render_lsv_marker(canvas, exon, exons[i + 1], pixel_factor, MARGIN);
                    }
                    if (exon.lsv_type === 2) {  // Single target
                        render_lsv_marker(canvas, exon, exons[i - 1], pixel_factor, MARGIN);
                    }
                    if (exon.lsv_type === 3) {  // Both single source and single target
                        exon.lsv_type = 1;
                        render_lsv_marker(canvas, exon, exons[i + 1], pixel_factor, MARGIN);
                        exon.lsv_type = 2;
                        render_lsv_marker(canvas, exon, exons[i - 1], pixel_factor, MARGIN);
                        exon.lsv_type = 3;
                    }
                }
            }

            return genomic_data;
        }
    }

    function renderLsvSpliceGraph(canvas, gene) {
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


            if (gene) {
                var exons_gene = gene.exons.filter(function (v) {
                    return v.exon_type < 3 && !v.intron_retention;
                });
            }

            if (exon_lsv_coords) {
                var regex_coords = /[\[|(](\d+?),\s*(\d+?)[\]|)]/.exec(exon_lsv_coords);
                exon_lsv_coords = [regex_coords[1], regex_coords[2]];

                // Find LSV exon number in the splice graph
                var lsv_exon = exons_gene.find(function (element) {
                    return element.start == exon_lsv_coords[0] && element.end == exon_lsv_coords[1]
                });
                exon_lsv_number = exons_gene.indexOf(lsv_exon) + 1;

                // If negative strand, complement the exon ordinal
                if (gene.strand === '-') {
                    exon_lsv_number = exons_gene.length - exon_lsv_number + 1;
                }
            }

            var lsv_data = canvas.getAttribute('data-lsv-string');
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

                if (lsvs_fields[1] == 0) continue;
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
                var number_exon = (direction > 0 ? i : i + 1 );
                number_exon = (number_exon == 0 || number_exon == num_exons ? exon_lsv_number : '');
                render_exon(canvas, exon, pixel_factor, margins, percentage_exon, number_exon);
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
                render_intron_retained(canvas, intron, pixel_factor, margins, percentage_exon, intron_ret_i - 1);
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

                var coords_x_start_e = coords[0] + ((exon_width / 2) / num_ss ) * (ss - 1);
                var coords_x_target_e = null;

                if (lsvs_fields[1] == 0) {
                    coords_x_target_e = coords_x_start_e;
                    count_starts_ends++;
                } else {
                    var target_exon = (direction > 0 ? exons[lsvs_fields[1][0]] : exons[lsvs_fields[1][0] - 1] );
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
                if (direction > 0 && ss != num_ss || direction < 0 && ss != 1) {

                    if (!rendered_exons[lsvs_fields[1][0]]) {  // Render ssites only if they haven't been rendered
                        ctx.strokeStyle = "rgba(0, 0, 0, 0.6)";
                        var offset_ss_aux = null,
                            coords_x_target_ss = null;
                        for (var ii = 1; ii < target_num_ss; ii++) {
                            coords_x_target_ss = coords_x_target_ref;
                            if (direction > 0) {
                                offset_ss_aux = ii * offset_ss_step;
                            } else {
                                offset_ss_aux = ( ss_reg[target_e] - ii) * offset_ss_step;
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
                if (direction > 0 && ss != num_ss || direction < 0 && ss != 1) {
                    // Check if is a special exon (starter or finisher)
                    drawDashedLine(ctx, Math.round(coords_x_start_e), Math.round(coords[1]), Math.round(coords_x_start_e), Math.round(coords[1] + exon_height), 2);
                }

                // splice sites dashed lines in target exon
                if (target_ss != 1 && direction > 0 || target_ss != ss_reg[target_e] && direction < 0) {
                    drawDashedLine(ctx, Math.round(coords_x_target_e), Math.round(coords[1]), Math.round(coords_x_target_e), Math.round(coords[1] + exon_height), 2);
                }

                var junc_h_pos = Math.round(margins[3] * 8 * ss); // Math.round((1 - (Math.abs(coords_x_start_e - coords_x_target_e)/canvas.group_width)) * (canvas.group_height*(1-percentage_exon)));
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

    }

    function renderLsvLegend(canvas) {
        if (canvas.getContext) {
            var i;

            // Render LSV representation from a string text representing the LSV i.e.: s|1,2e1.1,1.2|3
            var ctx = canvas.getContext("2d");

            // Clear previous draw
            ctx.clearRect(0, 0, canvas.width, canvas.height);

            var margins = [0, 0, 1, 1],
                pixel_factor = 1,
                percentage_exon = .2;

            var lsv_data = canvas.getAttribute('data-lsv-string');
            var lsvs = lsv_data.split('|');
            var sourceOrTarget = lsvs[0],
                exons_desc = lsvs[1],
                num_types = lsvs[2];
            var ss_coords = exons_desc.split('e');
            var ssFirstExon = ss_coords[0],
                restExons = ss_coords[1];

            restExons = restExons.split(',');
            var num_exons = parseInt(restExons[restExons.length - 1].split('.')[0]) + 1;
            var area = [canvas.width * 1.08 - margins[0] - margins[1], canvas.height - margins[2] - margins[3]];

            var exon_width = area[0] / (num_exons * 2 - 1);
            var start = -exon_width / 2;// margins[0];
            var direction = sourceOrTarget === 's' ? 1 : -1;

            // Render exons_obj
            var exons = [];
            for (i = 0; i < num_exons; i++) {
                var exon = {
                    'coords': [Math.round(start), Math.round(start + exon_width)],
                    'type': (direction > 0 && i === 0 || direction < 0 && i === num_exons - 1 ? 1 : 0)
                };
                exons.push(exon);
                render_exon(canvas, exon, pixel_factor, margins, percentage_exon, '');
                start += exon_width * 2;
            }

            // Render floating box
            var box_length = area[0] * 0.05,
                box_coords = [
                    (sourceOrTarget === 's' ? Math.round(margins[0] + exon_width - box_length / 2) : Math.round(-exon_width + margins[0] + exon_width * ((num_exons * 2 - 1) - 1) - box_length / 2)),
                    2.2 * percentage_exon * canvas.height
                ];

            ctx.strokeStyle = "rgba(0, 0, 0, 1)";
            ctx.fillStyle = "rgba(0, 0, 0, 0.65)";
            ctx.font = "8pt Arial";
            ctx.strokeRect(box_coords[0], box_coords[1], box_length, -box_length);
            ctx.fillRect(box_coords[0], box_coords[1], box_length, -box_length);
            renderNumReads(ctx, box_coords[0] + box_length / 2, box_coords[1] - box_length - 2, num_types);

            // Render junctions to black box
            var list_ss = ssFirstExon.split(',');
            var index_first = direction > 0 ? 0 : exons.length - 1;
            var coords = exons[index_first].coords.slice(0);
            var exon_height = percentage_exon * canvas.height;
            coords[0] += exon_width / 2 + direction * exon_width / 4;
            coords[1] = canvas.height - exon_height - margins[3];

            for (i = 0; i < list_ss.length; i++) {
                if (list_ss.length === 1) {
                    coords[0] += direction * (exon_width / 4 );
                }
                drawLine(ctx, Math.round(coords[0]), Math.round(coords[1]), Math.round(box_coords[0] + box_length / 2 - (direction * box_length / 2)), Math.round(box_coords[1]));
                if (i < list_ss.length - 1) {
                    drawDashedLine(ctx, Math.round(coords[0]), Math.round(coords[1]), Math.round(coords[0]), Math.round(coords[1] + exon_height), 2);
                }
                coords[0] += direction * ((exon_width / 4 ) / (list_ss.length - 1));
            }

            // Render junctions from black box
            for (i = 0; i < restExons.length; i++) {
                var index_ss = restExons[i].split('.');
                var x_dest = exons[index_first].coords[(direction > 0 ? 0 : 1)];
                x_dest += direction * index_ss[0] * 2 * exon_width;
                x_dest += direction * (exon_width / 8 * (index_ss[1] - 1));

                drawLine(ctx, Math.round(box_coords[0] + (direction > 0 ? box_length : 0)), Math.round(box_coords[1]), x_dest, Math.round(coords[1]));
                if (index_ss[1] > 1) {
                    drawDashedLine(ctx, x_dest, Math.round(coords[1]), x_dest, Math.round(coords[1] + exon_height), 2);
                }

            }
        }

    }


    /**
     * Main - Beginning of the execution
     * */


        // add sortable functionality to the table
    var exon_table = $('.exon_table');
    if (exon_table.length) {
        exon_table.tablesorter();
    }


    var canvases = $('.spliceGraph');

    function renderSpliceGraphZoomedPopUp(canvas) {
        /*
         * Adding extraPlots functionality
         */
        $(canvas).on("click", {exon_list: canvas.getAttribute('data-exon-list')}, function (e) {
            e.preventDefault();
            var my_window = window.open("", "Zoomed Splice Graph", "status=1,group_width=2000,group_height=420,scrollbars=yes", true);
            my_window.document.writeln(
                '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">' +
                '<html>' +
                '<head> ' +
                '<script type="text/javascript">' +
                '    CanvasRenderingContext2D.prototype.dashedLine = function (x1, y1, x2, y2, dashLen) {' +
                '        if (dashLen == undefined) dashLen = 2;' +
                '            this.moveTo(x1, y1);' +
                '            var dX = x2 - x1;' +
                '            var dY = y2 - y1;' +
                '            var dashes = Math.floor(Math.sqrt(dX * dX + dY * dY) / dashLen);' +
                '            var dashX = dX / dashes;' +
                '            var dashY = dY / dashes;' +
                '            var q = 0;' +
                '            while (q++ < dashes) {' +
                '                x1 += dashX;' +
                '                y1 += dashY;' +
                '                this[q % 2 == 0 ? \'moveTo\' : \'lineTo\'](x1, y1);' +
                '            }' +
                '            this[q % 2 == 0 ? \'moveTo\' : \'lineTo\'](x2, y2);' +
                '    };' +
                '</script>' +
                '<title>' + $(this)[0].id + '</title>' +
                '</head>' +
                '<body>' +
                '<canvas id="spliceGraph1" class="spliceGraph" group_width="4000px" group_height="400px" data-exon-list="' + e.data.exon_list + '">' +
                'This browser or document mode doesn\'t support canvas' +
                '</canvas>' +
                '</body>' +
                '</html>'
            );
            my_window.document.close();
            my_window.focus();
            renderSpliceGraph($('.spliceGraph', my_window.document.documentElement)[0]);

        });

    }


    return {
        renderLsvLegend: function (canvas) {
            return renderLsvLegend(canvas);
        },
        renderLsvSpliceGraph: function (canvas, gene) {
            return renderLsvSpliceGraph(canvas, gene);
        },
        renderSpliceGraph: function (canvas) {
            var MAX_INTRON = 300;   // Global
            var MAX_EXON = 300;   // Global
            return renderSpliceGraph(canvas);
        },
        renderSpliceGraphZoomedPopUp: function (canvas) {
            return renderSpliceGraphZoomedPopUp(canvas);
        },
        renderFloatingLegend: function (canvas) {
            return renderFloatingLegend(canvas);
        }
    }
};
