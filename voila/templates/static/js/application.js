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

function moveSpliceGraphDataToLsv(spliceGraphs, lsv) {
    var lsvStart = lsv.exons[0].start;
    var lsvEnd = lsv.exons[lsv.exons.length - 1].end;
    lsv.splice_graphs = [];

    delete lsv.exons;
    delete lsv.junctions;

    spliceGraphs.forEach(function (spliceGraph) {
        var foundJunctions = [];
        var starts = false;
        var ends = false;

        var foundExons = spliceGraph.exons.filter(function (exon) {
            return (exon.start >= lsvStart) && (exon.end <= lsvEnd)
        });

        spliceGraph.junctions.forEach(function (junction, juncIndex) {
            foundExons.forEach(function (exon) {
                if ((junction.start >= exon.start) && (junction.start <= exon.end))
                    starts = true;
                if ((junction.end >= exon.start) && (junction.end <= exon.end))
                    ends = true
            });

            if (starts && ends)
                foundJunctions.push(junction);

            starts = false;
            ends = false
        });

        foundExons.forEach(function (exon) {
            exon.a3 = [];
            exon.a5 = [];
            foundJunctions.forEach(function (junction, juncIndex) {
                if ((junction.start >= exon.start) && (junction.start <= exon.end))
                    exon.a5.push(juncIndex);
                if ((junction.end >= exon.start) && (junction.end <= exon.end))
                    exon.a3.push(juncIndex)
            })
        });

        lsv.splice_graphs.push({'exons': foundExons, 'junctions': foundJunctions})

    })


}

$(document).ready(function () {
    new Clipboard('.copy-to-clipboard', {
        text: function (trigger) {
            var lsv = $(trigger).parents('tr').find('.primer').attr('data-lsv');
            lsv = JSON.parse(lsv.replace(/'/g, '"'));

            var spliceDivs = $(trigger).parents('.gene-container').find('.spliceDiv').get();
            var spliceGraphs = spliceDivs.reduce(function (accu, currVal) {
                accu.push(JSON.parse(currVal.getAttribute('data-exon-list').replace(/'/g, '"')));
                return accu
            }, []);

            moveSpliceGraphDataToLsv(spliceGraphs, lsv);
            return JSON.stringify(lsv)
        }
    });

    // add sortable functionality to the table
    window.gene_objs = [];
    window.gene_obj_list = [];

    var $sgFilters = $('form#sg-filters');

    $("#simplified").on('change', function () {
        $sgFilters.slideToggle(100).submit();
    });

    /**
     * Process splice graph filters
     */
    $sgFilters.on('keyup submit', function (event) {
        // prevent form submit
        event.preventDefault();

        // variables
        var numReads = this.querySelector("[name='numReads']");
        var numReadsValue = parseInt(numReads.value);
        var simplified = document.querySelector("[name='simplified']").checked;
        var currentNumReads;

        // filter elements with these four classes
        ['.junction', '.readcounts', '.irlines', '.irreads'].forEach(function (value) {
            d3.selectAll(value).classed('sgfilter', function (d) {
                if (!simplified) return false;

                // check read counts toggle button for read counts state
                // var toggleReadCounts = this.parentNode.parentNode.parentNode.parentNode.querySelector('.readCounts');

                // For now.. we're going to comment this out and add per LSV clean reads after this release.
                // var displayNormReads = !toggleReadCounts.checked;
                var displayNormReads = true;

                // get read counts based on this specific splice graphs state
                var cleanReads = d.clean_reads ? d.clean_reads : 0;
                currentNumReads = displayNormReads ? d.reads : cleanReads;

                // return if this element should have the 'sgfilter' class
                return !isNaN(numReadsValue) && currentNumReads <= numReadsValue;
            });
        })
    });

    $('.floatingLegend').each(function () {
        if ($(this)[0].getContext) {
            splicegraph().renderFloatingLegend($(this)[0]);
        }
    });

    $('.lsvSelect').on('change', function () {
        splicegraph().selectLSV(this, this);
    });

    $('.tablesorter').each(function () {
        $(this).tablesorter({
            sortList: [
                [0, 0]
            ]
        }); // Disable sort function in column PDF  , headers: {3: {sorter: false}, 4: {sorter: false}, 5: {sorter: false}}
        $(this).tablesorterPager({
            widthFixed: true,
            widgets: ['zebra', 'renderCanvas'],
            container: $(this).parent().children(".pager")
        });
    });

    $('.lsvLegendThumb').each(function () {
        var collapsed = this.getAttribute('data-collapsed');
        if (collapsed != 0) {
            splicegraph().renderLsvLegend(this);
        } else {
            splicegraph().renderLsvSpliceGraph(this);
        }
        var can = this;

        function dlCanvas() {
            this.href = can.toDataURL('image/png');
        }

        var dl_canvas_link = $(this).parent().children(".lsv_type")[0];
        dl_canvas_link.addEventListener('click', dlCanvas, false);


    });

    /** Tooltip for barchart */
    var tooltips = $('.tooltip');
    if (tooltips.length) {
        tooltips.tooltipster({
            theme: 'tooltipster-light'
        });
    }

    /**
     * LSV Filters functionality
     * */
    var lsvFilters = $('#lsv-filters');
    if (lsvFilters.length) {
        lsvFilters.on('focusin', ":input", function () {
            $(this).data('oldValue', $(this).val());
        });
        lsvFilters.on("change", ":input", function () {
            var eventFired = $(this)[0];
            if (eventFired.name in {'ES': null, 'prime3': null, 'prime5': null, 'source': null, 'target': null}) {
                $('.' + eventFired.name).each(function () {
                    $(this).toggleClass(eventFired.name + 'hidden');
                });
            }
            else {
                var nval = parseInt(eventFired.value);
                if (nval > -1) {
                    $('.lsvrow').each(function () {
                        if (eventFired.name.indexOf('from') > -1) {
                            if (parseInt($(this)[0].getAttribute("data-" + eventFired.name.split('from')[0])) < nval) {
                                $(this).addClass(eventFired.name + 'hidden');
                            } else {
                                $(this).removeClass(eventFired.name + 'hidden');
                            }
                        }
                        if (eventFired.name.indexOf('to') > -1) {
                            if (parseInt($(this)[0].getAttribute("data-" + eventFired.name.split('to')[0])) > nval) {
                                $(this).addClass(eventFired.name + 'hidden');
                            } else {
                                $(this).removeClass(eventFired.name + 'hidden');
                            }
                        }
                    });
                } else {
                    eventFired.value = $(this).data('oldValue');
                }
            }
        });
    }

});

var initLargeCanvasSettings = function (num_bins, canvas) {

    // Calculate canvas drawable settings
    var settingsCanvas;

    return function () {
        if (settingsCanvas) {
            return settingsCanvas;
        }
        settingsCanvas = {};
        settingsCanvas.margin_y = [10, 25];
        settingsCanvas.margin_x = [30, 1];
        settingsCanvas.area_pixels = [canvas.width - settingsCanvas.margin_x.reduce(add, 0),
            canvas.height - settingsCanvas.margin_y.reduce(add, 0)];
        settingsCanvas.coords_origin = [settingsCanvas.margin_x[0], settingsCanvas.margin_y[0] + settingsCanvas.area_pixels[1]];
        settingsCanvas.bar_width = settingsCanvas.area_pixels[0] / num_bins;
        settingsCanvas.labels_steps = [Math.ceil((num_bins + 1) / 8), 4];
        return settingsCanvas;
    }();

};

function drawLine(contextO, startx, starty, endx, endy) {
    contextO.beginPath();
    contextO.moveTo(startx, starty);
    contextO.lineTo(endx, endy);
    contextO.closePath();
    contextO.stroke();
}

function drawDashedLine(contextO, startx, starty, endx, endy, dashLen) {
    contextO.beginPath();
    contextO.dashedLine(startx, starty, endx, endy, dashLen);
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

function roundRect(ctx, x, y, width, height, radius, fill, stroke) {
    if (typeof stroke == "undefined") {
        stroke = true;
    }
    if (typeof radius === "undefined") {
        radius = 5;
    }
    ctx.beginPath();
    ctx.moveTo(x + radius, y);
    ctx.lineTo(x + width - radius, y);
    ctx.quadraticCurveTo(x + width, y, x + width, y + radius);
    ctx.lineTo(x + width, y + height - radius);
    ctx.quadraticCurveTo(x + width, y + height, x + width - radius, y + height);
    ctx.lineTo(x + radius, y + height);
    ctx.quadraticCurveTo(x, y + height, x, y + height - radius);
    ctx.lineTo(x, y + radius);
    ctx.quadraticCurveTo(x, y, x + radius, y);
    ctx.closePath();
    if (stroke) {
        ctx.stroke();
    }
    if (fill) {
        ctx.fill();
    }
}


function createGradientSinglePlots(ctx, canvas, margins) {
    //create a gradient object from the canvas context
    var gradient = ctx.createLinearGradient(0, 0, canvas.width, 0);

    // Add the colors with fixed stops at 1/4 of the width.
//    gradient.addColorStop(0, "#0000FF");
//    gradient.addColorStop(.45, "cyan");
//    gradient.addColorStop(.50, "#00FF00");
//    gradient.addColorStop(.65, "yellow");
//    gradient.addColorStop(1, "#FF0000");
    if (margins) {
        gradient = ctx.createLinearGradient(margins[0], 0, canvas.width - margins[1], 0);
    }

    gradient.addColorStop(1, rgbToHex(213, 62, 79));
    gradient.addColorStop(.875, rgbToHex(244, 109, 67));
    gradient.addColorStop(.75, rgbToHex(253, 174, 97));
    gradient.addColorStop(.625, rgbToHex(254, 224, 139));
    gradient.addColorStop(.5, rgbToHex(255, 255, 191));
    gradient.addColorStop(.375, rgbToHex(230, 245, 152));
    gradient.addColorStop(.25, rgbToHex(171, 221, 164));
    gradient.addColorStop(.125, rgbToHex(102, 194, 165));
    gradient.addColorStop(0, rgbToHex(50, 136, 189));


    return gradient;
}


function drawBarchartWithCanvasId(eventOrCanvasid, zoomInc, canvasSettings) {
    // To cover the case when the reset button is clicked and when we just render the canvas passing a canvas id

    var canvasId;
    if (eventOrCanvasid.type == "click") {
        // This method was fired by a click event! We need to retrieve the params...
        zoomInc = eventOrCanvasid.data.zoomInc;
        canvasId = eventOrCanvasid.data.canvasId;
    } else {
        canvasId = eventOrCanvasid;
    }
    var canvas = $("#".concat(canvasId))[0];
    renderLargeCanvas(canvas, zoomInc, canvasSettings);
}

function attachResetZoomLink(canvas, drawFunction, canvasSettingsP) {
    var idResetImage = 'resetZoom'.concat(canvas.getAttribute('id'));
    //find the td in which we need to append the reset button
    $(canvas).parent().prepend('<img id="'.concat(idResetImage).concat('" src="static/img/undo.png" class="hidden resetZoom tooltip" title="Reset zoom" />'));
    $("#".concat(idResetImage)).on("click", {
        canvasId: canvas.id,
        zoomInc: 0,
        canvasSettings: canvasSettingsP
    }, drawFunction);
}


function drawBarChart(context, bins, settingsCanvasSingle, zoomLevel) {

    // Draw the x and y axes
    context.lineWidth = "1.0";

    // Y-axis
    drawLine(context,
        settingsCanvasSingle.coords_origin[0],
        settingsCanvasSingle.coords_origin[1] + 5, //Mark axis
        settingsCanvasSingle.coords_origin[0],
        settingsCanvasSingle.coords_origin[1] - settingsCanvasSingle.area_pixels[1]);
    // X-axis
    drawLine(context,
        settingsCanvasSingle.coords_origin[0], //Mark axis
        settingsCanvasSingle.coords_origin[1],
        settingsCanvasSingle.coords_origin[0] + settingsCanvasSingle.area_pixels[0],
        settingsCanvasSingle.coords_origin[1]);
//    drawLine(context, settingsCanvasSingle.coords_origin[0], canvasHeight-marginYaxis, canvasWidth, canvasHeight-marginYaxis);
    context.lineWidth = "0.0";

    var chartHeight = settingsCanvasSingle.area_pixels[1];
    chartHeight += zoomLevel;
//    var maxValue = chartHeight;

    //Preserve the gradient color style
    var gradientColorStyle = context.fillStyle;

    // Offset from origin to write the labels in x-axis
    var xLabelsCoordX = Math.min(15, settingsCanvasSingle.margin_y[1]);

    for (var i = 0; i < bins.length; i++) {
        var name = i / bins.length;
        var height = parseInt(bins[i] * chartHeight);
        //if (parseInt(height) > parseInt(maxValue)) height = maxValue;

        // Write the data to the chart
        //context.fillStyle = "#b90000";
        drawRectangle(context,
            settingsCanvasSingle.coords_origin[0] + (i * settingsCanvasSingle.bar_width),
            settingsCanvasSingle.coords_origin[1] - height,
            settingsCanvasSingle.bar_width,
            height,
            true);

        // Add the column title to the x-axis
        context.textAlign = "center";
        context.fillStyle = "black";
        name = name * 100;
        if (i % settingsCanvasSingle.labels_steps[0] == 0) {
            name = name.toFixed(1);
            context.fillText(name,
                settingsCanvasSingle.coords_origin[0] + (i * settingsCanvasSingle.bar_width),
                settingsCanvasSingle.coords_origin[1] + xLabelsCoordX
            );
            drawLine(context,
                settingsCanvasSingle.coords_origin[0] + (i * settingsCanvasSingle.bar_width),
                settingsCanvasSingle.coords_origin[1] + 1, // +1 to don't overlap
                settingsCanvasSingle.coords_origin[0] + (i * settingsCanvasSingle.bar_width),
                settingsCanvasSingle.coords_origin[1] + 4 // size of the delimiter=4
            )
        }

        context.fillStyle = gradientColorStyle;
    }

    // Add some data markers to the y-axis
    var convertedMarkDataIncrementsIn = settingsCanvasSingle.area_pixels[1] / settingsCanvasSingle.labels_steps[1]; // chartHeight/settingsCanvasSingle.labels_steps[1];

    context.textAlign = "right";
    context.fillStyle = "black";
    var markerValue = 0;
    var offset = 0;
    while (offset <= settingsCanvasSingle.area_pixels[1]) {
        context.fillText(markerValue.toFixed(1), settingsCanvasSingle.coords_origin[0] - 2, settingsCanvasSingle.coords_origin[1] - offset, 50);
        markerValue += (100 / settingsCanvasSingle.labels_steps[1]) * (settingsCanvasSingle.area_pixels[1] / chartHeight);
        offset += convertedMarkDataIncrementsIn;
    }
}

function renderLargeCanvas(canvas, zoomInc) {
    if (canvas.getContext) {  // check for support
        var ctx = canvas.getContext("2d");

        // Clear canvas
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        // Create gradient for coloured barchart
        ctx.fillStyle = createGradientSinglePlots(ctx, canvas);

        var bins = JSON.parse(canvas.getAttribute('bins'));

        // Compute the zoom
        canvas.zoom = Math.max(canvas.zoom + zoomInc * 15, 0); // Only positive zoom
        canvas.zoom = Math.abs(zoomInc) * canvas.zoom; // To reset zoom
        (canvas.zoom == 0) ? hideResetZoomLink(canvas) : showResetZoomLink(canvas);

        drawBarChart(ctx, bins, window.settingsCanvasSingle, canvas.zoom); // Canvas settings shared by all (in window)
    }

}

function showResetZoomLink(canvas) {
    $(canvas).parent().children('.resetZoom').removeClass('hidden');
}

function hideResetZoomLink(canvas) {
    $(canvas).parent().children('.resetZoom').addClass('hidden');
}
// Get Color from Brewer Palette
BREWER_PALETTE = [
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


function getColor(colorNumber, palette, hue) {
    colorNumber = colorNumber % 16;
    return "rgba(" + palette[colorNumber].toString() + ", " + hue + ")";
//        return rgbToHex(palette[colorNumber][0], palette[colorNumber][1], palette[colorNumber][2]);
}


function drawLSVCompactStackBars(canvas, fillMode) {

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
        // Clear canvas
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        var groups_str = '[' + canvas.getAttribute("data-lsv") + ']';
        var groups = JSON.parse(groups_str.replace(/'/g, '"'));

        // Calculate origins_coords
        var header_height = 0; // canvas.height*.1;
        var num_groups = groups.length;

        var sub_canvas_w = canvas.width / num_groups;
        var sub_canvas_h = canvas.height - header_height;
//        var sub_canvas_margins = [sub_canvas_w *.00, 0, header_height, sub_canvas_h *.21];
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

            for (var lsv_count = 0; lsv_count < group.bins.length; lsv_count++) {
                // Calculate the height of the accumulated mean
                acc_height += group.means[lsv_count];

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
                    origins_coords[count][1] - (acc_height - group.means[lsv_count]) * sub_canvas_pixels[1],
                    origins_coords[count][0] + sub_canvas_pixels[0],
                    origins_coords[count][1] - (acc_height) * sub_canvas_pixels[1]
                );
            }
        }
    }
}

function drawDeltaLSVCompactSVG(htmlElementId, lsv, threshold) {
    var width = 200,
        height = 20;
    var margin = {top: 1, bottom: 8, left: 2, right: 2};
    var domain = [-1, 1];
    var x = d3.scale.linear()
        .range([margin.left, margin.right])
        .domain(domain);
    var border_frame = 2,
        MIN_DELTAPSI = .05;

    var svgContainer = d3.select("#" + htmlElementId)
        .append("svg")
        .attr("class", "excl-incl-rect")
        .attr("width", width)
        .attr("height", height);

    var markerWidth = 6,
        markerHeight = 6,
        cRadius = 30, // play with the cRadius value
        refX = cRadius + (markerWidth * 2),
        refY = -Math.sqrt(cRadius);

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

}


function drawDeltaBox(canvas) {
    if (canvas.getContext) {  // check for support
        var ctx = canvas.getContext("2d");

        // Clear canvas
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        var excl_incl_set = JSON.parse(canvas.dataset.exclIncl);
        var halfCanvasWidth = parseInt(canvas.width / 2);
        var margin_top = 2;
        var margin_bottom = 10;
        // Create gradient exclusion/inclusion

        ctx.fillStyle = "grey";
        ctx.lineWidth = 0.5;
        drawRectangle(ctx, 1, 1, canvas.width - 2, canvas.height - margin_bottom, 0); // Canvas box
        // Fill with gradient
//        ctx.fillStyle=createGradientDeltaPlots(ctx, 0, canvas.width);

        // draw exclusion/inclusion Box
        ctx.fillStyle = "blue";
        ctx.fillStyle = rgbToHex(55, 126, 184); // blue
        var set18qual1 = "rgba(228,26,28)";
        var set18qual2 = "rgba(55,126,184)";
        ctx.fillRect(halfCanvasWidth, margin_top, -parseInt(excl_incl_set[0] * canvas.width / 2), canvas.height - margin_top - margin_bottom);
        ctx.fillStyle = "red";
        ctx.fillStyle = rgbToHex(228, 26, 28); // red
        ctx.fillRect(halfCanvasWidth, margin_top, parseInt(excl_incl_set[1] * canvas.width / 2), canvas.height - margin_top - margin_bottom);


        // draw -1, 0 and 1
        ctx.textAlign = "center";
        ctx.font = "bold 8pt Arial";
        ctx.fillStyle = "blue";
        ctx.fillStyle = rgbToHex(55, 126, 184); // blue
        ctx.fillText("-" + (excl_incl_set[0] * 100).toFixed(1) + " ", Math.max(halfCanvasWidth - Math.max(parseInt(excl_incl_set[0] * canvas.width / 2), 10), 20), canvas.height);
        ctx.fillStyle = "red";
        ctx.fillStyle = rgbToHex(228, 26, 28); // red
//        ctx.fillText("1", canvas.width - margin_top, canvas.height);
        ctx.fillText("+" + (excl_incl_set[1] * 100).toFixed(1), Math.min(halfCanvasWidth + Math.max(parseInt(excl_incl_set[1] * canvas.width / 2), 10), halfCanvasWidth * 2 - 20), canvas.height);
        ctx.fillStyle = "black";
//        ctx.fillText("0", halfCanvasWidth, canvas.height);


        // draw middle bar (0)
        ctx.lineWidth = 1.5;
        ctx.strokeStyle = "black";
        drawLine(ctx, halfCanvasWidth, parseInt(margin_top / 2), halfCanvasWidth, canvas.height - margin_bottom + margin_top / 2);
    }
}

function drawDeltaBarChart(context, binsArray, settingsCanvasDelta, zoomLevel, pThreshold, psi_junc) {
    // pThreshold is gonna be used to differentiate between single and delta PSI

    // Draw the x and y axes
    context.lineWidth = "1.0";

    var chartHeight = settingsCanvasDelta.area_pixels[1];
    chartHeight += zoomLevel;
    var maxValue = 0;

    for (var binsCount = 0; binsCount < 1; binsCount++) { //binsArray.length
        var bins = binsArray[binsCount];
        //Preserve the gradient color style
//        var gradientColorStyle = createGradientDeltaPlots(context, settingsCanvasDelta.coords_origin[0],
//                settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0], binsCount, binsCount, pThreshold);

        context.fillStyle = getColor(psi_junc, BREWER_PALETTE, .5); //gradientColorStyle;
        context.strokeStyle = getColor(psi_junc, BREWER_PALETTE, 1);

        // Offset from origin to write the labels in x-axis
        //    var xLabelsCoordX = Math.min(15, settingsCanvasDelta.margin_y[1]);

        for (var i = 0; i < bins.length; i++) {
            var name = i / bins.length;
            var height = Math.round(bins[i] * chartHeight);
            if (height > maxValue) maxValue = height;

            // Write the data to the chart
            drawRectangle(context,
                Math.round(settingsCanvasDelta.coords_origin[0] + (i * settingsCanvasDelta.bar_width)),
                Math.round(settingsCanvasDelta.coords_origin[1] - height),
                settingsCanvasDelta.bar_width,
                height,
                true);

        }

    } // End Vertical bars

    // Add the column title to the x-axis
    context.textAlign = "center";
    context.fillStyle = "black";
    context.strokeStyle = "black";

    // X-axis
    drawLine(context,
        settingsCanvasDelta.coords_origin[0], //Mark axis
        settingsCanvasDelta.coords_origin[1],
        settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0],
        settingsCanvasDelta.coords_origin[1]);

    // Y-axis
    var yaxis_x = settingsCanvasDelta.coords_origin[0],
        marks = [0, 0.5, 1];
    if (pThreshold) {
        yaxis_x += (settingsCanvasDelta.area_pixels[0]) / 2;
        marks = [-1, 0, 1];
    }

    drawLine(context,
        yaxis_x,
        settingsCanvasDelta.coords_origin[1] + 5, //Mark axis
        yaxis_x,
        settingsCanvasDelta.coords_origin[1] - settingsCanvasDelta.area_pixels[1]);
    name *= 100;
    context.fillText(marks[0], settingsCanvasDelta.coords_origin[0], settingsCanvasDelta.coords_origin[1] + 15);
    drawLine(context, settingsCanvasDelta.coords_origin[0],
        settingsCanvasDelta.coords_origin[1],
        settingsCanvasDelta.coords_origin[0],
        settingsCanvasDelta.coords_origin[1] + 2 // size of the delimiter=4
    );

    context.fillText(marks[1], settingsCanvasDelta.coords_origin[0] + (settingsCanvasDelta.area_pixels[0]) / 2, settingsCanvasDelta.coords_origin[1] + 15);
    drawLine(context, settingsCanvasDelta.coords_origin[0] + (settingsCanvasDelta.area_pixels[0]) / 2,
        settingsCanvasDelta.coords_origin[1],
        settingsCanvasDelta.coords_origin[0] + (settingsCanvasDelta.area_pixels[0]) / 2,
        settingsCanvasDelta.coords_origin[1] + 2 // size of the delimiter=4
    );

    context.textAlign = "right";
    context.fillText(marks[2], settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0], settingsCanvasDelta.coords_origin[1] + 15);
    drawLine(context, settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0],
        settingsCanvasDelta.coords_origin[1], // +1 to don't overlap
        settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0],
        settingsCanvasDelta.coords_origin[1] + 2 // size of the delimiter=4
    );

    // Vertical bars for percentages of inclusion/exclusion (by Default, in -20 and +20)
    if (pThreshold) {
        context.textAlign = "center";
        context.dashedLine(
            settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0] / 2 - parseInt(settingsCanvasDelta.area_pixels[0] / 2 * pThreshold),
            settingsCanvasDelta.coords_origin[1] + 5, //Mark axis
            settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0] / 2 - parseInt(settingsCanvasDelta.area_pixels[0] / 2 * pThreshold),
            settingsCanvasDelta.coords_origin[1] - settingsCanvasDelta.area_pixels[1]
        );
        context.stroke();
        context.fillText("-" + pThreshold + " ",
            settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0] / 2 - parseInt(settingsCanvasDelta.area_pixels[0] / 2 * pThreshold),
            settingsCanvasDelta.coords_origin[1] + 15 //Mark axis
        );

        context.dashedLine(
            settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0] / 2 + parseInt(settingsCanvasDelta.area_pixels[0] / 2 * pThreshold),
            settingsCanvasDelta.coords_origin[1] + 5, //Mark axis
            settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0] / 2 + parseInt(settingsCanvasDelta.area_pixels[0] / 2 * pThreshold),
            settingsCanvasDelta.coords_origin[1] - settingsCanvasDelta.area_pixels[1]
        );
        context.stroke();
        context.fillText(pThreshold,
            settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0] / 2 + parseInt(settingsCanvasDelta.area_pixels[0] / 2 * pThreshold),
            settingsCanvasDelta.coords_origin[1] + 15 //Mark axis
        );

    }
    // Add some data markers to the y-axis
    var convertedMarkDataIncrementsIn = settingsCanvasDelta.area_pixels[1] / settingsCanvasDelta.labels_steps[1]; // chartHeight/settingsCanvasDelta.labels_steps[1];

    context.textAlign = "right";
    context.fillStyle = "black";
    var markerValue = 0;
    var offset = 0;
    while (offset <= settingsCanvasDelta.area_pixels[1]) {
        context.fillText(markerValue.toFixed(1), settingsCanvasDelta.coords_origin[0] - 2, settingsCanvasDelta.coords_origin[1] - offset, 50);
        markerValue += (100 / settingsCanvasDelta.labels_steps[1]) * (settingsCanvasDelta.area_pixels[1] / chartHeight);
        offset += convertedMarkDataIncrementsIn;
    }


}

function drawExpDeltaWithCanvasId(eventOrCanvasid, zoomInc, canvasSettings) {
    // To cover the case when the reset button is clicked and when we just render the canvas passing a canvas id

    var canvasId;
    if (eventOrCanvasid.type == "click") {
        // This method was fired by a click event! We need to retrieve the params...
        zoomInc = eventOrCanvasid.data.zoomInc;
        canvasId = eventOrCanvasid.data.canvasId;
        canvasSettings = eventOrCanvasid.data.canvasSettings;
    } else {
        canvasId = eventOrCanvasid;
    }
    canvas = $("#".concat(canvasId))[0];
    renderExpandedDeltaCanvas(canvas, zoomInc, canvasSettings);
}

function renderExpandedDeltaCanvas(canvas, zoomInc, canvasSettings) {
    if (canvas.getContext) {  // check for support
        var ctx = canvas.getContext("2d");

        // Clear canvas
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        // Create gradient for coloured barchart
//        ctx.fillStyle = createGradientSinglePlots(ctx, canvas);

        var lsv = JSON.parse(canvas.getAttribute('data-lsv'));
//        var bins = lsv.bins;

        // Compute the zoom
        var zoom = parseInt(canvas.getAttribute('data-zoom'));
        zoom = Math.max(zoom + zoomInc * 15, 0); // Only positive zoom
        zoom = Math.abs(zoomInc) * zoom; // To reset zoom

//        canvas.zoom = Math.abs(zoomInc)*canvas.zoom; // To reset zoom
        (zoom == 0) ? hideResetZoomLink(canvas) : showResetZoomLink(canvas);

        // Canvas attributes to reset barchart from zoom view
        canvas.setAttribute('data-zoom', zoom.toString());

        var pThreshold = canvas.getAttribute('data-threshold');
        drawDeltaBarChart(ctx, lsv.bins, canvasSettings, zoom, pThreshold, lsv.psi_junction); // Canvas settings shared by all (in window)
    }

}

function initExpandedDeltaCanvas(canvas, canvasSettings) {
    canvas.setAttribute('data-zoom', '0');
    attachResetZoomLink(canvas, drawExpDeltaWithCanvasId, canvasSettings);
    renderExpandedDeltaCanvas(canvas, 0, canvasSettings);
}

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

function mul(a, b) {
    return a * b;
}

function sum(numbers) {
    return reduce(add, 0, numbers);
}

function map(func, array) {
    var result = [];
    forEach(array, function (element) {
        result.push(func(element));
    });
    return result;
}

function componentToHex(c) {
    var hex = c.toString(16);
    return hex.length == 1 ? "0" + hex : hex;
}

function rgbToHex(r, g, b) {
    return "#" + componentToHex(r) + componentToHex(g) + componentToHex(b);
}

function renderViolin(htmlElementId, results, tableId, params) {


    function addViolin(svg, results, height, width, domain, imposeMax, count, id_svg, id_table) {

        var data = d3.layout.histogram()
            .bins(params.num_bins)
            .range(domain)
            (results);


        var y = d3.scale.linear()
            .range([width / 2, 0])
            .domain([0, Math.max(imposeMax, d3.max(data, function (d) {
                return d.y;
            }))]);

        var x = d3.scale.linear()
            .range([height - margin.bottom, margin.top - 3]) //-margin.left, margin.right])
            .domain(domain)
            .nice();

        var area = d3.svg.area()
            .interpolate(interpolation)
            .x(function (d) {
                if (interpolation == "step-before")
                    return x(d.x + d.dx / 2);
                return x(d.x);
            })
            .y0(width / 2)
            .y1(function (d) {
                return y(d.y);
            });

        var line = d3.svg.line()
            .interpolate(interpolation)
            .x(function (d) {
                if (interpolation == "step-before")
                    return x(d.x + d.dx / 2);
                return x(d.x);
            })
            .y(function (d) {
                return y(d.y);
            });

        svg.append("linearGradient")
            .attr("id", "violin-gradient" + id_svg + count + id_table)
            .attr("gradientUnits", "userSpaceOnUse")
            .attr("x1", margin.top).attr("y1", 0)
            .attr("x2", height - margin.bottom).attr("y2", 0)
            .selectAll("stop")
            .data([
                {offset: "0%", color: getColor(count, BREWER_PALETTE, 1)},
                {offset: "100%", color: getColor(count, BREWER_PALETTE, 1)} //"steelblue" "gray"
            ])
            .enter().append("stop")
            .attr("offset", function (d) {
                return d.offset;
            })
            .attr("stop-color", function (d) {
                return d.color;
            });


        var gPlus = svg.append("g");
        var gMinus = svg.append("g");

        gPlus.append("path")
            .datum(data)
            .attr("class", "area")
            .attr("d", area);

        gMinus.append("path")
            .datum(data)
            .attr("class", "area")
            .attr("d", area);

        gPlus.attr("transform", "rotate(90,0,0)  translate(0,-" + width + ")");
        gMinus.attr("transform", "rotate(90,0,0) scale(1,-1)");

        gPlus.attr('fill', 'url(#violin-gradient' + id_svg + count + id_table + ')');
        gMinus.attr('fill', 'url(#violin-gradient' + id_svg + count + id_table + ')');


    }

    function addBoxPlot(svg, results, height, width, domain, boxPlotWidth) {
        var y = d3.scale.linear()
            .range([height - margin.bottom, margin.top])
            .domain(domain);

        var x = d3.scale.linear()
            .range([0, width]);

        var left = 0.5 - boxPlotWidth / 2;
        var right = 0.5 + boxPlotWidth / 2;

        var probs = [0.05, 0.25, 0.5, 0.75, 0.95];
        for (var i = 0; i < probs.length; i++) {
            probs[i] = y(d3.quantile(results, probs[i]));
        }

        svg.append("rect")
            .attr("class", "boxplot fill")
            .attr("x", x(left))
            .attr("width", x(right) - x(left))
            .attr("y", probs[3])
            .attr("height", -probs[3] + probs[1]);
//            .attr("y", probs[3])
//            .attr("height", probs[1]);

        svg.append("circle")
            .attr("class", "boxplot mean")
            .attr("cx", x(0.5))
            .attr("cy", y(d3.mean(results)))
            .attr("r", x(boxPlotWidth / 5));

        var iS = [0, 2, 4];
        var iSclass = ["", "median", ""];
        for (var i = 0; i < iS.length; i++) {
            svg.append("line")
                .attr("class", "boxplot " + iSclass[i])
                .attr("x1", x(left))
                .attr("x2", x(right))
                .attr("y1", probs[iS[i]])
                .attr("y2", probs[iS[i]]);
        }

        iS = [[0, 1], [3, 4]];
        for (var i = 0; i < iS.length; i++) {
            svg.append("line")
                .attr("class", "boxplot")
                .attr("x1", x(0.5))
                .attr("x2", x(0.5))
                .attr("y1", probs[iS[i][0]])
                .attr("y2", probs[iS[i][1]]);
        }

//        svg.append("rect")
//            .attr("class", "boxplot")
//            .attr("x", x(left))
//            .attr("width", x(right)-x(left))
//            .attr("y", probs[3])
//            .attr("height", -probs[3]+probs[1]);


    }

    function addExpectedPSI(svg, mean_value, height, boxWidth, count) {
        svg.append("text")
            .attr("x", boxWidth / 2)
            .attr("y", height)
            .attr("text-anchor", "middle")
            .attr("font-size", "12px")
            .attr("fill", getColor(count, BREWER_PALETTE, 1))
            .text(mean_value.toFixed(3));
    }


    var margin = {top: 10, bottom: 30, left: 30, right: 10};
    var width = 100 * results.length; //element_jq.getAttribute('width');
    var height = 200; //element_jq.getAttribute('height');
    var spacing_space = (width - margin.left - margin.right) * .05;
    var boxWidth = Math.round(((width - margin.left - margin.right) - spacing_space) / results.length);
    var boxSpacing = Math.round(spacing_space / results.length);

    var domain = [0, 1];
    if (params.delta) {
        domain = [-1, 1];
    }

    var resolution = params.num_bins;
    var interpolation = 'basis'; // 'step-before'; 'basis'

    var y = d3.scale.linear()
        .range([height - margin.bottom, margin.top])
        .domain(domain);

    var yAxis = d3.svg.axis()
        .scale(y)
        .ticks(6)
        .orient("left")
        .tickSize(5, 0, 5);


    var svg = d3.select("#" + htmlElementId)
        .append("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("class", "violin-boxplot");

    svg.append("line")
        .attr("class", "boxplot")
        .attr("x1", margin.left)
        .attr("x2", width - margin.right)
        .attr("y1", y(0))
        .attr("y2", y(0));

    for (var i = 0; i < results.length; i++) {
        var g = svg.append("g").attr("transform", "translate(" + (i * (boxWidth + boxSpacing) + margin.left) + ",0)");
        addViolin(g, results[i], height, boxWidth, domain, 0.25, i, htmlElementId, tableId);
        addBoxPlot(g, results[i], height, boxWidth, domain, .15);
        addExpectedPSI(g, d3.mean(results[i]), height, boxWidth, i);
    }

    svg.append("g")
        .attr('class', 'axis')
        .attr("transform", "translate(" + margin.left + ",0)")
        .call(yAxis);

    return svg[0];
}

function translate_lsv_bins(lsv_bins, num_samples) {
    var adjusted_bins = [];
    for (var lsv_way = 0; lsv_way < lsv_bins.length; lsv_way++) {
        var tmp_bins = [];
        var bins_size = lsv_bins[lsv_way].length;

        for (var ii = 1; ii < bins_size + 1; ii++) {

            var num_copies = Math.round(num_samples * lsv_bins[lsv_way][ii - 1]);

            for (var bins_i = 0; bins_i < num_copies; bins_i++) {
                tmp_bins.push((1 / bins_size) / 2 + ((ii - 1) / bins_size));
            }
        }
        adjusted_bins.push(tmp_bins);
    }
    return adjusted_bins
}


function translate_delta_lsv_bins(lsv_bins, num_samples) {
    var adjusted_bins = [];
    for (var lsv_way = 0; lsv_way < lsv_bins.length; lsv_way++) {
        var tmp_bins = [];
        var bins_size = lsv_bins[lsv_way].length;
        var start_offset = -1 + 1 / bins_size;
        for (var ii = 0; ii < bins_size; ii++) {
            var num_copies = Math.round(num_samples * lsv_bins[lsv_way][ii]);
            for (var bins_i = 0; bins_i < num_copies; bins_i++) {
                tmp_bins.push(start_offset + ii * 2 / bins_size);
            }
        }
        adjusted_bins.push(tmp_bins);
    }
    return adjusted_bins
}
