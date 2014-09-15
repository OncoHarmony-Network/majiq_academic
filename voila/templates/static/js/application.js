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


$( document ).ready(function(){

    // add sortable functionality to the table
    $('.tablesorter').each(function() {
        $(this).tablesorter({sortList: [
            [0, 0]
        ]}); // Disable sort function in column PDF  , headers: {3: {sorter: false}, 4: {sorter: false}, 5: {sorter: false}}
        $(this).tablesorterPager({widthFixed: true, widgets: ['zebra', 'renderCanvas'], container: $(this).parent().children(".pager")});
    });

    // Check if the html has psiPlots and largePsiPlots (render legend if so)
    if ($(document).find('.psiPlot').length){
        // Define canvas settings and store them as a document attribute (global)
        var psiPlot_list = $(document).find('.psiPlot');
        var largePsiPlot_list = $(document).find('.largePsiPlot');

        var binsSingle = JSON.parse(psiPlot_list[0].getAttribute('bins'));
        var canvasSingle = largePsiPlot_list[0];

        window.settingsCanvasSingle = initLargeCanvasSettings(binsSingle.length, canvasSingle);
    }

    if ($(document).find('.extendedDeltaPsi').length){
        // Define canvas settings and store them as a document attribute (global)
        var extendedDeltas = $(document).find('.extendedDeltaPsi');
        var binsDelta = JSON.parse(extendedDeltas[0].getAttribute('data-bins'));
        var canvasDelta = extendedDeltas[0];

        window.settingsCanvasDelta = initLargeCanvasSettings(binsDelta.length, canvasDelta);
    }

    function renderLegendCanvas() {
        var canvas = $('.legendCanvas')[0];
        if (canvas.getContext){
            $(canvas).css({"margin-top": "-6px", "margin-bottom": "-6px"});
            var ctx = canvas.getContext("2d");
            ctx.textAlign = "center";
            // TODO: remove hardcoded values
            ctx.fillText("0", 27, 20);
            ctx.fillText("1", 417, 20);
            ctx.fillText("0.5", 222, 20);
            ctx.beginPath();
            ctx.dashedLine(27, 10, 27, 0, 1);
            ctx.dashedLine(417, 10, 417, 0, 1);
            ctx.dashedLine(222, 10, 222, 0, 1);
            ctx.stroke();
        }
    }

    if ($('.legendCanvas').length)
        renderLegendCanvas();

    // Render all canvas
    //NOTE! Canvas with bins are rendered in the customized widgets section of jquery.tablesorter.js

    /**
     * Delta PSI rendering
     * */
    var deltaTables = $('.deltaTable');
    if (deltaTables.length){
        deltaTables.css({"width": "auto"});
    }


    /**
     * Adding GC Content plot
     */
    $('.gcContent').on("click", function(e){
        e.preventDefault();
        my_window = window.open("", "GC Content", "status=1,width=1024,height=800");
        my_window.document.write('<h1>GC content</h1>');
        my_window.document.write('<img src="' + $(this)[0].getAttribute('data-gc-src') + '" alt="GC Content" border="0">');
    });


    // Single LSVs - TODO: Move from here to jquery.tablesorter.js
    $('.lsvLegendThumb').each( function(){
        var collapsed = this.getAttribute('data-collapsed');
        if (collapsed){
            splicegraph().renderLsvLegend(this);
        } else {
            splicegraph().renderLsvSpliceGraph(this);
        }
        var can = this;
        function dlCanvas() {
            var dt = can.toDataURL('image/png');
            this.href = dt;
        };
        var dl_canvas_link = $(this).parent().children(".lsv_type")[0];
        dl_canvas_link.addEventListener('click', dlCanvas, false);


    });

    var tooltips = $('.tooltip');
    if (tooltips.length){
        $('.tooltip').tooltipster({
            theme: 'tooltipster-light'
        });
    }

});

var initLargeCanvasSettings = function (num_bins, canvas) {

    // Calculate canvas drawable settings
    var settingsCanvas;

    return function(){
        if (settingsCanvas){
            return settingsCanvas;
        }
        settingsCanvas = {};
        settingsCanvas.margin_y = [10, 25];
        settingsCanvas.margin_x = [30, 1];
        settingsCanvas.area_pixels = [canvas.width - settingsCanvas.margin_x.reduce(add, 0),
                canvas.height - settingsCanvas.margin_y.reduce(add, 0)];
        settingsCanvas.coords_origin = [settingsCanvas.margin_x[0], settingsCanvas.margin_y[0] + settingsCanvas.area_pixels[1]];
        settingsCanvas.bar_width = settingsCanvas.area_pixels[0]/num_bins;
        settingsCanvas.labels_steps = [Math.ceil((num_bins+1)/8), 4];
        return settingsCanvas;
    }();

};

function addExpandedViewFunction(){
    // detailed view buttons behaviour
    $('.detailedPsiLink').on("click", function(e){
        e.preventDefault();
        var iconExpand = $(this).get(0);

        iconExpand.setAttribute('isExpanded', (parseInt(iconExpand.getAttribute('isExpanded')) + 1) % 2) ;
        var expandIcon;
        if (iconExpand.getAttribute('isExpanded') == false){
            expandIcon = "../templates/static/img/arrow_down.png";
            hideResetZoomLink($(this));
        } else {
            expandIcon = "../templates/static/img/arrow_up.png";
        }

        // change icon
        iconExpand.setAttribute('src', expandIcon);
        var idDetailed = parseInt(this.id.match(/(\d+)$/)[0], 10);
        var targetPlot = $('#psiExtendedPlot'+idDetailed);
        targetPlot.toggle("show").focus();
    });

    $('.expandDelta').on("click", function(e){
        e.preventDefault();
        var iconExpand = $(this).get(0);

        iconExpand.setAttribute('data-is-expanded', (parseInt(iconExpand.dataset.isExpanded) + 1) % 2) ;
        var expandIcon;
        if (iconExpand.getAttribute('data-is-expanded') == false){
            expandIcon = "../templates/static/img/arrow_down.png";
            hideResetZoomLink($(this));
        } else {
            expandIcon = "../templates/static/img/arrow_up.png";
        }

        // change icon
        iconExpand.setAttribute('src', expandIcon);
        var idDetailed = parseInt(this.id.match(/(\d+)$/)[0], 10);
        $('#deltaPsi'+idDetailed).toggle("show");
        var targetElement = $('#extendedDeltaPsi'+idDetailed);
        targetElement.toggle("show", function(){
            $(this).get(0).focus();
        });


    });

    // detailed view buttons behaviour
    $('.zoomMatrix').on("click", function(e){
        e.preventDefault();
        var iconExpand = $(this).get(0);

        iconExpand.setAttribute('data-is-expanded', (parseInt(iconExpand.getAttribute('data-is-expanded')) + 1) % 2) ;
        var expandIcon;
        if (iconExpand.getAttribute('data-is-expanded') == false){
            expandIcon = "static/img/zoom_in.png";
        } else {
            expandIcon = "static/img/zoom_out.png";
        }

        // change icon
        iconExpand.setAttribute('src', expandIcon);
        var idDetailed = parseInt(this.id.match(/(\d+)$/)[0], 10);
        $('#matrixPlot'+idDetailed).toggle("show").focus();
    });

    $('.expandLsv').on("click", function(e){
        e.preventDefault();
        var iconExpand = $(this).get(0);

        iconExpand.setAttribute('data-is-expanded', (parseInt(iconExpand.dataset.isExpanded) + 1) % 2) ;
        var expandIcon;
        if (iconExpand.getAttribute('data-is-expanded') == false){
            expandIcon = "static/img/arrow_down.png";
        } else {
            expandIcon = "static/img/arrow_up.png";
        }

        // change icon
        iconExpand.setAttribute('src', expandIcon);
        var lsv_selected = $(this).siblings('.lsvSingleExtendedBoxplot').get(0);

        $(lsv_selected).toggle("show");
        drawLSEBoxplots(lsv_selected);
        var targetElement = $(this).siblings('.lsvSingleCompactPercentiles').get(0);
        $(targetElement).toggle("show", function(){
            $(this).get(0).focus();
        });


    });

}


function drawLine(contextO, startx, starty, endx, endy) {
  contextO.beginPath();
  contextO.moveTo(startx, starty);
  contextO.lineTo(endx, endy);
  contextO.closePath();
  contextO.stroke();
}

function drawDashedLine(contextO, startx, starty, endx, endy, dashLen) {
    contextO.beginPath();
    contextO.dashedLine(startx, starty,endx, endy, dashLen);
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
    if (typeof stroke == "undefined" ) {
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


function createGradientSinglePlots(ctx, canvas, margins){
    //create a gradient object from the canvas context
    var gradient = ctx.createLinearGradient(0, 0, canvas.width, 0);

    // Add the colors with fixed stops at 1/4 of the width.
//    gradient.addColorStop(0, "#0000FF");
//    gradient.addColorStop(.45, "cyan");
//    gradient.addColorStop(.50, "#00FF00");
//    gradient.addColorStop(.65, "yellow");
//    gradient.addColorStop(1, "#FF0000");
    if (margins){
        gradient = ctx.createLinearGradient(margins[0], 0, canvas.width - margins[1], 0);
    }

    gradient.addColorStop(1, rgbToHex(213,62,79));
    gradient.addColorStop(.875, rgbToHex(244,109,67));
    gradient.addColorStop(.75, rgbToHex(253,174,97));
    gradient.addColorStop(.625, rgbToHex(254,224,139));
    gradient.addColorStop(.5, rgbToHex(255,255,191));
    gradient.addColorStop(.375, rgbToHex(230,245,152));
    gradient.addColorStop(.25, rgbToHex(171,221,164));
    gradient.addColorStop(.125, rgbToHex(102,194,165));
    gradient.addColorStop(0, rgbToHex(50,136,189));


    return gradient;
}

function createGradientDeltaPlots(ctx, gradientFrom, gradientTo, fromColor, toColor, pThreshold){
    //create a gradient object from the canvas context
    var gradient = ctx.createLinearGradient(gradientFrom, 0, gradientTo, 0);

    // Add the colors with fixed stops.
    gradient.addColorStop(0, getColor(fromColor, BREWER_PALETTE,.5));
//    gradient.addColorStop(0.5-pThreshold/2,"white");
//    gradient.addColorStop(0.5-pThreshold/2,"grey");
////    gradient.addColorStop(0.5,"white");
//    gradient.addColorStop(0.5+pThreshold/2,"grey");
//    gradient.addColorStop(0.5+pThreshold/2,"white");
    gradient.addColorStop(1, getColor(toColor, BREWER_PALETTE,.5));

    return gradient;
}

function drawBarchartWithCanvasId(eventOrCanvasid, zoomInc, canvasSettings){
    // To cover the case when the reset button is clicked and when we just render the canvas passing a canvas id

    var canvasId;
    if (eventOrCanvasid.type == "click"){
        // This method was fired by a click event! We need to retrieve the params...
        zoomInc = eventOrCanvasid.data.zoomInc;
        canvasId = eventOrCanvasid.data.canvasId;
    } else {
        canvasId = eventOrCanvasid;
    }
    canvas = $("#".concat(canvasId))[0];
    renderLargeCanvas(canvas, zoomInc, canvasSettings);
}

function attachResetZoomLink(canvas, drawFunction, canvasSettingsP) {
    var idResetImage = 'resetZoom'.concat(canvas.getAttribute('id'));
    //find the td in which we need to append the reset button
    $(canvas).parent().prepend('<img id="'.concat(idResetImage).concat('" src="static/img/undo.png" class="hidden resetZoom tooltip" title="Reset zoom" />'));
    $("#".concat(idResetImage)).on("click", {canvasId: canvas.id, zoomInc: 0, canvasSettings: canvasSettingsP}, drawFunction);
}

function drawInitialBarplotCanvas(canvas){
    canvas.zoom = 0; // TODO: refactor, change local attributes with "element".setAttribute('att_name', 'att_value')

    attachResetZoomLink(canvas, drawBarchartWithCanvasId);
    renderLargeCanvas(canvas, 0);
}

function drawBarChart(context, bins, settingsCanvasSingle, zoomLevel) {

    // Draw the x and y axes
    context.lineWidth = "1.0";

    // Y-axis
    drawLine(context,
        settingsCanvasSingle.coords_origin[0],
        settingsCanvasSingle.coords_origin[1]+5, //Mark axis
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

    for (var i=0; i < bins.length; i++) {
        var name = i/bins.length;
        var height = parseInt(bins[i]*chartHeight);
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
        if ( i % settingsCanvasSingle.labels_steps[0] == 0 ){
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
    var convertedMarkDataIncrementsIn = settingsCanvasSingle.area_pixels[1]/settingsCanvasSingle.labels_steps[1]; // chartHeight/settingsCanvasSingle.labels_steps[1];

    context.textAlign = "right";
    context.fillStyle = "black";
    var markerValue = 0;
    var offset = 0;
    while (offset <= settingsCanvasSingle.area_pixels[1]) {
        context.fillText(markerValue.toFixed(1), settingsCanvasSingle.coords_origin[0] - 2, settingsCanvasSingle.coords_origin[1] - offset, 50);
        markerValue += (100/settingsCanvasSingle.labels_steps[1])*(settingsCanvasSingle.area_pixels[1]/chartHeight);
        offset += convertedMarkDataIncrementsIn;
    }
}

function renderLargeCanvas(canvas, zoomInc){
    if (canvas.getContext) {  // check for support
        var ctx = canvas.getContext("2d");

        // Clear canvas
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        // Create gradient for coloured barchart
        ctx.fillStyle = createGradientSinglePlots(ctx, canvas);

        var bins = JSON.parse(canvas.getAttribute('bins'));

        // Compute the zoom
        canvas.zoom = Math.max(canvas.zoom + zoomInc*15, 0); // Only positive zoom
        canvas.zoom = Math.abs(zoomInc)*canvas.zoom; // To reset zoom
        (canvas.zoom == 0) ? hideResetZoomLink(canvas) : showResetZoomLink(canvas);

        drawBarChart(ctx, bins, window.settingsCanvasSingle, canvas.zoom); // Canvas settings shared by all (in window)
    }

}

function showResetZoomLink(canvas){
    $(canvas).parent().children('.resetZoom').removeClass('hidden');
}

function hideResetZoomLink(canvas){
    $(canvas).parent().children('.resetZoom').addClass('hidden');
}

function drawBoxplotHeatmap(canvas) {
    if (canvas.getContext) {  // check for support
        var ctx = canvas.getContext("2d");

        // Clear canvas
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        //create a gradient object from the canvas context
        var gradient = createGradientSinglePlots(ctx, canvas);


        // TODO: Refactor this code: replace canvas.attributes
        var sizeBins = JSON.parse(canvas.attributes.bins.nodeValue).length;
        var mean_x = parseInt(canvas.attributes.mean_psi.nodeValue * canvas.width);

        // TODO: Recfactor this code
        var quartiles = JSON.parse(canvas.attributes.quartiles.nodeValue);
        var quartile_10 = parseInt(quartiles[0]/sizeBins * canvas.width);
        var quartile_25 = parseInt(quartiles[1]/sizeBins * canvas.width);
        var quartile_50 = parseInt(quartiles[2]/sizeBins * canvas.width);
        var quartile_75 = parseInt(quartiles[3]/sizeBins * canvas.width);
        var quartile_90 = parseInt(quartiles[4]/sizeBins * canvas.width);


        var previous_style = ctx.fillStyle;
        // Use the gradient to draw
        ctx.fillStyle = gradient;
        drawRectangle(ctx, quartile_25, 5, quartile_75 - quartile_25, canvas.height-10, 1); // Box
        ctx.fillStyle = previous_style;

//        drawRectangle(ctx, 0, 0, canvas.width, canvas.height, 0); // Canvas box
        drawLine(ctx, quartile_50, 5, quartile_50, canvas.height - 5) ;     // Median
        drawLine(ctx, quartile_10, canvas.height/2, quartile_25-1, canvas.height/2 );// Range left
        drawLine(ctx, quartile_75+1, canvas.height/2, quartile_90, canvas.height/2 ); // Range right

        ctx.strokeStyle = "#FF0000"; //Not working...
        ctx.lineWidth = 2;
        ctx.beginPath();
//            ctx.dashedLine(mean_x, 0, mean_x, canvas.height , 1); // Mean
        drawLine(ctx, mean_x, 0, mean_x, canvas.height); // Mean
        ctx.stroke();
        ctx.lineWidth = 1;
        ctx.strokeStyle = "#000000";

    }
}

// Get Color from Brewer Palette
BREWER_PALETTE = [
    [228,26,28],
    [55,126,184],
    [77,175,74],
    [152,78,163],
    [255,127,0],
//    [255,255,51],
    [166,86,40],
    [247,129,191],
    [153,153,153],

    [28,126,128],
    [155,226,29],
    [177,275,19],
    [252,178,8],
    [55,227,100],
//    [55,55,151],
    [266,186,140],
    [47,229,36],
    [253,253,253]
];


function getColor(colorNumber, palette, hue){
    return "rgba("+palette[colorNumber].toString()+", "+hue+")";
//        return rgbToHex(palette[colorNumber][0], palette[colorNumber][1], palette[colorNumber][2]);
}

function drawLSVGroups(canvas){

    function createGradientLSVGroups(coords, count) {
        //create a gradient object from the canvas context
        var gradient = ctx.createLinearGradient(coords.x1, coords.y1, coords.x2, coords.y2);

        // Add the colors with fixed stops.
        gradient.addColorStop(0,"white");
//        gradient.addColorStop(0,getColor(count));  // This suppress the color gradient
        gradient.addColorStop(1,getColor(count, BREWER_PALETTE, 1));

        return gradient;
    }

    if (canvas.getContext) {  // check for support
        var ctx = canvas.getContext("2d");

        // Clear canvas
        ctx.clearRect(0, 0, canvas.width, canvas.height);
//        var groups_str = '{"Brain1": {"quartiles": [0.0085869884077270017, 0.0185507490528867, 0.039331969099260386, 0.070499931202467433, 0.10776128732936643], "mean": 0.7040854774660984}, "Liver1": {"quartiles": [0.0167127755328524, 0.027718954397289025, 0.045679859838293688, 0.068044715851760751, 0.093933591236696534], "mean": 0.051347056627349189}, "CNS1": {"quartiles": [0.0050204656945765352, 0.0094307862661881618, 0.016861982876914608, 0.02711108795391242, 0.039047249315262544], "mean": 0.019890393118988429}, "Liver2": {"quartiles": [0.033018523374089559, 0.052006124437603635, 0.08022851891171659, 0.11349486090515559, 0.15063461121609234], "mean": 0.087625252655997701}, "Brain2": {"quartiles": [0.16696786454087925, 0.23725294831529714, 0.31983500541741855, 0.42557043358099123, 0.52059433186890247], "mean": 0.33619105810705663}}';
        var groups_str = '{"Brain1": {"PSI1": {"quartiles": [0.58558305885339956, 0.64927950284874325, 0.7123383481909511, 0.76410640951310627, 0.81020384655356392], "mean": 0.7040854774660984}, "PSI2": {"quartiles": [0.1399576942378119, 0.18648402982087461, 0.23816476004292436, 0.29543137579707301, 0.36034427584126916], "mean": 0.2456880468921798}, "PSI3": {"quartiles": [0.0085869884077270017, 0.0185507490528867, 0.039331969099260386, 0.070499931202467433, 0.10776128732936643], "mean": 0.05022647564172205}}, "Liver1": {"PSI1": {"quartiles": [0.017312674221672222, 0.027306003129876962, 0.04387031495479439, 0.066852711784719368, 0.09095349084757412], "mean": 0.050121886445342538}, "PSI2": {"quartiles": [0.84268654065686244, 0.87360616269822611, 0.90371812976703758, 0.92839330731234337, 0.94836009037537949], "mean": 0.89853105692730806}, "PSI3": {"quartiles": [0.0167127755328524, 0.027718954397289025, 0.045679859838293688, 0.068044715851760751, 0.093933591236696534], "mean": 0.051347056627349189}}, "CNS1": {"PSI1": {"quartiles": [0.72634413748218041, 0.75492912648127097, 0.78489698300156308, 0.81318882833940087, 0.83444882266707976], "mean": 0.78196816665515079}, "PSI2": {"quartiles": [0.14808838181253248, 0.16819905374556013, 0.19378444116212087, 0.22541219253649442, 0.25354826411316617], "mean": 0.19814144022586128}, "PSI3": {"quartiles": [0.0050204656945765352, 0.0094307862661881618, 0.016861982876914608, 0.02711108795391242, 0.039047249315262544], "mean": 0.019890393118988429}}, "Liver2": {"PSI1": {"quartiles": [0.035286197919954999, 0.05253611937149797, 0.081574971529918078, 0.11692023224409571, 0.15648746333728336], "mean": 0.089254143288114138}, "PSI2": {"quartiles": [0.73424632923438948, 0.78441260247880984, 0.83017998625195932, 0.87020836193544171, 0.90181310236123868], "mean": 0.82312060405588838}, "PSI3": {"quartiles": [0.033018523374089559, 0.052006124437603635, 0.08022851891171659, 0.11349486090515559, 0.15063461121609234], "mean": 0.087625252655997701}}, "Brain2": {"PSI1": {"quartiles": [0.10474716073787305, 0.15545030362768197, 0.2268214264629328, 0.31877883386939282, 0.40346178787707565], "mean": 0.24247623462756193}, "PSI2": {"quartiles": [0.2307375733025599, 0.31663796698787738, 0.41937502332727283, 0.52272833287239551, 0.61401689633784939], "mean": 0.42133270726538152}, "PSI3": {"quartiles": [0.16696786454087925, 0.23725294831529714, 0.31983500541741855, 0.42557043358099123, 0.52059433186890247], "mean": 0.33619105810705663}}}'

//        var groups = JSON.parse(canvas.getAttribute('data-groups'));  // Not used at this stage
        var groups = JSON.parse(groups_str.replace(/"""/g, "'"));


        // Calculate origins_coords
        var header_height = canvas.height*.1;
        var num_groups = Object.keys(groups).length;
        for (var kk=0;kk<num_groups; kk++){
            var num_psis = Object.keys(groups[kk]).length;
            break;
        }
        var sub_canvas_w = canvas.width / num_groups;
        var sub_canvas_h = canvas.height-header_height;
        var sub_canvas_margins = [sub_canvas_w *.1, sub_canvas_w *.05, header_height+sub_canvas_h*.05, sub_canvas_h *.1];
        var sub_canvas_pixels = [sub_canvas_w-sub_canvas_margins[0] -sub_canvas_margins[1], canvas.height-sub_canvas_margins[3]-sub_canvas_margins[2]];

        // Draw sub-boxes separators and headers
        drawLine(ctx, 0, header_height, canvas.width, header_height);
        ctx.textAlign = "center";
        ctx.font = "10pt Arial";

        var i=0;
        for (var key in groups){
            drawLine(ctx, sub_canvas_w*i, 0, sub_canvas_w*i, canvas.height);
            i++;
            ctx.fillText(key,sub_canvas_w*(i-1/2), header_height/2+5);
        }

        var origins_coords = [];
        for (var count_groups=0; count_groups<num_groups; count_groups++){
            origins_coords[count_groups]=[sub_canvas_margins[0] + count_groups*sub_canvas_w, canvas.height-sub_canvas_margins[3]];
        }

        // Body of the plotting function
        var ii=0;


        for (var group_key in groups){
            // X-axis
            drawLine(ctx, origins_coords[ii][0]-5, origins_coords[ii][1], origins_coords[ii][0]+sub_canvas_w-(sub_canvas_margins[0] + sub_canvas_margins[1])+1, origins_coords[ii][1]);
            // Y-axis
            drawLine(ctx, origins_coords[ii][0], origins_coords[ii][1]+5, origins_coords[ii][0], sub_canvas_margins[2]);
            var prev_style = ctx.fillStyle;
            var prev_font = ctx.font;
            ctx.fillStyle = "#000000";
            ctx.font = "8pt Arial";
            ctx.fillText("0", origins_coords[ii][0]-6, origins_coords[ii][1]+10);
            ctx.fillText("1", origins_coords[ii][0]-6, origins_coords[ii][1]-sub_canvas_pixels[1]+5);
//            ctx.fillStyle = prev_style;
            ctx.font = prev_font;

            // Separators
            var offset = 0;
            var count = 0;
            for (var psi_name in groups[group_key]){
                var coords_gradient = {
                    'x1': origins_coords[ii][0]+offset, 'y1': origins_coords[ii][1],
                    'x2': origins_coords[ii][0]+offset+sub_canvas_pixels[0]/num_psis, 'y2': origins_coords[ii][1]-sub_canvas_pixels[1]};
                ctx.fillStyle = createGradientLSVGroups(coords_gradient, count);

                count++;
                var group = groups[group_key][psi_name];

                // 25 and 75 percentile box
                drawRectangle(ctx, origins_coords[ii][0]+offset,
                    Math.round(origins_coords[ii][1]-group.quartiles[1]*sub_canvas_pixels[1]),
                    sub_canvas_pixels[0]/num_psis,
                    -Math.round((group.quartiles[3]-group.quartiles[1])*sub_canvas_pixels[1]), true
                );
                ctx.strokeRect(origins_coords[ii][0]+offset,
                    Math.round(origins_coords[ii][1]-group.quartiles[1]*sub_canvas_pixels[1]),
                    sub_canvas_pixels[0]/num_psis,
                    -Math.round((group.quartiles[3]-group.quartiles[1])*sub_canvas_pixels[1]), true
                );
                // 10 and 75 tales
                drawLine(ctx,
                    Math.round(origins_coords[ii][0]+offset+sub_canvas_pixels[0]/num_psis/2),
                    Math.round(origins_coords[ii][1]-group.quartiles[0]*sub_canvas_pixels[1]),
                    Math.round(origins_coords[ii][0]+offset+sub_canvas_pixels[0]/num_psis/2),
                    Math.round(origins_coords[ii][1]-group.quartiles[1]*sub_canvas_pixels[1])
                );
                drawLine(ctx,
                    Math.round(origins_coords[ii][0]+offset+sub_canvas_pixels[0]/num_psis/2),
                    Math.round(origins_coords[ii][1]-group.quartiles[3]*sub_canvas_pixels[1]),
                    Math.round(origins_coords[ii][0]+offset+sub_canvas_pixels[0]/num_psis/2),
                    Math.round(origins_coords[ii][1]-group.quartiles[4]*sub_canvas_pixels[1])
                );

                // Median (50 percentile)
                drawLine(ctx, origins_coords[ii][0]+offset, Math.round(origins_coords[ii][1]-group.quartiles[2]*sub_canvas_pixels[1]),
                    origins_coords[ii][0]+offset+sub_canvas_pixels[0]/num_psis,
                    Math.round(origins_coords[ii][1]-group.quartiles[2]*sub_canvas_pixels[1]));

                // Mean
                ctx.lineWidth = 2;
                ctx.strokeStyle = rgbToHex(197,27,138);
                drawLine(ctx, origins_coords[ii][0]+offset, Math.round(origins_coords[ii][1]-group.mean*sub_canvas_pixels[1]), origins_coords[ii][0]+offset+sub_canvas_pixels[0]/num_psis, Math.round(origins_coords[ii][1]-group.mean*sub_canvas_pixels[1]));
                ctx.strokeStyle = "#000000";
                ctx.lineWidth = 1;

                offset += sub_canvas_pixels[0]/num_psis;
                if (count < num_psis){
                    ctx.lineWidth = .5;
                    drawDashedLine(ctx, origins_coords[ii][0]+offset, sub_canvas_margins[2]-sub_canvas_margins[2] *.1, origins_coords[ii][0]+offset, canvas.height-sub_canvas_margins[3] *.2, 1);
                }
                ctx.lineWidth = 1;
                ctx.fillStyle = '#000000';
                ctx.font = '8pt Arial';
                ctx.fillText(psi_name, origins_coords[ii][0]+offset-sub_canvas_pixels[0]/num_psis/2, canvas.height-sub_canvas_margins[3] *.5);

            }

            ii++;
        }
    }
}

function drawLSVCompactStackBars(canvas, fillMode){

    function getColor(colorNumber, palette, hue){
        return "rgba("+palette[colorNumber].toString()+", "+hue+")";
//        return rgbToHex(palette[colorNumber][0], palette[colorNumber][1], palette[colorNumber][2]);
    }

    function createGradientLSVGroupsCompact(coords, group, count, fillMode, hue) {
        //create a gradient object from the canvas context

        var gradient = ctx.createLinearGradient(coords.x1, coords.y1, coords.x2, coords.y2);
        // Add the colors with fixed stops.
        if (fillMode < 2 || fillMode === 4){
            gradient.addColorStop(0,getColor(count, BREWER_PALETTE, 1));  // No gradient
            gradient.addColorStop(1,getColor(count, BREWER_PALETTE, 1));

        } else if (fillMode >= 2){
            if (hue)
            {
                gradient.addColorStop(0,getColor(count, BREWER_PALETTE, hue));  // No gradient
                gradient.addColorStop(1,getColor(count, BREWER_PALETTE, hue));
            } else {
                var newHue = 1;
                if (fillMode === 2){
                    newHue = 1- (group.quartiles[count][4]-group.quartiles[count][0]);
                } else {
                    newHue = 1-group.variances[count];
                }
                gradient.addColorStop(0,getColor(count, BREWER_PALETTE, newHue));  // No gradient
                gradient.addColorStop(1,getColor(count, BREWER_PALETTE, newHue));
            }
        }
        return gradient;
    }

    var fillingPSIArea = function (ctx, fillMode, group, lsv_count, x1, y1, x2, y2) {
        var area = [];
        if (fillMode === 0){
            // Fill from center to both left and right
            // First approach, use 1 - the 25 to 75 percentile to fill the area
            area[0] = Math.round((1 - (group.quartiles[lsv_count][4] - group.quartiles[lsv_count][0])) * (x2 - x1));
            area[1] = y2 - y1;
            ctx.strokeStyle = ctx.fillStyle;
            drawRectangle(ctx, x1+(x2-x1-area[0])/2, y1, area[0], area[1], true);
//            ctx.strokeRect(x1+(x2-x1-area[0])/2, y1, area[0], area[1]);
        } else if (fillMode === 1){
            // Fill from center to both left and right
            // Second approach, use 1 - the variance
//            console.log(group.variances[lsv_count]);
//            console.log(d3.variance(translate_lsv_bins(group.bins[lsv_count],1000));

            area[0] = (1 - group.variances[lsv_count]) * (x2 - x1);
            area[1] = y2 - y1;
            ctx.strokeStyle = ctx.fillStyle;
            drawRectangle(ctx, Math.round(x1+(x2-x1-area[0])/2), Math.round(y1), Math.floor(area[0]), Math.floor(area[1]), true);
//            ctx.strokeRect(x1+(x2-x1-area[0])/2, y1, area[0], area[1]);
        } else if (fillMode === 2 || fillMode === 3 ){
            // Fill all area, use the hue to represent variance
            area[0] = x2 - x1;
            area[1] = y2 - y1;
            ctx.strokeStyle = ctx.fillStyle;
            drawRectangle(ctx, x1, y1, area[0], area[1], true);
//            ctx.strokeRect(x1, y1, area[0], area[1]);
        } else if (fillMode === 4){
            // Fill from left to right
            // Using 1 - the 25 to 75 percentile to fill the area
            area[0] = Math.round((1 - (group.quartiles[lsv_count][4] - group.quartiles[lsv_count][0])) * (x2 - x1));
            area[1] = y2 - y1;
            ctx.strokeStyle = ctx.fillStyle;
            drawRectangle(ctx, x1, y1, area[0], area[1], true);
        }
    };

    if (canvas.getContext) {  // check for support
        var ctx = canvas.getContext("2d");
        // Clear canvas
        ctx.clearRect(0, 0, canvas.width, canvas.height);
//        var groups_str = '{"Brain1": {"quartiles": [0.0085869884077270017, 0.0185507490528867, 0.039331969099260386, 0.070499931202467433, 0.10776128732936643], "mean": 0.7040854774660984}, "Liver1": {"quartiles": [0.0167127755328524, 0.027718954397289025, 0.045679859838293688, 0.068044715851760751, 0.093933591236696534], "mean": 0.051347056627349189}, "CNS1": {"quartiles": [0.0050204656945765352, 0.0094307862661881618, 0.016861982876914608, 0.02711108795391242, 0.039047249315262544], "mean": 0.019890393118988429}, "Liver2": {"quartiles": [0.033018523374089559, 0.052006124437603635, 0.08022851891171659, 0.11349486090515559, 0.15063461121609234], "mean": 0.087625252655997701}, "Brain2": {"quartiles": [0.16696786454087925, 0.23725294831529714, 0.31983500541741855, 0.42557043358099123, 0.52059433186890247], "mean": 0.33619105810705663}}';
//        var groups_str = '{"Brain1": {"PSI1": {"quartiles": [0.58558305885339956, 0.64927950284874325, 0.7123383481909511, 0.76410640951310627, 0.81020384655356392], "mean": 0.7040854774660984}, "PSI2": {"quartiles": [0.1399576942378119, 0.18648402982087461, 0.23816476004292436, 0.29543137579707301, 0.36034427584126916], "mean": 0.2456880468921798}, "PSI3": {"quartiles": [0.0085869884077270017, 0.0185507490528867, 0.039331969099260386, 0.070499931202467433, 0.10776128732936643], "mean": 0.05022647564172205}}, "Liver1": {"PSI1": {"quartiles": [0.017312674221672222, 0.027306003129876962, 0.04387031495479439, 0.066852711784719368, 0.09095349084757412], "mean": 0.050121886445342538}, "PSI2": {"quartiles": [0.84268654065686244, 0.87360616269822611, 0.90371812976703758, 0.92839330731234337, 0.94836009037537949], "mean": 0.89853105692730806}, "PSI3": {"quartiles": [0.0167127755328524, 0.027718954397289025, 0.045679859838293688, 0.068044715851760751, 0.093933591236696534], "mean": 0.051347056627349189}}, "CNS1": {"PSI1": {"quartiles": [0.72634413748218041, 0.75492912648127097, 0.78489698300156308, 0.81318882833940087, 0.83444882266707976], "mean": 0.78196816665515079}, "PSI2": {"quartiles": [0.14808838181253248, 0.16819905374556013, 0.19378444116212087, 0.22541219253649442, 0.25354826411316617], "mean": 0.19814144022586128}, "PSI3": {"quartiles": [0.0050204656945765352, 0.0094307862661881618, 0.016861982876914608, 0.02711108795391242, 0.039047249315262544], "mean": 0.019890393118988429}}, "Liver2": {"PSI1": {"quartiles": [0.035286197919954999, 0.05253611937149797, 0.081574971529918078, 0.11692023224409571, 0.15648746333728336], "mean": 0.089254143288114138}, "PSI2": {"quartiles": [0.73424632923438948, 0.78441260247880984, 0.83017998625195932, 0.87020836193544171, 0.90181310236123868], "mean": 0.82312060405588838}, "PSI3": {"quartiles": [0.033018523374089559, 0.052006124437603635, 0.08022851891171659, 0.11349486090515559, 0.15063461121609234], "mean": 0.087625252655997701}}, "Brain2": {"PSI1": {"quartiles": [0.10474716073787305, 0.15545030362768197, 0.2268214264629328, 0.31877883386939282, 0.40346178787707565], "mean": 0.24247623462756193}, "PSI2": {"quartiles": [0.2307375733025599, 0.31663796698787738, 0.41937502332727283, 0.52272833287239551, 0.61401689633784939], "mean": 0.42133270726538152}, "PSI3": {"quartiles": [0.16696786454087925, 0.23725294831529714, 0.31983500541741855, 0.42557043358099123, 0.52059433186890247], "mean": 0.33619105810705663}}}'

        var groups_str = canvas.getAttribute("data-lsv");
//        var groups_str = '{"Brain1": {"PSI1": {"var": 0.0074555154779834577, "quartiles": [0.58558305885339956, 0.64927950284874325, 0.7123383481909511, 0.76410640951310627, 0.81020384655356392], "mean": 0.7040854774660984}, "PSI2": {"var": 0.0068315643838000482, "quartiles": [0.1399576942378119, 0.18648402982087461, 0.23816476004292436, 0.29543137579707301, 0.36034427584126916], "mean": 0.2456880468921798}, "PSI3": {"var": 0.0016758508559471695, "quartiles": [0.0085869884077270017, 0.0185507490528867, 0.039331969099260386, 0.070499931202467433, 0.10776128732936643], "mean": 0.05022647564172205}}, "Liver1": {"PSI1": {"var": 0.0008928339892878025, "quartiles": [0.017312674221672222, 0.027306003129876962, 0.04387031495479439, 0.066852711784719368, 0.09095349084757412], "mean": 0.050121886445342538}, "PSI2": {"var": 0.0017652980221384569, "quartiles": [0.84268654065686244, 0.87360616269822611, 0.90371812976703758, 0.92839330731234337, 0.94836009037537949], "mean": 0.89853105692730806}, "PSI3": {"var": 0.00099690656783313412, "quartiles": [0.0167127755328524, 0.027718954397289025, 0.045679859838293688, 0.068044715851760751, 0.093933591236696534], "mean": 0.051347056627349189}}, "CNS1": {"PSI1": {"var": 0.001790118326605698, "quartiles": [0.72634413748218041, 0.75492912648127097, 0.78489698300156308, 0.81318882833940087, 0.83444882266707976], "mean": 0.78196816665515079}, "PSI2": {"var": 0.0017024220254318967, "quartiles": [0.14808838181253248, 0.16819905374556013, 0.19378444116212087, 0.22541219253649442, 0.25354826411316617], "mean": 0.19814144022586128}, "PSI3": {"var": 0.00019317205946046349, "quartiles": [0.0050204656945765352, 0.0094307862661881618, 0.016861982876914608, 0.02711108795391242, 0.039047249315262544], "mean": 0.019890393118988429}}, "Liver2": {"PSI1": {"var": 0.0023980329456464657, "quartiles": [0.035286197919954999, 0.05253611937149797, 0.081574971529918078, 0.11692023224409571, 0.15648746333728336], "mean": 0.089254143288114138}, "PSI2": {"var": 0.0043719250933379064, "quartiles": [0.73424632923438948, 0.78441260247880984, 0.83017998625195932, 0.87020836193544171, 0.90181310236123868], "mean": 0.82312060405588838}, "PSI3": {"var": 0.0023261733084575485, "quartiles": [0.033018523374089559, 0.052006124437603635, 0.08022851891171659, 0.11349486090515559, 0.15063461121609234], "mean": 0.087625252655997701}}, "Brain2": {"PSI1": {"var": 0.013894101212886185, "quartiles": [0.10474716073787305, 0.15545030362768197, 0.2268214264629328, 0.31877883386939282, 0.40346178787707565], "mean": 0.24247623462756193}, "PSI2": {"var": 0.020538079779782795, "quartiles": [0.2307375733025599, 0.31663796698787738, 0.41937502332727283, 0.52272833287239551, 0.61401689633784939], "mean": 0.42133270726538152}, "PSI3": {"var": 0.018431257833960251, "quartiles": [0.16696786454087925, 0.23725294831529714, 0.31983500541741855, 0.42557043358099123, 0.52059433186890247], "mean": 0.33619105810705663}}}'

        var groups = JSON.parse(groups_str.replace(/\\\"/g, "\'").replace(/\"/g,"").replace(/'/g, "\""));
//        "[{"bins": [[0.08094728896690942, 0.05535597311102183, 0.04829619922500768, 0.04446645580902593, 0.041721634801770185, 0.039453108348147456, 0.037437526254519314, 0.03558157768952228, 0.03384394113230546, 0.03220600926408134, 0.030658918210746445, 0.029197284429548734, 0.027816385809332732, 0.02651128461279872, 0.02527693196889278, 0.024108598311683816, 0.023002249012738767, 0.021954715867445938, 0.020963668171994174, 0.02002745738858422, 0.019144917971123068, 0.018315183228509967, 0.017537544028374458, 0.01681135483161965, 0.016135981597485458, 0.015510788363935867, 0.014935169318826096, 0.014408646215852305, 0.01393106451248005, 0.013502937175451535, 0.01312601059571582, 0.01280418133602641, 0.012545018321427463, 0.012362448642123864, 0.012281944003008824, 0.012351762318631469, 0.012671178471174774, 0.013477627771113658, 0.015525474200760477, 0.023793558710303378], [0.023793558710303135, 0.015525474200760453, 0.013477627771113655, 0.012671178471174774, 0.012351762318631474, 0.012281944003008822, 0.012362448642123864, 0.012545018321427464, 0.012804181336026415, 0.01312601059571582, 0.013502937175451535, 0.01393106451248005, 0.014408646215852305, 0.014935169318826096, 0.015510788363935865, 0.016135981597485458, 0.01681135483161965, 0.017537544028374458, 0.018315183228509974, 0.01914491797112307, 0.02002745738858422, 0.020963668171994174, 0.021954715867445976, 0.023002249012738767, 0.024108598311683816, 0.025276931968892743, 0.026511284612798725, 0.027816385809332732, 0.02919728442954874, 0.030658918210746452, 0.03220600926408134, 0.0338439411323057, 0.035581577689522303, 0.0374375262545194, 0.03945310834814751, 0.04172163480177015, 0.044466455809025945, 0.04829619922500767, 0.055355973111021854, 0.08094728896690963]], "conf_interval": [[0, 38], [1, 39]], "coords": "chr8:1000058-10000.12117003358", "id": "chr8:1000058-1000358", "matrix_link": "#", "means": [0.36253612941783353, 0.6374638705821889], "name": "PS1", "quartiles": [[1, 4, 11, 22, 33], [6, 17, 28, 35, 38]], "type": "Exon skipping"}]"

        // Calculate origins_coords
        var header_height = 0; // canvas.height*.1;
        var num_groups = groups.length;

        var sub_canvas_w = canvas.width / num_groups;
        var sub_canvas_h = canvas.height-header_height;
//        var sub_canvas_margins = [sub_canvas_w *.00, 0, header_height, sub_canvas_h *.21];
        var sub_canvas_margins = [0, 0, 0, 0];
        var sub_canvas_pixels = [sub_canvas_w-sub_canvas_margins[0] -sub_canvas_margins[1], canvas.height-sub_canvas_margins[3]-sub_canvas_margins[2]];

        // Draw sub-boxes separators and headers
//        drawLine(ctx, 0, sub_canvas_h-sub_canvas_margins[3], canvas.width, sub_canvas_h-sub_canvas_margins[3]);
        ctx.textAlign = "center";
        ctx.font = "8pt Arial";

        var tmpStyle = ctx.strokeStyle;
        ctx.strokeStyle = "rgba(0, 0, 0, .5)";
        var i=0;
        for (var ii=1; ii<num_groups; ii++){
            if (i){
                drawLine(ctx, sub_canvas_w*i, 0, sub_canvas_w*i, canvas.height);
            }
            i++;
            ctx.fillText(groups[ii].name, sub_canvas_w*(i-1/2), sub_canvas_h-sub_canvas_margins[3]+sub_canvas_h*.2 - 1);
        }
        ctx.strokeStyle = tmpStyle;

        var origins_coords = [];
        for (var count_groups=0; count_groups<num_groups; count_groups++){
            origins_coords[count_groups]=[sub_canvas_margins[0] + count_groups*sub_canvas_w, canvas.height-sub_canvas_margins[3]];
        }


        // Separators
        for (var count=0; count<num_groups; count++){
            var offset = 0;
            var acc_height = 0;
            var group = groups[count];

            for (var lsv_count=0; lsv_count<group.bins.length; lsv_count++) {
                // Calculate the height of the accumulated mean
                acc_height += group.means[lsv_count];

                var coords_gradient = {
                    'x1': origins_coords[count][0] + offset, 'y1': origins_coords[count][1],
                    'x2': origins_coords[count][0] + offset + sub_canvas_pixels[0], 'y2': origins_coords[count][1] - sub_canvas_pixels[1]};
                ctx.fillStyle   = createGradientLSVGroupsCompact(coords_gradient, group, lsv_count, fillMode, 1);
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

function drawDeltaLSVCompactSVG(htmlElementId, lsv) {
    /* TODO: Create SVG component; Draw canvas rectangle; Draw for each of the excl-incl pairs, the bars.*/
    var width = 200,
        height = 20;
    var margin = {top: 1, bottom: 8, left: 2, right: 2};
    var domain = [-1, 1];
    var x = d3.scale.linear()
        .range([margin.left, margin.right])
        .domain(domain);
    var border_frame = 2;

    var svgContainer = d3.select("#" + htmlElementId)
        .append("svg")
        .attr("class", "excl-incl-rect")
        .attr("width", width)
        .attr("height", height);

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

        if (last_excl_pos != width / 2 && (Math.round(width / 2 - margin.left) * lsv.excl_incl[ii][0]) >= 1) {
            console.log("This guys has more than 2 ways with exclusion >= 1: " + htmlElementId);
            console.log(last_excl_pos + ", " + (width / 2 - margin.left) * lsv.excl_incl[ii][0])
        }
        if (last_incl_pos != width / 2 && Math.round((width / 2 - margin.right) * lsv.excl_incl[ii][1]) >= 1) {
            console.log("This guys has more than 2 ways with inclusion >= 1: " + htmlElementId);
            console.log(last_incl_pos + ", " + (width / 2 - margin.right) * lsv.excl_incl[ii][1])
        }

        // Draw percentages text
        if (Math.round((width / 2 - margin.left) * lsv.excl_incl[ii][0]) >= 1){
            last_excl_pos -= Math.round((width / 2 - margin.left) * lsv.excl_incl[ii][0]);
            svgContainer.append("text")
                .attr("x", last_excl_pos)
                .attr("y", height)
                .attr("text-anchor", "middle")
                .attr("font-size", "9px")
                .attr("fill", getColor(ii, BREWER_PALETTE, 1))
                .text("-" + Math.round((width / 2 - margin.left) * lsv.excl_incl[ii][0]));

        }

        if (Math.round((width / 2 - margin.right) * lsv.excl_incl[ii][1]) >= 1){
            last_incl_pos += Math.round((width / 2 - margin.right) * lsv.excl_incl[ii][1]);
            svgContainer.append("text")
                .attr("x", last_incl_pos)
                .attr("y", height)
                .attr("text-anchor", "middle")
                .attr("font-size", "9px")
                .attr("fill", getColor(ii, BREWER_PALETTE, 1))
                .text("+" + Math.round((width / 2 - margin.right) * lsv.excl_incl[ii][1]));
        }

    }

    // Draw canvas frame
    svgContainer.append("rect")
        .attr("x", margin.left)
        .attr("y", margin.top)
        .attr("width", width - margin.right - margin.left)
        .attr("height", height - margin.bottom - margin.top)
        .style('stroke', 'black')
        .style('stroke-width', border_frame)
        .attr("stroke-opacity", .5)
        .style('fill', 'none');

    // Draw x-axis ticks
    svgContainer.append("text")
        .attr("x", width / 2)
        .attr("y", height)
        .attr("text-anchor", "middle")
        .attr("font-size", "10px")
        .attr("fill", "black")
        .text("0");
//    svgContainer.append("text")
//        .attr("x", 0)
//        .attr("y", height)
//        .attr("text-anchor", "start")
//        .attr("font-size", "10px")
//        .attr("fill", "black")
//        .text("-1");
//    svgContainer.append("text")
//        .attr("x", width)
//        .attr("y", height)
//        .attr("text-anchor", "end")
//        .attr("font-size", "10px")
//        .attr("fill", "black")
//        .text("+1");

    // Draw separator
    svgContainer.append("line")
        .attr("x1", width / 2)
        .attr("y1", margin.top)
        .attr("x2", width / 2)
        .attr("y2", height - margin.bottom)
        .attr("stroke-width", 2)
        .attr("stroke-opacity", .8)
        .attr("stroke", "black");

    return svgContainer;

}


function drawLSEBoxplots(canvas) {
    if (canvas.getContext) {  // check for support
        var ctx = canvas.getContext("2d");

        // Clear canvas
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        var lsv = JSON.parse(canvas.getAttribute("data-lsv").replace(/\\\"/g, "\'").replace(/\"/g,"").replace(/'/g, "\""));
        var arrayBins = lsv.bins;
        var arrayQuartiles = lsv.quartiles;
        var arrayMeansPsi = lsv.means;

        // Done: draw axis
        var marginsLseGraph = [30, 20, 10, 40]; //E, W, N, S
        var margin_boxes = .7;

        //create a gradient object from the canvas context
        var gradient = createGradientSinglePlots(ctx, canvas, marginsLseGraph); // TODO: Fix gradient too go from

        // Drawable area in pixels
        var area_pixels = [
            canvas.width-(marginsLseGraph[0]+marginsLseGraph[1]),
            canvas.height-(marginsLseGraph[2] + marginsLseGraph[3])
        ];
        // Y-axis
        drawLine(ctx,
            marginsLseGraph[0],
            marginsLseGraph[2]+area_pixels[1]+5, //Mark axis
            marginsLseGraph[0],
            marginsLseGraph[2]
        );
        // X-axis
        drawLine(ctx,
            marginsLseGraph[0]-5, //Mark axis
            marginsLseGraph[2]+area_pixels[1],
            marginsLseGraph[0]+area_pixels[0],
            marginsLseGraph[2]+area_pixels[1]
        );


        // Draw legends
        ctx.textAlign = "center";
        ctx.fillText("0", marginsLseGraph[0], marginsLseGraph[2]+area_pixels[1]+15);
        drawLine(ctx, marginsLseGraph[0]+area_pixels[0]/2, marginsLseGraph[2]+area_pixels[1], marginsLseGraph[0]+area_pixels[0]/2, marginsLseGraph[2]+area_pixels[1]+5);
        ctx.fillText("0.5", marginsLseGraph[0]+area_pixels[0]/2, marginsLseGraph[2]+area_pixels[1]+15);
        drawLine(ctx, marginsLseGraph[0]+area_pixels[0], marginsLseGraph[2]+area_pixels[1]-1, marginsLseGraph[0]+area_pixels[0], marginsLseGraph[2]+area_pixels[1]+6);
        ctx.fillText("1", marginsLseGraph[0]+area_pixels[0], marginsLseGraph[2]+area_pixels[1]+15);


        // Done: For each LSE, draw Box Plot
        var sizeBins = arrayBins[0].length; // Bins per LSE? Assuming all same size, check the first one..
        var sub_area = [area_pixels[0], Math.round(area_pixels[1]/arrayQuartiles.length)]; // new height to draw

        for (var i=0; i<arrayQuartiles.length; i++){

            var box_height = sub_area[1]*(1-margin_boxes);
            var mean_x = arrayMeansPsi[i]*sub_area[0];

            var quartiles = arrayQuartiles[i];
            var quartile_10 = parseInt(quartiles[0]/sizeBins * sub_area[0]);
            var quartile_25 = parseInt(quartiles[1]/sizeBins * sub_area[0]);
            var quartile_50 = parseInt(quartiles[2]/sizeBins * sub_area[0]);
            var quartile_75 = parseInt(quartiles[3]/sizeBins * sub_area[0]);
            var quartile_90 = parseInt(quartiles[4]/sizeBins * sub_area[0]);

            var previous_style = ctx.fillStyle;
            // Use the gradient to draw
            ctx.fillStyle = gradient;
            drawRectangle(ctx, marginsLseGraph[0] + quartile_25,
                canvas.height - (marginsLseGraph[3] + sub_area[1]*(i+1/2) + box_height/2),  // Is complicated ...
                quartile_75 - quartile_25,
                box_height, 1); // Box
            ctx.fillStyle = previous_style;

            if (i<arrayQuartiles.length-1){
                // Draw delimiter:
                ctx.lineWidth = .5;
                ctx.beginPath();
                ctx.dashedLine(marginsLseGraph[0],
                    canvas.height - (marginsLseGraph[3] + (i+1)*sub_area[1]),
                    marginsLseGraph[0] + sub_area[0],
                    canvas.height - (marginsLseGraph[3] + (i+1)*sub_area[1]),
                    5
                );
                ctx.closePath();
                ctx.stroke();
                ctx.lineWidth = 1;
            }


            // Median
            drawLine(ctx, marginsLseGraph[0] + quartile_50,
                canvas.height - (marginsLseGraph[3] + sub_area[1]*(i+1/2) + box_height/2),
                marginsLseGraph[0] + quartile_50,
                canvas.height - (marginsLseGraph[3] + sub_area[1]*(i+1/2))  + box_height/2
            );

            // Line 10-25 quartiles
            drawLine(ctx, marginsLseGraph[0] + quartile_10,
                canvas.height - (marginsLseGraph[3] + sub_area[1]*(i+1/2)),
                marginsLseGraph[0] + quartile_25-1,
                canvas.height - (marginsLseGraph[3] + sub_area[1]*(i+1/2))
            );

            // Line 75-90 quartiles
            drawLine(ctx, marginsLseGraph[0] + quartile_75+1,
                canvas.height - (marginsLseGraph[3] + sub_area[1]*(i+1/2)),
                marginsLseGraph[0] + quartile_90,
                canvas.height - (marginsLseGraph[3] + sub_area[1]*(i+1/2))
            );

            // Mean
            var prevStrokeStyle = ctx.strokeStyle;
            ctx.strokeStyle = "#006afc"; // "#FF00FF";
            ctx.lineWidth = 2;
            drawLine(ctx, marginsLseGraph[0] + mean_x,
                canvas.height - (marginsLseGraph[3] + sub_area[1]*(i+1/2) + box_height - 10),
                marginsLseGraph[0] + mean_x,
                canvas.height - (marginsLseGraph[3] + sub_area[1]*(i+1/2))  + box_height -10); // Mean
            ctx.lineWidth = 1;
            ctx.strokeStyle = prevStrokeStyle;

            // Name each Local Single Variant
            var marginNameBox = 10;
            ctx.font = "10pt Arial";
            ctx.textAlign = "left";
            ctx.fillText("PSI"+(i+1),
                marginsLseGraph[0] + marginNameBox + 5,
                canvas.height - (marginsLseGraph[3] + sub_area[1]*(i+1)) + marginNameBox + 15
            );

            var nameBoxDim = [("PSI"+(i+1)).length*9 + 5, 20];
            roundRect(ctx,
                marginsLseGraph[0] + marginNameBox,
                canvas.height - (marginsLseGraph[3] + sub_area[1]*(i+1)) + marginNameBox,
                nameBoxDim[0],
                nameBoxDim[1]
            );

        }
    }

}

function drawDeltaBox(canvas) {
    if (canvas.getContext) {  // check for support
        var ctx = canvas.getContext("2d");

        // Clear canvas
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        var excl_incl_set = JSON.parse(canvas.dataset.exclIncl);
        var halfCanvasWidth = parseInt(canvas.width/2);
        var margin_top = 2;
        var margin_bottom = 10;
        // Create gradient exclusion/inclusion

        ctx.fillStyle = "grey";
        ctx.lineWidth = 0.5;
        drawRectangle(ctx, 1, 1, canvas.width-2, canvas.height-margin_bottom, 0); // Canvas box
        // Fill with gradient
//        ctx.fillStyle=createGradientDeltaPlots(ctx, 0, canvas.width);

        // draw exclusion/inclusion Box
        ctx.fillStyle = "blue";
        ctx.fillStyle = rgbToHex(55, 126, 184); // blue
        var set18qual1 = "rgba(228,26,28)";
        var set18qual2 = "rgba(55,126,184)";
        ctx.fillRect(halfCanvasWidth, margin_top, -parseInt(excl_incl_set[0]*canvas.width/2), canvas.height-margin_top-margin_bottom);
        ctx.fillStyle = "red";
        ctx.fillStyle = rgbToHex(228, 26, 28); // red
        ctx.fillRect(halfCanvasWidth, margin_top, parseInt(excl_incl_set[1]*canvas.width/2), canvas.height-margin_top-margin_bottom);



        // draw -1, 0 and 1
        ctx.textAlign = "center";
        ctx.font = "bold 8pt Arial";
        ctx.fillStyle = "blue";
        ctx.fillStyle = rgbToHex(55, 126, 184); // blue
        ctx.fillText("-"+(excl_incl_set[0]*100).toFixed(1)+" ", Math.max(halfCanvasWidth-Math.max(parseInt(excl_incl_set[0]*canvas.width/2), 10), 20), canvas.height);
        ctx.fillStyle = "red";
        ctx.fillStyle = rgbToHex(228, 26, 28); // red
//        ctx.fillText("1", canvas.width - margin_top, canvas.height);
        ctx.fillText("+"+(excl_incl_set[1]*100).toFixed(1), Math.min(halfCanvasWidth+Math.max(parseInt(excl_incl_set[1]*canvas.width/2), 10), halfCanvasWidth*2-20), canvas.height);
        ctx.fillStyle = "black";
//        ctx.fillText("0", halfCanvasWidth, canvas.height);



        // draw middle bar (0)
        ctx.lineWidth = 1.5;
        ctx.strokeStyle = "black";
        drawLine(ctx, halfCanvasWidth, parseInt(margin_top/2), halfCanvasWidth, canvas.height-margin_bottom+margin_top/2);
    }
}

function drawDeltaBarChart(context, binsArray, settingsCanvasDelta, zoomLevel, pThreshold) {

    // Draw the x and y axes
    context.lineWidth = "1.0";



    var chartHeight = settingsCanvasDelta.area_pixels[1];
    chartHeight += zoomLevel;
    var maxValue = 0;

    for (var binsCount=0; binsCount<binsArray.length; binsCount++) {
        var bins = binsArray[binsCount];
        //Preserve the gradient color style
        var gradientColorStyle = createGradientDeltaPlots(context, settingsCanvasDelta.coords_origin[0],
                settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0], binsCount, binsCount, pThreshold);

        context.fillStyle = getColor(binsCount, BREWER_PALETTE,.5); //gradientColorStyle;
        context.strokeStyle = getColor(binsCount, BREWER_PALETTE,1);

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


        // End Vertical bars
    }

        // Add the column title to the x-axis
    context.textAlign = "center";
    context.fillStyle = "black";
    context.strokeStyle = "black";


    // Y-axis
    drawLine(context,
            settingsCanvasDelta.coords_origin[0] + (settingsCanvasDelta.area_pixels[0])/2,
            settingsCanvasDelta.coords_origin[1]+5, //Mark axis
            settingsCanvasDelta.coords_origin[0] + (settingsCanvasDelta.area_pixels[0])/2,
            settingsCanvasDelta.coords_origin[1] - settingsCanvasDelta.area_pixels[1]);
    // X-axis
    drawLine(context,
        settingsCanvasDelta.coords_origin[0], //Mark axis
        settingsCanvasDelta.coords_origin[1],
            settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0],
        settingsCanvasDelta.coords_origin[1]);

    name = name * 100;
    context.fillText("-1", settingsCanvasDelta.coords_origin[0], settingsCanvasDelta.coords_origin[1] + 15);
    drawLine(context, settingsCanvasDelta.coords_origin[0],
        settingsCanvasDelta.coords_origin[1],
        settingsCanvasDelta.coords_origin[0],
            settingsCanvasDelta.coords_origin[1] + 2 // size of the delimiter=4
    );

    context.fillText("0", settingsCanvasDelta.coords_origin[0] + (settingsCanvasDelta.area_pixels[0]) / 2, settingsCanvasDelta.coords_origin[1] + 15);
    drawLine(context, settingsCanvasDelta.coords_origin[0] + (settingsCanvasDelta.area_pixels[0]) / 2,
        settingsCanvasDelta.coords_origin[1],
            settingsCanvasDelta.coords_origin[0] + (settingsCanvasDelta.area_pixels[0]) / 2,
            settingsCanvasDelta.coords_origin[1] + 2 // size of the delimiter=4
    );

    context.textAlign = "right";
    context.fillText("1", settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0], settingsCanvasDelta.coords_origin[1] + 15);
    drawLine(context, settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0],
        settingsCanvasDelta.coords_origin[1], // +1 to don't overlap
            settingsCanvasDelta.coords_origin[0] + settingsCanvasDelta.area_pixels[0],
            settingsCanvasDelta.coords_origin[1] + 2 // size of the delimiter=4
    );

    // Vertical bars for percentages of inclusion/exclusion (by Default, in -20 and +20)
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
    // Add some data markers to the y-axis
    var convertedMarkDataIncrementsIn = settingsCanvasDelta.area_pixels[1]/settingsCanvasDelta.labels_steps[1]; // chartHeight/settingsCanvasDelta.labels_steps[1];

    context.textAlign = "right";
    context.fillStyle = "black";
    var markerValue = 0;
    var offset = 0;
    while (offset <= settingsCanvasDelta.area_pixels[1]) {
        context.fillText(markerValue.toFixed(1), settingsCanvasDelta.coords_origin[0] - 2, settingsCanvasDelta.coords_origin[1] - offset, 50);
        markerValue += (100/settingsCanvasDelta.labels_steps[1])*(settingsCanvasDelta.area_pixels[1]/chartHeight);
        offset += convertedMarkDataIncrementsIn;
    }


}

function drawExpDeltaWithCanvasId(eventOrCanvasid, zoomInc, canvasSettings){
    // To cover the case when the reset button is clicked and when we just render the canvas passing a canvas id

    var canvasId;
    if (eventOrCanvasid.type == "click"){
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

function renderExpandedDeltaCanvas(canvas, zoomInc, canvasSettings){
    if (canvas.getContext) {  // check for support
        var ctx = canvas.getContext("2d");

        // Clear canvas
        ctx.clearRect(0, 0, canvas.width, canvas.height);

        // Create gradient for coloured barchart
//        ctx.fillStyle = createGradientSinglePlots(ctx, canvas);

        var bins = JSON.parse(canvas.getAttribute('data-lsv')).bins;

        // Compute the zoom
        var zoom = parseInt(canvas.getAttribute('data-zoom'));
        zoom = Math.max(zoom + zoomInc*15, 0); // Only positive zoom
        zoom = Math.abs(zoomInc)*zoom; // To reset zoom

//        canvas.zoom = Math.abs(zoomInc)*canvas.zoom; // To reset zoom
        (zoom == 0) ? hideResetZoomLink(canvas) : showResetZoomLink(canvas);

        // Canvas attributes to reset barchart from zoom view
        canvas.setAttribute('data-zoom', zoom.toString());

        var pThreshold = canvas.getAttribute('data-threshold');
        drawDeltaBarChart(ctx, bins, canvasSettings, zoom, pThreshold); // Canvas settings shared by all (in window)
    }

}

function initExpandedDeltaCanvas(canvas, canvasSettings){
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

function mul(a, b){
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

function renderViolin(htmlElementId, results, tableId, params){


    function addViolin(svg, results, height, width, domain, imposeMax, count, id_svg, id_table){


        var data = d3.layout.histogram()
            .bins(42)
            .range(domain)
        (results);

        var y = d3.scale.linear()
            .range([width/2, 0])
            .domain([0, Math.max(imposeMax, d3.max(data, function(d) { return d.y; }))]);

        var x = d3.scale.linear()
            .range([height-margin.bottom, margin.top]) //-margin.left, margin.right])
            .domain(domain)
            .nice();

        var area = d3.svg.area()
            .interpolate(interpolation)
            .x(function(d) {
                if(interpolation=="step-before")
                    return x(d.x+d.dx/2)
                return x(d.x);
            })
            .y0(width/2)
            .y1(function(d) { return y(d.y); });

        var line=d3.svg.line()
            .interpolate(interpolation)
            .x(function(d) {
                if(interpolation=="step-before")
                    return x(d.x+d.dx/2)
                return x(d.x);
            })
            .y(function(d) { return y(d.y); });

        svg.append("linearGradient")
            .attr("id", "violin-gradient"+id_svg+count+id_table)
            .attr("gradientUnits", "userSpaceOnUse")
            .attr("x1", margin.top).attr("y1", 0 )
            .attr("x2", height-margin.bottom).attr("y2", 0)
            .selectAll("stop")
            .data([
                {offset: "0%", color: getColor(count, BREWER_PALETTE, 1)},
                {offset: "100%", color: getColor(count, BREWER_PALETTE, 1)} //"steelblue" "gray"
//                {offset: "50%", color: "gray"},

            ])
            .enter().append("stop")
            .attr("offset", function(d) { return d.offset; })
            .attr("stop-color", function(d) { return d.color; });



        var gPlus=svg.append("g");
        var gMinus=svg.append("g");

        gPlus.append("path")
            .datum(data)
            .attr("class", "area")
            .attr("d", area);

//        gPlus.append("path")
//            .datum(data)
//            .attr("class", "violin")
//            .attr("d", line);


        gMinus.append("path")
            .datum(data)
            .attr("class", "area")
            .attr("d", area);

//        gMinus.append("path")
//            .datum(data)
//            .attr("class", "violin")
//            .attr("d", line);


        gPlus.attr("transform", "rotate(90,0,0)  translate(0,-"+width+")");
        gMinus.attr("transform", "rotate(90,0,0) scale(1,-1)");

        gPlus.attr('fill','url(#violin-gradient' +id_svg+count+id_table +')');
        gMinus.attr('fill','url(#violin-gradient' +id_svg+count+id_table +')');



    }
    function addBoxPlot(svg, results, height, width, domain, boxPlotWidth){
        var y = d3.scale.linear()
            .range([height-margin.bottom, margin.top])
            .domain(domain);

        var x = d3.scale.linear()
            .range([0, width]);

        var left=0.5-boxPlotWidth/2;
        var right=0.5+boxPlotWidth/2;

        var probs=[0.1,0.25,0.5,0.75,0.9];
        for(var i=0; i<probs.length; i++){
            probs[i]=y(d3.quantile(results, probs[i]));
        }


        svg.append("rect")
            .attr("class", "boxplot fill")
            .attr("x", x(left))
            .attr("width", x(right)-x(left))
            .attr("y", probs[3])
            .attr("height", -probs[3]+probs[1]);
//            .attr("y", probs[3])
//            .attr("height", probs[1]);

        svg.append("circle")
            .attr("class", "boxplot mean")
            .attr("cx", x(0.5))
            .attr("cy", y(d3.mean(results)))
            .attr("r", x(boxPlotWidth/5));

        var iS=[0,2,4];
        var iSclass=["","median",""];
        for(var i=0; i<iS.length; i++){
            svg.append("line")
                .attr("class", "boxplot "+iSclass[i])
                .attr("x1", x(left))
                .attr("x2", x(right))
                .attr("y1", probs[iS[i]])
                .attr("y2", probs[iS[i]]);
        }

        iS=[[0,1],[3,4]];
        for(var i=0; i<iS.length; i++){
            svg.append("line")
                .attr("class", "boxplot")
                .attr("x1", x(0.5))
                .attr("x2", x(0.5))
                .attr("y1", probs[iS[i][0]])
                .attr("y2", probs[iS[i][1]]);
        }

        svg.append("rect")
            .attr("class", "boxplot")
            .attr("x", x(left))
            .attr("width", x(right)-x(left))
            .attr("y", probs[3])
            .attr("height", -probs[3]+probs[1]);


    }

    function addExpectedPSI(svg, mean_value, height, boxWidth, count){
        svg.append("text")
            .attr("x", boxWidth/2)
            .attr("y", height)
            .attr("text-anchor", "middle")
            .attr("font-size", "12px")
            .attr("fill", getColor(count, BREWER_PALETTE, 1))
            .text(mean_value.toFixed(3));
    }


    var margin={top:10, bottom:30, left:30, right:10};
    var width =  100 * results.length ; //element_jq.getAttribute('width');
    var height=  200; //element_jq.getAttribute('height');
    var spacing_space = (width - margin.left - margin.right)*.05;
    var boxWidth=Math.round(((width - margin.left - margin.right)-spacing_space)/results.length);
    var boxSpacing=Math.round(spacing_space/results.length);

    var domain=[0,1];
    if (params.delta){
        domain=[-1,1];
    }

    var resolution=42;
    var interpolation='basis'; // 'step-before'; 'basis'

    var y = d3.scale.linear()
        .range([height-margin.bottom, margin.top])
        .domain(domain);

    var yAxis = d3.svg.axis()
        .scale(y)
        .ticks(6)
        .orient("left")
        .tickSize(5,0,5);


    var svg = d3.select("#"+htmlElementId)
        .append("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("class", "violin-boxplot");

    svg.append("line")
        .attr("class", "boxplot")
        .attr("x1", margin.left)
        .attr("x2", width-margin.right)
        .attr("y1", y(0))
        .attr("y2", y(0));

    for(var i=0; i<results.length; i++){
        results[i]=results[i]; //.sort(d3.ascending);
        var g=svg.append("g").attr("transform", "translate("+(i*(boxWidth+boxSpacing)+margin.left)+",0)");
        addViolin(g, results[i], height, boxWidth, domain, 0.25, i, htmlElementId, tableId);
        addBoxPlot(g, results[i], height, boxWidth, domain, .15);
        addExpectedPSI(g, d3.mean(results[i]), height, boxWidth, i);

    }

    svg.append("g")
        .attr('class', 'axis')
        .attr("transform", "translate("+margin.left+",0)")
        .call(yAxis);

    return svg[0];
}

function translate_lsv_bins(lsv_bins, num_samples) {
    var adjusted_bins = [];
    for (var lsv_way=0; lsv_way<lsv_bins.length; lsv_way++){
        var tmp_bins = [];
        var bins_size = lsv_bins[lsv_way].length;
        for (var ii=1; ii< bins_size + 1; ii++) {
            var num_copies = Math.round(num_samples * lsv_bins[lsv_way][ii - 1]);
            for (var bins_i=0; bins_i<num_copies; bins_i++){
                tmp_bins.push((1/bins_size)/2 + ((ii-1) / bins_size));
            }
            console.log((1/bins_size)/2 + ((ii-1) / bins_size));
        }
        adjusted_bins.push(tmp_bins);
    }

    return adjusted_bins
}


function translate_delta_lsv_bins(lsv_bins, num_samples) {
    var adjusted_bins = [];
    for (var lsv_way=0; lsv_way<lsv_bins.length; lsv_way++){
        var tmp_bins = [];
        var bins_size = lsv_bins[lsv_way].length;
        for (var ii=1; ii< bins_size + 1; ii++) {
            var num_copies = Math.round(num_samples * lsv_bins[lsv_way][ii - 1]);
            for (var bins_i=0; bins_i<num_copies; bins_i++){
                tmp_bins.push(-.975 + ii * 2 / bins_size);
            }
        }
        adjusted_bins.push(tmp_bins);
    }
    return adjusted_bins
}