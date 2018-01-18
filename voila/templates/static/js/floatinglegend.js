var FloatingLegend = function () {

};

FloatingLegend.prototype.renderFloatingLegend = function (canvas) {
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
    this.renderNumReads(ctx, Math.round(x + (area_figures[0] / num_fig - SEP_FIG) / 2), MARGINS[2] + font_height, 32);
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

};

FloatingLegend.prototype.renderNumReads = function (ctx, x, y, num_reads) {
    if (parseInt(num_reads) === 0) return;
    ctx.fillStyle = "rgba(0, 0, 0, .8)";
    ctx.font = "9pt Arial";
    ctx.textAlign = "center";
    ctx.fillText(num_reads, x, y - 2);
};
