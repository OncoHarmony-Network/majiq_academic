$(document).on('mouseover', '.junction-grp', function () {
    $('.junction-grp').css('opacity', 0.1);
    $(this).css('opacity', '');
});

$(document).on('mouseout', '.junction-grp', function () {
    $('.junction-grp').css('opacity', '');
});

$(document).on('mouseenter', '.exon, .junction-grp, .half-exon', function () {
    var d = d3.select(this).datum();
    var gene_container = $(this).closest('.gene-container');

    if ([4, 5].includes(d.exon_type)) {
        if (d.exon_type === 4)
            gene_container.find('.coordinates').text('MISSING' + ' - ' + d.end);
        else
            gene_container.find('.coordinates').text(d.start + ' - ' + 'MISSING');
        gene_container.find('.length').text('UNKNOWN')
    } else {
        gene_container.find('.coordinates').text(d.start + ' - ' + d.end);
        gene_container.find('.length').text(d.end - d.start)
    }
});

$(document).on('change', '.splice-graph-selectors select', function () {
    var gene_container = $(this).closest('.gene-container');
    var sg_div = gene_container.find('.splice-graph');
    var experiment = gene_container.find('.splice-graph-selectors select option:selected').text().trim();
    sg.update(sg_div[0], experiment)
});

$(document).on('click', '.toggle-scale', function () {
    var gene_container = $(this).closest('.gene-container');
    var splice_graph = gene_container.find('.splice-graph');
    var experiment = gene_container.find('.splice-graph-selectors select option:selected').text().trim();
    splice_graph.toggleClass('default-view');
    sg.update(splice_graph[0], experiment)
});

var zoom = function (el, value, reset) {
    var gene_container = $(el).closest('.gene-container');
    var experiment = gene_container.find('.splice-graph-selectors select option:selected').text().trim();
    var sg_div = gene_container.find('.splice-graph');
    var zoom = sg_div.attr('data-zoom');
    if (reset) {
        sg_div.attr('data-zoom', value);
        sg.update(sg_div[0], experiment)
    }
    else {
        var zoom_value = parseFloat(zoom) + value;
        if (zoom_value > 0) {
            sg_div.attr('data-zoom', zoom_value);
            sg.update(sg_div[0], experiment)
        }
    }

};

$(document).on('click', '.zoom-in', function () {
    zoom(this, 1)
});


$(document).on('click', '.zoom-out', function () {
    zoom(this, -1)
});


$(document).on('click', '.zoom-reset', function () {
    zoom(this, 1, true)
});


