$(document).on('mouseover', '.junction-grp', function () {
    $('.junction-grp').css('opacity', 0.1);
    $(this).css('opacity', '');
});

$(document).on('mouseout', '.junction-grp', function () {
    $('.junction-grp').css('opacity', '');
});

$(document).on('mouseenter', '.exon, .junction-grp, .half-exon', function () {
    var d = d3.select(this).datum();
    var splice_graph_tools = $(this).closest('.gene-container').children('.splice-graph-tools');
    var splice_graph = $(this).closest('.splice-graph');
    var experiment = splice_graph.attr('data-experiment');
    var gene_id = splice_graph.attr('data-gene-id');

    db.get(gene_id).then(function (gene) {
        var exon_type;
        try {
            exon_type = gene.exon_types[d.start][d.end][experiment];
        } catch (TypeError) {
            exon_type = -1
        }
        if ([4, 5].includes(exon_type)) {
            if (exon_type === 4)
                splice_graph_tools.find('.coordinates').text('MISSING' + ' - ' + d.end);
            else
                splice_graph_tools.find('.coordinates').text(d.start + ' - ' + 'MISSING');
            splice_graph_tools.find('.length').text('UNKNOWN')
        } else {
            splice_graph_tools.find('.coordinates').text(d.start + ' - ' + d.end);
            splice_graph_tools.find('.length').text(d.end - d.start)
        }
    })
});

$(document).on('change', '.splice-graph-selectors select', function () {
    var $this = $(this);
    var group = $this.attr('data-group');
    var experiment = $this.find(':selected').attr('value');
    $this
        .closest('.gene-container')
        .find('.splice-graph.' + group)
        .attr('data-experiment', experiment)
        .each(function () {
            sg.update(this)
        })
});

$(document).on('click', '.toggle-scale', function () {
    $(this)
        .closest('.gene-container')
        .find('.splice-graph')
        .toggleClass('default-view')
        .each(function () {
            sg.update(this)
        });
});

var zoom = function (el, value, reset) {
    $(el)
        .closest('.gene-container')
        .find('.splice-graph')
        .each(function () {
            var zoom = this.getAttribute('data-zoom');
            if (reset) {
                this.setAttribute('data-zoom', value);
                sg.update(this)
            } else {
                var zoom_value = parseFloat(zoom) + value;
                if (zoom_value > 0) {
                    this.setAttribute('data-zoom', zoom_value);
                    sg.update(this)
                }
            }
        });
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


