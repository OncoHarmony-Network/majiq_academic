$(document).on('mouseover', '.junction-grp', function () {
    $('.junction-grp').css('opacity', 0.2);
    $(this).css('opacity', '');
});

$(document).on('mouseout', '.junction-grp', function () {
    $('.junction-grp').css('opacity', '');
});

$(document).on('mouseenter', '.exon, .junction-grp, .half-exon, .intron-retention', function () {
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
    var group = this.getAttribute('data-group');
    var experiment = this.options[this.selectedIndex].getAttribute('value');
    var gene_container = $(this).closest('.gene-container')[0];
    var splice_graph = gene_container.querySelector('.splice-graph[data-group=' + group + ']');
    var lsv_ids = get_lsv_ids(gene_container);

    splice_graph.setAttribute('data-experiment', experiment);
    sg.update(splice_graph, lsv_ids)
});

var get_lsv_ids = function (gene_container) {
    var highlights = gene_container.querySelectorAll('.highlight-btn:checked');
    return Array.from(highlights).reduce(function (acc, el) {
        var td = el.parentElement.parentElement.parentElement;
        var lsv_id = td.getAttribute('data-lsv-id');
        var weighted = td.querySelector('.weighted-btn:checked');
        acc.push([lsv_id, Boolean(weighted)]);
        return acc
    }, []);
};


$(document).on('click', '.toggle-scale', function () {
    var gene_container = $(this).closest('.gene-container')[0];
    var lsv_ids = get_lsv_ids(gene_container);
    var splice_graphs = gene_container.querySelectorAll('.splice-graph');
    for (var i = 0; i < splice_graphs.length; i++) {
        var splice_graph = splice_graphs[i];
        splice_graph.classList.toggle('default-view');
        sg.update(splice_graph, lsv_ids)
    }
});


var zoom = function (el, value, reset) {
    var gene_container = $(el).closest('.gene-container')[0];
    var lsv_ids = get_lsv_ids(gene_container);
    var splice_graphs = gene_container.querySelectorAll('.gene-container .splice-graph');
    for (var i = 0; i < splice_graphs.length; i++) {
        var splice_graph = splice_graphs[i];
        var zoom = splice_graph.getAttribute('data-zoom');
        if (reset) {
            splice_graph.setAttribute('data-zoom', value);
            sg.update(splice_graph, lsv_ids)
        } else {
            var zoom_value = parseFloat(zoom) + value;
            if (zoom_value > 0) {
                splice_graph.setAttribute('data-zoom', zoom_value);
                sg.update(splice_graph, lsv_ids)
            }
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


