$(document).on('mouseover', '.junction-grp', function () {
    var container = this.parentElement.parentElement.parentElement;
    var datum = d3.select(this).datum();
    d3.select(container)
        .selectAll('.junction-grp')
        .style('opacity', function (d) {
            if (datum.start !== d.start || datum.end !== d.end)
                return 0.2
        })
});

$(document).on('mouseout', '.junction-grp', function () {
    $('.junction-grp').css('opacity', '');
});

$(document).on('mouseenter', '.exon, .junction-grp, .half-exon, .intron-retention', function () {
    var d = d3.select(this).datum();
    var splice_graph_tools = $(this).closest('.gene-container').children('.splice-graph-tools');
    if (d.half_exon) {
        if (d.half_exon === 'start')
            splice_graph_tools.find('.coordinates').text('MISSING' + ' - ' + d.end);
        else
            splice_graph_tools.find('.coordinates').text(d.start + ' - ' + 'MISSING');
        splice_graph_tools.find('.length').text('UNKNOWN')
    } else {
        splice_graph_tools.find('.coordinates').text(d.start + ' - ' + d.end);
        splice_graph_tools.find('.length').text(d.end - d.start)
    }
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
    var gene_container = this.closest('.gene-container');
    var lsv_ids = get_lsv_ids(gene_container);
    var splice_graphs = gene_container.querySelectorAll('.splice-graph');

    gene_container.querySelector('.splice-graph-container').classList.toggle('default-view');

    for (var i = 0; i < splice_graphs.length; i++) {
        var splice_graph = splice_graphs[i];
        sg.update(splice_graph, lsv_ids)
    }
});


var zoom = function (el, value, reset) {
    var gene_container = el.closest('.gene-container');
    var lsv_ids = get_lsv_ids(gene_container);
    var splice_graphs = gene_container.querySelectorAll('.gene-container .splice-graph');
    var sg_container = gene_container.querySelector('.splice-graph-container');
    var zoom = sg_container.getAttribute('data-zoom');

    if (reset) {
        sg_container.setAttribute('data-zoom', value)
    } else {
        sg_container.setAttribute('data-zoom', Math.max(1, parseFloat(zoom) + value).toString());
    }

    for (var i = 0; i < splice_graphs.length; i++) {
        sg.update(splice_graphs[i], lsv_ids)
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


