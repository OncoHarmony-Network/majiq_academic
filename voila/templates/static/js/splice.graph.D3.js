/**
 * Created by abarrera on 8/27/14.
 */
/** Utils */

function clone(obj) {
    if(obj == null || typeof(obj) != 'object')
        return obj;

    var temp = obj.constructor(); // changed

    for(var key in obj) {
        if(obj.hasOwnProperty(key)) {
            temp[key] = clone(obj[key]);
        }
    }
    return temp;
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
    [255, 255, 51],
    [166, 86, 40],
    [247, 129, 191],
    [153, 153, 153],

    [28,126,128],
    [155,226,29],
    [177,275,19],
    [252,178,8],
    [55,227,100],
    [55,55,151],
    [266,186,140],
    [47,229,36],
    [253,253,253]
];

function getColor(colorNumber, palette, hue) {
    return "rgba(" + palette[colorNumber].toString() + ", " + hue + ")";
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
    var reshape_exons = false,
        exons_mapped_tmp = [],
        exon_tmp,
        index_exons_to_update = [];  // For exons_obj that need to be updated (shorted)

    // First offset (first exon can have long intron)
    var offset = reshape_intron({'coords': [0, 0]}, exons[0], reshape_exons);
    var acc_offset = offset;

    // Note: to account for overlapping exons_obj where exons_obj within a very large exon have long introns, we should
    // ^^^^^ store the last
    var coords_extra = [];
    for (var k = 0; k < exons[0].coords_extra.length; k++) {
        coords_extra.push(map(function (x) {
            return add(x, -acc_offset);
        }, exons[0].coords_extra[k]));
    }
    exon_tmp = {
        'coords': map(function (x) {
            return add(x, -acc_offset);
        }, resize_exon(exons[0].coords, reshape_exons)),
        'type_exon': exons[0].type_exon,
        'intron_retention': exons[0].intron_retention,
        'lsv_type': exons[0].lsv_type,
        'a3': exons[0].a3,
        'a5': exons[0].a5
    };
    exon_tmp.coords_extra = coords_extra;
//        exons_mapped.push({'coords': map(function(x) { return add(x, -acc_offset);}, resize_exon(exons_obj[0].coords, reshape_exons)), 'type': exons_obj[0].type_exon});

    exons_mapped_tmp[0] = exon_tmp;
    var last_end = exon_tmp.coords[1];

    for (var i = 0; i < exons[0]['a3'].length; i++) {
        junctions[exons[0]['a3'][i]].coords[1] -= offset;
    }

    for (var i = 0; i < exons[0]['a5'].length; i++) {
        junctions[exons[0]['a5'][i]].coords[0] -= offset;
    }

    for (var i = 1; i < exons.length; i++) {
        offset = reshape_intron(exons[i - 1], exons[i], reshape_exons);
        acc_offset += offset;

        // Check if there are exons_obj to make shorter (intron retention)
        coords_extra = [];
        for (var k = 0; k < exons[i].coords_extra.length; k++) {
            coords_extra.push(map(function (x) {
                return add(x, -acc_offset);
            }, exons[i].coords_extra[k]));
        }
        exon_tmp = {
            'coords': map(function (x) {
                return add(x, -acc_offset);
            }, resize_exon(exons[i].coords, reshape_exons)),
            'type_exon': exons[i].type_exon,
            'intron_retention': exons[i].intron_retention,
            'lsv_type': exons[i].lsv_type,
            'a3': exons[i].a3,
            'a5': exons[i].a5
        };
        exon_tmp.coords_extra = coords_extra;

        exons_mapped_tmp[i] = exon_tmp;

        // Check if any previous exon needs to be updated
        if (exons[i].coords[0] - exons[i - 1].coords[1] < 0) {
            index_exons_to_update.push(i - 1);
        }

        if (offset) {
            for (var k = 0; k < index_exons_to_update.length; k++) {
                if (exons[index_exons_to_update[k]].coords[1] > exons[i].coords[0]) {
                    exons_mapped_tmp[index_exons_to_update[k]].coords[1] -= offset;
                }
            }
        }

        for (var j = 0; j < exons[i]['a3'].length; j++) {
            junctions[exons[i]['a3'][j]].coords[1] -= acc_offset;
        }
        for (var j = 0; j < exons[i]['a5'].length; j++) {
            junctions[exons[i]['a5'][j]].coords[0] -= acc_offset;
        }
    }
    return exons_mapped_tmp;

}

var add_keys = function(elems){
    var elems_dkeys = [];
    for (var i=0; i<elems.length; i++){
        elems_dkeys.push({'key': i, 'value': elems[i]})
    }
    return elems_dkeys;
};

/** End utils */


var mostLeftCoord = function(exons){
    return d3.min(exons, function(d) {
        return d.value.coords[0];
    });
};
var mostRightCoord = function(exons){
    return d3.max(exons, function(d) {
        return d.value.coords[1];
    });
};

var longestJunc = function(junctions, scale){
    var longest=-1;
    for (var i=0; i<junctions.length; i++){
        if ((scale(junctions[i].coords[1]) - scale(junctions[i].coords[0])) > longest){
            longest = scale(junctions[i].coords[1]) - scale(junctions[i].coords[0]);
        }
    }
    return longest;
};

var exonKey = function(d){
    return d.key;
};


var renderExons = function(datap, keyf, scaleX){

    var exons = svgCanvas.selectAll("rect.exon").data(datap, keyf);

    exons.enter().append("rect");
    exons.attr("class", "exon")
        .transition()
        .duration(1000)
        .ease("linear")
        .attr("x", function(d){
            return scaleX(d.value.coords[0]);
        })
        .attr("y", h*JUNC_AREA)
        .attr("width", function(d){
            return Math.max(EXON_MIN_W, Math.round(scaleX(d.value.coords[1]) - scaleX(d.value.coords[0])));
        })
        .attr("height", EXON_H);

    return exons;
};


// Add dispersion factors to the junctions
var addDispersion = function(exons, junctions){
    exons.each(function(d){
        if (d.value.type_exon == 0) {d3.select(this).classed("found", true);}
        if (d.value.type_exon == 1) {d3.select(this).classed("novel", true);}
        if (d.value.type_exon == 2) {d3.select(this).classed("missing", true);}
        for (var i=0; i< d.value.a3.length; i++){
            junctions[d.value.a3[i]].dispFrom = junctions[d.value.a3[0]].coords[0];
            if (i == 0){
                junctions[d.value.a3[i]].dispersion = 1;
            }else{
                junctions[d.value.a3[i]].dispersion = d.value.a3.length - i + 1;
            }
        }
    });

};



var renderJunctions = function(data, scaleX){
    var juncs = svgCanvas.selectAll('ellipse').data(data);
    var maxJunc = longestJunc(data, scaleX);
    juncs.enter().append("ellipse");

    juncs.attr("class", "junction")
        .attr("clip-path", "url(#" + "cut-off-junctions)")
        .transition()
        .duration(1000)
        .ease("linear")

        .attr("ry", function(d){
            return Math.round(
                    (scaleX(d.coords[1]) - scaleX(d.coords[0]))/maxJunc * JUNC_AREA * (h - padding[0] -padding[2])
                    - ((scaleX(d.coords[1]) - scaleX(d.coords[0]))/maxJunc * JUNC_AREA * (h - padding[0] -padding[2])/d.dispersion)  * (d.dispersion - 1 ? 1 : 0)
            );
        })
        .attr("rx", function(d){
            return Math.round((scaleX(d.coords[1]) - scaleX(d.coords[0]))/2);
        })
        .attr("cy", h * JUNC_AREA)
        .attr("cx", function(d){
            return Math.round((scaleX(d.coords[1]) + scaleX(d.coords[0]))/2);
        });
    return juncs;
};



var addClassJuncs= function(juncs){
    juncs.each(function(d){
        if (d.type_junction == 0) {d3.select(this).classed("found", true);}
        if (d.type_junction == 1) {d3.select(this).classed("novel", true);}
        if (d.type_junction == 2) {d3.select(this).classed("missing", true);}
    });
};

var renderNumReads = function(junctions, scaleX) {
    var maxJunc = longestJunc(junctions, scaleX);
    var labels = svgCanvas.selectAll("text.readcounts")
        .data(junctions);
    labels.enter().append("text");
    labels.classed("readcounts", true)
        .transition()
        .duration(1000)
        .ease("linear")
        .text(function (d) {
            if (d.num_reads) {
                return d.num_reads;
            }
            return '';
        })
        .attr("x", function (d) {
            return Math.round((scaleX(d.coords[1]) + scaleX(d.coords[0])) / 2);
        })
        .attr("y", function (d) {
            return Math.round(h * JUNC_AREA - 2 -
                    (scaleX(d.coords[1]) - scaleX(d.coords[0])) / maxJunc * JUNC_AREA * (h - padding[0] - padding[2])
                    + ((scaleX(d.coords[1]) - scaleX(d.coords[0])) / maxJunc * JUNC_AREA * (h - padding[0] - padding[2]) / d.dispersion) * (d.dispersion - 1 ? 1 : 0)
            );
        });

};

var renderNumExons = function(exons, scaleX) {
    var labels = svgCanvas.selectAll("text.numexon")
        .data(exons);
    labels.enter().append("text");
    labels.classed("numexon", true)
        .transition()
        .duration(1000)
        .ease("linear")
        .text(function (d, i) {
            return (1+i).toString();
        })
        .attr("x", function (d, i) {
            return Math.round((scaleX(d.value.coords[1]) + scaleX(d.value.coords[0])) / 2);
        })
        .attr("y", h*JUNC_AREA + EXON_H /2);
};

var updateScale = function(datap) {
    return d3.scale.linear()
        .domain([mostLeftCoord(datap), mostRightCoord(datap)])
        .rangeRound([padding[3], w - padding[1]]);

};

d3.select('.toogleScale').on('click', function(){

    if (d3.select(this).classed('scaled')) {
        renderSpliceGraph(orig_objs.exons, orig_objs.junc);
        d3.select(this).classed('scaled', false);
    } else {
        renderSpliceGraph(exons_mapped, junctions_obj);
        d3.select(this).classed('scaled', true);
    }

});

var renderSpliceGraph = function (exonsp, junctionsp) {

    /** Compute the scale used */
    var scaleX = updateScale(exonsp);

    /** Render exons and compute disperison for junctions */
    var exons = renderExons(exonsp, exonKey, scaleX);
    addDispersion(exons, junctionsp);
    renderNumExons(exonsp, scaleX);

    /** Render junctions and read numbers */
    var junctions = renderJunctions(junctionsp, scaleX);
    addClassJuncs(junctions);
    renderNumReads(junctionsp, scaleX);

    /** Add interactivity for exons ... */
    exons.on('mouseover', function(d){  // Highlight exons when hovered
        d3.select(this).classed("hovered", true);
        // Tooltip
        //Get this x/y values, then augment for the tooltip
        var xPosition = parseFloat(d3.select(this).attr("x")) ;
        var yPosition = h + 120;

        //Update the tooltip position and value
        d3.select(".tooltipD3")
            .style("left", xPosition + "px")
            .style("top", yPosition + "px")
            .select("#coordsLabel")
            .text(function(){
                return exons_obj[d.key].coords[0] + "-" + exons_obj[d.key].coords[1];
            });

        //Show the tooltip
        d3.select(".tooltipD3").classed("hidden", false);

    })
        .on('mouseout', function(d){
            d3.select(this).classed("hovered", false);
            d3.select(".tooltipD3").classed("hidden", true);
        });

    /** ..and junctions! */
    junctions.on('mouseover', function(d, i){
            d3.selectAll('.junction').classed("blurred", true);
            d3.selectAll('.readcounts').classed("blurred", true);
            d3.select(d3.selectAll('.readcounts')[0][i]).classed("blurred", false).classed("highlighted", true);
            d3.select(this).classed("blurred", false);
            d3.select(this).classed("hovered", true);

//        .attr("x", function (d) {
//            return Math.round((scaleX(d.coords[1]) + scaleX(d.coords[0])) / 2);
//        })
//            .attr("y", function (d) {
//                return Math.round(h * JUNC_AREA - 2 -
//                        (scaleX(d.coords[1]) - scaleX(d.coords[0])) / maxJunc * JUNC_AREA * (h - padding[0] - padding[2])
//                        + ((scaleX(d.coords[1]) - scaleX(d.coords[0])) / maxJunc * JUNC_AREA * (h - padding[0] - padding[2]) / d.dispersion) * (d.dispersion - 1 ? 1 : 0)
//                );
//            })

            //Update the tooltip position and value
            d3.select(".tooltipD3")
                .style("left", Math.round((scaleX(d.coords[1]) + scaleX(d.coords[0])) / 2) + "px")
                .style("top", 50 + "px")
                .select("#coordsLabel")
                .text(function(){
                    return orig_objs.junc[i].coords[0] + "-" + orig_objs.junc[i].coords[1];
                });
            //Show the tooltip
            d3.select(".tooltipD3").classed("hidden", false);

        })
        .on('mouseout', function(){
            d3.selectAll('.junction').classed('blurred', false);
            d3.selectAll('.readcounts').classed("blurred", false).classed("highlighted", false);
            d3.select(this).classed("hovered", false);

            //Show the tooltip
            d3.select(".tooltipD3").classed("hidden", true);
        });


};



/**
 * D3 - SpliceGraph
 * */

//var genes_obj = JSON.parse('{\'exons\': [{\'a3\': [], \'a5\': [0], \'coords\': [8509651, 8509995], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [0], \'a5\': [1, 2], \'coords\': [8520289, 8520458], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [1], \'a5\': [3], \'coords\': [8526853, 8527465], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [2, 3], \'a5\': [4, 5], \'coords\': [8528381, 8528388], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [4, 5], \'a5\': [6, 7], \'coords\': [8528477, 8528736], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [7], \'a5\': [10, 8, 9], \'coords\': [8530208, 8530256], \'coords_extra\': [], \'intron_retention\': true, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [8, 10], \'a5\': [11], \'coords\': [8530364, 8530399], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [9, 11], \'a5\': [12], \'coords\': [8531119, 8531272], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [12], \'a5\': [13], \'coords\': [8532419, 8532468], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [13, 6], \'a5\': [14, 15], \'coords\': [8533658, 8533718], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [14], \'a5\': [16, 17], \'coords\': [8536210, 8536311], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [16], \'a5\': [18], \'coords\': [8538548, 8538592], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [17, 18], \'a5\': [19, 20], \'coords\': [8539051, 8539128], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [], \'a5\': [21], \'coords\': [8546353, 8546440], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [15, 19, 21], \'a5\': [22], \'coords\': [8547433, 8548095], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [20, 22], \'a5\': [23], \'coords\': [8550487, 8550608], \'coords_extra\': [], \'intron_retention\': true, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [23], \'a5\': [24], \'coords\': [8550746, 8551289], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [24], \'a5\': [25, 26], \'coords\': [8551908, 8552068], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [25], \'a5\': [27], \'coords\': [8552186, 8552950], \'coords_extra\': [], \'intron_retention\': true, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [26, 27], \'a5\': [], \'coords\': [8553575, 8553998], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}], \'junctions\': [{\'coords\': [8509995, 8520289], \'num_reads\': 219, \'type_junction\': 0}, {\'coords\': [8520458, 8527413], \'num_reads\': 163, \'type_junction\': 0}, {\'coords\': [8520458, 8528381], \'num_reads\': 42, \'type_junction\': 0}, {\'coords\': [8527465, 8528381], \'num_reads\': 320, \'type_junction\': 0}, {\'coords\': [8528384, 8528477], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [8528388, 8528477], \'num_reads\': 446, \'type_junction\': 0}, {\'coords\': [8528522, 8533700], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [8528570, 8530208], \'num_reads\': 530, \'type_junction\': 0}, {\'coords\': [8530246, 8530364], \'num_reads\': 229, \'type_junction\': 0}, {\'coords\': [8530246, 8531119], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [8530256, 8530364], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [8530399, 8531119], \'num_reads\': 196, \'type_junction\': 0}, {\'coords\': [8531272, 8532419], \'num_reads\': 203, \'type_junction\': 0}, {\'coords\': [8532468, 8533658], \'num_reads\': 103, \'type_junction\': 0}, {\'coords\': [8533718, 8536210], \'num_reads\': 234, \'type_junction\': 0}, {\'coords\': [8533718, 8548042], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [8536311, 8538548], \'num_reads\': 177, \'type_junction\': 0}, {\'coords\': [8536311, 8539051], \'num_reads\': 23, \'type_junction\': 0}, {\'coords\': [8538592, 8539051], \'num_reads\': 174, \'type_junction\': 0}, {\'coords\': [8539128, 8548042], \'num_reads\': 325, \'type_junction\': 0}, {\'coords\': [8539128, 8550487], \'num_reads\': 120, \'type_junction\': 0}, {\'coords\': [8546440, 8548042], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [8548095, 8550487], \'num_reads\': 269, \'type_junction\': 0}, {\'coords\': [8550608, 8550746], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [8551289, 8551908], \'num_reads\': 209, \'type_junction\': 0}, {\'coords\': [8551959, 8552893], \'num_reads\': 3, \'type_junction\': 0}, {\'coords\': [8551959, 8553575], \'num_reads\': 172, \'type_junction\': 0}, {\'coords\': [8552950, 8553575], \'num_reads\': 0, \'type_junction\': 2}], \'name\': \'ENST00000348943\', \'strand\': \'+\'}'.replace(/\'/g, "\"").replace(/'/g, ""));
//var genes_obj = JSON.parse('{\'exons\': [{\'a3\': [], \'a5\': [1, 0], \'coords\': [6291633, 6291769], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [], \'a5\': [3, 2], \'coords\': [6291985, 6292513], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [], \'a5\': [4], \'coords\': [6292605, 6292714], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [0, 1, 2, 3, 4], \'a5\': [5], \'coords\': [6296952, 6297677], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [5], \'a5\': [6, 7], \'coords\': [6313788, 6313979], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [6], \'a5\': [8, 9], \'coords\': [6316765, 6316855], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [7], \'a5\': [10], \'coords\': [6316978, 6317080], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [8], \'a5\': [11], \'coords\': [6324878, 6324922], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [9, 10, 11], \'a5\': [12], \'coords\': [6334533, 6334648], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [12], \'a5\': [13], \'coords\': [6338558, 6338712], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [13], \'a5\': [14], \'coords\': [6339109, 6339255], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [14], \'a5\': [15], \'coords\': [6339838, 6339928], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [15], \'a5\': [16], \'coords\': [6340443, 6340622], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [16], \'a5\': [17], \'coords\': [6342495, 6342623], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [17], \'a5\': [18], \'coords\': [6346515, 6346694], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [18], \'a5\': [19], \'coords\': [6347098, 6347250], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [19], \'a5\': [20], \'coords\': [6347829, 6347931], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [20], \'a5\': [21], \'coords\': [6348235, 6348363], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [21], \'a5\': [22], \'coords\': [6348658, 6348808], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [22], \'a5\': [23], \'coords\': [6349100, 6349227], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [23], \'a5\': [24], \'coords\': [6349315, 6349493], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [24], \'a5\': [25], \'coords\': [6349842, 6349913], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [25], \'a5\': [26], \'coords\': [6350606, 6350734], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [26], \'a5\': [27], \'coords\': [6352126, 6352198], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [27], \'a5\': [28], \'coords\': [6354929, 6355092], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [28], \'a5\': [29], \'coords\': [6355291, 6355445], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [29], \'a5\': [], \'coords\': [6355554, 6356642], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}], \'junctions\': [{\'coords\': [6291693, 6296952], \'num_reads\': 79, \'type_junction\': 0}, {\'coords\': [6291769, 6296952], \'num_reads\': 6, \'type_junction\': 0}, {\'coords\': [6292291, 6296952], \'num_reads\': 3, \'type_junction\': 0}, {\'coords\': [6292513, 6296952], \'num_reads\': 1, \'type_junction\': 0}, {\'coords\': [6292714, 6296952], \'num_reads\': 2, \'type_junction\': 0}, {\'coords\': [6297199, 6313788], \'num_reads\': 68, \'type_junction\': 0}, {\'coords\': [6313979, 6316765], \'num_reads\': 22, \'type_junction\': 0}, {\'coords\': [6313979, 6316978], \'num_reads\': 16, \'type_junction\': 0}, {\'coords\': [6316855, 6324878], \'num_reads\': 28, \'type_junction\': 0}, {\'coords\': [6316855, 6334533], \'num_reads\': 8, \'type_junction\': 0}, {\'coords\': [6317080, 6334533], \'num_reads\': 25, \'type_junction\': 0}, {\'coords\': [6324922, 6334533], \'num_reads\': 14, \'type_junction\': 0}, {\'coords\': [6334648, 6338558], \'num_reads\': 61, \'type_junction\': 0}, {\'coords\': [6338712, 6339109], \'num_reads\': 62, \'type_junction\': 0}, {\'coords\': [6339255, 6339838], \'num_reads\': 70, \'type_junction\': 0}, {\'coords\': [6339928, 6340443], \'num_reads\': 42, \'type_junction\': 0}, {\'coords\': [6340622, 6342495], \'num_reads\': 26, \'type_junction\': 0}, {\'coords\': [6342623, 6346515], \'num_reads\': 74, \'type_junction\': 0}, {\'coords\': [6346694, 6347098], \'num_reads\': 56, \'type_junction\': 0}, {\'coords\': [6347250, 6347829], \'num_reads\': 97, \'type_junction\': 0}, {\'coords\': [6347931, 6348235], \'num_reads\': 94, \'type_junction\': 0}, {\'coords\': [6348363, 6348658], \'num_reads\': 23, \'type_junction\': 0}, {\'coords\': [6348808, 6349100], \'num_reads\': 62, \'type_junction\': 0}, {\'coords\': [6349227, 6349315], \'num_reads\': 40, \'type_junction\': 0}, {\'coords\': [6349493, 6349842], \'num_reads\': 80, \'type_junction\': 0}, {\'coords\': [6349913, 6350606], \'num_reads\': 72, \'type_junction\': 0}, {\'coords\': [6350734, 6352126], \'num_reads\': 78, \'type_junction\': 0}, {\'coords\': [6352198, 6354929], \'num_reads\': 69, \'type_junction\': 0}, {\'coords\': [6355092, 6355291], \'num_reads\': 57, \'type_junction\': 0}, {\'coords\': [6355445, 6355554], \'num_reads\': 78, \'type_junction\': 0}], \'name\': \'ENSMUST00000003461\', \'strand\': \'+\'}'.replace(/\'/g, "\"").replace(/'/g, ""));
//var genes_obj = JSON.parse('{\'exons\': [{\'a3\': [], \'a5\': [0], \'coords\': [137475455, 137475855], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 2}, {\'a3\': [0], \'a5\': [1], \'coords\': [137476316, 137476571], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 2}, {\'a3\': [1], \'a5\': [2], \'coords\': [137480866, 137480934], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 2}, {\'a3\': [2], \'a5\': [3], \'coords\': [137481478, 137481567], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [3], \'a5\': [4], \'coords\': [137485329, 137485486], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [4], \'a5\': [5], \'coords\': [137486434, 137486697], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [5], \'a5\': [6], \'coords\': [137488171, 137488449], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [], \'a5\': [7, 8, 9], \'coords\': [137492571, 137492956], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [7], \'a5\': [10], \'coords\': [137493331, 137493583], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [6, 8, 10], \'a5\': [11], \'coords\': [137495244, 137495288], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [9, 11], \'a5\': [12], \'coords\': [137495758, 137495862], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [12], \'a5\': [13], \'coords\': [137496580, 137496924], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [13], \'a5\': [14, 15], \'coords\': [137497485, 137497553], \'coords_extra\': [], \'intron_retention\': true, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [14, 15], \'a5\': [16, 17], \'coords\': [137497741, 137497835], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [16], \'a5\': [18, 19], \'coords\': [137498819, 137499033], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [17, 18], \'a5\': [20, 21], \'coords\': [137499276, 137499332], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [19, 20], \'a5\': [22], \'coords\': [137499776, 137499822], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [21, 22], \'a5\': [23, 24], \'coords\': [137500009, 137500102], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [23], \'a5\': [25], \'coords\': [137500403, 137500641], \'coords_extra\': [], \'intron_retention\': true, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [25], \'a5\': [26], \'coords\': [137500711, 137500855], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [24, 26], \'a5\': [27], \'coords\': [137501517, 137501651], \'coords_extra\': [], \'intron_retention\': true, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [27], \'a5\': [28, 29], \'coords\': [137501688, 137501914], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [28], \'a5\': [31, 30], \'coords\': [137502207, 137502634], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [29, 30, 31], \'a5\': [32, 33], \'coords\': [137503332, 137503767], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [32], \'a5\': [34], \'coords\': [137504159, 137504377], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [33, 34], \'a5\': [35, 36], \'coords\': [137504911, 137505047], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [], \'a5\': [37], \'coords\': [137505052, 137505077], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [35, 36, 37], \'a5\': [38], \'coords\': [137506034, 137506098], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [38], \'a5\': [39], \'coords\': [137506419, 137506601], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [39], \'a5\': [40, 41], \'coords\': [137506727, 137506849], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [40, 41], \'a5\': [42, 43], \'coords\': [137507034, 137507099], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [42], \'a5\': [44, 45, 46], \'coords\': [137507754, 137507823], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [44], \'a5\': [47, 48], \'coords\': [137508409, 137508540], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [43, 46, 47, 48, 45], \'a5\': [49], \'coords\': [137512853, 137513356], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [49], \'a5\': [], \'coords\': [137514285, 137514675], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}], \'junctions\': [{\'coords\': [137475855, 137476394], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137476571, 137480866], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137480934, 137481478], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137481567, 137485329], \'num_reads\': 2, \'type_junction\': 0}, {\'coords\': [137485486, 137486434], \'num_reads\': 1, \'type_junction\': 0}, {\'coords\': [137486697, 137488171], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137488449, 137495244], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137492956, 137493331], \'num_reads\': 16, \'type_junction\': 0}, {\'coords\': [137492956, 137495244], \'num_reads\': 119, \'type_junction\': 0}, {\'coords\': [137492956, 137495758], \'num_reads\': 71, \'type_junction\': 0}, {\'coords\': [137493583, 137495244], \'num_reads\': 1, \'type_junction\': 0}, {\'coords\': [137495288, 137495758], \'num_reads\': 173, \'type_junction\': 0}, {\'coords\': [137495862, 137496580], \'num_reads\': 304, \'type_junction\': 0}, {\'coords\': [137496757, 137497485], \'num_reads\': 173, \'type_junction\': 0}, {\'coords\': [137497553, 137497741], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137497553, 137497743], \'num_reads\': 176, \'type_junction\': 0}, {\'coords\': [137497835, 137498819], \'num_reads\': 450, \'type_junction\': 0}, {\'coords\': [137497835, 137499276], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137499033, 137499276], \'num_reads\': 23, \'type_junction\': 0}, {\'coords\': [137499033, 137499776], \'num_reads\': 257, \'type_junction\': 0}, {\'coords\': [137499332, 137499776], \'num_reads\': 10, \'type_junction\': 0}, {\'coords\': [137499332, 137500009], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137499822, 137500009], \'num_reads\': 283, \'type_junction\': 0}, {\'coords\': [137500102, 137500403], \'num_reads\': 478, \'type_junction\': 0}, {\'coords\': [137500102, 137501517], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137500641, 137500711], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137500855, 137501517], \'num_reads\': 323, \'type_junction\': 0}, {\'coords\': [137501651, 137501688], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137501797, 137502207], \'num_reads\': 51, \'type_junction\': 0}, {\'coords\': [137501797, 137503623], \'num_reads\': 7, \'type_junction\': 0}, {\'coords\': [137502299, 137503623], \'num_reads\': 5, \'type_junction\': 0}, {\'coords\': [137502416, 137503623], \'num_reads\': 108, \'type_junction\': 0}, {\'coords\': [137503767, 137504159], \'num_reads\': 140, \'type_junction\': 0}, {\'coords\': [137503767, 137504911], \'num_reads\': 61, \'type_junction\': 0}, {\'coords\': [137504377, 137504911], \'num_reads\': 86, \'type_junction\': 0}, {\'coords\': [137505038, 137506034], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137505047, 137506034], \'num_reads\': 188, \'type_junction\': 0}, {\'coords\': [137505077, 137506034], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137506098, 137506521], \'num_reads\': 232, \'type_junction\': 0}, {\'coords\': [137506601, 137506727], \'num_reads\': 366, \'type_junction\': 0}, {\'coords\': [137506849, 137507034], \'num_reads\': 5, \'type_junction\': 0}, {\'coords\': [137506849, 137507050], \'num_reads\': 287, \'type_junction\': 0}, {\'coords\': [137507099, 137507754], \'num_reads\': 268, \'type_junction\': 0}, {\'coords\': [137507099, 137513260], \'num_reads\': 2, \'type_junction\': 0}, {\'coords\': [137507823, 137508409], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137507823, 137513255], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137507823, 137513260], \'num_reads\': 184, \'type_junction\': 0}, {\'coords\': [137508515, 137513260], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [137508540, 137513260], \'num_reads\': 1, \'type_junction\': 0}, {\'coords\': [137513356, 137514285], \'num_reads\': 498, \'type_junction\': 0}], \'name\': \'ENST00000254900\', \'strand\': \'-\'}'.replace(/\'/g, "\"").replace(/'/g, ""));
var genes_obj = JSON.parse('{\'exons\': [{\'a3\': [], \'a5\': [0, 1], \'coords\': [135502453, 135502674], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [0, 1], \'a5\': [2, 3], \'coords\': [135502940, 135507158], \'coords_extra\': [[135502940, 135507014]], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [2], \'a5\': [4], \'coords\': [135508972, 135509043], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [3, 4], \'a5\': [5, 6], \'coords\': [135510929, 135511021], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [5], \'a5\': [7], \'coords\': [135511265, 135511485], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [6, 7], \'a5\': [8], \'coords\': [135513462, 135513696], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [8], \'a5\': [9, 10], \'coords\': [135514976, 135515056], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [9], \'a5\': [11, 12, 13, 14, 15], \'coords\': [135515494, 135515824], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [11, 13], \'a5\': [16], \'coords\': [135516098, 135516219], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [10, 12, 14, 16], \'a5\': [21, 22, 23, 20, 17, 18, 19], \'coords\': [135516886, 135517140], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [17, 21], \'a5\': [24], \'coords\': [135517864, 135518046], \'coords_extra\': [], \'intron_retention\': true, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [18, 20, 22, 24], \'a5\': [25], \'coords\': [135518099, 135518461], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [15, 19, 23, 25], \'a5\': [26, 27], \'coords\': [135520046, 135520188], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [26], \'a5\': [28], \'coords\': [135520664, 135520719], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [27, 28], \'a5\': [29], \'coords\': [135521223, 135521337], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [29], \'a5\': [30, 31, 33, 32], \'coords\': [135521428, 135521812], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [30, 32, 33], \'a5\': [34, 35], \'coords\': [135522777, 135522887], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [34], \'a5\': [36, 37], \'coords\': [135523552, 135523807], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [36], \'a5\': [38], \'coords\': [135524086, 135524087], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 2}, {\'a3\': [31, 35, 37, 38], \'a5\': [39, 40], \'coords\': [135524355, 135524462], \'coords_extra\': [], \'intron_retention\': true, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [39], \'a5\': [], \'coords\': [135524854, 135525088], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [40], \'a5\': [], \'coords\': [135539002, 135540311], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}], \'junctions\': [{\'coords\': [135502674, 135502940], \'num_reads\': 3, \'type_junction\': 1}, {\'coords\': [135502674, 135507041], \'num_reads\': 534, \'type_junction\': 0}, {\'coords\': [135507158, 135508972], \'num_reads\': 487, \'type_junction\': 0}, {\'coords\': [135507158, 135510929], \'num_reads\': 249, \'type_junction\': 0}, {\'coords\': [135509043, 135510929], \'num_reads\': 1055, \'type_junction\': 0}, {\'coords\': [135511021, 135511265], \'num_reads\': 904, \'type_junction\': 0}, {\'coords\': [135511021, 135513462], \'num_reads\': 30, \'type_junction\': 0}, {\'coords\': [135511485, 135513462], \'num_reads\': 393, \'type_junction\': 0}, {\'coords\': [135513696, 135514976], \'num_reads\': 692, \'type_junction\': 0}, {\'coords\': [135515056, 135515494], \'num_reads\': 501, \'type_junction\': 0}, {\'coords\': [135515056, 135516886], \'num_reads\': 22, \'type_junction\': 0}, {\'coords\': [135515589, 135516098], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135515589, 135516886], \'num_reads\': 34, \'type_junction\': 0}, {\'coords\': [135515598, 135516098], \'num_reads\': 4, \'type_junction\': 0}, {\'coords\': [135515598, 135516886], \'num_reads\': 600, \'type_junction\': 0}, {\'coords\': [135515598, 135520046], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135516219, 135516886], \'num_reads\': 4, \'type_junction\': 0}, {\'coords\': [135517055, 135517864], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135517055, 135518099], \'num_reads\': 1, \'type_junction\': 0}, {\'coords\': [135517055, 135520046], \'num_reads\': 3, \'type_junction\': 0}, {\'coords\': [135517092, 135518099], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135517140, 135517864], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135517140, 135518099], \'num_reads\': 2, \'type_junction\': 0}, {\'coords\': [135517140, 135520046], \'num_reads\': 207, \'type_junction\': 0}, {\'coords\': [135518046, 135518099], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135518461, 135520046], \'num_reads\': 1, \'type_junction\': 0}, {\'coords\': [135520188, 135520664], \'num_reads\': 14, \'type_junction\': 0}, {\'coords\': [135520188, 135521223], \'num_reads\': 429, \'type_junction\': 0}, {\'coords\': [135520719, 135521223], \'num_reads\': 14, \'type_junction\': 0}, {\'coords\': [135521337, 135521428], \'num_reads\': 380, \'type_junction\': 0}, {\'coords\': [135521553, 135522777], \'num_reads\': 365, \'type_junction\': 0}, {\'coords\': [135521553, 135524355], \'num_reads\': 1, \'type_junction\': 0}, {\'coords\': [135521695, 135522777], \'num_reads\': 5, \'type_junction\': 0}, {\'coords\': [135521812, 135522777], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135522887, 135523552], \'num_reads\': 4, \'type_junction\': 0}, {\'coords\': [135522887, 135524355], \'num_reads\': 743, \'type_junction\': 0}, {\'coords\': [135523807, 135524086], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135523807, 135524355], \'num_reads\': 1, \'type_junction\': 0}, {\'coords\': [135524087, 135524355], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135524462, 135524854], \'num_reads\': 2, \'type_junction\': 0}, {\'coords\': [135524462, 135539002], \'num_reads\': 535, \'type_junction\': 0}], \'name\': \'ENST00000339290\', \'strand\': \'+\'}'.replace(/\'/g, "\"").replace(/'/g, ""));
var exons_obj = genes_obj.exons;
var junctions_obj = genes_obj.junctions;

var orig_objs = {'exons': add_keys(clone(exons_obj)), 'junc': clone(junctions_obj)};

var exons_mapped = map_exon_list(exons_obj, junctions_obj); //exons_obj; //
exons_mapped = add_keys(exons_mapped);

var padding = [10, 5, 5, 5];
var w = 1000,
    h = 200,
    JUNC_AREA=0.9,
    EXON_H = Math.round(h * (1-JUNC_AREA) - padding[2]),
    EXON_MIN_W= 2;

var svgCanvas = d3.select("#testDiv")
    .append("svg")
    .attr("width", w)
    .attr("height", h);

svgCanvas.append('clipPath')
    .attr("id", "cut-off-junctions")
    .append("rect")
    .attr("x", 0)
    .attr("y", 0)
    .attr("width", w)
    .attr("height", h * JUNC_AREA);

/** Render initial splice graph */
renderSpliceGraph(exons_mapped, junctions_obj);
