/**
 * Created by abarrera on 8/27/14.
 */
/** Utils */

function clone(obj) {
    if(obj === null || typeof(obj) != 'object')
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

    exons_mapped_tmp[0] = exon_tmp;
    var last_end = exon_tmp.coords[1];

    for (var i = 0; i < exons[0].a3.length; i++) {
        junctions[exons[0].a3[i]].coords[1] -= offset;
    }

    for (var i = 0; i < exons[0].a5.length; i++) {
        junctions[exons[0].a5[i]].coords[0] -= offset;
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

        for (var j = 0; j < exons[i].a3.length; j++) {
            junctions[exons[i].a3[j]].coords[1] -= acc_offset;
        }
        for (var j = 0; j < exons[i].a5.length; j++) {
            junctions[exons[i].a5[j]].coords[0] -= acc_offset;
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



function spliceGraphD3() {


    var width = 1000, // default width
        height = 120, // default height
        padding = [10, 5, 5, 5],
        JUNC_AREA=0.85;
    var EXON_H = Math.round(height * (1-JUNC_AREA) - padding[2]),
        EXON_MIN_W= 2,
        EXON_MAX_W=200;
    var orig_objs;

    function my(selection) {

        selection.each(function(d, i){
            // Select the svg element, if it exists.
            var svgCanvas = d3.select(this).selectAll("svg").data([d]);

            // Otherwise, create the skeletal chart.
            svgCanvas.enter().append("svg");

            svgCanvas.attr("width", width).attr("height", height);

            var clipP = svgCanvas.selectAll('#cut-off-junctions').data([1]);
            clipP.enter().append('svg:clipPath');
            clipP.attr("id", "cut-off-junctions");
            var rectClip = clipP.selectAll("rect").data([1]);
            rectClip.enter().append("rect");
            rectClip
                .attr("x", 0)
                .attr("y", 0)
                .attr("width", width)
                .attr("height", height * JUNC_AREA);

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

            var renderNumReads = function(junctions, scaleX) {
                var maxJunc = longestJunc(junctions, scaleX);
                var labels = svgCanvas.selectAll("text.readcounts")
                    .data(junctions);
                labels.enter().append("text");
                labels.classed("readcounts", true)
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .text(function (d) {
                        if (d.num_reads) {
                            return d.num_reads;
                        }
                        return '';
                    })
                    .attr("x", function (d) {
                        if (strand == '-')
                            return Math.round((2*width - scaleX(d.coords[0]) - scaleX(d.coords[1])) / 2);
                        return Math.round((scaleX(d.coords[1]) + scaleX(d.coords[0])) / 2);
                    })
                    .attr("y", function (d) {
                        return Math.round(height * JUNC_AREA - 2 -
                                (scaleX(d.coords[1]) - scaleX(d.coords[0])) / maxJunc * JUNC_AREA * (height - padding[0] - padding[2])
                                + ((scaleX(d.coords[1]) - scaleX(d.coords[0])) / maxJunc * JUNC_AREA * (height - padding[0] - padding[2]) / d.dispersion) * (d.dispersion - 1 ? 1 : 0)
                        );
                    });

            };

            var renderJunctions = function(data, scaleX){
                var juncs = svgCanvas.selectAll('ellipse').data(data);
                var maxJunc = longestJunc(data, scaleX);
                juncs.enter().append("ellipse");

                juncs.attr("class", "junction")
                    .attr("clip-path", "url(#" + "cut-off-junctions)")
                    .transition()
                    .duration(100)
                    .ease("linear")

                    .attr("ry", function(d){
                        return Math.round(
                                (scaleX(d.coords[1]) - scaleX(d.coords[0]))/maxJunc * JUNC_AREA * (height - padding[0] -padding[2])
                                - ((scaleX(d.coords[1]) - scaleX(d.coords[0]))/maxJunc * JUNC_AREA * (height - padding[0] -padding[2])/d.dispersion)  * (d.dispersion - 1 ? 1 : 0)
                        );
                    })
                    .attr("rx", function(d){
                        return Math.round((scaleX(d.coords[1]) - scaleX(d.coords[0]))/2);
                    })
                    .attr("cy", height * JUNC_AREA)
                    .attr("cx", function(d){
                        return Math.round((scaleX(d.coords[1]) + scaleX(d.coords[0]))/2);
                    });
//                    .attr("transform", function(d){
//                        return "translate(" + width + ",0) scale(-1 , 1)";
//                    });

                juncs.classed("found", function(d){ return d.type_junction == 0; });
                juncs.classed("novel", function(d){ return d.type_junction == 1; });
                juncs.classed("missing", function(d){ return d.type_junction == 2; });

                return juncs;
            };

            var spliceSites = function(junctions, scaleX) {
                var ssites3 = svgCanvas.selectAll("line.ssite3")
                    .data(junctions);
                ssites3.enter().append("line");
                ssites3.classed("ssite3", true)
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .attr("x1", function (d) {
                        return Math.round(scaleX(d.coords[0]));
                    })
                    .attr("y1", function (d) {
                        return Math.round(height * JUNC_AREA);
                    })
                    .attr("x2", function (d) {
                        return Math.round(scaleX(d.coords[0]));
                    })
                    .attr("y2", function (d) {
                        return Math.round(height - padding[2]);
                    });
                ssites3.classed("found", function(d){ return d.type_junction == 0; });
                ssites3.classed("novel", function(d){ return d.type_junction == 1; });
                ssites3.classed("missing", function(d){ return d.type_junction == 2; });

                var ssites5 = svgCanvas.selectAll("line.ssite5")
                    .data(junctions);
                ssites5.enter().append("line");
                ssites5.classed("ssite5", true)
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .attr("x1", function (d) {
                        return Math.round(scaleX(d.coords[1]));
                    })
                    .attr("y1", function (d) {
                        return Math.round(height * JUNC_AREA);
                    })
                    .attr("x2", function (d) {
                        return Math.round(scaleX(d.coords[1]));
                    })
                    .attr("y2", function (d) {
                        return Math.round(height - padding[2]);
                    });
                ssites5.classed("found", function(d){ return d.type_junction == 0; });
                ssites5.classed("novel", function(d){ return d.type_junction == 1; });
                ssites5.classed("missing", function(d){ return d.type_junction == 2; });

            };

            var renderExons = function(datap, keyf, scaleX){

                var exons = svgCanvas.selectAll("rect.exon").data(datap, keyf);

                exons.enter().append("rect");
                exons.attr("class", "exon")
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .attr("x", function(d){
                        return scaleX(d.value.coords[0]);
                    })
                    .attr("y", height*JUNC_AREA)
                    .attr("width", function(d){
                        return Math.max(EXON_MIN_W, Math.round(scaleX(d.value.coords[1]) - scaleX(d.value.coords[0])));
                    })
                    .attr("height", EXON_H);

                return exons;
            };

            var renderNumExons = function(exons, scaleX) {
                var labels = svgCanvas.selectAll("text.numexon")
                    .data(exons);
                labels.enter().append("text");
                labels.classed("numexon", true)
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .text(function (d, i) {
                        if (strand == '-')
                            return (exons.length - (i)).toString();
                        return (i + 1).toString();
                    })
                    .attr("x", function (d, i) {
                        if (strand == '-')
                            return Math.round((2*width - scaleX(d.value.coords[0]) - scaleX(d.value.coords[1])) / 2);
                        return Math.round((scaleX(d.value.coords[1]) + scaleX(d.value.coords[0])) / 2);
                    })
                    .attr("y", height*JUNC_AREA + EXON_H /2)
                    .attr("alignment-baseline", "central");

            };

            var renderIntronRetention = function(exons, scaleX) {
                var intronsRet = svgCanvas.selectAll("rect.intronret")
                    .data(exons.filter(function(v){ return v.value.intron_retention;}));  // Only exons with intron retention
                intronsRet.enter().append("rect");
                intronsRet.classed("intronret", true)
                    .classed('missing', function(d){
                        return d.value.type_exon == 2;
                    })
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .attr("x", function(d){
                        return scaleX(d.value.coords[1]);
                    })
                    .attr("y", Math.round(height*JUNC_AREA + height*(1-JUNC_AREA)/5))
                    .attr("width", function(d){
                        return Math.round(scaleX(exons[d.key+1].value.coords[0]) - scaleX(d.value.coords[1]));
                    })
                    .attr("height", Math.round((EXON_H)*2/5));
            };

            var renderCoordsExtra = function(exons, scaleX) {
                var coords_extra = [];
                exons.forEach(function(e){
                    if (e.value.coords_extra.length>0 && e.value.type_exon != 2) {
                        e.value.coords_extra.forEach(function (ce) {
                            coords_extra.push(ce);
                        });
                    }
                });

                var partialNewExons = svgCanvas.selectAll("rect.newpartialexon").data(coords_extra);  // Only exons with coords extra
                partialNewExons.enter().append("rect");
                partialNewExons.classed("newpartialexon", true)
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .attr("x", function(d){
                        return scaleX(d[0]);
                    })
                    .attr("y", Math.round(height*JUNC_AREA))
                    .attr("width", function(d){
                        return Math.round(scaleX(d[1]) - scaleX(d[0]));
                    })
                    .attr("height", Math.round(EXON_H));
                partialNewExons.exit().remove();
            };

            var updateScale = function(datap) {
                return d3.scale.linear()
                    .domain([mostLeftCoord(datap), mostRightCoord(datap)])
                    .rangeRound([padding[3], width - padding[1]]);

            };

            // generate chart here, using `w` and `h`
            var exonsp = d[0],
                junctionsp=d[1],
                strand=d[2];

            /** Compute the scale used */
            var scaleX = updateScale(exonsp);

            /** Render exons and compute disperison for junctions */
            var exons = renderExons(exonsp, exonKey, scaleX);
            addDispersion(exons, junctionsp);
            renderNumExons(exonsp, scaleX);
            renderIntronRetention(exonsp, scaleX);
            renderCoordsExtra(exonsp, scaleX);

            /** Render junctions and read numbers */
            var junctions = renderJunctions(junctionsp, scaleX);
            renderNumReads(junctionsp, scaleX);
            spliceSites(junctionsp, scaleX);

            /** Add interactivity for exons ... */
            exons.on('mouseover', function(d){  // Highlight exons when hovered
                d3.select(this).classed("hovered", true);
                // Tooltip
                //Get this x/y values, then augment for the tooltip
                var xPosition = parseFloat(d3.select(this).attr("x")) ;
                var yPosition = height + 120;

                //Update the tooltip position and value
                d3.select(this.parentNode.parentNode).select(".tooltipD3")
//                    .style("left", xPosition + "px")
//                    .style("top", yPosition + "px")
                    .select(".coordsLabel")
                    .text(function(){
                        return orig_objs.exons[d.key].value.coords[0] + "-" + orig_objs.exons[d.key].value.coords[1];
                    });

                //Show the tooltip
                d3.select(this.parentNode.parentNode).select(".tooltipD3").classed("hidden", false);

            })
                .on('mouseout', function(d){
                    d3.select(this).classed("hovered", false);
                    d3.select(this.parentNode.parentNode).select(".tooltipD3").classed("hidden", true);
                });

            /** ..and junctions! */
            junctions.on('mouseover', function(d, i){
                d3.select(this.parentNode).selectAll('.junction').classed("blurred", true);
                d3.select(this.parentNode).selectAll('.readcounts').classed("blurred", true);
                d3.select(d3.select(this.parentNode).selectAll('.ssite3')[0][i]).classed("highlighted", true);
                d3.select(d3.select(this.parentNode).selectAll('.ssite5')[0][i]).classed("highlighted", true);
                d3.select(d3.select(this.parentNode).selectAll('.readcounts')[0][i]).classed("blurred", false).classed("highlighted", true);
                d3.select(this).classed("blurred", false);
                d3.select(this).classed("hovered", true);

                //Update the tooltip position and value
                d3.select(this.parentNode.parentNode).select(".tooltipD3")
//                    .style("left", mouseCoords[0]+ "px")
//                    .style("top", mouseCoords[1]+ "px")
                    .classed("hidden", false)
                    .select(".coordsLabel")
                    .text(function(){
                        return orig_objs.junc[i].coords[0] + "-" + orig_objs.junc[i].coords[1];
                    });
                //Show the tooltip
//                d3.select(".tooltipD3").classed("hidden", false);

            })
            .on('mouseout', function(){
                d3.selectAll('.junction').classed('blurred', false);
                d3.selectAll('.readcounts').classed("blurred", false).classed("highlighted", false);
                d3.selectAll('.ssite3').classed("highlighted", false);
                d3.selectAll('.ssite5').classed("highlighted", false);
                d3.select(this).classed("hovered", false);

                //Show the tooltip
                d3.select(this.parentNode.parentNode).select(".tooltipD3").classed("hidden", true);
            });

//            d3.selectAll(this.childNodes[0].childNodes).attr("transform", "translate(" + width + ",0) scale(-1 , 1)");
//            d3.select(this.childNodes[0]).selectAll("text").attr("transform", "translate(" + 10 + ",0)");
            if (strand == '-')
                d3.select(this).selectAll(":not(text)").attr("transform", "translate(" + width + ",0) scale(-1 , 1)");
        });

    }

    my.width = function(value) {
        if (!arguments.length) return width;
        width = Math.max(value, 1000);
        return my;
    };

    my.height = function(value) {
        if (!arguments.length) return height;
        height = Math.max(Math.min(value, 600), 120);
        EXON_H = Math.round(height * (1-JUNC_AREA) - padding[2]);
        return my;
    };

    my.orig_objs = function(value) {
        if (!arguments.length) return orig_objs;
        orig_objs = value;
        return my;
    };

    return my;

}


/**
* D3 - SpliceGraph - ONLY ACTIVE WHILE DEBUGGING!!
* */

//var genes_obj = JSON.parse('{\'exons\': [{\'a3\': [], \'a5\': [0, 1], \'coords\': [135502453, 135502674], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [0, 1], \'a5\': [2, 3], \'coords\': [135502940, 135507158], \'coords_extra\': [[135502940, 135507014]], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [2], \'a5\': [4], \'coords\': [135508972, 135509043], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [3, 4], \'a5\': [5, 6], \'coords\': [135510929, 135511021], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [5], \'a5\': [7], \'coords\': [135511265, 135511485], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [6, 7], \'a5\': [8], \'coords\': [135513462, 135513696], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [8], \'a5\': [9, 10], \'coords\': [135514976, 135515056], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [9], \'a5\': [11, 12, 13, 14, 15], \'coords\': [135515494, 135515824], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [11, 13], \'a5\': [16], \'coords\': [135516098, 135516219], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [10, 12, 14, 16], \'a5\': [21, 22, 23, 20, 17, 18, 19], \'coords\': [135516886, 135517140], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [17, 21], \'a5\': [24], \'coords\': [135517864, 135518046], \'coords_extra\': [], \'intron_retention\': true, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [18, 20, 22, 24], \'a5\': [25], \'coords\': [135518099, 135518461], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [15, 19, 23, 25], \'a5\': [26, 27], \'coords\': [135520046, 135520188], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [26], \'a5\': [28], \'coords\': [135520664, 135520719], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [27, 28], \'a5\': [29], \'coords\': [135521223, 135521337], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [29], \'a5\': [30, 31, 33, 32], \'coords\': [135521428, 135521812], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [30, 32, 33], \'a5\': [34, 35], \'coords\': [135522777, 135522887], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [34], \'a5\': [36, 37], \'coords\': [135523552, 135523807], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [36], \'a5\': [38], \'coords\': [135524086, 135524087], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 2}, {\'a3\': [31, 35, 37, 38], \'a5\': [39, 40], \'coords\': [135524355, 135524462], \'coords_extra\': [], \'intron_retention\': true, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [39], \'a5\': [], \'coords\': [135524854, 135525088], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}, {\'a3\': [40], \'a5\': [], \'coords\': [135539002, 135540311], \'coords_extra\': [], \'intron_retention\': false, \'lsv_type\': 0, \'type_exon\': 0}], \'junctions\': [{\'coords\': [135502674, 135502940], \'num_reads\': 3, \'type_junction\': 1}, {\'coords\': [135502674, 135507041], \'num_reads\': 534, \'type_junction\': 0}, {\'coords\': [135507158, 135508972], \'num_reads\': 487, \'type_junction\': 0}, {\'coords\': [135507158, 135510929], \'num_reads\': 249, \'type_junction\': 0}, {\'coords\': [135509043, 135510929], \'num_reads\': 1055, \'type_junction\': 0}, {\'coords\': [135511021, 135511265], \'num_reads\': 904, \'type_junction\': 0}, {\'coords\': [135511021, 135513462], \'num_reads\': 30, \'type_junction\': 0}, {\'coords\': [135511485, 135513462], \'num_reads\': 393, \'type_junction\': 0}, {\'coords\': [135513696, 135514976], \'num_reads\': 692, \'type_junction\': 0}, {\'coords\': [135515056, 135515494], \'num_reads\': 501, \'type_junction\': 0}, {\'coords\': [135515056, 135516886], \'num_reads\': 22, \'type_junction\': 0}, {\'coords\': [135515589, 135516098], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135515589, 135516886], \'num_reads\': 34, \'type_junction\': 0}, {\'coords\': [135515598, 135516098], \'num_reads\': 4, \'type_junction\': 0}, {\'coords\': [135515598, 135516886], \'num_reads\': 600, \'type_junction\': 0}, {\'coords\': [135515598, 135520046], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135516219, 135516886], \'num_reads\': 4, \'type_junction\': 0}, {\'coords\': [135517055, 135517864], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135517055, 135518099], \'num_reads\': 1, \'type_junction\': 0}, {\'coords\': [135517055, 135520046], \'num_reads\': 3, \'type_junction\': 0}, {\'coords\': [135517092, 135518099], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135517140, 135517864], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135517140, 135518099], \'num_reads\': 2, \'type_junction\': 0}, {\'coords\': [135517140, 135520046], \'num_reads\': 207, \'type_junction\': 0}, {\'coords\': [135518046, 135518099], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135518461, 135520046], \'num_reads\': 1, \'type_junction\': 0}, {\'coords\': [135520188, 135520664], \'num_reads\': 14, \'type_junction\': 0}, {\'coords\': [135520188, 135521223], \'num_reads\': 429, \'type_junction\': 0}, {\'coords\': [135520719, 135521223], \'num_reads\': 14, \'type_junction\': 0}, {\'coords\': [135521337, 135521428], \'num_reads\': 380, \'type_junction\': 0}, {\'coords\': [135521553, 135522777], \'num_reads\': 365, \'type_junction\': 0}, {\'coords\': [135521553, 135524355], \'num_reads\': 1, \'type_junction\': 0}, {\'coords\': [135521695, 135522777], \'num_reads\': 5, \'type_junction\': 0}, {\'coords\': [135521812, 135522777], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135522887, 135523552], \'num_reads\': 4, \'type_junction\': 0}, {\'coords\': [135522887, 135524355], \'num_reads\': 743, \'type_junction\': 0}, {\'coords\': [135523807, 135524086], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135523807, 135524355], \'num_reads\': 1, \'type_junction\': 0}, {\'coords\': [135524087, 135524355], \'num_reads\': 0, \'type_junction\': 2}, {\'coords\': [135524462, 135524854], \'num_reads\': 2, \'type_junction\': 0}, {\'coords\': [135524462, 135539002], \'num_reads\': 535, \'type_junction\': 0}], \'name\': \'ENST00000339290\', \'strand\': \'+\'}'.replace(/\'/g, "\"").replace(/'/g, ""));
//var exons_obj = genes_obj.exons;
//var junctions_obj = genes_obj.junctions;
//var strand = "+";
//
//var orig_objs = {'exons': add_keys(clone(exons_obj)), 'junc': clone(junctions_obj)};
//
//var exons_mapped = map_exon_list(exons_obj, junctions_obj); //exons_obj; //
//exons_mapped = add_keys(exons_mapped);
//
//
///** Render initial splice graph */
//var chart = spliceGraphD3().orig_objs(orig_objs);
//var spliceg = d3.select("#testDiv")
//    .datum([exons_mapped, junctions_obj, strand])
//    .call(chart);
//
//d3.select('.zoomInSplice').on('click', function() {
//    chart.width(chart.width() + 600);
//    chart.height(chart.height() + 100);
//    spliceg.call(chart);
//});
//
//d3.select('.zoomOutSplice').on('click', function() {
//    chart.width(chart.width() - 600);
//    chart.height(chart.height() - 100);
//    spliceg.call(chart);
//});
//
//d3.select('.zoomResetSplice').on('click', function() {
//    chart.width(1000);
//    chart.height(120);
//    spliceg.call(chart);
//});
//
//d3.select('.toogleScale').on('click', function(){
////    chart.width(chart.width() + 500);
//    if (d3.select(this).classed('scaled')) {
//        spliceg.datum([orig_objs.exons, orig_objs.junc, strand])
//            .call(chart);
//        d3.select(this).classed('scaled', false);
//    } else {
//        spliceg.datum([exons_mapped, junctions_obj, strand])
//            .call(chart);
//        d3.select(this).classed('scaled', true);
//    }
//});

