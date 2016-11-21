/**
 * Created by abarrera on 8/27/14.
 */
/** Utils */

function clone(obj) {
    if (obj === null || typeof(obj) != 'object')
        return obj;

    var temp = obj.constructor(); // changed

    for (var key in obj) {
        if (obj.hasOwnProperty(key)) {
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

function add(a, b) {
    return a + b;
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

    var reshaped_intron = constant_size(exon2.start - exon1.end);
    if (reduce_exon) {
        if ((exon1.end - exon1.start) > MAX_INTRON) {
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
    var k;
    var i;
    var j;
    var reshape_exons = false;
    var exons_mapped_tmp = [];
    var exon_tmp;
    var index_exons_to_update = [];  // For exons_obj that need to be updated (shorted)

    // First offset (first exon can have long intron)

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

    var coords = map(function (x) {
        return add(x, -acc_offset);
    }, resize_exon([exons[0].start, exons[0].end], reshape_exons));

    exon_tmp = {
        'start': coords[0],
        'end': coords[1],
        'exon_type': exons[0].exon_type,
        'intron_retention': exons[0].intron_retention,
        'lsv_type': exons[0].lsv_type,
        'a3': exons[0].a3,
        'a5': exons[0].a5
    };
    exon_tmp.coords_extra = coords_extra;

    exons_mapped_tmp[0] = exon_tmp;

    for (i = 0; i < exons[0].a3.length; i++) {
        junctions[exons[0].a3[i]].end -= offset;
    }

    for (i = 0; i < exons[0].a5.length; i++) {
        junctions[exons[0].a5[i]].start -= offset;
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

        coords = map(function (x) {
            return add(x, -acc_offset);
        }, resize_exon([exons[i].start, exons[i].end], reshape_exons));

        exon_tmp = {
            'start': coords[0],
            'end': coords[1],
            'exon_type': exons[i].exon_type,
            'intron_retention': exons[i].intron_retention,
            'lsv_type': exons[i].lsv_type,
            'a3': exons[i].a3,
            'a5': exons[i].a5
        };
        exon_tmp.coords_extra = coords_extra;

        exons_mapped_tmp[i] = exon_tmp;

        // Check if any previous exon needs to be updated
        if (exons[i].start - exons[i - 1].end < 0) {
            index_exons_to_update.push(i - 1);
        }

        if (offset) {
            for (k = 0; k < index_exons_to_update.length; k++) {
                if (exons[index_exons_to_update[k]].end > exons[i].start) {
                    exons_mapped_tmp[index_exons_to_update[k]].end -= offset;
                }
            }
        }

        for (j = 0; j < exons[i].a3.length; j++) {
            junctions[exons[i].a3[j]].end -= acc_offset;
        }
        for (j = 0; j < exons[i].a5.length; j++) {
            junctions[exons[i].a5[j]].start -= acc_offset;
        }
    }
    return exons_mapped_tmp;

}

var add_keys = function (elems) {
    var elems_dkeys = [];
    for (var i = 0; i < elems.length; i++) {
        elems_dkeys.push({'key': i, 'value': elems[i]})
    }
    return elems_dkeys;
};

/** End utils */


var mostLeftCoord = function (exons) {
    return d3.min(exons, function (d) {
        return d.value.start;
    });
};
var mostRightCoord = function (exons) {
    return d3.max(exons, function (d) {
        return d.value.end;
    });
};

var longestJunc = function (junctions, scale) {
    var longest = -1;
    for (var i = 0; i < junctions.length; i++) {
        if ((scale(junctions[i].end) - scale(junctions[i].start)) > longest) {
            longest = scale(junctions[i].end) - scale(junctions[i].start);
        }
    }
    return longest;
};

var exonKey = function (d) {
    return d.key;
};


function toolTipD3(strand, begin, end, el) {
    // if anything is selected, don't change the tool tip
    if ($(el).closest('.splice-div-container').find('.selected').length)
        return;

    if (d3.select(el).classed('halfexon'))
        if (d3.select(el).classed('missingStart'))
            begin = 'MISSING';
        else
            end = 'MISSING';

    var length = (end - begin) + 1;
    var toolTip = d3.select(el.parentNode.parentNode.parentNode.parentNode).select('.tooltipD3');

    if (strand === '-')
        toolTip.select('.coordsLabel').text(end + ' - ' + begin);
    else
        toolTip.select('.coordsLabel').text(begin + ' - ' + end);

    toolTip.select('.lengthLabel').text(isNaN(length) ? 'UNKNOWN' : length);
}

function spliceGraphD3() {

    var width = 1000; // default width
    var height = 160; // default height
    var padding = [60, 5, 5, 5];
    var JUNC_AREA = 0.8;
    var EXON_H = Math.round(height * (1 - JUNC_AREA) - padding[2]);
    var EXON_MIN_W = 2;
    var orig_objs;

    //This is the accessor function we talked about above
    var lineFunction = d3.svg.line()
        .x(function (d) {
            return d.x;
        })
        .y(function (d) {
            return d.y;
        })
        .interpolate("linear");


    function my(selection) {

        selection.each(function (d) {

            // Select the svg element, if it exists.
            var svgCanvas = d3.select(this).selectAll("svg").data([d]);

            // Otherwise, create the skeletal chart.
            svgCanvas.enter().append("svg");

            svgCanvas.attr("width", width).attr("height", height);

            var clipP = svgCanvas.selectAll('.cut-off-junctions').data([1]);
            clipP.enter().append('svg:clipPath');
            clipP.classed("cut-off-junctions", true);
            clipP.attr("id", "cut-off-junctions-" + this.id);
            var rectClip = clipP.selectAll("rect").data([1]);
            rectClip.enter().append("rect");
            rectClip
                .attr("x", 0)
                .attr("y", 0)
                .attr("width", width)
                .attr("height", height * JUNC_AREA);


            // Add dispersion factors to the junctions
            var addDispersion = function (exons, junctions) {

                for (var ii = 0; ii < exons.length; ii++) {
                    var d = exons[ii];
                    for (var i = 0; i < d.value.a3.length; i++) {
                        junctions[d.value.a3[i]].dispFrom = junctions[d.value.a3[0]].start;
                        if (i == 0) {
                            junctions[d.value.a3[i]].dispersion = 1;
                        } else {
                            junctions[d.value.a3[i]].dispersion = d.value.a3.length - i + 1;
                        }
                    }
                }

            };

            var renderNumReads = function (junctions, scaleX, displayCorrectedReads) {

                var maxJunc = longestJunc(junctions, scaleX);

                var labels = svgCanvas.selectAll("text.readcounts")
                    .data(junctions.filter(function (v) {
                        return v.intron_retention < 1;
                    }));

                labels.enter().append("text");

                labels
                    .attr("class", "readcounts")
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .text(function (d) {
                        if (d.reads) {
                            if (displayCorrectedReads)
                                return d.clean_reads;
                            return d.reads;
                        }
                        return '';
                    })
                    .attr("x", function (d) {
                        if (strand == '-')
                            return Math.round((2 * width - scaleX(d.start) - scaleX(d.end)) / 2);
                        return Math.round((scaleX(d.end) + scaleX(d.start)) / 2);
                    })
                    .attr("y", function (d) {
                        var posY = (scaleX(d.end) - scaleX(d.start)) / maxJunc * JUNC_AREA * (height - padding[0] - padding[2]);
                        return Math.round(height * JUNC_AREA - 2 - posY + (posY / d.dispersion) * (d.dispersion - 1 ? 1 : 0));
                    })
                    .style("fill", displayCorrectedReads ? "red" : "black");
                return labels;
            };


            var renderJunctions = function (data, scaleX, docId) {
                var juncs = svgCanvas.selectAll('ellipse').data(data.filter(function (v) {
                    return v.intron_retention < 1;
                }));  // Skip intron retention junctions
                var maxJunc = longestJunc(data, scaleX);
                juncs.enter().append("ellipse");

                juncs.attr("class", "junction")
                    .attr("style", "")
                    .attr("clip-path", "url(#" + "cut-off-junctions-" + docId + ")")
                    .transition()
                    .duration(100)
                    .ease("linear")

                    .attr("ry", function (d) {
                        return Math.round(
                            (scaleX(d.end) - scaleX(d.start)) / maxJunc * JUNC_AREA * (height - padding[0] - padding[2])
                            - ((scaleX(d.end) - scaleX(d.start)) / maxJunc * JUNC_AREA * (height - padding[0] - padding[2]) / d.dispersion) * (d.dispersion - 1 ? 1 : 0)
                        );

                    })
                    .attr("rx", function (d) {
                        return Math.round((scaleX(d.end) - scaleX(d.start)) / 2);
                    })
                    .attr("cy", height * JUNC_AREA)
                    .attr("cx", function (d) {
                        return Math.round((scaleX(d.end) + scaleX(d.start)) / 2);
                    });

                juncs.classed("found", function (d) {
                    return d.junction_type == 0;
                });
                juncs.classed("novel", function (d) {
                    return d.junction_type == 1;
                });
                juncs.classed("missing", function (d) {
                    return d.junction_type == 2;
                });
                juncs.classed("foundInOthers", function (d) {
                    return d.junction_type == 3;
                });
                juncs.classed("novelInOthers", function (d) {
                    return d.junction_type == 1 && d.reads < 1;
                });

                return juncs;
            };


            var renderIntRetReads = function (data, scaleX) {
                var labels = svgCanvas.selectAll("text.irreads")
                    .data(data.filter(function (v) {
                        return v.intron_retention > 0;
                    }));
                labels.enter().append("text");
                labels.attr("class", "irreads")
                    .attr("text-anchor", function (d) {
                        return ((d.ir === 1 && strand == '-' || d.ir === 2 && strand == '+' ) ? "end" : "start" );
                    })
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .text(function (d) {
                        return d.num_reads;
                    })
                    .attr("x", function (d) {
                        var pos = Math.round(scaleX(d.coords[(d.ir == 1 ? 1 : 0)]) + (d.ir == 1 ? 2 : -2));
                        if (strand == '-')
                            return Math.round(width - pos);
                        else
                            return pos;
                    })
                    .attr("y", Math.round(height * JUNC_AREA + EXON_H / 5));
                return labels;
            };

            var renderIntRetLines = function (data, scaleX) {
                var irlines = svgCanvas.selectAll("line.irlines")
                    .data(data.filter(function (v) {
                        return v.intron_retention > 0;
                    }));
                irlines.enter().append("line");
                irlines.attr("class", "irlines")
                    .attr("style", "")
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .attr("x1", function (d) {
                        return Math.round(scaleX(d.coords[(d.ir == 1 ? 1 : 0)]));
                    })
                    .attr("y1", function () {
                        return Math.round(height * JUNC_AREA + EXON_H / 5 + 2);
                    })
                    .attr("x2", function (d) {
                        return Math.round(scaleX(d.coords[(d.ir == 1 ? 1 : 0)]) + ((d.ir == 1) ? 14 : -14));
                    })
                    .attr("y2", function () {
                        return Math.round(height * JUNC_AREA + EXON_H / 5 + 2);
                    });
                irlines.classed("found", function (d) {
                    return d.junction_type == 0;
                });
                irlines.classed("novel", function (d) {
                    return d.junction_type == 1;
                });
                return irlines;
            };

            var spliceSites = function (junctions, scaleX) {
                var ssites3 = svgCanvas
                    .selectAll("line.ssite3")
                    .data(junctions);

                ssites3
                    .enter()
                    .append("line");

                ssites3
                    .classed("ssite3", true)
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .attr("x1", function (d) {
                        return Math.round(scaleX(d.start));
                    })
                    .attr("y1", function () {
                        return Math.round(height * JUNC_AREA);
                    })
                    .attr("x2", function (d) {
                        return Math.round(scaleX(d.start));
                    })
                    .attr("y2", function () {
                        return Math.round(height - padding[2]);
                    });

                ssites3
                    .classed("found", function (d) {
                        return d.junction_type == 0;
                    })
                    .classed("novel", function (d) {
                        return d.junction_type == 1;
                    })
                    .classed("missing", function (d) {
                        return d.junction_type == 2;
                    });

                var ssites5 = svgCanvas
                    .selectAll("line.ssite5")
                    .data(junctions);

                ssites5
                    .enter()
                    .append("line");

                ssites5
                    .classed("ssite5", true)
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .attr("x1", function (d) {
                        return Math.round(scaleX(d.end));
                    })
                    .attr("y1", function () {
                        return Math.round(height * JUNC_AREA);
                    })
                    .attr("x2", function (d) {
                        return Math.round(scaleX(d.end));
                    })
                    .attr("y2", function () {
                        return Math.round(height - padding[2]);
                    });

                ssites5
                    .classed("found", function (d) {
                        return d.junction_type == 0;
                    })
                    .classed("novel", function (d) {
                        return d.junction_type == 1;
                    })
                    .classed("missing", function (d) {
                        return d.junction_type == 2;
                    });

            };

            var renderExons = function (datap, keyf, scaleX) {
                // Represent exons as rectangles - filter out half exons
                var exons = svgCanvas.selectAll("rect.exon").data(datap.filter(function (v) {
                    return v.value.exon_type < 3 && !v.value.intron_retention;
                }), keyf);

                exons.enter().append("rect");
                exons.attr("class", "exon")
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .attr("x", function (d) {
                        return scaleX(d.value.start);
                    })
                    .attr("y", height * JUNC_AREA)
                    .attr("width", function (d) {
                        return Math.max(EXON_MIN_W, Math.round(scaleX(d.value.end) - scaleX(d.value.start)));
                    })
                    .attr("height", EXON_H);

                exons.each(function (d) {
                    if (d.value.exon_type == 0) {
                        d3.select(this).classed("found", true);
                    }
                    if (d.value.exon_type == 1) {
                        d3.select(this).classed("novel", true);
                    }
                    if (d.value.exon_type == 2) {
                        d3.select(this).classed("missing", true);
                    }
                });

                return exons;
            };

            var renderHalfExons = function (datap, keyf, scaleX) {

                var halfExons = svgCanvas.selectAll("path.halfexon").data(datap.filter(function (v) {
                    return v.value.exon_type > 2;
                }), keyf);
                halfExons.enter().append("path");
                halfExons.attr("class", "halfexon")
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .attr("d", function (d) {
                        var index_missing = 0;  // missing start
                        if (d.value.exon_type == 4) //missing end
                            index_missing = 1;
                        var halfExonsPoints = [
                            {'x': scaleX(d.value.coords[(index_missing + 1) % 2]), 'y': height * JUNC_AREA},
                            {'x': scaleX(d.value.coords[index_missing]), 'y': height * JUNC_AREA},
                            {'x': scaleX(d.value.coords[index_missing]), 'y': height * JUNC_AREA + EXON_H},
                            {'x': scaleX(d.value.coords[(index_missing + 1) % 2]), 'y': height * JUNC_AREA + EXON_H}
                        ];
                        return lineFunction(halfExonsPoints);
                    });

                halfExons.each(function (d) {
                    if (d.value.exon_type == 3) {
                        d3.select(this).classed("missingStart", true);
                    }
                    if (d.value.exon_type == 4) {
                        d3.select(this).classed("missingEnd", true);
                    }
                });
                return halfExons;
            };

            var renderNumExons = function (exons, scaleX) {
                var exons_only = exons.filter(function (v) {
                    return v.value.exon_type < 3 && !v.value.intron_retention;
                });
                var labels = svgCanvas.selectAll("text.numexon")
                    .data(exons_only);
                labels.enter().append("text");
                labels.classed("numexon", true)
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .text(function (d, i) {
                        if (strand == '-')
                            return (exons_only.length - (i)).toString();
                        return (i + 1).toString();
                    })
                    .attr("x", function (exon) {
                        if (strand == '-')
                            return Math.round((2 * width - scaleX(exon.value.start) - scaleX(exon.value.end)) / 2);
                        return Math.round((scaleX(exon.value.end) + scaleX(exon.value.start)) / 2);
                    })
                    .attr("y", height * JUNC_AREA + EXON_H / 2)
                    .attr("alignment-baseline", "central");

            };

            var renderIntronRetention = function (exons, scaleX) {
                var intronsRet = svgCanvas.selectAll("rect.intronret")
                    .data(exons.filter(function (v) {
                        return v.value.intron_retention;
                    }));  // Only exons with intron retention
                intronsRet.enter().append("rect");
                intronsRet
                    .attr("class", "intronret")
                    .attr("style", "")
                    .classed('missing', function (d) {
                        return d.value.exon_type == 2;
                    })
                    .classed('novel', function (d) {
                        return d.value.exon_type == 1;
                    })
                    .transition()
                    .duration(100)
                    .ease("linear")
                    .attr("x", function (d) {
                        return scaleX(d.value.coords[0]);
                    })
                    .attr("y", Math.round(height * JUNC_AREA + EXON_H * 2 / 5))
                    .attr("width", function (d) {
                        return Math.round(scaleX(d.value.coords[1]) - scaleX(d.value.coords[0]));
                    })
                    .attr("height", Math.round(EXON_H * 2 / 5));
                return intronsRet;
            };

            var renderCoordsExtra = function (exons, scaleX) {
                var coords_extra = [];
                exons.forEach(function (e) {
                    if (e.value.coords_extra.length > 0 && e.value.exon_type != 2) {
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
                    .attr("x", function (d) {
                        return scaleX(d[0]);
                    })
                    .attr("y", Math.round(height * JUNC_AREA))
                    .attr("width", function (d) {
                        return Math.round(scaleX(d[1]) - scaleX(d[0]));
                    })
                    .attr("height", Math.round(EXON_H));
                partialNewExons.exit().remove();
            };

            var updateScale = function (datap) {
                return d3.scale.linear()
                    .domain([mostLeftCoord(datap), mostRightCoord(datap)])
                    .rangeRound([padding[3], width - padding[1]]);
            };

            var highlightLSV = function (spliceDiv, strand, d3Els) {
                var lsvID;
                var coords;
                var origCoords;
                var row;
                var halfExonsLength;
                var colorEls;
                var isSource;
                var wls;
                var irArray;
                var alreadySelected;
                var c;
                var w;
                var isWeighted;

                var highlightLSVs = spliceDiv.parentNode.parentNode.querySelectorAll('[value="highlight"]:checked');
                var colors = ['red', 'blue', 'green', 'purple', 'orange', 'brown', 'pink', 'grey'];

                var exons = orig_objs.exons.filter(function (d) {
                    return !d.value.intron_retention;
                });

                var junctions = orig_objs.junc.filter(function (d) {
                    return d.intron_retention == 0;
                });

                var intronRetentions = orig_objs.exons.filter(function (d) {
                    return d.value.intron_retention;
                });

                var irLines = orig_objs.junc.filter(function (d) {
                    return d.intron_retention > 0;
                });

                var weightedLines = function (row) {
                    var index = d3.select(spliceDiv).classed('exp1') ? 3 : 5;
                    var jsonData = row[index].querySelector('canvas').getAttribute('data-lsv');
                    var lsv_list = JSON.parse(jsonData.replace(/'/g, '"'));
                    var trans_lsv = translate_lsv_bins(lsv_list[0].bins, 1000);
                    return trans_lsv.reduce(function (prev, curr) {
                        return prev.concat((d3.mean(curr) * 3).toFixed(3));
                    }, []);
                };

                if (highlightLSVs.length) {
                    d3Els.junctions.classed('highlight-lsv-blurred', true);
                    d3Els.numReads.classed('highlight-lsv-blurred', true);
                    d3Els.irReads.classed('highlight-lsv-blurred', true);
                    d3Els.irLines.classed('highlight-lsv-blurred', true);
                }

                highlightLSVs.forEach(function (highlighLSV) {
                    row = highlighLSV.parentNode.parentNode.parentNode.parentNode.querySelectorAll('td');
                    isWeighted = row[0].querySelector('.weighted:checked');
                    lsvID = row[1].textContent.split(':');
                    coords = lsvID[1].split("-");
                    isSource = lsvID[2] === "source";
                    colorEls = [];
                    wls = weightedLines(row);

                    exons.forEach(function (exon, index) {
                        origCoords = [exon.value.start, exon.value.end];
                        if (origCoords[0] == coords[0] && origCoords[1] == coords[1]) {
                            // number of half exons before this exon
                            halfExonsLength = d3Els.halfExons.filter(function (d) {
                                return d.key < exon.key
                            })[0].length;
                            var exon_index = index - halfExonsLength;
                            d3.select(d3Els.exons[0][exon_index]).classed('highlight-lsv', true);
                        }
                    });

                    junctions.forEach(function (junction, index) {
                        origCoords = [junction.start, junction.end];
                        if (coords[0] <= origCoords[isSource ? 0 : 1] &&
                            coords[1] >= origCoords[isSource ? 0 : 1]) {
                            colorEls.push([d3Els.junctions[0][index]]);
                            d3.select(d3Els.numReads[0][index]).classed('highlight-lsv-blurred', false);
                        }
                    });

                    // junctions in genes with strand '-' are colored backwards
                    if (strand == '-')
                        colorEls.reverse();

                    // ir and ir lines are colored at the same time
                    irArray = [];

                    intronRetentions.forEach(function (ir, index) {
                        if (isSource && end == ir.value.start - 1 ||
                            !isSource && start == ir.value.end + 1) {
                            irArray.push(d3Els.intronRets[0][index]);
                        }
                    });

                    irLines.forEach(function (irLine, index) {
                        if (!isSource && start == irLine.end ||
                            isSource && end == irLine.start) {
                            d3.select(d3Els.irReads[0][index]).classed('highlight-lsv-blurred', false);
                            irArray.push(d3Els.irLines[0][index]);
                        }
                    });

                    colorEls.push(irArray);

                    colorEls.forEach(function (junction, index) {
                        d3.selectAll(junction)
                            .classed('highlight-lsv-blurred', false)
                            .attr("style", function () {
                                alreadySelected = d3.select(this).classed('highlight-lsv');
                                c = alreadySelected ? "black" : colors[index];
                                w = !alreadySelected && isWeighted ? ";stroke-width:" + wls[index] : "";
                                return 'stroke:' + c + ';fill:' + c + w + ";";
                            })
                            .classed('highlight-lsv', true);
                    })
                })
            };

            // are we displaying normalized reads?
            var toggleReadCounts = this.parentNode.parentNode.querySelector('.readCounts');
            var displayCorrectedReads = toggleReadCounts ? toggleReadCounts.checked : false;

            // generate chart here, using `w` and `h`
            var exonsp = d[0],
                junctionsp = d[1],
                strand = d[2];

            /** Compute the scale used */
            var scaleX = updateScale(exonsp);

            /** Render exons and compute disperison for junctions */
            var exons = renderExons(exonsp, exonKey, scaleX);
            addDispersion(exonsp, junctionsp);
            renderNumExons(exonsp, scaleX);
            var introns_ret = renderIntronRetention(exonsp, scaleX);
            renderCoordsExtra(exonsp, scaleX);

            /** Render half exons */
            var halfExons = renderHalfExons(exonsp, exonKey, scaleX);

            /** Render junctions and read numbers */
            var junctions = renderJunctions(junctionsp, scaleX, this.id);
            var numReads = renderNumReads(junctionsp, scaleX, displayCorrectedReads);
            var irReads = renderIntRetReads(junctionsp, scaleX);
            var irLines = renderIntRetLines(junctionsp, scaleX);
            spliceSites(junctionsp, scaleX);

            /** Add interactivity for ... */
            [exons, introns_ret, halfExons].forEach(function (el) {
                el
                    .on('mouseover', function (d) {
                        d3.select(this).classed("hovered", true);
                        toolTipD3(strand, orig_objs.exons[d.key].value.start, orig_objs.exons[d.key].value.end, this);
                    })
                    .on('mouseout', function () {
                        d3.select(this).classed("hovered", false);
                    });
            });

            var juncs_orig_filtered = orig_objs.junc.filter(function (v) {
                return v.intron_retention < 1;
            });

            /** ..and junctions! */
            junctions
                .on('mouseover', function (d, i) {
                    var d3svg = d3.select(this.parentNode);
                    var d3this = d3.select(this);
                    d3svg.selectAll('.junction').classed("blurred", true);
                    d3svg.selectAll('.readcounts').classed("blurred", true);
                    d3.select(d3svg.selectAll('.ssite3')[0][i]).classed("highlighted", true);
                    d3.select(d3svg.selectAll('.ssite5')[0][i]).classed("highlighted", true);
                    d3.select(d3svg.selectAll('.readcounts')[0][i]).classed("blurred", false).classed("highlighted", true);
                    d3this.classed("blurred", false);
                    d3this.classed("hovered", true);
                    toolTipD3(strand, juncs_orig_filtered[i].start, juncs_orig_filtered[i].end, this);
                })
                .on('mouseout', function () {
                    d3.selectAll('.junction').classed('blurred', false);
                    d3.selectAll('.readcounts').classed("blurred", false).classed("highlighted", false);
                    d3.selectAll('.ssite3').classed("highlighted", false);
                    d3.selectAll('.ssite5').classed("highlighted", false);
                    d3.select(this).classed("hovered", false);
                });

            /**
             * Add mouse click.
             */
            [junctions, exons, introns_ret, halfExons].forEach(function (els) {
                els.on('click', function () {
                    var className = 'selected';
                    var isThisSelected = d3.select(this).classed(className);
                    var spliceDivContainer = $(this).closest('.splice-div-container').get(0);
                    var event = document.createEvent('SVGEvents');

                    // un-select all elements
                    d3.select(spliceDivContainer).selectAll('.junction, .exon, .intronret, .halfexon').classed(className, false);

                    // trigger mouseover event while nothing's selected to update tool tip
                    event.initEvent('mouseover', true, true);
                    this.dispatchEvent(event);

                    // toggle selected on this element
                    d3.select(this).classed(className, !isThisSelected);
                });
            });

            if (strand == '-')
                d3.select(this).selectAll(":not(text)").attr("transform", "translate(" + width + ",0) scale(-1 , 1)");

            var d3Elements = {
                'exons': exons,
                'halfExons': halfExons,
                'junctions': junctions,
                'numReads': numReads,
                'intronRets': introns_ret,
                'irLines': irLines,
                'irReads': irReads
            };

            highlightLSV(this, strand, d3Elements);

        });


        $('#sg-filters').submit();

    }

    my.width = function (value) {
        if (!arguments.length) return width;
        width = Math.max(value, 1000);
        return my;
    };

    my.height = function (value) {
        if (!arguments.length) return height;
        height = Math.max(Math.min(value, 600), 120);
        EXON_H = Math.round(height * (1 - JUNC_AREA) - padding[2]);
        return my;
    };

    my.orig_objs = function (value) {
        if (!arguments.length) return orig_objs;
        orig_objs = value;
        return my;
    };

    return my;

}
