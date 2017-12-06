/*
 *
 * TableSorter 2.0 - Client-side table sorting with ease!
 * Version 2.0.5b
 * @requires jQuery v1.2.3
 *
 * Copyright (c) 2007 Christian Bach
 * Examples and docs at: http://tablesorter.com
 * Dual licensed under the MIT and GPL licenses:
 * http://www.opensource.org/licenses/mit-license.php
 * http://www.gnu.org/licenses/gpl.html
 *
 */
/**
 *
 * @description Create a sortable table with multi-column sorting capabilitys
 *
 * @example $('table').tablesorter();
 * @desc Create a simple tablesorter interface.
 *
 * @example $('table').tablesorter({ sortList:[[0,0],[1,0]] });
 * @desc Create a tablesorter interface and sort on the first and secound column column headers.
 *
 * @example $('table').tablesorter({ headers: { 0: { sorter: false}, 1: {sorter: false} } });
 *
 * @desc Create a tablesorter interface and disableing the first and second  column headers.
 *
 *
 * @example $('table').tablesorter({ headers: { 0: {sorter:"integer"}, 1: {sorter:"currency"} } });
 *
 * @desc Create a tablesorter interface and set a column parser for the first
 *       and second column.
 *
 *
 * @param Object
 *            settings An object literal containing key/value pairs to provide
 *            optional settings.
 *
 *
 * @option String cssHeader (optional) A string of the class name to be appended
 *         to sortable tr elements in the thead of the table. Default value:
 *         "header"
 *
 * @option String cssAsc (optional) A string of the class name to be appended to
 *         sortable tr elements in the thead on a ascending sort. Default value:
 *         "headerSortUp"
 *
 * @option String cssDesc (optional) A string of the class name to be appended
 *         to sortable tr elements in the thead on a descending sort. Default
 *         value: "headerSortDown"
 *
 * @option String sortInitialOrder (optional) A string of the inital sorting
 *         order can be asc or desc. Default value: "asc"
 *
 * @option String sortMultisortKey (optional) A string of the multi-column sort
 *         key. Default value: "shiftKey"
 *
 * @option String textExtraction (optional) A string of the text-extraction
 *         method to use. For complex html structures inside td cell set this
 *         option to "complex", on large tables the complex option can be slow.
 *         Default value: "simple"
 *
 * @option Object headers (optional) An array containing the forces sorting
 *         rules. This option let's you specify a default_view sorting rule. Default
 *         value: null
 *
 * @option Array sortList (optional) An array containing the forces sorting
 *         rules. This option let's you specify a default_view sorting rule. Default
 *         value: null
 *
 * @option Array sortForce (optional) An array containing forced sorting rules.
 *         This option let's you specify a default_view sorting rule, which is
 *         prepended to user-selected rules. Default value: null
 *
 * @option Boolean sortLocaleCompare (optional) Boolean flag indicating whatever
 *         to use String.localeCampare method or not. Default set to true.
 *
 *
 * @option Array sortAppend (optional) An array containing forced sorting rules.
 *         This option let's you specify a default_view sorting rule, which is
 *         appended to user-selected rules. Default value: null
 *
 * @option Boolean widthFixed (optional) Boolean flag indicating if tablesorter
 *         should apply fixed widths to the table columns. This is usefull when
 *         using the pager companion plugin. This options requires the dimension
 *         jquery plugin. Default value: false
 *
 * @option Boolean cancelSelection (optional) Boolean flag indicating if
 *         tablesorter should cancel selection of the table headers text.
 *         Default value: true
 *
 * @option Boolean debug (optional) Boolean flag indicating if tablesorter
 *         should display debuging information usefull for development.
 *
 * @type jQuery
 *
 * @name tablesorter
 *
 * @cat Plugins/Tablesorter
 *
 * @author Christian Bach/christian.bach@polyester.se
 */

(function ($) {
    $.extend({
        tablesorter: new function () {

            var parsers = [],
                widgets = [];

            this.defaults = {
                cssHeader: "header",
                cssAsc: "headerSortUp",
                cssDesc: "headerSortDown",
                cssChildRow: "expand-child",
                sortInitialOrder: "asc",
                sortMultiSortKey: "shiftKey",
                sortForce: null,
                sortAppend: null,
                sortLocaleCompare: true,
                textExtraction: "simple",
                parsers: {}, widgets: [],
                widgetZebra: {
                    css: ["even", "odd"]
                }, headers: {}, widthFixed: false,
                cancelSelection: true,
                sortList: [],
                headerList: [],
                dateFormat: "us",
                decimal: '/\.|\,/g',
                onRenderHeader: null,
                selectorHeaders: 'thead th',
                debug: false
            };

            /* debuging utils */

            function benchmark(s, d) {
                log(s + "," + (new Date().getTime() - d.getTime()) + "ms");
            }

            this.benchmark = benchmark;

            function log(s) {
                if (typeof console != "undefined" && typeof console.debug != "undefined") {
                    console.log(s);
                } else {
                    alert(s);
                }
            }

            /* parsers utils */

            function buildParserCache(table, $headers) {

                if (table.config.debug) {
                    var parsersDebug = "";
                }

                if (table.tBodies.length == 0) return; // In the case of empty tables
                var rows = table.tBodies[0].rows;

                if (rows[0]) {

                    var list = [],
                        cells = rows[0].cells,
                        l = cells.length;

                    for (var i = 0; i < l; i++) {

                        var p = false;

                        if ($.metadata && ($($headers[i]).metadata() && $($headers[i]).metadata().sorter)) {

                            p = getParserById($($headers[i]).metadata().sorter);

                        } else if ((table.config.headers[i] && table.config.headers[i].sorter)) {

                            p = getParserById(table.config.headers[i].sorter);
                        }
                        if (!p) {

                            p = detectParserForColumn(table, rows, -1, i);
                        }

                        if (table.config.debug) {
                            parsersDebug += "column:" + i + " parser:" + p.id + "\n";
                        }

                        list.push(p);
                    }
                }

                if (table.config.debug) {
                    log(parsersDebug);
                }

                return list;
            }

            function detectParserForColumn(table, rows, rowIndex, cellIndex) {
                var l = parsers.length,
                    node = false,
                    nodeValue = false,
                    keepLooking = true;
                while (nodeValue == '' && keepLooking) {
                    rowIndex++;
                    if (rows[rowIndex]) {
                        node = getNodeFromRowAndCellIndex(rows, rowIndex, cellIndex);
                        nodeValue = trimAndGetNodeText(table.config, node);
                        if (table.config.debug) {
                            log('Checking if value was empty on row:' + rowIndex);
                        }
                    } else {
                        keepLooking = false;
                    }
                }
                for (var i = 1; i < l; i++) {
                    if (parsers[i].is(nodeValue, table, node)) {
                        return parsers[i];
                    }
                }
                // 0 is always the generic parser (text)
                return parsers[0];
            }

            function getNodeFromRowAndCellIndex(rows, rowIndex, cellIndex) {
                return rows[rowIndex].cells[cellIndex];
            }

            function trimAndGetNodeText(config, node) {
                return $.trim(getElementText(config, node));
            }

            function getParserById(name) {
                var l = parsers.length;
                for (var i = 0; i < l; i++) {
                    if (parsers[i].id.toLowerCase() == name.toLowerCase()) {
                        return parsers[i];
                    }
                }
                return false;
            }

            /* utils */

            function buildCache(table) {

                if (table.config.debug) {
                    var cacheTime = new Date();
                }

                var totalRows = (table.tBodies[0] && table.tBodies[0].rows.length) || 0,
                    totalCells = (table.tBodies[0].rows[0] && table.tBodies[0].rows[0].cells.length) || 0,
                    parsers = table.config.parsers,
                    cache = {
                        row: [],
                        normalized: []
                    };

                for (var i = 0; i < totalRows; ++i) {

                    /** Add the table data to main data array */
                    var c = $(table.tBodies[0].rows[i]),
                        cols = [];

                    // if this is a child row, add it to the last row's children and
                    // continue to the next row
                    if (c.hasClass(table.config.cssChildRow)) {
                        cache.row[cache.row.length - 1] = cache.row[cache.row.length - 1].add(c);
                        // go to the next for loop
                        continue;
                    }

                    cache.row.push(c);

                    for (var j = 0; j < totalCells; ++j) {
                        cols.push(parsers[j].format(getElementText(table.config, c[0].cells[j]), table, c[0].cells[j]));
                    }

                    cols.push(cache.normalized.length); // add position for rowCache
                    cache.normalized.push(cols);
                    cols = null;
                }
                if (table.config.debug) {
                    benchmark("Building cache for " + totalRows + " rows:", cacheTime);
                }

                return cache;
            }

            function getElementText(config, node) {

                var text = "";

                if (!node) return "";

                if (!config.supportsTextContent) config.supportsTextContent = node.textContent || false;

                if (config.textExtraction == "simple") {
                    if (config.supportsTextContent) {
                        text = node.textContent;
                    } else {
                        if (node.childNodes[0] && node.childNodes[0].hasChildNodes()) {
                            text = node.childNodes[0].innerHTML;
                        } else {
                            text = node.innerHTML;
                        }
                    }
                } else {
                    if (typeof(config.textExtraction) == "function") {
                        text = config.textExtraction(node);
                    } else {
                        text = $(node).text();
                    }
                }
                return text;
            }

            function appendToTable(table, cache) {

                if (table.config.debug) {
                    var appendTime = new Date()
                }

                var c = cache,
                    r = c.row,
                    n = c.normalized,
                    totalRows = n.length,
                    checkCell = (n[0].length - 1),
                    tableBody = $(table.tBodies[0]),
                    rows = [];


                for (var i = 0; i < totalRows; i++) {
                    var pos = n[i][checkCell];

                    rows.push(r[pos]);

                    if (!table.config.appender) {

                        //var o = ;
                        var l = r[pos].length;
                        for (var j = 0; j < l; j++) {
                            tableBody[0].appendChild(r[pos][j]);
                        }

                        //
                    }
                }


                if (table.config.appender) {

                    table.config.appender(table, rows);
                }

                rows = null;

                if (table.config.debug) {
                    benchmark("Rebuilt table:", appendTime);
                }

                // apply table widgets
                //applyWidget(table); // DEACTIVATED BECAUSE SINCE tablesorter.pager RENDER WIDGETS BY DEFAULT!!!

                // trigger sortend
                setTimeout(function () {
                    $(table).trigger("sortEnd");
                }, 0);

            }

            function buildHeaders(table) {

                if (table.config.debug) {
                    var time = new Date();
                }

                var meta = ($.metadata) ? true : false;

                var header_index = computeTableHeaderCellIndexes(table);

                $tableHeaders = $(table.config.selectorHeaders, table).each(function (index) {

                    this.column = header_index[this.parentNode.rowIndex + "-" + this.cellIndex];
                    // this.column = index;
                    this.order = formatSortingOrder(table.config.sortInitialOrder);


                    this.count = this.order;

                    if (checkHeaderMetadata(this) || checkHeaderOptions(table, index)) this.sortDisabled = true;
                    if (checkHeaderOptionsSortingLocked(table, index)) this.order = this.lockedOrder = checkHeaderOptionsSortingLocked(table, index);

                    if (!this.sortDisabled) {
                        var $th = $(this).addClass(table.config.cssHeader);
                        if (table.config.onRenderHeader) table.config.onRenderHeader.apply($th);
                    }

                    // add cell to headerList
                    table.config.headerList[index] = this;
                });

                if (table.config.debug) {
                    benchmark("Built headers:", time);
                    log($tableHeaders);
                }

                return $tableHeaders;

            }

            // from:
            // http://www.javascripttoolbox.com/lib/table/examples.php
            // http://www.javascripttoolbox.com/temp/table_cellindex.html


            function computeTableHeaderCellIndexes(t) {
                var matrix = [];
                var lookup = {};
                var thead = t.getElementsByTagName('THEAD')[0];
                var trs = thead.getElementsByTagName('TR');
                var k;

                for (var i = 0; i < trs.length; i++) {
                    var cells = trs[i].cells;
                    for (var j = 0; j < cells.length; j++) {
                        var c = cells[j];

                        var rowIndex = c.parentNode.rowIndex;
                        var cellId = rowIndex + "-" + c.cellIndex;
                        var rowSpan = c.rowSpan || 1;
                        var colSpan = c.colSpan || 1;
                        var firstAvailCol;
                        if (typeof(matrix[rowIndex]) == "undefined") {
                            matrix[rowIndex] = [];
                        }
                        // Find first available column in the first row
                        for (k = 0; k < matrix[rowIndex].length + 1; k++) {
                            if (typeof(matrix[rowIndex][k]) == "undefined") {
                                firstAvailCol = k;
                                break;
                            }
                        }
                        lookup[cellId] = firstAvailCol;
                        for (k = rowIndex; k < rowIndex + rowSpan; k++) {
                            if (typeof(matrix[k]) == "undefined") {
                                matrix[k] = [];
                            }
                            var matrixrow = matrix[k];
                            for (var l = firstAvailCol; l < firstAvailCol + colSpan; l++) {
                                matrixrow[l] = "x";
                            }
                        }
                    }
                }
                return lookup;
            }

            function checkCellColSpan(table, rows, row) {
                var arr = [],
                    r = table.tHead.rows,
                    c = r[row].cells;

                for (var i = 0; i < c.length; i++) {
                    var cell = c[i];

                    if (cell.colSpan > 1) {
                        arr = arr.concat(checkCellColSpan(table, headerArr, row++));
                    } else {
                        if (table.tHead.length == 1 || (cell.rowSpan > 1 || !r[row + 1])) {
                            arr.push(cell);
                        }
                    }
                }
                return arr;
            }

            function checkHeaderMetadata(cell) {
                return !!(($.metadata) && ($(cell).metadata().sorter === false));
            }

            function checkHeaderOptions(table, i) {
                return !!((table.config.headers[i]) && (table.config.headers[i].sorter === false));
            }

            function checkHeaderOptionsSortingLocked(table, i) {
                if ((table.config.headers[i]) && (table.config.headers[i].lockedOrder)) return table.config.headers[i].lockedOrder;
                return false;
            }

            function applyWidget(table) {
                var c = table.config.widgets;
                var l = c.length;
                for (var i = 0; i < l; i++) {

                    getWidgetById(c[i]).format(table);
                }

            }

            function getWidgetById(name) {
                var l = widgets.length;
                for (var i = 0; i < l; i++) {
                    if (widgets[i].id.toLowerCase() == name.toLowerCase()) {
                        return widgets[i];
                    }
                }
            }

            function formatSortingOrder(v) {
                if (typeof(v) !== "number") {
                    return (v.toLowerCase() == "desc") ? 1 : 0;
                } else {
                    return (v == 1) ? 1 : 0;
                }
            }

            function isValueInArray(v, a) {
                var l = a.length;
                for (var i = 0; i < l; i++) {
                    if (a[i][0] == v) {
                        return true;
                    }
                }
                return false;
            }

            function setHeadersCss(table, $headers, list, css) {
                // remove all header information
                $headers.removeClass(css[0]).removeClass(css[1]);

                var h = [];
                $headers.each(function (offset) {
                    if (!this.sortDisabled) {
                        h[this.column] = $(this);
                    }
                });

                var l = list.length;
                for (var i = 0; i < l; i++) {
                    h[list[i][0]].addClass(css[list[i][1]]);
                }
            }

            function fixColumnWidth(table, $headers) {
                var c = table.config;
                if (c.widthFixed) {
                    var colgroup = $('<colgroup>');
                    $("tr:first td", table.tBodies[0]).each(function () {
                        colgroup.append($('<col>').css('width', $(this).width()));
                    });
                    $(table).prepend(colgroup);
                }
            }

            function updateHeaderSortCount(table, sortList) {
                var c = table.config,
                    l = sortList.length;
                for (var i = 0; i < l; i++) {
                    var s = sortList[i],
                        o = c.headerList[s[0]];
                    o.count = s[1];
                    o.count++;
                }
            }

            /* sorting methods */
            function multisort(table, sortList, cache) {
                var i;

                if (table.config.debug) {
                    var sortTime = new Date();
                }

                var dynamicExp = "var sortWrapper = function(a,b) {",
                    l = sortList.length;


                for (i = 0; i < l; i++) {
                    var c = sortList[i][0];
                    var order = sortList[i][1];
                    var s = (table.config.parsers[c].type == "text") ? ((order == 0) ? makeSortFunction("text", "asc", c) : makeSortFunction("text", "desc", c)) : ((order == 0) ? makeSortFunction("numeric", "asc", c) : makeSortFunction("numeric", "desc", c));
                    var e = "e" + i;

                    dynamicExp += "var " + e + " = " + s;
                    dynamicExp += "if(" + e + ") { return " + e + "; } ";
                    dynamicExp += "else { ";
                }

                // if value is the same keep orignal order
                var orgOrderCol = cache.normalized[0].length - 1;
                dynamicExp += "return a[" + orgOrderCol + "]-b[" + orgOrderCol + "];";

                for (i = 0; i < l; i++) {
                    dynamicExp += "}; ";
                }

                dynamicExp += "return 0; ";
                dynamicExp += "}; ";

                if (table.config.debug) {
                    benchmark("Evaling expression:" + dynamicExp, new Date());
                }

                eval(dynamicExp);

                cache.normalized.sort(sortWrapper);

                if (table.config.debug) {
                    benchmark("Sorting on " + sortList.toString() + " and dir " + order + " time:", sortTime);
                }

                return cache;
            }

            function makeSortFunction(type, direction, index) {
                var a = "a[" + index + "]",
                    b = "b[" + index + "]";
                if (type == 'text' && direction == 'asc') {
                    return "(" + a + " == " + b + " ? 0 : (" + a + " === null ? Number.POSITIVE_INFINITY : (" + b + " === null ? Number.NEGATIVE_INFINITY : (" + a + " < " + b + ") ? -1 : 1 )));";
                } else if (type == 'text' && direction == 'desc') {
                    return "(" + a + " == " + b + " ? 0 : (" + a + " === null ? Number.POSITIVE_INFINITY : (" + b + " === null ? Number.NEGATIVE_INFINITY : (" + b + " < " + a + ") ? -1 : 1 )));";
                } else if (type == 'numeric' && direction == 'asc') {
                    return "(" + a + " === null && " + b + " === null) ? 0 :(" + a + " === null ? Number.POSITIVE_INFINITY : (" + b + " === null ? Number.NEGATIVE_INFINITY : " + a + " - " + b + "));";
                } else if (type == 'numeric' && direction == 'desc') {
                    return "(" + a + " === null && " + b + " === null) ? 0 :(" + a + " === null ? Number.POSITIVE_INFINITY : (" + b + " === null ? Number.NEGATIVE_INFINITY : " + b + " - " + a + "));";
                }
            }

            /* public methods */
            this.construct = function (settings) {
                return this.each(function () {
                    var j;
                    // if no thead or tbody quit.
                    if (!this.tHead || !this.tBodies) return;
                    // declare
                    var $this, $document, $headers, cache, config, shiftDown = 0,
                        sortOrder;
                    // new blank config object
                    this.config = {};
                    // merge and extend.
                    config = $.extend(this.config, $.tablesorter.defaults, settings);
                    // store common expression for speed
                    $this = $(this);
                    // save the settings where they read
                    $.data(this, "tablesorter", config);
                    // build headers
                    $headers = buildHeaders(this);
                    // try to auto detect column type, and store in tables config
                    this.config.parsers = buildParserCache(this, $headers);
                    // build the cache for the tbody cells
                    cache = buildCache(this);
                    // get the css class names, could be done else where.
                    var sortCSS = [config.cssDesc, config.cssAsc];
                    // fixate columns if the users supplies the fixedWidth option
                    fixColumnWidth(this);
                    // apply event handling to headers
                    // this is to big, perhaps break it out?
                    $headers.click(
                        function (e) {
                            var totalRows = ($this[0].tBodies[0] && $this[0].tBodies[0].rows.length) || 0;
                            if (!this.sortDisabled && totalRows > 0) {
                                // Only call sortStart if sorting is
                                // enabled.
                                $this.trigger("sortStart");
                                // store exp, for speed
                                var $cell = $(this);
                                // get current column index
                                var i = this.column;
                                // get current column sort order
                                this.order = this.count++ % 2;
                                // always sort on the locked order.
                                if (this.lockedOrder) this.order = this.lockedOrder;

                                // user only whants to sort on one
                                // column
                                if (!e[config.sortMultiSortKey]) {
                                    // flush the sort list
                                    config.sortList = [];
                                    if (config.sortForce != null) {
                                        var a = config.sortForce;
                                        for (j = 0; j < a.length; j++) {
                                            if (a[j][0] != i) {
                                                config.sortList.push(a[j]);
                                            }
                                        }
                                    }
                                    // add column to sort list
                                    config.sortList.push([i, this.order]);
                                    // multi column sorting
                                } else {
                                    // the user has clicked on an all
                                    // ready sortet column.
                                    if (isValueInArray(i, config.sortList)) {
                                        // revers the sorting direction
                                        // for all tables.
                                        for (j = 0; j < config.sortList.length; j++) {
                                            var s = config.sortList[j],
                                                o = config.headerList[s[0]];
                                            if (s[0] == i) {
                                                o.count = s[1];
                                                o.count++;
                                                s[1] = o.count % 2;
                                            }
                                        }
                                    } else {
                                        // add column to sort list array
                                        config.sortList.push([i, this.order]);
                                    }
                                }
                                setTimeout(function () {
                                    // set css for headers
                                    setHeadersCss($this[0], $headers, config.sortList, sortCSS);
                                    appendToTable(
                                        $this[0], multisort(
                                            $this[0], config.sortList, cache)
                                    );
                                }, 1);
                                // stop normal event by returning false
                                return false;
                            }
                            // cancel selection
                        }).mousedown(function () {
                        if (config.cancelSelection) {
                            this.onselectstart = function () {
                                return false
                            };
                            return false;
                        }
                    });
                    // apply easy methods that trigger binded events
                    $this.bind("update", function () {
                        var me = this;
                        setTimeout(function () {
                            // rebuild parsers.
                            me.config.parsers = buildParserCache(
                                me, $headers);
                            // rebuild the cache map
                            cache = buildCache(me);
                        }, 1);
                    }).bind("updateCell", function (e, cell) {
                        var config = this.config;
                        // get position from the dom.
                        var pos = [(cell.parentNode.rowIndex - 1), cell.cellIndex];
                        // update cache
                        cache.normalized[pos[0]][pos[1]] = config.parsers[pos[1]].format(
                            getElementText(config, cell), cell);
                    }).bind("sorton", function (e, list) {
                        $(this).trigger("sortStart");
                        config.sortList = list;
                        // update and store the sortlist
                        var sortList = config.sortList;
                        // update header count index
                        updateHeaderSortCount(this, sortList);
                        // set css for headers
                        setHeadersCss(this, $headers, sortList, sortCSS);
                        // sort the table and append it to the dom
                        appendToTable(this, multisort(this, sortList, cache));
                    }).bind("appendCache", function () {
                        appendToTable(this, cache);
                    }).bind("applyWidgetId", function (e, id) {
                        getWidgetById(id).format(this);
                    }).bind("applyWidgets", function () {
                        // apply widgets
                        applyWidget(this);
                    });
                    if ($.metadata && ($(this).metadata() && $(this).metadata().sortlist)) {
                        config.sortList = $(this).metadata().sortlist;
                    }
                    // if user has supplied a sort list to constructor.
                    if (config.sortList.length > 0) {
                        $this.trigger("sorton", [config.sortList]);
                    }
                    // apply widgets
                    applyWidget(this);
                });
            };
            this.addParser = function (parser) {
                var l = parsers.length,
                    a = true;
                for (var i = 0; i < l; i++) {
                    if (parsers[i].id.toLowerCase() == parser.id.toLowerCase()) {
                        a = false;
                    }
                }
                if (a) {
                    parsers.push(parser);
                }

            };
            this.addWidget = function (widget) {
                widgets.push(widget);
            };
            this.formatFloat = function (s) {
                var i = parseFloat(s);
                return (isNaN(i)) ? 0 : i;
            };
            this.formatInt = function (s) {
                var i = parseInt(s);
                return (isNaN(i)) ? 0 : i;
            };
            this.isDigit = function (s, config) {
                // replace all an wanted chars and match.
                return /^[-+]?\d*$/.test($.trim(s.replace(/[,.']/g, '')));
            };
            this.clearTableBody = function (table) {
                $('tbody', table).empty();
            };
        }
    });

    // extend plugin scope
    $.fn.extend({
        tablesorter: $.tablesorter.construct
    });

    // make shortcut
    var ts = $.tablesorter;

    // add default_view parsers
    ts.addParser({
        id: "text",
        is: function (s) {
            return true;
        }, format: function (s) {
            return $.trim(s.toLocaleLowerCase());
        }, type: "text"
    });

    ts.addParser({
        id: "digit",
        is: function (s, table) {
            var c = table.config;
            return $.tablesorter.isDigit(s, c);
        }, format: function (s) {
            return $.tablesorter.formatFloat(s);
        }, type: "numeric"
    });

    ts.addParser({
        id: "currency",
        is: function (s) {
            return /^[£$€?.]/.test(s);
        }, format: function (s) {
            return $.tablesorter.formatFloat(s.replace(new RegExp(/[£$€]/g), ""));
        }, type: "numeric"
    });

    ts.addParser({
        id: "ipAddress",
        is: function (s) {
            return /^\d{2,3}[\.]\d{2,3}[\.]\d{2,3}[\.]\d{2,3}$/.test(s);
        }, format: function (s) {
            var a = s.split("."),
                r = "",
                l = a.length;
            for (var i = 0; i < l; i++) {
                var item = a[i];
                if (item.length == 2) {
                    r += "0" + item;
                } else {
                    r += item;
                }
            }
            return $.tablesorter.formatFloat(r);
        }, type: "numeric"
    });

    ts.addParser({
        id: "url",
        is: function (s) {
            return /^(https?|ftp|file):\/\/$/.test(s);
        }, format: function (s) {
            return jQuery.trim(s.replace(new RegExp(/(https?|ftp|file):\/\//), ''));
        }, type: "text"
    });

    ts.addParser({
        id: "isoDate",
        is: function (s) {
            return /^\d{4}[\/-]\d{1,2}[\/-]\d{1,2}$/.test(s);
        }, format: function (s) {
            return $.tablesorter.formatFloat((s != "") ? new Date(s.replace(
                new RegExp(/-/g), "/")).getTime() : "0");
        }, type: "numeric"
    });

    ts.addParser({
        id: "percent",
        is: function (s) {
            return /\%$/.test($.trim(s));
        }, format: function (s) {
            return $.tablesorter.formatFloat(s.replace(new RegExp(/%/g), ""));
        }, type: "numeric"
    });

    ts.addParser({
        id: "usLongDate",
        is: function (s) {
            return s.match(new RegExp(/^[A-Za-z]{3,10}\.? [0-9]{1,2}, ([0-9]{4}|'?[0-9]{2}) (([0-2]?[0-9]:[0-5][0-9])|([0-1]?[0-9]:[0-5][0-9]\s(AM|PM)))$/));
        }, format: function (s) {
            return $.tablesorter.formatFloat(new Date(s).getTime());
        }, type: "numeric"
    });

    ts.addParser({
        id: "shortDate",
        is: function (s) {
            return /\d{1,2}[\/\-]\d{1,2}[\/\-]\d{2,4}/.test(s);
        }, format: function (s, table) {
            var c = table.config;
            s = s.replace(/\-/g, "/");
            if (c.dateFormat == "us") {
                // reformat the string in ISO format
                s = s.replace(/(\d{1,2})[\/\-](\d{1,2})[\/\-](\d{4})/, "$3/$1/$2");
            } else if (c.dateFormat == "uk") {
                // reformat the string in ISO format
                s = s.replace(/(\d{1,2})[\/\-](\d{1,2})[\/\-](\d{4})/, "$3/$2/$1");
            } else if (c.dateFormat == "dd/mm/yy" || c.dateFormat == "dd-mm-yy") {
                s = s.replace(/(\d{1,2})[\/\-](\d{1,2})[\/\-](\d{2})/, "$1/$2/$3");
            }
            return $.tablesorter.formatFloat(new Date(s).getTime());
        }, type: "numeric"
    });
    ts.addParser({
        id: "time",
        is: function (s) {
            return /^(([0-2]?[0-9]:[0-5][0-9])|([0-1]?[0-9]:[0-5][0-9]\s(am|pm)))$/.test(s);
        }, format: function (s) {
            return $.tablesorter.formatFloat(new Date("2000/01/01 " + s).getTime());
        }, type: "numeric"
    });
    ts.addParser({
        id: "metadata",
        is: function (s) {
            return false;
        }, format: function (s, table, cell) {
            var c = table.config,
                p = (!c.parserMetadataName) ? 'sortValue' : c.parserMetadataName;
            return $(cell).metadata()[p];
        }, type: "numeric"
    });
    // add default_view widgets
    ts.addWidget({
        id: "zebra",
        format: function (table) {
            if (table.config.debug) {
                var time = new Date();
            }
            var $tr, row = -1,
                odd;
            // loop through the visible rows
            $("tr:visible", table.tBodies[0]).each(function (i) {
                $tr = $(this);
                // style children rows the same way the parent
                // row was styled
                if (!$tr.hasClass(table.config.cssChildRow)) row++;
                odd = (row % 2 == 0);
                $tr.removeClass(
                    table.config.widgetZebra.css[odd ? 0 : 1]).addClass(
                    table.config.widgetZebra.css[odd ? 1 : 0])
            });
            if (table.config.debug) {
                $.tablesorter.benchmark("Applying Zebra widget", time);
            }
        }
    });

    // add customized widgets
    ts.addWidget({
        id: "renderCanvas",
        format: function (table) {
            d3.select(table.parentNode).selectAll('.spliceDiv').each(function () {
                /**
                 * D3 - SpliceGraph
                 * */

                var genes_obj = JSON.parse($(this)[0].getAttribute('data-exon-list').replace(/'/g, '"'));


                var exons_obj = genes_obj.exons;
                var junctions_obj = genes_obj.junctions;


                var orig_objs = {'exons': add_keys(clone(exons_obj)), 'junc': clone(junctions_obj)};

                var exons_mapped = map_exon_list(exons_obj, junctions_obj);
                exons_mapped = add_keys(exons_mapped);

                var gene_obj_cpy = {
                    'orig': orig_objs,
                    'mapped': [exons_mapped, junctions_obj],
                    'strand': genes_obj.strand
                };

                /** Render initial splice graph */
                var chart = spliceGraphD3().orig_objs(orig_objs);


                var spliceg = d3.select("#" + this.id)
                    .datum([exons_mapped, junctions_obj, genes_obj.strand])
                    .call(chart);

                gene_objs.push(gene_obj_cpy);
                if (!$(this).hasClass('exp2')) {
                    gene_obj_list.push(genes_obj);
                }

                // toggle norm flag on read counts button and redraw splicegraph
                d3.select(this.parentNode.parentNode).select('.readCounts').on('click', function () {
                    var spliceDivs = this.parentNode.parentNode.parentNode.querySelectorAll('.spliceDiv');
                    d3.selectAll(spliceDivs).call(chart);
                });

                d3.select(this.parentNode.parentNode).select('.toogleScale').on('click', function () {
                    var toogleScale = this;
                    $(this).toggleClass('scaled');

                    d3.select(this.parentNode.parentNode).selectAll('.spliceDiv').each(function () {
                        var index_gene = parseInt(this.id.split('_')[1]);

                        if (d3.select(this).classed('delta')) {
                            index_gene *= 2;
                            if (d3.select(this).classed('exp2')) {
                                index_gene++
                            }
                        }

                        var gene_datum = geneDatum(toogleScale, gene_objs[index_gene]);
                        d3.select(this).datum(gene_datum).call(chart);
                    });

                });


                d3.select(this.parentNode.parentNode).select('.zoomInSplice').on('click', function () {
                    chart.width(chart.width() + 600);
                    chart.height(chart.height() + 100);
                    d3.select(this.parentNode.parentNode).selectAll('.spliceDiv').call(chart);
                });

                d3.select(this.parentNode.parentNode).select('.zoomOutSplice').on('click', function () {
                    chart.width(chart.width() - 600);
                    chart.height(chart.height() - 100);
                    d3.select(this.parentNode.parentNode).selectAll('.spliceDiv').call(chart);
                });

                d3.select(this.parentNode.parentNode).select('.zoomResetSplice').on('click', function () {
                    chart.width(1000);
                    chart.height(160);
                    d3.select(this.parentNode.parentNode).selectAll('.spliceDiv').call(chart);
                });

                d3.select(this.parentNode.parentNode).selectAll('.weighted').on('change', function () {
                    var highlight = $(this).closest('form').find('.highlight').get(0);
                    var geneContainer = $(this).closest('.gene-container').get(0);
                    if (this.checked && !highlight.checked)
                        highlight.checked = true;
                    d3.select(geneContainer).selectAll('.spliceDiv').call(chart);
                });

                d3.select(this.parentNode.parentNode).selectAll('.highlight').on('change', function () {
                    var geneContainer = $(this).closest('.gene-container').get(0);
                    d3.select(geneContainer).selectAll('.spliceDiv').call(chart);
                });

                /**
                 * Splice Graph selector
                 * */
                var spliceGraphSelector = d3.select(this.parentNode.parentNode).selectAll('.spliceGraphSelector');
                spliceGraphSelector.on('change', function () {
                    var parent = this.parentNode.parentNode.parentNode;
                    var toogleScale = parent.querySelector('.toogleScale');
                    var spliceDiv;

                    if (d3.select(this).classed('exp2'))
                        spliceDiv = parent.querySelector('.spliceDiv.exp2');
                    else if (d3.select(this).classed('exp1'))
                        spliceDiv = parent.querySelector('.spliceDiv.exp1');

                    var index_gene = 2 * spliceDiv.id.split("_")[1];
                    if ($(parent).hasClass('exp2'))
                        index_gene++;

                    var genes_obj = JSON.parse(this.value.replace(/'/g, '"'));

                    var exons_obj = genes_obj.exons;
                    var junctions_obj = genes_obj.junctions;

                    var orig_objs = {'exons': add_keys(clone(exons_obj)), 'junc': clone(junctions_obj)};

                    var exons_mapped = map_exon_list(exons_obj, junctions_obj);
                    exons_mapped = add_keys(exons_mapped);

                    gene_objs[index_gene] = {
                        'orig': orig_objs,
                        'mapped': [exons_mapped, junctions_obj],
                        'strand': genes_obj.strand
                    };

                    var gene_datum = geneDatum(toogleScale, gene_objs[index_gene]);
                    d3.select(spliceDiv).datum(gene_datum).call(chart);

                });

                function geneDatum(toogleScale, gene) {
                    if (d3.select(toogleScale).classed('scaled'))
                        return [gene.mapped[0], gene.mapped[1], gene.strand];
                    else
                        return [gene.orig.exons, gene.orig.junc, gene.strand];
                }
            });

            $('.lsvLegend', table).each(function () {
                splicegraph().renderLsvSpliceGraph(this, gene_obj_list[this.id]);
            });

            /**
             * Single LSV visualization
             */
            $('.lsvSingleCompactPercentiles', table).each(function () {
                drawLSVCompactStackBars($(this)[0], 1);

                $(this).on("click", function (e) {

                    $(this).toggle("show");
                    var svg_children = $(this).parent().children("svg");
                    if (svg_children.length) {
                        $(svg_children[0]).toggle();
                        return;
                    }

                    // NOTE: lsv_data is an array to support groups
                    var lsv_list = JSON.parse(
                        $(this)[0]
                            .getAttribute("data-lsv")
                            .replace(/'/g, "\"")
                    );


                    var sampled_bins = translate_lsv_bins(lsv_list.bins, 1000);

                    var svg = renderViolin($(this).parent()[0].id, sampled_bins, table.id, {
                        'delta': 0,
                        'num_bins': lsv_list.bins[0].length
                    });

                    $(svg).on("click", function (e) {
                        $(this).toggle("show");
                        var lsvCompact = $(this).parent().children('.lsvSingleCompactPercentiles');
                        if (lsvCompact.length) {
                            $(lsvCompact[0]).toggle();
                        }
                    });

                });

            });


            /**
             * Delta PSI LSV visualization
             */

            $('.lsvDeltaCompact', table).each(function () {
                var lsv = JSON.parse($(this)[0].getAttribute("data-lsv").replace(/'/g, '"'));
                var threshold = $(this)[0].getAttribute("data-threshold");

                var svg_children = $(this).children(".excl-incl-rect"),
                    svgExclInclRect;

                if (svg_children.length) {
                    svgExclInclRect = svg_children[0];
                } else {
                    svgExclInclRect = drawDeltaLSVCompactSVG(this.id, lsv, threshold)[0];
                }

                $(svgExclInclRect).on("click", function (e) {
                    e.preventDefault();
                    $(this).toggle("show");


                    // For now: if LSV has 2 ways, show zoomable barchart, otherwise violin boxplots
                    if (1) { //(lsv.bins.length > 2){
                        var svgViolin;

                        var svg_children = $(this).parent().children("violin-boxplot");

                        if (svg_children.length) {
                            svgViolin = svg_children[0];
                            $(svgViolin).toggle();
                        } else {
                            var sampled_bins = translate_delta_lsv_bins(lsv.bins, 1000);
                            svgViolin = renderViolin($(this).parent()[0].id, sampled_bins, table.id, {
                                'delta': 1,
                                'num_bins': lsv.bins[0].length
                            })[0];
                        }
                        $(svgViolin).on("click", function (e) {
                            e.preventDefault();
                            $(this).toggle("show");
                            var lsvCompact = $(this).parent().children(".excl-incl-rect");
                            if (lsvCompact.length) {
                                $(lsvCompact[0]).toggle();
                            }
                        });
                    } else {
                        var parentTd = $(this).parent()[0];
                        var canvasChildren = $(parentTd).children('.extendedDeltaPsi');
                        var canvasBarchart,
                            canvasSettings;
                        if (canvasChildren.length) {
                            canvasBarchart = canvasChildren[0];
                            canvasSettings = initLargeCanvasSettings(lsv.bins[0].length, canvasBarchart);
                            $(canvasBarchart).toggle();
                        } else {
                            canvasBarchart = $('<canvas/>', {
                                'id': "barchart_" + $(this).closest('td')[0].id,
                                'class': 'extendedDeltaPsi tooltip',
                                'Title': 'Mousewheel up/down to zoom in/out'
                            })[0];
                            canvasBarchart.width = 419;
                            canvasBarchart.height = 400;
                            canvasBarchart.setAttribute('data-lsv', JSON.stringify(lsv));
                            canvasBarchart.setAttribute('data-threshold', threshold);

                            canvasSettings = initLargeCanvasSettings(lsv.bins[0].length, canvasBarchart);
                            $(parentTd).append(canvasBarchart);
                            initExpandedDeltaCanvas(canvasBarchart, canvasSettings);

                            $(canvasBarchart).on("click", function (e) {
                                e.preventDefault();
                                $(this).toggle("show");
                                var lsvCompact = $(this).parent().children(".excl-incl-rect");
                                if (lsvCompact.length) {
                                    $(lsvCompact[0]).toggle();
                                }
                            });

                            $(canvasBarchart).on('mousewheel', function (event) {
                                event.preventDefault();
                                if (event.deltaY > 0) {
                                    drawExpDeltaWithCanvasId($(this)[0].id, 1, canvasSettings);
                                } else if (event.deltaY < 0) {
                                    drawExpDeltaWithCanvasId($(this)[0].id, -1, canvasSettings);
                                }
                            });

                            $('.tooltip').tooltipster({
                                theme: 'tooltipster-shadow'
                            });
                        }
                    }
                });
            });

            $(".violin-boxplot", table).on("click", function (e) {
                e.preventDefault();
                $(this).toggle("show");
                var lsvCompact = $(this).parent().children(".excl-incl-rect");
                if (lsvCompact.length) {
                    $(lsvCompact[0]).toggle();
                }
            });

            var canvasBarcharts = $('.extendedDeltaPsi', table);
            $(canvasBarcharts).on("click", function (e) {
                e.preventDefault();
                $(this).toggle("show");
                var lsvCompact = $(this).parent().children(".excl-incl-rect");
                if (lsvCompact.length) {
                    $(lsvCompact[0]).toggle();
                }
            });

            $(canvasBarcharts).on('mousewheel', function (event) {
                event.preventDefault();

                var lsv = JSON.parse($(this)[0].getAttribute('data-lsv'));
                var canvasSettings = initLargeCanvasSettings(lsv.bins[0].length, this);

                if (event.deltaY > 0) {
                    drawExpDeltaWithCanvasId($(this)[0].id, 1, canvasSettings);
                } else if (event.deltaY < 0) {
                    drawExpDeltaWithCanvasId($(this)[0].id, -1, canvasSettings);
                }
            });

            $('.resetZoom').on("click", function () {
                var canvas = $(this).parent().children('.extendedDeltaPsi')[0];
                var lsv = JSON.parse(canvas.getAttribute('data-lsv'));
                var canvasSettings = initLargeCanvasSettings(lsv.bins[0].length, canvas);
                drawExpDeltaWithCanvasId(canvas.id, 0, canvasSettings);
            });

            var tooltip = $('.tooltip');
            if (tooltip.length) {
                tooltip.tooltipster({
                    theme: 'tooltipster-shadow'
                });
            }

        }
    });
})(jQuery);
