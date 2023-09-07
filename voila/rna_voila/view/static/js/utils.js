function d3_parent(){ return this.parentNode; }

function isInArr(item, array){
    return array.indexOf(item) !== -1;
}

function attrDeltaInt(elem, attr, delta){
    elem.attr(attr, parseInt(elem.attr(attr)) + delta);
}

function attrDelta(elem, attr, delta){
    elem.attr(attr, elem.attr(attr) + delta);
}

function removeFromArray(item, array){
    array.splice(array.indexOf(item), 1);
}

$.fn.sortClass = function sortDivs(_class, _attr) {
    $("> ." + _class, this[0]).sort(dec_sort).appendTo(this[0]);
    function dec_sort(a, b){ return ($(b).data(_attr)) < ($(a).data(_attr)) ? 1 : -1; }
}

function dispFadeAlert(text){
    $('body').append(`<div class="tmp-alert">${text}</div>`);
    $('.tmp-alert').fadeOut(2000, function(){
        $(this).remove();
    });
}

function isEmpty(obj){
    return obj === undefined ? true : Object.keys(obj).length === 0;
}

// from https://bl.ocks.org/tophtucker/62f93a4658387bb61e4510c37e2e97cf
// set to 'sans-serif' currently
function measureText(string, fontSize = 10) {
  const widths = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0546875,0.4,0.6,0.8,0.8,1.1,0.9,0.4,0.6,0.5,0.6,0.8,0.4,0.5,0.4,0.5,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.4,0.4,0.8,0.8,0.8,0.8,1.2,0.9,0.9,0.9,0.9,0.9,0.8,1,0.9,0.4,0.7,0.9,0.8,1,0.9,1,0.9,1,0.9,0.9,0.8,0.9,0.9,1.2,0.9,0.9,0.8,0.5,0.5,0.5,0.7,0.9,0.5,0.8,0.8,0.7,0.7,0.8,0.5,0.7,0.7,0.4,0.5,0.8,0.4,1,0.7,0.8,0.8,0.7,0.6,0.7,0.5,0.7,0.7,1.1,0.7,0.7,0.7,0.6,0.4,0.6,0.8]
  const avg = 0.7342598684210524
  return string
    .split('')
    .map(c => c.charCodeAt(0) < widths.length ? widths[c.charCodeAt(0)] : avg)
    .reduce((cur, acc) => acc + cur) * fontSize
}

/**
 * Uses canvas.measureText to compute and return the width of the given text of given font in pixels.
 *
 * @param {String} text The text to be rendered.
 * @param {String} font The css font descriptor that text is to be rendered with (e.g. "bold 14px verdana").
 *
 * @see https://stackoverflow.com/questions/118241/calculate-text-width-with-javascript/21015393#21015393
 * Best usage: getTextWidth(text, getCanvasFontSize(myEl))
 */
function getTextWidth(text, font = getCanvasFontSize()) {
    // re-use canvas object for better performance
    const canvas = getTextWidth.canvas || (getTextWidth.canvas = document.createElement("canvas"));
    const context = canvas.getContext("2d");
    context.font = font;
    const metrics = context.measureText(text);
    return metrics.width;
}

function getCssStyle(element, prop) {
    return window.getComputedStyle(element, null).getPropertyValue(prop);
}

function getCanvasFontSize(el = document.body) {
    const fontWeight = getCssStyle(el, 'font-weight') || 'normal';
    const fontSize = getCssStyle(el, 'font-size') || '16px';
    const fontFamily = getCssStyle(el, 'font-family') || 'Times New Roman';

    return `${fontWeight} ${fontSize} ${fontFamily}`;
}

function save_data_as_file(mimetype, content, filename){
    const element = document.createElement('a');
    element.setAttribute('href', mimetype + encodeURIComponent(content));
    element.setAttribute('download', filename);

    element.style.display = 'none';
    document.body.appendChild(element);

    element.click();

    document.body.removeChild(element);
}

function save_html_file(html_content, filename){
    save_data_as_file('data:text/html;charset=utf-8,', html_content, filename);
}

function download_svg_elem(svg, filename){
    save_data_as_file('data:image/svg+xml;charset=utf-8,', svg, filename);
}