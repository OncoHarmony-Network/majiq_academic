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

function download_svg_elem(svg, filename){
    const element = document.createElement('a');
    element.setAttribute('href', 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svg));
    element.setAttribute('download', filename);

    element.style.display = 'none';
    document.body.appendChild(element);

    element.click();

    document.body.removeChild(element);
}