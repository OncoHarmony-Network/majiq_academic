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