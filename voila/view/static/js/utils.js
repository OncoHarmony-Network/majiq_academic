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