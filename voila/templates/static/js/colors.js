var Colors = function () {
    this.BREWER_PALETTE = [
        [228, 26, 28],
        [55, 126, 184],
        [77, 175, 74],
        [152, 78, 163],
        [255, 127, 0],
//    [255,255,51],
        [166, 86, 40],
        [247, 129, 191],
        [153, 153, 153],

        [28, 126, 128],
        [155, 226, 29],
        [177, 275, 19],
        [252, 178, 8],
        [55, 227, 100],
//    [55,55,151],
        [11, 186, 140],
        [47, 229, 36],
        [253, 253, 253]
    ];

    this.colors = ['red', 'blue', 'green', 'purple', 'orange', 'brown', 'pink']
};

Colors.prototype.brewer = function (color_index) {
    return 'rgb(' + this.BREWER_PALETTE[color_index % this.BREWER_PALETTE.length].join(',') + ')'
};

Colors.prototype.brewer_het = function (color_index) {

    return function () {
        return 'rgb(' + this.BREWER_PALETTE[color_index % this.BREWER_PALETTE.length].join(',') + ')'
    };
};

Colors.prototype.color = function (color_index) {
    return this.colors[color_index % this.colors.length]
};