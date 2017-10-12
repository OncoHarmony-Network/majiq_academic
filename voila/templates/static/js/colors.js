var Colors = function (opts) {
    this.color_arr = BREWER_PALETTE
};

Colors.prototype.toArray = function () {
    return this.color_arr;
};

Colors.prototype.toRGBArray = function () {
    rgb_arr = [];
    this.color_arr.forEach(function (color, color_index) {
        rgb_arr.push('rgb(' + BREWER_PALETTE[color_index % BREWER_PALETTE.length].join(',') + ')')
    });
    return rgb_arr
};