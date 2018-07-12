window.addEventListener('load', () => {
    document.querySelector('.lsv-filters').onchange = (event) => {
        load_lsvs()
    }
});