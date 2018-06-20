window.addEventListener('load', () => {
    // get query string
    const urlParams = new URLSearchParams(document.location.search);

    // splice_graph_update sg container dataset
    const sg_cont = document.querySelector('.splice-graph-container');
    sg_cont.dataset.geneId = urlParams.get('gene_id');
    sg_cont.dataset.zoom = '1';

    // menu drop downs
    document.querySelector('.tools.pure-menu-link').onclick = (event) => {
        event.preventDefault();
        document.querySelector('.tools.submenu').classList.toggle('hide-submenu');
    };

    // header info
    db_gene.createIndex({
        index: {fields: ['_id']}
    }).then(() => {
        return db_gene.find({
            selector: {
                _id: urlParams.get('gene_id')
            },
            fields: ['name', 'chromosome', 'strand', 'start', 'end']
        })
    }).then(results => {
        const gene = results.docs[0];
        document.querySelector('.gene-header .gene-name').textContent = `Gene name: ${gene.name}; ${ gene.chromosome }:${ gene.strand }:${ gene.start }-${ gene.end };`;
        document.querySelector('.gene-header .gene-id').textContent = `Gene ID: ${urlParams.get('gene_id')};`;
    });

    // populate splice graph selector groups
    db_gene.createIndex({
        index: {fields: ['_id']}
    }).then(() => {
        return db_gene.find({
            selector: {
                _id: 'metadata'
            },
        })
    }).then(results => {
        const meta = results.docs[0];
        d3.select('.groups select')
            .selectAll('option')
            .data(meta.group_names)
            .enter()
            .append('option')
            .text(d => {
                return d
            });
    });

    // populate splice graph experiments when group is changed
    document.querySelector('#groups').onchange = (event) => {
        db_gene.createIndex({
            index: {fields: ['_id']}
        }).then(() => {
            return db_gene.find({
                selector: {
                    _id: 'metadata'
                },
            })
        }).then(results => {
            const sgs = JSON.parse(localStorage.getItem('splice_graphs'));
            let exps = results.docs[0].experiment_names[event.target.selectedIndex];
            if (sgs)
                exps = exps.filter(exp => {
                    return !sgs.some(sg => sg[1] === exp)
                });

            const s = d3.select('.experiments select')
                .selectAll('option')
                .data(exps);

            s.text(d => {
                return d
            });

            s.enter()
                .append('option')
                .text(d => {
                    return d
                });

            s.exit().remove();

        });
    };

    // force change event to populate experiments on initial page load
    document.querySelector('#groups').dispatchEvent(new Event('change'));

    // submit event for splice graph selector
    document.querySelector('.splice-graph-form').onsubmit = (event) => {
        event.preventDefault();
        const f = event.target;
        const g = f.querySelector('#groups').value;
        const e = f.querySelector('#experiments').value;
        if (g && e)
            sgs.create(g, e)
    };

    document.querySelector('.toggle-scale').onclick = () => {
        document.querySelector('.splice-graph-container').classList.toggle('default-view');
        sgs.update();
    }
});
