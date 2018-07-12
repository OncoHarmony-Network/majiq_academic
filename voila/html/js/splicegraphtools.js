window.addEventListener('load', () => {
    // get query string
    const urlParams = new URLSearchParams(document.location.search);

    // splice_graph_update sg container dataset
    const sg_cont = document.querySelector('.splice-graph-container');
    sg_cont.dataset.geneId = urlParams.get('gene_id');
    sg_cont.dataset.zoom = '1';

    // menu drop downs
    document.querySelector('.splice-graph-tools.pure-menu-link').onclick = (event) => {
        event.preventDefault();
        document.querySelector('.splice-graph-tools.submenu').classList.toggle('hide-submenu');
    };

    document.querySelector('.lsv-tools.pure-menu-link').onclick = (event) => {
        event.preventDefault();
        document.querySelector('.lsv-tools.submenu').classList.toggle('hide-submenu');
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
        if (!gene)
            first_gene();
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

    // toggle splice graph scale
    document.querySelector('.toggle-scale').onclick = () => {
        document.querySelector('.splice-graph-container').classList.toggle('default-view');
        sgs.update(250);
    };

    // zoom in on splice graph
    document.querySelector('.zoom-in').onclick = () => {
        let zoom = parseInt(document.querySelector('.splice-graph-container').dataset.zoom);
        zoom++;
        document.querySelector('.splice-graph-container').dataset.zoom = zoom;
        sgs.update(250);
    };

    // zoom out on splice graph
    document.querySelector('.zoom-out').onclick = () => {
        let zoom = parseInt(document.querySelector('.splice-graph-container').dataset.zoom);
        if (zoom !== 1) {
            zoom--;
            document.querySelector('.splice-graph-container').dataset.zoom = zoom;
            sgs.update(250);
        }
    };

    // reset zoom for splice graph
    document.querySelector('.zoom-reset').onclick = () => {
        document.querySelector('.splice-graph-container').dataset.zoom = 1;
        sgs.update(250);
    };

    // activate/deactivate junction reads filter
    document.querySelector('#junction-reads-filter').onchange = (event) => {
        document.querySelectorAll('#reads-greater-than, #reads-less-than').forEach(el => el.disabled = !el.disabled);
        if (event.target.checked) {
            document.querySelectorAll('#reads-greater-than, #reads-less-than').forEach(e => e.dispatchEvent(new Event('input')));
        } else {
            d3.selectAll('.junction-grp')
                .attr('opacity', null)
        }
    };

    // adjust greater than and less than fields in junction filter
    const junctions_filter = () => {
        const gt = parseInt(document.querySelector('#reads-greater-than').value);
        const lt = parseInt(document.querySelector('#reads-less-than').value);
        d3.selectAll('.junction-grp')
            .attr('opacity', (d, i, a) => {
                let r = parseInt(a[i].querySelector('.junction-reads').textContent);
                r = isNaN(r) ? 0 : r;
                return (!isNaN(gt) && !isNaN(lt) && r <= gt || r >= lt) || (!isNaN(gt) && r <= gt) || (!isNaN(lt) && r >= lt) ? .2 : null;
            })
    };
    document.querySelector('#reads-greater-than').oninput = junctions_filter;
    document.querySelector('#reads-less-than').oninput = junctions_filter;

    // // highlight junctions when you mouse over them
    // new MutationObserver(mutation_list => {
    //     mutation_list.forEach(m => {
    //         m.addedNodes.forEach(n => {
    //             if (n.classList && n.classList.contains('junction-grp')) {
    //                 n.addEventListener('mouseover', (event) => {
    //                     const datum = d3.select(event.target).datum();
    //                     d3.selectAll('.junction-grp')
    //                         .each((d, i, a) => {
    //                             const el = d3.select(a[i]);
    //                             if (datum.start !== d.start || datum.end !== d.end) {
    //                                 el.attr('opacity', 0.2)
    //                             }
    //                         })
    //                 });
    //
    //                 n.addEventListener('mouseout', () => {
    //                     d3.selectAll('.junction-grp')
    //                         .attr('opacity', null);
    //                 })
    //             }
    //         })
    //     });
    // })
    //     .observe(document.querySelector('.splice-graph-container'), {
    //         childList: true,
    //         subtree: true
    //     });
});
