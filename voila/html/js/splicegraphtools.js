class SpliceGraphTools {
    constructor(sgs) {
        this.sgs = sgs;
        this.init();
        this.splice_graph_onload();
    }

    init() {
        // menu drop downs
        document.querySelector('.splice-graph-tools.tools-menu-btn').onclick = (event) => {
            event.preventDefault();
            d3.selectAll('.lsv-tools.tools-menu')
                .classed('hide-tools-menu', true);
            document.querySelector('.splice-graph-tools.tools-menu').classList.toggle('hide-tools-menu');
        };

        document.querySelector('.lsv-tools.tools-menu-btn').onclick = (event) => {
            event.preventDefault();
            d3.selectAll('.splice-graph-tools.tools-menu')
                .classed('hide-tools-menu', true);
            document.querySelector('.lsv-tools.tools-menu').classList.toggle('hide-tools-menu');
        };

        // header info
        this.sgs.gene.then(gene => {
            document.querySelector('.gene-header .gene-name').textContent = `Gene name: ${gene.name}; ${ gene.chromosome }:${ gene.strand }:${ gene.start }-${ gene.end };`;
            document.querySelector('.gene-header .gene-id').textContent = `Gene ID: ${urlParams.get('gene_id')};`;
            document.querySelector('.ucsc-gene').setAttribute('href', `http://genome.ucsc.edu/cgi-bin/hgTracks?db=${gene.genome}&position=${ gene.chromosome }:${ gene.start }-${ gene.end }`)
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

                s.text(d => d);

                s.enter()
                    .append('option')
                    .text(d => d);

                s.exit().remove();

            });
        };

        // force change event to populate experiments on initial page load
        SpliceGraphTools._populate_sg_form();

        // submit event for splice graph selector
        document.querySelector('.splice-graph-form').onsubmit = (event) => {
            event.preventDefault();
            const f = event.target;
            const g = f.querySelector('#groups').value;
            const e = f.querySelector('#experiments').value;
            if (g && e) {
                this.sgs.create(g, e);
                SpliceGraphTools._populate_sg_form()
            }
        };

        // toggle splice graph scale
        document.querySelector('.toggle-scale').onclick = () => {
            document.querySelector('.splice-graph-container').classList.toggle('default-view');
            this.sgs.update(250);
        };

        // zoom in on splice graph
        document.querySelector('.zoom-in').onclick = () => {
            let zoom = parseInt(document.querySelector('.splice-graph-container').dataset.zoom);
            zoom++;
            document.querySelector('.splice-graph-container').dataset.zoom = zoom;
            this.sgs.update(250);
        };

        // zoom out on splice graph
        document.querySelector('.zoom-out').onclick = () => {
            let zoom = parseInt(document.querySelector('.splice-graph-container').dataset.zoom);
            if (zoom !== 1) {
                zoom--;
                document.querySelector('.splice-graph-container').dataset.zoom = zoom;
                this.sgs.update(250);
            }
        };

        // reset zoom for splice graph
        document.querySelector('.zoom-reset').onclick = () => {
            document.querySelector('.splice-graph-container').dataset.zoom = 1;
            this.sgs.update(250);
        };

        // activate/deactivate junction reads filter
        document.querySelector('#junction-reads-filter').onchange = (event) => {
            document.querySelectorAll('#reads-greater-than, #reads-less-than').forEach(el => el.disabled = !el.disabled);
            if (event.target.checked) {
                document.querySelectorAll('#reads-greater-than, #reads-less-than').forEach(e => e.dispatchEvent(new Event('input')));
            } else {
                d3.selectAll('.junction-grp')
                    .classed('reads-filter', false)
            }
        };

        // adjust greater than and less than fields in junction filter
        const junctions_filter = () => {
            const gt = parseInt(document.querySelector('#reads-greater-than').value);
            const lt = parseInt(document.querySelector('#reads-less-than').value);
            d3.selectAll('.junction-grp')
                .classed('reads-filter', (d, i, a) => {
                    let r = parseInt(a[i].querySelector('.junction-reads').textContent);
                    r = isNaN(r) ? 0 : r;
                    return (!isNaN(gt) && !isNaN(lt) && r <= gt || r >= lt) || (!isNaN(gt) && r <= gt) || (!isNaN(lt) && r >= lt);
                })
        };
        document.querySelector('#reads-greater-than').oninput = junctions_filter;
        document.querySelector('#reads-less-than').oninput = junctions_filter;


    }

    /**
     * Add events to splice graph after they're added to page.
     */
    splice_graph_onload() {
        new MutationObserver(mutation_list => {

            const added_nodes = Array.from(mutation_list)
                .reduce((acc, curr) => acc.concat(Array.from(curr.addedNodes)), []);

            // highlight junctions and intron retentions when you mouse over them
            added_nodes
                .filter(el => el.classList && (el.classList.contains('junction-grp') || el.classList.contains('intron-retention-grp')))
                .forEach(el => {
                    const datum = d3.select(el).datum();
                    el.onmouseover = () => {
                        const coords = document.querySelector('.coordinates');
                        if (!coords.classList.contains('select')) {
                            coords.innerHTML = `Coordinates: ${datum.start}-${datum.end}; Length: ${datum.end - datum.start}`;

                            el.classList.add('mouseover');
                            document.querySelectorAll('.junction-grp, .intron-retention-grp').forEach(el => {
                                d3.select(el)
                                    .classed('mouseover-filter', d => d.start !== datum.start || d.end !== datum.end)
                            });
                        }
                    };

                    el.onmouseout = () => {
                        const coords = document.querySelector('.coordinates');
                        if (!coords.classList.contains('select'))
                            coords.innerHTML = null;

                        el.classList.remove('mouseover');
                        document.querySelectorAll('.junction-grp, .intron-retention-grp').forEach(el => el.classList.remove('mouseover-filter'));
                    };

                    el.onclick = () => {
                        const click_new = !el.classList.contains('select');

                        document.querySelectorAll('.select-filter, .select').forEach(x => {
                            x.classList.remove('select-filter');
                            x.classList.remove('select')
                        });

                        if (click_new) {
                            el.dispatchEvent(new Event('mouseover'));
                            document.querySelectorAll('.mouseover-filter').forEach(el => el.classList.add('select-filter'));
                            el.classList.add('select');
                            document.querySelector('.coordinates').classList.add('select');
                            this.copy_select_ucsc(el)
                        } else {
                            el.dispatchEvent(new Event('mouseout'));
                        }
                    }
                });

            // add click event to remove icon
            added_nodes
                .filter(el => el.classList && el.classList.contains('splice-graph-remove'))
                .forEach(el => {
                    el.onclick = (event) => {
                        const s = event.target.closest('.splice-graph');
                        this.sgs.remove_localstorage(s);
                        s.remove();
                        SpliceGraphTools._populate_sg_form()
                    }
                })

        })
            .observe(document.querySelector('.splice-graph-container'), {
                childList: true,
                subtree: true
            });
    }

    copy_select_ucsc(el) {
        this.sgs.gene.then(gene => {
            const d = d3.select(el).datum();
            const textArea = document.createElement("textarea");
            textArea.value = `${gene.chromosome}:${d.start}-${d.end}`;
            document.body.appendChild(textArea);
            textArea.focus();
            textArea.select();
            document.execCommand('copy');
            document.body.removeChild(textArea);
        })
    }

    static _populate_sg_form() {
        document.querySelector('#groups').dispatchEvent(new Event('change'));
    }

}

