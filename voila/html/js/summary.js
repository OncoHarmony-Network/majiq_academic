window.addEventListener('load', () => {
    const urlParams = new URLSearchParams(document.location.search);
    const gene_id = urlParams.get('gene_id');

    const gene_ids = Array.from(new Set(lsvs_arr.map(x => x.gene_id)));
    gene_ids.sort((a, b) => a.localeCompare(b));

    const gene_idx = gene_ids.indexOf(gene_id);

    const first_gene = () => {
        window.location.href = `summary.html?gene_id=${gene_ids[0]}`;
    };

    const last_gene = () => {
        window.location.href = `summary.html?gene_id=${gene_ids[gene_ids.length - 1]}`;
    };

    if (!gene_id)
        first_gene();

    document.querySelector('.prev-gene').addEventListener('click', () => {
        const prev_gene_id = gene_ids[gene_idx - 1];
        if (prev_gene_id)
            window.location.href = `summary.html?gene_id=${prev_gene_id}`;
        else
            last_gene()
    });

    document.querySelector('.next-gene').addEventListener('click', () => {
        const next_gene_id = gene_ids[gene_idx + 1];
        if (next_gene_id)
            window.location.href = `summary.html?gene_id=${next_gene_id}`;
        else
            first_gene()
    });
});

