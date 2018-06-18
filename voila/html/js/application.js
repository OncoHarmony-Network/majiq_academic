const first_gene = () => {
    db_gene.createIndex({
        index: {fields: ['_id']}
    }).then(() => {
        return db_gene.find({
            selector: {_id: {$ne: 'metadata'}},
            limit: 1,
            sort: ['_id'],
            fields: ['_id'],
        })
    }).then(results => {
        window.location.href = `summary.html?gene_id=${results.docs[0]._id}`;
    })
};

const last_gene = () => {
    db_gene.createIndex({
        index: {fields: ['_id']}
    }).then(() => {
        return db_gene.find({
            selector: {_id: {$ne: 'metadata'}},
            limit: 1,
            sort: [{'_id': 'desc'}],
            fields: ['_id'],
        })
    }).then(results => {
        window.location.href = `summary.html?gene_id=${results.docs[0]._id}`;
    })
};

window.addEventListener('load', () => {
    const gene_id = urlParams.get('gene_id');
    if (!gene_id)
        first_gene();

    document.querySelector('.prev-gene').addEventListener('click', () => {
        db_gene.createIndex({
            index: {fields: ['_id']}
        }).then(() => {
            return db_gene.find({
                selector: {_id: {$lt: urlParams.get('gene_id')}},
                limit: 2,
                sort: [{_id: 'desc'}],
                fields: ['_id'],
            })
        }).then(results => {
            if (results.docs.length)
                window.location.href = `summary.html?gene_id=${results.docs[0]._id}`;
            else
                last_gene()
        })
    });
    document.querySelector('.next-gene').addEventListener('click', () => {
        db_gene.createIndex({
            index: {fields: ['_id']}
        }).then(() => {
            return db_gene.find({
                selector: {_id: {$gt: urlParams.get('gene_id')}},
                limit: 2,
                sort: ['_id'],
                fields: ['_id'],
            })
        }).then(results => {
            if (results.docs.length)
                window.location.href = `summary.html?gene_id=${results.docs[0]._id}`;
            else
                first_gene()
        })
    });
});

