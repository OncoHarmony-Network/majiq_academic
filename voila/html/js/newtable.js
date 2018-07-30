class NewTable {
    constructor(table_selector, db_gene, db_lsv, opts) {
        this.table_selector = table_selector;
        this.get_data = opts.get_data;
        this.curate_data = opts.curate_data;
        this.show_data = opts.show_data;

        this.update();

        this.table.onclick = (event) => {
            const col = event.target;
            if (col.tagName === 'TH') {
                const col_idx = Array.from(this.table.querySelectorAll('th')).indexOf(col);
                col.asc = !col.asc;
                this.sort(col_idx, col.asc);
            }
        }
    }

    get table() {
        return document.querySelector(this.table_selector)
    }

    get table_body() {
        return document.querySelector(`${this.table_selector} tbody`)
    }

    static col_value(tr, col_idx) {
        const cell = tr.children[col_idx];
        return 'sortValue' in cell.dataset ? cell.dataset.sortValue : cell.textContent
    }

    update() {
        this.get_data(db_gene, db_lsv).then(this.curate_data).then(data => this.show_data(data, this.table_body));
    }

    comparator(col_idx, asc) {
        return (a, b) => {
            const a1 = NewTable.col_value(a, col_idx);
            const b1 = NewTable.col_value(b, col_idx);
            return asc ? a1.toString().localeCompare(b1) : b1.toString().localeCompare(a1);
        };
    }

    sort(col_idx, asc) {
        Array.from(this.table.querySelectorAll('tbody tr'))
            .sort(this.comparator(col_idx, asc))
            .forEach(tr => this.table.querySelector('tbody').appendChild(tr));
    }
}