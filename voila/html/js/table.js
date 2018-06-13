class Table {
    constructor(table) {
        this.table = table;
        table.onclick = (event) => {
            const col = event.target;
            if (col.tagName === 'TH') {
                const col_idx = Array.from(table.querySelectorAll('th')).indexOf(col);
                col.asc = !col.asc;
                this.sort(col_idx, col.asc);
            }
        }
    }

    comparator(col_idx, asc) {
        return (a, b) => {
            const a1 = Table.col_value(a, col_idx);
            const b1 = Table.col_value(b, col_idx);
            return asc ? a1.toString().localeCompare(b1) : b1.toString().localeCompare(a1);
        };
    }

    sort(col_idx, asc) {
        Array.from(this.table.querySelectorAll('tbody tr'))
            .sort(this.comparator(col_idx, asc))
            .forEach(tr => this.table.querySelector('tbody').appendChild(tr));
    }

    static col_value(tr, col_idx) {
        return tr.children[col_idx].textContent
    }
}