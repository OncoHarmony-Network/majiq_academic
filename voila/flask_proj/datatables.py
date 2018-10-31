from flask import request


class DataTables:
    def __init__(self, records, callback=None):
        self._form = request.form
        self._records = list(records)

        # record values
        self.records_len = len(self._records)
        self.draw = self._form['draw']
        self.start = int(self._form['start'])
        self.length = int(self._form['length'])
        self.column_sort = int(self._form['order[0][column]'])
        self.sort_direction = self._form['order[0][dir]']
        self.search_value = self._form['search[value]'].lower().strip()

        # Sort records
        self._records.sort(key=lambda x: x[self.column_sort], reverse=self.sort_direction == 'desc')

        # Filter records
        if self.search_value:
            self.filtered_sorted_records = list(filter(self._filter, self._records))
        else:
            self.filtered_sorted_records = self._records

        # Get length of all filtered records
        self.filtered_len = len(self.filtered_sorted_records)

        # Slicing records to fit current table view
        if self.length != -1:
            self.filtered_sorted_records = self.filtered_sorted_records[self.start:self.start + self.length]

        if callback:
            callback(self.filtered_sorted_records)

    def __iter__(self):
        yield 'data', self.filtered_sorted_records
        yield 'draw', self.draw
        yield 'recordsTotal', self.records_len
        yield 'recordsFiltered', self.filtered_len

    def _filter(self, vs):
        def parse(v):
            try:
                v = v.get('display', '')
            except AttributeError:
                pass

            return self.search_value in str(v).lower()

        return any(parse(v) for v in vs)
