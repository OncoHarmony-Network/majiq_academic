from flask import request


class DataTables:
    def __init__(self, records, record_fn=None):
        self._form = request.form
        self._records = records
        self._records_fn = record_fn

        # record values
        self.records_len = len(records)
        self.draw = self._form['draw']
        self.start = int(self._form['start'])
        self.length = int(self._form['length'])
        self.column_sort = int(self._form['order[0][column]'])
        self.sort_direction = self._form['order[0][dir]']
        self.search_value = self._form['search[value]'].lower()

        # Sort records
        self._records.sort(key=lambda d: d[self.column_sort], reverse=self.sort_direction == 'desc')

        # Filter records
        self.filtered_sorted_records = list(
            filter(lambda x: any(self.search_value in y.lower() for y in x), self._records))

        # Get length of all filtered records
        self.filtered_len = len(self.filtered_sorted_records)

        # Slicing records to fit current table view
        self.filtered_sorted_records = self.filtered_sorted_records[self.start:self.start + self.length]

        # Running records function if there is one
        if self._records_fn is not None:
            self._records_fn(self.filtered_sorted_records)

    def __iter__(self):
        yield 'data', self.filtered_sorted_records
        yield 'draw', self.draw
        yield 'recordsTotal', self.records_len
        yield 'recordsFiltered', self.filtered_len
