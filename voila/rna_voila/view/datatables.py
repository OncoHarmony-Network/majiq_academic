import json

from flask import request, session

from rna_voila.exceptions import SortFunctionNotFound


class DataTables:
    def __init__(self, records, sort_columns, sort=True, slice=True):
        records = list(records)
        self._form = request.form

        self.sort_columns = sort_columns
        self.filtered_len = 0

        # record values
        self.records_len = len(records)
        self.draw = self._form.get('draw', '')
        self.start = int(self._form.get('start', 0))
        self.length = int(self._form.get('length', 0))
        self.column_sort = int(self._form.get('order[0][column]', 0))
        self.sort_direction = self._form.get('order[0][dir]', '')
        self.search_value = self._form.get('search[value]', '').lower().strip()
        self.extra_sort = {}

        self._records = records
        self._filter()

        if sort:
            self.sort()

        if slice:
            self.slice()

    def __iter__(self):
        yield 'data', self._records
        yield 'draw', self.draw
        yield 'recordsTotal', self.records_len
        try:
            yield 'recordsFiltered', self.filtered_len
        except AttributeError:
            yield 'recordsFiltered', 0

    def _search_value_filter(self, vs):
        return any(self.search_value in str(v).lower() for v in vs.values())

    def _filter(self):
        lsv_filters = list(self.filter_values('lsv_filter').keys())

        if self.search_value:
            self._records = list(r for r in self._records if self._search_value_filter(r))

        for f in lsv_filters:
            self._records = list(r for r in self._records if r[f])

    def sort(self):
        col_name = self.sort_columns[self.column_sort]

        try:
            if col_name in self.extra_sort:
                self._records.sort(
                    key=lambda x: self.extra_sort[col_name](x),
                    reverse=self.sort_direction == 'desc')
            else:
                self._records = list(self._records)
                self._records.sort(key=lambda x: x[col_name],
                                   reverse=self.sort_direction == 'desc')
        except KeyError:
            raise SortFunctionNotFound(self.column_sort)

    def extra_filter(self, filter_fn):
        self._records = list(r for r in self._records if filter_fn(r))

    def slice(self):
        self._records = list(self._records)

        # Get length of all filtered records
        self.filtered_len = len(self._records)

        # Slicing records to fit current table view
        if self.length != -1:
            self._records = self._records[self.start:self.start + self.length]

    def callback(self):
        for idx in range(len(self._records)):
            yield idx, self._records[idx], self._records

    def filter_values(self, filter_name):
        return {k.split('[')[1][:-1]: v for k, v in self._form.items() if k.startswith(filter_name)}

    def heterogen_filters(self):
        het_filters = self.filter_values('het_filter')

        dpsi_thresh = het_filters['dpsi_threshold']
        stat_thresh = het_filters['stat_threshold']
        stat_type_i = int(het_filters['stat_type'])

        try:
            dpsi_thresh = float(dpsi_thresh)
        except ValueError:
            dpsi_thresh = 0

        try:
            stat_thresh = float(stat_thresh)
        except ValueError:
            stat_thresh = 0

        #
        # for rs in self._records:
        #     print(stat_thresh, json.loads(rs['stat_threshold']))

        self.extra_filter(lambda rs: json.loads(rs['dpsi_threshold']) >= dpsi_thresh)
        self.extra_filter(lambda rs: json.loads(rs['stat_threshold'])[stat_type_i] <= stat_thresh)

        return stat_type_i

    def delta_psi_filters(self):
        dpsi_filter_values = self.filter_values('dpsi_filter')

        dpsi_thresh = dpsi_filter_values['dpsi_threshold']
        confid_thresh = dpsi_filter_values['confidence_threshold']

        try:
            dpsi_thresh = float(dpsi_thresh)
        except ValueError:
            dpsi_thresh = 0

        try:
            confid_thresh = float(confid_thresh)
        except ValueError:
            confid_thresh = 0

        confid_idx = int(dpsi_thresh * 10)

        self.extra_filter(lambda rs: any(r >= dpsi_thresh for r in json.loads(rs['dpsi_threshold'])))
        self.extra_filter(lambda rs: confid_thresh <= json.loads(rs['confidence_threshold'])[confid_idx])

    def add_sort(self, sort_col, sort_fn):
        self.extra_sort[sort_col] = sort_fn

    class Sort:
        @staticmethod
        def highlight(row):
            lsv_id = row['lsv_id'].decode('utf-8')
            try:
                return session['highlight'][lsv_id]
            except KeyError:
                return [False, False]

        @staticmethod
        def lsv_id(row):
            ref_exon = row['lsv_id'].decode('utf-8')
            ref_exon = ref_exon.split(':')
            ref_exon = ref_exon[-1]
            ref_exon = ref_exon.split('-')
            try:
                ref_exon = list(map(int, ref_exon))
            except ValueError:
                if ref_exon[0] == 'na':
                    ref_exon[1] = int(ref_exon[1])
                    ref_exon[0] = ref_exon[1] - 10
                if ref_exon[1] == 'na':
                    ref_exon[0] = int(ref_exon[0])
                    ref_exon[1] = ref_exon[0] + 10
            return ref_exon
