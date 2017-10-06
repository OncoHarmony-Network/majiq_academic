import h5py
import numpy as np


class PropertyDoesNotExist(Exception):
    def __init__(self, k, sgt):
        msg = '"{0}" is not a property of {1}. Expecting: {2}'.format(k, sgt.__class__.__name__, ', '.join(sgt._props))
        super(PropertyDoesNotExist, self).__init__(msg)


class SpliceGraphType:
    def __init__(self, hdf5_grp):
        self._hdf5_grp = hdf5_grp
        self._props = ()
        self._process = None

    def __iter__(self):
        for k, v in self._hdf5_grp.attrs.items():
            yield k, v

    def __setattr__(self, name, value):
        if name[0] == '_':
            return super().__setattr__(name, value)
        return self.parse_attrs(**{name: value})

    def __getattr__(self, item):
        if item in self._props:
            try:
                value = self._hdf5_grp.attrs[item]
                if isinstance(value, np.bytes_):
                    value = value.decode('utf-8')
                return value
            except KeyError:
                return None

        raise PropertyDoesNotExist(item, self)

    @property
    def id(self):
        return self._hdf5_grp.name.split('/')[-1]

    def parse_attrs(self, **kwargs):
        a = self._hdf5_grp.attrs
        for k, v in kwargs.items():
            if k in self._props:

                if self._process and k in self._process:
                    v = self._process[k](v)

                # filter out values that are empty list/strings and False booleans.
                if hasattr(v, '__iter__') and not len(v):
                    continue
                elif v is False or v is np.False_:
                    continue

                if isinstance(v, str):
                    try:
                        v = np.string_(v)
                    except UnicodeEncodeError:
                        pass
                elif isinstance(v, (bool, np.bool_)):
                    v = int(v)

                try:
                    a[k] = v
                except TypeError:
                    print(v)
                    raise

            else:
                raise PropertyDoesNotExist(k, self)


class Junction(SpliceGraphType):
    def __init__(self, hdf5_grp, **kwargs):
        super().__init__(hdf5_grp)
        self._props = {'start', 'end', 'junction_type_list', 'reads_list', 'transcripts', 'intron_retention'}
        self.parse_attrs(**kwargs)

    def update_reads(self, experiment_name, reads):
        try:
            exp_names = self._hdf5_grp.file.attrs['experiment_names']
        except KeyError:
            raise KeyError("Experiment names have not been set.")

        try:
            idx = np.where(exp_names == experiment_name)[0][0]
        except IndexError:
            raise IndexError('Experiment name could not be found.')

        if self.reads_list is None:
            reads_list = [0] * len(exp_names)
        else:
            reads_list = self.reads_list

        reads_list[idx] = reads
        self.reads_list = reads_list


class Exon(SpliceGraphType):
    def __init__(self, hdf5_grp, **kwargs):
        super().__init__(hdf5_grp)
        self._props = {'end', 'start', 'a3', 'a5', 'exon_type_list', 'coords_extra', 'intron_retention', 'lsv_type',
                       'alt_starts', 'alt_ends'}
        self.parse_attrs(**kwargs)

    @property
    def coords(self):
        return self.start, self.end


class Gene(SpliceGraphType):
    def __init__(self, hdf5_grp, **kwargs):
        super().__init__(hdf5_grp)
        self._props = {'name', 'strand', 'chromosome', 'junctions', 'exons'}
        # self._process = {'junctions': self._references_juncs, 'exons': self._references_exons}
        self._process = {'junctions': self._ids, 'exons': self._ids}
        self.parse_attrs(**kwargs)

    @staticmethod
    def _ids(vs):
        return ','.join(v.id for v in vs)
        # return np.array(tuple(v.id for v in vs), dtype=h5py.special_dtype(vlen=np.unicode))

    def _references(self, vs):
        return np.array(tuple(v._hdf5_grp.ref for v in vs), dtype=h5py.special_dtype(ref=h5py.Reference))

    def _references_juncs(self, vs):
        for j in vs:
            if not isinstance(j, Junction):
                raise Exception('Wrong!')
        return self._references(vs)

    def _references_exons(self, vs):
        for e in vs:
            if not isinstance(e, Exon):
                raise Exception('Wrong!')
        return self._references(vs)

    @property
    def start(self):
        return next(self.exons).start

    @property
    def end(self):
        return tuple(self.exons)[-1].end

    @property
    def junctions(self):
        for junc_id in self._hdf5_grp.attrs['junctions'].decode('utf-8').split(','):
            yield Junction(self._hdf5_grp.file['Junctions'][junc_id])

    @property
    def exons(self):
        for exon_id in self._hdf5_grp.attrs['exons'].decode('utf-8').split(','):
            yield Exon(self._hdf5_grp.file['Exons'][exon_id])

    def get_experiment(self, experiment_index):
        d = dict(self)
        d['exons'] = [dict(e) for e in self.exons]
        d['junctions'] = [dict(j) for j in self.junctions]
        d['start'] = d['exons'][0]['start']
        d['end'] = d['exons'][-1]['end']

        for j in d['junctions']:
            j['reads'] = j['reads_list'][experiment_index]
            del j['reads_list']
            if 'junction_type_list' in j:
                del j['junction_type_list']
            j['junction_type'] = 0
            if 'intron_retention' not in j:
                j['intron_retention'] = 0

        for e in d['exons']:
            if 'exon_type_list' in e:
                del e['exon_type_list']
            e['exon_type'] = 0
            if 'intron_retention' not in e:
                e['intron_retention'] = 0

        return d

    def combine(self, experiment_index, gene_dict=None):
        if not gene_dict:
            return self.get_experiment(experiment_index)

        # gene_dict['junctions'] = [s.combine(experiment_index, d) for s, d in
        #                           zip(self.junctions, gene_dict['junctions'])]
        #
        # gene_dict['exons'] = [s.combine(experiment_index, d) for s, d in
        #                       zip(self.exons, gene_dict['exons'])]

        for j, s in zip(gene_dict['junctions'], self.junctions):
            j['reads'] += s.reads_list[experiment_index]

        return gene_dict
