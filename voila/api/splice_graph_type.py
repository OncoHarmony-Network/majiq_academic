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
                return self._hdf5_grp.attrs[item]
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
        self._process = {'junctions': self._references, 'exons': self._references}
        self.parse_attrs(**kwargs)

    def _references(self, vs):
        return np.array(tuple(v._hdf5_grp.ref for v in vs), dtype=h5py.special_dtype(ref=h5py.Reference))

    @property
    def junctions(self):
        return tuple(Junction(self._hdf5_grp[ref]) for ref in self._hdf5_grp.attrs['junctions'])

    @property
    def exons(self):
        return tuple(Exon(self._hdf5_grp[ref]) for ref in self._hdf5_grp.attrs['exons'])
