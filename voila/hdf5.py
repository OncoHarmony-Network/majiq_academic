from os.path import join

import numpy
from h5py.h5r import RegionReference

from voila.constants import EXPERIMENT_NAMES
from voila.utils.voila_log import voila_log

VOILA_FILE_VERSION = '/voila_file_version'


class HDF5VersionException(Exception):
    def __init__(self):
        super(HDF5VersionException, self).__init__('The hdf5 file version does not match the current version of Voila.')
        voila_log().error(self.message)


class HDF5(object):
    def __init__(self):
        """
        Move data in this class to and from HDF5 files.
        """
        pass

    @staticmethod
    def copy_group(src, dst):
        """
        Copy HDF5 group data to another HDF5 group.  This can be within the same file or across two files.  Sometimes
        there are attributes on the source group, these are also copied to the destination group.  This will copy
        references to the destination, but if the destination is in a different file, the referenced data will not be
        copied to the new file.
        :param src: source HDF5 group
        :param dst: destination HDF5 group
        :return: None
        """
        for key in src:
            try:
                src.copy(key, dst)
            except ValueError:
                del dst[key]
                src.copy(key, dst)

        for key in src.attrs:
            dst.attrs[key] = src.attrs[key]

    @staticmethod
    def copy_lsv_graphic(lsv_src, lsv_dst):
        """
        Copy LsvGraphic data from one VoilaLsv object to another.
        :param lsv_src: source VoilaLsv
        :param lsv_dst: destination VoilaLsv
        :return: None
        """
        lsv_graphic = 'lsv_graphic'
        try:
            dst = lsv_dst[lsv_graphic]
        except KeyError:
            dst = lsv_dst.create_group(lsv_graphic)
        HDF5.copy_group(lsv_src[lsv_graphic], dst)

    def exclude(self):
        """
        Exclude class attributes from being processed.  This infers that these attributes will be processed by the
        subclass in their to_hdf5 and from_hdf5 methods.
        :return: list of attributes that aren't processed.
        """
        return ()

    def cls_list(self):
        """
        Identify class attributes that are a list of a specific class.  e.g [VoilaLsv(), VoilaLsv(), ... ].  The
        returned dictionary should contain the name of the attribute that points to the list of classes and how
        instantiate the classes in the from_hdf.  See VoilaInput for an example.
        :return: Dictionary of classes
        """
        return {}

    def to_hdf5(self, h, use_id=True):
        """
        Adds attributes from this class to HDF5 file.
        :param use_id:
        :param h: HDF5 file object
        :return: None
        """

        attrs_dict = self.__dict__.copy()

        for key in self.exclude():
            del attrs_dict[key]

        for cls in self.cls_list():
            for index, cls_obj in enumerate(attrs_dict[cls]):
                cls_obj.to_hdf5(h.create_group(join(cls, str(index))), use_id)
            del attrs_dict[cls]

        for key in attrs_dict:
            try:
                h.attrs[key] = attrs_dict[key]
            except TypeError as e:
                # Where the value stored is None, skip it.  This assumes that the default value for the attribute in
                # the class is sufficient and that 'None' doesn't have some extra meaning for this attribute.
                if attrs_dict[key] is None:
                    pass
                else:
                    raise type(e)('There was an issue with key "{0}" with value "{1}"'.format(key, attrs_dict[key]))

    def from_hdf5(self, h):
        """
        Populate this class with values from HDF5 file object.
        :param h: HDF5 file object
        :return: self
        """

        cls_dict = self.cls_list()

        for cls in cls_dict:
            if cls in h:
                cls_grp = h[cls]
                self.__dict__[cls] = [None] * len(cls_grp)
                for index in cls_grp:
                    new_class = cls_dict[cls]
                    self.__dict__[cls][int(index)] = new_class.easy_from_hdf5(cls_grp[index])

        for key in h.attrs:
            if key not in self.exclude():
                # H5py stores attributes as numpy objects where it can.  Numpy objects
                # cause issues with the json conversion, therefore corner cases are handled below.
                value = h.attrs[key]
                if type(value) is numpy.ndarray:
                    value = value.tolist()
                elif type(value) is numpy.bool_:
                    value = value.item()

                self.__dict__[key] = value

        return self

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return self.__str__()

    @classmethod
    def easy_from_hdf5(cls, h):
        """
        Create object from hdf5 file without knowing this class's arguments.
        :param h: hdf5 file pointer
        :return:
        """
        return cls().from_hdf5(h)


class DataSet(object):
    def __init__(self, h, ds_name, objs, dtype=numpy.float64):
        """
        Store data as data sets in hdf5 file.
        :param h: hdf5 file
        :param ds_name: name of data set
        :param objs: objects needing to be stored
        """
        self.dtype = dtype
        self.h = h
        self.ds_name = ds_name
        self.objs = objs
        self.width = None

        self.validate()

    def validate(self):
        assert False, 'Validate has not been implemented for this class.'

    def unsigned_bits(self, bits):
        if self.objs is not None:
            return all(int(x).bit_length() <= bits and x >= 0 for x in self.objs)

        return True

    def array_2d(self):
        if self.objs is not None:
            return numpy.array(self.objs).ndim == 2
        return True

    def encode(self):
        """
        Encode a list of data.
        :return: None
        """
        if self.objs is not None and list(self.objs):
            if all(x == self.objs[0] for x in self.objs):
                self.h.attrs[self.ds_name] = self.objs[0]
            else:
                self.objs = [self.objs]
                self.encode_list()

    def decode(self):
        """
        Decode a list of data.
        :return: list of data
        """
        l = self.decode_list()
        if l:
            return l[0]

    def encode_list(self):
        """
        Encode a list of lists.
        :return: None
        """
        if self.objs is None and not list(self.objs):
            return

        self.width = len(self.objs[0])
        self.resize()
        ds = self.dataset()
        index = ds.attrs['index']
        start_index = index

        for obj in self.objs:
            ds[index] = obj
            index += 1

        self.h.attrs[self.ds_name] = ds.regionref[start_index:index]

        ds.attrs['index'] = index

    def decode_list(self):
        """
        Decode a list of lists.
        :return: list of stored data
        """
        if self.ds_name in self.h.attrs:
            if isinstance(self.h.attrs[self.ds_name], RegionReference):
                ref = self.h.attrs[self.ds_name]
                if ref:
                    return self.dataset()[ref].tolist()
            else:
                experiment_length = len(self.h['/'].attrs[EXPERIMENT_NAMES].tolist())
                return [[self.h.attrs[self.ds_name]] * experiment_length]

    def resize(self):
        ds = self.dataset()
        curr_length = len(ds)
        ds.resize((curr_length + len(self.objs), self.width))

    def dataset(self):
        datasets_folder = '/voila_datasets'
        try:
            dataset = self.h[join(datasets_folder, self.ds_name)]
        except KeyError:
            dataset = self.h.create_dataset(
                name=join(datasets_folder, self.ds_name),
                shape=(0, self.width),
                dtype=self.dtype,
                chunks=(1, self.width),
                maxshape=(None, self.width),
                compression="gzip",
                compression_opts=9,
                shuffle=True,
            )

            # store how many objects we've worked with in the HDF5 file
            dataset.attrs['index'] = 0

        return dataset


class BinsDataSet(DataSet):
    def __init__(self, h, objs=None):
        """
        VoilaLsv bins dataset.
        :param h: HDF5 file object
        :param objs: objects to be stored as Bins Data
        """
        super(BinsDataSet, self).__init__(h, 'bins', objs)

    def validate(self):
        assert self.array_2d(), 'Bins data must be 2d array.'


class Psi1DataSet(DataSet):
    def __init__(self, h, objs=None):
        """
        VoilaLsv PSI1 dataset
        :param h: HDF5 file object
        :param objs: objects to be stored as Psi1 data
        """
        super(Psi1DataSet, self).__init__(h, 'psi1', objs)

    def validate(self):
        assert self.array_2d(), 'Psi1 data must be 2d array.'


class Psi2DataSet(DataSet):
    def __init__(self, h, objs=None):
        """
        VoilaLsv PSI2 dataset
        :param h:  HDF5 file object
        :param objs: Objects to be stored as Psi2 data
        """
        super(Psi2DataSet, self).__init__(h, 'psi2', objs)

    def validate(self):
        assert self.array_2d(), 'Psi2 data must be 2d array.'


class ExonTypeDataSet(DataSet):
    def __init__(self, h, objs=None):
        """
        ExonGraphic's type list
        :param h: HDF5 file pointer
        :param objs: List to be stored
        """
        # sanity check
        super(ExonTypeDataSet, self).__init__(h, 'exon_type', objs, dtype=numpy.uint8)

    def validate(self):
        assert self.unsigned_bits(8), 'Exon Type must be unsigned 8 bit integer.'


class JunctionTypeDataSet(DataSet):
    def __init__(self, h, objs=None):
        """
        JunctionGraphic's type list
        :param h: HDF5 file pointer
        :param objs: junction type list
        """
        # sanity check
        super(JunctionTypeDataSet, self).__init__(h, 'junction_type', objs, dtype=numpy.uint8)

    def validate(self):
        assert self.unsigned_bits(8), 'Junction Type must be unsigned 8 bit integer'


class ReadsDataSet(DataSet):
    def __init__(self, h, objs=None):
        """
        JuctionGraphic's reads list
        :param h:  HDF5 file pointer
        :param objs: reads list
        """
        super(ReadsDataSet, self).__init__(h, 'reads', objs, dtype=numpy.uint32)

    def validate(self):
        assert self.unsigned_bits(32), 'Reads data must be unsigned 32 bit integer.'
