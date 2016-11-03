import numpy

from voila import constants
from voila.utils.voilaLog import voilaLog

VOILA_FILE_VERSION = '/voila_file_version'


class HDF5VersionException(Exception):
    def __init__(self):
        super(HDF5VersionException, self).__init__('The hdf5 file version does not match the current version of Voila.')
        voilaLog().error(self.message)


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
        try:
            h[VOILA_FILE_VERSION] = constants.FILE_VERSION
        except RuntimeError:
            pass

        attrs_dict = self.__dict__.copy()

        for key in self.exclude():
            del attrs_dict[key]

        for cls in self.cls_list():
            cls_grp = h.create_group(cls)
            for index, cls_obj in enumerate(attrs_dict[cls]):
                cls_obj.to_hdf5(cls_grp.create_group(str(index)), use_id)
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

        if h[VOILA_FILE_VERSION].value != constants.FILE_VERSION:
            raise HDF5VersionException()

        cls_dict = self.cls_list()

        for cls in cls_dict:
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

    def encode(self):
        self.objs = [self.objs]
        self.encode_list()

    def decode(self):
        return self.decode_list()[0]

    def encode_list(self):
        """
        Encode attribute data into datasets.
        :return: None
        """
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
        Decode stored data
        :return: list of stored data
        """
        ref = self.h.attrs[self.ds_name]
        return self.dataset()[ref].tolist()

    def resize(self):
        ds = self.dataset()
        curr_length = len(ds)
        ds.resize((curr_length + len(self.objs), self.width))

    def dataset(self):
        try:
            dataset = self.h['/datasets/' + self.ds_name]
        except KeyError:
            dataset = self.h.create_dataset(
                name='/datasets/' + self.ds_name,
                shape=(0, self.width),
                dtype=self.dtype,
                chunks=(1, self.width),
                maxshape=(None, self.width),
                compression="lzf",
                shuffle=True,
            )

            # store how many objects we've worked with in the HDF5 file
            dataset.attrs['index'] = 0

        return dataset


class BinsDataSet(DataSet):
    def __init__(self, h, objs=((),)):
        """
        VoilaLsv bins dataset.
        :param h: HDF5 file object
        :param objs: objects to be stored as Bins Data
        """
        super(BinsDataSet, self).__init__(h, 'bins', objs)


class Psi1DataSet(DataSet):
    def __init__(self, h, objs=((),)):
        """
        VoilaLsv PSI1 dataset
        :param h: HDF5 file object
        :param objs: objects to be stored as Psi1 data
        """

        super(Psi1DataSet, self).__init__(h, 'psi1', objs)


class Psi2DataSet(DataSet):
    def __init__(self, h, objs=((),)):
        """
        VoilaLsv PSI2 dataset
        :param h:  HDF5 file object
        :param objs: Objects to be stored as Psi2 data
        """
        super(Psi2DataSet, self).__init__(h, 'psi2', objs)


class ExonTypeDataSet(DataSet):
    def __init__(self, h, objs=()):
        """
        ExonGraphic's type list
        :param h: HDF5 file pointer
        :param objs: List to be stored
        """
        super(ExonTypeDataSet, self).__init__(h, 'exon_type', objs, dtype=numpy.int8)


class JunctionTypeDataSet(DataSet):
    def __init__(self, h, objs=()):
        """
        JunctionGraphic's type list
        :param h: HDF5 file pointer
        :param objs:
        """
        super(JunctionTypeDataSet, self).__init__(h, 'junction_type', objs, dtype=numpy.int8)


class ReadsDataSet(DataSet):
    def __init__(self, h, objs=()):
        super(ReadsDataSet, self).__init__(h, 'reads', objs, dtype=numpy.int32)


class CleanReadsDataSet(DataSet):
    def __init__(self, h, objs=()):
        super(CleanReadsDataSet, self).__init__(h, 'clean_reads', objs, dtype=numpy.int32)
