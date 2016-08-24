import numpy


class HDF5(object):
    def __init__(self):
        """
        Move data in this class to and from HDF5 files.
        """
        pass

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
        :return:
        """
        return {}

    def to_hdf5(self, h):
        """
        Adds attributes from this class to HDF5 file.
        :param h: HDF5 file object
        :return: None
        """
        attrs_dict = self.__dict__.copy()

        for key in self.exclude():
            del attrs_dict[key]

        for cls in self.cls_list():
            cls_grp = h.create_group(cls)
            for index, cls_obj in enumerate(attrs_dict[cls]):
                cls_obj.to_hdf5(cls_grp.create_group(str(index)))
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
        attrs_dict = dict(h.attrs)
        cls_dict = self.cls_list()

        for cls in cls_dict:
            cls_grp = h[cls]
            self.__dict__[cls] = [None] * len(cls_grp)
            for index in cls_grp:
                new_class = cls_dict[cls]['class']
                args = cls_dict[cls]['args']
                self.__dict__[cls][int(index)] = new_class(*args).from_hdf5(cls_grp[index])

        for key in self.exclude():
            # Some of the exclude keys will be in attributes and others will be
            # in groups, we want to make sure we don't write over these values.
            try:
                del attrs_dict[key]
            except KeyError:
                pass

        for key in attrs_dict:
            # H5py stores attributes as numpy objects where it can.  Numpy objects
            # cause issues with the json conversion, therefore corner cases are handled below.
            value = attrs_dict[key]
            if type(value) is numpy.ndarray:
                value = value.tolist()
            elif type(value) is numpy.bool_:
                value = value.item()

            self.__dict__[key] = value

        return self


class DataSet(object):
    def __init__(self, h, ds_name, shape):
        """
        Encode and Decode class attributes in to a HDF5 dataset.
        :param h: HDF5 file object
        :param ds_name: name of datset
        :param shape: shape of data stored in dataset
        """
        self.h = h
        self.ds_name = ds_name

        try:
            self.ds = h['/' + ds_name]
        except KeyError:
            self.ds = h.create_dataset('/' + ds_name, shape, dtype=numpy.float64, chunks=(1, shape[1]))
            # store how many objects we've worked with in the HDF5 file
            self.ds.attrs['index'] = 0

    def encode_list(self, objs):
        """
        Encode attribute data into datasets.
        :param objs: list of values to store in dataset.
        :return: None
        """
        index = self.ds.attrs['index']
        start_index = index
        for obj in objs:
            self.ds[index] = obj
            index += 1
        self.h.attrs[self.ds_name] = self.ds.regionref[start_index:index]
        self.ds.attrs['index'] = index

    def decode_list(self):
        """
        Decode stored data
        :return: list of stored data
        """
        ref = self.h.attrs[self.ds_name]
        return self.ds[ref].tolist()


class BinsDataSet(DataSet):
    def __init__(self, h, length=None):
        """
        VoilaLsv bins dataset.
        :param h: HDF5 file object
        :param length: total number of bins' rows in VoilaInput
        """
        super(BinsDataSet, self).__init__(h, 'bins', (length, 39))


class Psi1DataSet(DataSet):
    def __init__(self, h, length=None):
        """
        VoilaLsv PSI1 dataset
        :param h: HDF5 file object
        :param length: total number PSI1 rows in VoilaInput
        """
        super(Psi1DataSet, self).__init__(h, 'psi1', (length, 20))


class Psi2DataSet(DataSet):
    def __init__(self, h, length=None):
        """
        VoilaLsv PSI2 dataset
        :param h:  HDF5 file object
        :param length: total number of PSI2 rows in VoilaInput
        """
        super(Psi2DataSet, self).__init__(h, 'psi2', (length, 20))
