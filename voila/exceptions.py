import argparse


class VoilaException(Exception):
    pass


class NoLsvsFound(VoilaException):
    def __init__(self):
        """
        No LSVs found.
        """
        m = "There are no LSVs found.  It could be the threshold is too high or the filters are incorrectly set."
        super(NoLsvsFound, self).__init__(m)


class GeneIdNotFoundInVoilaFile(VoilaException):
    def __init__(self, gene_id):
        """
        Error thrown when Gene ID cannot be foudn in Voila file.
        :param gene_id: 
        """
        m = 'Gene ID "{0}" was not found in Voila file'.format(gene_id)
        super(GeneIdNotFoundInVoilaFile, self).__init__(m)


class LsvIdNotFoundInVoilaFile(VoilaException):
    def __init__(self, lsv_id):
        """
        Error thrown when Gene ID cannot be foudn in Voila file.
        :param lsv_id: 
        """
        m = 'LSV ID "{0}" was not found in Voila file'.format(lsv_id)
        super(LsvIdNotFoundInVoilaFile, self).__init__(m)


class AttributeNotFoundInVoilaFile(VoilaException):
    def __init__(self, lsv_id, attr):
        """
        Error thrown when Gene ID cannot be foudn in Voila file.
        :param lsv_id:
        """
        m = 'LSV ID {0} the attribute "{1}" was not found in Voila file'.format(lsv_id, attr)
        super(AttributeNotFoundInVoilaFile, self).__init__(m)


class NoExonsInGene(VoilaException):
    def __init__(self, gene_id):
        """
        Thrown when gene, for some reason, doesn't have any exons...
        """
        super(NoExonsInGene, self).__init__('There are no exons in gene {0}'.format(gene_id))


class CouldNotBeConverted(VoilaException):
    def __init__(self, values):
        super(CouldNotBeConverted, self).__init__(
            'Meta Data {0} could not be converted to indices'.format(', '.join(values)))


class GeneIdNotFoundInSpliceGraphFile(VoilaException):
    def __init__(self, gene_id):
        """
        Error thrown when Gene ID cannot be found in Splice Graph file.
        :param gene_id: 
        """
        m = 'Gene ID {0} was not found in Splice Graph file'.format(gene_id)
        super(GeneIdNotFoundInSpliceGraphFile, self).__init__(m)


class ExperimentIndexError(VoilaException):
    def __init__(self):
        """
        Error thrown when experiment index is out of range.
        """
        super(ExperimentIndexError, self).__init__(
            'Attempted to access an out of range experiment.')


class ExperimentUnknowField(VoilaException):
    def __init__(self, field):
        """
        Error thrown when attempting to access a class attribute that hasn't been marked as an "experiment" attribute.
        :param field:
        """
        super(ExperimentUnknowField, self).__init__(
            '"{0}" might not contain data that has been sorted into a list by experiment.'.format(field))


class InValidAnalysisType(VoilaException):
    def __init__(self):
        m = 'Analysis type is not valid for the Voila file.'
        super(InValidAnalysisType, self).__init__(m)


class NotNumpyObject(VoilaException):
    def __init__(self, obj):
        m = 'Must be a numpy object: {0} - {1}'.format(type(obj), obj)
        super(NotNumpyObject, self).__init__(m)


class NotCorrectUnsignedBits(VoilaException):
    pass


class VoilaCantFindFile(argparse.ArgumentTypeError):
    def __init__(self, value):
        super(VoilaCantFindFile, self).__init__('cannot find file "{0}"'.format(value))


class NotPsiVoilaFile(VoilaException):
    def __init__(self, file_name):
        m = 'Voila file has not been quantified using PSI: ' + file_name
        super().__init__(m)

class NotDeltaPsiVoilaFile(VoilaException):
    def __init__(self, file_name):
        m = 'Voila file has not been quantified using DeltaPSI: ' + file_name
        super().__init__(m)