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
        m = 'Gene ID {0} was not found in Voila file'.format(gene_id)
        super(GeneIdNotFoundInVoilaFile, self).__init__(m)


class LsvIdNotFoundInVoilaFile(VoilaException):
    def __init__(self, gene_id):
        """
        Error thrown when Gene ID cannot be foudn in Voila file.
        :param gene_id:
        """
        m = 'LSV ID {0} was not found in Voila file'.format(gene_id)
        super(LsvIdNotFoundInVoilaFile, self).__init__(m)


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


class ExperimentUnknownField(VoilaException):
    def __init__(self, field):
        """
        Error thrown when attempting to access a class attribute that hasn't been marked as an "experiment" attribute.
        :param field:
        """
        super(ExperimentUnknownField, self).__init__(
            '"{0}" might not contain data that has been sorted into a list by experiment.'.format(field))


class DirectoryDoesNotExist(argparse.ArgumentTypeError):
    def __init__(self, directory):
        super(DirectoryDoesNotExist, self).__init__('Directory does not exist: {0}'.format(directory))


class CannotFindFile(argparse.ArgumentTypeError):
    def __init__(self, value):
        super(CannotFindFile, self).__init__('Cannot find file "{0}"'.format(value))


class EmptyFilters(VoilaException):
    def __init__(self):
        super(EmptyFilters, self).__init__('All elements in filters could not be found')


class IntronRetentionNotFound(VoilaException):
    def __init__(self, ir):
        super().__init__('Intron Retention not found: {} {}-{}'.format(ir.gene_id, ir.start, ir.end))


class JunctionNotFound(VoilaException):
    def __init__(self, junc):
        super().__init__('Junction not found: {} {}-{}'.format(junc.gene_id, junc.start, junc.end))
