import argparse


class VoilaException(Exception):
    pass


class NoLsvsFound(VoilaException):
    def __str__(self):
        return "There are no LSVs found.  It could be the threshold is too high or the filters are incorrectly set."


class GeneIdNotFoundInVoilaFile(VoilaException):
    def __init__(self, filename, gene_id):
        """
        Error thrown when Gene ID cannot be foudn in Voila file.
        :param gene_id:
        """
        self.filename = filename
        self.gene_id = gene_id

    def __str__(self):
        return '{1}: Gene ID "{0}" was not found in Voila file'.format(self.gene_id, self.filename)


class LsvIdNotFoundInVoilaFile(VoilaException):
    def __init__(self, filename, lsv_id):
        """
        Error thrown when Gene ID cannot be found in Voila file.
        :param lsv_id:
        """
        self.filename = filename
        self.lsv_id = lsv_id

    def __str__(self):
        return '{}: LSV ID "{}" was not found in Voila file'.format(self.filename, self.lsv_id)


class CanNotFindFile(argparse.ArgumentTypeError):
    def __init__(self, value):
        """
        Specified file could not be found.
        :param value:
        """
        super(CanNotFindFile, self).__init__('cannot find "{0}"'.format(value))


class UnknownAnalysisType(VoilaException):
    def __init__(self, analysis_type):
        """
        Supplied analysis type (psi, delta psi, etc.) was not expected.
        :param analysis_type:
        """
        super().__init__('Unknown analysis type: ' + str(analysis_type))


class IndexNotFound(VoilaException):
    def __init__(self):
        """
        Somehow, index was not found when it should be there.
        """
        super().__init__('Index was not found.')


class SortFunctionNotFound(VoilaException):
    def __init__(self, column_idx):
        """
        While attempting to sort a column for a datatable, a sort function could not be found.  This could be that the
        column was no listed as a sortable column and it is (or visa versa).  Also this could be a column that should
        listed in the extra sort dictionary and isn't. 
        """
        super().__init__(
            'Column {} in datatables list requires specified sort function that was not found.'.format(column_idx))


class FoundNoSpliceGraphFile(VoilaException):
    def __init__(self):
        """
        Within the directories and files supplied, there were no splice graphs found.
        """
        super().__init__('Splice Graph file was not found.')


class FoundMoreThanOneSpliceGraph(VoilaException):
    def __init__(self):
        """
        Within the directories and files supplied, more then on splice graph was found.
        """
        super().__init__('Found more then one Splice Graph file.')


class MixedAnalysisTypeVoilaFiles(VoilaException):
    def __init__(self):
        """
        Within the directories and files supplied, more then on voila file was found and they contain more then on
        analysis type.
        """
        super().__init__('Found Voila files have more then one analysis type.')


class FoundMoreThanOneVoilaFile(VoilaException):
    def __init__(self):
        """
        Found more then one voila file when there should only be one voila file.
        """
        super().__init__('In the files or directories supplied, there was more than on Voila file found.')


class AnalysisTypeNotFound(VoilaException):
    def __init__(self):
        """
        There was no analysis type found, which probably means there were no voila files found.
        """
        super().__init__('No Voila files were found in the files or directories provided.')
