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
        super(CanNotFindFile, self).__init__('cannot find "{0}"'.format(value))


class UnknownAnalysisType(VoilaException):
    def __init__(self, analysis_type):
        super().__init__('Unknown analysis type: ' + str(analysis_type))


class IndexNotFound(VoilaException):
    pass


class SortFunctionNotFound(VoilaException):
    pass


class FoundNoSpliceGraphFile(VoilaException):
    pass


class FoundMoreThanOneSpliceGraph(VoilaException):
    pass


class MixedAnalysisTypeVoilaFiles(VoilaException):
    pass


class FoundMoreThanOneVoilaFile(VoilaException):
    def __init__(self):
        super().__init__('In the files or directories supplied, there was more than on Voila file found.')


class AnalysisTypeNotFound(VoilaException):
    def __init__(self):
        super().__init__('No Voila files were found in the files or directories provided.')
