import os
import configparser
from scipy import interpolate
import numpy as np
from majiq.src.constants import *


class SingletonMetaClass(type):
    def __init__(cls, name, bases, dict):
        super(SingletonMetaClass, cls).__init__(name,bases, dict)
        original_new = cls.__new__

        def my_new(cls, *args, **kwds):

            if cls.instance:
                cls.instance = original_new(cls, *args, **kwds)
            return cls.instance
        cls.instance = None
        cls.__new__ = staticmethod(my_new)


class Config(object):

    instance = None

    def __new__(cls, *argv):
        if not Config.instance:
            Config.instance = Config.__Config(*argv)
        return Config.instance

    def __getattr__(self, name):
        return getattr(self.instance, name)

    def __setattr__(self, name):
        return setattr(self.instance, name)

    class __Config(object):

        def _set_strandness(self, experiment_name, val):
            self.strand_specific[experiment_name] = val

        def __init__(self, filename, params, only_db=False):

            self.__dict__.update(params.__dict__)

            if not os.path.exists(self.outDir):
                os.makedirs(self.outDir)

            config = configparser.ConfigParser()
            config.read(filename)

            general = Config.config_section_map(config, "info")
            self.tissue_repl = {}
            self.exp_list = []
            count = 0
            if not os.path.exists(self.outDir):
                os.makedirs(self.outDir)
            self.sam_dir = general['samdir']
            self.genome = general['genome']
            self.readLen = int(general['readlen'])

            # self.simplify_threshold = 0.0
            # self.simplify_type = SIMPLIFY_ALL
            # if self.simplify:
            #
            #     self.simplify_threshold = float(params.simplify[1])
            #     self.simplify_type = self.simplify[0]
            #     if self.simplify_type not in (
            #     SIMPLIFY_ALL, SIMPLIFY_DB, SIMPLIFY_DENOVO) or not 0 <= self.simplify_threshold <= 1:
            #         raise RuntimeError(
            #             'Error in simplify option, first argument should be "all|denovo|annotated" and second'
            #             ' a float between 0..1')

            exps = Config.config_section_map(config, "experiments")
            self.juncfile_list = []
            for exp_idx, lstnames in exps.items():
                self.tissue_repl[exp_idx] = []
                # self.juncfile_list[exp_idx] = []
                elist = lstnames.split(',')
                for exp in elist:
                    self.exp_list.append(exp)
                    self.tissue_repl[exp_idx].append(count)
                    count += 1

            self.num_experiments = len(self.exp_list)
            self.samfile_name_list = []
            self.sam_list = []

            for name, ind_list in self.tissue_repl.items():
                for exp_idx in ind_list:
                    exp = self.exp_list[exp_idx]
                # for exp_idx, exp in enumerate(self.exp_list):
                    samfile = "%s/%s.%s" % (self.sam_dir, exp, SEQ_FILE_FORMAT)
                    if not os.path.exists(samfile):
                        raise RuntimeError("Skipping %s.... not found" % samfile)
                    baifile = "%s/%s.%s" % (self.sam_dir, exp, SEQ_INDEX_FILE_FORMAT)
                    if not os.path.exists(baifile):
                        raise RuntimeError("Skipping %s.... not found ( index file for bam file is required)" % baifile)
                    juncfile = "%s/%s.%s" % (self.sam_dir, exp, JUNC_FILE_FORMAT)
                    self.sam_list.append((exp, os.path.exists(juncfile), name))

            opt_dict = {'strandness': self._set_strandness}
            strandness = {'forward': FWD_STRANDED, 'reverse': REV_STRANDED, 'none': UNSTRANDED}
            if 'strandness' in general:
                try:
                    global_strand = strandness[general['strandness'].lower()]
                except:
                    raise RuntimeError('Incorrect Strand-specific option [forward, reverse, none]')
            else:
                global_strand = strandness['none']
            self.strand_specific = {xx: global_strand for xx in self.exp_list}

            opt = Config.config_section_map(config, "opts")
            for exp_id, opts_list in opt.items():
                elist = opts_list.split(',')
                for opt in elist:
                    op_id, op_val = opt.split(':')
                    try:
                        opt_dict[op_id](exp_id, op_val)
                    except KeyError:
                        raise RuntimeError('Option %s do not exist. The options available '
                                           'are %s' % (op_id, ','.join(opt_dict.keys())))

            return

        def __str__(self):
            return repr(self) + self.val

    @staticmethod
    def config_section_map(config_d, section):
        dict1 = {}
        try:
            options = config_d.options(section)
        except configparser.NoSectionError:
            return dict1
        for option in options:
            try:
                dict1[option] = config_d.get(section, option)
                if dict1[option] == -1:
                    print("skip: %s" % option)
            except:
                print("exception on %s!" % option)
                dict1[option] = None
        return dict1
