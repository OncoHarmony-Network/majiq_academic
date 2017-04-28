import os
import ConfigParser
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
        def __init__(self, filename, params, only_db=False):

            self.__dict__.update(params.__dict__)

            if not os.path.exists(self.outDir):
                os.makedirs(self.outDir)

            config = ConfigParser.ConfigParser()
            config.read(filename)
            # TODO: check if filename exists
            exps = Config.config_section_map(config, "experiments")
            general = Config.config_section_map(config, "info")
            self.tissue_repl = {}
            self.exp_list = []
            count = 0
            if not os.path.exists(self.outDir):
                os.makedirs(self.outDir)

            self.sam_dir = general['samdir']
            self.genome = general['genome']
            self.genome_path = general['genome_path']
            self.readLen = int(general['readlen'])

            if 'type' in general:
                self.strand_specific = (general['type'] == 'strand-specific')
            else:
                self.strand_specific = False

            self.simplify_threshold = 0.0
            self.simplify_type = SIMPLIFY_ALL
            if self.simplify:

                self.simplify_threshold = float(params.simplify[1])
                self.simplify_type = self.simplify[0]
                if self.simplify_type not in (
                SIMPLIFY_ALL, SIMPLIFY_DB, SIMPLIFY_DENOVO) or not 0 <= self.simplify_threshold <= 1:
                    raise RuntimeError(
                        'Error in simplify option, first argument should be "all|denovo|annotated" and second'
                        ' a float between 0..1')
            #self.simplify = self.simplify is not None

            for exp_idx, lstnames in exps.items():
                self.tissue_repl[exp_idx] = []
                elist = lstnames.split(',')
                for exp in elist:
                    self.exp_list.append(exp)
                    self.tissue_repl[exp_idx].append(count)
                    count += 1

            self.num_experiments = len(self.exp_list)
            self.gene_tlb = {}

            self.samfile_name_list = []
            self.sam_list = []
            for exp_idx, exp in enumerate(self.exp_list):
                samfile = "%s/%s.bam" % (self.sam_dir, exp)
                if not os.path.exists(samfile):
                    raise RuntimeError("Skipping %s.... not found" % samfile)
                baifile = "%s/%s.bam.bai" % (self.sam_dir, exp)
                if not os.path.exists(baifile):
                    raise RuntimeError("Skipping %s.... not found ( index file for bam file is required)" % baifile)
                # self.samfile_name_list.append(exp)
                self.sam_list.append(exp)
#                self.exp_list[exp_idx] = os.path.split(exp)[1]

            return

        def __str__(self):
            return `self` + self.val

    @staticmethod
    def config_section_map(config_d, section):
        dict1 = {}
        options = config_d.options(section)
        for option in options:
            try:
                dict1[option] = config_d.get(section, option)
                if dict1[option] == -1:
                    print("skip: %s" % option)
            except:
                print("exception on %s!" % option)
                dict1[option] = None
        return dict1
