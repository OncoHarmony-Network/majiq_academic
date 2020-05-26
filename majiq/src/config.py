"""
config.py

Class for parsing configuration files for MAJIQ build wrapping configparser
"""

import os
import configparser
import numpy as np
import majiq.src.constants as constants
import warnings


class SingletonMetaClass(type):
    def __init__(cls, name, bases, dict):
        super(SingletonMetaClass, cls).__init__(name, bases, dict)
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

        def __init__(self, filename, params):

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

            try:
                sam_dirlist = general["bamdirs"].split(",")
            except KeyError:
                if "samdir" in general:
                    raise UserWarning(
                        "samdir is a deprecated value, please use bamdirs instead"
                    )
                sam_dirlist = ["."]
                warnings.warn(
                    'bamdirs parameter not found in config file, using "./" instead'
                )

            try:
                junc_dirlist = general["sjdirs"].split(",")
            except KeyError:
                junc_dirlist = ["."]
                warnings.warn(
                    'sjdirs parameter not found in config file, using "./" instead'
                )

            self.genome = general["genome"]

            if "readlen" in general:
                # removing it from function/class definitions/initializations
                warnings.warn(
                    '"readlen" parameter is deprecated and will not be used.'
                    " MAJIQ now detects the maximum read length of each"
                    " experiment automatically. "
                )

            exps = Config.config_section_map(config, "experiments")
            self.juncfile_list = []
            for exp_idx, lstnames in exps.items():
                self.tissue_repl[exp_idx] = []
                elist = lstnames.split(",")
                for exp in elist:
                    self.exp_list.append(exp)
                    self.tissue_repl[exp_idx].append(count)
                    count += 1

            self.num_experiments = len(self.exp_list)
            self.samfile_name_list = []
            self.sam_list = []

            self.min_experiments = {}
            for name, ind_list in self.tissue_repl.items():
                if self.min_exp < 1:
                    mexp = np.ceil(len(ind_list) * self.min_exp)
                else:
                    mexp = self.min_exp
                self.min_experiments[name] = int(min(len(ind_list), mexp))
                for exp_idx in ind_list:
                    # try to find input file with this prefix from input folders
                    prefix = self.exp_list[exp_idx]
                    found = False  # update once found
                    if self.aggregate:
                        for j_dir in junc_dirlist:
                            juncfile = f"{j_dir}/{prefix}.{constants.JUNC_FILE_FORMAT}"
                            if os.path.isfile(juncfile):
                                found = True
                                self.sam_list.append(
                                    (self.exp_list[exp_idx], juncfile, True)
                                )
                                break
                    if found:
                        continue
                    for s_dir in sam_dirlist:
                        bamfile = f"{s_dir}/{prefix}.{constants.SEQ_FILE_FORMAT}"
                        baifile = f"{s_dir}/{prefix}.{constants.SEQ_INDEX_FILE_FORMAT}"

                        if os.path.isfile(bamfile) and os.path.isfile(baifile):
                            found = True
                            self.sam_list.append(
                                (self.exp_list[exp_idx], bamfile, False)
                            )
                            break

                    if not found:
                        raise RuntimeError(
                            f"No matching {constants.SEQ_FILE_FORMAT}"
                            f" (and {constants.SEQ_INDEX_FILE_FORMAT})"
                            f" or {constants.JUNC_FILE_FORMAT} file was found"
                            f" matching the experiment {prefix} in any of"
                            " the provided paths in the build configuration"
                        )

            opt_dict = {"strandness": self._set_strandness}
            strandness = {
                "forward": constants.FWD_STRANDED,
                "reverse": constants.REV_STRANDED,
                "none": constants.UNSTRANDED,
            }
            if "strandness" in general:
                try:
                    global_strand = strandness[general["strandness"].lower()]
                except Exception:
                    raise RuntimeError(
                        "Incorrect Strand-specific option [forward, reverse, none]"
                    )
            else:
                global_strand = strandness["none"]
            self.strand_specific = {xx: global_strand for xx in self.exp_list}

            opt = Config.config_section_map(config, "opts")
            for exp_id, opts_list in opt.items():
                elist = opts_list.split(",")
                for opt in elist:
                    op_id, op_val = opt.split(":")
                    try:
                        opt_dict[op_id](exp_id, op_val)
                    except KeyError:
                        raise RuntimeError(
                            "Option %s do not exist. The options available "
                            "are %s" % (op_id, ",".join(opt_dict.keys()))
                        )

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
            except Exception:
                print("exception on %s!" % option)
                dict1[option] = None
        return dict1
