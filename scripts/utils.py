from collections import defaultdict
import fnmatch
import os
import pickle
from matplotlib import pyplot
import numpy
import numpy as np

__author__ = 'abarrera'


def get_mean_step(bins):
    bins = numpy.array(bins)
    step = 1.0 / bins.size
    projection_prod = bins * numpy.arange(step / 2, 1, step)
    return numpy.sum(projection_prod)


def save_or_show(plotpath=None, name=None, exten='png'):
    if plotpath:
        if os.path.isdir(plotpath):
            plot_base_path, plot_name = plotpath, name
        else:
            plot_base_path, plot_name = os.path.split(plotpath)
            if not os.path.exists(plot_base_path):
                os.makedirs(plot_base_path)
            if not plot_name:
                plot_name = name
        pyplot.savefig("%s/%s.%s"%(plot_base_path, plot_name, exten), width=300, height=300, dpi=100)
        print "Saved in:\n%s/%s.%s" % (plot_base_path, plot_name, exten)

        pyplot.clf()
    else:
        pyplot.show()


def list_files_or_dir(file_or_dir_list, suffix='*', containing='*'):
    if type(file_or_dir_list) != list: return file_or_dir_list

    files = []
    for file_or_dir in file_or_dir_list:
        if os.path.isdir(file_or_dir):
            for root, dirnames, filenames in os.walk(file_or_dir):
                for filename in fnmatch.filter(filenames, '*%s*%s' % (containing, suffix)):
                    files.append(os.path.join(root, filename))
            # for file in os.listdir(file_or_dir):
            #     if not suffix or file.endswith(suffix):
            #         files.append(file_or_dir+'/'+file)
        else:
            files.append(file_or_dir)
    return files


def miso_delta_reader(path, dofilter=False, complex_lsvs=False, result_dict=None, from_majiq=True):
    """Read delta psi calculations from MISO."""
    ret = []
    for line in open(path):
        sline = line.split('\t')
        if sline[0] != "event_name":
            #event_name      sample1_posterior_mean  sample1_ci_low  sample1_ci_high sample2_posterior_mean  sample2_ci_low  sample2_ci_high diff    bayes_factor    isoforms        sample1_counts  sample1_assigned_counts sample2_counts  sample2_assigned_counts chrom   strand  mRNA_starts     mRNA_ends
            #ENST00000301332:145697325-145697429:target      0.98    0.91    1.00    0.98    0.92    1.00    0.00    0.05    'ENST00000301332:145697325-145697429:target.0.ex_ENST00000301332:145697325-145697429:target.0.lsv','ENST00000301332:145697325-145697429:target.1.ex_ENST00000301332:145697325-145697429:target.1.lsv'   (0,0):36,(1,0):37       0:37    (0,0):37,(1,0):47,(1,1):3       0:50    chr8    +       145694873,145695965     145697429,145697429
            try:
                transcripts = sline[9].split(",")
            except:
                print line
                raise
            delta_psi = sline[7].split(",")
            bayes_factor = sline[8].split(",")
            event_name = sline[0]
            
            if not from_majiq:
                if result_dict is not None: 
                    result_dict.append(float(delta_psi[0]))
                else:
                    ret.append([event_name, float(delta_psi[0]), float(bayes_factor[0]), 1])
            else: 
                if complex_lsvs or len(transcripts) == 2:  # only interested in 2 transcripts events for now
                    if result_dict is not None:
                        for i, dpsi in enumerate(delta_psi):
                            result_dict["%s#%d" % (event_name, i)].append(float(dpsi))
                        if len(delta_psi) == 1:
                            result_dict["%s#1" % event_name].append(-float(delta_psi[0]))
                    else:
                        max_dpsi = 0
                        max_junc = 0
                        for i, dpsi in enumerate(delta_psi):
                            if abs(float(delta_psi[i])) > max_dpsi:
                                max_dpsi = abs(float(delta_psi[i]))
                                max_junc = i
                        ret.append([event_name, float(delta_psi[max_junc]), float(bayes_factor[max_junc]), 1])
    return ret


def coverage_from_file(file, cov_suffix):
    result_dict = defaultdict(list)  # First the num. reads, second the positions
    with open(file) as majiq_file:
        majiq_builder = pickle.load(majiq_file)
        for lsv in majiq_builder[1]:
            for i, junc in enumerate(lsv.junction_list):
                result_dict[lsv.id+"#"+str(i)].append([np.sum(junc.data[0]), junc.nnz])
    # Save results
    pickle.dump(result_dict, open(file+cov_suffix, 'w'))
    print "Coverage saved in: %s" % file+cov_suffix
