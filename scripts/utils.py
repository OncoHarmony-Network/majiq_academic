import fnmatch
import os
from matplotlib import pyplot
import numpy

__author__ = 'abarrera'


def get_mean_step(bins):
    bins = numpy.array(bins)
    step = 1.0 / bins.size
    projection_prod = bins * numpy.arange(step / 2, 1, step)
    return numpy.sum(projection_prod)


def _save_or_show(plotpath=None, name=None):
    if plotpath:
        if os.path.isdir(plotpath):
            plot_base_path, plot_name = plotpath, name
        else:
            plot_base_path, plot_name = os.path.split(plotpath)
            if not os.path.exists(plot_base_path):
                os.makedirs(plot_base_path)
            if not plot_name:
                plot_name = name
        pyplot.savefig("%s/%s.png"%(plot_base_path, plot_name), width=300, height=300, dpi=100)
        print "Saved in:\n%s/%s" % (plot_base_path, plot_name)

        pyplot.clf()
    else:
        pyplot.show()


def list_files_or_dir(file_or_dir_list, suffix='*'):
    if type(file_or_dir_list) != list: return file_or_dir_list

    files = []
    for file_or_dir in file_or_dir_list:
        if os.path.isdir(file_or_dir):
            for root, dirnames, filenames in os.walk(file_or_dir):
                for filename in fnmatch.filter(filenames, '*%s' % suffix):
                    files.append(os.path.join(root, filename))
            # for file in os.listdir(file_or_dir):
            #     if not suffix or file.endswith(suffix):
            #         files.append(file_or_dir+'/'+file)
        else:
            files.append(file_or_dir)
    return files


def miso_delta_reader(path, filter_complex=True, result_dict=None):
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
            bayes_factor = sline[8]
            event_name = sline[0]

            if not filter_complex or len(transcripts) == 2:  # only interested in 2 transcripts events for now
                if result_dict is not None:
                    for i, dpsi in enumerate(delta_psi):
                        result_dict["%s#%d" % (event_name, i)].append(float(dpsi))
                else:
                    ret.append([event_name, float(delta_psi[0]), float(bayes_factor)])
    return ret
