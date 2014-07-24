import os
from matplotlib import pyplot
import numpy

__author__ = 'abarrera'


def get_mean_step(bins):
    bins = numpy.array(bins)
    step = 1.0 / bins.size
    projection_prod = bins * numpy.arange(step / 2, 1, step)
    return numpy.sum(projection_prod)


def _save_or_show(plotpath, name):
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