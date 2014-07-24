import os
import pickle
from matplotlib import pyplot
import argparse
import utils
__author__ = 'abarrera'

def barchart_expected(expected_psis, plotpath, mfile):
     f = pyplot.figure()
     pyplot.hist(expected_psis, bins=40)
     name = os.path.basename(mfile)
     utils._save_or_show(plotpath, name +'_expected_dist')


def plot_from_file(mfile, plotpath):
    expected_psis = []
    with open(mfile) as mfile_open:
        mpickle = pickle.load(mfile_open)
        for i, lsv_info in enumerate(mpickle[0]):
            if len(mpickle[0][i])<3:
                expected_psis.append(utils.get_mean_step(mpickle[0][i][0]))
    barchart_expected(expected_psis, plotpath, mfile)


def main():
    parser = argparse.ArgumentParser(description="Plot distribution of expected PSIs.")
    parser.add_argument("filesOrDir", nargs='+', help="MAJIQ PSIs quantification files.")
    parser.add_argument('--plotpath', default='output')
    args = parser.parse_args()

    for mfileOrDir in args.filesOrDir:
        if os.path.isdir(mfileOrDir):
            for file in os.listdir(mfileOrDir):
                if file.endswith(".majiq_psi.pickle"):
                    plot_from_file(mfileOrDir + "/" +file, args.plotpath)
        else:
            plot_from_file(mfileOrDir, args.plotpath)

if __name__ == '__main__':
    main()