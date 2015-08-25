"""
Create the intersection between the MAJIQ, MISO and MATS events
"""
import argparse

from pylab import *
from scipy.sparse import lil_matrix


def load_names(path, prev_dict):
    data = pickle.load(open(path))
    grimoire_obj = data[1][:,0]
    for i, junction in enumerate(grimoire_obj):
        if type(junction.coverage[0]) == lil_matrix:
            prev_dict[junction.name] = junction.id

    return prev_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--majiqfiles', nargs='+')
    parser.add_argument('--majiqnames', nargs='+')
    args = parser.parse_args()
    
    print "Building translation dictionary..."
    newname = {}
    for path in args.majiqfiles:
        print path
        newname = load_names(path, newname)

    print "Reading event names..."
    corrected_paths = [open(path+'.correct', 'w') for path in args.majiqnames] #paths to write the selected matrices
    for i, path in enumerate(args.majiqnames):
        corrected_names = []
        for name in pickle.load(open(path)):
            try:
                corrected_names.append(newname[name])
            except KeyError:
                corrected_names.append(name)

        print "Dumping..."
        pickle.dump(corrected_names, corrected_paths[i])

if __name__ == '__main__':
    main()