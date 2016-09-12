import h5py
import pickle
from voila.splice_graphics import ExonGraphic, LsvGraphic, JunctionGraphic

def visual_hdf5_to_pickle(visual):

    junc_list =[]
    for jid in visual['junctions'].keys():
        junc = visual['junctions/%s'% jid]
        jj = JunctionGraphic(junc.attrs['coords'], junc.attrs['type'], junc.attrs['num_reads'],
                                         transcripts=junc.attrs['transcripts'], ir=junc.attrs['ir'])
        junc_list.append(jj)

    exon_list = []
    for exid in visual['exons'].keys():
        ex = visual['exons/%s'% exid]
        eg = ExonGraphic(ex.attrs['a3'], ex.attrs['a5'], ex.attrs['coords'], type_exon=ex.attrs['lsv_type'],
                         coords_extra=ex.attrs['coords_extra'], intron_retention=ex.attrs['ir'],
                         alt_starts=ex.attrs['alt_starts'], alt_ends=ex.attrs['alt_ends'])
        exon_list.append(eg)

    splice_lsv = LsvGraphic(type_lsv=visual.attrs['type'], coords=visual.attrs['coords'], id=visual.attrs['id'],
                            name=visual.attrs['name'],
                            strand=visual.attrs['strand'], exons=exon_list, junctions=junc_list,
                            chrom=visual.attrs['chrom'])
    return splice_lsv



if __name__ == '__main__':

    import sys
    import os

    filename = sys.argv[1]

    output_filename = sys.argv[2]

    if not os.path.exists(filename):
        print "ERROR INPUT FILE DOESN'T EXISTS"


    dict_lsv = {}
    data = h5py.File(filename)
    for lsvid in data['LSVs'].keys():
        lsv = data['LSVs/%s' % lsvid]
        dict_lsv[lsv.attrs['id']] = visual_hdf5_to_pickle(lsv['visuals'])

    pickle.dump(dict_lsv, open(output_filename,'w+'))








