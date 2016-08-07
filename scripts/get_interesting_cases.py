import sys
grimoire_path = '/Volumes/jordi/working/old_majiq/'
if grimoire_path not in sys.path:
    sys.path.append(grimoire_path)

from scripts.pair_rank import rank_majiq
import pickle

import numpy as np


def read_ensamble_strand():
    fp = open('./ensembl.mm10.gff3')
    lines = fp.readlines()
    
    dd_strand = {}
    chromosome = {}
    for li in lines:
        if li.startswith('#'):
            continue
        tab = li.strip().split()
        if tab[2] != 'gene':
            continue
        for attr in tab[8].split(';'):
            tab2 = attr.split('=')
            if tab2[0] == 'ID':
                genename = tab2[1]
                break
        
        dd_strand[genename] = tab[6]
        chromosome[genename] = tab[0]
        
#    print dd_strand
    print chromosome
    return chromosome, dd_strand


def correct_type(strand_dict, gne, info_ev):

#    gene = info_ev[1].split(':')[0]
    
    if strand_dict[gne] == '-':
        if info_ev[2][0] == 't':
            new_ev = 's'+info_ev[1:]
        else:
            new_ev = 't'+info_ev[1:]
            
        info_ev = new_ev

    return info_ev


def get_read_count(list1, list2):

    dict_ev = [{}, {}]
    res_dic = {} 
    ev = set()

    for fl in list1:
        raw_reads1 = pickle.load(open(fl))
        for nm1, nr1 in raw_reads1:
            ev.add(nm1)
            if nm1 not in dict_ev:
                dict_ev[0][nm1] = []
            dict_ev[0][nm1].append(nr1.sum())

    for fl in list2:
        raw_reads2 = pickle.load(open(fl))
        for nm1, nr1 in raw_reads2:
            ev.add(nm1)
            if nm1 not in dict_ev:
                dict_ev[1][nm1] = []
            dict_ev[1][nm1].append(nr1)
    
    for nm in ev:
        avg1 = np.mean(dict_ev[0][nm])
        avg2 = np.mean(dict_ev[1][nm])
    
        res_dic[nm] = max(avg1, avg2)

    return res_dic

if __name__ == '__main__':

    reads1 = ['Cer_Liver/clean_reads.Cer_28_40_52.pkl']
    reads2 = ['Cer_Liver/clean_reads.Liv_28_40_52.pkl']

    filename = 'Cer_Liver/Cer_28_40_52_Liv_28_40_52.deltapsi.pickle'
    #file_list = [['/Volumes/Data/ucsc/test2/output_Hip1_Liv1/toJuan.Hippocampus1.old_majiq',
                    # /Volumes/Data/ucsc/test2/output_Hip1_Liv1/toJuan.Heart1.old_majiq']]
    label = ['Cer vs Liver1']
    fidx = 0
    N = 633
    matrices, matched_info, meta, psi1, psi2 = pickle.load(open(filename))
    #
    tmp = rank_majiq(matrices, matched_info, E=True)
    
    tmp.sort(key=lambda x: -x[1])
    relevant = set([xx[0][1] for xx in tmp[:N]])

    print len(relevant)

    chrom = {}
    strand = {}
    for tt in tmp:
        chrom[tt[0][4].id] = tt[0][4].chrom
        strand[tt[0][4].id] = tt[0][4].strand
                                
    dic = get_read_count(reads1, reads2)
    
    chrom, strand = read_ensamble_strand()

    bins = [range(0, 15), range(15, 20), range(20, 40), range(40, 100), range(100, 90000)]
    bins2 = [[0, 15], [15, 20], [20, 40], [40, 100], [100, 90000]]
    read_bins = [[], [], [], [], []]

    BINNAMES = ['0-15', '15-20', '20-40', '40-100', '100-xx']
    
    for name, nr in dic.items():
        if not name in relevant:
            continue
        for idx, el in enumerate(bins2):
            print el, nr
            if el[0] <= nr < el[1]:
                print "INTO"
                read_bins[idx].append(name)
                break
                    
    interest = []           
    for ll in read_bins:
        n = min(len(ll), 20)
        #interest.append(random.sample(ll,n))
        interest.append(ll)

    result = {}

    for bb in BINNAMES:
        result[bb] = []

    fall = 'Cer_Liver/Cer_CT28.mm10.sorted.old_majiq'
    alllsv = pickle.load(open(fall))
    outf = open('./event_set.tab', 'w+')
    header = '#event_id\t coverage_bin\t exon_coords\t type\t Junc1_coord\t junc2_coord\n'
    outf.write(header)
    for lidx, lsv in enumerate(alllsv[1]):
        for bidx, bn in enumerate(interest):
    #            if lsv.id in oldones: continue
            if lsv.id in bn:
                gene = lsv.id.split(':')[0]
                #typen =correct_type( strand, gene, lsv.type)
                typen = lsv.type
                event = "%s\t%s\t%s-%s\t%s" % (lsv.id, BINNAMES[bidx], lsv.coords[0], lsv.coords[1], typen)

                data = [lsv.id, BINNAMES[bidx], chrom[gene], strand[gene], lsv.coords, lsv.type]

                for junc in lsv.junction_id:
                    event += "\t%s" % (junc.split(':')[1])
                    data.append(junc.split(':')[1])

                print data
                outf.write('%s\n' % event)
                result[BINNAMES[bidx]].append(data)

    outf.close()
    pickle.dump(result, open('./events_for_mat.pkl', 'wb'))
