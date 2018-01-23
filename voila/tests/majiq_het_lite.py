import multiprocessing
import os

import numpy as np
import scipy.stats

from voila.api import Voila

NPROCS = 4


def build_psi_dic(voila_psi_file_list, voila_suffix='.psi_psi.tsv', voila_dir_path='./'):
    samples = []
    dic = {}
    all_juncs_count = 0
    for sample in voila_psi_file_list:
        if sample in samples:
            raise RuntimeError('ERROR! samples should have unique names...%s was already processed' % sample)
        print('processing PSI values for sample %s...' % sample)
        samples.append(sample)
        fd = open(voila_dir_path + sample + voila_suffix, 'r')
        for line in fd:
            if line.startswith('#'):
                continue
            l = line.strip().split('\t')
            ID = l[2]  # lsv id
            psis = [float(x) for x in l[3].split(';')]  # expected psi list
            juncs = l[14].split(';')  # junction coords list
            try:
                dic[ID]
            except:
                dic[ID] = {}
                dic[ID]['juncs'] = juncs
                num_juncs = len(juncs)
                all_juncs_count += num_juncs
            dic[ID][sample] = psis
        fd.close()
    print('finished all files. %s unique LSVs found' % len(dic.keys()))
    print('\tand from those %s junctions found in 1 or more sample' % all_juncs_count)
    return dic, samples, all_juncs_count


def get_lsv_data(voila_file, sample, lsv_ids):
    dic = {}
    with Voila(voila_file, 'r') as v:
        for gene_id, lsv_id in lsv_ids:
            lsv = v.get_voila_lsv(gene_id, lsv_id)
            ID = lsv_id
            juncs = tuple(j.coords() for j in lsv.junctions)
            psis = lsv.means
            dic[ID] = {sample: psis, 'juncs': juncs}

        return dic


def chunkify(lst, n):
    lst = tuple(lst)
    for i in range(n):
        yield lst[i::n]


def build_psi_dic_voila_file(voila_psi_file_list, voila_suffix, voila_dir_path, pool):
    samples = []
    dic = {}
    all_juncs_count = 0
    for sample in voila_psi_file_list:

        if sample in samples:
            raise RuntimeError('ERROR! samples should have unique names...%s was already processed' % sample)
        print('processing PSI values for sample %s...' % sample)
        samples.append(sample)

        voila_file = os.path.join(voila_dir_path, '{0}.{1}'.format(sample, voila_suffix))
        with Voila(voila_file, 'r') as v:
            lsvs = tuple(v.get_lsvs())

        multiple_results = [pool.apply_async(get_lsv_data, (voila_file, sample, x)) for x in chunkify(lsvs, NPROCS)]

        for res in multiple_results:
            for ID, value in res.get().items():
                for sample, psis in value.items():
                    try:
                        dic[ID][sample] = psis
                    except KeyError:
                        dic[ID] = value

        all_juncs_count = sum(len(value['juncs']) for key, value in dic.items())

        get_lsv_data(voila_file, sample, lsvs)

        print('finished all files. %s unique LSVs found' % len(dic.keys()))
        print('\tand from those %s junctions found in 1 or more sample' % all_juncs_count)

    return dic, samples, all_juncs_count


def create_psi_matrix(psi_dic, samples, all_juncs, outname_file='PSI_summary_table.tsv', write_output_matrix=False):
    # psi_dic, samples, all_juncs =  build_psi_dic(['AACB482_alpha','ABDG032_alpha'], dir_path = './psi_tsvs/')
    all_quant = 0
    junc_IDs = []
    psi_data = np.empty(
        (len(samples), all_juncs))  # rows are the samples and columns reprsent each junction (matrix[:,j] slice)
    psi_data.fill(np.nan)
    junc_idx = 0
    for lsv in psi_dic.keys():
        if len(psi_dic[lsv]) + 1 == len(samples):
            all_quant += 1
        for j in range(len(psi_dic[lsv]['juncs'])):
            junc = psi_dic[lsv]['juncs'][j]
            junc_IDs.append('%s_%s' % (lsv, junc))
            for s in range(len(samples)):
                samp = samples[s]
                try:
                    psi = psi_dic[lsv][samp][j]
                    psi_data[s][junc_idx] = psi
                except:
                    continue
            junc_idx += 1
            # return psi_data[:,0], samples, junc_IDs[0]
    print('Number of LSVs found in all %s samples: %s' % (len(samples), all_quant))
    print('shape of PSI array is ')
    print(np.shape(psi_data))
    print('Finished creating PSI matrix')
    return psi_data, junc_IDs


def reduce_psi_matrix(name_list1, name_list2, column_names, psi_matrix):
    rows1 = [n for n in range(len(column_names)) if column_names[n] in name_list1]
    rows2 = [n for n in range(len(column_names)) if column_names[n] in name_list2]
    set1_data = psi_matrix[rows1]
    set2_data = psi_matrix[rows2]

    return set1_data, set2_data


def run_dpsi_test(set1_data, set2_data, junc_IDs, out_prefix, stat='RankSum', dPSI_metric='median', dPSI_thresh=0.1,
                  pval_thresh=0.05):
    fw = open(out_prefix + '_%sdPSI_%s_%s_pval_%s.txt' % (dPSI_metric, int(dPSI_thresh * 100), stat, pval_thresh), 'w')
    fw.write('#LSV_junc\tpval\tdPSI\tSet1_PSIs\tSet2_PSIs\tNumQuant_Set1\tNumQuant_Set2\n')
    sig_juncs = 0
    sig_LSVs = set([])

    for n in range(len(junc_IDs)):
        junc = junc_IDs[n]
        psis1 = set1_data[:, n]
        psis2 = set2_data[:, n]
        non_NaN_psis1 = psis1[~np.isnan(psis1)]
        non_NaN_psis2 = psis2[~np.isnan(psis2)]

        if len(non_NaN_psis1) == 0:
            continue
        if len(non_NaN_psis2) == 0:
            continue
        if stat == 'RankSum':
            test_stat, pval = scipy.stats.ranksums(non_NaN_psis1, non_NaN_psis2)
        if pval > pval_thresh:
            continue

        if dPSI_metric == 'median':
            PSI1 = np.median(non_NaN_psis1)
            PSI2 = np.median(non_NaN_psis2)
        dPSI = PSI2 - PSI1
        if abs(dPSI) < dPSI_thresh:
            continue
        LSV = junc.split('_')[0]
        sig_LSVs.add(LSV)
        sig_juncs += 1
        fw.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
            junc, pval, dPSI, ';'.join([str(x) for x in non_NaN_psis1]), ';'.join([str(x) for x in non_NaN_psis2]),
            len(non_NaN_psis1), len(non_NaN_psis2)))
    fw.close()
    print('Finished analysis...\tFound %s changing junctions in %s unique LSVs' % (sig_juncs, len(sig_LSVs)))


def run_all(name_list1, name_list2, out_prefix, stat='RankSum', dPSI_metric='median', dPSI_thresh=0.1, pval_thresh=0.05,
            voila_dir_path='./', voila_suffix='.psi_psi.tsv', is_file=False, hdf5=False, pool=None):
    '''
    name_list1 and 2: List of strings containing Voila PSI file prefixes for group1 and group2 you want to compare
        if is_file = True [default_view: False], then name_list1 and name_list2 should be text files containing the file names
        seperated by newlines
    out_prefix: For naming the results file
    stat: Which statistical test to use [Currently only Rank Sum]
    dPSI_metric: Which "dPSI" comparison to use [Currently only median differences]
    dPSI_thresh: Minimum dPSI for a junction to be reported [Default: 0.1]
    pval_thresh: Maximum p-value for a junction to be reported [Default: 0.05]
    voila_dir_path: Path where all the Voila PSI text files are located
    voila_suffix: Voila PSI text file suffix that should be shared by all files [Default: 'psi_psi.tsv']       
    '''

    if is_file:
        name_list1 = open(name_list1, 'r').read().splitlines()
        name_list2 = open(name_list2, 'r').read().splitlines()
    set1 = set(name_list1)
    set2 = set(name_list2)
    if not len(set1 & set2) == 0:
        raise RuntimeError(
            'ERROR! The list of file names 1 and 2 contain overlaping samples (%s)...they should be unique' % (
                set1 & set2))
    all_files = name_list1 + name_list2

    if hdf5:
        psi_dic, samples, all_juncs_count = build_psi_dic_voila_file(all_files, voila_dir_path=voila_dir_path,
                                                                     voila_suffix=voila_suffix, pool=pool)
    else:
        psi_dic, samples, all_juncs_count = build_psi_dic(all_files, voila_dir_path=voila_dir_path,
                                                          voila_suffix=voila_suffix)

    psi_matrix, junc_IDs = create_psi_matrix(psi_dic, samples, all_juncs_count)
    data1, data2 = reduce_psi_matrix(name_list1, name_list2, samples, psi_matrix)
    run_dpsi_test(data1, data2, junc_IDs, out_prefix, stat=stat, dPSI_metric=dPSI_metric, dPSI_thresh=dPSI_thresh,
                  pval_thresh=pval_thresh)


if __name__ == "__main__":
    pool = multiprocessing.Pool(processes=NPROCS)
    run_all(['SRR1782688'], ['SRR1782690'], 'IBD', hdf5=True,
            voila_dir_path='/Users/cjgreen/Downloads/majiq_psi',
            voila_suffix='psi.voila', pool=pool)
