import numpy as np
from majiq.src.constants import *
import h5py
from majiq.src.sample import sample_from_junctions

class InvalidLSV(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


class LSV():

    def __init__(self, gene_id, gene_strand, ex, ss):

        self.junctions = []
        if ss:
            jncs = list(ex.ob)
        else:
            jncs = list(ex.ib)

        self.type = self.set_type(jncs, ex, gene_strand, ss)
        self.id = "%s:%s:%s-%s" % (gene_id, self.type[0], ex.start, ex.end)
        self.junctions.sort(key=lambda jj: (jj.start, jj.end))
        if len(self.junctions) < 2:
            raise InvalidLSV("not enougth junctions")

    def to_hdf5(self, hdf5grp, lsv_idx):
        njunc = len(self.junctions)
        h_lsv = hdf5grp.create_group("LSVs/%s" % self.id)
        h_lsv.attrs['id'] = self.id
        h_lsv.attrs['type'] = self.type
        h_lsv.attrs['coverage'] = [lsv_idx, lsv_idx + njunc]

        vh_lsv = h_lsv.create_group('visual')
        # self.get_visual_lsv(self.junctions, exp_idx).to_hdf5(vh_lsv)
        return lsv_idx + njunc

    def sample_lsvs(self, exp_idx, junc_mtrx, fitfunc_r, majiq_config):

        ex_index = sorted([xx.index for xx in self.junctions])
        cover = junc_mtrx[exp_idx, ex_index]

        s_lsv = sample_from_junctions(junction_list=cover,
                                      m=majiq_config.m,
                                      k=majiq_config.k,
                                      fitted_one_over_r=fitfunc_r,
                                      debug=majiq_config.debug)

        lsv_trs = np.array([cover.sum(axis=1), np.count_nonzero(cover, axis=1)]).T

        return s_lsv, lsv_trs

    def set_type(self, jlist, ref_exon, gene_strand, ss):

        # rev = (gene_strand == '-')
        sp_list = []
        ref_exon_id = "%s-%s" %(ref_exon.start, ref_exon.end)
        for junc in jlist:
            if ss:
                ex = junc.acceptor
                coord = junc.end
                ref_coord = junc.start
            else:
                ex = junc.donor
                coord = junc.start
                ref_coord = junc.end
            sp_list.append((ex.intron, coord, ref_coord, ex, junc))

        if gene_strand == '-':
            sp_list.sort(key=lambda xx: (-xx[0], xx[1]), reverse=True)
        else:
            sp_list.sort(key=lambda xx: (xx[0], xx[1]), reverse=False)

        ext_type = 's' if (ss and gene_strand != '-') or (not ss and gene_strand == '-') else 't'

        prev_ex = "%s-%s" % (sp_list[0][3].start, sp_list[0][3].end)
        excount = 1

        ref_ss = [xx.start for xx in ref_exon.ob] if ss else [xx.end for xx in ref_exon.ib]

        for intron, coord, ref_coord, ex, junc in sp_list:

            exid = "%s-%s" % (ex.start, ex.end)

            jidx = ref_ss.index(ref_coord)
            jidx += 1

            if exid == ref_exon_id:
                continue

            if intron:
                ext_type += "|i"
                continue

            if ex is None:
                ext_type += "|%se0" % jidx
                continue

            if ss:
                ss_list = sorted([xx.end for xx in ex.ib])
            else:
                ss_list = sorted([xx.start for xx in ex.ob])

            if prev_ex != exid:
                prev_ex = exid
                excount += 1
            ext_type += "|%se%s.%so%s" % (jidx, excount, ss_list.index(coord)+1, len(ss_list))

            self.junctions.append(junc)
        return ext_type

    @staticmethod
    def junc_cov_to_hdf5(hdf5grp, boots, cov):

        njunc = boots.shape[0]
        lsv_idx = hdf5grp.attrs['lsv_idx']
        if lsv_idx + njunc > 2:
            shp = hdf5grp[JUNCTIONS_DATASET_NAME].shape
            shp_new = shp[0] + 5000
            hdf5grp[JUNCTIONS_DATASET_NAME].resize((shp_new, shp[1]))
            hdf5grp['junc_cov'].resize((shp_new, 2))

        hdf5grp[JUNCTIONS_DATASET_NAME][lsv_idx:lsv_idx + njunc] = boots
        hdf5grp['junc_cov'][lsv_idx:lsv_idx + njunc] = cov

        hdf5grp.attrs['lsv_idx'] = lsv_idx + njunc


def detect_lsvs(list_exons, junc_mtrx, fitfunc_r, locks, gid, gstrand, majiq_config):

    count = 0

    sum_trx = junc_mtrx.sum(axis=2)
    pos_trx = np.count_nonzero(junc_mtrx, axis=2)
    for name, ind_list in majiq_config.tissue_repl.items():

        lsv_list = [[], []]

        for ex in list_exons:
            for ii, jjlist in enumerate([ex.ib, ex.ob]):
                if len(jjlist) < 2:
                    continue
                ex_index = sorted([xx.index for xx in jjlist])
                ex_mtrx_s = sum_trx[ind_list][:, ex_index]
                ex_mtrx_p = pos_trx[ind_list][:, ex_index]

                if np.any(np.logical_and(ex_mtrx_s > majiq_config.minpos, ex_mtrx_p > majiq_config.minreads)):
                    try:
                        lsv_list[ii].append(LSV(gid, gstrand, ex, ss=(ii == 1)))
                    except InvalidLSV:
                        continue

        inlist = []
        for ss in lsv_list[0]:
            for st in lsv_list[1]:
                if set(ss.junctions).issubset(set(st.junctions)) and not set(ss.junctions).issuperset(set(st.junctions)):
                    break
            else:
                inlist.append(ss)
        for st in lsv_list[1]:
            for ss in lsv_list[0]:
                if set(st.junctions).issubset(set(ss.junctions)):
                    break
            else:
                inlist.append(st)

        for exp_idx in majiq_config.tissue_repl[name]:
            np_jjlist = []
            attrs_list = []

            for lsvobj in inlist:
                b, c = lsvobj.sample_lsvs(exp_idx, junc_mtrx, fitfunc_r=fitfunc_r[exp_idx],
                                          majiq_config=majiq_config)
                np_jjlist.append(b)
                attrs_list.append(c)

            if len(np_jjlist) > 0:
                mtrx = np.concatenate(np_jjlist, axis=0)
                mtrx_attrs = np.concatenate(attrs_list, axis=0)

                locks[exp_idx].acquire()
                with h5py.File('%s/%s.majiq' % (majiq_config.outDir, majiq_config.sam_list[exp_idx]), 'r+') as f:
                    lsv_idx = f.attrs['lsv_idx']
                    for lsv in inlist:
                        lsv_idx = lsv.to_hdf5(f, lsv_idx)
                    LSV.junc_cov_to_hdf5(f, mtrx, mtrx_attrs)
                    f.attrs['num_lsvs'] = f.attrs['num_lsvs'] + len(inlist)
                locks[exp_idx].release()
                count += len(inlist)

    return count
