from majiq.src.constants import *
import h5py
import numpy as np
import os
import pickle
from majiq.src.config import Config

#######
# HDF5 API
#######


def _prepare_and_dump(filename, logging=None):
    majiq_config = Config()
    if logging is not None:
        logging.debug("Number of Genes in DB", len(majiq_config.gene_tlb))
    db_f = h5py.File(filename, 'w')
    for gidx, gn in enumerate(majiq_config.gene_tlb.values()):
        gn.to_hdf5(db_f)


def get_const_junctions(filename, logging=None):
    if not os.path.exists(filename):
        if logging is not None:
            logging.error('File % doesn\'t exists' % filename)
            raise UserWarning
    else:
        db_f = h5py.File(filename, 'r')
        cc = db_f[JUNCTIONS_DATASET_NAME][()]
        db_f.close()
        return np.array(cc)

def load_data_lsv(path, group_name, logger=None):
    """Load data from the preprocess step. Could change to a DDBB someday"""
    data = h5py.File(path, 'r')
    lsv_cov_list = []
    lsv_info = []

    meta_info = dict()
    meta_info['group'] = group_name
    meta_info['sample_id'] = data.attrs['sample_id']
    meta_info['fitfunc'] = data.attrs['fitfunc']

    try:
        for lsvid in data['LSVs'].keys():
            lsv = data['LSVs/%s' % lsvid]
            lsv_info.append([lsv.attrs['id'], lsv.attrs['type'], lsv['visual']])
            lsv_cov_list.append(data[JUNCTIONS_DATASET_NAME][lsv.attrs['coverage']])
            sh = data[JUNCTIONS_DATASET_NAME][lsv.attrs['coverage']].shape
            if sh[0] < 2:
                print("WRONG LSV %s" % lsvid)

    except KeyError:
        logger.info("No LSVs in file")
        raise

    return meta_info, [lsv_cov_list, lsv_info]


def load_lsvgraphic_from_majiq(h5df_grp, lsv_id):
    try:
        return h5df_grp['/LSVs/%s/visual' % lsv_id]
    except KeyError:
        return None


def read_meta_info(list_of_files):
    meta = {'experiments': []}
    for fl in list_of_files:
        with h5py.File(fl, 'r') as fp :
            meta['experiments'].append(fp.attrs['sample_id'])
            try:
                if meta['genome'] != fp.attrs['genome']:
                    raise RuntimeError('Combining experiments from different genome assemblies. Exiting')
            except KeyError:
                meta['genome'] = fp.attrs['genome']
                continue

    return meta


def open_hdf5_file(filename, **kwargs):
    return h5py.File(filename, 'r')


def close_hdf5_file(fp):
    return fp.close()



# bootstrap files

from majiq.grimoire.lsv import quant_lsv


def get_extract_lsv_list(list_of_lsv_id, file_list):
    result = []

    for lsv_id in list_of_lsv_id:
        lsv_cov = []
        lsv_type = None
        for fidx, fname in enumerate(file_list):

            with open_hdf5_file(fname) as data:
                try:
                    if lsv_type is None:
                        lsv_type = data['LSVs/%s' % lsv_id].attrs['type']

                    assert data['LSVs/%s' % lsv_id].attrs['type'] == lsv_type, "ERROR lsv_type doesn't match for %s" % lsv_id
                    cov = data['LSVs/%s' % lsv_id].attrs['coverage']
                    lsv_cov.append(data[JUNCTIONS_DATASET_NAME][cov[0]:cov[1]])
                except KeyError:
                    lsv_cov.append(None)

#        lsv_cov = np.array(lsv_cov)
        qq = quant_lsv(lsv_id, lsv_type, lsv_cov)
        result.append(qq)
    return result


def add_lsv_to_bootstrapfile_with_lock(lsv_id, lsv_type, samples, num_exp, lock_per_file, outdir, name):

    for ii in range(num_exp):
        vals = {'samples': samples[ii], 'id': lsv_id, 'type': lsv_type}
        file_name = '%s/%s.%d.boots.hdf5' % (outdir, name, ii)
        lock_per_file[ii].acquire()
        with h5py.File(file_name, 'r+') as f:
            lsv_idx = f.attrs['lsv_idx']
            lsv_idx = boots_write(f, vals, lsv_idx)
            f.attrs['lsv_idx'] = lsv_idx
        lock_per_file[ii].release()


def add_lsv_to_bootstrapfile(f, vals):
    lsv_idx = f.attrs['lsv_idx']
    lsv_idx = boots_write(f, vals, lsv_idx)
    f.attrs['lsv_idx'] = lsv_idx
    f.attrs['num_lsvs'] += 1


def create_bootstrap_file(file_list, outdir, name, m=100):
    import datetime
    for ii, ff in enumerate(file_list):
        f = h5py.File('%s/%s.%d.boots.hdf5' % (outdir, name, ii), 'w')
        f.create_dataset('junctions', (5000, m), maxshape=(None, m))
        # fill meta info
        f.attrs['sample_id'] = ff
        f.attrs['date'] = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
        f.attrs['VERSION'] = VERSION
        f.attrs['lsv_idx'] = 0
        f.close()


def load_bootstrap_samples(lsv_id, file_list, weight=True):
    lsv_samples = []
    lsv_type = None
    for ii, f in enumerate(file_list):
        samples = f['junctions'][f['LSVs/%s' % lsv_id].attrs['coverage']]
        lsv_type = f['LSVs/%s' % lsv_id].attrs['type']
        if weight:
            #wght = f['weights'][f["LSVs/%s" % lsv_id].attrs['weight_idx']]
            wght = f["LSVs/%s" % lsv_id].attrs['weight']
            samples *= wght
        lsv_samples.append(samples)
    return lsv_samples, lsv_type


def store_weights_bootstrap(lsv_list, wgts, file_list, outdir, name):
    for ii, ff in enumerate(file_list):
        file_name = '%s/%s.%d.boots.hdf5' % (outdir, name, ii)
        with h5py.File(ff, 'r+') as f:

            for idx, lsv in lsv_list.items():
                f["LSVs/%s" % lsv.id].attrs['weight'] = wgts[idx, ii]
            f.close()


def open_bootstrap_samples(num_exp, directory, name):
    result = []
    for ii in range(num_exp):
        file_name = '%s/%s.%d.boots.hdf5' % (directory, name, ii)
        result.append(h5py.File(file_name, 'r+'))
    return result


def close_bootstrap_file(file_list, outdir, name, m=100):
    for ii, ff in enumerate(file_list):
        file_name = '%s/%s.%d.boots.hdf5' % (outdir, name, ii)
        f = h5py.File(file_name, 'r+')
        lsv_idx = f.attrs['lsv_idx']
        f['junctions'].resize((lsv_idx, m))
        f.close()


def boots_write(hg_grp, vals, lsv_idx, dpsi=False):

    njunc = vals['samples'].shape[0]
    if lsv_idx + njunc > 2:
        shp = hg_grp['junctions'].shape
        shp_new = shp[0] + 5000
        hg_grp['junctions'].resize((shp_new, shp[1]))
        hg_grp['junc_cov'].resize((shp_new, 2))

    hg_grp['junctions'][lsv_idx:lsv_idx+njunc] = vals['samples']
    #print('%s::%s' % (vals['junc_attr'], type(vals['junc_attr'])))
    hg_grp['junc_cov'][lsv_idx:lsv_idx + njunc] = vals['junc_attr']

    h_lsv = hg_grp.create_group("LSVs/%s" % vals['id'])
    h_lsv.attrs['id'] = vals['id']
    h_lsv.attrs['type'] = vals['type']
    h_lsv.attrs['coverage'] = [lsv_idx, lsv_idx + njunc]
    #TODO: CHECK
    vh_lsv = h_lsv.create_group('visual')
    vals['lsv_graphic'].to_hdf5(vh_lsv)

    return lsv_idx + njunc