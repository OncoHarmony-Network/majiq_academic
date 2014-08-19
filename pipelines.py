import os
from collections import defaultdict
import abc
import pickle
from multiprocessing import Pool, Manager, current_process

from pylab import *
from numpy.ma import masked_less
from scipy.stats import norm
from scipy.stats import scoreatpercentile
from grimoire.utils.utils import create_if_not_exists, get_logger

from analysis.polyfitnb import fit_nb 
from analysis.sample import sample_from_junctions, mean_junction
from analysis.psi import calc_psi, mean_psi, simple_psi, DirichletCalc, reads_given_psi, BINS_CENTER, lsv_psi
from analysis.adjustdelta import adjustdelta
from analysis.weight import local_weight_eta_nu, global_weight_ro
from analysis.matrix import rank_deltas_lsv, collapse_matrix, rank_empirical_delta_lsv
import analysis.filter as majiq_filter
import analysis.io as majiq_io
import analysis.psi as  majiq_psi
import analysis.sample as majiq_sample
#from analysis.filter import norm_junctions, discardlow, discardhigh, discardminreads, discardmaxreads, discardminreads_and, mark_stacks
import pipe as pipe
################################
# Data loading and Boilerplate #
################################




def get_clean_raw_reads( matched_info, matched_lsv, outdir, names, num_exp ):

    res = []

    import ipdb
    ipdb.set_trace()
    for eidx in xrange(num_exp):
        for ldx, lsv in enumerate(matched_info):
            num = matched_lsv[eidx][ldx].sum()
            res.append([lsv[1],num])

        with open('%s/clean_reads.%s%d.pkl'%(outdir, names,eidx),'wb') as fp:
            pickle.dump(res, fp)




def _pipeline_run(pipeline, lsv=False, logger=None):
    "Exception catching for all the pipelines"
    try:
        return pipeline.run(lsv)

    except KeyboardInterrupt:
        if pipeline.logger: pipeline.logger.info("MAJIQ manually interrupted. Avada kedavra...")

# PLOTTING STUFF. All this should eventually go to a plotting module 
def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        savefig("%s%s.png"%(plotpath, plotname), bbox_inches='tight') 
        clf()
    else:
        show()


def preprocess(args):
    raise NotImplemented

class BasicPipeline:
    def __init__(self, args):
        """Basic configuration shared by all pipelines"""
        #trick to dump argparse arguments into self
        self.__dict__.update(args.__dict__)
        create_if_not_exists(self.output)
        if self.plotpath:
            create_if_not_exists(self.plotpath)
        #grab logger
        logger_path = self.logger
        if not logger_path:
            logger_path = self.output

        self.logger = get_logger("%smajiq.log"%logger_path, silent=self.silent, debug=self.debug)
        self.lsv = args.lsv
        self.nthreads = args.nthreads
        self.nz = args.nz
        self.psi_paths = []
        try:
            self.replica_len = [len(self.files1), len(self.files2)]
        except AttributeError:
            pass

    @abc.abstractmethod
    def run(self, lsv):
        """This is the entry point for all pipelines"""
        return

    def gc_content_norm_lsv( self, lsv_list, const_list ) :
        "Normalize the matrix using the gc content"
        self.logger.info("GC content normalization...")
        if self.gcnorm:
            for lidx, lsv in enumerate(lsv_list[0]):
                lsv = lsv * lsv_list[2][lidx]
            conts_list[0] = const_list[0] * const_list[2]
        return lsv_list, const_list

    def gc_content_norm(self, all_junctions):
        "Normalize the matrix using the gc content"
        self.logger.info("GC content normalization...")
        if self.gcnorm:
            for junc_set in all_junctions.keys():
                all_junctions[junc_set] = majiq_filter.norm_junctions(all_junctions[junc_set]["junctions"], all_junctions[junc_set]["gc_content"])

        return all_junctions

    def fitfunc(self, const_junctions):
        "Calculate the Negative Binomial function to sample from using the Constitutive events"
        if self.debug:
            self.logger.debug("Skipping fitfunc because --debug!")
            return poly1d([1, 0])
        else:
            self.logger.info("Fitting NB function with constitutive events...")
            return fit_nb(const_junctions, "%s_nbfit"%self.output, self.plotpath, nbdisp=self.nbdisp, logger=self.logger, discardb=True, bval=True)

    def mark_stacks_lsv(self, lsv_list, fitfunc):
        if self.markstacks >= 0:
            self.logger.info("Marking and masking stacks for...")
            lsv_list = majiq_filter.lsv_mark_stacks(lsv_list, fitfunc, self.markstacks, self.nbdisp, self.logger)

        return lsv_list

    def mark_stacks(self, all_junctions, fitfunc):
        if self.markstacks >= 0:
            self.logger.info("Marking and masking stacks for...")
            for junc_set in all_junctions.keys():
                if junc_set.find("const") == -1:
                    self.logger.info("... %s"%junc_set)
                    all_junctions[junc_set] = majiq_filter.mark_stacks(all_junctions[junc_set], fitfunc, self.markstacks, self.nbdisp, self.logger)

                all_junctions[junc_set] = masked_less(all_junctions[junc_set], 0) #remask the stacks

        return all_junctions






def parallel_calcpsi ( conf, lsv_junc, fitfunc, name, chunk, tempdir):
    __parallel_calcpsi_lsv( conf, lsv_junc, fitfunc, name, chunk, tempdir)
    return "FIN"


################################
# PSI calculation pipeline     #
################################

def calcpsi(args):
    return _pipeline_run(CalcPsi(args), args.lsv)

def __parallel_calcpsi_lsv( conf , lsv_junc, fitfunc, name, chunk, tempdir):

    try:
        if not os.path.isdir(tempdir):
            os.mkdir(tempdir)
        thread_logger = get_logger("%s/majiq.w%s.log"%(tempdir,chunk), silent=False, debug=conf['debug'])
        thread_logger.info( "[Th %s]: START child,%s"%(chunk,current_process().name))
        thread_logger.info('[Th %s]: Filtering ...'%(chunk))


        num_lsv = len(lsv_junc[0])
        ''' 
            fon[0] = False deactivates num_reads >= 10
            fon[1] = False deactivates npos >= 5
        '''
        fon = [True, True]
        lsv_junc = majiq_filter.lsv_quantifiable( lsv_junc, conf['minnonzero'], conf['minreads'], thread_logger , fon)
        thread_logger.info('[Th %s]: %s/%s lsv remaining'%(chunk, len(lsv_junc[0]), num_lsv))

        thread_logger.info("[Th %s]: Bootstrapping samples..."%(chunk)) 
        lsv_sample = []
        for ii in lsv_junc[0]:

            m_lsv, var_lsv, s_lsv = sample_from_junctions(  ii,
                                                            conf['m'],
                                                            conf['k'],
                                                            discardzeros= conf['discardzeros'],
                                                            trimborder  = conf['trimborder'],
                                                            fitted_func = fitfunc,
                                                            debug       = conf['debug'],
                                                            Nz          = conf['Nz'])
            lsv_sample.append( s_lsv )

        thread_logger.info("[Th %s]: Calculating PSI for %s ..."%(chunk, name))
        psi = lsv_psi(lsv_sample, conf['alpha'], conf['n'], conf['debug'])

        thread_logger.info("[Th %s]: Saving PSI..."%chunk)
        output = open("%s/%s_th%s.psi.pickle"%(tempdir, name, chunk), 'w')
        pickle.dump((psi, lsv_junc[1]), output)
        thread_logger.info("[Th %s]: PSI calculation for %s ended succesfully! Result can be found at %s"%(chunk, name, output.name))
    except Exception as e:
        print "%s"%sys.exc_traceback.tb_lineno, e
        sys.stdout.flush()

class CalcPsi(BasicPipeline):

    def run(self, lsv=False):
        ret = []
        if self.newmodel:
            ret = self.calcpsi_newmodel()
        else:
            ret = self.calcpsi_lsv(path)
        return ret


    def calcpsi_newmodel(self):
        """
        Given a file path with the junctions, return psi distributions. 
        write_pickle indicates if a .pickle should be saved in disk
        """
        self.logger.info("")
        self.logger.info("Running PSI groups...")
        self.logger.info("GROUP : %s"%self.files)

        num_exp = len(self.files)

        filtered_lsv = [None]*num_exp
        fitfunc = [None]*num_exp
        for ii, file in enumerate(self.files):
            lsv_junc, const = majiq_io.load_data_lsv(file, self.logger) 

            #fitting the function
#            lsv_junc = self.gc_content_norm_lsv( lsv_junc, const )
            fitfunc[ii] = self.fitfunc(const[0])
            filtered_lsv[ii] = self.mark_stacks_lsv( lsv_junc, fitfunc[ii])
        matched_lsv, matched_info = majiq_filter.quantifiable_in_group( filtered_lsv, self.minnonzero, self.minreads, self.logger , 0.10 )

        conf = { 'minnonzero':self.minnonzero,
                 'minreads': self.minreads,
                 'm':self.m,
                 'k':self.k,
                 'discardzeros':self.discardzeros,
                 'trimborder':self.trimborder,
                 'debug':self.debug,
                 'alpha':self.alpha,
                 'n':self.n, 
                 'nbins':40,
                 'nz':self.nz}

        if self.nthreads == 1:
            posterior_matrix, names = pipe.calcpsi( matched_lsv, matched_info, num_exp, conf, fitfunc, self.logger)
        else:
            pool = Pool(processes=self.nthreads)
            csize = len(matched_lsv) / int(self.nthreads)
            self.logger.info("CREATING THREADS %s with <= %s lsv"%(self.nthreads,csize))
            jobs = []

            for nt in xrange(self.nthreads):
                lb = nt * csize
                ub = min( (nt+1) * csize, len(matched_lsv) )
                if nt == self.nthreads - 1 : ub = len(matched_lsv) 
                lsv_list = matched_lsv[lb:ub]
                lsv_info = matched_info[lb:ub]
                print nt, ub, lb
                pool.apply_async( pipe.parallel_lsv_child_calculation, [ pipe.calcpsi,
                                                                         [ lsv_list, lsv_info, num_exp, conf, fitfunc ],
                                                                         lsv_info,
                                                                         '%s/tmp'%os.path.dirname(self.output),
                                                                         self.name,
                                                                         nt] )
            pool.close()
            pool.join()

            posterior_matrix = []
            names = []
            self.logger.info("GATHER pickles")
            for nt in xrange(self.nthreads):
                tempfile = open("%s/tmp/%s_th%s.calcpsi.pickle"%(os.path.dirname(self.output), self.name, nt))
                ptempt = pickle.load( tempfile )
                posterior_matrix.extend( ptempt[0] )
                names.extend(ptempt[1])


        pickle_path = "%s/%s_psigroup.pickle"%(self.output, self.name)
        pickle.dump([posterior_matrix, names], open(pickle_path, 'w'))
        self.logger.info("PSI calculation for %s ended succesfully! Result can be found at %s"%(self.name, self.output))
        self.logger.info("Alakazam! Done.")


    def calcpsi_lsv(self, path, write_pickle=True):
        """
        Given a file path with the junctions, return psi distributions. 
        write_pickle indicates if a .pickle should be saved in disk
        """
        name = os.path.basename(path)
        self.logger.info("")
        self.logger.info("Loading %s..."%path)
        lsv_junc, const = majiq_io.load_data_lsv(path, self.logger) 
        self.logger.debug("SHAPES for lsv %s,  constitutive %s"%(len(lsv_junc[0]), const[0].shape))
        self.logger.info("Loaded.")

        all_junctions = self.gc_content_norm_lsv( lsv_junc, const )

        fitfunc = self.fitfunc(const[0])
        filter_lsv = self.mark_stacks_lsv( lsv_junc, fitfunc)

        conf = { 'minnonzero':self.minnonzero,
                 'minreads': self.minreads,
                 'm':self.m,
                 'k':self.k,
                 'discardzeros':self.discardzeros,
                 'trimborder':self.trimborder,
                 'debug':self.debug,
                 'alpha':self.alpha,
                 'n':self.n, 
                 'Nz':self.nz}

        if self.nthreads == 1: 
            parallel_calcpsi( conf, filter_lsv, fitfunc, name, 0, '%s/tmp'%os.path.dirname(self.output))
            tempfile = open("%s/tmp/%s_th0.psi.pickle"%(self.output, name))
            ptempt = pickle.load( tempfile )
            psi = ptempt[0] 
            info =  ptempt[1] 
        else:
            try:
                pool = Pool(processes=self.nthreads)
                csize = len(filter_lsv[0]) / int(self.nthreads)
                self.logger.info("CREATING THREADS %s"%self.nthreads)
                jobs = []
    
                for nt in xrange(self.nthreads):
                    lb = nt * csize
                    ub = min( (nt+1) * csize, len(filter_lsv[0]) )
                    lsv_list = [filter_lsv[0][lb:ub],filter_lsv[1][lb:ub]]
                    jobs.append(pool.apply_async( parallel_calcpsi, [conf, lsv_list, fitfunc, name, nt, '%s/tmp'%os.path.dirname(self.output)] ))
                pool.close()
                pool.join()
            except Exception as e:
                print e

            psi = []
            info = []
            self.logger.info("GATHER pickles")
            for nt in xrange(self.nthreads):
                tempfile = open("%s/tmp/%s_th%s.psi.pickle"%(os.path.dirname(self.output), name, nt))
                ptempt = pickle.load( tempfile )
                psi.extend( ptempt[0] )
                info.extend( ptempt[1] )

        self.logger.info("Saving PSI...")
        if write_pickle:
            output = open("%s%s_psi.pickle"%(self.output, name), 'w')
            pickle.dump((psi, info), output)
            self.logger.info("PSI calculation for %s ended succesfully! Result can be found at %s"%(name, output.name))

        if self.debug > 0:
            return psi, info[:self.debug]
        return psi, info



################################
#          Delta PSI           #
################################

def deltapair(args):
    _pipeline_run( DeltaPair(args), args.lsv )




def deltapsi_calc( matched_list, matched_info, fitfunc, conf, chunk, prior_matrix, logr ):
    logr.info("[Th %s]: Bootstrapping for all samples..."%chunk)
    lsv_samples = [[],[]]
    for idx_exp, experiment in enumerate(matched_list):
        for idx, ii in enumerate(experiment):
            m_lsv, var_lsv, s_lsv = sample_from_junctions(  junction_list = ii,
                                                            m = conf['m'],
                                                            k = conf['k'],
                                                            discardzeros= conf['discardzeros'],
                                                            trimborder  = conf['trimborder'],
                                                            fitted_func = fitfunc[idx_exp],
                                                            debug       = conf['debug'],
                                                            Nz          = conf['Nz'])
            lsv_samples[idx_exp].append( s_lsv )

    logr.info("[Th %s]: Calculating P(Data | PSI_i, PSI_j)..."%chunk)
    #P(Data | PSI_i, PSI_j) = P(vector_i | PSI_i) * P(vector_j | PSI_j)
    numbins= 20
    data_given_psi = []
    try:
        for lsv_idx, info in enumerate(matched_info):
            if lsv_idx % 50 == 0:
                print "%s...."%lsv_idx,
                sys.stdout.flush()
            data_given_psi1 = majiq_psi.reads_given_psi_lsv( lsv_samples[0][lsv_idx], conf['psi_space'] )
            data_given_psi2 = majiq_psi.reads_given_psi_lsv( lsv_samples[1][lsv_idx], conf['psi_space'] )
            data_psi = []
            for psi, data1 in enumerate(data_given_psi1) :
            #TODO Tensor product is calculated with scipy.stats.kron. Probably faster, have to make sure I am using it correctly.
                sys.stdout.flush()

                data_psi.append(data_given_psi1[psi].reshape(-1, numbins) * data_given_psi2[psi].reshape(numbins, -1))
                if psi == 0: 
                    majiq_psi.plot_matrix(  data_psi[-1],
                                            "P(Data | PSI 1, PSI 2) Event %s.%s (Psi1: %s Psi2: %s)"%(lsv_idx,psi, sum(data_given_psi1[psi]), sum(data_given_psi2[psi])), 
                                            "datagpsi_%s.%s"%(info[1], psi),
                                            conf['plotpath'] )

            data_given_psi.append(data_psi)
        print 
        sys.stdout.flush()
    except Exception as e:
        print "%s"%sys.exc_traceback.tb_lineno, e
        sys.stdout.flush()


    #Finally, P(PSI_i, PSI_j | Data) equivalent to P(PSI_i, PSI_j)* P(Data | PSI_i, PSI_j) 
    logr.info("[Th %s]: Calculate Posterior Delta Matrices..."%chunk)
    posterior_matrix = []
    for lidx, lsv in enumerate(matched_info) :
        if lsv_idx % 50 == 0: 
            print "%s...."%lsv_idx,
            sys.stdout.flush()
        lsv_psi_matrix = []
        for psi in range(len(data_given_psi[lidx])) :
            pm = (prior_matrix * data_given_psi[lidx][psi])
            psi_mat = (pm / sum(pm))
            lsv_psi_matrix.append( psi_mat )
            if psi == 0: 
                majiq_psi.plot_matrix(  psi_mat,
                                    "Posterior Delta Event %s.%s (Psi1: %s Psi2: %s)"%(lidx,psi, sum(data_given_psi1[psi]), sum(data_given_psi2[psi])),
                                    "posterior_dpsi.%s.%s"%(lsv[1],psi),
                                    conf['plotpath'] )
        posterior_matrix.append(lsv_psi_matrix)

    return posterior_matrix

def delta_calculation(matched_info, lsv_samples, conf, prior_matrix, logr=None, chunk=0):
    numbins= 20
    data_given_psi = []
    for lsv_idx, info in enumerate(matched_info):
        if lsv_idx % 50 == 0:
            print "%s...."%lsv_idx,
            sys.stdout.flush()
        data_given_psi1 = majiq_psi.reads_given_psi_lsv( lsv_samples[0][lsv_idx], conf['psi_space'] )
        data_given_psi2 = majiq_psi.reads_given_psi_lsv( lsv_samples[1][lsv_idx], conf['psi_space'] )
        data_psi = []
        for psi, data1 in enumerate(data_given_psi1) :
        #TODO Tensor product is calculated with scipy.stats.kron. Probably faster, have to make sure I am using it correctly.
            sys.stdout.flush()

            data_psi.append(data_given_psi1[psi].reshape(-1, numbins) * data_given_psi2[psi].reshape(numbins, -1))
            if psi == 0: 
                majiq_psi.plot_matrix(  data_psi[-1],
                                        "P(Data | PSI 1, PSI 2) Event %s.%s (Psi1: %s Psi2: %s)"%(lsv_idx,psi, sum(data_given_psi1[psi]), sum(data_given_psi2[psi])), 
                                        "datagpsi_%s.%s"%(info[1], psi),
                                        conf['plotpath'] )

        data_given_psi.append(data_psi)
    print
    sys.stdout.flush()

    #Finally, P(PSI_i, PSI_j | Data) equivalent to P(PSI_i, PSI_j)* P(Data | PSI_i, PSI_j) 
    logr.info("[Th %s]: Calculate Posterior Delta Matrices..."%chunk)
    posterior_matrix = []
    for lidx, lsv in enumerate(matched_info) :
        if lsv_idx % 50 == 0: 
            print "%s...."%lsv_idx,
            sys.stdout.flush()
        lsv_psi_matrix = []
        for psi in range(len(data_given_psi[lidx])) :
            pm = (prior_matrix * data_given_psi[lidx][psi])
            psi_mat = (pm / sum(pm))
            lsv_psi_matrix.append( psi_mat )
            if psi == 0: 
                majiq_psi.plot_matrix(  psi_mat,
                                    "Posterior Delta Event %s.%s (Psi1: %s Psi2: %s)"%(lidx,psi, sum(data_given_psi1[psi]), sum(data_given_psi2[psi])),
                                    "posterior_dpsi.%s.%s"%(lsv[1],psi),
                                    conf['plotpath'] )
        posterior_matrix.append(lsv_psi_matrix)

    return posterior_matrix

def parallel_delta_calculation ( matched_info, lsv_samples, conf, prior_matrix, tempdir, chunk):

    try:
        if not os.path.isdir(tempdir):
                os.mkdir(tempdir)
        thread_logger = get_logger("%s/majiq.w%s.log"%(tempdir,chunk), silent=False, debug=conf['debug'])
        thread_logger.info( "[Th %s]: START child,%s"%(chunk,current_process().name))
        thread_logger.info('[Th %s]: Filtering ...'%(chunk))
        post_matrix = delta_calculation( matched_info, lsv_samples, conf, prior_matrix, logr=thread_logger, chunk=chunk)
        print "%s/%s_%s_th%s.deltapsi.pickle"%(tempdir, conf['names'][0], conf['names'][1], chunk)
        sys.stdout.flush()
        thread_logger.info("[Th %s]: Saving DeltaPSI..."%chunk)
        output = open("%s/%s_%s_th%s.deltapsi.pickle"%(tempdir, conf['names'][0], conf['names'][1], chunk), 'w')
        pickle.dump((post_matrix, matched_info), output)
        thread_logger.info("[Th %s]: DeltaPSI calculation for XX ended succesfully! Result can be found at %s"%(chunk, output.name))
    except Exception as e:
        print "%s"%sys.exc_traceback.tb_lineno, e
        sys.stdout.flush()

    return

def parallel_delta_psi_wrapper ( matched_list, matched_info, fitfunc, conf, prior_matrix, tempdir, chunk):

    if not os.path.isdir(tempdir):
            os.mkdir(tempdir)
    thread_logger = get_logger("%s/majiq.w%s.log"%(tempdir,chunk), silent=False, debug=conf['debug'])
    thread_logger.info( "[Th %s]: START child,%s"%(chunk,current_process().name))
    thread_logger.info('[Th %s]: Filtering ...'%(chunk))

    post_matrix = deltapsi_calc(matched_list, matched_info, fitfunc, conf, chunk, prior_matrix, thread_logger)
    

    print "%s/%s_%s_th%s.deltapsi.pickle"%(tempdir, conf['names'][0], conf['names'][1], chunk)
    sys.stdout.flush()
    thread_logger.info("[Th %s]: Saving DeltaPSI..."%chunk)
    output = open("%s/%s_%s_th%s.deltapsi.pickle"%(tempdir, conf['names'][0], conf['names'][1], chunk), 'w')
    pickle.dump((post_matrix, matched_info), output)
    thread_logger.info("[Th %s]: DeltaPSI calculation for %s ended succesfully! Result can be found at %s"%(chunk, name, output.name))
    return

def delta_precomputed( info, name1, name2, outDir, chunk=0, logger=None ):
    matrices = None
    matrix_path = "%s/%s_%s_th%s.deltapsi.pickle"%(outDir, name1, name2, chunk)
    print "matrix_path", matrix_path
    if os.path.exists(matrix_path):
        logger.info("%s exists! Loading..."%matrix_path)
        tlb, mats = pickle.load(open(matrix_path))
        matrices = []
        for lidx, lsv in enumerate(info):
            new_index = tlb[lsv[1]]
            matrices.append( mats[ new_index] ) 
        
    return matrices

def store_delta_psi(posterior, info, name1, name2, outDir, chunk=0, thread_logger=None):

    thread_logger.info("[Th %s]: Saving DeltaPSI..."%chunk)
    output = open("%s/%s_%s_th%s.deltapsi.pickle"%(outDir, name1, name2, chunk), 'w')

    tlb = {}
    for lidx, lsv in enumerate(info):
        tlb[lsv[1]] = lidx

    pickle.dump([tlb,posterior], output)
    thread_logger.info("[Th %s]: DeltaPSI calculation for %s|%s ended succesfully! Result can be found at %s"%(chunk, name1, name2, outDir))


class DeltaPair(BasicPipeline):

    def _get_delta_info(self, newdist, norm_space):
        """
        Print out some information on the delta PSI
        """
        delta_limit = 0.2
        for i in xrange(len(norm_space)):
            if norm_space[i] > -delta_limit:
                pos_index_neg = i
                break

        for i in xrange(len(norm_space)-1, 0, -1):
            if norm_space[i] < delta_limit:
                pos_index_pos = i
                break

        self.logger.info("Synthetic prior: Delta PSI > %s (%.4f)"%(delta_limit, newdist[0:pos_index_neg].sum()))
        self.logger.info("Delta PSI < -%s (%.4f)"%(delta_limit, newdist[pos_index_pos+1:].sum()))
        message = "norm(mean=0, std=%s) + uniform Delta PSI > 0.2 = %.2f Delta PSI < -0.2 = %.2f"%(self.priorstd, newdist[pos_index_pos+1:].sum(), newdist[0:pos_index_neg].sum())
        self.logger.info(message)
        title(message)

    def run( self, lsv=False ):
        if lsv :
            self.pairdelta_lsv(self.file1, self.file2, self.output)
        else:
            self.pairdelta(self.file1, self.file2, self.output)

    def pairdelta_lsv(self, file1, file2, output):
        self.logger.info("")
        self.logger.info("Processing pair %s - %s..."%(file1, file2))

        lsv_junc1, const1 = majiq_io.load_data_lsv(file1, self.logger) 
        lsv_junc2, const2 = majiq_io.load_data_lsv(file2, self.logger) 

        #fitting the function
        fitfunc1 = self.fitfunc(const1[0])
        fitfunc2 = self.fitfunc(const2[0])

        filtered_lsv1 = self.mark_stacks_lsv( lsv_junc1, fitfunc1)
        filtered_lsv2 = self.mark_stacks_lsv( lsv_junc2, fitfunc2)



        #Quantifiable junctions filter
        ''' Quantify and unify '''
        self.logger.info('Filtering ...')
        num_lsv1 = len(filtered_lsv1[0])
        num_lsv2 = len(filtered_lsv2[0])

        ''' 
            fon[0] = False deactivates num_reads >= 20
            fon[1] = False deactivates npos >= 5
        '''
        fon = [True, True]
        filtered_lsv1 = majiq_filter.lsv_quantifiable( filtered_lsv1, self.minnonzero, self.minreads, self.logger , fon)
        self.logger.info('%s/%s lsv remaining'%(len(filtered_lsv1[0]),num_lsv1))
        filtered_lsv2 = majiq_filter.lsv_quantifiable( filtered_lsv2, self.minnonzero, self.minreads, self.logger , fon)
        self.logger.info('%s/%s lsv remaining'%(len(filtered_lsv2[0]),num_lsv2))


        psi_space, prior_matrix = majiq_psi.gen_prior_matrix(self, filtered_lsv1, filtered_lsv2, output)

        matched_lsv, matched_info = majiq_filter.lsv_intersection( filtered_lsv1, filtered_lsv2 )

        get_clean_raw_reads( matched_info, matched_lsv, output, self.names )


        conf = { 'minnonzero':self.minnonzero,
                 'minreads': self.minreads,
                 'm':self.m,
                 'k':self.k,
                 'discardzeros':self.discardzeros,
                 'trimborder':self.trimborder,
                 'debug':self.debug,
                 'alpha':self.alpha,
                 'n':self.n, 
                 'plotpath':self.plotpath,
                 'Nz':self.nz,
                 'names':self.names,
                 'psi_space':psi_space}

        if self.nthreads == 1:
            posterior_matrix = deltapsi_calc(matched_lsv, matched_info, [fitfunc1,fitfunc2], conf, 'master' , prior_matrix, self.logger)
        else:
            try:
                pool = Pool(processes=self.nthreads)
                csize = len(matched_lsv[0]) / int(self.nthreads)
                self.logger.info("CREATING THREADS %s with <= %s lsv"%(self.nthreads,csize))
                jobs = []
    
                for nt in xrange(self.nthreads):
                    lb = nt * csize
                    ub = min( (nt+1) * csize, len(matched_lsv[0]) )
                    if nt == self.nthreads - 1 : ub = len(matched_lsv[0]) 
                    lsv_list = [matched_lsv[0][lb:ub],matched_lsv[1][lb:ub]]
                    lsv_info = matched_info[lb:ub]
                    jobs.append(pool.apply_async( parallel_delta_psi_wrapper, [ lsv_list,
                                                                                lsv_info,
                                                                                [fitfunc1,fitfunc2], 
                                                                                conf, 
                                                                                prior_matrix, 
                                                                                '%s/tmp'%os.path.dirname(self.output), 
                                                                                nt ] ) )
                pool.close()
                pool.join()
            except Exception as e:
                print "e", e
                sys.stdout.flush()

            posterior_matrix = []
            self.logger.info("GATHER pickles")
            for nt in xrange(self.nthreads):
                tempfile = open("%s/tmp/%s_%s_th%s.deltapsi.pickle"%(os.path.dirname(self.output), self.names[0], self.names[1], nt))
                ptempt = pickle.load( tempfile )
                posterior_matrix.extend( ptempt[0] )

        pickle_path = "%s%s_%s_deltamatrix.pickle"%(output, self.names[0], self.names[1])
        pickle.dump([posterior_matrix,matched_info], open(pickle_path, 'w'))
        self.logger.info("Done!")
        return posterior_matrix, matched_info



def deltagroup(args):
    _pipeline_run(DeltaGroup(args), args.lsv)


class DeltaGroup(DeltaPair, CalcPsi):

    def _load_psis(self):
        "Load the already calculated PSI paired, and their corresponding names"
        psis = []
        names = []
        for i in xrange(len(self.files1)):
            for j in xrange(len(self.files2)):
                psi_path = self._psipaired_path(i, j)
                events_path = self._events_path(i, j)
                self.logger.info("Loading PSI %s for weights..."%psi_path)
                self.logger.info("... and corresponding event names %s"%events_path)
                psis.append(array(pickle.load(open(psi_path))))
                names.append(array(pickle.load(open(events_path))))

        return array(psis), array(names)


    def get_local_and_median_psi(self, filtered_psis, group=0):
        """Given a dictionary event names as key and the PSI distributions as values, calculate the median distribution""" 
        #get the events_paired from filtered_psis 
        events_paired = defaultdict(list) #pairs the events by event_name, disregarding experiment
        for experiment_num, event_dict in filtered_psis.items():
            for event_name, event in event_dict.items():
                events_paired[event_name].append(event) #for calculating the median psi per event
                #self.logger.debug("EVENT %s PAIRED: #events: %s Event: %s"%(event_name, len(events_paired[event_name]), event))

        #calculates the median PSI and local_values from events_paired
        median_ref = [] 
        local_values = [[] for x in xrange(self.replica_len[group])]
        for event_name, events in events_paired.items():
            if len(events) == self.replica_len[group]: #only if event is in all experiments
                median_ref.append(median(array(events), axis=0))
                median_ref[-1] /= sum(median_ref[-1])
                for i, event in enumerate(events):
                    local_values[i].append(event)
                #print "#%s EVENTS"%len(events), events
                #print "MEDIAN", median_ref[-1], sum(median_ref[-1]), '\n'

        local_values = [array(experiment) for experiment in local_values]
        #print "LOCAL VALUES", local_values
        return local_values, array(median_ref)

    def calc_weights_lsv(self, lsv_samples, info, relevant_names, nexp, ro_type=3, group=0 ):


        self.logger.info("WEIGHTS: Calculating local weights...")
        eta_wgt, nu_wgt = local_weight_eta_nu( lsv_samples, nexp )
        self.logger.info("WEIGHTS: Calculating global weights...")

        relevants = np.zeros(shape=(len(relevant_names), nexp), dtype=np.dtype('object'))
        jdx = 0
        for idx, lsv_info in enumerate(info):
            if lsv_info[1] in relevant_names:
                relevants[jdx] = lsv_samples[idx]
                jdx +=1
        ro_wgt = global_weight_ro( relevants, nexp )
        print "GLOBALS",ro_wgt
        gweights_path = "%s%s_all_weights.pickle"%(self.output, self.names[group])
        pickle.dump([ ro_wgt, eta_wgt, nu_wgt ], open(gweights_path, 'w'))
        
        gweights = np.zeros(shape=(len(lsv_samples),nexp), dtype=np.dtype('object'))
        for ne in range(nexp):
            for lidx in xrange(len(info)):
                gweights[lidx,ne] = []
                for jinc in xrange(len(eta_wgt[lidx,ne])):
                    try:
#                        gweights[lidx,ne].append( ro_wgt[ne] )
                        gweights[lidx,ne].append( ro_wgt[ne] * eta_wgt[lidx,ne][jinc] * nu_wgt[lidx,ne][jinc] )
                    except:
                        pdb.set_trace()
#            for ii in range(len(lsv_samples)):
#                gweights[ii,ne] = ro_wgt[ne]
        self.logger.info("WEIGHTS: Done")

        gweights_path = "%s%s_weights.pickle"%(self.output, self.names[group])
        pickle.dump(gweights, open(gweights_path, 'w'))
        return gweights

    def calc_weights(self, relevant, group=0):
        """
        With relevant set, calculate weigths from the PSIs between experiments (kind of delta PSI)
        Group: 0 means get weights for self.files1, 1 gets weights for self.files2 
        All this conversion between defaultdicts is probably easier using Pandas...
        """
        self.logger.info("WEIGHTS: Loading PSIs...")
        psis, event_names = self._load_psis()
        filtered_psis_dict = defaultdict(lambda : defaultdict(list)) #to hold the filtered PSI values in the relevant "best changing" set
        
        self.logger.info("WEIGHTS: Filtering PSIS with relevant set...")
        # Union of entries that paired with a changing event in the opposite group
        for k, experiment_pair in enumerate(psis):
            i, j = self.k_ref[k] #i and j are pairs for the different permutations of self.files1 and self.files2
            if group == 0: g = i
            else: g = j
            experiment = experiment_pair[0] #PSI values come in pairs of all the permutations of experiments. We only need to analyze the first one to cover everything (the ones in i)
            my_names = event_names[k] #picks up the right sequence of names 
            for event_num, event in enumerate(experiment):
                event_name = my_names[event_num]
                if event_name in relevant: #all the psis in the relevant set
                    #filtered_psis_paired[event_name].append(event)
                    filtered_psis_dict[g][event_name] = event #this way we join all duplicated PSI entries into one

        self.logger.info("WEIGHTS: Calculating intersection of events and median...")
        local_values, median_ref = self.get_local_and_median_psi(filtered_psis_dict, group)
        self.logger.info("WEIGHTS: Calculating local weights...")
        filter_lw = local_weights(local_values, self.weightsL1, median_ref)
        self.logger.info("WEIGHTS: Calculating global weights...")
        gweights = global_weights(locweights=filter_lw)
        gweights_path = "%s%s_weights.pickle"%(self.output, self.names[group])
        pickle.dump(gweights, open(gweights_path, 'w'))
        self.logger.info("WEIGHTS: Done")
        return gweights

    def equal_if_not(self, weights):
        "If weigths havent been set, automatically fix them equally"
        if not weights.size:
            numpairs = len(self.k_ref)
            return [1./numpairs]*numpairs

        return weights

    def comb_replicas_lsv(self, data_posterior_lsv, matched_info, weights1=None, weights2=None):
        "Combine all replicas per event into a single average replica. Weights per experiment can be provided"
#        weights1 = self.equal_if_not(weights1) 
#        weights2 = self.equal_if_not(weights2) 
        comb_matrix = []
        comb_names = []
        #FILTERMIN = len(self.k_ref)-1 #filter events that show on all replicas

        for lidx, lsv in enumerate(matched_info):
            comb_matrix.append([])
            njunc = len(data_posterior_lsv[lidx][0][0])
            comb = [ np.zeros(shape=(20,20),dtype=np.float) for xx in range(njunc)]
            for idx, matrices in enumerate(data_posterior_lsv[lidx]):
                for jdx, junc in enumerate(matrices):
                    for k, matrix in enumerate(junc):
                        comb[k] += matrix*weights1[lidx,idx][k]*weights2[lidx,jdx][k]
            for junc_idx, cmb in enumerate(comb):
                cmb /= sum(cmb) #renormalize so it sums 1
                comb_matrix[lidx].append(cmb)
            comb_names.append(lsv)

        return comb_matrix, comb_names

    def _pair_path(self, i, j):
        return "%s%s_%s_"%(self.output, i, j)

    def _psipaired_path(self, i, j):
        return "%s%s_%s_psipaired.pickle"%(self._pair_path(i, j), self.names[0], self.names[1])

    def _events_path(self, i, j):
        return "%s%s_%s_eventnames.pickle"%(self._pair_path(i, j), self.names[0], self.names[1]) 

    def _matrix_path(self, i, j):
        return "%s%s_%s_deltamatrix.pickle"%(self._pair_path(i, j), self.names[0], self.names[1])


    def run(self, lsv=False):

        if lsv:
            self.delta_groups_opt()
        elif self.newmodel :
            self.delta_groups_newmodel()
        else:
            self.delta_groups()

    def delta_groups( self ):
        self.logger.info("")
        self.logger.info("Running deltagroups...")
        self.logger.info("GROUP 1: %s"%self.files1)
        self.logger.info("GROUP 2: %s"%self.files2)

        num_exp = [len(self.files1), len(self.files2)]

        filtered_lsv1 = [None]*num_exp[0]
        fitfunc1 = [None]*num_exp[0]
        for ii, file in enumerate(self.files1):
            lsv_junc, const = majiq_io.load_data_lsv(file, self.logger) 

            #fitting the function
            fitfunc1[ii] = self.fitfunc(const[0])
            filtered_lsv1[ii] = self.mark_stacks_lsv( lsv_junc, fitfunc1[ii])
        filtered_lsv1 = majiq_filter.quantifiable_in_group( filtered_lsv1, self.minnonzero, self.minreads, self.logger , 0.10 )

        filtered_lsv2 = [None]*num_exp[1]
        fitfunc2 = [None]*num_exp[1]
        for ii, file in enumerate(self.files2):
            lsv_junc, const = majiq_io.load_data_lsv(file, self.logger) 

            #fitting the function
            fitfunc2[ii] = self.fitfunc(const[0])
            filtered_lsv2[ii] = self.mark_stacks_lsv( lsv_junc, fitfunc2[ii])
        filtered_lsv2 = majiq_filter.quantifiable_in_group( filtered_lsv2, self.minnonzero, self.minreads, self.logger , 0.10 )

        matched_lsv, matched_info = majiq_filter.lsv_intersection( filtered_lsv1, filtered_lsv2 )

        lsv_samples1 = [[] for xx in self.files1]
        lsv_samples1 = np.zeros(shape=(len(matched_info), len(self.files1)), dtype=np.dtype('object'))
        lsv_samples2 = np.zeros(shape=(len(matched_info), len(self.files2)), dtype=np.dtype('object'))


        fitfunc = [fitfunc1, fitfunc2]
        self.logger.info("Bootstrapping for all samples...")
        for grp_idx, group in enumerate(matched_lsv):
            for lidx, lsv_all in enumerate(group):
                for eidx, lsv in enumerate(lsv_all):
                    m_lsv, var_lsv, s_lsv = sample_from_junctions(  junction_list = lsv,
                                                                    m = self.m,
                                                                    k = self.k,
                                                                    discardzeros= self.discardzeros,
                                                                    trimborder  = self.trimborder,
                                                                    fitted_func = fitfunc[grp_idx][eidx],
                                                                    debug       = self.debug,
                                                                    Nz          = self.nz)
                    if grp_idx == 0:
                        lsv_samples1[lidx,eidx] = s_lsv
                    else:
                        lsv_samples2[lidx,eidx] = s_lsv




        self.logger.info("Calculating pairs...")

        conf = { 'minnonzero':self.minnonzero,
                 'minreads': self.minreads,
                 'm':self.m,
                 'k':self.k,
                 'discardzeros':self.discardzeros,
                 'trimborder':self.trimborder,
                 'debug':self.debug,
                 'alpha':self.alpha,
                 'n':self.n, 
                 'plotpath':self.plotpath,
                 'Nz':self.nz,
                 'names':self.names}

        if self.newmodel :
            group1, group2 = pipe.combine_for_priormatrix( matched_lsv[0], matched_lsv[1], matched_info, num_exp)
            psi_space, prior_matrix = majiq_psi.gen_prior_matrix( self, group1, group2, self.output)
            comb_matrix = pipe.model2(lsv_samples1, lsv_samples2, matched_info, num_exp, prior_matrix, psi_space, M = self.m)
            comb_names = matched_info
        else:
            relevant_events = []
            pairs_posteriors = np.zeros(shape=(len(matched_info), len(self.files1), len(self.files2)), dtype=np.dtype('object'))

            for idx in xrange(num_exp[0]):
                for jdx in xrange(num_exp[1]):
                    name1 = "%s%d"%(self.names[0],idx)
                    name2 = "%s%d"%(self.names[1],jdx)
                    self.logger.info("Calculating %s | %s"%(name1, name2))
                    matrices = delta_precomputed( matched_info, name1, name2, self.output, logger=self.logger)

                    if matrices is None:
                        lsv_exp1 = [ xx[idx] for xx in matched_lsv[0] ]
                        lsv_exp2 = [ xx[jdx] for xx in matched_lsv[1] ]
                        for xidx, xx in enumerate(matched_lsv[0]):
                            if len(xx[idx]) != len(matched_lsv[1][xidx][jdx]) : 
                                import pdb
                                pdb.set_trace()
                        psi_space, prior_matrix = majiq_psi.gen_prior_matrix(   self,
                                                                                [lsv_exp1, matched_info],
                                                                                [lsv_exp2, matched_info],
                                                                                self.output)
                        conf['psi_space'] =psi_space

                        if self.nthreads == 1:
                            Dt1_Dt2 = [lsv_samples1[:,idx],lsv_samples2[:,jdx]]
                            matrices = delta_calculation( matched_info, Dt1_Dt2, conf, prior_matrix[ii], self.logger )
                        else:
                            pool = Pool(processes=self.nthreads)
                            csize = len(matched_lsv[0]) / int(self.nthreads)
                            self.logger.info("CREATING THREADS %s with <= %s lsv"%(self.nthreads,csize))
                            jobs = []

                            for nt in xrange(self.nthreads):
                                lb = nt * csize
                                ub = min( (nt+1) * csize, len(matched_lsv[0]) )
                                if nt == self.nthreads - 1 : ub = len(matched_lsv[0]) 
                                Dt1_Dt2 = [lsv_samples1[lb:ub,idx],lsv_samples2[lb:ub,jdx]]
                                lsv_info = matched_info[lb:ub]
                                jobs.append(pool.apply_async( parallel_delta_calculation, [ lsv_info,
                                                                                            Dt1_Dt2,
                                                                                            conf,
                                                                                            prior_matrix[ii],
                                                                                            '%s/tmp'%os.path.dirname(self.output),
                                                                                            nt]))
                            pool.close()
                            pool.join()
                            posterior_matrix = []
                            self.logger.info("GATHER pickles")
                            matrices = []
                            for nt in xrange(self.nthreads):
                                tempfile = open("%s/tmp/%s_%s_th%s.deltapsi.pickle"%(os.path.dirname(self.output), self.names[0], self.names[1], nt))
                                ptempt = pickle.load( tempfile )
                                matrices.extend( ptempt[0] )

                        store_delta_psi( matrices, matched_info, name1, name2, self.output, thread_logger=self.logger)

                    if not self.fixweights1: #get relevant events for weights calculation
                        tmp = rank_deltas_lsv(matrices, matched_info, E=True)[:self.numbestchanging]
                        relevant_events.extend(tmp)

                    for lidx, matrix_lsv in enumerate(matrices):
#                        pairs_posteriors[lidx][idx].append( matrix_lsv )
                        pairs_posteriors[lidx,idx,jdx] = matrix_lsv


            self.logger.info("All pairs calculated, calculating weights...")
            self.logger.info("Obtaining weights from Delta PSI...")
            self.logger.info("... obtaining 'relevant set'")
            #sort again the combined ranks
            relevant_events.sort(key=lambda x: -x[1])
            
            rel_events = set([xx[0] for xx in relevant_events])


            self.logger.info("Obtaining weights for relevant set...")
            weights1 = self.calc_weights_lsv( lsv_samples1, matched_info, rel_events, num_exp[0], ro_type = 3, group=0 )
            self.logger.info("Weigths for %s are (respectively) %s"%(self.files1, weights1))
            weights2 = self.calc_weights_lsv( lsv_samples2, matched_info, rel_events, num_exp[1], ro_type = 3, group=1 )
            self.logger.info("Weigths for %s are (respectively) %s"%(self.files2, weights2))

            self.logger.info("Normalizing with weights...")
            comb_matrix, comb_names = self.comb_replicas_lsv( pairs_posteriors, matched_info, weights1=weights1, weights2=weights2 )
            self.logger.info("%s events matrices calculated"%len(comb_names))
        pickle_path = "%s%s_%s_deltagroup.pickle"%(self.output, self.names[0], self.names[1])
        pickle.dump([comb_matrix, comb_names], open(pickle_path, 'w'))
        self.logger.info("Alakazam! Done.")
##########


    def delta_groups_newmodel( self ):
        self.logger.info("")
        self.logger.info("Running deltagroups new model ...")
        self.logger.info("GROUP 1: %s"%self.files1)
        self.logger.info("GROUP 2: %s"%self.files2)

        exec_id = '%s_%s'%(self.names[0],self.names[1])
        tempfile = '%s/%s_temp_mid_exec.pickle'%(self.output,exec_id)
        num_exp = [len(self.files1), len(self.files2)]
        if not os.path.exists(tempfile):


            filtered_lsv1 = [None]*num_exp[0]
            fitfunc = [ [None]*num_exp[0], [None]*num_exp[1] ]
            for ii, file in enumerate(self.files1):
                lsv_junc, const = majiq_io.load_data_lsv(file, self.logger) 

                #fitting the function
                fitfunc[0][ii] = self.fitfunc(const[0])
                filtered_lsv1[ii] = self.mark_stacks_lsv( lsv_junc, fitfunc[0][ii])
            filtered_lsv1 = majiq_filter.quantifiable_in_group( filtered_lsv1, self.minnonzero, self.minreads, self.logger , 0.10 )

            filtered_lsv2 = [None]*num_exp[1]
            for ii, file in enumerate(self.files2):
                lsv_junc, const = majiq_io.load_data_lsv(file, self.logger) 

                #fitting the function
                fitfunc[1][ii] = self.fitfunc(const[0])
                filtered_lsv2[ii] = self.mark_stacks_lsv( lsv_junc, fitfunc[1][ii])
            filtered_lsv2 = majiq_filter.quantifiable_in_group( filtered_lsv2, self.minnonzero, self.minreads, self.logger , 0.10 )

            matched_lsv, matched_info = majiq_filter.lsv_intersection( filtered_lsv1, filtered_lsv2 )


            group1, group2 = pipe.combine_for_priormatrix( matched_lsv[0], matched_lsv[1], matched_info, num_exp)
            psi_space, prior_matrix = majiq_psi.gen_prior_matrix( self, group1, group2, self.output)


            #TEMP 
            tout = open(tempfile, 'w+')
            pickle.dump([matched_info, matched_lsv, psi_space, prior_matrix, fitfunc], tout)
            tout.close()
            #END TEMP

        else:
            matched_info, matched_lsv, psi_space, prior_matrix, fitfunc = pickle.load(open(tempfile))

        conf = { 'minnonzero':self.minnonzero,
                 'minreads': self.minreads,
                 'm':self.m,
                 'k':self.k,
                 'discardzeros':self.discardzeros,
                 'trimborder':self.trimborder,
                 'debug':self.debug,
                 'alpha':self.alpha,
                 'n':self.n, 
                 'plotpath':self.plotpath,
                 'nz':self.nz,
                 'names':self.names}

        get_clean_raw_reads( matched_info, matched_lsv[0], self.output, self.names[0], num_exp[0] )
        get_clean_raw_reads( matched_info, matched_lsv[1], self.output, self.names[1], num_exp[1] )

        if self.nthreads == 1:
            posterior_matrix, names = pipe.model2( matched_lsv, matched_info, num_exp, conf, prior_matrix, fitfunc, psi_space, self.logger)
        else:


            pool = Pool(processes=self.nthreads)
            csize = len(matched_lsv[0]) / int(self.nthreads)
            self.logger.info("CREATING THREADS %s with <= %s lsv"%(self.nthreads,csize))
            jobs = []

            for nt in xrange(self.nthreads):
                lb = nt * csize
                ub = min( (nt+1) * csize, len(matched_lsv[0]) )
                if nt == self.nthreads - 1 : ub = len(matched_lsv[0]) 
                lsv_list = [matched_lsv[0][lb:ub],matched_lsv[1][lb:ub]]
                lsv_info = matched_info[lb:ub]
                pool.apply_async( pipe.parallel_lsv_child_calculation, [ pipe.model2, 
                                                                         [ lsv_list, lsv_info, num_exp, conf, prior_matrix, fitfunc, psi_space ],
#                pipe.parallel_lsv_child_calculation ( [ matched_lsv, matched_info, num_exp, conf, prior_matrix, fitfunc, psi_space ],
                                                                         matched_info,
                                                                         '%s/tmp'%os.path.dirname(self.output),
                                                                         '%s_%s'%(self.names[0], self.names[1]),
#                                                                         nt )
                                                                         nt] )
            pool.close()
            pool.join()

            posterior_matrix = []
            names = []
            self.logger.info("GATHER pickles")
            for nt in xrange(self.nthreads):
                tempfile = open("%s/tmp/%s_%s_th%s.deltapsi.pickle"%(os.path.dirname(self.output), self.names[0], self.names[1], nt))
                ptempt = pickle.load( tempfile )
                posterior_matrix.extend( ptempt[0] )
                names.extend(ptempt[1])


        pickle_path = "%s%s_%s_deltagroup.pickle"%(self.output, self.names[0], self.names[1])
        pickle.dump([posterior_matrix, names], open(pickle_path, 'w'))
        self.logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s"%(self.names[0],self.names[1], self.output))
        self.logger.info("Alakazam! Done.")

##########



    def empirical_relevant_set( self, matched_lsv, matched_info, num_exp, logger=None ):

        relevant_events = []
        for iidx in xrange(num_exp[0]):
            for jidx in xrange(num_exp[1]):
                lsv_exp1 = [ xx[iidx] for xx in matched_lsv[0] ]
                lsv_exp2 = [ xx[jidx] for xx in matched_lsv[1] ]

                filtered_lsv1 = majiq_filter.lsv_quantifiable([lsv_exp1, matched_info], minnonzero=10, min_reads=20, logger=logger)
                filtered_lsv2 = majiq_filter.lsv_quantifiable([lsv_exp2, matched_info], minnonzero=10, min_reads=20, logger=logger)

                ids1 = set([xx[1] for xx in filtered_lsv1[1]])
                ids2 = set([xx[1] for xx in filtered_lsv2[1]])
                matched_names = ids1.intersection(ids2)
                best_set_mean1 = [[],[]]
                best_set_mean2 = [[],[]]

                for ii in matched_names:
                    for idx, nm in enumerate(filtered_lsv1[1]):
                        if nm[1] == ii:
                            nz = np.count_nonzero(filtered_lsv1[0][idx])
                            best_set_mean1[0].append(nz * majiq_sample.mean_junction(filtered_lsv1[0][idx]))
                            best_set_mean1[1].append(filtered_lsv1[1][idx])
                            break
                    for idx, nm in enumerate(filtered_lsv2[1]):
                        if nm[1] == ii:
                            nz = np.count_nonzero(filtered_lsv2[0][idx])
                            best_set_mean2[0].append(nz * majiq_sample.mean_junction(filtered_lsv2[0][idx]))
                            best_set_mean2[1].append(filtered_lsv2[1][idx])
                            break

                best_delta_psi = majiq_psi.empirical_delta_psi(best_set_mean1[0], best_set_mean2[0])

#                import pdb
#                pdb.set_trace()
                relevant_events.extend(rank_empirical_delta_lsv(best_delta_psi, matched_info)[:self.numbestchanging])

        relevant_events.sort(key=lambda x: -x[1])
        
        rel_events = set([xx[0] for xx in relevant_events])
        return rel_events


