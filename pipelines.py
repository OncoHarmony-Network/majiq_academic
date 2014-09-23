import os
from collections import defaultdict
import abc
import pickle
from multiprocessing import Pool, current_process

from pylab import *
from numpy.ma import masked_less
from grimoire.utils.utils import create_if_not_exists, get_logger
from analysis.polyfitnb import fit_nb
from analysis.weight import local_weight_eta_nu, global_weight_ro
from analysis.matrix import rank_deltas_lsv, rank_empirical_delta_lsv
import analysis.filter as majiq_filter
import analysis.io as majiq_io
import analysis.psi as  majiq_psi
import analysis.sample as majiq_sample

import pipe as pipe
################################
# Data loading and Boilerplate #
################################


def get_clean_raw_reads(matched_info, matched_lsv, outdir, names, num_exp):

    res = []
    for eidx in xrange(num_exp):
        for ldx, lsv in enumerate(matched_info):
            num = matched_lsv[ldx][eidx].sum()
            res.append([lsv[1],num])

        with open('%s/clean_reads.%s%d.pkl' % (outdir, names, eidx), 'wb') as fp:
            pickle.dump(res, fp)


def _pipeline_run(pipeline, lsv=False, logger=None):
    """ Exception catching for all the pipelines """
    try:
        return pipeline.run(lsv)
    except KeyboardInterrupt:
        if pipeline.logger: pipeline.logger.info("MAJIQ manually interrupted. Avada kedavra...")

# PLOTTING STUFF. All this should eventually go to a plotting module 
def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        savefig("%s%s.png" % (plotpath, plotname), bbox_inches='tight')
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
        logger_path = self.logger
        if not logger_path:
            logger_path = self.output

        self.logger = get_logger("%smajiq.log" % logger_path, silent=self.silent, debug=self.debug)
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

    def gc_content_norm_lsv(self, lsv_list, const_list):
        "Normalize the matrix using the gc content"
        self.logger.info("GC content normalization...")
        if self.gcnorm:
            for lidx, lsv in enumerate(lsv_list[0]):
                lsv = lsv * lsv_list[2][lidx]
            const_list[0] = const_list[0] * const_list[2]
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

def __parallel_calcpsi_lsv(conf, lsv_junc, fitfunc, name, chunk, tempdir):

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

            m_lsv, var_lsv, s_lsv = majiq_sample.sample_from_junctions(  ii,
                                                            conf['m'],
                                                            conf['k'],
                                                            discardzeros= conf['discardzeros'],
                                                            trimborder  = conf['trimborder'],
                                                            fitted_func = fitfunc,
                                                            debug       = conf['debug'],
                                                            Nz          = conf['Nz'])
            lsv_sample.append( s_lsv )

        thread_logger.info("[Th %s]: Calculating PSI for %s ..."%(chunk, name))
        psi = majiq_psi.lsv_psi(lsv_sample, conf['alpha'], conf['n'], conf['debug'])

        thread_logger.info("[Th %s]: Saving PSI..."%chunk)
        output = open("%s/%s_th%s.psi.pickle"%(tempdir, name, chunk), 'w')
        pickle.dump((psi, lsv_junc[1]), output)
        thread_logger.info("[Th %s]: PSI calculation for %s ended succesfully! Result can be found at %s"%(chunk, name, output.name))
    except Exception as e:
        print "%s"%sys.exc_traceback.tb_lineno, e
        sys.stdout.flush()


class CalcPsi(BasicPipeline):

    def run(self, lsv=False):
        self.calcpsi()

    def calcpsi(self):
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

        conf = {'minnonzero': self.minnonzero,
                'minreads': self.minreads,
                'm': self.m,
                'k': self.k,
                'discardzeros': self.discardzeros,
                'trimborder': self.trimborder,
                'debug': self.debug,
                'alpha': self.alpha,
                'n': self.n,
                'nbins': 40,
                'nz': self.nz}

        if self.nthreads == 1:
            posterior_matrix, names = pipe.calcpsi(matched_lsv, matched_info, num_exp, conf, fitfunc, self.logger)
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



################################
#          Delta PSI           #
################################

def deltapair(args):
    _pipeline_run(DeltaPair(args), args.lsv)


class DeltaPair(BasicPipeline):

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



    def equal_if_not(self, weights):
        "If weigths havent been set, automatically fix them equally"
        if not weights.size:
            numpairs = len(self.k_ref)
            return [1./numpairs]*numpairs

        return weights

    def _pair_path(self, i, j):
        return "%s%s_%s_"%(self.output, i, j)

    def _psipaired_path(self, i, j):
        return "%s%s_%s_psipaired.pickle"%(self._pair_path(i, j), self.names[0], self.names[1])

    def _events_path(self, i, j):
        return "%s%s_%s_eventnames.pickle"%(self._pair_path(i, j), self.names[0], self.names[1]) 

    def _matrix_path(self, i, j):
        return "%s%s_%s_deltamatrix.pickle"%(self._pair_path(i, j), self.names[0], self.names[1])


    def run(self):
        self.delta_groups()

    def delta_groups(self):
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
                                                                         matched_info,
                                                                         '%s/tmp'%os.path.dirname(self.output),
                                                                         '%s_%s'%(self.names[0], self.names[1]),
                                                                         nt] )
            pool.close()
            pool.join()

            posterior_matrix = []
            names = []
            self.logger.info("GATHER pickles")
            for nt in xrange(self.nthreads):
                tempfile = open("%s/tmp/%s_%s_th%s.%s.pickle"%(os.path.dirname(self.output), self.names[0], self.names[1], nt, pipe.model2.__name__))
                ptempt = pickle.load( tempfile )
                posterior_matrix.extend( ptempt[0] )
                names.extend(ptempt[1])


        pickle_path = "%s%s_%s.%s.pickle"%(self.output, self.names[0], self.names[1],pipe.model2.__name__)
        pickle.dump([posterior_matrix, names], open(pickle_path, 'w'))
        self.logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s"%(self.names[0],self.names[1], self.output))
        self.logger.info("Alakazam! Done.")
