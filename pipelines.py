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
from analysis.weight import local_weights, global_weights
from analysis.matrix import rank_deltas, collapse_matrix
import analysis.filter as majiq_filter
import analysis.io as majiq_io
import analysis.psi as  majiq_psi
#from analysis.filter import norm_junctions, discardlow, discardhigh, discardminreads, discardmaxreads, discardminreads_and, mark_stacks

################################
# Data loading and Boilerplate #
################################



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

        if not os.path.isdir(tempdir):
            os.mkdir(tempdir)
        thread_logger = get_logger("%s/majiq.w%s.log"%(tempdir,chunk), silent=False, debug=conf['debug'])
        thread_logger.info( "[Th %s]: START child,%s"%(chunk,current_process().name))
        thread_logger.info('[Th %s]: Filtering ...'%(chunk))


        num_lsv = len(lsv_junc[0])
        ''' 
            fon[0] = False deactivates num_reads >= 20
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

class CalcPsi(BasicPipeline):

    def run(self, lsv=False):
        ret = []
        for path in self.files:
            if lsv :
                ret.append(self.calcpsi_lsv(path))
            else:
                self.calcpsi(path)
        return ret


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




def deltapsi_calc( matched_list, matched_info, fitfunc, conf, chunk, prior_matrix, logr):

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
                                    "Posterior Delta Event %s.%s (Psi1: %s Psi2: %s)"%(lsv_idx,psi, sum(data_given_psi1[psi]), sum(data_given_psi2[psi])), 
                                    "posterior_dpsi.%s.%s"%(lsv[1],psi),
                                    conf['plotpath'] )
        posterior_matrix.append(lsv_psi_matrix)

    return posterior_matrix


def parallel_delta_psi_wrapper ( matched_list, matched_info, fitfunc, conf, prior_matrix, tempdir, chunk):

    if not os.path.isdir(tempdir):
            os.mkdir(tempdir)
    thread_logger = get_logger("%s/majiq.w%s.log"%(tempdir,chunk), silent=False, debug=conf['debug'])
    thread_logger.info( "[Th %s]: START child,%s"%(chunk,current_process().name))
    thread_logger.info('[Th %s]: Filtering ...'%(chunk))

    post_matrix = deltapsi_calc(matched_list, matched_info, fitfunc, conf, chunk, prior_matrix, thread_logger)
    

    print "%s/%s_%s_th%s.deltapsi.pickle"%(tempdir, conf['names'][0], conf['names'][1], chunk)
    sys.stdout.flush()
    thread_logger.info("[Th %s]: Saving PSI..."%chunk)
    output = open("%s/%s_%s_th%s.deltapsi.pickle"%(tempdir, conf['names'][0], conf['names'][1], chunk), 'w')
    pickle.dump((post_matrix, matched_info), output)
    thread_logger.info("[Th %s]: PSI calculation for %s ended succesfully! Result can be found at %s"%(chunk, name, output.name))
    return

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
    _pipeline_run(DeltaGroup(args))


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

    def calc_weights_lsv(self, group, releveant, type=3 ):


        eta_wgt = local_weight_eta( group )


        nu_wgt = local_weight_nu( group, lsv_medians )

        ro_wgt = global_weight_ro( group, lsv_medians )

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

    def comb_replicas_lsv(self, delta_posterior_lsv, weights1=None, weights2=None):
        "Combine all replicas per event into a single average replica. Weights per experiment can be provided"
        weights1 = self.equal_if_not(weights1) 
        weights2 = self.equal_if_not(weights2) 
        comb_matrix = []
        comb_names = []
        FILTERMIN = len(self.k_ref)-1 #filter events that show on all replicas

        for l_idx,lsv in enumerate( delta_posterior_lsv):
            comb_matrix.append([])
#            name = 
            for  matrices in lsv.items():
                for k, matrix in enumerate(matrices):
                    #i = k / len(self.files1) #infers the weight1 index given the position of the matrix
                    #j = k % len(self.files2)
                    i, j = self.k_ref[k]
                    comb += matrix*weights1[l_idx][i]*weights2[l_idx][j]

                if k == FILTERMIN:
                    comb /= sum(comb) #renormalize so it sums 1
                    comb_matrix[l_idx].append(comb)
                    comb_names.append(name)

        return comb_matrix, comb_names


    def comb_replicas(self, pairs_posteriors, weights1=None, weights2=None):
        "Combine all replicas per event into a single average replica. Weights per experiment can be provided"
        weights1 = self.equal_if_not(weights1) 
        weights2 = self.equal_if_not(weights2) 
        comb_matrix = []
        comb_names = []
        FILTERMIN = len(self.k_ref)-1 #filter events that show on all replicas
        for name, matrices in pairs_posteriors.items():
            for k, matrix in enumerate(matrices):
                #i = k / len(self.files1) #infers the weight1 index given the position of the matrix
                #j = k % len(self.files2)
                i, j = self.k_ref[k]
                if k == 0: #first iteration
                    comb = matrix*weights1[i]*weights2[j]
                else:
                    comb += matrix*weights1[i]*weights2[j]

            if k == FILTERMIN:
                comb /= sum(comb) #renormalize so it sums 1
                comb_matrix.append(comb)
                comb_names.append(name)

        return comb_matrix, comb_names

    def _pair_path(self, i, j):
        return "%s%s_%s_"%(self.output, i, j)

    def _psipaired_path(self, i, j):
        return "%s%s_%s_psipaired.pickle"%(self._pair_path(i, j), self.names[0], self.names[1])

    def _events_path(self, i, j):
        return "%s%s_%s_eventnames.pickle"%(self._pair_path(i, j), self.names[0], self.names[1]) 

    def _matrix_path(self, i, j):
        return "%s%s_%s_deltamatrix.pickle"%(self._pair_path(i, j), self.names[0], self.names[1])

    def run_lsv(self):
        self.logger.info("")
        self.logger.info("Running deltagroups...")
        self.logger.info("GROUP 1: %s"%self.files1)
        self.logger.info("GROUP 2: %s"%self.files2)
        self.logger.info("Calculating pairs...")
        pairs_posteriors = defaultdict(list)
        self.k_ref = [] #maps the i, j reference for every matrix (this could be do in less lines with div and mod, but this is safer) 
        relevant_events = []

        for i in xrange(len(self.files1)):
            for j in xrange(len(self.files2)):
                self.k_ref.append([i, j])
                path = self._pair_path(i, j)
                matrix_path = self._matrix_path(i, j)
                if os.path.exists(matrix_path):
                    self.logger.info("%s exists! Loading..."%path)
                    matrices[i,j], info[i,j] = pickle.load(open(matrix_path))
                else:
                    self.logger.info("Calculating pair %s"%path)
                    matrices[i,j], info[i,j] = self.pairdelta_lsv(self.files1[i], self.files2[j], path)
                    self.logger.info("Saving pair posterior for %s, %s"%(i, j))

                if not self.fixweights1: #get relevant events for weights calculation
                    relevant_events.extend(rank_deltas_lsv(matrices, info, E=True)[:self.numbestchanging])
                    pairs_posteriors[name].append(matrices[k]) #pairing all runs events

        self.logger.info("All pairs calculated, calculating weights...")
        self.logger.info("Obtaining weights from Delta PSI...")
        self.logger.info("... obtaining 'relevant set'")
        #sort again the combined ranks
        relevant_events.sort(key=lambda x: -x[1])
        
        self.logger.info("Obtaining weights for relevant set...")
        weights1 = self.calc_weights_lsv(relevant_events, ro_type = self.ro_type )
        self.logger.info("Weigths for %s are (respectively) %s"%(self.files1, weights1))
        weights2 = self.calc_weights_lsv(relevant_events, ro_type = self.ro_type )
        self.logger.info("Weigths for %s are (respectively) %s"%(self.files2, weights2))

        self.logger.info("Normalizing with weights...")
        comb_matrix, comb_names = self.comb_replicas_lsv( pairs_posteriors, weights1=weights1, weights2=weights2 )
        self.logger.info("%s events matrices calculated"%len(comb_names))
        pickle_path = "%s%s_%s_deltacombmatrix.pickle"%(self.output, self.names[0], self.names[1])
        name_path = "%s%s_%s_combeventnames.pickle"%(self.output, self.names[0], self.names[1])
        pickle.dump(comb_matrix, open(pickle_path, 'w'))
        pickle.dump(comb_names, open(name_path, 'w'))
        self.logger.info("Alakazam! Done.")

    def run(self):
        self.logger.info("")
        self.logger.info("Running deltagroups...")
        self.logger.info("GROUP 1: %s"%self.files1)
        self.logger.info("GROUP 2: %s"%self.files2)
        self.logger.info("Calculating pairs...")
        pairs_posteriors = defaultdict(list)
        self.k_ref = [] #maps the i, j reference for every matrix (this could be do in less lines with div and mod, but this is safer) 
        relevant_events = []
        for i in xrange(len(self.files1)):
            for j in xrange(len(self.files2)):
                self.k_ref.append([i, j])
                path = self._pair_path(i, j)
                matrix_path = self._matrix_path(i, j)
                events_path = self._events_path(i, j)
                if os.path.exists(matrix_path):
                    self.logger.info("%s exists! Loading..."%path)
                    matrices = pickle.load(open(matrix_path))
                    names = pickle.load(open(events_path))
                else:
                    self.logger.info("Calculating pair %s"%path)
                    matrices, names = self.pairdelta(self.files1[i], self.files2[j], path)
                    self.logger.info("Saving pair posterior for %s, %s"%(i, j))

                if not self.fixweights1: #get relevant events for weights calculation
                    relevant_events.extend(rank_deltas_lsv(matrices, names, E=True)[:self.numbestchanging])

                for k, name in enumerate(names):
                    pairs_posteriors[name].append(matrices[k]) #pairing all runs events

        self.logger.info("All pairs calculated, calculating weights...")
        if self.fixweights1: #Weights are manually fixed
            weights1 = self.fixweights1
            weights2 = self.fixweights2
        else:
            self.logger.info("Obtaining weights from Delta PSI...")
            self.logger.info("... obtaining 'relevant set'")
            #sort again the combined ranks
            relevant_events.sort(key=lambda x: -x[1])

            #gather the first NUMEVENTS names
            relevant = []
            for name, matrix in relevant_events:
                if name not in relevant_events: #ignore duplicated entries
                    relevant.append(name)

                if len(relevant) == self.numbestchanging:
                    break #we got enough elements
                    
            self.logger.info("Obtaining weights for relevant set...")
            weights1 = self.calc_weights(relevant, group=0)
            self.logger.info("Weigths for %s are (respectively) %s"%(self.files1, weights1))
            weights2 = self.calc_weights(relevant, group=1)
            self.logger.info("Weigths for %s are (respectively) %s"%(self.files2, weights2))

        self.logger.info("Normalizing with weights...")
        comb_matrix, comb_names = self.comb_replicas(pairs_posteriors, weights1=weights1, weights2=weights2)        
        self.logger.info("%s events matrices calculated"%len(comb_names))
        pickle_path = "%s%s_%s_deltacombmatrix.pickle"%(self.output, self.names[0], self.names[1])
        name_path = "%s%s_%s_combeventnames.pickle"%(self.output, self.names[0], self.names[1])
        pickle.dump(comb_matrix, open(pickle_path, 'w'))
        pickle.dump(comb_names, open(name_path, 'w'))
        self.logger.info("Alakazam! Done.")



