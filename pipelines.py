import os
from collections import defaultdict
import abc
import pickle

from pylab import *
from numpy.ma import masked_less
from scipy.io import loadmat
from scipy.stats import norm
from scipy.sparse import lil_matrix
from scipy.stats import scoreatpercentile

from grimoire.utils.utils import create_if_not_exists, get_logger
from analysis.polyfitnb import fit_nb 
from analysis.filter import norm_junctions, discardlow, discardhigh, discardminreads, discardmaxreads, discardminreads_and, mark_stacks
from analysis.sample import sample_from_junctions, mean_junction
from analysis.psi import calc_psi, mean_psi, simple_psi, DirichletCalc, reads_given_psi, BINS_CENTER
from analysis.adjustdelta import adjustdelta
from analysis.weight import local_weights, global_weights
from analysis.matrix import rank_deltas, collapse_matrix


################################
# Data loading and Boilerplate #
################################

DELTA_RATIO = 0.2 #TODO to parameters


def _find_len(grimoire_obj, logger):
    junc_len = 666
    for junction in grimoire_obj:
        if hasattr(junction[0], 'coverage'):
            if type(junction[0].coverage[0]) == lil_matrix:
                event = junction[0].coverage[0].toarray()
                junc_len = event[0].shape[0]
            else:
                junc_len = junction[0].coverage[0].shape[0]

            if junc_len != 1: 
                if logger: logger.debug("Junction length is %s"%junc_len)
                return junc_len  

def _find_len2(grimoire_obj, logger):
    junc_len = 666
    for junction in grimoire_obj:
        if hasattr(junction, 'coverage'):
            if type(junction.coverage[0]) == lil_matrix:
                event = junction.coverage[0].toarray()
                junc_len = event[0].shape[0]
            else:
                junc_len = junction.coverage[0].shape[0]

            if junc_len != 1: 
                if logger: logger.debug("Junction length is %s"%junc_len)
                return junc_len  



def _load_data_const(grimoire_obj, logger=None):
    CONST_MAX = 2000 # we don't really need more than 2000 constitutive exons 
    """
    Overriding Jordis objects. 
    Should be deleted at some point as the majiq.analysis should instead read them and the object should be extended

    Plus, this should change when the gc_content becomes available again.
    """
    ret = []
    #first get junction length, should be almost always the first entry
    my_len = _find_len(grimoire_obj, logger)

    for i, junction in enumerate(grimoire_obj):
        if i == CONST_MAX: # we don't need more constitutive junctions
            return array(ret)
            
        if type(junction[0].coverage[0]) == lil_matrix:
            event = junction[0].coverage[0].toarray()
            event = event[0]
        else:
            event = junction[0].coverage[0]

        if event.shape[0] == 1:
            event = [0]*my_len


        ret.append(list(event))

    return array(ret)

def _load_data2(grimoire_obj, logger=None, getnames = False):
    """
    Same as above
    """
    ret = []
    #first get junction length, should be almost always the first entry
    my_len = _find_len2(grimoire_obj, logger)
    for i, junction in enumerate(grimoire_obj):  
        if type(junction.coverage[0]) == lil_matrix:
            event = junction.coverage[0].toarray()
            event = event[0]
        else:
            event = junction.coverage[0]

        if event.shape[0] == 1:
            event = [0]*my_len


        ret.append(list(event))

    return array(ret)

def load_data(path, logger=None):
    "Load data from the preprocess step. Could change to a DDBB someday"
    data = pickle.load(open(path))
    return _load_data2(data[1][:,0], logger), _load_data2(data[1][:,1], logger), _load_data_const(data[2], logger) #inc, exc, const

def load_data_n(paths, logger=None, const=False):
    "Load and pair n replicas (for the groups). Doesn't get the constitutive data"
    datas = []
    for path in paths:
        if logger: logger.info("Loading %s..."%path)
        datas.append(pickle.load(open(path))) 

    if logger: logger.info("Combining all matrices into 'paired_samples'...")
    paired_samples = defaultdict(list)
    for data in datas:
        _load_data_pair2(data[1], paired_samples, logger)

    pair_length = len(paths) * 2 # the length of an event that has all replicas is double the number of experiments (each experiment has both inc and exc) 

    #print "PAIRED SAMPLES", len(paired_samples)
    event_names = [] 
    ret = []
    for event_name, junctions in paired_samples.items():
        if len(junctions) == pair_length:
            ret.append(junctions)
            event_names.append(event_name)

    return ret, event_names


def _load_data_pair2(data, paired_samples, logger=None):
    "Load the information in data inside paired_samples defaultdict"
    my_len = _find_len(data, logger)
    for i, junction in enumerate(data): #iterate junctions in 1 and add to the inc and exc to the dictionary
        if hasattr(junction[0], 'coverage') and hasattr(junction[1], 'coverage'):

            if type(junction[0].coverage[0]) == lil_matrix:  
                my_inc = junction[0].coverage[0].toarray()
                my_inc = my_inc[0]
            else:
                my_inc = junction[0].coverage[0]

            if type(junction[1].coverage[0]) == lil_matrix:  
                my_exc = junction[1].coverage[0].toarray()
                my_exc = my_exc[0]
            else:
                my_exc = junction[1].coverage[0]

            if my_inc.shape[0] == 1:
                my_inc = [0]*my_len

            if my_exc.shape[0] == 1:
                my_exc = [0]*my_len

            my_name = junction[0].name
            if not my_name:
                my_name = junction[1].name

            paired_samples[my_name].extend([list(my_inc), list(my_exc)])    

def load_data_pair(path1, path2, logger=None):
    """Pairing functionality should be extracted of this function"""
    #const extracting doesnt change
    data1 = pickle.load(open(path1))
    data2 = pickle.load(open(path2))
    const1 = _load_data_const(data1[2], logger)
    const2 = _load_data_const(data2[2], logger)

    paired_samples = defaultdict(list)
    _load_data_pair2(data1[1], paired_samples, logger)
    _load_data_pair2(data2[1], paired_samples, logger)

    #pair the events: Only keep the paired events
    ret = []
    event_names = [] 
    for event_name, junctions in paired_samples.items():
        if len(junctions) == 4:
            ret.append(junctions)
            event_names.append(event_name)

    ret = array(ret)
    return ret[:,0], ret[:,1], const1, ret[:,2], ret[:,3], const2, event_names

def _pipeline_run(pipeline, logger=None):
    "Exception catching for all the pipelines"
    try:
        pipeline.run()

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

def plot_matrix(matrix, my_title, plotname, plotpath):
    clf()
    ax = subplot(1,1,1)
    title(my_title)
    imshow(matrix)
    xlabel(u"PSI i")
    ylabel(u"PSI j")
    ax.set_xticklabels([0, 0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels([0, 0, 0.25, 0.5, 0.75, 1])

    _save_or_show(plotpath, plotname=plotname)

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

    @abc.abstractmethod
    def run(self):
        """This is the entry point for all pipelines"""
        return

    def gc_content_norm(self, all_junctions):
        "Normalize the matrix using the gc content"
        self.logger.info("GC content normalization...")
        if self.gcnorm:
            for junc_set in all_junctions.keys():
                all_junctions[junc_set] = norm_junctions(all_junctions[junc_set]["junctions"], all_junctions[junc_set]["gc_content"])

        return all_junctions

    def fitfunc(self, const_junctions):
        "Calculate the Negative Binomial function to sample from using the Constitutive events"
        if self.debug:
            self.logger.debug("Skipping fitfunc because --debug!")
            return poly1d([1, 0])
        else:
            self.logger.info("Fitting NB function with constitutive events...")
            return fit_nb(const_junctions, "%s_nbfit"%self.output, self.plotpath, nbdisp=self.nbdisp, logger=self.logger, discardb=True)

    def mark_stacks(self, all_junctions, fitfunc):
        if self.markstacks >= 0:
            self.logger.info("Marking and masking stacks...")
            for junc_set in all_junctions.keys():
                if junc_set.find("const") == -1:
                    print "... %s"%junc_set
                    all_junctions[junc_set] = mark_stacks(all_junctions[junc_set], fitfunc, self.markstacks, self.nbdisp)

                all_junctions[junc_set] = masked_less(all_junctions[junc_set], 0) #remask the stacks

        return all_junctions

################################
# PSI calculation pipeline     #
################################

def calcpsi(args):
    _pipeline_run(CalcPsi(args))

class CalcPsi(BasicPipeline):

    def run(self):
        for path in self.files:
            self.calcpsi(path)


    def calcpsi(self, path, write_pickle=True):
        """
        Given a file path with the junctions, return psi distributions. 
        write_pickle indicates if a .pickle should be saved in disk
        """
        name = os.path.basename(path)
        self.logger.info("")
        self.logger.info("Loading %s..."%path)
        inc, exc, const = load_data(path, self.logger) 
        self.logger.debug("SHAPES for inclusion, %s exclusion, %s constitutive %s"%(inc.shape, exc.shape, const.shape))
        self.logger.info("Loaded.")
        all_junctions = {"inc": inc, "exc": exc, "const": const }
        all_junctions = self.gc_content_norm(all_junctions)

        self.logger.info("Masking non unique...")
        for junc_set in all_junctions.keys():
            all_junctions[junc_set] = masked_less(all_junctions[junc_set], 0) 

        fitfunc = self.fitfunc(all_junctions["const"])
        all_junctions = self.mark_stacks(all_junctions, fitfunc)
        #FILTER_JUNCTIONS?
        self.logger.info('Filtering ...')
        filter_junctions = defaultdict(array)
        filter_junctions["exc"], filter_junctions["inc"] = discardlow(self.minnonzero, True, self.logger, all_junctions["exc"], all_junctions["inc"])
        filter_junctions["exc"], filter_junctions["inc"] = discardminreads(self.minreads, True, self.logger, False, filter_junctions["exc"], filter_junctions["inc"])
        
        self.logger.info("Bootstrapping samples...") 
        mean_exc, var_exc, exc_samples = sample_from_junctions(filter_junctions["exc"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc, debug=self.debug)
        mean_inc, var_inc, inc_samples = sample_from_junctions(filter_junctions["inc"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc, debug=self.debug)      

        self.logger.info("\nCalculating PSI for %s ..."%(name))
        psi = calc_psi(inc_samples, exc_samples, name, self.alpha, self.n, self.debug, self.psiparam)

        self.logger.info("Saving PSI...")
        if write_pickle:
            output = open("%s%s_psi.pickle"%(self.output, name), 'w')
            pickle.dump(psi, output)
            self.logger.info("PSI calculation for %s ended succesfully! Result can be found at %s"%(name, output.name))

        return psi



################################
#          Delta PSI           #
################################

def deltapair(args):
    _pipeline_run(DeltaPair(args))

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

    def run(self):
        self.pairdelta(self.file1, self.file2, self.output)

    def pairdelta(self, file1, file2, output):
        self.logger.info("")
        self.logger.info("Processing pair %s - %s..."%(file1, file2))
        inc1, exc1, const1, inc2, exc2, const2, event_names = load_data_pair(file1, file2, self.logger) 
        self.logger.debug("Shapes for inclusion, %s exclusion, %s constitutive %s"%(inc1.shape, exc1.shape, const1.shape))
        all_junctions = {"inc1": inc1, "exc1": exc1, "const1": const1, "inc2": inc2, "exc2": exc2, "const2": const2 }
        all_junctions = self.gc_content_norm(all_junctions)
        self.logger.info("Masking non unique...")
        for junc_set in all_junctions.keys():
            #print junc_set, all_junctions[junc_set], all_junctions[junc_set].shape
            all_junctions[junc_set] = masked_less(array(all_junctions[junc_set]), 0) 

        fitfunc1 = self.fitfunc(all_junctions["const1"])
        fitfunc2 = self.fitfunc(all_junctions["const2"])

        if self.markstacks >= 0:
            self.logger.info("Marking and masking stacks...")
            for junc_set in all_junctions.keys():
                if junc_set.find("const") == -1:
                    if junc_set.endswith("1"):
                        f = fitfunc1
                    else:
                        f = fitfunc2

                    print "... %s"%junc_set
                    all_junctions[junc_set] = mark_stacks(all_junctions[junc_set], f, self.markstacks, self.nbdisp)

        #Start prior matrix
        self.logger.info("Calculating prior matrix...")
        numbins = 20 #half the delta bins TODO to parameters

        self.logger.info("Calculate jefferies matrix...")
        dircalc = DirichletCalc() 
        #Adjust prior matrix with Jefferies prior        
        jefferies = []
        psi_space = linspace(0, 1-self.binsize, num=numbins) + self.binsize/2
        for i in psi_space:
            jefferies.append([])
            for j in psi_space:
                jefferies[-1].append(dircalc.pdf([i, 1-i, j, 1-j], [self.alpha, self.alpha, self.alpha, self.alpha]))

        #jefferies = array([dircalc.pdf([x, 1-x], [0.5, 0.5]) for x in psi_space])
        jefferies = array(jefferies)
        jefferies /= sum(jefferies)
        plot_matrix(jefferies, "Jefferies Matrix", "jefferies_matrix", self.plotpath)
 
        if self.synthprior:
            #Use a synthetic matrix to generate the values
            prior_matrix = [] 
            uniform = self.prioruniform/numbins 
            mydist = norm(loc=0, scale=self.priorstd)
            norm_space = linspace(-1, 1-self.binsize, num=numbins*2) + self.binsize/2
            pdfnorm = mydist.pdf(norm_space)

            newdist = (pdfnorm+uniform)/(pdfnorm+uniform).sum()

            plot(linspace(-1, 1, num=len(list(pdfnorm))), pdfnorm)
            _save_or_show(self.plotpath, plotname="prior_distribution")
            #generate the matrix
            for i in xrange(numbins):
                prior_matrix.append(list(newdist[numbins-i:(numbins*2)-i]))

            prior_matrix = array(prior_matrix)
            prior_matrix /= sum(prior_matrix) #renormalize so it sums 1
            self._get_delta_info(newdist, norm_space)
            plot_matrix(prior_matrix, "Prior Matrix (before Jefferies)", "prior_matrix_no_jefferies", self.plotpath)

        elif not self.jefferiesprior:
            #Using the empirical data to get the prior matrix
            self.logger.info('Filtering to obtain "best set"...')
            best_set = defaultdict(array)
            best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardminreads_and(incexcpairs=[[all_junctions["exc1"], all_junctions["inc1"]], [all_junctions["exc2"], all_junctions["inc2"]]], minreads=self.minandreads, logger=self.logger)
            best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardlow(self.minnonzero, True, self.logger, best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"])
            best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardminreads(self.minreads, True, self.logger, False, best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"])
            best_exc_reads1 = mean_junction(best_set['exc1'])
            best_inc_reads1 = mean_junction(best_set['inc1'])
            best_exc_reads2 = mean_junction(best_set['exc2'])
            best_inc_reads2 = mean_junction(best_set['inc2'])
            self.logger.info("'Best set' is %s events (out of %s)"%(best_set["inc1"].shape[0], all_junctions["inc1"].shape[0]))
            self.logger.info("Calculating PSI for 'best set'...")
            best_psi1 = simple_psi(best_inc_reads1, best_exc_reads1)
            best_psi2 = simple_psi(best_inc_reads2, best_exc_reads2)
            self.logger.info("Calculating delta PSI for 'best set'...")
            best_delta_psi = array(best_psi1 - best_psi2)

            self.logger.info("Parametrizing 'best set'...")
            mixture_pdf = adjustdelta(best_delta_psi, output, plotpath=self.plotpath, title=" ".join(self.names), numiter=self.iter, breakiter=self.breakiter, V=self.V)

            pickle.dump(mixture_pdf, open("%s%s_%s_bestset.pickle"%(output, self.names[0], self.names[1]), 'w'))
        
            prior_matrix = []
            for i in xrange(numbins):
                prior_matrix.extend(mixture_pdf[numbins-i:(numbins*2)-i])
            prior_matrix = array(prior_matrix).reshape(numbins, -1)

        #some info for later analysis
        pickle.dump(event_names, open("%s%s_%s_eventnames.pickle"%(output, self.names[0], self.names[1]), 'w')) 
        if not self.jefferiesprior:
            plot_matrix(prior_matrix, "Prior Matrix (before Jefferies)", "prior_matrix_no_jefferies", self.plotpath)

        #Calculate prior matrix
        self.logger.info("Adding a Jefferies prior to prior (alpha=%s)..."%(self.alpha))
        #Normalize prior with jefferies
        if self.jefferiesprior:
            self.logger.info("Using the Uniform distribution + Jefferies...")
            prior_matrix = jefferies + (self.prioruniform/numbins)
        else: 
            prior_matrix *= jefferies 

        prior_matrix /= sum(prior_matrix) #renormalize so it sums 1
        
        plot_matrix(prior_matrix, "Prior Matrix", "prior_matrix", self.plotpath)

        self.logger.info("Saving prior matrix for %s..."%(self.names))
        pickle.dump(prior_matrix, open("%s%s_%s_priormatrix.pickle"%(output, self.names[0], self.names[1]), 'w'))

        self.logger.info("Bootstrapping for all samples...")
        mean_exc1, var_exc1, exc_samples1 = sample_from_junctions(all_junctions["exc1"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc1, debug=self.debug)
        mean_inc1, var_inc1, inc_samples1 = sample_from_junctions(all_junctions["inc1"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc1, debug=self.debug)      
        mean_exc2, var_exc2, exc_samples2 = sample_from_junctions(all_junctions["exc2"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc2, debug=self.debug)
        mean_inc2, var_inc2, inc_samples2 = sample_from_junctions(all_junctions["inc2"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc2, debug=self.debug)

        self.logger.info("\nCalculating PSI (just for reference) for %s and %s..."%(self.names[0], self.names[1]))
        psi1 = calc_psi(inc_samples1, exc_samples1, self.names[0], self.alpha, self.n, self.debug, self.psiparam)
        psi2 = calc_psi(inc_samples2, exc_samples2, self.names[1], self.alpha, self.n, self.debug, self.psiparam)        
        pickle.dump(psi1, open("%s%s_psipaired.pickle"%(output, self.names[0]), 'w'))
        pickle.dump(psi2, open("%s%s_psipaired.pickle"%(output, self.names[1]), 'w'))

        self.logger.info("Calculating P(Data | PSI_i, PSI_j)...")
        #P(Data | PSI_i, PSI_j) = P(vector_i | PSI_i) * P(vector_j | PSI_j)
        data_given_psi1 = reads_given_psi(inc_samples1, exc_samples1, psi_space)
        data_given_psi2 = reads_given_psi(inc_samples2, exc_samples2, psi_space)

        data_given_psi = []
        for sample in xrange(data_given_psi1.shape[0]):
            #TODO Tensor product is calculated with scipy.stats.kron. Probably faster, have to make sure I am using it correctly.
            data_given_psi.append(data_given_psi1[sample].reshape(-1, numbins) * data_given_psi2[sample].reshape(numbins, -1))
            if self.debug: plot_matrix(data_given_psi[sample], "P(Data | PSI 1, PSI 2) Event %s (Inc1: %s, Exc1: %s Inc2: %s Exc2: %s)"%(sample, sum(inc_samples1[sample]), sum(exc_samples1[sample]), sum(inc_samples2[sample]), sum(exc_samples2[sample])), "datagpsi_sample%s"%sample, self.plotpath)

        #Finally, P(PSI_i, PSI_j | Data) equivalent to P(PSI_i, PSI_j)* P(Data | PSI_i, PSI_j) 
        self.logger.info("Calculate Posterior Delta Matrices...")
        posterior_matrix = []
        for sample in xrange(len(data_given_psi)):
            pm = (prior_matrix * data_given_psi[sample])
            posterior_matrix.append(pm / sum(pm))
            if self.debug: plot_matrix(posterior_matrix[-1], "Posterior Delta Event %s (Inc1: %s, Exc1: %s Inc2: %s Exc2: %s)"%(sample, sum(inc_samples1[sample]), sum(exc_samples1[sample]), sum(inc_samples2[sample]), sum(exc_samples2[sample])), "postdelta_sample%s"%sample, self.plotpath)

        pickle_path = "%s%s_%s_deltamatrix.pickle"%(output, self.names[0], self.names[1])
        pickle.dump(posterior_matrix, open(pickle_path, 'w'))
        self.logger.info("Done!")
        return posterior_matrix, event_names

def deltagroup(args):
    _pipeline_run(DeltaGroup(args))


class DeltaGroup(DeltaPair, CalcPsi):

    def calc_weights(self, replicas, relevant=None):
        self.logger.info("Loading data...")
        paired_replicas, event_names = load_data_n(replicas, logger=self.logger)
        self.logger.info("Calculating weights for %s..."%replicas)        

        paired_replicas = array(paired_replicas)
        paired_replicas_norm = []
        #gc content norm and masking
        for replica_num in xrange(0, len(replicas)*2, 2): 
            inc, exc = array(paired_replicas[:,replica_num]), array(paired_replicas[:,replica_num+1])
            self.logger.info("WEIGHTS: Calculating PSI for 'High coverage'...")
            all_junctions = {"inc": inc, "exc": exc }
            all_junctions = self.gc_content_norm(all_junctions)
            all_junctions = self.mark_stacks(all_junctions, fitfunc)
            self.logger.info("Masking non unique...")
            for junc_set in all_junctions.keys():
                all_junctions[junc_set] = masked_less(all_junctions[junc_set], 0) 

            paired_replicas_norm.append([all_junctions["exc"], all_junctions["inc"]])

        paired_replicas_norm = array(paired_replicas_norm)
        print "PAIRED", paired_replicas_norm.shape
        #filtering
        self.logger.info("WEIGHTS: Calculating local weights...")
        psis = []
        for exc, inc in paired_replicas_norm: 
            #self.logger.info("WEIGHTS: 'High coverage' is %s events (out of %s)"%(high_inc.shape[0], inc.shape[0]))
            mean_exc, var_exc, exc_samples = sample_from_junctions(exc, self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=poly1d([1,0]), debug=self.debug)
            mean_inc, var_inc, inc_samples = sample_from_junctions(inc, self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=poly1d([1,0]), debug=self.debug)      
            psis.append(calc_psi(exc_samples, inc_samples, None, self.alpha, self.n, self.debug, self.psiparam))
            
        if relevant: #if we have a relevant set, calculate weigths from the delta
            filtered_psis = defaultdict(list) #to hold the filtered PSI values in the relevant "best changing" set
            filtered_psis_paired = defaultdict(list) #to calculate the median PSI

            for i, experiment in enumerate(psis):
                for event_num, event in enumerate(experiment):
                    event_name = event_names[event_num]
                    if event_name in relevant: #all the psis in the relevant set
                        filtered_psis_paired[event_name].append(event)
                        filtered_psis[i].append(list(event))

            #calculate median PSI
            median_ref = [] #will hold the median values for each median reference   
            for event_name, events in filtered_psis_paired.items():
                median_ref.append(median(array(events), axis=0))
                median_ref[-1] /= sum(median_ref[-1])
                #print "#%s EVENTS"%len(events), events
                #print "MEDIAN", median_ref[-1], sum(median_ref[-1]), '\n'

            filter_lw = local_weights(array(filtered_psis.values()), self.weightsL1, array(median_ref))

        else: 
            lw = local_weights(psis)
            self.logger.info("WEIGHTS: Filtering local weights...")
            max_diffs = []
            #calculate the maximum difference between the weights
            #Also filter out the events that are not present in all experiments
            for event_weights in rollaxis(lw, 1):
                max_diff = 0
                #print event_weights
                num_experiments = 0
                max_experiments = len(event_weights)**2 #maximum number of matrices per event
                for i in xrange(len(event_weights)):
                    for j in xrange(len(event_weights)):
                        num_experiments += 1
                        max_diff = max(max_diff, event_weights[i] - event_weights[j] )

                if num_experiments == max_experiments:
                    max_diffs.append(max_diff)
            
            #filter the weights with the greatest changing the percentile of distance length, but leaving the ones that change the most (as they could be stacks)
            mindiff = scoreatpercentile(max_diffs, self.changsetpercentile)  
            self.logger.info("WEIGHTS: Total number of events is %s..."%lw.shape[1])  
            filter_lw = []
            for k, diff in enumerate(max_diffs):
                if diff >= mindiff:
                    filter_lw.append(lw[:,k])

            filter_lw = array(filter_lw)
            filter_lw = filter_lw.transpose(1, 0)
            self.logger.info("... total used for weights calculation is %s (%s %%) events that change the most"%(filter_lw.shape[1], self.changsetpercentile))
            
        gweights = global_weights(locweights=filter_lw)
        lweights_path = "%s_%s_%s_localweights.pickle"%(self.output, self.names[0], self.names[1])
        gweights_path = "%s_%s_%s_globalweights.pickle"%(self.output, self.names[0], self.names[1])
        pickle.dump(filter_lw,  open(lweights_path, 'w'))
        pickle.dump(gweights, open(gweights_path, 'w'))
        return gweights


    def _sub_calcweights(self, files, fixweights, relevant):
        if fixweights:
            weights = fixweights
        else:
            weights = self.calc_weights(files, relevant=relevant)

        self.logger.info("Weights for %s are %s (respectively)"%(files, weights))
        return weights

    def equal_if_not(self, weights):
        "If weigths havent been set, automatically fix them equally"
        if not weights.size:
            numpairs = len(self.k_ref)
            return [1./numpairs]*numpairs

        return weights

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

    def run(self):
        self.logger.info("")
        #calculate weights for both sets
        if self.replicaweights or self.fixweights1:
            weights1 = self._sub_calcweights(self.files1, self.fixweights1)
            weights2 = self._sub_calcweights(self.files2, self.fixweights2)

        #NOTE: global_weights better
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
                path = "%s%s_%s_"%(self.output, i, j)
                matrix_path = "%s%s_%s_deltamatrix.pickle"%(path, self.names[0], self.names[1])
                events_path = "%s%s_%s_eventnames.pickle"%(path, self.names[0], self.names[1]) 
                if os.path.exists(matrix_path):
                    self.logger.info("%s exists! Loading..."%path)
                    matrices = pickle.load(open(matrix_path))
                    names = pickle.load(open(events_path))
                else:
                    self.logger.info("Calculating pair %s"%path)
                    matrices, names = self.pairdelta(self.files1[i], self.files2[j], path)
                    self.logger.info("Saving pair posterior for %s, %s"%(i, j))

                if not self.replicaweights and not self.fixweights1: #get relevant events for weights calculation
                    relevant_events.extend(rank_deltas(matrices, names, E=True)[:self.numbestchanging])

                for k, name in enumerate(names):
                    pairs_posteriors[name].append(matrices[k]) #pairing all runs events

        #sort again the combined ranks
        relevant_events.sort(key=lambda x: -x[1])

        #gather the first NUMEVENTS names
        relevant = []
        for name, matrix in relevant_events:
            if name not in relevant_events: #ignore duplicated entries
                relevant.append(name)

            if len(relevant) == self.numbestchanging:
                break #we got enough elements

        #print relevant[:20], "...", relevant[-20:]
        if not self.replicaweights and not self.fixweights1:
            #relevant = self.get_relevant(pairs_posteriors)
            self.logger.info("Obtaining weights from Delta PSI...")
            weights1 = self._sub_calcweights(self.files1, None, relevant=relevant)
            weights2 = self._sub_calcweights(self.files2, None, relevant=relevant)

        self.logger.info("Normalizing with weights...")
        comb_matrix, comb_names = self.comb_replicas(pairs_posteriors, weights1=weights1, weights2=weights2)        
        self.logger.info("%s events matrices calculated"%len(comb_names))
        pickle_path = "%s%s_%s_deltacombmatrix.pickle"%(self.output, self.names[0], self.names[1])
        name_path = "%s%s_%s_combeventnames.pickle"%(self.output, self.names[0], self.names[1])
        pickle.dump(comb_matrix, open(pickle_path, 'w'))
        pickle.dump(comb_names, open(name_path, 'w'))
        self.logger.info("Alakazam! Done.")



