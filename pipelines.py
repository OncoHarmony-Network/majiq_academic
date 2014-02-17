import os
from collections import defaultdict
import abc
import pickle

from pylab import *
from numpy.ma import masked_less
from scipy.io import loadmat
from scipy.stats import norm
from scipy.sparse import lil_matrix

from grimoire.utils.utils import create_if_not_exists, get_logger
from analysis.polyfitnb import fit_nb 
from analysis.filter import norm_junctions, discardlow, discardhigh, discardminreads, discardmaxreads, discardminreads_and, mark_stacks
from analysis.sample import sample_from_junctions, mean_junction
from analysis.psi import calc_psi, mean_psi, simple_psi, DirichletCalc, reads_given_psi, BINS_CENTER
from analysis.adjustdelta import adjustdelta

################################
# Data loading and Boilerplate #
################################

DELTA_RATIO = 0.2 #TODO to parameters


def _find_len(grimoire_obj, logger):
    junc_len = 666
    for junction in grimoire_obj:
        if hasattr(junction[0], 'coverage'):
            junc_len = junction[0].coverage[0].shape[0]
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

        ret.append(list(junction[0].coverage[0]))

    return array(ret)

def _load_data2(grimoire_obj, logger=None, getnames = False):
    """
    Same as above
    """
    ret = []
    #first get junction length, should be almost always the first entry
    _find_len(grimoire_obj, logger)

    for junction in grimoire_obj:   
        if hasattr(junction, 'coverage'):
            ret.append(list(junction.coverage[0]))
            #print "Empirical", len(ret[-1])
        else:
            ret.append([0]*junc_len)
            #print "Theoretical", len(ret[-1])

    return array(ret)

def load_data(path, logger=None):
    "Load data from the preprocess step. Could change to a DDBB someday"
    data = pickle.load(open(path))
    return _load_data2(data[1][:,0], logger), _load_data2(data[1][:,1], logger), _load_data_const(data[2], logger) #inc, exc, const

def load_data_pair(path1, path2, logger=None):
    """Pairing functionality should be extracted of this function"""
    #const extracting doesnt change
    data1 = pickle.load(open(path1))
    data2 = pickle.load(open(path2))
    const1 = _load_data_const(data1[2], logger)
    const2 = _load_data_const(data2[2], logger)
    paired_samples = defaultdict(list)

    my_len = _find_len(data1[1], logger)

    for i, junction in enumerate(data1[1]): #iterate junctions in 1 and add to the inc and exc to the dictionary
        if hasattr(junction[0], 'coverage') and hasattr(junction[1], 'coverage'):
            my_inc = list(junction[0].coverage[0])
            my_exc = list(junction[1].coverage[0])

            if type(my_inc[0]) == lil_matrix:
                my_inc = [0]*my_len

            if type(my_exc[0]) == lil_matrix:
                my_exc = [0]*my_len

            my_name = junction[0].name
            if not my_name:
                my_name = junction[1].name

            paired_samples[my_name].extend([my_inc, my_exc])

    #same for the second pair
    for i, junction in enumerate(data2[1]): #iterate junctions in 1 and add to the inc and exc to the dictionary 
        if hasattr(junction[0], 'coverage') and hasattr(junction[1], 'coverage'):
            my_inc = list(junction[0].coverage[0])
            my_exc = list(junction[1].coverage[0])
            if type(my_inc[0]) == lil_matrix:
                my_inc = [0]*my_len

            if type(my_exc[0]) == lil_matrix:
                my_exc = [0]*my_len

            my_name = junction[0].name
            if not my_name:
                my_name = junction[1].name

            paired_samples[my_name].extend([my_inc, my_exc])

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
        if pipeline.logger: pipeline.logger.info("MAJIQ manually interrupted. Exiting...")

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
        self.logger.info("GC content normalization...")
        if self.gcnorm:
            for junc_set in all_junctions.keys():
                all_junctions[junc_set] = norm_junctions(all_junctions[junc_set]["junctions"], all_junctions[junc_set]["gc_content"])

        return all_junctions

################################
# PSI calculation pipeline     #
################################

def calcpsi(args):
    _pipeline_run(CalcPsi(args))

class CalcPsi(BasicPipeline):

    def run(self):
        for path in self.files:
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

            if self.debug:
                self.logger.debug("Skipping fitfunc because --debug!")
                fitfunc = poly1d([1, 0])
            else:
                self.logger.info("Fitting NB function with constitutive events...")
                fitfunc = fit_nb(all_junctions["const"], "%s_nbfit"%self.output, self.plotpath, nbdisp=self.nbdisp, logger=self.logger, discardb=True)

            if self.markstacks >= 0:
                self.logger.info("Marking and masking stacks...")
                for junc_set in all_junctions.keys():
                    if junc_set.find("const") == -1:
                        print "... %s"%junc_set
                        all_junctions[junc_set] = mark_stacks(all_junctions[junc_set], fitfunc, self.markstacks, self.nbdisp)

                    all_junctions[junc_set] = masked_less(all_junctions[junc_set], 0) #remask the stacks

            self.logger.info('Filtering ...')
            filter_junctions = defaultdict(array)
            filter_junctions["exc"], filter_junctions["inc"] = discardlow(self.minnonzero, True, self.logger, all_junctions["exc"], all_junctions["inc"])
            filter_junctions["exc"], filter_junctions["inc"] = discardminreads(self.minreads, True, self.logger, False, filter_junctions["exc"], filter_junctions["inc"])
             
            self.logger.info("Bootstrapping samples...")
            mean_exc, var_exc, exc_samples = sample_from_junctions(filter_junctions["exc"], self.m, self.k, discardzeros=False, trimborder=self.trimborder, fitted_func=fitfunc, debug=self.debug)
            mean_inc, var_inc, inc_samples = sample_from_junctions(filter_junctions["inc"], self.m, self.k, discardzeros=False, trimborder=self.trimborder, fitted_func=fitfunc, debug=self.debug)      

            self.logger.info("\nCalculating PSI for %s ..."%(name))
            psi = calc_psi(inc_samples, exc_samples, name, self.alpha, self.n, self.debug, self.psiparam)

            self.logger.info("Saving PSI...")
            output = open("%s%s_psi.pickle"%(self.output, name), 'w')
            pickle.dump(psi, output)
            self.logger.info("PSI calculation for %s ended succesfully! Result can be found at %s"%(name, output.name))


################################
#          Delta PSI           #
################################

def deltapair(args):
    _pipeline_run(DeltaPair(args))

class DeltaPair(BasicPipeline):

    def run(self):
        self.logger.info("")
        self.logger.info("Processing pair %s..."%self.file1)
        inc1, exc1, const1, inc2, exc2, const2, event_names= load_data_pair(self.file1, self.file2, self.logger) 
        self.logger.debug("SHAPES for inclusion, %s exclusion, %s constitutive %s"%(inc1.shape, exc1.shape, const1.shape))
        all_junctions = {"inc1": inc1, "exc1": exc1, "const1": const1, "inc2": inc2, "exc2": exc2, "const2": const2 }
        all_junctions = self.gc_content_norm(all_junctions)
        self.logger.info("Masking non unique...")
        for junc_set in all_junctions.keys():
            all_junctions[junc_set] = masked_less(all_junctions[junc_set], 0) 

        if self.debug:
            self.logger.debug("Skipping fitfunc because --debug!!")
            fitfunc1 = poly1d([1, 0])
            fitfunc2 = poly1d([1, 0])
        else:
            self.logger.info("Fitting NB function with constitutive events...")
            fitfunc1 = fit_nb(all_junctions["const1"], "%s_nbfit1"%self.output, self.plotpath, nbdisp=self.nbdisp, discardb=True)
            fitfunc2 = fit_nb(all_junctions["const2"], "%s_nbfit2"%self.output, self.plotpath, nbdisp=self.nbdisp, discardb=True)

            self.logger.info("Fitfunc 1: %s %s Fitfunc2: %s %s"%(fitfunc1.c[0], fitfunc1.c[1], fitfunc2.c[0], fitfunc2.c[1]))

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

        self.logger.info("Calculating prior matrix...")
        numbins = 20 #half the delta bins
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
            mydist = norm(loc=0, scale=self.priorstd)
            norm_space = linspace(-1, 1-self.binsize, num=numbins*2) + self.binsize/2
            pdfnorm = mydist.pdf(norm_space)
            uniform = self.prioruniform/numbins
            newdist = (pdfnorm+uniform)/(pdfnorm+uniform).sum()
            pdfnorm /= sum(newdist)
            print pdfnorm
            print pdfnorm.shape
            plot(linspace(-1, 1, num=len(list(pdfnorm))), pdfnorm)
            _save_or_show(self.plotpath, plotname="prior_distribution")
            for i in xrange(numbins):
                prior_matrix.append(list(pdfnorm[numbins-i:(numbins*2)-i]))

            prior_matrix = array(prior_matrix)
            print prior_matrix.shape
            print prior_matrix
            prior_matrix /= sum(prior_matrix) #renormalize so it sums 1
            plot_matrix(prior_matrix, "Prior Matrix (before Jefferies)", "prior_matrix_no_jefferies", self.plotpath)

        else:
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
            mixture_pdf = adjustdelta(best_delta_psi, self.output, plotpath=self.plotpath, title=" ".join(self.names), numiter=self.iter, breakiter=self.breakiter, V=self.V)

            pickle.dump(mixture_pdf, open("%s%s_%s_bestset.pickle"%(self.output, self.names[0], self.names[1]), 'w'))
        
            prior_matrix = []
            for i in xrange(numbins):
                prior_matrix.extend(mixture_pdf[numbins-i:(numbins*2)-i])

            prior_matrix = array(prior_matrix).reshape(numbins, -1)

        #some info for later analysis
        pickle.dump(event_names, open("%s%s_%s_eventnames.pickle"%(self.output, self.names[0], self.names[1]), 'w')) 
        plot_matrix(prior_matrix, "Prior Matrix (before Jefferies)", "prior_matrix_no_jefferies", self.plotpath)
        self.logger.info("Adding a uniform to prior (%.3f)..."%self.prioruniform)
        #prior_matrix /= sum(prior_matrix) #renormalize so it sums 1
        self.logger.info("Adding a Jefferies prior to prior (alpha=%s)..."%(self.alpha))

        prior_matrix *= jefferies #Normalize PSI with jefferies
        prior_matrix /= sum(prior_matrix) #renormalize so it sums 1
        
        plot_matrix(prior_matrix, "Prior Matrix", "prior_matrix", self.plotpath)

        self.logger.info("Saving prior matrix for %s..."%(self.names))
        pickle.dump(prior_matrix, open("%s%s_%s_priormatrix.pickle"%(self.output, self.names[0], self.names[1]), 'w'))

        self.logger.info("Bootstrapping for all samples...")
        mean_exc1, var_exc1, exc_samples1 = sample_from_junctions(all_junctions["exc1"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc1, debug=self.debug)
        mean_inc1, var_inc1, inc_samples1 = sample_from_junctions(all_junctions["inc1"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc1, debug=self.debug)      
        mean_exc2, var_exc2, exc_samples2 = sample_from_junctions(all_junctions["exc2"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc2, debug=self.debug)
        mean_inc2, var_inc2, inc_samples2 = sample_from_junctions(all_junctions["inc2"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc2, debug=self.debug)

        self.logger.info("\nCalculating PSI (just for reference) for %s and %s..."%(self.names[0], self.names[1]))
        psi1 = calc_psi(inc_samples1, exc_samples1, self.names[0], self.alpha, self.n, self.debug, self.psiparam)
        psi2 = calc_psi(inc_samples2, exc_samples2, self.names[1], self.alpha, self.n, self.debug, self.psiparam)        
        pickle.dump(prior_matrix, open("%s%s_%s_psipaired.pickle"%(self.output, self.names[0], self.names[1]), 'w'))
        pickle.dump(prior_matrix, open("%s%s_%s_psipaired.pickle"%(self.output, self.names[0], self.names[1]), 'w'))

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

        pickle.dump(posterior_matrix, open("%s%s_%s_deltamatrix.pickle"%(self.output, self.names[0], self.names[1]), 'w'))
        self.logger.info("Done!")



def deltagroup(args):
    pass