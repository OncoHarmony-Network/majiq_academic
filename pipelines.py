import os
from collections import defaultdict
import abc
import pickle

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
#from analysis.filter import norm_junctions, discardlow, discardhigh, discardminreads, discardmaxreads, discardminreads_and, mark_stacks

################################
# Data loading and Boilerplate #
################################





def _pipeline_run(pipeline, lsv=False, logger=None):
    "Exception catching for all the pipelines"
    try:
        pipeline.run(lsv)

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
        self.lsv = args.lsv
        self.psi_paths = []
        
        try:
            self.replica_len = [len(self.files1), len(self.files2)]
        except AttributeError:
            pass


    @abc.abstractmethod
    def run(self, lsv):
        """This is the entry point for all pipelines"""
        return

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
            return fit_nb(const_junctions, "%s_nbfit"%self.output, self.plotpath, nbdisp=self.nbdisp, logger=self.logger, discardb=True, bval=False)

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

################################
# PSI calculation pipeline     #
################################

def calcpsi(args):
    _pipeline_run(CalcPsi(args), args.lsv)

class CalcPsi(BasicPipeline):

    def run(self, lsv=False):
        for path in self.files:
            print "LSV",lsv
            if lsv :
                self.calcpsi_lsv(path)
            else:
                self.calcpsi(path)

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
#        all_junctions = {"inc": inc, "exc": exc, "const": const }
#        all_junctions = self.gc_content_norm(all_junctions)

       # self.logger.info("Masking non unique...")


        #for junc_set in all_junctions.keys():
        #    all_junctions[junc_set] = masked_less(all_junctions[junc_set], 0) 

        fitfunc = self.fitfunc(const[0])
        filter_lsv = self.mark_stacks_lsv( lsv_junc[0], fitfunc)
        #FILTER_JUNCTIONS?
        self.logger.info('Filtering ...')
        lsv_junc = majiq_filter.lsv_quantifiable( filter_lsv, self.minnonzero, self.minreads, self.logger )

        self.logger.info("Bootstrapping samples...") 
        lsv_sample = []
        for ii in lsv_junc[0]:
            m_lsv, var_lsv, s_lsv = sample_from_junctions(ii, self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc, debug=self.debug) 
            lsv_sample.append( s_lsv )

        self.logger.info("\nCalculating PSI for %s ..."%(name))
        psi = lsv_psi(lsv_sample, name, self.alpha, self.n, self.debug)

        self.logger.info("Saving PSI...")
        if write_pickle:
            output = open("%s%s_psi.pickle"%(self.output, name), 'w')
            pickle.dump(psi, output)
            self.logger.info("PSI calculation for %s ended succesfully! Result can be found at %s"%(name, output.name))

        return psi


    def calcpsi(self, path, write_pickle=True):
        """
        Given a file path with the junctions, return psi distributions. 
        write_pickle indicates if a .pickle should be saved in disk
        """
        name = os.path.basename(path)
        self.logger.info("")
        self.logger.info("Loading %s..."%path)
        inc, exc, const = majiq_io.load_data(path, self.logger) 
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
        filter_junctions["exc"], filter_junctions["inc"] = majiq_filter.discardlow(self.minnonzero, True, self.logger, None, all_junctions["exc"], all_junctions["inc"])
        filter_junctions["exc"], filter_junctions["inc"] = majiq_filter.discardminreads(self.minreads, True, self.logger, False, None, filter_junctions["exc"], filter_junctions["inc"])

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
    _pipeline_run( DeltaPair(args), args.lsv )

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

        lsv_junc1, const1 = majiq_io.load_data_lsv(path, self.logger) 
        lsv_junc2, const2 = majiq_io.load_data_lsv(path, self.logger) 

#        self.logger.debug("Shapes for inclusion, %s exclusion, %s constitutive %s"%(inc1.shape, exc1.shape, const1.shape))
#        all_junctions = {"inc1": inc1, "exc1": exc1, "const1": const1, "inc2": inc2, "exc2": exc2, "const2": const2 }
#        all_junctions = self.gc_content_norm(all_junctions)
        self.logger.info("Masking non unique...")
        for junc_set in all_junctions.keys():
            #print junc_set, all_junctions[junc_set], all_junctions[junc_set].shape
            all_junctions[junc_set] = masked_less(array(all_junctions[junc_set]), 0) 

        #fitting the function
        fitfunc1 = self.fitfunc(const1[0])
        fitfunc2 = self.fitfunc(const2[0])
        filtered_lsv1 = self.mark_stacks_lsv( lsv_junc1[0], fitfunc1)
        filtered_lsv2 = self.mark_stacks_lsv( lsv_junc1[0], fitfunc1)

        #Quantifiable junctions filter
        ''' Quantify and unify '''
        self.logger.info('Filtering ...')
#        all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"], event_names = majiq_filter.discardlow(self.minnonzero, True, self.logger, event_names, all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"],)
#        all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"], event_names = majiq_filter.discardminreads(self.minreads, True, self.logger, False, event_names, all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"],)

        psi_space, prior_matrix = majiq_psi.gen_prior_matrix()

        self.logger.info("Bootstrapping for all samples...")
        lsv_sample1 = []
        for ii in lsv_junc1[0]:
            m_lsv, var_lsv, s_lsv = sample_from_junctions(ii, self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc, debug=self.debug) 
            lsv_sample1.append( s_lsv )
        lsv_sample2 = []
        for ii in lsv_junc2[0]:
            m_lsv, var_lsv, s_lsv = sample_from_junctions(ii, self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc, debug=self.debug) 
            lsv_sample2.append( s_lsv )

        ''' This should be debug code '''
#        self.logger.info("\nCalculating PSI (just for reference) for %s and %s..."%(self.names[0], self.names[1]))
#        psi1 = calc_psi(inc_samples1, exc_samples1, self.names[0], self.alpha, self.n, self.debug, self.psiparam)
#        psi2 = calc_psi(inc_samples2, exc_samples2, self.names[1], self.alpha, self.n, self.debug, self.psiparam) 
#        psi_path = "%s%s_%s_psipaired.pickle"%(output, self.names[0], self.names[1])
#        pickle.dump([psi1, psi2], open(psi_path, 'w'))

        self.logger.info("Calculating P(Data | PSI_i, PSI_j)...")
        #P(Data | PSI_i, PSI_j) = P(vector_i | PSI_i) * P(vector_j | PSI_j)

        #for num_psi:
        #    data_given_psi1 = reads_given_psi(inc_samples1, exc_samples1, psi_space)
        #    data_given_psi2 = reads_given_psi(inc_samples2, exc_samples2, psi_space)

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








    def pairdelta(self, file1, file2, output):
        self.logger.info("")
        self.logger.info("Processing pair %s - %s..."%(file1, file2))
        inc1, exc1, const1, inc2, exc2, const2, event_names = majiq_io.load_data_pair(file1, file2, self.logger, self.tracklist) 
        self.logger.debug("Shapes for inclusion, %s exclusion, %s constitutive %s"%(inc1.shape, exc1.shape, const1.shape))
        all_junctions = {"inc1": inc1, "exc1": exc1, "const1": const1, "inc2": inc2, "exc2": exc2, "const2": const2 }
        all_junctions = self.gc_content_norm(all_junctions)
        self.logger.info("Masking non unique...")
        for junc_set in all_junctions.keys():
            #print junc_set, all_junctions[junc_set], all_junctions[junc_set].shape
            all_junctions[junc_set] = masked_less(array(all_junctions[junc_set]), 0) 

        #Quantifiable junctions filter
        self.logger.info('Filtering ...')
        all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"], event_names = majiq_filter.discardlow(self.minnonzero, True, self.logger, event_names, all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"],)
        all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"], event_names = majiq_filter.discardminreads(self.minreads, True, self.logger, False, event_names, all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"],)

        #fitting the function
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

                    if self.logger: self.logger.info("... %s"%junc_set)
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

        if self.tracklist:
            self.logger.info("TRACKLIST: After filters")
            for i, event_name in enumerate(event_names):
                if event_name in self.tracklist:
                    logger.info("TRACKLIST After filters (%s): %s"%(event_name, all_junctions['inc1'][i], all_junctions['exc1'][i], all_junctions['exc2'][i], all_junctions['inc2'][i]))

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
            best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = majiq_filter.discardminreads_and(incexcpairs=[[all_junctions["exc1"], all_junctions["inc1"]], [all_junctions["exc2"], all_junctions["inc2"]]], minreads=self.priorminandreads, logger=self.logger)
            best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = majiq_filter.discardlow(self.priorminnonzero, True, self.logger, None, best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"])
            best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = majiq_filter.discardminreads(self.priorminreads, True, self.logger, False, None, best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"])
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
            mixture_pdf = adjustdelta(best_delta_psi, output, plotpath=self.plotpath, title=" ".join(self.names), numiter=self.iter, breakiter=self.breakiter, V=self.V, logger=self.logger)

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
        mean_exc1, var_exc1, exc_samples1 = sample_from_junctions(all_junctions["exc1"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc1, debug=self.debug, tracklist=self.tracklist, names=event_names)
        mean_inc1, var_inc1, inc_samples1 = sample_from_junctions(all_junctions["inc1"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc1, debug=self.debug, tracklist=self.tracklist, names=event_names)      
        mean_exc2, var_exc2, exc_samples2 = sample_from_junctions(all_junctions["exc2"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc2, debug=self.debug, tracklist=self.tracklist, names=event_names)
        mean_inc2, var_inc2, inc_samples2 = sample_from_junctions(all_junctions["inc2"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc2, debug=self.debug, tracklist=self.tracklist, names=event_names)
        self.logger.info("\nCalculating PSI (just for reference) for %s and %s..."%(self.names[0], self.names[1]))
        psi1 = calc_psi(inc_samples1, exc_samples1, self.names[0], self.alpha, self.n, self.debug, self.psiparam)
        psi2 = calc_psi(inc_samples2, exc_samples2, self.names[1], self.alpha, self.n, self.debug, self.psiparam) 
        psi_path = "%s%s_%s_psipaired.pickle"%(output, self.names[0], self.names[1])
        pickle.dump([psi1, psi2], open(psi_path, 'w'))

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
                    relevant_events.extend(rank_deltas(matrices, names, E=True)[:self.numbestchanging])

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



