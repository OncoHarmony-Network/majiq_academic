import os
from collections import defaultdict
import abc
import pickle

from pylab import *
from numpy.ma import masked_less
from scipy.io import loadmat
from grimoire.utils.utils import create_if_not_exists, get_logger

from analysis.polyfitnb import fit_nb 
from analysis.filter import norm_junctions, discardlow, discardhigh, discardminreads, discardmaxreads, mark_stacks
from analysis.sample import sample_from_junctions
from analysis.psi import calc_psi, mean_psi, DirichletCalc, reads_given_psi, BINS_CENTER


################################
# Data loading and Boilerplate #
################################

def _load_data_delta3(my_mat, name, my_type):
    return {"junctions":my_mat[name][my_type][0, 0][0, 0]['cov'], "gc_content": my_mat[name][my_type][0, 0][0, 0]['gc_val']}

def _load_data_delta2(my_mat, name):
    inc = _load_data_delta3(my_mat, name, 'Inc')
    exc = _load_data_delta3(my_mat, name, 'Exc')
    const = _load_data_delta3(my_mat, name, 'rand10k')
    return inc, exc, const

def load_data_delta(path, name1, name2):
    "Load data from a matlab file. Should be deprecated soon"
    my_mat = loadmat(path)
    inc1, exc1, const1 = _load_data_delta2(my_mat, name1)
    inc2, exc2, const2 = _load_data_delta2(my_mat, name2)
    return inc1, exc1, const1, inc2, exc2, const2


def _load_data_psi2(my_mat, my_type):
    return {"junctions":my_mat[my_type][0, 0]['cov'], "gc_content": my_mat[my_type][0, 0]['gc_val']}


def load_data_psi(path):
    my_mat = loadmat(path)
    return _load_data_psi2(my_mat, 'Inc'), _load_data_psi2(my_mat, 'Exc'), _load_data_psi2(my_mat, 'rand10k')


def _pipeline_run(pipeline):
    "Exception catching for all the pipelines"
    try:
        pipeline.run()

    except KeyboardInterrupt:
        if pipeline.logger: pipeline.logger.info("MAJIQ manually interrupted. Exiting...")


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
            self.logger.info("Processing %s..."%path)
            inc, exc, const = load_data_psi(path) 
            all_junctions = {"inc": inc, "exc": exc, "const": const }

            self.logger.info("GC content normalization...")
            for junc_set in all_junctions.keys():
                all_junctions[junc_set] = norm_junctions(all_junctions[junc_set]["junctions"], 
                                                         all_junctions[junc_set]["gc_content"])

            self.logger.info("Masking non unique...")
            for junc_set in all_junctions.keys():
                all_junctions[junc_set] = masked_less(all_junctions[junc_set], 0) 

            if self.debug:
                self.logger.debug("Skipping fitfunc because --debug!")
                fitfunc = poly1d([1, 0])

            else:
                self.logger.info("Fitting NB function with constitutive events...")
                fitfunc = fit_nb(all_junctions["const"], "%s_nbfit"%self.output, self.plotpath, nbdisp=self.nbdisp, logger=self.logger)

            if self.markstacks >= 0:
                self.logger.info("Marking and masking stacks...")
                for junc_set in all_junctions.keys():
                    if junc_set.find("const") == -1:
                        print "... %s"%junc_set
                        all_junctions[junc_set] = mark_stacks(all_junctions[junc_set], fitfunc, self.markstacks, self.nbdisp)


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
            output_path = open("%s%s_psi.pickle"%(self.output, name), 'w')
            pickle.dump(psi, output_path)
            self.logger.info("PSI calculation for %s ended succesfully! Result "%output_path)
            return psi


################################
#          Delta PSI           #
################################


def deltapair(args):
    _pipeline_run(DeltaPair(args))

class DeltaPair(BasicPipeline):

    def run(self):
        name1, name2 = ((self.file1.split("/")[-1]).split('.mat'))[0].split('_') #Esto asin no, cambiar cuando venga lo nuevo
        self.logger.info("")
        self.logger.info("Processing pair %s..."%self.file1)
        inc1, exc1, const1, inc2, exc2, const2 = load_data_delta(self.file1, name1, name2) #loading the paired matrixes
        all_junctions = {"inc1": inc1, "exc1": exc1, "const1": const1, 
                         "inc2": inc2, "exc2": exc2, "const2": const2 }

        self.logger.info("GC content normalization...")
        for junc_set in all_junctions.keys():
            all_junctions[junc_set] = norm_junctions(all_junctions[junc_set]["junctions"], 
                                                     all_junctions[junc_set]["gc_content"])

        self.logger.info("Masking non unique...")
        for junc_set in all_junctions.keys():
            all_junctions[junc_set] = masked_less(all_junctions[junc_set], 0) 

        if self.debug:
            logger.debug("Skipping fitfunc because --debug!!")
            fitfunc1 = poly1d([1, 0])
            fitfunc2 = poly1d([1, 0])
        else:
            self.logger.info("Fitting NB function with constitutive events...")
            fitfunc1 = fit_nb(all_junctions["const1"], "%s_nbfit1"%self.output, self.plotpath, nbdisp=self.nbdisp)
            fitfunc2 = fit_nb(all_junctions["const2"], "%s_nbfit2"%self.output, self.plotpath, nbdisp=self.nbdisp)

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


        #TODO This should be happening for every pair of N/M experiments (for N in exp: for M in exp:)
        self.logger.info('Filtering to obtain "best set"...')
        best_set = defaultdict(array)

        best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardminreads_and(incexcpairs=[[all_junctions["inc1"], all_junctions["exc1"]], [all_junctions["inc2"], all_junctions["exc2"]]], minreads=minandreads, logger=False)
        best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardlow(self.minnonzero, True, logger, best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"])
        best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardminreads(self.minreads, True, logger, best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"])
         
        self.logger.info("Bootstrapping for 'best set'...")
        mean_exc1, var_exc1, exc_best_samples1 = sample_from_junctions(best_set["exc1"], self.m, self.k, discardzeros=False, trimborder=self.trimborder, fitted_func=fitfunc1, debug=self.debug)
        mean_inc1, var_inc1, inc_best_samples1 = sample_from_junctions(best_set["inc1"], self.m, self.k, discardzeros=False, trimborder=self.trimborder, fitted_func=fitfunc1, debug=self.debug)      
        mean_exc2, var_exc2, exc_best_samples2 = sample_from_junctions(best_set["exc2"], self.m, self.k, discardzeros=False, trimborder=self.trimborder, fitted_func=fitfunc2, debug=self.debug)
        mean_inc2, var_inc2, inc_best_samples2 = sample_from_junctions(best_set["inc2"], self.m, self.k, discardzeros=False, trimborder=self.trimborder, fitted_func=fitfunc2, debug=self.debug)

        self.logger.info("\nCalculating PSI for 'best set' %s ..."%(self.names))
        best_psi1 = calc_psi(inc_best_samples1, exc_best_samples1, name1, self.alpha, self.n, self.debug, self.psiparam)
        best_psi2 = calc_psi(inc_best_samples2, exc_best_samples2, name2, self.alpha, self.n, self.debug, self.psiparam)
        
        self.logger.info("\nCalculating delta PSI for 'best set' %s ..."%(self.names))
        best_delta_psi = mean_psi(best_psi1) - mean_psi(best_psi2)

        self.logger.info("Obtaning prior matrix for 'best set'...")
        mixture_pdf = adjustdelta(best_delta_psi, self.output, plotpath=self.plotpath, title=" ".join(self.names), numiter=self.iter, breakiter=self.breakiter, V=self.V)

        self.logger.info("Calculating prior matrix...")
        numbins = int(round(len(mixture_pdf)/2)) #half the delta bins
        dircalc = DirichletCalc() 
        #Calculate prior matrix
        prior_matrix = []
        for i in xrange(numbins):
            prior_matrix.extend(mixture_pdf[numbins-i:(numbins*2)-i])

        prior_matrix = array(prior_matrix).reshape(numbins, -1)

        clf()
        title("Prior Matrix")
        imshow(prior_matrix)
        xlabel("PSI_i")
        ylabel("PSI_j")
        _save_or_show(self.plotpath, plotname="prior_matrix_no_jefferies")

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

        clf()
        title("Jefferies Matrix")
        imshow(jefferies)
        xlabel("PSI_i")
        ylabel("PSI_j")
        _save_or_show(self.plotpath, plotname="jefferies_matrix")

        prior_matrix *= jefferies #Normalize PSI with jefferies
        prior_matrix /= sum(prior_matrix) #renormalize so it sums 1

        clf()
        title("Prior Matrix")
        imshow(prior_matrix)
        xlabel("PSI 1")
        ylabel("PSI 2")
        _save_or_show(self.plotpath, plotname="prior_matrix")

        self.logger.info("Saving prior matrix for %s..."%(self.names))
        pickle.dump(prior_matrix, open("%s%s_%s_priormatrix.pickle"%(self.output, name1, name2), 'w'))

        self.logger.info("Bootstrapping for all samples...")
        mean_exc1, var_exc1, exc_samples1 = sample_from_junctions(all_junctions["exc1"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc1, debug=self.debug)
        mean_inc1, var_inc1, inc_samples1 = sample_from_junctions(all_junctions["inc1"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc1, debug=self.debug)      
        mean_exc2, var_exc2, exc_samples2 = sample_from_junctions(all_junctions["exc2"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc2, debug=self.debug)
        mean_inc2, var_inc2, inc_samples2 = sample_from_junctions(all_junctions["inc2"], self.m, self.k, discardzeros=self.discardzeros, trimborder=self.trimborder, fitted_func=fitfunc2, debug=self.debug)
  
        self.logger.info("Writing samples in disk...")
        __write_samples(exc_samples1, self.output, self.names, 0)
        __write_samples(inc_samples1, self.output, self.names, 1)
        __write_samples(exc_samples2, self.output, self.names, 0)
        __write_samples(inc_samples2, self.output, self.names, 1)

        self.logger.info("\nCalculating PSI %s ..."%(self.names))
        psi1 = calc_psi(inc_samples1, exc_samples1, name1, self.alpha, self.n, self.debug, self.psiparam)
        psi2 = calc_psi(inc_samples2, exc_samples2, name2, self.alpha, self.n, self.debug, self.psiparam)

        pickle.dump(psi1, open("%s%s_psi.pickle"%(self.output, name1), 'w'))
        pickle.dump(psi1, open("%s%s_psi.pickle"%(self.output, name2), 'w'))

        self.logger.info("Calculating P(Data | PSI_i, PSI_j)...")
        #P(Data | PSI_i, PSI_j) = P(vector_i | PSI_i) * P(vector_j | PSI_j)
        data_given_psi1 = reads_given_psi(inc_samples1, exc_samples1, psi_space)
        data_given_psi2 = reads_given_psi(inc_samples2, exc_samples2, psi_space)


        data_given_psi = []
        for sample in xrange(data_given_psi1.shape[0]):
            #TODO Tensor product is calculated with scipy.stats.kron. Probably faster, have to make sure I am using it correctly.
            data_given_psi.append(data_given_psi1[sample].reshape(-1, numbins) * data_given_psi2[sample].reshape(numbins, -1))
            print "Event %s"%sample
            print inc_samples1[sample], exc_samples1[sample]
            print
            print sum(data_given_psi[sample])
            clf()
            title("P(Data | PSI 1, PSI 2) Event %s (Inc1: %s, Exc1: %s Inc2: %s Exc2: %s)"%(sample, sum(inc_samples1[sample]), sum(exc_samples1[sample]), sum(inc_samples2[sample]), sum(exc_samples2[sample])))
            imshow(data_given_psi[sample])
            xlabel("PSI 1")
            ylabel("PSI 2")
            _save_or_show(self.plotpath, plotname="datagpsi_sample%s"%sample)


        #Finally, P(PSI_i, PSI_j | Data) proportional P(PSI_i, PSI_j)* P(Data | PSI_i, PSI_j) 
        self.logger.info("Calculate Posterior Delta Matrices...")
        posterior_matrix = []
        for sample in xrange(len(data_given_psi)):
            pm = (prior_matrix * data_given_psi[sample])
            posterior_matrix.append(pm / sum(pm))

            clf()
            title("Posterior Delta Event %s ()"%sample)
            imshow(posterior_matrix[-1])
            xlabel("PSI i")
            ylabel("PSI j")
            _save_or_show(self.plotpath, plotname="postdelta_sample%s"%sample)

        pickle.dump(posterior_matrix, open("%s%s_%s_deltamatrix.pickle"%(self.output, name1, name2), 'w'))
        self.logger.info("Done!")






def deltagroup(args):
    pass