from __future__ import division
import scipy.io
import scipy
from collections import defaultdict
import numpy as np
from pylab import *
from matplotlib import pyplot
from scipy.stats import poisson
from scipy.stats.mstats import mquantiles

tissues = {'Heart':[0,1,2,3],'Hippocampus':[4,5,6,7],'Lung':[8,9,10,11],'Spleen':[12,13,14,15],'Thymus':[16,17,18,19,20]}
#tissues = {'Heart':[0,1],'Hippocampus':[2,3],'Lung':[4,5],'Spleen':[6,7],'Thymus':[8,9]}
tissues = {'Heart':[0]}
#files=['Heart1chr1']
#files = ['Heart1','Heart3', 'Hippocampus1','Hippocampus2', 'Lung1','Lung2', 'Spleen1','Spleen2', 'Thymus1','Thymus2']
files = ['Heart1','Heart3','Heart5','Heart6',
            'Hippocampus1','Hippocampus2','Hippocampus4','Hippocampus6',
            'Lung2','Lung3','Lung5','Lung6',
            'Spleen1','Spleen2','Spleen3','Spleen5',
            'Thymus1','Thymus2','Thymus3','Thymus4','Thymus6']

dr = '/data/ucsc/reads/test_1k/new_mat/'
fidx = 0
RPKMlist = [[] for ii in range(len(files))]
junxReads = [[] for ii in range(len(files))]
junxPos = [[0]*77 for ii in range(len(files))]
junxPos_nu = [[0]*77 for ii in range(len(files))]
pvalues = [np.array([],np.dtype('float')) for ii in range(len(files))]
lmbd = [[] for ii in range(len(files))]
global_junxr = {}
junxposhist = [[] for ii in range(len(files))]
junx_readatpos = [[0]*77 for ii in range(len(files))]
junx_readatpos_nu =  [[0]*77 for ii in range(len(files))]
counter = [0]*len(files)


for fl_idx, fl in enumerate(files):
    fil = dr+'test_'+fl
    print fl
    heart1 = scipy.io.loadmat( fil)
    readxpos_mean = [0]*77
    readxpos_var = [0]*77
    id2 = 0
    list_p = []
    p2 = []
    xpos = [[] for ii in range(77)]
    total_gen = len(heart1['RNASEQ_INFO'][0])
    for ii in heart1['RNASEQ_INFO'][0]:
        if id2 % 100 ==0: print "%d/%d"%(id2,total_gen)
#        print id2
        #print "GET RPKM ", id2
        if ii[0,0] == 0 : break
        #print "KK"
        RPKMlist[fl_idx].append(ii['RPKM'][0,0][0,0])
        dct_junc = {}
        dct_junc2 = {}

        read_matrix = ii['reads'][0,0]
        num_reads,num_ss = read_matrix.shape
        if num_reads == 0 or num_ss == 0: 
            id2+=1
            continue
      
        for jdx in range(num_reads):
            
            jname = ""
            cnt = 0
            for idx in range(num_ss):
                zz = read_matrix[jdx,idx]
#                print idx, zz
                if zz == 0: continue
                jname += ii['header'][0,0][0,idx][0]
                cnt +=1
                if jname == '': continue
                if cnt== 2: break
                if not jname in dct_junc:
                    dct_junc[jname] = 0
                    dct_junc2[jname] = [0]*77
                #print jname
                dct_junc[jname] += 1 
                
                if cnt == 1:
                    startpoint = ii['junc'][0,0][0,idx] - ii['coordinates'][0,0][jdx,0]
                    if startpoint >76 :
                        startpoint = 77
                    elif startpoint > 68 or startpoint <8:
                        continue
                    if zz >0:
                        junxPos[fl_idx][77-startpoint] += zz
                        dct_junc2[jname][77-startpoint] += zz
                    elif zz < 0:
                        junxPos_nu[fl_idx][77-startpoint] +=zz

        for jj in dct_junc.values():
            junxReads[fl_idx].append(jj)
            gname = "%s_%s"%(jdx,jname)
            if not gname in global_junxr:
                global_junxr[gname]=[0]*len(files)
            global_junxr[gname][fl_idx]=jj
        for indx in range(77):
            for jj in dct_junc2.values():
                if jj[indx] > 0:
                    xpos[indx].append(jj[indx])                     

        jj_aux = []
        for jj in dct_junc2.values():
            tp = 0
            for jidx,jji in enumerate(jj):
                if jji > 0:
                    junx_readatpos[fl_idx][jidx] +=1
                    tp +=1
                    jj_aux.append(jji)
            list_p.append( poisson.sf(jj_aux,np.mean(jj_aux)) )
            p2.append(np.mean(jj_aux))
            junxposhist[fl_idx].append(tp)
        id2+=1
#    if len(list_p) > 0:
#        pvalues[fl_idx] = np.concatenate((list_p))
#
#    for nn in pvalues[fl_idx]:
#        if float(nn) < float(10 **-3):
#            counter[fl_idx] +=1
#    percent = float(counter[fl_idx]*100)/len(pvalues[fl_idx])
#    print "%s Loosing : %d/%d (%.2f %%)"%(fl,counter[fl_idx], len(pvalues[fl_idx]), percent)


#import pickle
#fp = open('all_data.dat','wb')
#pickle.dump([RPKMlist,junxReads,junxposhist,junxPos,junx_readatpos,pvalues,counter,global_junxr],fp)


print "GENERATE HISTOGRAM"    

#''' RPKM x Genes plot'''
#pyplot.figure(fidx)
#for lb,ii in enumerate(RPKMlist):
#    pyplot.hist(ii, 250, range=(0,250),normed=1, histtype='step', cumulative=True, label=files[lb])
#pyplot.title("#Genes vs RPKM")                                                                                                                                                                                                                                                                                                                                   
#pyplot.xlabel("RPKM") 
#pyplot.legend(loc='lower right') 
#pyplot.ylabel("#Genes")
#pyplot.ylim((0,1.0))
#pyplot.grid()
#pyplot.show()
#fidx +=1
#
#''' #junctions vs #reads'''
#pyplot.figure(fidx)
#for lb,ii in enumerate(junxReads):
#    pyplot.hist(ii,200, range=(0,200),normed=1, histtype='step', cumulative=True,label=files[lb])
#pyplot.title("#Junctions vs #reads")                                                                                                                                                                                                                                                                                                                                   
#pyplot.xlabel("#reads") 
#pyplot.ylabel("#junctions") 
#pyplot.xlim((0,100))
#pyplot.ylim((0,1.0))
#pyplot.legend(loc='lower right') 
#pyplot.grid()
#pyplot.show()
#fidx +=1
#
#''' #junctions vs #pos for #reads # zero'''
#pyplot.figure(fidx)
#for lb,ii in enumerate(junxposhist):
#    pyplot.hist(ii,77, range=(0,77),normed=1, histtype='step', cumulative=True,label=files[lb])
#pyplot.title("junctions that have any read in that position")                                                                                                                                                                                                                                                                                                                                   
#pyplot.xlabel("position") 
#pyplot.grid()
#pyplot.ylabel("#junctions") 
#pyplot.xlim((0,68))
#pyplot.ylim((0,1.0))
#pyplot.legend(loc='lower right') 
#pyplot.show()
#fidx +=1
#
#
##''' subplots with additive count of reads'''
##ncol = 2
##nrow = max(2,len(tissues))
##f,splt = pyplot.subplots(nrow, ncol, sharex='col',sharey='row')
##ind = ['out'] + range(-76,0)
##for lb,ii in enumerate(junx_readatpos):
##    row = lb/ncol
##    col = lb%ncol
##    splt[row][col].bar(range(-77,0),ii)
##    splt[row][col].grid()                                                                                                                                                                                                                                                                                                               
##    splt[row][col].set_title(files[lb])
##    splt[row][col].set_xlabel("read start position")
##    splt[row][col].set_ylabel("#junctions") 
###pyplot.legend(loc='lower right') 
##pyplot.title("#Junctions vs start position")
##pyplot.show()
##fidx +=1
#
#
#

bins = [0]*len(files)
coef = np.zeros(shape=(len(files),60))
for lb, ii in enumerate(junxPos):
    jj = junx_readatpos[lb]
    for idx in range(8,68):
        if ii[idx] ==0 : coef[lb,idx-8]= 0
        else:
            coef[lb,idx-8]=float(ii[idx])/float(jj[idx])
#    print coef
    bins[lb] = mquantiles(coef[lb])[2]
    
mn = np.mean(bins) 
print mn

final= [0]*76
for ii in range(76):
    if ii<8 or ii >=68: continue
    pp = np.mean(coef[:,ii-8])
    if pp >= mn : 
        final[ii] = 1.0
    else:
        final[ii] = (pp) /mn

#''' subplots with additive count of reads'''
#
#for lb,ii in enumerate(junxPos):
#    nrow = 3
#    f,(a,b,c) = pyplot.subplots(nrow, 1, sharex=True)
#    a.bar(range(-77,0),ii)
#    a.grid()
#    a.set_title("A:Total number of reads per position")
#    a.set_ylabel("#reads") 
#    b.bar(range(-77,0),junx_readatpos[lb])
#    b.set_title("B:Number of junctions with reads in that position")
#    b.grid()
#    b.set_ylabel("#junctions") 
#    coef = [0.0]*77; idx=0
#    for x,y in zip(ii, junx_readatpos[lb]):
#        if x ==0 : coef[idx]= 0
#        else : coef[idx]=float(x)/float(y)
#        idx += 1
#    c.bar(range(-77,0),coef,color='g')
#    c.grid()
#    c.set_xticks(np.arange(-75,0,5))
#    c.set_xlabel("read start position")
#    c.set_ylabel("#reads/#junctions") 
#
##pyplot.title("#Junctions vs start position")
#pyplot.show()
#fidx +=1
#
#
##f,splt = pyplot.subplots(nrow, ncol, sharex='col',sharey='row')
##ind = ['out'] + range(-76,0)
##for lb,ii in enumerate(junxPos):
##    row = lb/ncol
##    col = lb%ncol
##    splt[row][col].bar(range(-77,0),junxPos_nu[lb],color='g')
##
##    splt[row][col].grid()                                                                                                                                                                                                                                                                                                               
##    splt[row][col].set_title(files[lb])
##    splt[row][col].set_xlabel("read start position")
##    splt[row][col].set_ylabel("#junctions") 
###pyplot.legend(loc='lower right') 
##pyplot.title("#Junctions vs start position")
##pyplot.show()
##fidx +=1
#
##print pvalues, len(pvalues)    
#a =arange(0,1,0.001)
#pyplot.figure(fidx)
#pyplot.plot(a,a,'k--',label="rand Poisson")
#for lb,ii in enumerate(pvalues):
#    pyplot.hist(ii,10000,range=(0,1), normed=1,histtype='step', cumulative=True,label=files[lb])
#pyplot.title("pvalues for junction positions") 
#pyplot.ylim((0,1.0))
#pyplot.grid()
#pyplot.legend(loc='lower right') 
#pyplot.show()
#fidx +=1
#
#pyplot.figure(fidx)
#pyplot.plot(a,a,'k--',label="rand Poisson")
#for lb,ii in enumerate(pvalues):
#    pyplot.hist(ii,10000,range=(0,1), normed=1,histtype='step', cumulative=True,label=files[lb])
#pyplot.title("pvalues for junction positions logscale") 
#pyplot.ylim((0,0.1))
#pyplot.xscale('log')
##pyplot.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
#pyplot.grid()
#pyplot.legend(loc='upper left') 
#pyplot.show()
#fidx +=1
#
#
#
#pyplot.figure(fidx)
#pyplot.plot(a,a,'k--',label="rand Poisson")
#for lb,ii in enumerate(pvalues):
#    pyplot.hist(ii,10000,range=(0,1), normed=1,histtype='step', cumulative=True,label=files[lb])
#pyplot.title("pvalues for junction positions logscale") 
#pyplot.ylim((0,0.1))
#pyplot.yticks(np.arange(0,0.1,0.005))
##pyplot.xscale('log')
##pyplot.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
#pyplot.grid()
#pyplot.legend(loc='upper left') 
#pyplot.show()
#fidx +=1
#
#
##f,splt = pyplot.subplots(len(tissues), 2, sharex='col')
#for kk,vv in tissues.items():
#    pyplot.figure(fidx)
#    mu = []
#    sigma = []
#    for junc,lst in global_junxr.items():
#        tmp = []
#        for idx in vv: 
#            if lst[idx]>0: tmp.append(lst[idx])
#        if len(tmp)>1:
#            mu.append(np.mean(tmp))
#            sigma.append(np.std(tmp))
#    if mu == []: 
#        print "There is no junc > 0"
#        continue
#    pyplot.scatter(mu,sigma)
#    pyplot.plot(a,a,'k--')
#    pyplot.grid()
#    pyplot.xlabel("mu")
#    pyplot.ylabel("sigma")  
#    pyplot.title("%s"%kk) 
#    pyplot.show()
#    fidx +=1
