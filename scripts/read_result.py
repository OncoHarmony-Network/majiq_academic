import scipy.io
from collections import defaultdict
import numpy as np
from pylab import *
from matplotlib import pyplot


file = '/Volumes/data/ucsc/reads/test_1k/test.mat'
heart1 = scipy.io.loadmat(file)


RPKMlist = []
junxReads = []
junxPos = []
id2 = 0
for ii in [heart1['RNASEQ_INFO'][15]]:
    
    #print "GET RPKM ", id2
    if ii[0][0][0] == 0 : break
    #print "KK"
    RPKMlist.append(ii[0]['RPKM'][0,0][0,0])
    dct_junc = {}
    print "pre reads"
    for jj in ii[0]['reads']:
        stop()
        jname = ""
        cnt = 0
        if len(jj[0])==0: break
        print jj[0][0]
        for idx,zz in enumerate(jj[0][0]):
            if zz == 0: continue
            jname += ii[0]['header'][0,0][idx,0][0]
            cnt +=1
            if cnt== 2: break
        if jname == '': continue
        if not jname in dct_junc:
            dct_junc[jname] = 0             
        print jname
        dct_junc[jname] += 1 
    
    print len(dct_junc)
    for jj in dct_junc.values():
        print "KKKKK"
        junxReads.append(jj)
    id2+=1

print "GENERATE HISTOGRAM"

pyplot.figure(1)
pyplot.hist(RPKMlist, 1000, range=(0,2000))
#pyplot.hist(RPKMlist,range=(0,2000), bins=1000)
pyplot.show()
pyplot.figure(2)
pyplot.hist(junxReads)
pyplot.show()
