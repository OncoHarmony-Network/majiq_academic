import numpy as np
import math
from gene import Gene, Transcript
from exon import Exon
from junction import Junction
from utils import utils, stats
import globals

def analize_genes(gene_list, out_file, pos_list, neg_list, ASvsConst = 'AS'):
 #   of = open("%s_%s"%(out_file,ASvsConst), 'w+')
    ok_detect = 0
    wrong_detect = 0
    total = 0
    out_list=set()

    for strand,glist in gene_list.items():
        for gn in glist:
            as_db_list = []
            gn.get_exon_list()
            if len(gn.transcript_list ) == 1 and ASvsConst == 'AS': continue
            if len(gn.exons) < 3: continue
            mat = np.ndarray(shape=(len(gn.transcript_list),len(gn.exons)),dtype='bool')
            idx_t = 0
            for tpt in gn.transcript_list:
                for g_ex in gn.exons:
                    if set(tpt.exon_list).intersection( set(g_ex.exonTx_list)) :
                        mat[idx_t,g_ex.id-1]= 1
                    else:
                        mat[idx_t,g_ex.id-1]= 0

                idx_t +=1
                
            if ASvsConst == 'AS':
                out = analize_gene_bin_matrix(mat)
            else:
                out = analize_bin_matrix_const_firstlast(mat)
            
            if gn in pos_list:
                for (ex1,ex2,ex3) in pos_list[gn]:
                    as_db_list.append((ex1.get_id(),ex2.get_id(), ex3.get_id()))
            merge_as = set(out).union(set(as_db_list))
            
            gn.add_transcript_AS_candidates(list(merge_as))
#            gn.add_transcript_AS_candidates(out)
            total += len(out)
##            for extrip in out:
#
#                exon2 =  gn.exons[extrip[2]]
#                if len(exon2.exonTx_list) > 1 :
#                    continue
#                out_list.add(exon2)
#                
##                tsc = gn.transcript_list[extrip[0]]
#                exon1 =  gn.exons[extrip[1]]
#                if len(exon1.exonTx_list) > 1 :
#                    continue
##                txex = exon1.get_exon_definition(tsc)
##                ex1 = ("%s-%s"%txex.get_coordinates())
#
##                txex = exon2.get_exon_definition(tsc)
##                ex2 = ("%s-%s"%txex.get_coordinates())
#                
#                exon3 =  gn.exons[extrip[3]]
##                txex = exon3.get_exon_definition(tsc)
##                ex3 = ("%s-%s"%txex.get_coordinates())
#                if len(exon3.exonTx_list) > 1 :
#                    continue
#
#               # out_list.append(exon2)
#               # "%s:%s"%(name,extrip[2]))
#
#                found = False
#                for pp in pos_list:
#                    if exon1 == pp[0] and exon2 == pp[1] and exon3 == pp[2]:
#                        ok_detect +=1
#                        found = True
#                         del pp
#                        break
#                if found : continue
#                for pp in neg_list:
#                    if exon1 == pp[0] and exon2 == pp[1] and exon3 == pp[2]:
#                        wrong_detect +=1
#                        found = True
#                        del pp
#                        break
##                if not found:
##                    of.write("%s:%s:%s:%s:%s\n"%(chr,gn.strand,ex1,ex2,ex3 ))
##    of.close()
    print "%s \tCORRECT: %d, FALSE: %d, TOTAL: %d"%(ASvsConst,ok_detect, wrong_detect,total)
    return out_list

def analize_bin_matrix_const_firstlast(mat):
    
    crrct = []
    incrrct = []
    trans = []
    out = []
    for jj in range(1,mat.shape[1]-1):
        bords = (jj-1,jj,jj+1)
        for ii in range(mat.shape[0]):
            if mat[ii, jj-1] != mat[ii,jj] or mat[ii,jj+1] != mat[ii,jj]: 
                if not bords in incrrct:
                    incrrct.append(bords)
                break
#            i_pre = jj -2
#            while i_pre >= 0 :
#                if mat[ii,i_pre] == 1: break
#                i_pre -=1
#
#            i_post = jj +2
#            while i_post < mat.shape[1] :
#                if mat[ii,i_post] == 1: break
#                i_post +=1

#            if i_pre >= 0 and i_post < mat.shape[1]: 
#                continue

            if mat[ii,jj] == 1:
                if not bords in crrct:
                    crrct.append(bords)
                    trans.append(ii)
    for idx,ii in enumerate(crrct) :
        if not ii in incrrct:
            #out.append((trans[idx],ii[0],ii[1],ii[2]))
            out.append((ii[0],ii[1],ii[2]))
    return out

def analize_gene_bin_matrix(mat):
    out = []
    v_101 = []
    v_111 = []
    trans = []
    for jj in range(1,mat.shape[1]-1):
        for ii in range(mat.shape[0]):
            i_pre = jj -1
            while i_pre >= 0 :
                if mat[ii,i_pre] == 1: break
                i_pre -=1
            if i_pre< 0 : continue

            i_post = jj +1
            while i_post < mat.shape[1] :
                if mat[ii,i_post] == 1: break
                i_post +=1
            if i_post >= mat.shape[1]: continue

            bords = (i_pre,jj,i_post)
            if mat[ii,jj] == 1:
                if not bords in v_111:
                    v_111.append(bords)
                    trans.append(ii)
            else:
                if not bords in v_101:
                    v_101.append(bords)
    for idx,ii in enumerate(v_111) :
        if ii in v_101:
            #out.append((trans[idx],ii[0],ii[1],ii[2]))
            out.append(ii)
    return out


def analize_junction_reads_with_DB(gene_list, alt_list, const_list, b_list ): 

    ok_detect = [0,0]
    wrong_detect = [0,0]
    total = [0,0]
    ASvsConst = ['ALT','CONST']
    list = [alt_list,const_list]
#    out_files = ("./altern.ok.events","./const.ok.events")
    #fp_list = [0,0]
    conj_list = [0,0]
    for idxn in (0,1):
        conj_list[idxn] = set()
#        fp_list[idxn] = open(out_files[idxn],'w+')
#        fp_list[idxn].write("#EVENTS from AVISPA DB detected on junctions detection methods\n")

    for chr in gene_list.keys():
        for name,gn in gene_list[chr].items():

            ex_list = gn.get_exon_list()
            mat = np.zeros(shape=(len(ex_list),len(ex_list)),dtype='int')
            junc_list = gn.get_all_junctions()
            for junc in junc_list:
                if junc.donor is None or junc.acceptor is None or junc.donor.id is None or junc.acceptor.id is None: continue
                if junc.donor.id >= len(ex_list) or junc.acceptor.id >= len(ex_list): continue
                #if junc.donor.id > junc.acceptor.id : 
#                    print "KAKAKAKAKAKAKAKAKAK",junc.donor,junc.acceptor,junc.donor.id,junc.acceptor.id, junc.donor.get_coordinates(), junc.acceptor.get_coordinates()
                mat[junc.donor.id-1, junc.acceptor.id-1] = junc.readN

            (alt,const) = rnaSeq_const_detection(mat, b_list)
            idxn = 0
            prmat = False
            for out in (alt,const):
                total[idxn] += len(out)
                pos_list = list[idxn]
                neg_list = list[1-idxn]
                for ii in out:
                    exon1 =  gn.exons[ii-1]
                    if len(exon1.exonTx_list) > 1 :
                        continue
                    exon2 =  gn.exons[ii]
                    if len(exon2.exonTx_list) > 1 :
                        continue
                    exon3 =  gn.exons[ii+1]
                    if len(exon3.exonTx_list) > 1 :
                        continue
                    found = False
                    conj_list[idxn].add((exon1,exon2,exon3))
                    #conj_list[idxn].append(exon2)
#"%s:%s"%(gn.id, ii))
                    for pp in pos_list:
                        if exon1 == pp[0] and exon2 == pp[1] and exon3 == pp[2]:
#                            fp_list[idxn].write("%s.%s.%s.%s\n"%(gn.id,pp[0].id, pp[1].id, pp[2].id))
                            ok_detect[idxn] +=1
                            found = True
                            break
                    if found : continue
                    for pp in neg_list:
                        if exon1 == pp[0] and exon2 == pp[1] and exon3 == pp[2]:
                            if not prmat :
                                prmat = True
#                                print ASvsConst[idxn],pp[0].id, pp[1].id, pp[2].id
#                                print mat
                            wrong_detect[idxn] +=1
                            found = True
                            del pp
                            break
                idxn +=1

    if not False in b_list:
        print "RESULT ALL ENABLED"

    elif not b_list[0]:
        print "RESULT I"

    elif not b_list[1]:
        print "RESULT II"

    elif not b_list[2]:
        print "RESULT III"

    for idxn in range(2):
        print "%s \tCORRECT: %d, FALSE: %d, TOTAL: %d"%(ASvsConst[idxn],ok_detect[idxn], wrong_detect[idxn],total[idxn])
    
    return conj_list




def __junction_filter_check( junc): 

    ''' Check total_read >= 10 '''
    filter = False
#    if junc.readN[exp_idx] < 10 : return False 
    count = 0
    for exp_idx in range(globals.num_experiments) :
        cov_cnt = 0
        if junc.readN[exp_idx] >= 10: 
            for val in junc.coverage[exp_idx]:
                if val >0 : cov_cnt += 1
            if cov_cnt < 3 : continue
            count +=1
        if count > (0.1 * (globals.num_experiments)) :
            filter = True
            break

    return (filter)


def __get_enabled_junction(con):
    print "KK"
    for jrow in con:
        for jj in jrow:
            print jj
            if jj is None: continue
            break
        else: 
            continue
        break
    return jj

def analize_junction_reads( gene_list,chr ):


    num_discard = 0
    junc_set = [[]]*globals.num_experiments

    tab_out = [0]*globals.num_experiments
    for exp_idx in range(globals.num_experiments) :
        name = globals.exp_list[exp_idx]
        tab_out[exp_idx] = open("./%s.corr.tab"%name,"w+")

    total = 0
    total_cisfrm = 0
    total_aisfrm = 0
    some_none = 0
    overlp = 0
    notjunc = 0
    ss_variant = 0
    for strand, glist  in gene_list.items():
        for gn in glist:
            count = 0
            for ll in gn.get_RNAread_list():
               count += len(ll)
            if count == 0: continue
#            print "gene",name
            ex_list = gn.get_exon_list()
            ss3_l = []
            ss5_l = []
            tlb = {}
            exidx = 0
            for ex in ex_list:
                if ex.id is None: continue
                s3 = [ss3 for ss3 in set(ex.ss_3p_list)]
                s5 = [ss5 for ss5 in set(ex.ss_5p_list)]
                s3 = sorted(list(set(s3)))
                s5 = sorted(list(set(s5)))
                st3 = len(ss3_l)
                st5 = len(ss5_l)
                ss3_l += s3
                ss5_l += s5
                tlb[exidx] = [range(st3,len(ss3_l)),range(st5,len(ss5_l))]
                exidx += 1
#            print tlb
            mat  = np.zeros(shape=(len(ss5_l),len(ss3_l)),dtype='int')
            jmat = np.zeros(shape=(len(ss5_l),len(ss3_l)),dtype='object')
            jmat.fill(None)
            junc_list = gn.get_all_junctions()
#            print junc_list
            for junc in junc_list:
#                print "PREE" , junc.get_coordinates(), junc.readN
                st,end = junc.get_coordinates()
                if not st in ss5_l or not end in ss3_l:
                    notjunc += 1
                    continue
                overlp +=1
                x = ss5_l.index(st)
                y = ss3_l.index(end)
                mat [ x, y ] = junc.readN.sum()
                jmat[ x, y ] = junc

            print "KAKAKAKAKA",len(gn.get_transcript_AS_candidates())
            print gn.get_transcript_AS_candidates()
            (alt, cisfrm,aisfrm) = rnaSeq_const_detection(mat, tlb, (True,True,True),gn.get_transcript_AS_candidates())

#            if len(alt) > 0:
#                print "VALUES : ",gn.id, "CANDIDATES:",alt
#                utils.print_junc_matrices(mat,tlb,1)
#                print "END VALUES"
            total += len(alt) 
            total_cisfrm += len(cisfrm)
            total_aisfrm += len(aisfrm)
            for ii in alt:
                a  = tlb[ii]
                c1 = tlb[ii-1]
                c2 = tlb[ii+1]

                c1_a = jmat[ c1[1][0] : c1[1][-1]+1,  a[0][0] :  a[0][-1]+1 ]
                a_c2 = jmat[  a[1][0] :  a[1][-1]+1, c2[0][0] : c2[0][-1]+1 ]
                c1c2 = jmat[ c1[1][0] : c1[1][-1]+1, c2[0][0] : c2[0][-1]+1 ]
               
#                if c1a is None  or ac2 is None or c1c2 is None: 
#                    some_none +=1
#                    continue
#                print jj,c1a,ac2,c1c2
    
                pass_filter = (True)
                variant = False
                for junc_idx, junc in enumerate((c1_a,c1c2)):
                    pre = 0
                    junc_or = False
                    for junc_row in junc:
                        for jj in junc_row:
#                            print "KK",jj
                            if jj is None: continue
                            pre +=1
                            print "START", junc_idx
                            print jj.readN, jj.coverage
                            print "END"
                            junc_or = junc_or or __junction_filter_check(jj)
                    pass_filter = pass_filter or junc_or
                    if pre >2: 
                        variant = True 
                if variant:
                    ss_variant +=1
                    continue
                if not pass_filter: 
                    num_discard +=1
                    continue
                    
                elif pre < 2:
                    some_none +=1
#                    continue
            
                jc1a  = __get_enabled_junction(c1_a)
                jc1c2 = __get_enabled_junction(c1c2)
                for exp_idx in range(globals.num_experiments):
                    junc_set[exp_idx].append((jc1a,jc1c2))
                    
                    ''' Correlation check: Createing tab files for experiment '''
                    
                    x,y = jc1c2.coverage.shape
                    if jc1a is None:
                        n1 = 0
                    else:
                        n1 =  jc1a.coverage[exp_idx][y-38]
                    n2 = jc1c2.coverage[exp_idx][y-38]
                    n3 = jc1c2.coverage[exp_idx][y-60]

#                    print "####", jc1a.coverage[exp_idx],jc1c2.coverage[exp_idx]
                    print "#######", n1, n2 ,n3
#                    exit()
                    tab_out[exp_idx].write("%s.%s\t%s\t%s\t%s\n"%(name,ii+1,n1,n2,n3))


                #STEP 3: Delete stacked

                #STEP 4:  multinomial distribution
                stats.junction_expression(jc1c2.coverage[exp_idx],jc1c2.get_gc_content()[exp_idx],globals.gc_factor[exp_idx],0.8)

#print len(junc_set[0])
    [x.close() for x in tab_out]
    print "AS %s DISCARDED JUNCTIONS PER experiment"%chr,num_discard,"/",total, some_none
    print "AS %s constitutive isoform"%total_cisfrm
    print "AS %s AS isoform present in transcript analysis"%total_aisfrm
    print "AS %s SKIPPEDO junction"%chr, notjunc,"/",overlp
    print "AS %s How many events with ss variants"%chr, ss_variant
    return

def rnaSeq_const_detection(mat, exon_to_ss, b_list, pre_list=None):
    '''
       Rules for const are:
        1. All the junction from A should go to C1 or C2
        2. All the junction from C1 should go to A
        3. All the junction to C2 should come from A
        4. Nomber of reads from C1-A should be equivalent to number of reads from A-C2
    '''
    alt_list = []
    alt_Cisfrm = []
    alt_Aisfrm = []
    const = abs(math.log(1.5/1.0,2))

    #change bucle for iterate by exons
    for ii in range( 1, len(exon_to_ss) -1 ) :

        a = exon_to_ss[ii]
        c1 = exon_to_ss[ii-1]
        c2 = exon_to_ss[ii+1]

#        print c1[1][-1]
        c1_a = mat[ c1[1][0] : c1[1][-1]+1,  a[0][0] :  a[0][-1]+1 ]
        a_c2 = mat[  a[1][0] :  a[1][-1]+1, c2[0][0] : c2[0][-1]+1 ]
        c1c2 = mat[ c1[1][0] : c1[1][-1]+1, c2[0][0] : c2[0][-1]+1 ]
       
        toA    = mat[          : c1[1][0],     a[0][0]    :  a[0][-1]+1 ]
        toC2   = mat[          : c1[1][0],    c2[0][0]    : c2[0][-1]+1 ]
        fromA  = mat[  a[1][0] :  a[1][-1]+1, c2[0][-1]+1 :             ]
        fromC1 = mat[ c1[1][0] : c1[1][-1]+1, c2[0][-1]+1 :             ]
        
#        print "C1_A"
#        utils.print_junc_matrices(c1_a,1)
#        print "A_C2"
#        utils.print_junc_matrices(a_c2,1)
#        print "C1_C2"
#        utils.print_junc_matrices(c1c2,1)
#        print "END CCC"

        pre_read = c1_a.sum()
        post_read = a_c2.sum()
        x = abs(math.log((float(pre_read)+0.1)/(float(post_read)+0.1),2))
        if b_list[2] and x>const : continue
        if b_list[1] and np.count_nonzero(toA) >0 : continue #TO A
        if b_list[1] and ii+2 < mat.shape[1] and np.count_nonzero(fromA)  >0 : continue #FROM A
        if b_list[0] and ii+2 < mat.shape[1] and np.count_nonzero(fromC1) >0 : continue #FROM C1
        if b_list[0] and np.count_nonzero(toC2) >0 : continue # TO C2 

        if c1c2.sum() > 0 :
            #print "ALT"
            alt_list.append(ii)
            if (not pre_list is None) and (ii-1,ii,ii+1) in pre_list:
                alt_Aisfrm.append((ii-1,ii,ii+1))
        else:
            if (not pre_list is None) and (ii-1,ii,ii+1) in pre_list:
                print "LLL"
                alt_Cisfrm.append((ii-1,ii,ii+1))


    return (alt_list,alt_Cisfrm,alt_Aisfrm)

