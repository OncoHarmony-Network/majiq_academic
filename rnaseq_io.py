from gene import Gene, Transcript
from exon import Exon, detect_exons
from junction import  Junction
from rna_reads import RNA_read, parse_read
from RNAexperiment import RNAexperiment
import gc
import numpy as np



def reads_for_junc_coverage(filename, gene_list, readlen, exp_index):

    ''' Read and parse reads by junctions '''
    junctions = {}
    g_chromosome = ''
    counter = [0]*6
    total_reads = 0
    idx = {'+':0,'-':0}
    junctions = {'+':[],'-':[]}
    print "START READING ",filename
    sorted_gene_list = gene_list
    BUFFER = int(10E8) #1 gigabyte buffer
    file = open(filename, 'r')
    text = file.readlines(BUFFER)
    while text != []:
        print "New Buffer"
        for ii in text:
            if ii.startswith('@'): continue
            tab = ii.strip().split()
            chrom = tab[2]
            flag = int(tab[1])
            strand = '+' if flag & 0x10 == 0 else '-'
            (isRead,read) = parse_read( tab, readlen )
            if not isRead :
                r_start = int(tab[3])
                n_reads = read
            else:
                n_reads = read.get_read_count()
                r_start, r_end = read.get_coordinates()

            while idx[strand] < len (sorted_gene_list[strand]):
                gene = sorted_gene_list[strand][idx[strand]]
#                gene = sorted_gene_list[idx[strand]]
                (g_start, g_end) = gene.get_coordinates()
#                if isRead:
#                    print "###", r_start, g_start, read.get_junctions_info()
                if r_start < g_start :
                    ''' That means that we found new junctions out of any transcript'''
                    '''mark as non-transcript one and next junction'''
#                   print ii 
                    counter[1] += 1
                    break
                elif r_start > g_end or (strand != gene.strand):
                    '''next_gene'''
                    idx[strand] += 1
                    continue
                #end elif junc_start ...
                '''
                If the read is not a junction read, we add the count to the gene 
                and continue with the next read 
                '''
                total_reads += n_reads
                if not isRead : 
                    gene.add_read_count(n_reads,exp_index)
                    break
                gene.add_read(read,exp_index)
                counter[2] += 1
                
                j_list = gene.get_all_junctions()
                for (junc_start,junc_end) in read.get_junctions_info():
                    if junc_start - r_start > readlen : 
                        r_start = junc_start - (readlen - 16) -1
                    elif junc_start - r_start >= readlen -8 or junc_start -r_start <= 8: continue
                    found = False
                    if junc_end - junc_start < 10 :
                        counter[0] +=1
                        continue
                    for jj in j_list:
                        ( j_st,j_ed ) = jj.get_coordinates()
                        if   j_st > junc_start  or (j_st == junc_start and j_ed > junc_end): break
                        elif j_st < junc_start  or (j_st == junc_start and j_ed < junc_end): continue
                        elif junc_start == j_st and junc_end == j_ed:
                            ''' update junction and add to list'''
                            counter[3] +=1
                            jj.update_junction_read(exp_index,n_reads,r_start, read.get_GCcontent(), read.get_unique())
                            if not (junc_start,'5prime',jj) in junctions[strand]:
                                junctions[strand].append((junc_start,'5prime',jj))
                                junctions[strand].append((junc_end,'3prime',jj))
                            found = True
                            break
                        #end elif junc_start == ...
                    #end for jj in j_list
                    if not found :
                        ''' update junction and add to list'''
                        junc = None
                        for (coord,t,j) in junctions[gene.strand]:
                            if (j.start == junc_start and j.end == junc_end):
                                junc = j
                                junc.update_junction_read(exp_index, n_reads, r_start, read.get_GCcontent(), read.get_unique())
                                break
                            #end if (j.start) == ...
                        #end for (coord,t,j) ...
                        if junc is None:
                            '''mark a new junction '''
                            counter[4] += 1
                            junc = Junction( junc_start, junc_end, None, None, gene, readN=n_reads)
                            junc.update_junction_read(exp_index, n_reads, r_start, read.get_GCcontent(), read.get_unique())
                            junctions[strand].append((junc_start,'5prime',junc))
                            junctions[strand].append((junc_end,'3prime',junc))
                    #end if not found ...
                #end for junc ...
                break
            #end while ...
        #end for t in text
        gc.collect()
        text = file.readlines(BUFFER)
    #END WHILE
    file.close()

    print "INVALID JUNC",    counter[0]
    print "READ WRONG GENE", counter[1]
    print "READ IN GENE",    counter[2]
    print "READ FOUND JUNC", counter[3]
    print "READ NEW JUNC",   counter[4]
    print "READ ALL JUNC",   counter[5]

    for str, ll2 in junctions.items():
        junctions[str].sort()
    detect_exons(gene_list, junctions,None)
    skip_gene = 0
    non_skip = 0
    for strand, glist in gene_list.items():
        for gene in glist:
            gene.prepare_exons()
            re = gene.calculate_RPKM(exp_index, total_reads)
            if re == 0 :
                #print "SIN READS??? :",gene.id
                skip_gene +=1
            else:
                non_skip +=1
    print "Skipped genes without exons",skip_gene
    print " Non skipped",non_skip
    return 

def read_transcript_ucsc(filename, refSeq = False):
    print "Reading transcript ucsc format"
    BUFFER = int(10E6) #10 megabyte buffer
    file = open(filename, 'r')
    text = file.readlines(BUFFER)
    all_genes = {}
    while text != []:
        print "New Buffer"
        for t in text:
            t = t.strip()
            if t.startswith('#') or t == '':
                continue
            tab = t.split('\t')
            if refSeq :
                chrom = tab[2]
                start= int(tab[4])
                end = int(tab[5])
                strand = tab[3]
                gene_id = tab[12]
                transcript_id = tab[1]
                ex_start = tab[9].split(',')
                ex_end = tab[10].split(',')
                nblocks = int(tab[8])
            else:
                chrom = tab[0]
                ''' BED format file is 0 based, so the first base of the chromosome is 0. That means that the first real base of an exon is start+1'''
                start = int(tab[1])+1
                end = int(tab[2])
                strand = tab[5]
                gene_id = tab[3]
                transcript_id = tab[3]
                off_start = tab[11].split(',')
                off_len = tab[10].split(',')
                nblocks = int(tab[9])

                ex_start = [0]*nblocks
                ex_end = [0]*nblocks
                for ii in range(nblocks):
                    ex_start[ii] = int(off_start[ii])+start
                    ex_end[ii] = ex_start[ii] + int(off_len[ii])-1

#            print "New Transcript %s"%transcript_id
#            if not gene_id in all_genes :
            if not chrom in all_genes:
                all_genes[chrom] = {'+':[],'-':[]}
            gn = Gene(gene_id,chrom,strand,start,end)
            gene = gn.is_gene_in_list(all_genes[chrom][strand], gene_id)
            if not gene is None:
                del gn
                gn = gene
            else:
                all_genes[chrom][strand].append(gn)
            trcpt = Transcript(transcript_id, gn,start,end )
            gn.add_transcript(trcpt)
            pre_end = None
            pre_ex = None
            for ii in xrange(nblocks):
                start = ex_start[ii]
                end = ex_end[ii]
                ex = gn.exist_exon(start,end)
                if ex is None :
                    ex = Exon(start,end,gn,strand)
                    gn.add_exon(ex)
                txex = ex.add_new_definition(start, end, trcpt)
                trcpt.add_exon(txex)
                if not pre_end is None:
                    junc = gn.exist_junction(pre_end,start)
                    if junc is None:
                        junc = Junction(pre_end,start, pre_ex, ex,gn)
                    trcpt.add_junction(junc)
                pre_end = end
                pre_ex = ex
            #END for ii in range(nblocks)
            trcpt._sort_in_list(strand)
        #end for t in text
        gc.collect()
        text = file.readlines(BUFFER)
    #END WHILE
    file.close()
    print "SORTING GENES"
    n_genes = 0
    for chr in all_genes.keys():
        for strand, gg in all_genes[chr].items():
            n_genes += len(gg)
            all_genes[chr][strand] = sorted(gg)
            for gene in all_genes[chr][strand]:
                gene.prepare_exons()
    print "NUM_GENES",n_genes
    return all_genes



def read_test_file(filename):
    fp = open(filename,'r+')
    data = fp.readlines()
    fp.close()

    solution = []
    name = ''
#    comment = ''
    idx = 0
    for line in data:
        line = line.strip()
        if line.startswith('#'):
            line = line[1:]
            elm = line.split(':')
            shape = elm[1].split(",")
            name = elm[0]
            solution = []
            if elm[2] != '':
                solution = elm[2].split(",")
            idx = 0
            mat = np.ndarray(shape=(int(shape[0]),int(shape[1])),dtype='bool')
            continue
        elif line.startswith('%'):
#            comment = line[1:]
            continue
        elif not line:
            out =  analize.analize_gene_bin_matrix(mat)
            if len(solution) == 0 and len(out) == 0:
                sol = "OK"
            elif set(solution) == set(out):
                sol = "OK"
            else:
                sol = "FAIL, output,",out,"expected:",solution
            print "Test %s ..%s"%(name,sol)
            continue
        elm = line.split()
        for jj in range(len(elm)):
            mat[idx,jj] = (elm[jj] == '1')
        idx +=1
    return


def get_junctions_STAR( filename, gene_list_per_chrom ):
    input_f = open(filename,'r')
    junctions = {}

    sorted_gene_list = []
    gene_list = {}
    g_chromosome = ''

    wrong_gene = 0
    in_gene = 0
    new_junc = 0
    found_junc = 0
    all_junc = 0
    ppp = 0
    for line in input_f:
        line = line.strip()
        tab = line.split()
        chrom = tab[0]
        junc_start = int(tab[1]) - 1
        junc_end = int(tab[2]) #+ 1
        strn_type = int(tab[3])
        if junc_end - junc_start < 10 :
            ppp +=1
            continue
        if strn_type >0:
            str_list = ('+','-')
        else:
            str_list = '+' if strn_type == 1 else '-'
        all_junc +=1
        if chrom != g_chromosome :
            if not chrom in gene_list_per_chrom: continue
            gene_list = gene_list_per_chrom[chrom]
            g_chromosome = chrom
            sorted_gene_list = sorted(gene_list.values())
            '''first_gene'''
            idx = {'+':0,'-':0}
        for strand in str_list:
            found = False
            while idx[strand] < len (sorted_gene_list):
                gene = sorted_gene_list[idx[strand]]
                (g_start, g_end) = gene.get_coordinates()
                if junc_end < g_start :
                    ''' That means that we found new junctions out of any transcript'''
                    '''mark as non-transcript one and next junction'''
                    wrong_gene += 1
#                    junc = Junction( ss5p, ss3p, None, None, None, readN=tab[6] )
#                    junctions[chrom][strand].append(junc)
                    break
                elif junc_start > g_end or (strand != gene.strand):
                    '''next_gene'''
                    idx[strand]  += 1
                    continue
                j_list = gene.get_all_junctions()
                if not chrom in junctions :
                    junctions[chrom] = {'+':[],'-':[]}
                in_gene +=1
                for jj in j_list:
                    ( j_st,j_ed ) = jj.get_coordinates()
                    if j_st > junc_start  or (j_st == junc_start and j_ed > junc_end):
                        break
                    elif j_st < junc_start or (j_st == junc_start and j_ed < junc_end):
                        continue
                    elif junc_start == j_st and junc_end == j_ed:
                        found_junc +=1
                        ''' update junction and add to list'''
                        jj.add_read_number(int(tab[6]),exp_index)
                        if not (junc_start,'5prime',jj) in junctions[chrom][strand]:
                            junctions[chrom][strand].append((junc_start,'5prime',jj))
                            junctions[chrom][strand].append((junc_end,'3prime',jj))
                        found = True
                        break
                if not found :
                    junc = None
                    for (coord,t,j) in junctions[chrom][gene.strand]:
                        if (j.start == junc_start and j.end == junc_end):
                            junc = j
                            junc.add_read_number(int(tab[6]),exp_index)
                            break
                    if junc is None:
                        new_junc +=1
                        '''mark a new junction '''
                        junc = Junction( junc_start, junc_end, None, None, gene, readN=int(tab[6]))
                        junctions[chrom][strand].append((junc_start,'5prime',junc))
                        junctions[chrom][strand].append((junc_end,'3prime',junc))
                break

    print "JUNC WRONG GENE", wrong_gene
    print "JUNC IN GENE", in_gene
    print "JUNC NEW JUNC", new_junc
    print "JUNC FOUND JUNC", found_junc
    print "JUNC ALL JUNC", all_junc
    print "INVALID JUNC", ppp
    for chr, ll in junctions.items():
        for str, ll2 in ll.items():
            junctions[chr][str].sort()
            print ll2
#            print chr, str, len(ll2)
    detect_exons(gene_list_per_chrom, junctions,None)


def read_triplets_bed(filename,all_genes):
    dict = {}
    info = {}
    inp = open(filename)
    for line in inp:
        line = line.strip()
        tab = line.split()
        id = tab[3].split('_')
        reg = id[-1]
        name = '_'.join(id[:-1])
        if not tab[0] in info:
            info [tab[0]] = {}
            info [tab[0]]['+'] = []
            info [tab[0]]['-'] = []
        
        if not name in dict:
            info[tab[0]][tab[5]].append(name)
            dict[name]={}
        dict[name][reg]=(int(tab[1]), int(tab[2]))

    inp.close()
    notfound = []
    ex_found = {}
    pepepe = []
    found_in_gene = 0
    incomplete_events = 0
    found_ex = 0
    exons_w_var = {'C1':0,'A':0,'C2':0}

    for chr,it in info.items():
        if not chr in all_genes :
            for strand, list in it.items():
                for name in list:
                    notfound.append(name)
                    notfound.append(name)
                    notfound.append(name)
            continue
        for strand, list in it.items():
            if strand == '+':
                start_reg = 'C1'
                end_reg = 'C2'
            else:
                start_reg = 'C2'
                end_reg = 'C1'
            gene_list =(all_genes[chr][strand])
            gene_list.sort()
            jj = {'+':0,'-':0}
            for name in list:
                vv = dict[name]
                #jj =0
                while jj[strand] < len(gene_list):
                    gg = gene_list[jj[strand]]
                    (gg_start,gg_end) = gg.get_coordinates()
                    if gg_end< vv[start_reg][0] or gg.strand != strand:
                        jj[strand]+=1
                        continue
                    elif gg_start> vv[end_reg][1]:
#                        jj+=1
#                        pass
                        
                        notfound.append(vv)
                        notfound.append(vv)
                        notfound.append(vv)
                        break
                    else:
                        found_in_gene += 1
                        ex_list = gg.get_exon_list()
                        triplet = []
                        idx = 0
                        for reg in (start_reg,'A',end_reg):
                            while idx< len(ex_list):
                    #check exon
                                (start,end) = ex_list[idx].get_coordinates()
                                #print name, gg.id, start,end, vv[reg]
                                if end < vv[reg][0] :
                                    #print "NEXT"
                                    idx +=1
                                    continue
                                elif start > vv[reg][1]:
                                    #print "WRONG"
                                    notfound.append(vv[reg])
#                                    idx +=1
                                    break
                                else:
                                    #print "FOUND",ex_list[idx].get_coordinates()
                                    found_ex +=1
                                    if len(ex_list[idx].exonTx_list ) >1 :
                                        exons_w_var[reg] += 1
                                     #   if reg != 'A':
                                        triplet.append(ex_list[idx]) 
                                    else:
                                        triplet.append(ex_list[idx]) 
                                    break
                        if len(triplet) == 3:
                            #print triplet[0].id, triplet[0].get_coordinates()
                            #print triplet[1].id, triplet[1].get_coordinates()
                            #print triplet[2].id, triplet[2].get_coordinates()
                            if gg not in ex_found:
                                ex_found[gg] = []
                            ex_found[gg].append(triplet)
                            pepepe.append((name, "%s.%s"%(gg.id,triplet[1].id)))
                        else:
                            incomplete_events +=1
                        break
                if jj == len(gene_list):
                    notfound.append(vv)
                    notfound.append(vv)
                    notfound.append(vv)

    total_ex = 0
    total_ex += sum([len(xx) for xx in ex_found.values()])

    print "Number of found exons :%d"%found_ex
    print "\twith splice variants C1:%d, A:%d, C2:%d"%(exons_w_var['C1'],exons_w_var['A'],exons_w_var['C2'])
    print "Number of non-found exons :%d"%len(notfound)
    print "---\nNumber of complete triplets: %d"%len(ex_found), total_ex
    print "Number of incomplete triplets: %d"%incomplete_events
    print " found in gene %d"%found_in_gene

#    fp = open ('./test','w')
#    for ii in pepepe:
#        fp.write("%s %s\n"%ii)

    return  ex_found



