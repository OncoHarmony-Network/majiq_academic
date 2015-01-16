from gene import Gene, Transcript
from exon import detect_exons
from grimoire.junction import Junction
import grimoire.exon as majiq_exons
import pysam
import gc
import os
import mglobals
from collections import namedtuple
import gzip
import urllib
import cPickle as pickle


def create_if_not_exists(my_dir, logger=False):
    """Create a directory path if it does not exist"""
    try:
        if logger:
            logger.info("\nCreating directory %s..." % my_dir)
        os.makedirs(my_dir)
    except OSError:
        if logger:
            logger.info("\nDirectory %s already exists..." % my_dir)


def load_bin_file(filename, logger=None):

    if not os.path.exists(filename):
        if logger:
            logger.error('Path %s for loading does not exist' % filename)
        return

    fop = open(filename, 'rb')

    fast_pickler = pickle.Unpickler(fop)
#    fast_pickler.fast = 1
    data = fast_pickler.load()
    fop.close()
    return data


def dump_bin_file(data, filename):

    with open(filename, 'wb') as ofp:
       fast_pickler = pickle.Pickler(ofp, protocol=2)
       #fast_pickler.fast = 1
       fast_pickler.dump(data)
#         pickle.dump(data, protocol=2)


def __cross_junctions(read):

    """
      This part will parse the jI from STAR
      # Column 17: jI:B:I,Start1,End1,Start2,End2,... Start and End of introns for all junctions (1- based)
      # jI:B:I,-1 means that the read doesn't cross any junction
    """

    jlist = []
    cross = False
    try:
        list_junc = read.opt('jI')
        if len(list_junc) > 1 and list_junc[0] != -1:
            for idx in range(0, len(list_junc), 2):
                junc_start = int(list_junc[idx])-1
                junc_end = int(list_junc[idx+1])+1
                jlist.append((junc_start, junc_end))
            #end for idx ...
            cross = True
        #end else ...
    except KeyError:
#    if len(jlist) != 0: print "STAR ::",jlist
#    print"THIS IS NOT a WELL defined STAR output"
        off = 0
        for op, num in read.cigar:
            if op in [0, 5, 6, 7, 8]:
                off += num
            elif op in [1, 5]:
                off += 0
            elif op == 2:
                off += num
            elif op == 3:
                jlist.append((read.pos+off, read.pos+off+num+1))
                off += num
#    if len(jlist) !=0 : print "NOSTAR:", jlist, read.cigar

    return cross, jlist


def __is_unique(read):
    unique = True
    try:
        if read.opt('NH') > 1:
            unique = False
    except KeyError:
        if read.flag & 0x100 == 1:
            unique = False
    return unique


def __get_num_reads(read):
    try:
        nreads = int(read.opt('HI'))
    except KeyError:
        nreads = 1
#    return 1
    return nreads


def count_mapped_reads(filename, exp_idx):
    stats = pysam.flagstat(filename)
    mapped_reads = int(stats[2].split()[0])
    mglobals.num_mapped_reads[exp_idx] = mapped_reads


def is_neg_strand(read):


    res = False
    if read.flag & 0x10 == 0x10:
#        print "FLAG",read.flag
        res = True

    if mglobals.strand_specific:
        res = not res

    return res


def get_junc_from_list(coords, list_elem):
    res = None
    for xx in list_elem:
        newcoord = xx.get_coordinates()
        if newcoord[0] == coords[0] and newcoord[1] == coords[1]:
            res = xx
            break
    return res


def rnaseq_intron_retention(filenames, gene_list, readlen, chrom, logging=None):

    samfile = [pysam.Samfile(xx, "rb") for xx in filenames]
    intron_ret_list = []
    for strand in ('+', '-'):
        for gne in gene_list[strand]:
            intron_list = gne.get_all_introns()
            virtua_juncs = []
            for exp_index in range(len(filenames)):
                for exon1, exon2 in intron_list:
                    ex1_start, ex1_end = exon1.get_coordinates()
                    ex2_start, ex2_end = exon2.get_coordinates()
                    v_junc1 = [ex1_end, ex1_end+1]
                    v_junc2 = [ex2_start-1, ex2_start]
                    new_junc = True
                    offset = readlen - 8
                    try:
                        read_iter = samfile[exp_index].fetch(chrom, ex1_end + 1 + 8, ex1_end + 1 + offset)
                    except ValueError:
                        #logging.info('There are no reads in %s:%d-%d' % (chrom, ex1_end, ex1_end+1))
                        continue

                    junc1 = get_junc_from_list(v_junc1, virtua_juncs)
                    if junc1 is None:
                        junc1 = Junction(v_junc1[0], v_junc1[1], exon1, None, gne, readN=0)
                        virtua_juncs.append(junc1)
                    else:
                        new_junc &= False

                    for read in read_iter:
                        is_cross, junc_list = __cross_junctions(read)
                        r_start = read.pos
                        if is_cross or r_start > ex1_end - 8:
                            continue

                        strand_read = '+' if not is_neg_strand(read) else '-'
                        unique = __is_unique(read)
                        if strand_read != strand or not unique:
                            continue

                        nreads = __get_num_reads(read)
                        gne.add_read_count(nreads, exp_index)


                        #intron_ju
                        nc = read.seq.count('C') + read.seq.count('c')
                        ng = read.seq.count('g') + read.seq.count('G')
                        gc_content = float(nc + ng) / float(len(read.seq))

                        junc1.update_junction_read(exp_index, nreads, r_start, gc_content, unique)

                    try:
                        read_iter = samfile[exp_index].fetch(chrom, ex2_start - offset - 1, ex2_start - 1 - 8)
                    except ValueError:
                        #logging.info('There are no reads in %s:%d-%d' % (chrom, ex1_end, ex1_end+1))
                        continue

                    junc2 = get_junc_from_list(v_junc2, virtua_juncs)
                    if junc2 is None:
                        junc2 = Junction(v_junc2[0], v_junc2[1], None, exon2, gne, readN=0)
                        virtua_juncs.append(junc2)
                    else:
                        new_junc &= False

                    for read in read_iter:
                        is_cross, junc_list = __cross_junctions(read)
                        r_start = read.pos
                        if is_cross or r_start < (ex2_start - 1 - offset):
                            continue
                        strand_read = '+' if not is_neg_strand(read) else '-'
                        unique = __is_unique(read)
                        if strand_read != strand or not unique:
                            continue

                        nreads = __get_num_reads(read)
                        gne.add_read_count(nreads, exp_index)
                        r_start = read.pos

                        #intron_ju
                        nc = read.seq.count('C') + read.seq.count('c')
                        ng = read.seq.count('g') + read.seq.count('G')
                        gc_content = float(nc + ng) / float(len(read.seq))
                        junc2.update_junction_read(exp_index, nreads, r_start, gc_content, unique)

                    if new_junc:
                        majiq_exons.new_exon_definition(v_junc1[1], v_junc2[0], None, junc1, junc2, gne, isintron=True)
                        logging.info("NEW INTRON RETENTION EVENT %s, %d-%d" % (gne.get_name(), v_junc1[0], v_junc2[1]))


def read_sam_or_bam(filenames, gene_list, readlen, chrom, nondenovo=False, logging=None):

    counter = [0] * 6
    samfile = [pysam.Samfile(xx, "rb") for xx in filenames]
    temp_ex = []
    non_unique_num = 0
    skip_gene = 0
    non_skip = 0

    for strand in ('+', '-'):
        for gne in gene_list[strand]:
            junctions = []
            strt, end = gne.get_coordinates()
            j_list = gne.get_all_junctions()
            ex_list = gne.get_exon_list()

            for exp_index in range(len(filenames)):
                try:
                    read_iter = samfile[exp_index].fetch(chrom, strt, end)
                except ValueError:
                    logging.info('There are no reads in %s:%d-%d' % (chrom, strt, end))
                    continue
                for read in read_iter:
                    strand_read = '+' if not is_neg_strand(read) else '-'

                    if strand_read != strand:
                        continue
                    unique = __is_unique(read)
                    if not unique:
                        non_unique_num += 1
                        continue
                    nreads = __get_num_reads(read)
                    gne.add_read_count(nreads, exp_index)
                    is_cross, junc_list = __cross_junctions(read)
                    r_start = read.pos
                    if r_start < strt or r_start > end:
                        continue

                    for ex_idx in range(len(ex_list)):
                        ex_start, ex_end = ex_list[ex_idx].get_coordinates()
                        if ex_start <= r_start <= ex_end:
                            ex_list[ex_idx].update_coverage(exp_index, nreads)
                            temp_ex.append(ex_list[ex_idx])
                            break
                    #else:
                    #    break
    
                    if not is_cross:
                        continue
                    nc = read.seq.count('C') + read.seq.count('c') 
                    ng = read.seq.count('g') + read.seq.count('G') 
                    gc_content = float(nc + ng) / float(len(read.seq))
                    for (junc_start, junc_end) in junc_list:
                        if junc_start - r_start > readlen:
                            r_start = junc_start - (readlen - 16) - 1
                        elif junc_start - r_start >= readlen - 8 or junc_start - r_start <= 8:
                            continue
                        found = False
                        if junc_end - junc_start < 10:
                            counter[0] += 1
                            continue

                        for jj in j_list:
                            (j_st, j_ed) = jj.get_coordinates()
                            if j_st > junc_start or (j_st == junc_start and j_ed > junc_end):
                                break
                            elif j_st < junc_start or (j_st == junc_start and j_ed < junc_end):
                                continue
                            elif junc_start == j_st and junc_end == j_ed:
                                ''' update junction and add to list'''
                                counter[3] += 1
                                jj.update_junction_read(exp_index, nreads, r_start, gc_content, unique)
                                if not (junc_start, '5prime', jj) in junctions:
                                    junctions.append((junc_start, '5prime', jj))
                                    junctions.append((junc_end, '3prime', jj))
                                found = True
                                break
                            #end elif junc_start == ...
                        #end for jj in j_list

                        if not found:
                            ''' update junction and add to list'''
                            junc = None
                            for (coord, t, junc) in junctions:
                                if junc.start == junc_start and junc.end == junc_end:
                                    junc.update_junction_read(exp_index, nreads, r_start, gc_content, unique)
                                    if not (junc_start, '5prime', junc) in junctions:
                                        junctions.append((junc_start, '5prime', junc))
                                        junctions.append((junc_end, '3prime', junc))
                                    break
                                #end if (j.start) == ...
                            #end for (coord,t,j) ...
                            if junc is None:
                                '''mark a new junction '''
                                counter[4] += 1
                                junc = Junction(junc_start, junc_end, None, None, gne, readN=nreads)
                                junc.update_junction_read(exp_index, nreads, r_start, gc_content, unique)
                                junctions.append((junc_start, '5prime', junc))
                                junctions.append((junc_end, '3prime', junc))
                        #end if not found ...
                    #end for junc ...
    #            print "JJJunctions", junctions
            if len(junctions) > 0:
                detect_exons(gne, junctions, None)
            gne.prepare_exons()

    for ss in samfile:
        ss.close()
    gc.collect()

    logging.debug("INVALID JUNC", counter[0])
    logging.debug("READ WRONG GENE", counter[1])
    logging.debug("READ IN GENE",    counter[2])
    logging.debug("READ FOUND JUNC", counter[3])
    logging.debug("READ NEW JUNC",   counter[4])
    logging.debug("READ ALL JUNC",   counter[5])
    logging.debug("Non Unique", non_unique_num)
            
    logging.debug("Skipped genes without exons", skip_gene)
    logging.debug(" Non skipped", non_skip)
    return 


def read_transcript_ucsc(filename):
    print "Reading transcript ucsc format"
    buffer_size = int(10E6)
    file_p = open(filename, 'r')
    text = file_p.readlines(buffer_size)
    all_genes = {}
    while text != []:
        print "New Buffer"
        for t in text:
            t = t.strip()
            if t.startswith('#') or t == '':
                continue
            tab = t.split('\t')
            chrom = tab[0]

            #BED format file is 0 based, so the first base of the chromosome is 0.
            # That means that the first real base of an exon is start+1'''
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

            if not chrom in all_genes:
                all_genes[chrom] = {'+': [], '-': []}
            gn = Gene(gene_id, chrom, strand, start, end)
            gene = gn.is_gene_in_list(all_genes[chrom][strand], gene_id)
            if not gene is None:
                del gn
                gn = gene
            else:
                all_genes[chrom][strand].append(gn)

            trcpt = Transcript(transcript_id, gn, start, end)
            gn.add_transcript(trcpt)
            pre_end = None
            pre_txex = None
            for ii in xrange(nblocks):
                start = ex_start[ii]
                end = ex_end[ii]

                txex = gn.new_annotated_exon(start, end, trcpt)
                trcpt.add_exon(txex)
                junc = gn.exist_junction(pre_end, start)
                if junc is None:
                    junc = Junction(pre_end, start, None, None, gn, annotated=True)
                trcpt.add_junction(junc)
                txex.add_3prime_junc(junc)
                if not pre_end is None:
                    pre_txex.add_5prime_junc(junc)
                pre_end = end
                pre_txex = txex
            #END for ii in range(nblocks)
            if not pre_end is None: 
                junc = Junction(pre_end, None, None, None, gn, annotated=True)
                trcpt.add_junction(junc)
                pre_txex.add_5prime_junc(junc)
            trcpt.sort_in_list()
        #end for t in text
        gc.collect()
        text = file_p.readlines(buffer_size)
    #END WHILE
    file_p.close()

    n_genes = 0
    for chrom in all_genes.keys():
        temp_ex = []
        for strand, gg in all_genes[chrom].items():
            n_genes += len(gg)
            all_genes[chrom][strand] = sorted(gg)
            for gene in all_genes[chrom][strand]:
#                print "COLLAPSING GENE", gene.get_id()
                gene.collapse_exons()
                temp_ex.extend(gene.get_exon_list())
#                print "NUM IR:", len(gene.ir_list)
        print "Calculating gc_content.........",
        majiq_exons.set_exons_gc_content(chrom, temp_ex)
        print "Done."

    print "NUM_GENES", n_genes
    return all_genes


def read_bed_pcr(filename, list_genes):
    
    input_f = open(filename, 'r')
    readlines = input_f.readlines()
    alt_exon = []
    pre_chrom = ''
    gene_list = {}
    lnum = 0
    while lnum < len(readlines):
        more = True
        event = {}
        while more:
            rl = readlines[lnum]
            tab = rl.strip().split()
            t = tab[3].split('_')
            event['name'] = '_'.join(t[:-1])
            reg = t[-1]
            if reg == 'C2':
                more = False
            event['chrom'] = tab[0]
            event['strand'] = tab[5]
            event[reg] = [int(tab[1]), int(tab[2])]
            score = tab[4].split('|')[0]
            if score == '?':
                score = 0
            score = float(score)
            if reg in ['A', 'A2']:
                alt_exon = event[reg]
            lnum += 1

        chrom = event['chrom']
        strand = event['strand']

        if chrom != pre_chrom:
            try:
                gene_list = list_genes[chrom]
            except KeyError:
                continue

            pre_chrom = chrom
            idx = {'+': 0, '-': 0}
        name = event['name']

        if strand == '-':
            region_list = ('C2', 'C1')
        else:
            region_list = ('C1', 'C2')

        for reg in region_list:
            exon_start = event[reg][0]
            exon_end = event[reg][1]

            while idx[strand] < len(gene_list[strand]):
                gn = gene_list[strand][idx[strand]]
                (g_start, g_end) = gn.get_coordinates()
                if exon_end < g_start:
                    break
                elif exon_start > g_end:
                    idx[strand] += 1
                    continue
                ex = gn.exist_exon(exon_start, exon_end)
                if ex is None:
                    break
                ex.set_pcr_score(name, score, alt_exon)

                break


gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)


def __parse_gff_attributes(attribute_string):
    """Parse the GFF3 attribute column and return a dict
    :param attribute_string:
    """  #
    if attribute_string == ".":
        return {}
    ret = {}
    for attribute in attribute_string.split(";"):
        key, value = attribute.split("=")
        key = urllib.unquote(key)
        if key in ret:
            key = 'extra_%s' % key
            if not key in ret:
                ret[key] = []
            ret[key].append(urllib.unquote(value))
        else:
            ret[key] = urllib.unquote(value)
    return ret


def __parse_gff3(filename):
    """
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.

    Supports transparent gzip decompression.
    """
    #Parse with transparent decompression
    open_func = gzip.open if filename.endswith(".gz") else open
    with open_func(filename) as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            #If this fails, the file format is not standard-compatible
            assert len(parts) == len(gffInfoFields)
            #Normalize data
            normalized_info = {
                "seqid": None if parts[0] == "." else urllib.unquote(parts[0]),
                "source": None if parts[1] == "." else urllib.unquote(parts[0]),
                "type": None if parts[2] == "." else urllib.unquote(parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else urllib.unquote(parts[6]),
                "phase": None if parts[7] == "." else urllib.unquote(parts[7]),
                "attributes": __parse_gff_attributes(parts[8])
            }
            #Alternatively, you can emit the dictionary here, if you need mutabwility:
            #    yield normalized_info
            yield GFFRecord(**normalized_info)


def _prepare_and_dump(genes, logging=None):
    n_genes = 0
    for chrom in genes.keys():
        temp_ex = []
        for strand, gg in genes[chrom].items():
            n_genes += len(gg)
            genes[chrom][strand] = sorted(gg)
            for gene in genes[chrom][strand]:
                gene.collapse_exons()
                temp_ex.extend(gene.get_exon_list())
        if not logging is None:
            logging.info("Calculating gc_content chromosome %s........." % chrom)
        majiq_exons.set_exons_gc_content(chrom, temp_ex)
        gc.collect()
        temp_dir = "%s/tmp/%s" % (mglobals.outDir, chrom)
        create_if_not_exists(temp_dir)
        # ipdb.set_trace()
        # objgraph.show_most_common_types(limit=20)
        if not logging is None:
            logging.info("Creating temporal annotation %s" % chrom)
        fname = '%s/annot_genes.pkl' % temp_dir
        dump_bin_file(genes[chrom], fname)

    tmp_chrom = "%s/tmp/chromlist.pkl" % mglobals.outDir
    dump_bin_file(genes.keys(), tmp_chrom)
    if not logging is None:
        logging.debug("Number of Genes", n_genes)



def read_gff(filename, pcr_filename, logging=None):
    """

    :param filename: GFF input filename
    :param pcr_filename: BED file name with the PCR validations
    :param logging: logger object
    :return: :raise RuntimeError:
    """
    all_genes = {}
    gene_id_dict = {}
    trcpt_id_dict = {}
    last_end = {}
    for record in __parse_gff3(filename):
        chrom = record.seqid
        strand = record.strand
        start = record.start
        end = record.end

        if record.type == 'gene':
            gene_name = record.attributes['Name']
            gene_id = record.attributes['ID']
            if not chrom in all_genes:
                all_genes[chrom] = {'+': [], '-': []}
            gn = Gene(gene_id, gene_name, chrom, strand, start, end)

            if gene_id in mglobals.gene_tlb and gn != mglobals.gene_tlb[gene_id]:
                raise RuntimeError('Two Different Genes with the same name %s' % gene_name)
            mglobals.gene_tlb[gene_id] = gn
            all_genes[chrom][strand].append(gn)
            gene_id_dict[record.attributes['ID']] = gn

        elif record.type == 'mRNA' or record.type == 'transcript':
            transcript_name = record.attributes['ID']
            parent = record.attributes['Parent']
            try:
                gn = gene_id_dict[parent]
                trcpt = Transcript(transcript_name, gn, start, end)
                gn.add_transcript(trcpt)
                trcpt_id_dict[record.attributes['ID']] = trcpt
                last_end[record.attributes['ID']] = (None, None)
            except KeyError:
                if not logging is None:
                    logging.Error("Error, incorrect gff. mRNA %s doesn't have valid gene %s"
                                  % (transcript_name, parent))
                raise

        elif record.type == 'exon':
            parent_tx_id = record.attributes['Parent']
            try:
                parent_tx = trcpt_id_dict[parent_tx_id]
                gn = parent_tx.get_gene()
                txex = gn.new_annotated_exon(start, end, parent_tx)
                parent_tx.add_exon(txex)

            except KeyError:
                if not logging is None:
                    logging.info("Error, incorrect gff. exon %s doesn't have valid mRNA %s" % (record.attributes['ID'],
                                                                                               parent_tx_id))
        #end elif
    #end for
    for kk, trcpt in trcpt_id_dict.items():

        exon_list = trcpt.prepare_exon_list()
        gn = trcpt.get_gene()
        pre_end = None
        pre_txex = None
        for ex in exon_list:
            start, end = ex.get_coordinates()
            junc = gn.new_annotated_junctions(pre_end, start, trcpt)
            ex.add_3prime_junc(junc)
            if not pre_txex is None:
                pre_txex.add_5prime_junc(junc)
            pre_end = end
            pre_txex = ex

        junc = gn.new_annotated_junctions(pre_end, None, trcpt)
        pre_txex.add_5prime_junc(junc)

    #end for
    _prepare_and_dump(all_genes, logging)
    if pcr_filename is not None:
        read_bed_pcr(pcr_filename, all_genes)

    chr_list = all_genes.keys()
    del all_genes
    return chr_list
