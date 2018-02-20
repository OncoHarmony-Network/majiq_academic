from majiq.src.internals.io_bam cimport IOBam
from majiq.src.multiproc import QueueMessage
from majiq.src.constants import *
from libcpp.string cimport string
from libcpp.set cimport set
from majiq.src.internals.junction cimport Junction

cdef class SeqParse:
    cdef IOBam c_iobam      # hold a C++ instance which we're wrapping
    def __cinit__(self, str bam1, int strandness1, str regions=None):
        cdef bytes py_bytes = bam1.encode()
        cdef char* c_bam1 = py_bytes
        if regions is None:
            self.c_iobam = IOBam(c_bam1, strandness1)
        # else:
        #     self.c_iobam = IOBam(c_bam1, strandness1, c_regions)

    def find_junctions(self):
        self.c_iobam.find_junctions()

    def check_junctions(self, dict_gtrees, dict_genes, junctions, stranded, dict_exons, queue, gname, logger):
        # cdef map[string, Junction].iterator it

        cdef set[string] set_junctions
        cdef set[string] set_chrom  #= set[string]({xx['chromosomes'] for xx in dict_genes.values()})
        # cdef set[string].iterator set_it = set_junctions.begin()
        cdef string cs1


        for gene_id, xx in junctions.items() :
            cs1 = dict_genes[gene_id]['chromosome'].encode('utf-8')
            set_chrom.insert(cs1)
            if stranded != UNSTRANDED:
                for yy in xx.values():
                    cs1 = ('%s:%s:%s-%s' % (dict_genes[gene_id]['chromosome'], dict_genes[gene_id]['strand'],
                                           yy.start, yy.end)).encode('utf-8')
                    set_junctions.insert(cs1)

            else:
                for yy in xx.values():
                    cs1 = ('%s:.:%s-%s' % (dict_genes[gene_id]['chromosome'], yy.start, yy.end)).encode('utf-8')
                    set_junctions.insert(cs1)
        logger.info("START READING")
        self.c_iobam.set_filters(set_chrom, set_junctions)
        self.c_iobam.find_junctions()
        logger.info("DONE READING...")

        #for it in self.c_iobam.get_dict():
        for it in self.c_iobam.junc_dict_:
            chrom = it.second.chrom.decode('utf-8')
            if chrom not in dict_gtrees: continue
            junc = [it.second.start, it.second.end, it.second.strand]

            found = False
            #for junc in sorted(jj_set):
            possible_genes = []
            for gobj in dict_gtrees[chrom].search(junc[0], junc[1]):

                if (stranded and gobj.data[1] != junc[2]) or (gobj.end < junc[1] and gobj.start > junc[0]):
                    continue

                gid = gobj.data[0]

                start_sp = [jj.start for ex in dict_exons[gid] for jj in ex.ob if jj.start > 0 and jj.end > 0]
                end_sp = [jj.end for ex in dict_exons[gid] for jj in ex.ib if jj.start > 0 and jj.end > 0]

                if junc[0] in start_sp or junc[1] in end_sp:
                    found = True
                    # print((gid, junc[0], junc[1], gname))
                    qm = QueueMessage(QUEUE_MESSAGE_BUILD_JUNCTION, (gid, junc[0], junc[1], gname), 0)
                    queue.put(qm, block=True)
                else:
                    possible_genes.append(gid)

            if not found and len(possible_genes) > 0:
                for gid in possible_genes:
                    # print((gid, junc[0], junc[1], gname))
                    qm = QueueMessage(QUEUE_MESSAGE_BUILD_JUNCTION, (gid, junc[0], junc[1], gname), 0)
                    queue.put(qm, block=True)



