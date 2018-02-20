from majiq.src.internals.io_bam cimport IOBam
from majiq.src.multiproc import QueueMessage
from majiq.src.constants import *
from libcpp.string cimport string
from libcpp.set cimport set
from libcpp.map cimport map
from majiq.src.sample import create_lt
import majiq.src.io as majiq_io
from majiq.grimoire.exon import detect_exons
from majiq.src.config import Config
from majiq.src.internals.junction cimport Junction

cdef _find_new_junctions(list file_list, int chunk, object conf, object logger):
    cdef str gne_id
    cdef dict gene_obj
    cdef string chrom, cs1
    cdef map[string, set[string]] set_junctions

    majiq_config = Config()
    list_exons = {}

    logger.info('Reading DB')


    for gne_id, gene_obj in conf.genes_dict.items():
        list_exons[gne_id] = []
        dict_junctions = {}
        chrom = conf.genes_dict[gne_id]['chromosome'].encode('utf-8')

        majiq_io.from_matrix_to_objects(gne_id, conf.elem_dict[gne_id], dict_junctions, list_exons[gne_id])
        # if chrom not in set_junctions:
        #
        for yy in dict_junctions.keys():
            for strnd in ['+', '-', '.']:
                cs1 = ('%s:%s:%s-%s' % (chrom, '+', yy[0], yy[1])).encode('utf-8')
                set_junctions[chrom].insert(cs1)

        detect_exons(dict_junctions, list_exons[gne_id])

    td = create_lt(conf.genes_dict)
    for exp_name, is_junc_file, name in file_list:
        fname = '%s/%s.%s' % (majiq_config.sam_dir, exp_name, SEQ_FILE_FORMAT)
        logger.info('READ JUNCS from %s' % fname)
        bbb = SeqParse(fname, majiq_config.strand_specific[exp_name])
        bbb.check_junctions(td, set_junctions, majiq_config.strand_specific[exp_name], list_exons,
                            conf.queue, name, logger)

        # read_juncs(fname, is_junc_file, list_exons, conf.genes_dict, td, dict_junctions,
        #            majiq_config.strand_specific[exp_name], conf.queue, gname=name, logger=logger)



cdef class SeqParse:
    cdef IOBam c_iobam      # hold a C++ instance which we're wrapping
    def __cinit__(self, str bam1, int strandness1, str regions=None):
        cdef bytes py_bytes = bam1.encode()
        cdef char* c_bam1 = py_bytes
        if regions is None:
            self.c_iobam = IOBam(c_bam1, strandness1)
        # else:
        #     self.c_iobam = IOBam(c_bam1, strandness1, c_regions)

    cdef find_junctions(self):
        self.c_iobam.find_junctions()

    cdef check_junctions(self, dict_gtrees, map[string, set[string]] set_junctions, stranded, dict_exons,
                         queue, gname, logger):
        #cdef string chrom
        cdef list junc

        logger.info("START READING")
        with nogil:
            self.c_iobam.set_filters(set_junctions)
            self.c_iobam.find_junctions()
        logger.info("DONE READING...")

        #for it in self.c_iobam.get_dict():
        for it in self.c_iobam.junc_dict_:
            chrom = it.second.chrom.decode('utf-8')
            if chrom not in dict_gtrees: continue
            found = False
            #for junc in sorted(jj_set):
            possible_genes = []
            for gobj in dict_gtrees[chrom].search(it.second.start, it.second.end):

                if (stranded and gobj.data[1] != it.second.strand) or (gobj.end < it.second.end and
                                                                       gobj.start > it.second.start):
                    continue

                gid = gobj.data[0]

                start_sp = [jj.start for ex in dict_exons[gid] for jj in ex.ob if jj.start > 0 and jj.end > 0]
                end_sp = [jj.end for ex in dict_exons[gid] for jj in ex.ib if jj.start > 0 and jj.end > 0]

                if it.second.start in start_sp or it.second.end in end_sp:
                    found = True
                    # print((gid, junc[0], junc[1], gname))
                    qm = QueueMessage(QUEUE_MESSAGE_BUILD_JUNCTION, (gid, it.second.start, it.second.end, gname), 0)
                    queue.put(qm, block=True)
                else:
                    possible_genes.append(gid)

            if not found and len(possible_genes) > 0:
                for gid in possible_genes:
                    # print((gid, junc[0], junc[1], gname))
                    qm = QueueMessage(QUEUE_MESSAGE_BUILD_JUNCTION, (gid, it.second.start, it.second.end, gname), 0)
                    queue.put(qm, block=True)



## OPEN API FOR PYTHON

def find_new_junctions(list file_list, int chunk, object conf, object logger):
    _find_new_junctions(file_list, chunk, conf, logger)

