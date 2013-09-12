import argparse

class RNA_read:

    def __init__ ( self, start, end, junc_list, readlen, count, gc,unique):
        self.start = start
        self.end = end
        self.readlen = readlen
        self.count = count
        self.junc_list = junc_list
        self.GCcontent = gc
        self.unique = unique

    def get_junctions_info(self):
        return self.junc_list

    def get_unique(self):
        return self.unique

    def get_read_count(self):
        return self.count

    def get_coordinates(self):
        return (self.start,self.end)

    def get_GCcontent(self):
        return self.GCcontent

    def get_gene(self):
        return self.gene

    def set_gene(self,gene):
        self.gene = gene
        return

def parse_CIGAR (string):
    num = 0
    lst = []
    for ii in string:
        if ii.isdigit():
            num = 10*num +int(ii)
        else:
            lst.append( (ii,num) )
            num = 0
    return lst


'''
    @parse_read
    @param read_string: This string has the SAM line that define
        the RNAseq read

    @output RNA_read object
'''

def parse_read(read_lst, readlen):

    jlist = []
    start = int(read_lst[3])

    cigar = parse_CIGAR(read_lst[5])

    options = read_lst[11:]
    next_read = False
    n_reads = 1
    unique = True
    for op in options:
        op_tab = op.split(':')
        if op_tab[0]== 'NH' and int(op_tab[2]) > 1:
            #NH determines if the read is uniquely mapped or not
            unique = False
#            next_read = True
#            break 
        elif op_tab[0] == 'HI':
            n_reads = int(op_tab[2])
        elif op_tab[0] == 'jI':
            '''
              This part will parse the jI from STAR
                # Column 17: jI:B:I,Start1,End1,Start2,End2,... Start and End of introns for all junctions (1- based)
                # jI:B:I,-1 means that the read doesn't cross any junction
            '''
            jc_coord = op_tab[2].split(',')
            if jc_coord[1] == '-1':
                next_read = True
                break 
            else:
                #TODO: check if the read pass more than one junction

                for idx in range(1,len(jc_coord),2):
                    junc_start = int(jc_coord[idx])-1
                    junc_end = int(jc_coord[idx+1])+1
#                    if junc_end == 46024065 :
#                        print "STAT 1",read_lst
                    jlist.append((junc_start,junc_end))
                #end for idx ...
            #end else ...
        #end elif jI ...
    #end for op ...
    #old_start,old_end = prev_read.get_coordinates()

    if next_read :
        read = n_reads
        is_read = False
    else:
        is_read = True
        for ii,jj in reversed(cigar):
            if ii == 'M':
                m = jj
                break
        #end for
        '''end is calculated by adding the last matching set on CIGAR string, to the last 3'ss specified on jI'''
        end = jlist[-1][1]+m
        gc = float(read_lst[9].count('C') + read_lst[9].count('G'))/float(len(read_lst[9]))
        read = RNA_read( start, end, jlist, readlen, n_reads, gc,unique )
        #end else ...

    return (is_read,read)
