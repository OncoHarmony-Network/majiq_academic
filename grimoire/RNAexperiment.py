import copy


class RNAexperiment(object):

    def __init__(self, tissue,replica_num,genome, gene_list):
        self.tissue = tissue
        self.replica = replica_num
        self.genome = genome
        self.read_list = []
        self.gene_list = copy.deepcopy(gene_list)


    def add_gene(self, gene):
        chr = gene.get_chromosome()
        if not chr in self.gene_list:
            self.gene_list[chr] = []
        
        self.gene_list[chr].append(gene)
        return gene


    def get_gene_list(self):
        return self.gene_list

   


