from voila.api import Voila


def create_expected_psi(experiment_names, junction_ids):
    expected_psis = []

    for _ in experiment_names:
        expected_psi = []

        for _ in junction_ids:
            expected_psi.append(.3)

        expected_psis.append(expected_psi)

    return expected_psis


def create_median():
    return [.2 for _ in range(40)]


if __name__ == "__main__":
    # with SpliceGraphs('/Users/cjgreen/Development/small_test/majiq_build/splicegraph.hdf5', 'r') as sg:
    #     for gene in sg.get_genes():
    #         for exon in gene.exons:
    #             if gene.gene_id == 'ENSMUSG00000067629':
    #                 print(gene.gene_id, exon.exon_type)



    with Voila('/Users/cjgreen/Development/small_test/majiq_deltapsi/Adr_Cer.deltapsi.voila', 'r') as v:
        for lsv in v.get_voila_lsvs():
            # print(lsv.trunc_means_psi1)
            print(lsv.means_psi1)
            print(len(lsv.junctions))
            print('----------------------')
            exit(1)