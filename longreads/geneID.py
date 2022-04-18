import pandas as pd


#loading the first dataset
df_tsv = pd.read_csv('summary.tsv', delimiter='\t', skiprows=4)
#loading the second dataset
df_gtf = pd.read_csv('collapse3.isoforms.gtf', sep='\t|;',
                        header=None, index_col=False, engine='python')

def extract_gene(df_tsv, df_gtf):

    #extracting the 'gene_id' column
    gene_id_tsv = df_tsv['gene_id']
    gene_id_tsv = set(gene_id_tsv)

    #extracting and formatting the 'gene_id' column
    df_gtf[8] = df_gtf[8].apply(lambda x : x.strip('gene_id').strip()[1:-1])
    gene_id_gtf = set(df_gtf[8])

    #finding the common gene_ids of both dataframes
    common_genes = gene_id_gtf.intersection(gene_id_tsv)
    common_genes = pd.Series(list(common_genes))

    #checking if the intersection is good
    for s in common_genes:
        assert s in gene_id_gtf
        assert s in gene_id_tsv

    #saving all the common gene_ids to a .txt file
    common_genes.to_csv('common_genes.tsv', sep='\t', header=None, index=False)



if __name__== "__main__":
    extract_gene(df_tsv=df_tsv, df_gtf=df_gtf)