import argparse
import flairParser
import pickle

parser = argparse.ArgumentParser(description='Converter from various long read non-majiq output softwares to voila.lr file')
parser.add_argument('-i', '--input-file', type=str, required=False,
                    help='please provide an input text file from a long read program, for example, a gtf')
parser.add_argument('-o', '--output-file', type=str, required=False,
                    help='the path to write the resulting voila file to')

args = parser.parse_args()

print('~~~Parsing Flair~~~')
flairreader = flairParser.FlairReader(args.input_file)
print('~~~Done Parsing Flair~~~')

all_genes = {}
for gene_id in flairreader.gene_ids:
    all_genes[gene_id] = []
    for transcript in flairreader.get_exons(gene_id):
        all_genes[gene_id].append([(abs(x.start), abs(x.end)) for x in transcript])


with open(args.output_file, 'wb') as f:
    pickle.dump(all_genes, f)

