import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Figure creator for majiq transcripts comparison')
    parser.add_argument('--gene-id', type=str, required=False,
                        help='process only this Gene Id for Majiq gene')
    parser.add_argument('--gene-id-flair', type=str, required=False,
                        help='Gene id to use for flair input file, if omitted, will use the majiq gene id')
    parser.add_argument('--gene-ids-file', type=str, required=False,
                        help='Path to a file to read a list of specific gene_ids from')
    parser.add_argument('--majiq-splicegraph-path', type=str, required=True,
                        help='path to majiq splicegraph output file')
    parser.add_argument('--flair-gtf-path', type=str, required=True,
                        help='path to flair output file')
    parser.add_argument('--output-path', type=str, required=True,
                        help='path to place output files in, will be created if it does not exist')
    parser.add_argument('-v', '--verbose', action='count', default=0)
    parser.add_argument('--per-module', action='store_true',
                        help='Show additional verbose output')
    parser.add_argument('--fuzziness', type=int, default=0,
                        help='Reduce fuzziness of long-read sequencing')
    args = parser.parse_args()
    return args