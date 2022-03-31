import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Figure creator for majiq transcripts comparison')
    parser.add_argument('--gene-id', type=str, required=True,
                        help='Gene Id for Majiq gene')
    parser.add_argument('--gene-id-flair', type=str, required=False,
                        help='Gene id to use for flair input file, if omitted, will use the majiq gene id')
    parser.add_argument('--majiq-splicegraph-path', type=str, required=True,
                        help='path to majiq splicegraph output file')
    parser.add_argument('--flair-gtf-path', type=str, required=True,
                        help='path to flair output file')
    parser.add_argument('--output-path', type=str, required=True,
                        help='path to place output files in, will be created if it does not exist')

    args = parser.parse_args()
    return args