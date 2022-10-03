import os

import matplotlib.pyplot as plt
import matplotlib.patches as patches

import majiqv2, flairParser
from graph import exon
from config import get_args
from tool_comparison import ToolComparer

args = get_args()

#majiq_splicegraph_path = '/slowdata/lrdata/majiq/splicegraph.sql'
#majiq_gene_id="gene:ENSG00000109534"
majiq_splicegraph_path = args.majiq_splicegraph_path
majiq_gene_id = args.gene_id

majiq_parser = majiqv2.MajiqV2Reader(majiq_splicegraph_path)
majiq_parser.parse_splicegraph(majiq_gene_id)

#flair_gtf_path = '/slowdata/lrdata/flair/flair_filter_transcripts.gtf'
#flair_gene_id = 'ENSG00000109534.16'
flair_gtf_path = args.flair_gtf_path

if args.gene_id_flair:
    flair_gene_id = args.gene_id_flair
else:
    flair_gene_id = majiq_gene_id

save_path = args.output_path
os.makedirs(save_path, exist_ok=True)


padding = 100
transcript_height = 20
transcript_height_padding = 5
majiqColor='r'
flairColor='g'
bothColor='b'
moduleColor='#c3b40840'

print('~~~Parsing Flair~~~')
flairreader = flairParser.FlairReader(flair_gtf_path)
print('~~~Done Parsing Flair~~~')



def plot(only_in_flair, only_in_majiq, in_flair_and_majiq, filename, module_extent=None, module_extents=None):

    fig, ax = plt.subplots(1)

    if not module_extent:
        gene_start = float('inf')
        gene_end = float('-inf')
        for transcript_type in (only_in_flair, only_in_majiq,):
            for transcript in transcript_type:
                for _exon in transcript:
                    gene_start = min(gene_start, abs(_exon.start))
                    gene_end = max(gene_end, abs(_exon.end))
        for transcript_type in (in_flair_and_majiq,):
            for m_transcript, f_transcript in transcript_type:
                for _exon in m_transcript:
                    gene_start = min(gene_start, abs(_exon.start))
                    gene_end = max(gene_end, abs(_exon.end))
                for _exon in f_transcript:
                    gene_start = min(gene_start, abs(_exon.start))
                    gene_end = max(gene_end, abs(_exon.end))
    else:
        gene_start, gene_end = module_extent[0], module_extent[1]

    import pprint
    print(gene_start, gene_end)
    print(only_in_flair)
    print(only_in_majiq)
    print(in_flair_and_majiq)
    print('---------------------')

    delayed_patches = []
    y = transcript_height_padding
    colors = (flairColor, bothColor, majiqColor)
    for i, transcripts in enumerate((only_in_flair, in_flair_and_majiq, only_in_majiq)):
        color = colors[i]

        if i == 1:
            transcripts = (x[1] for x in transcripts)  # use the majiq transcripts when both are provided

        for exons in transcripts:

            for _exon in exons:

                if _exon.start == -1:
                    _exon = exon(_exon.end, _exon.end)

                if _exon.end == -1:
                    _exon = exon(_exon.start, _exon.start)

                if _exon.start < 0:
                    _exon = exon(gene_start, _exon.end)

                if _exon.end < 0:
                    _exon = exon(_exon.start, gene_end)

                x = _exon.start

                gene_start = min(gene_start, _exon.start)
                gene_end = max(gene_end, _exon.end)

                width = _exon.end - _exon.start

                rect = patches.Rectangle((x, y), width, transcript_height, linewidth=1,
                                         edgecolor=color, facecolor="none")
                delayed_patches.append(rect)

            y += (transcript_height_padding + transcript_height)

    if module_extent:
        rect = patches.Rectangle((module_extent[0], transcript_height_padding/2.0), module_extent[1]-module_extent[0], y-(transcript_height_padding/2.0), linewidth=1,
                                 edgecolor='none', facecolor=moduleColor)
        ax.add_patch(rect)
    else:
        for _module_extent in module_extents:
            rect = patches.Rectangle((_module_extent.start, transcript_height_padding/2.0), _module_extent.end-_module_extent.start, y-(transcript_height_padding/2.0), linewidth=1,
                                     edgecolor='none', facecolor=moduleColor)
            ax.add_patch(rect)

    for rect in delayed_patches:
        ax.add_patch(rect)

    plt.xlim([gene_start-padding, gene_end+padding])
    plt.ylim([0, y])
    plt.title(f"{majiq_gene_id} / {flair_gene_id}")

    #plt.show()
    plt.savefig(os.path.join(save_path, filename))



flair_exons = flairreader.get_exons(flair_gene_id, majiq_module_extent=None, modules=None)
print('!', flair_exons)

modules_list = [majiq_parser.moduleExtent(i) for i in range(majiq_parser.getNumModules())]
modules_list = flairreader.extend_modules(modules_list, flairreader.get_exons(flair_gene_id))

print('modules', modules_list)

majiq_exons, majiq_denovo, majiq_has_reads = majiq_parser.allpaths_data(
    modules=None,
    module_idx=None,
    max_paths=args.max_paths,
    majiq_module_extent=None,
)

tc = ToolComparer(args)
only_in_flair, only_in_majiq, in_flair_and_majiq = tc.compare_fuzzy(flair_exons, majiq_exons, args.fuzziness5, args.fuzziness3)
# Plot for the entire gene
# here, we match transcripts exactly between majiq and flair

plot(only_in_flair, only_in_majiq, in_flair_and_majiq, f'{majiq_gene_id}_module_combined.png', module_extents=modules_list)






for module_idx, module in enumerate(modules_list):



    """
    for in-module, by default the exons we receive from majiq start/end are technically not part of the module
    UNLESS, they have different start/end. As a simple way to deal with this, for matching purposes, we will trim all
    exon coordinates at the start/end of the module to match the module coordinates
    
    """

    flair_exons = flairreader.get_exons(flair_gene_id, majiq_module_extent=module, modules=None)

    majiq_exons, majiq_denovo, majiq_has_reads = majiq_parser.allpaths_data(
        modules=None,
        module_idx=module_idx,
        max_paths=args.max_paths,
        majiq_module_extent=module
    )
    only_in_flair, only_in_majiq, in_flair_and_majiq = tc.compare_fuzzy(flair_exons, majiq_exons, args.fuzziness5, args.fuzziness3)


    plot(only_in_flair, only_in_majiq, in_flair_and_majiq, f'{majiq_gene_id}_module_{module_idx}.png', module_extent=module)






