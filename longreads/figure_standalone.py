import os

import matplotlib.pyplot as plt
import matplotlib.patches as patches

import majiqv2, flair
from graph import exon


majiq_splicegraph_path = '/slowdata/lrdata/majiq/splicegraph.sql'
majiq_gene_id="gene:ENSG00000109534"
parser = majiqv2.MajiqV2Reader(majiq_splicegraph_path)
parser.parse_splicegraph(majiq_gene_id)

flair_gtf_path = '/slowdata/lrdata/flair/flair_filter_transcripts.gtf'
flair_gene_id = 'ENSG00000109534.16'

save_path = '/tmp/lr_o'
os.makedirs(save_path, exist_ok=True)


padding = 100
transcript_height = 20
transcript_height_padding = 5
majiqColor='r'
flairColor='g'
bothColor='b'
moduleColor='#c3b40840'

print('~~~Parsing Flair~~~')
flairreader = flair.FlairReader(flair_gtf_path)
print('~~~Done Parsing Flair~~~')



def plot(only_in_flair, only_in_majiq, in_flair_and_majiq, filename, module_extent=None):

    fig, ax = plt.subplots(1)

    if not module_extent:
        gene_start, gene_end = parser.extent(majiq_gene_id)
    else:
        gene_start, gene_end = module_extent[0], module_extent[1]


    print(only_in_flair)
    print(only_in_majiq)
    print(in_flair_and_majiq)
    print('---------------------')



    delayed_patches = []
    y = transcript_height_padding
    colors = (flairColor, majiqColor, bothColor)
    for i, transcripts in enumerate((only_in_flair, only_in_majiq, in_flair_and_majiq)):
        color = colors[i]
        for exons in transcripts:

            for exon in exons:
                x = exon.start

                gene_start = min(gene_start, exon.start)
                gene_end = max(gene_end, exon.end)

                width = exon.end - exon.start

                rect = patches.Rectangle((x, y), width, transcript_height, linewidth=1,
                                         edgecolor=color, facecolor="none")
                delayed_patches.append(rect)

            y += (transcript_height_padding + transcript_height)

    if module_extent:
        rect = patches.Rectangle((module_extent[0], transcript_height_padding/2.0), module_extent[1]-module_extent[0], y-(transcript_height_padding/2.0), linewidth=1,
                                 edgecolor='none', facecolor=moduleColor)
        ax.add_patch(rect)
    else:
        for module_idx in range(parser.getNumModules()):
            module_extent = parser.moduleExtent(module_idx)
            rect = patches.Rectangle((module_extent[0], transcript_height_padding/2.0), module_extent[1]-module_extent[0], y-(transcript_height_padding/2.0), linewidth=1,
                                     edgecolor='none', facecolor=moduleColor)
            ax.add_patch(rect)

    for rect in delayed_patches:
        ax.add_patch(rect)

    plt.xlim([gene_start-padding, gene_end+padding])
    plt.ylim([0, y])
    plt.title(f"{majiq_gene_id} / {flair_gene_id}")

    #plt.show()
    plt.savefig(os.path.join(save_path, filename))



# Plot for the entire gene
# here, we match transcripts exactly between majiq and flair
flair_exons = set(x[0] for x in flairreader.gene(flair_gene_id))
majiq_exons = set(x[0] for x in parser.getAllPaths())
only_in_flair = flair_exons.difference(majiq_exons)
only_in_majiq = majiq_exons.difference(flair_exons)
in_flair_and_majiq = flair_exons.intersection(majiq_exons)
plot(only_in_flair, only_in_majiq, in_flair_and_majiq, f'{majiq_gene_id}_module_combined.png')

for module_idx in range(parser.getNumModules()):

    majiq_module_extent = parser.moduleExtent(module_idx)

    """
    for in-module, by default the exons we receive from majiq start/end are technically not part of the module
    UNLESS, they have different start/end. As a simple way to deal with this, for matching purposes, we will trim all
    exon coordinates at the start/end of the module to match the module coordinates
    
    """

    ord_flair_exons = tuple(x[0] for x in flairreader.gene(flair_gene_id, extent=majiq_module_extent))

    ord_majiq_exons = tuple(x[0] for x in parser.getAllPaths(module_idx=module_idx))
    flair_exons = set()
    for transcript in ord_flair_exons:
        flair_exons.add((exon(max(majiq_module_extent[0], e.start), min(majiq_module_extent[1], e.end)) for e in transcript))
    majiq_exons = set()
    for transcript in ord_majiq_exons:
        majiq_exons.add((exon(max(majiq_module_extent[0], e.start), min(majiq_module_extent[1], e.end)) for e in transcript))


    only_in_flair = flair_exons.difference(majiq_exons)
    only_in_majiq = majiq_exons.difference(flair_exons)
    in_flair_and_majiq = flair_exons.intersection(majiq_exons)


    plot(only_in_flair, only_in_majiq, in_flair_and_majiq, f'{majiq_gene_id}_module_{module_idx}.png', module_extent=majiq_module_extent)






