import os

import matplotlib.pyplot as plt
import matplotlib.patches as patches

import majiqv2, flair



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





def plot(flair_exons, majiq_exons, filename, module_extent=None):

    fig, ax = plt.subplots(1)

    if not module_extent:
        gene_start, gene_end = parser.extent(majiq_gene_id)
    else:
        gene_start, gene_end = module_extent[0], module_extent[1]


    only_in_flair = flair_exons.difference(majiq_exons)
    only_in_majiq = majiq_exons.difference(flair_exons)
    in_flair_and_majiq = flair_exons.intersection(majiq_exons)



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




flair_exons = set(x[0] for x in flair.FlairReader.parse_gtf(flair_gtf_path, flair_gene_id))
majiq_exons = set(x[0] for x in parser.getAllPaths())
plot(flair_exons, majiq_exons, f'{majiq_gene_id}_module_combined.png')

for module_idx in range(parser.getNumModules()):

    majiq_module_extent = parser.moduleExtent(module_idx)

    flair_exons = set(x[0] for x in flair.FlairReader.parse_gtf(flair_gtf_path, flair_gene_id, extent=majiq_module_extent))
    majiq_exons = set(x[0] for x in parser.getAllPaths(module_idx=module_idx))

    plot(flair_exons, majiq_exons, f'{majiq_gene_id}_module_{module_idx}.png', module_extent=majiq_module_extent)






