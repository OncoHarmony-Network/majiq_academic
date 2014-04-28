#!/usr/bin/python
import os
import utils.utils_voila as utils_voila

__author__ = 'abarrera'


def table_marks_set(size):
    # TODO: Customizable by script options
    ideal_set = (10, 20, 50, 100)
    index = 0
    for mark in ideal_set:
        if size < ideal_set[index]:
            break
        index += 1

    return ideal_set[0:index]


def mocked_get_array_bins():
    pass


def _render_template(output_dir, output_html, majiq_output, type_summary, threshold):
    """
    Rendering the summary template to create the HTML file.

    @param output_dir: directory where the HTML will be created
    @param majiq_output: event list from MAJIQ output used as input in Voila
    @param type_summary: defines summary template used
    """
    from jinja2 import Environment, FileSystemLoader
    env = Environment(loader=FileSystemLoader("templates/"))
    sum_template = env.get_template(type_summary + "_summary_template.html")

    if not output_dir.endswith('/'):
        output_dir += '/'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    voila_output = open(output_dir+output_html, 'w+')

    if type_summary == 'single':
        voila_output.write(sum_template.render(eventList=majiq_output['event_list'],
                                               tableMarks=table_marks_set(len(majiq_output['event_list'])),
                                               metadata=majiq_output['metadata_pre']))
    # TODO: Add experiment 1 and 2 info for the template
    elif type_summary == 'delta':
        voila_output.write(sum_template.render(eventList=majiq_output['event_list'],
                                               tableMarks=table_marks_set(len(majiq_output['event_list'])),
                                               expList=majiq_output['experiments_info'],
                                               threshold=threshold
                                               ))

    elif type_summary == 'lsv_single':
        voila_output.write(sum_template.render(eventList=majiq_output['event_list'][:1],
                                               tableMarks=table_marks_set(len(majiq_output['event_list'])),
                                               arrayBins=[majiq_output['event_list'][0].bins,
                                                          majiq_output['event_list'][1].bins,
                                                          majiq_output['event_list'][2].bins],
                                               arrayQuartiles=[majiq_output['event_list'][0].quartiles,
                                                               majiq_output['event_list'][1].quartiles,
                                                               majiq_output['event_list'][2].quartiles],
                                               arrayMeans=[majiq_output['event_list'][0].mean_psi,
                                                           majiq_output['event_list'][1].mean_psi,
                                                           majiq_output['event_list'][2].mean_psi],
                                               arrayNames=['PSI1', 'PSI2', 'PSI3']
                                               ))
    else:
        print "summary type not recognized %s." % type_summary
        import sys
        sys.exit(1)

    voila_output.close()


def create_summary(majiq_bins_file, output_dir, meta_preprocess, meta_postprocess, type_summary, threshold, confidence=.95):
    """This method generates an html summary from a majiq output file"""
    output_html = os.path.splitext(os.path.split(majiq_bins_file)[1])[0] + "_" + type_summary + "_" + str(threshold) + '.html'

    if type_summary == 'single':
        majiq_output = utils_voila.get_single_exp_data(majiq_bins_file, meta_preprocess, meta_postprocess, confidence)
    elif type_summary == 'delta':
        majiq_output = utils_voila.get_delta_exp_data(majiq_bins_file, meta_postprocess, confidence, threshold)
    elif type_summary == 'lsv_single':
        majiq_output = utils_voila.get_lsv_single_exp_data(majiq_bins_file, meta_preprocess, confidence)

    _render_template(output_dir, output_html, majiq_output, type_summary, threshold)
    return


def main():

    import argparse
    parser = argparse.ArgumentParser(description='Generate an HTML summary of Majiq output.')
    parser.add_argument('-b', '--bins', metavar='majiq_output.pickle', dest='majiq_bins', type=str, required=True, help='Pickle file with the bins produced by Majiq.')
    parser.add_argument('-o', '--output', metavar='output_dir', dest='output_dir', type=str, required=True, help='Output directory where the files will be placed.')
    parser.add_argument('--meta-pre', metavar='metadata_pre.majiq', dest='meta_preprocess', type=str, help='Metadata preprocess.')
    parser.add_argument('--event-names', metavar='event_names.majiq', dest='event_names', type=str, help='Event names.')
    parser.add_argument('--key-plots', metavar='keysplots.pickle', dest='keys_plots', type=str, help='Heatmap plots.')
    parser.add_argument('-t', '--type', type=str, choices=['single', 'delta', 'lsv_single'], dest='type_summary', default='single', help='Type of summary generated.')
    parser.add_argument('--threshold', type=float, dest='threshold', default=0.2, help='Probability threshold used to sum the accumulative probability of inclusion/exclusion.')

    # parser.add_argument('-c', '--confidence', metavar=0.95, dest='confidence', type=float,
    #                     default=0.95, help='Percentage of confidence required (by default, 0.95).')
    # TODO: add options to filter the output summary file by: genes, events, etc.

    args = parser.parse_args()
    create_summary(args.majiq_bins, args.output_dir, args.meta_preprocess, {'names': args.event_names, 'keys_plots': args.keys_plots}, args.type_summary, args.threshold)


if __name__ == '__main__':
    main()