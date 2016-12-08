from os import path

from voila.run_voila_utils import table_marks_set

from voila import io_voila
from voila.utils.run_voila_utils import get_env, get_summary_template, copy_static
from voila.utils.voila_log import voila_log


def collect_sample_names(sample_names):
    max_length = 10
    rtn_list = []

    for sample_name in sample_names:
        if len(sample_name) > max_length:
            trunc_sample_name = sample_name[:max_length - 3] + '...'
        else:
            trunc_sample_name = sample_name

        rtn_list.append({'truncated': trunc_sample_name, 'full': sample_name})

    return rtn_list


def conditional_table(args):
    log = voila_log()
    output_html = "%s_%s_comp_table_%.2f.html" % (args.cond_pair[0], args.cond_pair[1], args.threshold_change)

    lsvs_dict = io_voila.load_dpsi_tab(args)
    log.info("LSVs added to the table: %d" % len(lsvs_dict.keys()))

    env = get_env()
    sum_template = get_summary_template(args, env)

    sample_names = collect_sample_names(args.sample_names)

    with open(path.join(args.output, output_html), 'w') as voila_output:
        voila_output.write(
            sum_template.render(
                lsvs=lsvs_dict,
                sample_names=sample_names,
                tableMarks=table_marks_set(len(lsvs_dict)),
                cond_pair=args.cond_pair,
                thres=args.threshold_change
            )
        )

    copy_static(args)
