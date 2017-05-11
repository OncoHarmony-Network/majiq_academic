from os import path

from voila.utils.run_voila_utils import get_env, copy_static, get_output_html


def lsv_thumbnails(args):
    """
    Render html output for lsv thumbnails.
    :param args: command line arguments
    :return: None
    """
    output_dir = args.output
    output_html = get_output_html(args)

    env = get_env()
    template_file_name = args.type_analysis.replace("-", "_") + "_summary_template.html"
    sum_template = env.get_template(template_file_name)

    with open(path.join(output_dir, output_html), 'w') as voila_output:
        voila_output.write(
            sum_template.render(
                lsv_list=args.lsv_types,
                collapsed=args.collapsed
            )
        )

    copy_static(args)
