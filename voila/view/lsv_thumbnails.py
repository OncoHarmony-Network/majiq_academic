from os import path

from voila.view.html import Html


class LsvThumbnails():
    def __init__(self, args):
        lsv_thumbnails(args)

    # @classmethod
    # def arg_parents(cls):
    #     parser = cls.get_parser()
    #
    #     parser.add_argument('--collapsed',
    #                         action='store_true',
    #                         default=False,
    #                         help='Collapsed LSVs thumbnails in the HTML summary.')
    #     return (
    #         cls.base_args(), cls.lsv_type_search_args(), cls.output_args(), parser
    #     )


def lsv_thumbnails(args):
    """
    Render html output for lsv thumbnails.
    :param args: command line arguments
    :return: None
    """
    output_dir = args.output
    output_html = Html.get_output_html(args)

    env = Html.get_env()
    template_file_name = args.type_analysis.replace("-", "_") + "_summary_template.html"
    sum_template = env.get_template(template_file_name)

    with open(path.join(output_dir, output_html), 'w') as voila_output:
        voila_output.write(
            sum_template.render(
                lsv_list=args.lsv_types,
                collapsed=args.collapsed
            )
        )

    Html.copy_static(args)
