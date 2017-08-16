import pdb

from voila.tools import Tool
from voila.tools.utils import io_caleb


# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'


class ThisisLookup(Tool):
    help = 'Given a directory and Gene Name, or Gene ID, or LSV ID, prettily-print all LSVs'

    def arguments(self):
        parser = self.get_parser()
        parser.add_argument('directory',
                            type=str,
                            help='Directory or file list where voila texts are listed.')
        help_mes = "dPSI threshold by which to call junctions as changing"
        parser.add_argument('--dpsi_thresh',
                            '--dpsi',
                            type=float,
                            help=help_mes,
                            default=0.2)
        help_mes = "Prob(dPSI) threshold by which to call junctions as changing"
        parser.add_argument('--prob_dpsi_thresh',
                            '--prob',
                            type=float,
                            help=help_mes,
                            default=0.95)
        help_mes = 'Optional pattern matching to identify the voila text files'
        parser.add_argument('-p',
                            '--pattern',
                            default="*tsv",
                            type=str,
                            help=help_mes)
        help_mes = 'Which comparisons or samples to lookup ID in? Single space or comma separated please.'
        parser.add_argument('--names',
                            type=str,
                            help=help_mes)
        help_mes = "Flag: don't import IR LSVs"
        parser.add_argument('--no_ir',
                            action='store_true',
                            help=help_mes,
                            default=True)
        return parser

    def run(self, args):
        # parse the comparisons argument
        if args.names:
            if "," in args.names or " " in args.names:
                args.names.replace(" ", ",")
                to_lookup = args.names.split(",")
            else:
                to_lookup = [args.names]
            dont_remove_dups = False
        else:
            to_lookup = None
            dont_remove_dups=True
        imported = io_caleb.quick_import(input=args.directory,
                                         cutoff_d_psi=args.dpsi_thresh,
                                         cutoff_prob=args.prob_dpsi_thresh,
                                         pattern=args.pattern,
                                         keep_ir=args.no_ir,
                                         comparisons=to_lookup)
        io_caleb.check_is_ignant(imported, args.dpsi_thresh)
        pdb.set_trace()
