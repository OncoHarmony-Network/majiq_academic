import argparse
from .tools import *
import importlib




def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('l_res_files', nargs='+')
    parser.add_argument('l_tools', nargs='+')
    parser.add_argument('-t', '--threshold', default=0.2, type=float)
    parser.add_argument('--use-prob', default=False, action="store_true",
                        help="For MAJIQ, calculate sum_v P(deltaPSI > V)")

    args = parser.parse_args()

    results = []
    for idx, fl in enumerate(args.l_tools):
        importlib.import_module (".parser.%s.parse" % tools_dict[fl])
        results = parse(args.l_res_files[idx])

    for method_name, ranks_pair in ranks.items():
        print "Ranking %s...." % method_name
        rank1, rank2 = ranks_pair
        print "Num events", len(rank1), len(rank2)
        print "Calculating the ratios..."
        # calculate the ratios
        ratios = []
        events = []

        max_events = min(args.max, len(rank1))

        fdr = []
            import sys

            for i in xrange(max_events):
                chunk1 = list(rank1[0:i + 1])
                chunk2 = list(rank2[0:i + 1])
                # check if event1 is into chunk2
                found = 0
                for event1 in chunk1:
                    found += _is_in_chunk(event1, chunk2)
                    if i == max_events - 1:
                        events.append([event1, _is_in_chunk(event1, chunk2)])
                ratios.append(float(found) / max_events)
                if i % 20 == 0:
                    print "%s..." % i,
                    sys.stdout.flush()

            ratios = array(ratios)

        print "RESULT:", ratios[0:10], "...", ratios[-10:], "length", ratios.shape
        print "Saving in %s" % args.output

        if not os.path.exists(args.output):
            os.makedirs(args.output)

        pickle.dump(ratios,
                    open(args.output + "/ratios.%s.%s.pickle" % (str(args.type_rank).replace('-', '_'), method_name),
                         'w'))

        print "Saving events... in %s " % args.output
        pickle.dump(events,
                    open(args.output + "/events.%s.%s.pickle" % (str(args.type_rank).replace('-', '_'), method_name),
                         'w'))

        print "Saving N1 size... in %s " % args.output
        pickle.dump(n1[method_name],
                    open(args.output + "/n1.%s.%s.pickle" % (str(args.type_rank).replace('-', '_'), method_name), 'w'))

        if args.fdr:
            # print "FDR:", fdr[0:10], "...", fdr[-10:], "length", fdr.shape
            pickle.dump(fdr,
                        open("%s/fdr.%s.%s.pickle" % (args.output, method_name, str(args.type_rank).replace('-', '_')),
                             'w'))
            pickle.dump(v_values, open(
                "%s/fdr.%s.%s_v.pickle" % (args.output, method_name, str(args.type_rank).replace('-', '_')), 'w'))
            plot_fdr(args.output, method_name, fdr)

        if "old_majiq" in method_name:
            only_exp1_ranks.append(ratios)

    if args.create_restrict_plot:
        create_restrict_plot(ranks)

    print "Done!"
