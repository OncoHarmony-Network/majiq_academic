import os
import pdb
from voila.tools import Tool
from voila.tools.utils import io_caleb
from voila.utils.voila_log import voila_log



# Created by Matthew Gazzara, adapted by Caleb
# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'mgazzara'

LOG = voila_log()


class ThisisLsvVenns(Tool):
    help = ''

    def arguments(self):
        parser = self.get_parser()
        parser.add_argument('voila_files',
                            type=str,
                            help='Path to a file with line-by-line listing of voila text files.')
        parser.add_argument('set_names',
                            type=str,
                            help='List of 2 or 3 set names separated by commas and no spaces')
        help_mes = "dPSI threshold to use for calling sig LSVs)"
        parser.add_argument('dpsi_thresh',
                            type=float,
                            help=help_mes)
        help_mes = "P(dPSI) threshold to use for calling sig LSVs)"
        parser.add_argument('prob_thresh',
                            type=float,
                            help=help_mes)
        help_mes = "Outfile path and name. Include .pdf in the file path name if you want vector graphics."
        parser.add_argument('outfile',
                            type=str,
                            help=help_mes)
        help_mes = "Flag: only consider LSVs that were quantifiable in all comparisons."
        parser.add_argument('--remove_non_shared',
                            action='store_true',
                            default=False,
                            help=help_mes)
        help_mes = "Flag: Write to file the shared LSVs?"
        parser.add_argument('--write_shared',
                            action='store_true',
                            default=False,
                            help=help_mes)
        help_mes = "Flag: Turn off plot text so you can add your own?"
        parser.add_argument('--clean_plot',
                            action='store_true',
                            default=False,
                            help=help_mes)
        help_mes = "Flag: Consider intron retention events?"
        parser.add_argument('--consider_ir_events',
                            action='store_true',
                            default=False,
                            help=help_mes)
        help_mes = "Optional. Additional pylab.savefig() arguments."
        parser.add_argument('--savefig_args',
                            nargs='*',
                            help=help_mes)
        return parser

    def run(self, args):
        is_valid, the_res = io_caleb.is_likely_list_of_txtfiles(args.voila_files)
        if not is_valid:
            raise ValueError("%s Wasn't a path to a valid list of voila files.." % args.voila_files)
        set_names = args.set_names.split(",")
        if len(set_names) != 2 and len(set_names) != 3:
            raise RuntimeError("%s is not a valid list of set names. Please provide 2 or 3 sep by commas, no spaces."
                               % set_names)
        if args.write_shared:
            out_dir = os.path.dirname(args.outfile)
            outfile = os.path.basename(args.outfile)[0:-4] + "_OverlapOutput.txt"
        else:
            out_dir = "./"
            outfile = 'OverlapOutput'
        the_plot = make_venn(the_res,
                             set_names,
                             thresh=args.dpsi_thresh,
                             prob_thresh=args.prob_thresh,
                             remove_non_shared=args.remove_non_shared,
                             write_shared=args.write_shared,
                             out_dir=out_dir,
                             out_file=outfile,
                             no_txt=args.clean_plot,
                             consider_ir_events=args.consider_ir_events,
                             interactive_plotting=False)  # for non-gui systems, need to turn off interactive plotting.
        if args.savefig_args:
            # TODO this..
            raise RuntimeError("I haven't implemented this yet... - Caleb")
            #the_plot.savefig(args.outfile, args.savefig_args)
        else:
            the_plot.savefig(args.outfile)


def make_venn(voila_list,
              set_names,
              thresh=0.2,
              prob_thresh=0,
              no_txt=False,
              write_shared=False,
              out_file='OverlapOutput',
              remove_non_shared=False,
              title_prefix='',
              out_dir='./',
              compare_genes=False,
              interactive_plotting=False,
              consider_ir_events=False):
    """
    voila_list should be a list of strings that are the paths to the 2 or 3 dPSI voila text files you wish to overlap
    set_names should be a list of strings used as names in the Venn diagram
    thresh: minimum E(dPSI) threshold. An LSV should have 1 or more junctions >= thresh to be considered changing
    prob_thresh: minimum value of Prob(|dPSI| >= v) value, where v is the E(dPSI) threshold. LSV should have 1 or more
    junctions meeting this to be considered changing

    Other Options:
        no_txt: do not display the numbers nor the labels in the Venn diagram [default: False]
        compare_genes: Do overlaps / counts at gene level and not LSV level [default: False]
        remove_non_shared: Remove from consideration LSVs that were not quantified in all dPSI files given to make the
        comparison more fair
        consider_ir_events: should IR LSVs be ignored? (IR is ignored if its dPSI is >= 0.05 in a given LSV)
    """
    if not interactive_plotting:
        # Must do this before importing other matplotlib stuff
        # This turns of the display capabilities of matplotlib, which is necessary when working on a dumb linux machine
        import matplotlib
        matplotlib.use('Agg')
    import pylab as pyl
    from matplotlib_venn import venn2, venn3
    if not len(voila_list) == len(set_names):
        raise RuntimeError('number of voila paths input should equal number of set names...')
    if len(voila_list) > 3: raise RuntimeError \
        ('too many dPSI sets to compare...can only handle 2 or 3 for venn diagram')
    if len(voila_list) < 2: raise RuntimeError('too few dPSI sets to compare...can only handle 2 or 3 for venn diagram')
    nSets = len(voila_list)
    set_list = []
    all_sets = []
    all_irs = list()
    for n in range(len(set_names)):
        v_file = voila_list[n]
        evs = get_events(v_file,
                         thresh=thresh,
                         prob_thresh=prob_thresh,
                         only_genes=compare_genes,
                         ignore_ir=not consider_ir_events)
        if not consider_ir_events:
            evs, these_ir = evs
            all_irs.extend(these_ir)
        set_list.append(evs)
        if remove_non_shared == True:
            allEvents = get_events(v_file, thresh=0, prob_thresh=0, only_genes=compare_genes)
            all_sets.append(allEvents)
        LOG.info('Expt %s had %s events at E(dPSI) thresh of %s and prob thresh of %s' %
              (set_names[n], len(evs), thresh, prob_thresh))
    if not consider_ir_events:
        all_irs = set(all_irs)
        set_list = [thisset - all_irs for thisset in set_list]
    pyl.figure(0)
    if remove_non_shared == False:
        if nSets == 2:
            venn = venn2(subsets=(set_list[0], set_list[1]), set_labels=set_names)
            if no_txt:
                for vv in ['10', '11', '01']:
                    venn.get_label_by_id(vv).set_text('')
                for text in venn.set_labels:
                    text.set_text('')
        else:
            venn = venn3(subsets=(set_list[0], set_list[1], set_list[2]), set_labels=set_names)
            if no_txt:
                for vv in ['100', '010', '001', '110', '101', '011', '111']:
                    venn.get_label_by_id(vv).set_text('')
                for text in venn.set_labels:
                    text.set_text('')
        pyl.title('%s overlaps:\n|dPSI|>=%s at prob %s' % (title_prefix, thresh, prob_thresh))
    else:
        if nSets == 2:
            all_evs = all_sets[0] & all_sets[1]
            venn = venn2(subsets=(set_list[0] & all_evs, set_list[1] & all_evs), set_labels=set_names)
            if no_txt == True:
                for vv in ['10', '11', '01']:
                    venn.get_label_by_id(vv).set_text('')
                for text in venn.set_labels:
                    text.set_text('')
        else:
            all_evs = all_sets[0] & all_sets[1] & all_sets[2]
            venn = venn3(subsets=(set_list[0] & all_evs, set_list[1] & all_evs, set_list[2] & all_evs),
                         set_labels=set_names, set_colors=['#3182bd', '#de2d26', '#31a354'], alpha=0.6)
            if no_txt == True:
                for vv in ['100', '010', '001', '110', '101', '011', '111']:
                    venn.get_label_by_id(vv).set_text('')
                for text in venn.set_labels:
                    text.set_text('')
        pyl.title('%s overlaps:\n|dPSI|>=%s at prob %s' % (title_prefix, thresh, prob_thresh))

    if nSets == 2:
        all_changes = set_list[0] | set_list[1]
        shared2 = set_list[0] & set_list[1]
    else:
        all_changes = set_list[0] | set_list[1] | set_list[2]
        shared2 = (set_list[0] & set_list[1]) | (set_list[0] & set_list[2]) | (set_list[1] & set_list[2])
    LOG.info("Total number of significantly changing LSVs: %s" % len(all_changes))
    LOG.info("Total number of shared changing LSVs: %s" % len(shared2))

    if write_shared == True:
        if compare_genes == True:
            fw = open(out_dir + '/' + out_file + 'sharedGenes.txt', 'w')
            fw.write('#GeneID')
            for n in range(len(set_names)):
                fw.write('\t%s_hit' % set_names[n])
            fw.write('\n')
            for e in shared2:
                fw.write(e)
                for n in range(len(set_list)):
                    chg_set = set_list[n]
                    exp_set = all_sets[n]
                    if e in chg_set:
                        fw.write('\t1')
                    else:
                        if e in exp_set:
                            fw.write('\t0')
                        else:
                            fw.write('\tNaN')
                fw.write('\n')
            fw.close()
        else:
            dic = {}
            for n in range(len(set_names)):
                expt = set_names[n]
                fd = open(voila_list[n], 'r')
                for line in fd:
                    if line.startswith('#'):
                        continue
                    l = line.strip().split('\t')
                    ID = l[2]
                    if not ID in all_changes: continue
                    try:
                        dic[ID]
                    except:
                        dic[ID] = {}
                        dic[ID]['info'] = {}
                        dic[ID]['info']['name'] = l[0]
                        dic[ID]['info']['rest'] = l[7:]
                    dic[ID][expt] = {}
                    if ID in set_list[n]:
                        dic[ID][expt]['hit'] = 1
                    else:
                        dic[ID][expt]['hit'] = 0
                    dic[ID][expt]['dPSIs'] = l[3]
                    dic[ID][expt]['psi1'] = l[5]
                    dic[ID][expt]['psi2'] = l[6]
                fd.close()
            fw = open(out_dir + '/' + out_file + ".txt", 'w')
            fw.write('#Name\tLSV')
            for n in range(len(set_names)):
                fw.write('\t%s_hit' % set_names[n])
            for n in range(len(set_names)):
                fw.write('\t%s_dPSI' % set_names[n])
            fw.write('\tLSVtype\tA5SS\tA3SS\tES\tNumJunc\tNumExon\tDeNovoJunc\tchr\tstrand\t'
                     'JuncCoord\tExonCoord\tExonAltStart\tExonAltEnd\tIRcoords\tLink\n')
            for e in shared2:
                fw.write('%s\t%s' % (dic[e]['info']['name'], e))
                for exp in set_names:
                    try:
                        fw.write('\t%s' % (dic[e][exp]['hit']))
                    except:
                        fw.write('\tNaN')
                for exp in set_names:
                    try:
                        fw.write('\t%s' % (dic[e][exp]['dPSIs']))
                    except:
                        fw.write('\tNaN')
                for ii in dic[e]['info']['rest']:
                    fw.write('\t%s' % ii)
                fw.write('\n')
            fw.close()
    return pyl


def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb


def get_events(voila_all_path,
               thresh=0.2,
               prob_thresh=0,
               write_names=False,
               fname='ChangingLSVs',
               only_genes=False,
               ignore_ir=False):

    events = list()
    ir_event = list()
    LOG.info(voila_all_path)
    fd = open(voila_all_path, 'r')
    if only_genes == True:
        idx = 1
    else:
        idx = 2
    for line in fd:
        l = line.strip().split('\t')
        if line.startswith('#'):
            lsv_type_index = l.index("LSV Type")
            if prob_thresh > 0.0:
                voila_prob = float(l[4].split('=')[-1].split(')')[0])
                if not voila_prob == thresh:
                    raise RuntimeError \
                        ('WARNING! For dPSI file %s the voila probability threshold was run for %s ... '
                         'not your threshold of %s' %
                         (voila_all_path, voila_prob, thresh))
            continue
        dPSIs = [abs(float(x)) for x in l[3].split(';')]
        probs = [float(x) for x in l[4].split(';')]
        if max(dPSIs) >= thresh:
            if max(probs) >= prob_thresh:
                if ignore_ir:
                    is_intron_ret_ev = True if l[lsv_type_index][-1] == "i" else False
                    if is_intron_ret_ev:
                        # My threshold for an intron changing event being ignore-worthy
                        if dPSIs[-1] >= 0.05:
                            ir_event.append(l[idx])
                            continue
                events.append(l[idx])
    fd.close()
    if write_names == True:
        fw = open('%s_dPSI%s_Prob%s.txt' % (fname, str(int(thresh * 100)), str(int(prob_thresh * 100))), 'w')
        evs = set(events)
        for e in evs: fw.write('%s\n' % e)
        fw.close()
    if ignore_ir:
        return set(events), set(ir_event)
    return set(events)


def write_all_changing_data(voila_list, set_names, out_file, thresh=0.2, prob_thresh=0):
    all_changes = set()
    set_list = [[] for x in range(len(voila_list))]
    for n in range(len(voila_list)):
        voila_file = voila_list[n]
        chg_evs = get_events(voila_file, thresh=thresh, prob_thresh=prob_thresh)
        for ev in chg_evs:
            all_changes.add(ev)
            set_list[n].append(ev)
    dic = {}
    for ss in set_list:
        LOG.info(len(ss))
    for n in range(len(set_names)):
        expt = set_names[n]
        LOG.info('###%s' % expt)
        fd = open(voila_list[n], 'r')
        for line in fd:
            if line.startswith('#'):
                continue
            l = line.strip().split('\t')
            ID = l[2]
            if not ID in all_changes:
                continue
            try:
                dic[ID]
            except:
                dic[ID] = {}
                dic[ID]['info'] = {}
                dic[ID]['info']['name'] = l[0]
                dic[ID]['info']['rest'] = l[7:]
            dic[ID][expt] = {}
            if ID in set_list[n]:
                dic[ID][expt]['hit'] = 1
            else:
                dic[ID][expt]['hit'] = 0
            dic[ID][expt]['dPSIs'] = l[3]
            dic[ID][expt]['psi1'] = l[5]
            dic[ID][expt]['psi2'] = l[6]
        fd.close()
    fw = open(out_file, 'w')
    fw.write('#Name\tLSV')
    for n in range(len(set_names)):
        fw.write('\t%s_hit' % set_names[n])
    for n in range(len(set_names)):
        fw.write('\t%s_dPSI' % set_names[n])
    fw.write('\tLSVtype\tA5SS\tA3SS\tES\tNumJunc\tNumExon\tDeNovoJunc\t'
             'chr\tstrand\tJuncCoord\tExonCoord\tExonAltStart\tExonAltEnd\tIRcoords\tLink\n')
    for e in all_changes:
        fw.write('%s\t%s' % (dic[e]['info']['name'], e))
        for exp in set_names:
            try:
                fw.write('\t%s' % (dic[e][exp]['hit']))
            except:
                fw.write('\tNaN')
        for exp in set_names:
            try:
                fw.write('\t%s' % (dic[e][exp]['dPSIs']))
            except:
                fw.write('\tNaN')
        for ii in dic[e]['info']['rest']:
            fw.write('\t%s' % ii)
        fw.write('\n')
    fw.close()
