from voila.utils.voilaLog import voilaLog

log = voilaLog()


class AsciiSpliceGraph(object):
    def __init__(self, gene, experiment):
        """
        Show splice graph for gene object
        :param gene: gene object
        """
        self.experiment = experiment
        self.gene = gene

    def make_display(self):
        """
        Make splice graph display
        :return:
        """
        display = []
        display_length = 1
        while not display:
            display_length += 1
            display = self.make_junctions_display(self.make_exons_display(display_length))

        return display

    def display(self):
        """
        Return ascii string for splice graph.
        :return: String
        """
        display = self.make_display()
        return '\n'.join([''.join(d) for d in display])

    def make_exons_display(self, display_length):
        """
        Generate display for exons
        :param display_length: how long display will be
        :return: list of strings
        """
        empty_char = ' '
        display = [empty_char for _ in range(int(display_length))]

        for exon in self.gene.exons:
            exon_start_index, exon_end_index = self.get_indexes((exon.start, exon.end), display_length)

            if display[exon_start_index] != empty_char:
                return []
            display[exon_start_index] = '('

            if display[exon_end_index] != empty_char:
                return []
            display[exon_end_index] = ')'

            for index in range(exon_start_index + 1, exon_end_index):
                display[index] = '='

            if exon.intron_retention:
                if display[exon_start_index + 1] != '=':
                    return []
                display[exon_start_index + 1] = 'i'

            if display[exon_end_index - 1] != '=':
                return []
            display[exon_end_index - 1] = str(exon.get('exon_type', self.experiment))

        self.validate(display)
        return display

    def get_indexes(self, coords, display_length):
        """
        Get indexes for elements.
        :param coords: coordinates
        :param display_length: how long the display is
        :return:
        """
        return self.get_index(coords[0], display_length), self.get_index(coords[1], display_length)

    def get_index(self, coord, display_length):
        """
        Get index for a coordinate.
        :param coord: coordinate
        :param display_length: length of display
        :return: int
        """
        gene_length = float(self.gene.end() - self.gene.start())
        start_percent_from_start = (coord - self.gene.start()) / gene_length
        return int(start_percent_from_start * (display_length - 1))

    def make_junctions_display(self, exons_display):
        """
        Make display for juctions
        :param exons_display: already generated exon display
        :return: list of lists
        """
        if not exons_display:
            return []

        display_length = len(exons_display)
        junctions_display = []

        for junction in self.gene.junctions:
            level_index = 0

            junc_start_index, junc_end_index = self.get_indexes((junction.start, junction.end), display_length)

            try:
                while [junctions_display[level_index][junc_start_index],
                       junctions_display[level_index][junc_end_index]] != [' ', ' ']:
                    level_index += 1
            except IndexError:
                junctions_display.append([' ' for _ in range(len(exons_display))])

            jd = junctions_display[level_index]

            for index in range(junc_start_index + 1, junc_end_index):
                jd[index] = '~'

            jd[junc_start_index] = '/'

            if jd[junc_end_index] != ' ':
                return []

            jd[junc_end_index] = '\\'

            junction_type = list('t' + str(junction.get('junction_type', self.experiment)))
            junction_type.reverse()
            for ti, t in enumerate(junction_type):
                if jd[junc_end_index - 1 - ti] != '~':
                    return []
                jd[junc_end_index - 1 - ti] = t

            for ri, r in enumerate('r' + str(junction.get('reads', self.experiment))):
                if jd[junc_start_index + 1 + ri] != '~':
                    return []
                jd[junc_start_index + 1 + ri] = r

        junctions_display.reverse()
        junctions_display.append(exons_display)
        junctions_display.append(['-' for _ in range(display_length)])
        return junctions_display

    def validate(self, exons_display):
        """
        Validate exon display
        :param exons_display: exon display
        :return: None
        """
        start = True
        for index, item in enumerate(exons_display):
            if (item == '(' and not start) or (item == ')' and start):
                raise Exception("NOPE! Not a valid splicegraph. At index: {0}".format(index))
            if item == '(' and start:
                start = False
            if item == ')' and not start:
                start = True
