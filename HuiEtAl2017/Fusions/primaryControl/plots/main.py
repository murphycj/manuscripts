import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.pyplot
from itertools import cycle
import pyensembl

from agfusion import AGFusionDB, Fusion

HORIZONTAL_LEVELS = [1,2,3,4]

class _Plot(object):
    def __init__(self, filename='', height=0, width=0, dpi=0, fontsize=12,
                 scale=0):

        self.filename = filename
        self.scale = scale
        self.width = width
        self.height = height
        self.dpi = dpi
        self.fontsize = fontsize

        self.fig = plt.figure(
            figsize=(self.width, self.height), dpi=self.dpi, frameon=False
        )
        self.ax = self.fig.add_subplot(111)
        self.rr = self.fig.canvas.get_renderer()

    def save(self):

        self.fig.savefig(
            self.filename,
            dpi=self.dpi,
            bbox_inches='tight'
        )

        plt.close(self.fig)
        plt.clf()

    def _scale(self, seq_length):
        """
        scale the sequence (protein or DNA)
        """

        if self.scale is None or self.scale < seq_length:
            self.normalize = seq_length
        else:
            self.normalize = self.scale

        self.offset = 0.05 + (1.0 - float(seq_length)/self.normalize)*0.45

        assert self.normalize >= seq_length, "length normalization should be >= protein length"


class _PlotProtein(_Plot):

    def __init__(self, transcript=None, colors=None,
                 rename=None, no_domain_labels=False,
                 exclude = [], *args, **kwargs):
        super(_PlotProtein, self).__init__(*args, **kwargs)
        self.transcript = transcript
        self.colors = colors
        self.rename = rename
        self.no_domain_labels = no_domain_labels
        self.exclude = exclude
        self.vertical_offset = 0.55

    def _draw_domains(self, domains):
        # plot domains

        domain_labels = {i:[] for i in HORIZONTAL_LEVELS}
        domain_labels_levels = {}
        domain_label_boxes = {i:[] for i in HORIZONTAL_LEVELS}
        lowest_level_plotted = HORIZONTAL_LEVELS[0]

        domains.sort(key=lambda x: x[3])

        domain_count = 0

        for domain in domains:

            # use domain name if available, otherwise use its ID

            if domain[1] is None:
                domain_name = str(domain[0])
            else:
                domain_name = str(domain[1])

            if domain_name in self.exclude:
                continue

            if domain_name in self.rename:
                domain_name = self.rename[domain_name]

            domain_start = (int(domain[3])/float(self.normalize))*0.9 + self.offset
            domain_end = (int(domain[4])/float(self.normalize))*0.9 + self.offset
            domain_center = (domain_end-domain_start)/2. + domain_start

            if not self.no_domain_labels:

                # for each newly plotted domain, loop through all previous
                # plotted domains and calculated the extent of overlap
                # then horizontally adjust the domain label to be plotted
                # closest to the protein domain structure or be
                # on the same level as the label it overlaps with the least

                #domain_stack_level = cycle()
                overlaps = {i:0.0 for i in HORIZONTAL_LEVELS}
                overlaps_all_levels = True
                min_overlap = [float("inf"),HORIZONTAL_LEVELS[0]]
                # plot domain at 1st level

                for level in HORIZONTAL_LEVELS:

                    level_pos = self.vertical_offset - 0.15 - (level-1.0)*0.1

                    tmp_domain_label = self.ax.text(
                        domain_center,
                        level_pos,
                        domain_name,
                        horizontalalignment='center',
                        verticalalignment='center',
                        fontsize=self.fontsize
                    )
                    tmp_domain_label_box = tmp_domain_label.get_window_extent(renderer=self.rr)

                    #check to see if it overlaps with anything

                    if len(domain_label_boxes[level])>0:
                        max_overlap = max([i.x1 - tmp_domain_label_box.x0 for i in domain_label_boxes[level]])

                        if max_overlap > 0.0:
                            overlaps[level] = max_overlap
                            tmp_domain_label.remove()

                            if max_overlap <= min_overlap[0]:
                                min_overlap = [max_overlap, level]
                        else:
                            domain_labels[level].append(tmp_domain_label)
                            domain_label_boxes[level].append(tmp_domain_label_box)
                            overlaps_all_levels = False

                            domain_labels_levels[domain_count] = level

                            if level > lowest_level_plotted:
                                lowest_level_plotted = level

                            break
                    else:
                        domain_labels[level].append(tmp_domain_label)
                        domain_label_boxes[level].append(tmp_domain_label_box)
                        overlaps_all_levels = False

                        domain_labels_levels[domain_count] = level

                        if level > lowest_level_plotted:
                            lowest_level_plotted = level

                        break

                # if the domain label overlaps with something on all levels
                # then plot it on the level where is overlaps the least

                if overlaps_all_levels:

                    level_pos = self.vertical_offset - 0.15 - (min_overlap[1]-1.0)*0.1

                    tmp_domain_label = self.ax.text(
                        domain_center,
                        level_pos,
                        domain_name,
                        horizontalalignment='center',
                        verticalalignment='center',
                        fontsize=self.fontsize
                    )
                    tmp_domain_label_box = tmp_domain_label.get_window_extent(renderer=self.rr)

                    domain_labels[min_overlap[1]].append(tmp_domain_label)
                    domain_label_boxes[min_overlap[1]].append(tmp_domain_label_box)

                    domain_labels_levels[domain_count] = level

                    if min_overlap[1] > lowest_level_plotted:
                        lowest_level_plotted = min_overlap[1]
            domain_count += 1

        # now we know how many levels of domains labels are needed, so
        # remove all levels, make the correction to self.vertical_offset
        # and replot all labels.

        for level, label in list(domain_labels.items()):
            for ll in label:
                ll.remove()

        self.levels_plotted = HORIZONTAL_LEVELS.index(lowest_level_plotted)
        self.vertical_offset += (0.05 * self.levels_plotted)

        domain_count = 0

        for domain in domains:

            if domain[1] is None:
                domain_name = str(domain[0])
            else:
                domain_name = str(domain[1])

            if domain_name in self.exclude:
                continue

            color = '#3385ff'
            if domain_name in self.colors:
                color = self.colors[domain_name]

            domain_start = (int(domain[3])/float(self.normalize))*0.9 + self.offset
            domain_end = (int(domain[4])/float(self.normalize))*0.9 + self.offset
            domain_center = (domain_end-domain_start)/2. + domain_start

            self.ax.add_patch(
                patches.Rectangle(
                    (
                        domain_start,
                        self.vertical_offset,
                    ),
                    domain_end - domain_start,
                    0.1,
                    color=color
                )
            )

            # fetch the level the domain label was determined it was to be
            # plotted on

            level = domain_labels_levels[domain_count]

            level_pos = self.vertical_offset - 0.15 - (level-1.0)*0.1

            tmp_domain_label = self.ax.text(
                domain_center,
                level_pos,
                domain_name,
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=self.fontsize
            )

            domain_count += 1


    def _draw_protein_length_markers(self, protein_length):
        # plot protein length markers

        self.line_end = protein_length/float(self.normalize)*0.9 + self.offset

        self.ax.text(
            0.5,
            self.vertical_offset - (0.5 + self.levels_plotted * 0.1),
            "Amino acid position",
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=self.fontsize
        )

        self.ax.add_line(plt.Line2D(
            (
                self.offset,
                self.offset+self.protein_frame_length
            ),
            (
                self.vertical_offset - (0.35 + self.levels_plotted * 0.05),
                self.vertical_offset - (0.35 + self.levels_plotted * 0.05)
            ),
            color='black'
            )
        )

        # left marker

        self.left_marker_line = self.ax.add_line(plt.Line2D(
            (
                self.offset,
                self.offset
            ),
            (
                self.vertical_offset - (0.38 + self.levels_plotted * 0.05),
                self.vertical_offset - (0.35 + self.levels_plotted * 0.05)
            ),
            color='black'
            )
        )
        self.left_marker_text = self.ax.text(
            self.offset,
            self.vertical_offset - (0.43 + self.levels_plotted * 0.05),
            "0",
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=self.fontsize
        )

        # draw markers for increments of 100 amino acids

        for i in range(1, protein_length+1):
            if (i % 100) == 0:
                self.left_marker_line = self.ax.add_line(plt.Line2D(
                    (
                        self.offset+(i/float(self.normalize)*0.9),
                        self.offset+(i/float(self.normalize)*0.9)
                    ),
                    (
                        self.vertical_offset - (0.38 + self.levels_plotted * 0.05),
                        self.vertical_offset - (0.35 + self.levels_plotted * 0.05)
                    ),
                    color='black'
                    )
                )

        # right marker

        self.right_marker_line = self.ax.add_line(plt.Line2D(
            (
                self.offset+self.protein_frame_length,
                self.offset+self.protein_frame_length
            ),
            (
                self.vertical_offset - (0.38 + self.levels_plotted * 0.05),
                self.vertical_offset - (0.35 + self.levels_plotted * 0.05)
            ),
            color='black'
            )
        )

        self.right_marker_text = self.ax.text(
            self.offset+self.protein_frame_length,
            self.vertical_offset - (0.43 + self.levels_plotted * 0.05),
            str(protein_length),
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=self.fontsize
        )

    def _draw_main_body(self):
        """
        main protein frame
        """

        self.ax.add_patch(
            patches.Rectangle(
                (self.offset, self.vertical_offset),
                self.protein_frame_length,
                0.1,
                fill=False
            )
        )

class PlotFusionProtein(_PlotProtein):
    def __init__(self, *args, **kwargs):
        super(PlotFusionProtein, self).__init__(*args, **kwargs)

    def _draw_junction(self):
        # add the junction

        self.ax.add_line(plt.Line2D(
            (
                (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset,
                (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset
            ),
            (
                self.vertical_offset - 0.05,
                self.vertical_offset + 0.15
            ),
            color='black'
            )
        )

        # middle marker, loop until it does not overlap with right marker

        overlaps = True
        line_offset = (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset
        text_offset = (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset
        junction_label_vertical_offset = 0.0

        right_marker_text_box = self.right_marker_text.get_window_extent(renderer=self.rr)
        left_marker_text_box = self.left_marker_text.get_window_extent(renderer=self.rr)

        while overlaps:

            # middle_marker_line_1/2/3 are to draw angled line

            middle_marker_line_1 = self.ax.add_line(plt.Line2D(
                (
                    (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset,
                    (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset
                ),
                (
                    self.vertical_offset - (0.37 + self.levels_plotted * 0.05) - junction_label_vertical_offset,
                    self.vertical_offset - (0.35 + self.levels_plotted * 0.05)
                ),
                color='black'
                )
            )

            middle_marker_line_2 = self.ax.add_line(plt.Line2D(
                (
                    line_offset,
                    (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset
                ),
                (
                    self.vertical_offset - (0.37 + self.levels_plotted * 0.05) - junction_label_vertical_offset,
                    self.vertical_offset - (0.37 + self.levels_plotted * 0.05) - junction_label_vertical_offset
                ),
                color='black'
                )
            )

            middle_marker_line_3 = self.ax.add_line(plt.Line2D(
                (
                    line_offset,
                    line_offset
                ),
                (
                    self.vertical_offset - (0.4 + self.levels_plotted * 0.05) - junction_label_vertical_offset,
                    self.vertical_offset - (0.37 + self.levels_plotted * 0.05) - junction_label_vertical_offset
                ),
                color='black'
                )
            )

            middle_marker_text = self.ax.text(
                text_offset,
                self.vertical_offset - (0.45 + self.levels_plotted * 0.05) - junction_label_vertical_offset,
                str(self.transcript.transcript_protein_junction_5prime),
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=self.fontsize
            )

            # detect if text overlaps

            middle_marker_text_box = middle_marker_text.get_window_extent(
                renderer=self.rr
            )

            # if overlaps then offset the junction text to the left

            if (right_marker_text_box.fully_overlaps(middle_marker_text_box)) and (left_marker_text_box.fully_overlaps(middle_marker_text_box)):
                junction_label_vertical_offset = junction_label_vertical_offset + 0.01

                middle_marker_line_1.remove()
                middle_marker_line_2.remove()
                middle_marker_line_3.remove()
                middle_marker_text.remove()

            elif right_marker_text_box.fully_overlaps(middle_marker_text_box):
                line_offset = line_offset - 0.01
                text_offset = text_offset - 0.01

                middle_marker_line_1.remove()
                middle_marker_line_2.remove()
                middle_marker_line_3.remove()
                middle_marker_text.remove()
            elif left_marker_text_box.fully_overlaps(middle_marker_text_box):
                line_offset = line_offset + 0.01
                text_offset = text_offset + 0.01

                middle_marker_line_1.remove()
                middle_marker_line_2.remove()
                middle_marker_line_3.remove()
                middle_marker_text.remove()
            else:
                overlaps = False

    def draw(self):
        self._scale(self.transcript.protein_length)
        self.protein_frame_length = self.transcript.protein_length/float(self.normalize)*0.9
        self._draw_domains(self.transcript.domains['fusion'])
        self._draw_protein_length_markers(self.transcript.protein_length)
        self._draw_junction()

        self._draw_main_body()

        self.ax.axis('off')
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)


peptides = [
    [60,24],
    [84,9],
    [103,11],
    [119,11],
    [171,9],
    [298,13],
    [445,20],
    [528,8],
    [620,11],
    [631,9],
    [649,11],
    [660,8],
    [668,10],
    [687,10],
    [740,16],
    [770,8],
    [2330,8],
    [2476,13]
]

agfusion_db = AGFusionDB('/Users/charlesmurphy/Desktop/Research/0914_hui/results/Fusions/plots/agfusion.mus_musculus.84.db', debug=False)
agfusion_db.build = 'mus_musculus' + '_' + str(84)
pyensembl_data = pyensembl.EnsemblRelease(84, 'mus_musculus')

fusion = Fusion(
    gene5prime=['ENSMUSG00000030849'],
    gene5primejunction=130167703,
    gene3prime=['ENSMUSG00000055322'],
    gene3primejunction=74016186,
    db=agfusion_db,
    pyensembl_data=pyensembl_data,
    protein_databases=['pfam','tmhmm'],
    noncanonical=False
)

pplot = PlotFusionProtein(
    filename='FGFR2-TNS1.png',
    width=10,
    height=3,
    dpi=90,
    scale=0,
    fontsize=10,
    colors={'TMhelix':'black','Pkinase_Tyr':'red'},
    rename=[],
    no_domain_labels=False,
    transcript=fusion.transcripts['ENSMUST00000122054-ENSMUST00000169786'],
    exclude=[]
)
pplot.draw()

for peptide in peptides:

    domain_start = (int(peptide[0])/float(pplot.normalize))*0.9 + pplot.offset
    domain_end = (int(peptide[0]+peptide[1])/float(pplot.normalize))*0.9 + pplot.offset
    domain_center = (domain_end-domain_start)/2. + domain_start

    pplot.ax.add_patch(
        patches.Rectangle(
            (
                domain_start,
                pplot.vertical_offset+0.2,
            ),
            domain_end - domain_start,
            0.05,
            color='grey'
        )
    )

pplot.save()
