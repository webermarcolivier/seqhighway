import itertools
import matplotlib.colors as mcolors
import seaborn as sns
from copy import copy
import re
from Bio.SeqFeature import SeqFeature, FeatureLocation
from .config import config_default as CONFIG



VERBOSE = 2

restriction_dict = {
    # EcoRI
    # 5'---G     AATTC---3'
    # 3'---CTTAA     G---5'
    'EcoRI':{'motif':'GAATTC', 'cut_position':(+1, -1)},
    # BamHI
    # 5'---G     GATCC---3'
    # 3'---CCTAG     G---5'
    'BamHI':{'motif':'GGATCC', 'cut_position':(+1, -1)},
    # BglII
    # 5'---A     GATCT---3'
    # 3'---TCTAG     A---5'
    'BamHI':{'motif':'AGATCT', 'cut_position':(+1, -1)},
    # EcoRV
    # 5'---GAT  ATC---3'
    # 3'---CTA  TAG---5'
    'EcoRV':{'motif':'GATATC', 'cut_position':(+3, -3)},
    # HindIII
    # 5'---A     AGCTT---3'
    # 3'---TTCGA     A---5'
    'HindIII':{'motif':'AAGCTT', 'cut_position':(+1, -1)},
    # BglII
    # 5'---A     GATCT---3'
    # 3'---TCTAG     A---5'
    'BglII':{'motif':'AGATCT', 'cut_position':(+1, -1)},
    # PstI
    # 5'---CTGCA  G---3'
    # 3'---G  ACGTC---5'
    'PstI':{'motif':'CTGCAG', 'cut_position':(+1, -1)},
    # SalI
    # 5'---G  TCGAC---3'
    # 3'---CAGCT  G---5'
    'SalI':{'motif':'GTCGAC', 'cut_position':(+1, -1)},
}

def get_color_map(x="Blues"):
    color_list = None
    if x is None:
        return None
    elif type(x) is str:
        name = x
        color_list = sns.color_palette(name, 61)
    elif len(x[0]) > 1:
        color_list = x
    else:
        raise ValueError("argument x should be either the name of a seaborn color palette, or a list of rgb colors.")
    colormap = sns.blend_palette(color_list, as_cmap=True, input='rgb')
    return colormap


class Graph:
    """Adapted from dnafeaturesviewer, see https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer
    Author: Zulko

    Minimal implementation of non-directional graphs.
    Parameters
    ----------
    nodes
      A list of objects. They must be hashable
    edges
      A list of the form [(n1,n2), (n3,n4)...] where (n1, n2) represents
      an edge between nodes n1 and n2
    """
    def __init__(self, nodes, edges):
        self.nodes = nodes
        self.neighbors = {n: [] for n in nodes}
        for n1, n2 in edges:
            self.neighbors[n1].append(n2)
            self.neighbors[n2].append(n1)


def compute_features_levels(features, base_level=0):
    """Adapted from dnafeaturesviewer, see https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer
    Author: Zulko

    Compute the vertical levels on which the features should be displayed
    in order to avoid collisions.
    `features` must be a list of `dna_features_viewer.GraphicFeature`.
    The method used is basically a graph coloring:
    - The nodes of the graph are features and they will be colored with a level
    - Two nodes are neighbors iff their features's locations overlap
    - Levels are attributed to nodes iteratively starting with the nodes
      corresponding to the largest features.
    - A node receives the lowest level (starting at 0) that is not already
      the level of one of its neighbors.
    """
    edges = [
        (f1, f2)
        for f1, f2 in itertools.combinations(features, 2)
        if f1.overlaps_with(f2)
    ]
    graph = Graph(features, edges)
    levels = {n: None for n in graph.nodes}

    def collision(base_level, node):
        """Return True if the node placed at base_level collides with
        its neighbors in the graph."""
        for neighbor in graph.neighbors[node]:
            level = levels[neighbor]
            if level is None:
                continue
            # nlines ????
#             if 'nlines' in neighbor.data:
#                 top = numpy.ceil(level + 0.5 * neighbor.data['nlines'])
#                 if level <= base_level < top:
#                     return True
                
#                 top = numpy.ceil(base_level + 0.5 * node.data['nlines'])
#                 if base_level <= level < top:
#                     return True
#             else:
            if level == base_level:
                return True
        return False

    for node in sorted(graph.nodes, key=lambda f: -f.length):
        while collision(base_level, node):
            base_level += 1
        levels[node] = base_level
    return levels



def extract_feature_name(feature, qualifier_priority_list=['id', 'name', 'locus_tag']):
    quals = feature.qualifiers
    name = None
    for qual in qualifier_priority_list[::-1]:
        if qual in quals.keys():
            q = quals[qual]
            if type(q) is list: q = q[0]
            name = q
    return name


def convert_strand_numeric_to_strand_string(strand_num):
    # Convert strand +1/-1 to +/-
    if strand_num == +1:
        strand = '+'
    elif strand_num == -1:
        strand = '-'
    elif strand_num == '+' or strand_num == '-':
        strand = strand_num
    else:
        strand = None
    return strand


class TrackFeature():
    
    newid = itertools.count()
    
    def __init__(self, start, end, strand, xmin, xmax, feature_bio=None, name=None,
                 level=0, zorder=0, annot_class=None, track=None):
        self.start = start
        self.end = end
        self.start_uncropped = start
        self.end_uncropped = end
        self.strand = strand
        self.xmin = xmin
        self.xmax = xmax
        self.feature_bio = feature_bio
        self.length = None
        self.track = track

        self.visible = True
        self.cropped_left = False
        self.cropped_right = False

        if self.start is not None and self.end is not None:
            if self.xmin > self.start:
                self.cropped_left = True
                self.start = self.xmin
            if self.xmax < self.end:
                self.cropped_right = True
                self.end = self.xmax
            if self.xmin > self.end - 1 or self.xmax - 1 < self.start:
                self.visible = False
            
            if self.visible:
                self.length = self.end - self.start
        
        self.annot_class = annot_class
        self.name = name
        self.id_ = next(self.newid)
        self.level = level
        self.zorder = zorder
        self.format_attr = CONFIG.get(annot_class)
        if self.format_attr is None:
            self.format_attr = dict()
            
    @property
    def jstart(self):
        return self.start - self.xmin

    @property
    def jend(self):
        return self.end - self.xmin

    def __len__(self):
        return len(self.end - self.start)
        
    # from dnafeaturesviewer
    def overlaps_with(self, other):
        """Return True iff the feature's location overlaps with feature `other`
        """
        loc1, loc2 = (self.start, self.end), (other.start, other.end)
        loc1, loc2 = sorted(loc1), sorted(loc2)
        loc1, loc2 = sorted([loc1, loc2], key=lambda loc: loc[0])
        return loc1[1] > loc2[0]
    
    def get_css_attr(self):
        css_att = dict()
        
        # non-css attributes first
        for k, v in self.format_attr.items():
            if k == 'backgroundcolor':
                # Convert matplotlib colors
                css_att['background-color'] = mcolors.to_hex(v, keep_alpha=True)
                if VERBOSE >= 2: print("css_att['background-color']", css_att['background-color'])

            if k == 'fontcolor':
                # Convert matplotlib colors
                css_att['color'] = mcolors.to_hex(v, keep_alpha=True)

        # Compute color of the feature
        bg_score_name = self.format_attr.get('background_score_name')
        bg_score_cmap = get_color_map(self.format_attr.get('background_score_cmap'))
        vmin = self.format_attr.get('background_score_vmin')
        vmax = self.format_attr.get('background_score_vmax')
        if bg_score_name is not None and bg_score_cmap is not None:
            score = self.feature_bio.qualifiers.get(bg_score_name)
            if score is None:
                score = 1
            norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
            css_att['background-color'] = mcolors.to_hex(mcolors.to_rgb(bg_score_cmap(norm(score))), keep_alpha=True)

        # The css attribute will override the more general attributes defined above
        if 'css_attributes' in self.format_attr.keys():
            css_att.update(self.format_attr['css_attributes'])

        return css_att
    
    def get_css(self, track_height=None, css_units=None, colwidth=None, base_fontsize=None):
        indent = '    '
        css_class = indent + '#featId{:d}.SH_{} {{\n'.format(self.id_, self.annot_class)
        css_class += indent + indent + 'width: 100%;\n'
        css_class += indent + indent + 'height: 100%;\n'
        for att, val in self.get_css_attr().items():
            css_class += indent + indent + '{}: {};\n'.format(att, val)
        css_class += indent + '}\n'
        return css_class
    
    def get_css_tag(self, j=None):
        tag_open = '<div id="featId{:d}" class="SH_{}">'.format(self.id_, self.annot_class)
        tag_close = '</div>'
        class_list = ['SH_{}'.format(self.annot_class)]
        return tag_open, tag_close, class_list, self.id_

    def get_css_class_and_id(self, j=None):
        return 'id="featId{:d}" class="SH_{}"'.format(self.id_, self.annot_class)
    
    def get_text(self, text='', j=None):
        return text


class FeatureOverlay(TrackFeature):
    
    def __init__(self, start, end, strand, xmin, xmax, **kwargs):
        super().__init__(start=start, end=end, strand=strand, xmin=xmin, xmax=xmax, **kwargs)
        self.modifies_text = False
            
    def get_text(self, text, j=None):
        return text
        
        
class FeatureFloating(TrackFeature):
    
    def __init__(self, start, end, strand, xmin, xmax, **kwargs):
        super().__init__(start=start, end=end, strand=strand, xmin=xmin, xmax=xmax, **kwargs)
        self.modifies_text = True


class FeatureFixed(TrackFeature):
    """The FeatureFixed is fixed on the track."""
    
    def __init__(self, start, end, strand, xmin, xmax, level, **kwargs):
        super().__init__(start=start, end=end, strand=strand, xmin=xmin, xmax=xmax, level=level, **kwargs)


class FeatureArrow(FeatureFloating):
    
    def __init__(self, start, end, strand, xmin, xmax, **kwargs):
        super().__init__(start=start, end=end, strand=strand, xmin=xmin, xmax=xmax, **kwargs)
        self.show_label = True
        self.label_rel_pos = 1    # position of the annotation label with respect to the strand-specific start
     
    # Simple reprensentation based on text arrow ------->
#     def get_text(self, text, j):
#         if self.strand == '+':
#             if self.jstart <= j < self.jend - 1:
#                 return '-'
#             else:
#                 return '>'
#         elif self.strand == '-':
#             if self.jstart < j <= self.jend - 1:
#                 return '-'
#             else:
#                 return '<'

    def get_text(self, text='', j=None):
        new_text = text
#         new_text = u'⭪'
#         new_text = u'⭬'
        # show the label, and also an truncated ellipsis showing that the feature
        # has been cropped
        if self.show_label:
            name = self.name if self.name is not None else ''
            if self.strand == '+':
                if j == self.jstart + self.label_rel_pos:
                    new_text = name
                elif j == self.jstart and self.cropped_left:
                    new_text = u'⋯'
                elif j == self.jend - 2 and self.cropped_right:
                    new_text = u'⋯'
            elif self.strand == '-':
                if j == self.jend - 1 - self.label_rel_pos:
                    new_text = name
                elif j == self.jstart + 1 and self.cropped_left:
                    new_text = u'⋯'
                elif j == self.jend - 1 and self.cropped_right:
                    new_text = u'⋯'
        return new_text
            
    def get_css(self, track_height, css_units, colwidth, base_fontsize):
        css_class = super().get_css(track_height, css_units, colwidth, base_fontsize)
        css_class += """    #featId{:d}.SH_triangle-left {{
        width: 0; 
        height: 0; 
        border-top: {}{} solid transparent;
        border-bottom: {}{} solid transparent;   
        border-right: {}{} solid {};
        background: none;
    }}
    #featId{:d}.SH_triangle-right {{
        width: 0; 
        height: 0; 
        border-top: {}{} solid transparent;
        border-bottom: {}{} solid transparent;   
        border-left: {}{} solid {};
        background: none;
    }}\n""".format(self.id_,
                   0.5*track_height, css_units,
                   0.5*track_height, css_units,
                   # we use a small overshoot to improve pixel rounding
                   1.1*colwidth, css_units,
                   self.get_css_attr()['background-color'],
                   self.id_,
                   0.5*track_height, css_units,
                   0.5*track_height, css_units,
                   1.1*colwidth, css_units,
                   self.get_css_attr()['background-color'])
        if self.show_label:
            css_class += """    #featId{:d}.SH_label {{
        position: relative;
        z-index: 10;
        overflow-x: visible;
        float: {};
    }}\n""".format(self.id_, "left" if self.strand == '+' else "right")
        return css_class
    
    def get_css_tag(self, j):
        class_list = []
        shape = 'normal'
        if self.strand == '+':
            if self.jstart <= j < self.jend - 1:
                shape = 'normal'
            else:
                shape = 'arrow-right'
        elif self.strand == '-':
            if self.jstart < j <= self.jend - 1:
                shape = 'normal'
            else:
                shape = 'arrow-left'
        if shape == 'normal':
            class_list = [self.annot_class]
        elif shape == 'arrow-left':
            class_list = [self.annot_class, "triangle-left"]
        elif shape == 'arrow-right':
            class_list = [self.annot_class, "triangle-right"]
                    
        class_list = ['SH_' + c for c in class_list]
        tag_open = '<div id="featId{:d}" class="{}">'.format(self.id_, ' '.join(class_list))
        tag_close = '</div>'
        
        # We have to add one more div level, in order to be able to float the label to the left or right
        # in the inner div, and keep the outer div to the 100% width with its background.
        if self.show_label:
            if ((self.strand == '+' and j == self.jstart + self.label_rel_pos) or
                    (self.strand == '-' and j == self.jend - 1 - self.label_rel_pos)):
                tag_open = tag_open + '<div id="featId{:d}" class="SH_label">'.format(self.id_)
                tag_close = '</div>' + tag_close

        is_cropped_pos = ((self.strand == '+' and j == self.jstart and self.cropped_left) or
                          (self.strand == '+' and j == self.jend - 2 and self.cropped_right) or
                          (self.strand == '-' and j == self.jend - 1 and self.cropped_right) or
                          (self.strand == '-' and j == self.jstart + 1 and self.cropped_left))
        if is_cropped_pos:
            tag_open = tag_open + '<div id="featId{:d}" class="SH_label">'.format(self.id_)
            tag_close = '</div>' + tag_close

        return tag_open, tag_close, class_list, self.id_
    
    
class FeatureTranslation(FeatureFloating):
    
    def __init__(self, start, end, strand, xmin, xmax, feature_bio, seqrecord=None, codon_table=1,
                 cds=True, to_stop=False, **kwargs):
        super().__init__(start=start, end=end, strand=strand, xmin=xmin, xmax=xmax, feature_bio=feature_bio, **kwargs)
        self.seqrecord = seqrecord
        self.protein_seq = (self.feature_bio.extract(self.seqrecord).seq
                            .translate(table=codon_table, cds=cds, to_stop=to_stop))
        self.protein_seq = str(self.protein_seq)
        if cds:
            # add the stop codon
            self.protein_seq = self.protein_seq + '*'
        self.show_label = True

    def get_text(self, text, j):
        j1 = j + self.xmin
        if self.strand == '+':
            j_dna = j1 - self.start_uncropped
            if VERBOSE >= 2: print("FeatureTranslation:", self.protein_seq)
            if VERBOSE >= 2: print("strand", self.strand, "j", j, "j1", j1, "j_dna", j_dna, "j_dna // 3", j_dna // 3, "aa", self.protein_seq[j_dna // 3])
            if j_dna % 3 == 1:
                return self.protein_seq[j_dna // 3]
            else:
                return ''
        elif self.strand == '-':
            j_dna = self.end_uncropped - 1 - j1
            if VERBOSE >= 2: print("strand", self.strand, "j", j, "j1", j1, "j_dna", j_dna, "j_dna // 3", j_dna // 3, "aa", self.protein_seq[j_dna // 3])
            if j_dna % 3 == 1:
                return self.protein_seq[j_dna // 3]
            else:
                return ''


class FeatureRestrictionSite(FeatureOverlay):
    def __init__(self, start, end, cut_position=None, strand='+', **kwargs):
        self.cut_position = cut_position
        super().__init__(start=start, end=end, strand=strand, **kwargs)
        self.show_label = True


    def get_css(self, track_height, css_units, colwidth, base_fontsize):
        css_class = super().get_css(track_height, css_units, colwidth, base_fontsize)
        css_class = ''
        return css_class
    
    def get_css_tag(self, j):
        j1 = j + self.xmin
        class_list = []
        style = 'normal'
        cut_left = self.cut_position[0]
        cut_right = self.cut_position[1]
        print("FeatureRestrictionSite", self.cut_position, cut_left, cut_right, j1)
        class_list = ['SH_restriction_site']
        if style == 'normal':
            if j1 == cut_left:
                class1 = 'SH_restriction_site_strandp_left'
            elif cut_left < j1 and j1 < cut_right:
                class1 = 'SH_restriction_site_strandp_middle'
            elif j1 == cut_right:
                class1 = 'SH_restriction_site_strandp_right'
            else:
                class1 = None
        if class1 is not None:
            class_list.append(class1)

        tag_open = '<div id="featId{:d}" class="{}">'.format(self.id_, ' '.join(class_list))
        tag_close = '</div>'
        
        return tag_open, tag_close, class_list, self.id_


class FeatureTick(TrackFeature):
    
    def __init__(self, position, **kwargs):
        super().__init__(start=position, end=position + 1, strand='+', name='tick', annot_class='coordinates', **kwargs)
        self.modifies_text = True
        self.show_label = True

    def get_text(self, text='', j=None):
        new_text = text
        if self.show_label:
            x = j + self.xmin
            new_text = '{:d}'.format(x)
        return new_text
            
    def get_css(self, track_height, css_units, colwidth, base_fontsize):
        css_class = ''
        return css_class
    
    def get_css_tag(self, j):
        tag_open = ''
        tag_close = ''
        # We have to add one more div level, in order to be able to float the label to the left or right
        # in the inner div, and keep the outer div to the 100% width with its background.
        if self.show_label:
            tag_open = tag_open + '<div id="featId{:d}" class="SH_ticklabel">'.format(self.id_)
            tag_close = '</div>' + tag_close
        return tag_open, tag_close, ['SH_ticklabel'], self.id_


class FeatureTickLine(FeatureTick):
    
    def __init__(self, position, **kwargs):
        super().__init__(position=position, **kwargs)
        
    def get_text(self, text='', j=None):
        new_text = text
        if self.show_label:
            new_text = '|'
        return new_text


class FeatureSequenceLabel(FeatureFixed):
    
    def __init__(self, start, end, level, **kwargs):
        super().__init__(start=start, end=end, strand='+', level=level, annot_class='sequence_label', **kwargs)
        self.show_label = True

    def get_text(self, text, j):
        if self.strand == '+':
            # write the label at the start position
            if j == self.jstart:
                return self.name
            else:
                return text
            
    def get_css(self, track_height, css_units, colwidth, base_fontsize):
        css_class = ''
        return css_class
    
    def get_css_tag(self, j):
        tag_open = ''
        tag_close = ''
        # We have to add one more div level, in order to be able to float the label to the left or right
        # in the inner div, and keep the outer div to the 100% width with its background.
        if self.show_label:
            tag_open = tag_open + '<div id="featId{:d}" class="SH_sequence_label">'.format(self.id_)
            tag_close = '</div>' + tag_close
        return tag_open, tag_close, ['SH_sequence_label'], self.id_


class FeatureRestrictionSiteLabel(FeatureSequenceLabel):

    def __init__(self, start, end, level, **kwargs):
        super().__init__(start=start, end=end, level=level, **kwargs)
        self.show_triangle = True

    def get_css_tag(self, j):

        class_list = []
        tag_open = '<div id="featId{:d}" class="SH_sequence_label">'.format(self.id_)
        tag_close = '</div>'
        
        if ((self.strand == '+' and j == self.jstart) or (self.strand == '-' and j == self.jend - 1)):
            if self.show_triangle:
                content_triangle = ('<div id="featId{:d}" class="SH_restriction_triangle">'.format(self.id_) +
                                    '<svg width="8" height="8"><polygon points="0 0, 8 0, 4 8"/></svg>' +
                                    '</div>')
                tag_open = tag_open + content_triangle

            if self.show_label:
                tag_open = tag_open + '<div id="featId{:d}" class="SH_restriction_label">'.format(self.id_)
                tag_close = '</div>' + tag_close

        return tag_open, tag_close, class_list, self.id_



class Track():
    def __init__(self, name=None, seq=None, track_class='default', fixed=True):
        self.seq = seq
        self.name = name
        self.features = []
        self.track_class = track_class
        self.fixed = fixed


class TrackGroup():
    """A TrackGroup can only have one main sequence.
    """
    
    newid = itertools.count()
    
    def __init__(self, xmin, xmax, seqrecord=None, seqrecord_start=None, name=None, show_translation=False,
                 codon_table=1):
        if name is None:
            self.name = '{:d}'.format(next(self.newid))
        else:
            self.name = name
        self.tracks = dict()
        self.length = None
        self.start = None
        self.end = None
        self.xmin = xmin
        self.xmax = xmax
        self.features_overlay = []
        self.features_floating = []
        self.features_fixed = []
        self.codon_table = codon_table
        if seqrecord is not None:
            self.read_seqrecord(seqrecord, seqrecord_start)
        self.show_translation = show_translation
        self.has_sequence_label_feat = False
        
    def get_main_track(self):
        return self.tracks[0]

    def get_minus_strand_track(self):
        ts = [track for track in self.tracks.values() if track.track_class == 'minus_strand']
        if len(ts) > 0:
            return ts[0]
        else:
            return None
    
    def remove_all_features_from_track(self, track):
        for feat in self.features_overlay:
            if feat.track == track:
                feat.track = None
        for feat in self.features_floating:
            if feat.track == track:
                feat.track = None
        for feat in self.features_fixed:
            if feat.track == track:
                feat.track = None

    def read_seqrecord(self, seqrecord, seqrecord_start):
        # We read the seqrecord and build a list of tracks, with main track being the sequence
        # and with additional tracks for floating annotations

        start = seqrecord_start
        length = len(seqrecord)
        end = start + length
        # Crop the seqrecord if needed
        self.seqrecord_cropped = seqrecord[max(self.xmin - start, 0):
                                           min(self.xmax - start, len(seqrecord))]
        self.seqrecord = seqrecord
        self.start = max(start, self.xmin)
        self.end = min(end, self.xmax)
        self.length = self.end - self.start
        self.name = self.seqrecord.id
        seq = str(self.seqrecord_cropped.seq)
        
        # Create main track
        self.tracks = {0:Track(name=self.seqrecord.id, seq=seq, fixed=True)}
        
        # We iterate over *all* the features of the SeqRecord object. Note that we do not slice the SeqRecord object,
        # because the BioPython slice will drop all the features that are not fully contained in the slice. We want to
        # keep and draw cropped features.
        nfeat = len(seqrecord.features)
        for i, feature_bio in enumerate(seqrecord.features):
            if VERBOSE >= 1: print("TrackGroup:read_seqrecord, parsing feature", i, "/", nfeat)
            
            # Adjust position of the features to the global coordinates
            # we have to shift by seqrecord_start, *not by self.start*
            start = int(feature_bio.location.start) + seqrecord_start
            end = int(feature_bio.location.end) + seqrecord_start
            strand = convert_strand_numeric_to_strand_string(feature_bio.location.strand)
            if strand is None:
                strand = '+'
            feats = []
            name = extract_feature_name(feature_bio)
            zorder = 0
            if 'z_order' in feature_bio.qualifiers:
                zorder = feature_bio.qualifiers['z_order']
            if 'zorder' in feature_bio.qualifiers:
                zorder = feature_bio.qualifiers['zorder']
            if feature_bio.type in ['CDS']:
                feats = [FeatureArrow(feature_bio=feature_bio, start=start, end=end, strand=strand,
                                      xmin=self.xmin, xmax=self.xmax, name=name,
                                      annot_class='CDS', zorder=zorder)]
            elif feature_bio.type in ['ncRNA', 'sRNA', 'tRNA', 'siRNA']:
                feats = [FeatureArrow(feature_bio=feature_bio, start=start, end=end, strand=strand,
                                      xmin=self.xmin, xmax=self.xmax, name=name,
                                      annot_class='RNA', zorder=zorder)]
            elif feature_bio.type in ['start_codon']:
                feats = [FeatureOverlay(feature_bio=feature_bio, start=start, end=end, strand=strand,
                                        xmin=self.xmin, xmax=self.xmax, name=name,
                                        annot_class='start_codon', zorder=zorder)]
            elif feature_bio.type in ['region']:
                feats = [FeatureOverlay(feature_bio=feature_bio, start=start, end=end, strand=strand,
                                        xmin=self.xmin, xmax=self.xmax, name=name,
                                        annot_class='region', zorder=zorder)]
            elif feature_bio.type in ['restriction_site']:
                cut_position = None
                if 'enzyme' in feature_bio.qualifiers.keys():
                    enzyme = feature_bio.qualifiers['enzyme'][0]
                    if enzyme not in restriction_dict.keys():
                        raise ValueError("Enzyme name", enzyme, "is not in the restriction enzyme dictionary.")
                    else:
                        cut_shift_5p = restriction_dict[enzyme]['cut_position'][0]
                        cut_shift_3p = restriction_dict[enzyme]['cut_position'][1]
                        cut_position = (start + cut_shift_5p - 1, end + cut_shift_3p)
                else:
                    enzyme = None
                    # Enzyme name can be omitted if the exact cleavage site is known
                    if 'cut_position' in feature_bio.qualifiers.keys():
                        cut_position = feature_bio.qualifiers['cut_position']
                    else:
                        raise ValueError("Either the enzyme name or the exact cut_position has to be given.")
                feats = [FeatureRestrictionSite(feature_bio=feature_bio, start=start, end=end,
                                                xmin=self.xmin, xmax=self.xmax, name=enzyme,
                                                annot_class='restriction_site', zorder=zorder,
                                                cut_position=cut_position),
                         FeatureRestrictionSite(feature_bio=feature_bio, start=start, end=end,
                                                xmin=self.xmin, xmax=self.xmax, name=enzyme,
                                                annot_class='restriction_site', zorder=zorder,
                                                cut_position=cut_position, strand='-'),
                         FeatureRestrictionSiteLabel(start=cut_position[0] + 1, end=end,
                                                     xmin=self.xmin, xmax=self.xmax, name=enzyme,
                                                     zorder=zorder, level=1)]
            else:
                feats = [FeatureOverlay(feature_bio=feature_bio, start=start, end=end, strand=strand,
                                        xmin=self.xmin, xmax=self.xmax, name=name,
                                        annot_class=feature_bio.type, zorder=zorder)]

            if feats is not []:
                for feat in feats:
                    if feat.visible:
                        if issubclass(type(feat), FeatureFloating):
                            self.features_floating.append(feat)
                        elif issubclass(type(feat), FeatureOverlay):
                            self.features_overlay.append(feat)
                        elif issubclass(type(feat), FeatureSequenceLabel):
                            self.features_fixed.append(feat)
                            self.has_sequence_label_feat = True

    def add_minus_strand(self):
        negative_track_keys = [k for k in self.tracks.keys() if k < 0]
        if len(negative_track_keys) > 0:
            for k in negative_track_keys:
                self.tracks[k - 1] = self.tracks.pop(k)
        self.tracks[-1] = Track(name='', track_class='minus_strand', fixed=True,
                                seq=str(self.seqrecord_cropped.reverse_complement().seq)[::-1])

    def add_sequence_label_track(self):
        positive_track_keys = [k for k in self.tracks.keys() if k > 0]
        if len(positive_track_keys) > 0:
            # Shift all upper tracks one unit up
            for k in sorted(positive_track_keys, reverse=True):
                self.tracks[k + 1] = self.tracks.pop(k)
        self.tracks[1] = Track(name='', track_class='sequence_label', fixed=True)

    def add_restriction_sites(self):
        if self.tracks[1].track_class != 'sequence_label':
            self.add_sequence_label_track()

        # Search for restriction sites
        seq = self.get_main_track().seq

        for enzyme, restrict_dict in restriction_dict.items():
            re_motif = re.compile(restrict_dict['motif'])
            for m in re_motif.finditer(seq):
                if VERBOSE >= 2: print(enzyme, re_motif, m.start(), m.end(), m.group())
                start = m.start()
                end = m.end()
                if VERBOSE >= 2: print("restriction site motif found", enzyme, start, end)
                # Test if the feature has already been added
                exists = any([feat.start == start and issubclass(type(feat), FeatureRestrictionSite)
                              for feat in self.features_overlay])
                if not exists:
                    if VERBOSE >= 2: print("adding restriction site feature", enzyme, start, end)
                    cut_shift_5p = restrict_dict['cut_position'][0]
                    cut_shift_3p = restrict_dict['cut_position'][1]
                    cut_position = (start + cut_shift_5p - 1, end + cut_shift_3p)
                    feat = FeatureRestrictionSite(feature_bio=None, start=start, end=end,
                                                  xmin=self.xmin, xmax=self.xmax, name=enzyme,
                                                  annot_class='restriction_site',
                                                  cut_position=cut_position)
                    self.features_overlay.append(feat)
                    feat = FeatureRestrictionSite(feature_bio=None, start=start, end=end,
                                                  xmin=self.xmin, xmax=self.xmax, name=enzyme,
                                                  annot_class='restriction_site',
                                                  cut_position=cut_position, strand='-')
                    self.features_overlay.append(feat)
                    feat = FeatureRestrictionSiteLabel(start=cut_position[0] + 1, end=end,
                                                       xmin=self.xmin, xmax=self.xmax, name=enzyme,
                                                       level=1)
                    self.features_fixed.append(feat)
                    self.has_sequence_label_feat = True
        self.organize_tracks()


    def organize_tracks(self):
        if VERBOSE >= 1: print("TrackGroup:organize_tracks")

        # Remove sequence label track
        for i in list(self.tracks.keys()):
            if self.tracks[i].track_class == 'sequence_label':
                self.remove_all_features_from_track(self.tracks[i])
                self.tracks.pop(i)
                if VERBOSE >= 2: print("organize_tracks: removing track", i)

        # Create the fixed sequence label track and add features
        # Remark: For some reason self.has_sequence_label_feat here is set back to FALSE while it was True just before
        self.has_sequence_label_feat = any([issubclass(type(feat), FeatureSequenceLabel) for feat
                                            in self.features_fixed])
        if self.has_sequence_label_feat:
            if VERBOSE >= 2: print("adding sequence label track")
            self.add_sequence_label_track()
            for feat in self.features_fixed:
                if issubclass(type(feat), FeatureSequenceLabel):
                    self.tracks[1].features.append(feat)
                    feat.track = self.tracks[1]

        # Remove all floating tracks
        for i in list(self.tracks.keys()):
            if VERBOSE >= 2: print("organize_tracks: self.tracks[i].fixed", i, self.tracks[i].fixed)
            if not self.tracks[i].fixed:
                self.remove_all_features_from_track(self.tracks[i])
                self.tracks.pop(i)
                if VERBOSE >= 2: print("organize_tracks: removing track", i)

        # Organize floating features into multiple tracks, computing smart layout
        if len(self.features_floating) > 0:
            feat_len_list = [feat.length for feat in self.features_floating]
            # Sort features with longer one first
            self.features_floating = sorted(self.features_floating, key=lambda feat: feat.length)[::-1]
            
            base_level = max(self.tracks.keys()) + 1
            features_levels = compute_features_levels(self.features_floating, base_level=base_level)
            for feat, level in features_levels.items():
                feat.level = level
            for feat in self.features_floating:
                level = feat.level
                if level not in self.tracks.keys():
                    self.tracks[level] = Track(fixed=False)
                # Add feature to track
                self.tracks[level].features.append(feat)
                feat.track = self.tracks[level]

            if self.show_translation:
                # Add additional tracks for features that need several tracks together,
                # for example, translation sequence below the arrow for CDS.
                levels = set([feat.level for feat in self.features_floating if feat.annot_class == 'CDS'])
                tracks_copy = copy(self.tracks)
                self.tracks = dict()
                new_level = min(tracks_copy.keys())
                for old_level in sorted(tracks_copy.keys()):
                    if old_level in levels:
                        # add new empty track below
                        self.tracks[new_level] = Track(fixed=False)
                        self.tracks[new_level + 1] = tracks_copy[old_level]
                        for feat in self.tracks[new_level + 1].features:
                            feat.level = new_level + 1
                            feat.track = self.tracks[new_level + 1]
                        # add translation
                        for feat in self.tracks[new_level + 1].features:
                            if feat.annot_class == 'CDS':
                                feat2 = FeatureTranslation(feature_bio=feat.feature_bio,
                                                           start=feat.start, end=feat.end, strand=feat.strand,
                                                           xmin=feat.xmin, xmax=feat.xmax,
                                                           name=feat.name, annot_class='translation',
                                                           zorder=feat.zorder, seqrecord=self.seqrecord,
                                                           codon_table=self.codon_table)
                                feat2.level = new_level
                                feat2.track = self.tracks[new_level]
                                self.tracks[new_level].features.append(feat2)

                        new_level += 2
                    else:
                        # copy track
                        self.tracks[new_level] = tracks_copy[old_level]
                        for feat in self.tracks[new_level].features:
                            feat.level = new_level
                            feat.track = self.tracks[new_level]
                        new_level += 1

        # Remove all overlay features
        for feat in self.features_overlay:
            if feat.track is not None:
                feat.track.features.remove(feat)
                if VERBOSE >= 2: print("organize_tracks: removing overlay feature", feat,
                                       "from track #", feat.track)

        # if len(self.features_overlay) > 0:
        #     feat_len_list = [feat.length for feat in self.features_overlay]
        #     # Sort features with longest one first, which will have the lowest zorder
        #     self.features_overlay = sorted(self.features_overlay, key=lambda feat: feat.length)[::-1]
        #     base_level = 0
        #     features_levels = compute_features_levels(self.features_overlay, base_level=base_level)
        #     for feat, level in features_levels.items():
        #         # Not for features that have a specific z_order
        #         if feat.zorder == 0:
        #             feat.zorder = level
        #         else:
        #             # Keep the user-defined zorder
        #             feat.zorder = feat.zorder
        #     for feat in self.features_overlay:
        #         feat.level = 0
        #         assert 0 in self.tracks.keys()
        #         # Add feature to track
        #         self.tracks[0].features.append(feat)
        #         feat.track = self.tracks[0]

        # Add overlay features to the main sequence track
        if self.get_main_track() is not None:
            features = [feat for feat in self.features_overlay if feat.strand == '+']
            if len(features) > 0:
                feat_len_list = [feat.length for feat in features]
                # Sort features with longest one first, which will have the lowest zorder
                features = sorted(features, key=lambda feat: feat.length)[::-1]
                base_level = 0
                features_levels = compute_features_levels(features, base_level=base_level)
                for feat, level in features_levels.items():
                    # Not for features that have a specific z_order
                    if feat.zorder == 0:
                        feat.zorder = level
                    else:
                        # Keep the user-defined zorder
                        feat.zorder = feat.zorder
                for feat in features:
                    # Add feature to track
                    self.get_main_track().features.append(feat)
                    feat.track = self.get_main_track()

        # Add overlay features to the minus sequence track
        if self.get_minus_strand_track() is not None:
            features = [feat for feat in self.features_overlay if feat.strand == '-']
            if len(features) > 0:
                feat_len_list = [feat.length for feat in features]
                # Sort features with longest one first, which will have the lowest zorder
                features = sorted(features, key=lambda feat: feat.length)[::-1]
                base_level = 0
                features_levels = compute_features_levels(features, base_level=base_level)
                for feat, level in features_levels.items():
                    # Not for features that have a specific z_order
                    if feat.zorder == 0:
                        feat.zorder = level
                    else:
                        # Keep the user-defined zorder
                        feat.zorder = feat.zorder
                for feat in features:
                    # Add feature to track
                    self.get_minus_strand_track().features.append(feat)
                    feat.track = self.get_minus_strand_track()

        for feat in self.features_overlay:
            print("features_overlay", feat.track, feat)
        for feat in self.features_floating:
            print("features_floating", feat.track, feat)
        for feat in self.features_fixed:
            print("features_fixed", feat.track, feat)
                
    def get_dimensions(self):
        height = max(self.tracks.keys()) - min(self.tracks.keys()) + 1
        if self.start is not None and self.end is not None:
            self.length = self.end - self.start
        else:
            self.length = None
        return self.length, height
                

class Element():
    
    def __init__(self, text='', x=None, i=None, j=None, features=list(),
                 track_group_id=None, track_id=None, trackgroup=None, track=None):
        self.original_text = text
        self.text = text
        self.i = i
        self.j = j
        self.x = x
        self.features = features
        self.annot_class_list = list()
        self.track_group_id = track_group_id
        self.track_id = track_id
        self.trackgroup = trackgroup
        self.track = track
    
    def apply_features(self):
        # resolve zorder. last one, highest zorder, is the most important.
        self.features = sorted(self.features, key=lambda feat: feat.zorder)
        for feat in self.features:
            self.text = feat.get_text(text=self.text, j=self.j)
            self.annot_class_list.append(feat.annot_class)
            
    def __str__(self):
        return ('Element: x={:d}, i={:d}, j={:d}, track_id={:d}, track_group_id={:d}, text:{}, n features:{:d}, features_id:{:d}, features:{}'
                .format(self.x, self.i, self.j, self.track_id, self.track_group_id, self.text,
                        len(self.features), id(self.features), [f.name for f in self.features]))
        

class TrackHighway():
    
    def __init__(self, seqrecord, seqrecord_start=0, xmin=None, xmax=None, show_translation=False, codon_table=1):
        # xmin, xmax are the global coordinates system and are **FIXED**
        self.trackgroup_list = []
        self.xmin = xmin
        self.xmax = xmax
        if self.xmin is None:
            self.xmin = 0
        if xmax is None:
            self.xmax = len(seqrecord)
        self.add_seqrecord(seqrecord=seqrecord, seqrecord_start=seqrecord_start, show_translation=show_translation,
                           codon_table=codon_table)
        self.hspace = 1
        self.array = None
        self.length = 0
        self.height = 0
        self.header_list = []
        
    def add_seqrecord(self, seqrecord, seqrecord_start=0, show_translation=False, codon_table=1):
        assert seqrecord_start >= 0
        trackgroup = TrackGroup(xmin=self.xmin, xmax=self.xmax, seqrecord=seqrecord,
                                seqrecord_start=seqrecord_start, show_translation=show_translation, codon_table=codon_table)
        trackgroup.organize_tracks()
        self.trackgroup_list.append(trackgroup)
        
    def add_coordinates_track(self, ticks=True, step=10):
        trackgroup = TrackGroup(xmin=self.xmin, xmax=self.xmax)
        trackgroup.length = self.xmax - self.xmin
        trackgroup.tracks = {0:Track(name='', track_class='coordinates')}
        if ticks:
            trackgroup.tracks[1] = Track(track_class='coordinates')
        
        for x in range(self.xmin, self.xmax):
            j = x - self.xmin
            if x % step == 0:
                feat = FeatureTick(position=x, xmin=self.xmin, xmax=self.xmax)
                trackgroup.features_floating.append(feat)
                trackgroup.tracks[0].features.append(feat)
                if ticks:
                    feat = FeatureTickLine(position=x, xmin=self.xmin, xmax=self.xmax)
                    trackgroup.features_floating.append(feat)
                    trackgroup.tracks[1].features.append(feat)

        self.trackgroup_list.append(trackgroup)
        
    def add_minus_strand(self, trackgroup_name=None):
        if trackgroup_name is not None:
            trackgroup = [tg for tg in self.trackgroup_list if tg.name == trackgroup_name]
            if len(trackgroup) != 1:
                raise ValueError("trackgroup_name was not found in the list of existing trackGroups.")
            trackgroup = trackgroup[0]
        else:
            trackgroup = self.trackgroup_list[0]
        trackgroup.add_minus_strand()

    def add_restriction_sites(self):
        self.trackgroup_list[0].add_restriction_sites()
        
    def update_array(self):
        # generate an array representation of the tracks and annotations
        self.length = 0
        self.height = 0
        self.header_list = []
        
        # compute the xmin and xmax of all the trackgroups
        xmin_all = self.trackgroup_list[0].xmin
        xmax_all = self.trackgroup_list[0].xmax
        for track_group_id, trackgroup in enumerate(self.trackgroup_list):
            length, height = trackgroup.get_dimensions()
            self.height += height
            xmin_all = min(xmin_all, trackgroup.xmin)
            xmax_all = max(xmax_all, trackgroup.xmax)
        
        self.length = self.xmax - self.xmin
        self.array = [[Element(x=self.xmin + j, i=i, j=j, features=list())
                       for j in range(self.length)]
                      for i in range(self.height)]
        self.header_list = self.height*['']
        
        iOffset = 0
        for track_group_id, trackgroup in enumerate(self.trackgroup_list):
            ntracks = len(trackgroup.tracks)
            length, height = trackgroup.get_dimensions()
            level_min = min(trackgroup.tracks.keys())
            level_max = max(trackgroup.tracks.keys())
            if VERBOSE >= 2:
                print("TrackHighway:update_array ntracks", ntracks, "length", length, "height", height,
                      "level_min", level_min, "level_max", level_max)
            
            for level, track in trackgroup.tracks.items():
                i = level_max - level + iOffset
                
                # Add header
                if track.name is not None:
                    self.header_list[i] = track.name

                # Fill position of empty elements
                for j in range(self.length):
                    self.array[i][j].track_group_id = track_group_id
                    self.array[i][j].track_id = level
                    self.array[i][j].trackgroup = trackgroup
                    self.array[i][j].track = track
                
                if track.seq is not None:
                    for j0 in range(length):
                        jOffset = trackgroup.start - self.xmin
                        j = j0 + jOffset
                        if VERBOSE >= 2: print("TrackHighway:update_array level", level, "track", track, "i", i, 'j', j)
                        if j >= 0:
                            if VERBOSE >= 2: print("track.seq[j0]",track.seq[j0])
                            char = track.seq[j0]
                            self.array[i][j] = Element(text=char, x=self.xmin + j, i=i, j=j, features=[],
                                                       track_group_id=track_group_id, track_id=level,
                                                       trackgroup=trackgroup, track=track)

                for feat in track.features:
                    # Crop the feature
                    for j in range(feat.jstart, feat.jend):
                        ((self.array[i][j]).features).append(feat)
            
            iOffset += height
            
        # Apply the features to each element, taking into account the zorders
        for i in range(self.height):
            for j in range(self.length):
                self.array[i][j].apply_features()
                if VERBOSE >= 2: print(i, j, self.array[i][j].text, self.array[i][j].annot_class_list)
            
    def to_html(self, **kwargs):
        writer = HTMLWriter(self, **kwargs)
        return writer.to_html()


class HTMLWriter():
    
    def __init__(self, trackhighway, mode='fixed_table', table_width=70,
                 show_header=True, show_header_every_line=False,
                 css_track_height=1, css_seq_colwidth=0.67, css_base_fontsize=1, css_units='rem'):
        self.mode = mode
        self.trackhighway = trackhighway
        self.trackhighway.update_array()
        # Collect the list of all features
        self.all_features = set([feat for row in self.trackhighway.array for e in row for feat in e.features])
        self.css_base_fontsize = css_base_fontsize    # all sizes will be scaled based on this size
        self.css_units = css_units
        
        # For chrome, standalone html, css_track_height 1, css_seq_colwidth 0.6
        # For embedded jupyter on firefox, css_track_height 1, css_seq_colwidth 0.6???
        self.css_track_height = css_track_height
        self.css_seq_colwidth = css_seq_colwidth
        self.show_header = show_header
        self.show_header_every_line = show_header_every_line
        self.css_header_colwidth = 20
        self.table_width = table_width    # in columns (=sequence characters)

    def seqhighway_html(self):
        array = self.trackhighway.array
        indent = '    '


        ### FIXED TABLE MODE
        if self.mode == 'fixed_table':
            s = ''
            indent_container = ''
            s += indent_container + '<div class="SH_container">\n'
            indent_table = indent
            s += indent_table + '<table class="SH_table">\n'
            # colOpen = '<div class="SH_col">'
            # colClose = '</div>'
            nrows = self.trackhighway.height
            ncols = self.trackhighway.length
            n_chunk = ncols // self.table_width + 1
            for k_chunk in range(n_chunk):
                for i in range(nrows):
                    i_table = k_chunk*(nrows + 1) + i

                    indent_row = indent_table + indent
                    # Apply CSS of the track to the table row
                    element = array[i][0]
                    track_class = element.track.track_class
                    if track_class == 'default':
                        track_class_css = 'SH_track_default'
                    elif track_class == 'coordinates':
                        track_class_css = 'SH_track_default SH_track_coordinates'
                    elif track_class == 'minus_strand':
                        track_class_css = 'SH_track_default SH_track_minus_strand'
                    elif track_class == 'sequence_label':
                        track_class_css = 'SH_track_default SH_track_sequence_label'
                    s += (indent_row +
                          '<tr class="SH_tr SH_trackgroup{:d} SH_track{:d} {}">\n'
                          .format(element.track_group_id, element.track_id, track_class_css))
                    
                    indent_element = indent_row + indent
                    element_open = '<div class="SH_element">'
                    element_close = '</div>'
                    if self.show_header:
                        # write header column element
                        text = self.trackhighway.header_list[i]
                        element = array[i][0]
                        track_class = element.track.track_class
                        features = element.features
                        ss = text
                        ss = element_open + ss + element_close
                        ss = (indent_element +
                              '<td class="SH_td SH_header_col">' +
                              ss + '</td>\n')
                        s += ss

                    for j_table in range(self.table_width):

                        j = k_chunk*self.table_width + j_table
                        if j < ncols:
                            element = array[i][j]
                            features = element.features
                            if VERBOSE >= 2: print("HTMLWriter:seqhighway_html i, j, features", i, j, [f.name for f in features])
                            ss = element.text
                            # The most internal tag wins over the other tags. We use this rule in order to set the
                            # priority of the feature annotations.
                            # Features are ordered following zorder: last one, highest zorder, is the most important.
                            for feat in features[::-1]:
                                tag_open, tag_close, class_list, id_ = feat.get_css_tag(j)
                                ss = tag_open + ss + tag_close
                            ss = element_open + ss + element_close
                            # We also apply the css classes to the td, because width, height and padding
                            # can only be effective for the td element in the table.
                            css_classes_td = []
                            for feat in features[::-1]:
                                tag_open, tag_close, class_list, id_ = feat.get_css_tag(j)
                                css_classes_td = css_classes_td + class_list
                            css_classes_td = ' '.join(css_classes_td)
                            ss = (indent_element +
                                  '<td class="SH_td {}">'.format(css_classes_td) +
                                  ss + '</td>\n')
                            s += ss
                    
                    s += indent_row + '</tr>\n'
                # add empty row
                s += indent_row + '<tr class="SH_spacing_row"></tr>\n'


            s += indent_table + '</table>\n'
            s += indent_container + '</div>\n'


        ### INLINE-BLOCK MODE
        elif self.mode == 'inline-block':
            s = ''
            s += '<span class="SH_highway">\n'
            colOpen = '<div class="SH_col">'
            colClose = '</div>'
            nrows = self.trackhighway.height
            ncols = self.trackhighway.length
            for j in range(ncols):
                indentPrefix = ' -->' if j > 0 else '    '
                s += indentPrefix + colOpen + '\n'
                for i in range(nrows):
                    element = array[i][j]
                    features = element.features
                    track_class = element.track.track_class
                    if VERBOSE >= 2: print("HTMLWriter:seqhighway_html i, j, features", i, j, [f.name for f in features])
                    ss = element.text
                    # The most internal tag wins over the other tags. We use this rule in order to set the
                    # priority of the feature annotations.
                    # Features are ordered following zorder: last one, highest zorder, is the most important.
                    for feat in features[::-1]:
                        tag_open, tag_close, class_list, id_ = feat.get_css_tag(j)
                        ss = tag_open + ss + tag_close
                    element_open = '<div class="SH_element">'
                    element_close = '</div>'
                    ss = element_open + ss + element_close
                    
                    # Apply trackgroup and track classes
                    if track_class == 'default':
                        track_class_css = 'SH_track_default'
                    elif track_class == 'coordinates':
                        track_class_css = 'SH_track_default SH_track_coordinates'
                    elif track_class == 'minus_strand':
                        track_class_css = 'SH_track_default SH_track_minus_strand'
                    elif track_class == 'sequence_label':
                        track_class_css = 'SH_track_default SH_track_sequence_label'
                    ss = (indent + indent +
                          '<div class="SH_trackgroup{:d} SH_track{:d} {}">'
                          .format(element.track_group_id, element.track_id, track_class_css) +
                          ss + '</div>' + '\n')
                    s += ss
                
                suffix = '<!--' if j < ncols - 1 else ''
                s += indent + colClose + suffix + '\n'
            s += '</span>'
            
            if self.show_header:
                s_header = ''
                if self.show_header_every_line:
                    # we will put all the header in the left column,
                    # and all the sequence in the right column.
                    n_header = self.seq_max_HTML_rows
                    s_header += '<div class="SH_container">\n'
                    indent1 = indent
                else:
                    n_header = 1
                    indent1 = ''
                s_header += indent1 + '<div class="SH_header_col">\n'
                
                for i_header in range(n_header):
                    s_row = ''
                    for i in range(nrows):
                        text = self.trackhighway.header_list[i]
                        element = array[i][0]
                        track_class = element.track.track_class
                        if VERBOSE >= 2: print("HTMLWriter:seqhighway_html,header i, j, features", i, 0, self.trackhighway.header_list[i])
                        element_open = '<div class="SH_element">'
                        element_close = '</div>'
                        s_element = element_open + text + element_close
                        # Apply trackgroup and track classes
                        if track_class == 'default':
                            track_class_css = 'SH_track_default'
                        elif track_class == 'coordinates':
                            track_class_css = 'SH_track_default SH_track_coordinates'
                        elif track_class == 'minus_strand':
                            track_class_css = 'SH_track_default SH_track_minus_strand'
                        elif track_class == 'sequence_label':
                            track_class_css = 'SH_track_default SH_track_sequence_label'

                        s_element = (indent1 + indent + indent + indent +
                                     '<div class="SH_trackgroup{:d} SH_track{:d} {}">'
                                     .format(element.track_group_id, element.track_id, track_class_css) +
                                     s_element + '</div>' + '\n')
                        s_row += s_element
                    
                    s_row = (indent1 + indent + '<div class="SH_header_row{:d}">\n'.format(i_header) +
                             s_row +
                             indent1 + indent + '</div>\n')
                    s_header += s_row
                
                s_header += indent1 + '</div>' + '\n'
                
                s = s_header + '<div class="SH_seqhighway_col">\n' + s + '</div>\n'
                if self.show_header_every_line:
                    s += '</div>\n'

        return s


    def css_style(self):
        self.css_header_colwidth = self.css_seq_colwidth*max([len(header) for header in self.trackhighway.header_list])

        style = ''
        if self.mode == 'fixed_table':
            style += """
<style>
    .SH_container {{
        font-size: {}{};
        overflow-x: auto;
    }}""".format(self.css_base_fontsize, self.css_units)
            style += """
    .SH_table {
        table-layout: fixed;
        border-width: 0;
        border-collapse: collapse;
        border-spacing: 0;
        height: 0;    <!-- This is needed in order to set the height of 100% for rows -->
    }
    .rendered_html tbody tr.SH_tr {
        background: transparent;    <!-- set back the background to transparent,
                                         not to interfere with annotation colors -->
    }"""
            style += """
    td.SH_header_col {{
        min-width: {}{};
        text-align: right;
        padding-right: {}{};
    }}""".format(self.css_base_fontsize*self.css_header_colwidth, self.css_units,
                 2*self.css_base_fontsize*self.css_seq_colwidth, self.css_units)
            style += """
    .SH_track_default {{
        font-family:monospace;
        font-size: {}{};
    }}""".format(self.css_base_fontsize, self.css_units,)
            style += """
    tr.SH_spacing_row {{
        height: {}{};
    }}""".format(2*self.css_base_fontsize*self.css_track_height, self.css_units,)
            style += """
    .SH_track_default > td.SH_td {{
        height: {}{};
        padding: 0;
        padding-top: {}{};
        padding-bottom: {}{};
    }}""".format(# adding the padding value to the height of the row
                 (1 + 0.175 + 0.175)*self.css_base_fontsize*self.css_track_height, self.css_units,
                 0.175*self.css_base_fontsize*self.css_track_height, self.css_units,
                 0.175*self.css_base_fontsize*self.css_track_height, self.css_units,)
            style += """
    .SH_track0 > td.SH_td {{
        height: {}{};
        padding-top: {}{};
        padding-bottom: {}{};
    }}""".format(# padding on the top and bottom
                 (1 + 0.25 + 0.25)*self.css_base_fontsize*self.css_track_height, self.css_units,
                 0.25*self.css_base_fontsize*self.css_track_height, self.css_units,
                 0.25*self.css_base_fontsize*self.css_track_height, self.css_units,)
            style += """
    .SH_track_minus_strand > td.SH_td {{
        height: {}{};
        padding-bottom: {}{};
    }}""".format(# no padding
                 (1 + 0.2)*self.css_base_fontsize*self.css_track_height, self.css_units,
                 0.2*self.css_base_fontsize*self.css_track_height, self.css_units)
            style += """
    .SH_track_coordinates > td.SH_td {{
        height: {}{};
        padding-top: 0;
        padding-bottom: 0;
    }}""".format(0.8*self.css_base_fontsize*self.css_track_height, self.css_units,)
            style += """
    .SH_track_sequence_label > td.SH_td {{
        height: {}{};
        padding-bottom: 0;
    }}""".format(1*self.css_base_fontsize*self.css_track_height, self.css_units,)
            style += """
    div.SH_sequence_label {{
        position: relative;
        overflow-x: visible;
        float: left;
    }}""".format()
            style += """
    .SH_restriction_triangle {{
        display: inline-block;
        position: absolute;
        bottom: 0;
        left: -4;
    }}
    .SH_restriction_label {{
        display: inline-block;
        position: relative;
        left: 5;
    }}""".format()
            style += """
    td.SH_td.SH_CDS {{
    }}""".format()
            style += """
    td.SH_td.SH_translation {{
        padding-top: 0;
    }}""".format()
            style += """
    .SH_element {{
        white-space: nowrap;
        padding: 0px;
        margin-top: 0px;
        margin-bottom: 0px;
        width: {}{};
        height:100%;
        text-align: center;
    }}""".format(self.css_base_fontsize*self.css_seq_colwidth, self.css_units,)
            style += """
    .SH_label {{
        font-size: {}{};
    }}""".format(# IMPORTANT: the table row will adjust their height depending on the larger
                 # cell. Even if we set manually the height of the cell and internal div,
                 # it seems that the row adapts. Thus, it is important that no element,
                 # including the text labels, do not overcome the height limit.
                 0.8*self.css_base_fontsize, self.css_units,)
            style += """
    .SH_ticklabel {{
        font-size: {}{};
        z-index: 10;
        overflow-x: visible;
        text-align: center;
        margin-left: -100%;
        margin-right: -100%;
    }}\n""".format(0.6*self.css_base_fontsize*self.css_track_height, self.css_units)
            style += """
    .SH_track0 > td.SH_restriction_site {{
        border-right: 0;
        border-left: 0;
        border-top: 0;
        border-bottom: 0;
        border-right-style: solid;
        border-left-style: solid;
        border-top-style: solid;
        border-bottom-style: solid;
    }}
    .SH_track0 > td.SH_restriction_site_strandp_left {{
        border-right: 1px;
        border-right-style: solid;
    }}
    .SH_track0 > td.SH_restriction_site_strandp_middle {{
        border-bottom: 1px;
        border-bottom-style: solid;
    }}
    .SH_track0 > td.SH_restriction_site_strandp_right {{
    }}
    .SH_track_minus_strand > td.SH_restriction_site_strandp_left {{
    }}
    .SH_track_minus_strand > td.SH_restriction_site_strandp_middle {{
    }}
    .SH_track_minus_strand > td.SH_restriction_site_strandp_right {{
        border-left: 1px;
        border-left-style: solid;
    }}\n""".format()



        elif self.mode == 'inline-block':
            style += ("""
<style>
    .SH_container {{
        font-size: {}{};
    }}
    .SH_header_col {{
        width: {}{};
        line-height: {}{};
        float: left;
        display: inline-block;
        text-align: right;
        padding: 0px;
        margin-left: 0px;
        margin-right: {}{};
        border: 0px;
        border-style: solid;
        vertical-align: bottom;
        }}"""
                     .format(self.css_base_fontsize, self.css_units,
                             self.css_base_fontsize*self.css_header_colwidth, self.css_units,
                             self.css_base_fontsize*self.css_track_height, self.css_units,
                             self.css_base_fontsize*1, self.css_units
                             ))
                    
            seq_total_length = self.css_base_fontsize*self.css_seq_colwidth*(self.trackhighway.length + 2)
            # minimum width of the seq div
            seq_div_min_width = 10
            seq_max_HTML_rows = int((seq_total_length // seq_div_min_width) + 1)
            # maximum width of the seq div
            seq_div_max_width = 200
            seq_min_HTML_rows = int((seq_total_length // seq_div_max_width) + 1)
            self.seq_max_HTML_rows = seq_max_HTML_rows
            self.seq_min_HTML_rows = seq_min_HTML_rows
            if VERBOSE >= 2: print("HTMLWriter:seqhighway_html seq_min_HTML_rows, seq_max_HTML_rows", seq_min_HTML_rows, seq_max_HTML_rows)
            break_points = {i_header:seq_total_length/i_header + self.css_base_fontsize*(self.css_header_colwidth + 1)
                           for i_header in range(seq_min_HTML_rows, seq_max_HTML_rows + 1)}
            if VERBOSE >= 2: print("HTMLWriter:seqhighway_html seq_total_length, break_points", seq_total_length, break_points)
            # Simplify break points: we merge together break points that are very close
            break_points2 = []
            break_last = 1e9
            break0 = None
            i_list = None
            i_first = min(break_points.keys())
            for i, break_point in break_points.items():
                if break_last - break_point < 2:
                    # append to previous group
                    i_list.append(i)
                else:
                    # push previous group
                    if i > i_first:
                        if VERBOSE >= 2: print(i, break_point, (i_list, break0))
                        break_points2.append((i_list, break0))
                    # new group
                    i_list = [i]
                    break0 = break_point
                break_last = break_point
            break_points2.append((i_list, break0))
            if VERBOSE >= 2: print("HTMLWriter:seqhighway_html break_points2", break_points2)
        
            # First header row
            style += ("""
    .SH_header_row0 {{
        display: block;
        margin-bottom: {}{};
    }}\n"""
                      .format(i, (1.5 + 0.25)*self.css_base_fontsize*self.css_track_height, self.css_units))
            # From 2nd header row
            if self.show_header_every_line:
                # We always show the first group of header rows
                for i in range(1, self.seq_min_HTML_rows):
                    style += ("""
    .SH_header_row{:d} {{
        display: block;
        margin-bottom: {}{};
    }}\n"""
                              .format(i, (1.5 + 0.25)*self.css_base_fontsize*self.css_track_height, self.css_units))
                
                for i_list, break_point in break_points2:
                    for i in i_list:
                        style += "    .SH_header_row{:d} {{\n        display: none;\n    }}\n".format(i)
                    style += """    @media (max-width: {}{}) {{\n""".format(break_point, self.css_units)
                    for i in i_list:
                        style += ("""        .SH_header_row{:d} {{
            display: block;
            margin-bottom: {}{};
        }}\n"""
                                  .format(i, (1.5 + 0.25)*self.css_base_fontsize*self.css_track_height, self.css_units))
                    style += "    }\n"
                     
            style += ("""
    .SH_seqhighway_col {{
        width: 100%;
        color: "black";
    }}
    .SH_col {{
        display: inline-block;
        line-height: {}{};
        width: {}{};
        padding: 0px;
        margin-left: 0px;
        margin-right: 0px;
        margin-bottom: {}{};
        border: 0px;
        border-style: solid;
        text-align: center;
        vertical-align: bottom;
    }}
    .SH_highway {{
        word-break: break-all;
        text-align: left;
    }}
    .SH_track_default {{
        font-family:monospace;
        font-size: {}{};
        height: {}{};
        margin-bottom: {}{};
        position: relative;
    }}
    .SH_track0 {{
        height: {}{};
        margin-top: {}{};
    }}
    .SH_track_sequence_label {{
        margin-top: {}{};
        margin-bottom: {}{};
    }}
    .SH_restriction_triangle {{
        display: inline-block;
        position: absolute;
        bottom: 0;
        left: -4 px;
        margin-right: 6 px;
    }}
    .SH_track_coordinates {{
        height: {}{};
        margin-top: {}{};
        margin-bottom: {}{};
    }}
    .SH_track_minus_strand {{
        margin-top: {}{};
    }}
    .SH_element {{
        white-space: nowrap;
        padding: 0px;
        margin-top: 0px;
        margin-bottom: 0px;
        position:absolute;
        width:100%;
        height:100%;
        bottom:0;
    }}
    .SH_ticklabel {{
        position: relative;
        font-size: {}{};
        z-index: 10;
        overflow-x: visible;
        float: none;
        margin-left: -100%;
        margin-right: -100%;
        text-align: center;
    }}
"""
                .format(#SH_col
                        self.css_base_fontsize*self.css_track_height, self.css_units,
                        self.css_base_fontsize*self.css_seq_colwidth, self.css_units,
                        1.5*self.css_base_fontsize*self.css_track_height, self.css_units,
                        #SH_track_default
                        self.css_base_fontsize*1, self.css_units,
                        self.css_base_fontsize*self.css_track_height, self.css_units,
                        0.25*self.css_base_fontsize*self.css_track_height, self.css_units,
                        #SH_track0
                        self.css_base_fontsize*self.css_track_height, self.css_units,
                        0.5*self.css_base_fontsize*self.css_track_height, self.css_units,
                        #SH_track_sequence_label
                        0*self.css_base_fontsize*self.css_track_height, self.css_units,
                        0*self.css_base_fontsize*self.css_track_height, self.css_units,
                        #SH_trackCoordinates
                        0.8*self.css_base_fontsize*self.css_track_height, self.css_units,
                        0*self.css_base_fontsize*self.css_track_height, self.css_units,
                        0.25*self.css_base_fontsize*self.css_track_height, self.css_units,
                        #SH_track_minus_strand
                        0.25*self.css_base_fontsize*self.css_track_height, self.css_units,
                        #SH_ticklabel
                        0.6*self.css_base_fontsize*self.css_track_height, self.css_units))


        for feat in self.all_features:
            style += feat.get_css(track_height=self.css_base_fontsize*self.css_track_height,
                                  css_units=self.css_units,
                                  colwidth=self.css_base_fontsize*self.css_seq_colwidth,
                                  base_fontsize=self.css_base_fontsize)
        style += "\n</style>\n"

        return style


    def to_html(self):
        s = ''
        s += '<meta charset="utf-8"/>\n'    # improves UTF-8 comptability on Google chrome, although UTF8 should be default
        s += '<html>\n'
        s += self.css_style()
        s += self.seqhighway_html()
        s += '\n</html>'
        return s
