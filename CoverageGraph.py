from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord

import matplotlib.patches as mpatches
from matplotlib.pyplot import *

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.it'] = 'Helvetica:italic'
rcParams['mathtext.rm'] = 'Helvetica'

class MyCustomTranslator(BiopythonTranslator):
    def compute_feature_color(self, features):
        if features.type == "V1":
            return "red"
        elif features.type == "V2":
            return "blue"
        elif features.type == "V3":
            return "green"
        elif features.type == "V4":
            return "yellow"
        elif features.type == "V5":
            return "orange"
        elif features.type == "V6":
            return "purple"
        elif features.type == "V7":
            return "brown"
        else:
            return "gold"


    def compute_feature_box_color(self, feature):
        """Compute a box_color for this feature."""
        return "auto"

    def compute_filtered_features(self, features):
        """Do not display promoters. Just because."""
        return [
            feature for feature in features
            if (feature.type != "CDS")
            if (feature.type != "misc_feature")
            or ("BamHI" in str(feature.qualifiers.get("label", '')))
        ]

    def compute_feature_label(self, feature):
        return None

    def compute_feature_linewidth(self, feature):
        """Compute the edge width of the feature's arrow/rectangle."""
        return 1.0

    def compute_feature_box_color(self, feature):
        """Compute a box_color for this feature."""
        return "auto"

    def compute_feature_box_linewidth(self, feature):
        """Compute a box_linewidth for this feature."""
        return 0.3

graphic_record = MyCustomTranslator().translate_record("/Users/administrator/Desktop/SS14_tprD_new_2.gb")

zoom_start, zoom_end = 597, 4796
graphic_record = graphic_record.crop((zoom_start, zoom_end))

bx, _ = graphic_record.plot(figure_width=7,figure_height=1.5, strand_in_label_threshold=7)

red_patch = mpatches.Patch(color='red', label='V1')
blue_patch = mpatches.Patch(color='blue', label='V2')
green_patch = mpatches.Patch(color='green', label='V3')
yellow_patch = mpatches.Patch(color='yellow', label='V4')
orange_patch = mpatches.Patch(color='orange', label='V5')
purple_patch = mpatches.Patch(color='purple', label='V6')
brown_patch = mpatches.Patch(color='brown', label='V7')

bx.legend(handles=[red_patch,blue_patch,green_patch,yellow_patch,orange_patch,purple_patch,brown_patch], ncol=7, handlelength=0.5, handleheight=0.5, loc='center', frameon=False, bbox_to_anchor=(0.5, .8), fontsize=12, handletextpad=.35, columnspacing = 1)

labels = [item.get_text() for item in bx.get_xticklabels()]

bx.set_xticklabels(labels)

bx.set_xticks([])

bamfile = "/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix_New_030723/GenomeCov_Deletion/P2220168.depth"
bamfile2 = "/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix_New_030723/GenomeCov_Deletion/P2120135.depth"
bamfile3 = "/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix_New_030723/GenomeCov_Deletion/UW15970L.depth"

import matplotlib.pyplot as plt
import csv

x1 = []
y1 = []

x2 = []
y2 = []

x3 = []
y3 = []

with open(bamfile, 'r') as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter='\t')
    for row in tsvreader:
        x1.append(float(row[1]))
        y1.append(float(row[2]))

with open(bamfile2, 'r') as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter='\t')
    for row in tsvreader:
        x2.append(float(row[1]))
        y2.append(float(row[2]))

with open(bamfile3, 'r') as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter='\t')
    for row in tsvreader:
        x3.append(float(row[1]))
        y3.append(float(row[2]))

x1 = x1[148300:-950000]
y1 = y1[148300:-950000]
x2 = x2[148300:-950000]
y2 = y2[148300:-950000]
x3 = x3[148300:-950000]
y3 = y3[148300:-950000]

fig, axs = plt.subplots(nrows=3, sharex=True, figsize=(7, 5))

for i, ax in enumerate(axs):
    if i == 0:
        x, y = x1, y1
        x4=151315
        x5=151700
        ax.axvline(x=x4, color='red', linestyle='--', linewidth=1)
        ax.axvline(x=x5, color='red', linestyle='--', linewidth=1)
        ax.axvline(x=149519, color='black', linestyle='--', linewidth=1)
        ax.plot(x, y, color='#0E84B4')
        ax.fill_between(x, y, color='#0E84B4', alpha=0.5)
        ax.text(148420, 10, "P-22-20168 (Nichols)", fontsize=9.5, color='black',fontname='Helvetica')
        ax.fill_between([151315, 151700], [300, 300], color='gray', alpha=0.5)
        ax.tick_params(axis='both', which='major', labelsize=9)

    elif i == 1:
        x, y = x2, y2
        ax.axvline(x=149519, color='black', linestyle='--', linewidth=1)
        ax.text(148420, 3.33, "P-21-20135 (Nichols)", fontsize=9.5, color='black',fontname='Helvetica')
        ax.plot(x, y, color='#0E84B4')
        ax.fill_between(x, y, color='#0E84B4', alpha=0.5)
        ax.tick_params(axis='both', which='major', labelsize=9)

    else:

        x, y = x3, y3
        ax.plot(x, y, color='#B50A2A')
        ax.fill_between(x, y, color='#B50A2A', alpha=0.5)
        ax.text(148420, 10, "UW15970L (SS14)", fontsize=9.5, color='black',fontname='Helvetica')
        ax.set_xlabel("SS14 NC_021508 Genomic Coordinates",fontname='Helvetica')
        ax.tick_params(axis='both', which='major', labelsize=9)

    ax.set_xlim([min(x), max(x)])
    ax.spines['top'].set_visible(False)

axs[0].set_ylim([0, 300])
axs[1].set_ylim([0, 100])
axs[2].set_ylim([0, 300])

axs[0].set_xlim([148400, 152600])
axs[1].set_xlim([148400, 152600])
axs[2].set_xlim([148400, 152600])

#plt.xlabel("Genome Position")
fig.text(0.04, 0.5, 'Read Depth', va='center', rotation='vertical', fontname='Helvetica')

# Save the plots to files

fig.savefig("Cov_Graph.png", dpi=300)
bx.figure.savefig("TPRD.png", dpi=300)

fig.savefig("Cov_Graph.pdf", dpi=300)
bx.figure.savefig("TPRD.pdf", dpi=300)

fig.savefig("Cov_Graph.svg", transparent=True)
bx.figure.savefig("TPRD.svg", transparent=True)

plt.show()