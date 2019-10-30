"""
Methods to visualize IBD locations and groupings. NOTE: this assumes you have a
"figs" directory, and below that a directory for each chromosome. i.e. figs/chr2
Authors: Kelly Finke, Mikey Kourakos, Sara Mathieson
Date: 10/27/19
"""
# DONE

# for displaying figures over ssh
import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt

from tqdm import tqdm
from matplotlib import collections as mc

def group_conflicts_diagram(groups, conflicts, indv, CHROM, chrom_len, flag):
    """Plot each group in a different color"""
    ibd_counts = {}
    ibds = []
    orig_color_lst = [group.color for group in groups]
    additional_colors = [group.secondary_colors for group in groups]
    print("constructing diagram...")
    for group in groups: # can tqdm
        for ibd in group.ibds:
            if ibd in ibd_counts:
                ibd_counts[ibd] += 1
            else:
                ibd_counts[ibd] = 0
            ibds.append((ibd, group))
        ibds.append("next")

    colors = []
    segs = []
    lines = []
    y_val = 1
    for item in ibds:
        if item == "next":
            y_val += 20 # TODO change to 10?
        else:
            ibd, group = item
            colors.append(group.color)
            segs.append([(ibd.start, y_val), (ibd.end, y_val)])
            lines.append(1)
            for con_item in conflicts[ibd]:
                conflict, color = con_item
                if color not in orig_color_lst:
                    idx = find_color(color, additional_colors)

                    # TODO should deal with this by removing outdated conflicts,
                    # but is not actually a huge problem for now
                    if idx == None:
                        color = None
                    else:
                        color = orig_color_lst[idx]
                if color != None:
                    segs.append([(int(conflict)-100, y_val), (int(conflict)+ \
                        100, y_val)])
                    colors.append(color)
                    lines.append(10)
            y_val += 1

    # set up plot
    lc = mc.LineCollection(segs, colors=colors, linewidths=lines)
    fig, ax = plt.subplots()
    ax.add_collection(lc)
    ax.margins(0.1)
    plt.title("\""+indv+"\" Chrom " + str(CHROM) + " grouped IBD segments")
    plt.xlabel("IBD segment start and end points")
    plt.xlim(chrom_len[0], chrom_len[1])

    # save according to flag (assumes a directory i.e. figs/chr21)
    if flag == "init":
        figname = "figs/chr" + str(CHROM) + "/" +indv+"_pathLpng"
    else:
        figname = "figs/chr" + str(CHROM) + "/" + flag + "/" +indv+"_pathL.png"
    if flag != "geno":
        plt.savefig(figname, format='png', dpi=400, bbox_inches='tight', \
            pad_inches=0)

    # this is essential for not causing memory problems
    plt.close()

def find_color(color, lol):
    """Find the appropriate conflict color."""
    for i in range(len(lol)):
        if color in lol[i]:
            return i
    # groups that are removed are not getting conflicts removed (TODO fix this)
    return None
