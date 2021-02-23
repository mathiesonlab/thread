"""
Methods for determining the content of an IBD (actual SNP values), i.e. IBD
grouping algorithms.
Authors: Kelly Finke, Michael Kourakos, Sara Mathieson
Date: 10/27/19
"""
# DONE

# python imports
import json
from matplotlib import colors as mcolors
import random
import time
#from tqdm import tqdm

# our imports
from IBD_Group import Chrom
from IBD_Group import IBD_Group
import recombination_points as reco

# thresholds for combining IBDs and later groups
OVERLAP_MIN = 250
GROUP_OVERLAP = 100
CONFLICT_MAX = 5
MIN_LEN = 300000
MIN_SNPS = 100
MIN_IBDS = 5

def pairwise_conflicts(group, ibd, ibd_alleles):
    """
    Compares bases in ibd (from json file) to bases in template (passed in as a
    dictionary). Returns the overlap and conflicts between template and ibd.
    """

    overlap = 0
    conflicts = {}
    template = group.template
    for base in ibd_alleles:
        if base in template:
            if(template[base] == ibd_alleles[base]):
                overlap += 1
            else:
                conflicts[base] = (template[base], ibd_alleles[base])

    return (overlap, conflicts)


def pairwise_group_conflicts(group1, group2,pairwise_groups_data):
    """
    Compares bases in two groups (from json file). Returns the number of
    overlapping SNPs (overlap) and the NUMBER of conflicts. Right now this does
    not add conflicts in colors to IBDs within the groups, but it could.
    """
    # check in conflict object first (reduce pairwise group conflicts calls)

    group1_id, group1_version = group1.get_group_template_id()
    group2_id, group2_version = group2.get_group_template_id()

    if group1_id > group2_id:
        groups_id = (group1_id, group2_id)
        groups_version = (group1_version, group2_version)
    else:
        groups_id = (group2_id, group1_id)
        groups_version = (group2_version, group1_version)

    # if the number of conflicts hasn't been counted for this group, or
    # if the dictionary store the conflict for the previous versions of the group's template
    if not groups_id in pairwise_groups_data or pairwise_groups_data[groups_id][0] != groups_version:
        overlap = 0
        conflicts = 0

        for base in group1.template:
            if base in group2.template:
                if group1.template[base] == group2.template[base]:
                    overlap += 1
                else:
                    conflicts += 1

        pairwise_groups_data[groups_id] = (groups_version, (overlap, conflicts))

    return pairwise_groups_data[groups_id][1] # return (overlap, conflicts) part of the value

def group_IBDs(ibds, indv, ibd_dict,pairwise_groups_data):
    """
    Groups the given IBDs into groups based on overlaps.
    Start with one group (two if we already have haploid data) with first IBD
    for IBDs in overlap, add IBDs until number of conflicts reaches threshold
    if conflicts between new ibd and group > threshold, start new group with new
    IBDs for all subsequent IBDs, check against both groups, etc
    """

    # sort by length (longest to shortest)
    ibds = sorted(list(ibds), key=lambda ibd: len(ibd), reverse=True)

    # list of IBD_Group objects
    groups = []
    all_conflicts = {}

    # initialize group colors
    group_colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    group_colors = list(group_colors.keys())

    # find stretches of homozygosity and keep homozygous groups separate
    homozygous = find_homozygous(ibds, ibd_dict)
    homoz_groups, remaining_ibds = make_homoz_groups(ibds, ibd_dict, \
                                                     homozygous, all_conflicts, group_colors)

    # create initial heterozygous groups based on IBDs
    steady_state = False
    passes = 0
    while not steady_state:
        before = len(remaining_ibds)
        remaining_ibds = group_IBDs_helper(remaining_ibds, ibd_dict, \
                                           homoz_groups, groups, all_conflicts, group_colors, True, pairwise_groups_data)
        passes += 1

        # stop when we are no longer incorporating more IBDs
        if len(remaining_ibds) == before:
            steady_state = True

    # using remaining IBDs to merge groups
    all_groups, remaining_ibds = merge_groups(remaining_ibds, ibd_dict, homoz_groups, groups, all_conflicts,
                                              group_colors, indv, pairwise_groups_data)

    # try once more to add IBDs
    remaining_ibds = group_IBDs_helper(remaining_ibds, ibd_dict, [], \
                                       all_groups, all_conflicts, group_colors, False, pairwise_groups_data)

    # final pass to merge groups (not using IBDs, based on group overlaps)
    steady_state = False
    while not steady_state:
        prev_num = len(all_groups)
        i = 0
        # look at pairs of groups
        for group1 in all_groups:

            # find potential groups to merge with group1
            potentials = []
            for j in range(i+1, len(all_groups)):
                group2 = all_groups[j]
                has_possible_overlap = group1.start <= group2.end and group1.end >= group2.start
                if has_possible_overlap:
                    overlap,conflicts = pairwise_group_conflicts(group1, group2, pairwise_groups_data)
                else:
                    overlap, conflicts = 0,0
                if overlap > GROUP_OVERLAP and conflicts < CONFLICT_MAX:
                    potentials.append([group2, overlap, j])

            if len(potentials) >= 1:
                max_overlap = max(potentials, key=lambda item: item[1])
                result = group1.combine_with(max_overlap[0])
                if result:
                    idx = max_overlap[2]
                    all_groups.pop(idx)

            i += 1
        if prev_num == len(all_groups):
            steady_state = True

    # final step to merge non-overlapping groups into haplotypes
    all_groups = condense_groups(all_groups, pairwise_groups_data)

    # hopefully we have used up all IBDs by this point
    if len(remaining_ibds) > 0:
        print('HUGE PROBLEM: IBDs REMAIN!')

    # return all groups and conflicts between them
    return all_groups, all_conflicts

def merge_groups(remaining_ibds, ibd_dict, homoz_groups, groups, all_conflicts,group_colors, indv,
                 pairwise_groups_data):
    """Use remaining IBDs (that hopefully overlap 2 groups) to merge groups
    """

    # first duplicate homoz groups (need two of them)
    all_groups = homoz_groups + groups
    for homoz_g in homoz_groups:
        ibd = random.sample(homoz_g.ibds, 1)[0]
        ibd_alleles = ibd_dict[str(ibd)]
        new_g = homoz_g.duplicate(ibd, ibd_alleles, group_colors.pop(0))
        all_groups.append(new_g)

    passes = 0
    steady_state = False
    while len(remaining_ibds) > 0 and not steady_state:
        remaining = []
        remaining_ibds = sorted(remaining_ibds, key=lambda ibd: len(ibd), \
                                reverse=True)

        for ibd in remaining_ibds:

            ibd_alleles = ibd_dict[str(ibd)]
            if ibd not in all_conflicts.keys():
                all_conflicts[ibd] = []

            matching_groups = []

            for group in all_groups:
                has_possible_overlap = group.start <= ibd.end and group.end >= ibd.start
                if has_possible_overlap:
                    overlap, conflicts = pairwise_conflicts(group, ibd, ibd_alleles)
                else:
                    overlap, conflicts = 0,{}
                for conflict in conflicts.keys():
                    all_conflicts[ibd].append((conflict, group.color))

                if len(conflicts.keys()) < CONFLICT_MAX:
                    if overlap > OVERLAP_MIN:
                        matching_groups.append(group)

            placed = False
            if len(matching_groups) == 1:
                matching_groups[0].add(ibd, ibd_alleles)
                placed = True

            # case for 3 (often 2 could be homozygous)
            if passes >= 1 and len(matching_groups) == 3:
                if matching_groups[0] == matching_groups[1]:
                    matching_groups = [matching_groups[0], matching_groups[2]]
                elif matching_groups[0] == matching_groups[2]:
                    matching_groups = [matching_groups[0], matching_groups[1]]
                elif matching_groups[1] == matching_groups[2]:
                    matching_groups = [matching_groups[0], matching_groups[1]]

            if passes >= 1 and len(matching_groups) == 2 and matching_groups[0]\
                    != matching_groups[1]:
                # TODO don't merge if we have too many conflicts
                has_possible_overlap = group1.start <= group2.end and group1.end >= group2.start
                if has_possible_overlap:
                    overlap,conflicts = pairwise_group_conflicts( \
                    matching_groups[0], matching_groups[1], pairwise_groups_data)
                else:
                    overlap, conflicts = 0,0

                # overlap may not be that high since IBD might span
                if conflicts < CONFLICT_MAX:

                    # do this before adding IBD so we don't accidentially merge
                    # two homozygous groups
                    result = matching_groups[0].combine_with(matching_groups[1])
                    if result:
                        matching_groups[0].add(ibd, ibd_alleles)
                        # duplicated groups (i.e. homoz) will be equal
                        all_groups.remove(matching_groups[1])
                        placed = True

            if not placed:
                remaining.append(ibd)

        passes += 1
        if len(remaining_ibds) == len(remaining):
            steady_state = True
        remaining_ibds = remaining

    return all_groups, remaining_ibds

def initial_groups(ibds, indv, ibd_dict):
    """
    If already genotyped, we can use the divided haplotypes to initialize the
    groups.
    """
    all_conflicts = {}
    group_colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    group_colors = list(group_colors.keys())

    # initialize group1 and group2
    groups = [None, None]
    for ibd in ibds:
        # do this only once!
        all_conflicts[ibd] = []
        ibd_alleles = ibd_dict[str(ibd)]
        hap = ibd.get_indvs()[indv]

        # start groups
        if groups[0] == None and (hap == 0 or hap == 1):
            groups[0] = IBD_Group(ibd, ibd_alleles, group_colors.pop(0))
        if groups[1] == None and (hap == 0 or hap == 2):
            groups[1] = IBD_Group(ibd, ibd_alleles, group_colors.pop(0))

        for group in groups:
            if group != None:
                has_possible_overlap = group.start <= ibd.end and group.end >= ibd.start
                if has_possible_overlap:
                    overlap, conflicts = pairwise_conflicts(group, ibd, ibd_alleles)
                else:
                    overlap, conflicts = 0,{}
                for conflict in conflicts.keys():
                    all_conflicts[ibd].append((conflict, group.color))
        if hap == 0 or hap == 1:
            if groups[0] != None:
                groups[0].add(ibd, ibd_alleles)
        if hap == 0 or hap == 2:
            if groups[1] != None:
                groups[1].add(ibd, ibd_alleles)

    # in rare circumstances the entire chrom is homozygous
    # TODO we should probably duplicate in these circumstances
    if groups[0] == None:
        groups = [groups[1]]
    elif groups[1] == None:
        groups = [groups[0]]

    return groups, all_conflicts

def group_IBDs_helper(ibds, ibd_dict, homoz_groups, groups, all_conflicts, \
                      group_colors, ignore_homoz, pairwise_groups_data):
    """
    Runs one pass of trying to add IBDs to groups, saving those it hasn't added
    yet. Continue running this function until all IBDs are added. IBDs to be
    careful about: those that overlap with homozygous regions (these could create
    false haplotypes: wrong combination) and those that can be added to multiple groups.
    """
    remaining_ibds = []
    ibds = sorted(ibds, key=lambda ibd: len(ibd), reverse=True)
    for ibd in ibds: #can tqdm

        ibd_alleles = ibd_dict[str(ibd)]
        if ibd not in all_conflicts.keys():
            all_conflicts[ibd] = []

        if ignore_homoz and overlap_homoz(ibd, homoz_groups):
            remaining_ibds.append(ibd)

        else:
            placed = False
            # create first group if none exist
            if len(groups) == 0:
                # first group
                new_group = IBD_Group(ibd, ibd_alleles, group_colors.pop(0))
                groups.append(new_group)
                placed = True
            else:
                # overlap and conflicts between IBD and each group
                matches = 0
                potentials = []
                for group in groups:
                    has_possible_overlap = group.start <= ibd.end and group.end >= ibd.start
                    if has_possible_overlap:
                        overlap, conflicts = pairwise_conflicts(group, ibd, ibd_alleles)
                    else:
                        overlap, conflicts = 0,{}

                    for conflict in conflicts.keys():
                        all_conflicts[ibd].append((conflict, group.color))

                    if len(conflicts.keys()) < CONFLICT_MAX:
                        if overlap > OVERLAP_MIN:
                            matches += 1
                            potentials.append([group, overlap])
                            placed = True

                        # reminder that there was insufficient overlap
                        else:
                            pass

                # check if no matching groups found, create new group if so
                if matches == 0:
                    # TODO don't make a new group if IBD was just too short
                    # hack to reset group colors
                    # TODO: after merge, return extra group color to the pool
                    if len(group_colors) == 0:
                        group_colors = dict(mcolors.BASE_COLORS, \
                                            **mcolors.CSS4_COLORS)
                        group_colors = list(group_colors.keys())

                    new_group = IBD_Group(ibd, ibd_alleles, group_colors.pop(0))
                    groups.append(new_group)
                    placed = True

                # right now, if multiple potentials, add to one with max overlap
                elif matches >= 1:
                    max_overlap = max(potentials, key=lambda item: item[1])
                    max_overlap[0].add(ibd, ibd_alleles)

            if not placed:
                print("didn't place for some reason???")
                remaining_ibds.append(ibd)

    return remaining_ibds

def overlap_homoz(ibd, homoz_groups):
    for homoz_group in homoz_groups:
        start = homoz_group.start
        end = homoz_group.end
        # straddle start or end of homoz group
        if (ibd.start <= start <= ibd.end) or (ibd.start <= end <= ibd.end):
            return True
    return False

def find_homozygous(ibds, ibd_dict):
    """ibds is list of IBDs in this individual, ibd_dict contains info about
    the sequences of the ibds. This function finds stretches of homozygosity,
    i.e. all no IBDs in a long region have any conflicts
    """

    # first create one large Chrom with all alleles
    ibd_alleles = ibd_dict[str(ibds[0])]
    whole_chrom = Chrom(ibds[0], ibd_alleles)
    for ibd in ibds:
        ibd_alleles = ibd_dict[str(ibd)]
        whole_chrom.add(ibd, ibd_alleles)

    # then see what regions of this chromosome have no conflicts
    homoz_stretches = []
    stretch = []
    # note json must store as string, so convert to int just for this step!
    for base in sorted([int(b) for b in whole_chrom.template.keys()]):

        if len(whole_chrom.template[str(base)]) == 1: # homozygous
            stretch.append(base)
        else:
            # add to list of stretches and restart stretch
            homoz_stretches.append(stretch)
            stretch = []

    # finally, return long stretches
    long_stretches = []
    for stretch in homoz_stretches:
        if len(stretch) > 0:
            first_base = int(stretch[0])
            last_base = int(stretch[-1])

            # difference between first and last base
            length = last_base - first_base
            if length >= MIN_LEN and len(stretch) >= MIN_SNPS:
                long_stretches.append([first_base, last_base])

    return long_stretches

def make_homoz_groups(ibds, ibd_dict, homozygous, all_conflicts, group_colors):
    """Create homozygous groups"""
    homoz_groups = [None]*len(homozygous)
    remaining_ibds = []
    for ibd in ibds:
        placed = False

        # see if IBD belongs to any homozygous stretch
        for i in range(len(homozygous)):
            stretch = homozygous[i]
            start = stretch[0]
            end = stretch[1]

            # see if IBD inside stretch
            if ibd.start >= start and ibd.end <= end:
                if homoz_groups[i] == None:
                    homoz_groups[i] = [ibd]
                else:
                    homoz_groups[i].append(ibd)
                placed = True
        if not placed:
            remaining_ibds.append(ibd)

    # if no IBDs assigned, abandon homoz stretch
    filtered = []
    for g in homoz_groups:
        # case when we have enough IBDs
        if g != None and len(g) > MIN_IBDS:
            ibd = g[0]
            ibd_alleles = ibd_dict[str(ibd)]
            new_group = IBD_Group(ibd, ibd_alleles, group_colors.pop(0))
            for ibd in g:
                ibd_alleles = ibd_dict[str(ibd)]
                new_group.add(ibd, ibd_alleles)
                all_conflicts[ibd] = [] # no conflicts
            filtered.append(new_group)

        # if we have just a few, return them to remaining
        elif g != None:
            for ibd in g:
                remaining_ibds.append(ibd)

    # make sure to set homoz_stretches
    for g in filtered:
        g.homoz_stretches = [[g.start, g.end]]
    return filtered, remaining_ibds

def condense_groups(all_groups, pairwise_groups_data):
    """
    Combines groups that have no conflicts. Important that we only combine one
    pair per iteration!
    """
    # sort by length first, then we'll combine long groups first
    all_groups = sorted(all_groups, key=lambda group: len(group), reverse=True)

    steady_state = False
    while not steady_state:
        num_groups = len(all_groups)
        combined = False
        i = 0

        for group1 in all_groups:
            potentials = []
            for j in range(i+1, len(all_groups)):
                group2 = all_groups[j]
                has_possible_overlap = group1.start <= group2.end and group1.end >= group2.start
                if has_possible_overlap:
                    overlap,conflicts = pairwise_group_conflicts(group1, group2, pairwise_groups_data)
                else:
                    overlap, conflicts = 0,0
                if conflicts == 0:
                    potentials.append(group2)

            # combine with first (longest) group
            while len(potentials) > 0 and not combined:
                best_group = best_combine(group1, potentials)
                # special case at the end with 3 groups
                if len(all_groups) == 3:
                    override = True
                else:
                    override = False

                # try to merge
                result = group1.combine_with(best_group, override)
                # success
                if result:
                    all_groups.remove(best_group)
                    combined = True
                # if bad merge, remove and try again (made while loop)
                elif len(potentials) > 0:
                    potentials.remove(best_group) # remove bad try

            i += 1

        # check if number of groups has stabilized
        if num_groups == len(all_groups):
            steady_state = True

    return all_groups

def best_combine(group1, potentials):
    """
    From a list of potential groups (no conflicts), find the best one to
    combine with the given group. This could be based on the groups "lining up"
    one right after the other.
    """
    dist_group = []
    for group2 in potentials:
        if group1.end <= group2.start:
            distance = group2.start - group1.end
        elif group2.end <= group1.start:
            distance = group1.start - group2.end
        elif group1.start <= group2.start <= group1.end:
            distance = group2.start - group1.end # negative distance (overlap)
        elif group2.start <= group1.start <= group2.end:
            distance = group1.start - group2.end # negative distance (overlap)
        else:
            print('HUGE PROBLEM: not disjoint or overlapping!')
        dist_group.append([distance, group2])

    # return group with smallest distance
    return min(dist_group, key=lambda item: item[0])[1]

def add_IBDs(id, groups, ibds, ibd_dict, all_conflicts):
    """
    For a "good" individual where we already have two strong groups, add in new
    IBDs (should not conflict because we should have checked this when rooting).
    """
    if len(groups) != 2:
        print("HUGE PROBLEM! reconstructed indv", id, "with", len(groups), \
              "groups.")

    for ibd in ibds:

        # haven't added it yet
        if not ibd_in_groups(ibd, groups):
            matches = 0
            ibd_alleles = ibd_dict[str(ibd)]
            if ibd not in all_conflicts.keys():
                all_conflicts[ibd] = []

            for group in groups:
                has_possible_overlap = group.start <= ibd.end and group.end >= ibd.start
                if has_possible_overlap:
                    overlap, conflicts = pairwise_conflicts(group, ibd, ibd_alleles)
                else:
                    overlap, conflicts = 0,{}
                for conflict in conflicts.keys():
                    all_conflicts[ibd].append((conflict, group.color))
                if len(conflicts.keys()) < CONFLICT_MAX:
                    if overlap > OVERLAP_MIN:
                        matches += 1
                        group.add(ibd, ibd_alleles)

            if matches != 1:
                # TODO what to do in this case?
                #print("Unexpected number of matches!", matches)
                pass

def ibd_in_groups(ibd, groups):
    """Return True if IBD is in at least one group, False otherwise"""
    for group in groups:
        if ibd in group.ibds:
            return True
    return False
