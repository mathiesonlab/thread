"""
IBD class, as well as functions to give IBDs to individuals.
Authors: Kelly Finke, Michael Kourakos, Sara Mathieson
Date: 10/27/19
"""
# DONE

# python imports
import json
from tqdm import tqdm

class IBD:
    """
    A class to hold all the information about a certain IBD segment. Keeps track
    of the chromosome it is located on, the start and end points, the number of
    SNP's within the IBD segment, and a dictionary of individuals who have the
    segment. The dictionary has Individual IDs as keys, and a 0, 1, or 2 as each
    value (int). A 0 indicates homozygousity, and 1 or 2 denote the haplotype
    the IBD was found on.
    """
    def __init__(self, chromosome, start, end, total_SNPs, genetic_dist, \
        indv = None):
        # should not change
        self.chromosome = int(chromosome)
        self.start = int(start)
        self.end = int(end)
        self.SNPs = int(total_SNPs)
        self.genetic_distance = float(genetic_dist)
        self.id_indv = indv
        self.id = self.get_id(self.id_indv)

        # could change
        self._indvs = {}    # indv_id   -> haplotype number, 0 = both haplotypes
        self._sources = {}  # source_id -> ancestors
        self._curr_source = None # currently assigned source (id which is a str)
        self._vote_plus = 0 # positive support for current source
        self._vote_minus = 0 # negative support for current source

    # GETTERS
    def get_id(self, id_indv):
        """Note: this is used in str, eq/neq, and hash"""
        return str(self.chromosome) + " " + str(self.start) + " " + \
            str(self.end) + " " + ("skeleton" if id_indv == None else id_indv)

    def get_indvs(self): return self._indvs
    def get_sources(self): return self._sources
    def get_curr_source(self): return self._curr_source
    def get_vote_plus(self): return self._vote_plus
    def get_vote_minus(self): return self._vote_minus

    def get_hap(self, indv_id):
        assert indv_id in self._indvs
        return self._indvs[indv_id]

    # SETTERS
    def set_hap(self, indv_id, hap):
        """Make sure ID is a str and hap is an int"""
        assert str(indv_id) == indv_id
        assert int(hap) == hap
        self._indvs[indv_id] = hap

    def set_sources(self, ancestors):
        self._sources = ancestors

    def remove_source(self, source):
        # TODO should we consolidate remove_source and change_source?
        del self._sources[source]

    def change_source(self, new_source):
        """
        If the current source of this IBD creates a conflict, remove and
        update the source (and update votes).
        """
        # make sure we removed previous source
        assert self._curr_source not in self._sources
        self._curr_source = new_source
        # reset votes
        self._vote_plus = 0
        self._vote_minus = 0

    def remove_indv(self, indv_id):
        """Remove an individual from this IBD based on their id"""
        del self._indvs[indv_id]

    def add_vote_plus(self):
        self._vote_plus += 1

    def add_vote_minus(self):
        self._vote_minus += 1

    # OVERRIDE
    def __len__(self):
        """Length of the IBD segment"""
        return self.end - self.start

    def __str__(self):
        return self.id

    def __eq__(self, other):
        if isinstance(other, IBD):
            return (self.id == other.id)
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.id)

def get_IBDs(germ_filename, left_out):
    """
    From each GERMLINE file, read IBDs and return a list of all IBDs. Can
    optionally leave out some individuals.
    """
    # TODO no longer using multiple GERMLINE files at once, can remove function
    #IBDs = []
    #for germ_file in germ_files: # can tqdm
    return read_germline(germ_filename, left_out)

def read_germline(germ_file, left_out):
    """
    Reads a germline .match file, and creates a list of IBD instances.
    """
    IBDs = []
    # IBD -> lists of pairs of individual ids (ie 52.1) sharing the IBDs
    IBD_set = {}

    g_file = open(germ_file, "r")
    for line in g_file:
        tokens = line.strip().split()
        # "skeleton" IBD, contains all details but no individuals, since
        # multiple IBDs can have the same skeleton
        chrom = int(tokens[4])
        start = int(tokens[5])
        end = int(tokens[6])
        num_snps = int(tokens[9])
        genetic_dist = float(tokens[10])
        skeleton = IBD(chrom, start, end, num_snps, genetic_dist)

        # find pair of individuals
        indv1 = tokens[1] # str
        indv2 = tokens[3] # str

        if skeleton not in IBD_set:
            IBD_set[skeleton] = [[indv1, indv2]]
        else:
            IBD_set[skeleton].append([indv1, indv2])

    # go through each IBD and merge pairs into a group
    for skeleton, indv_pairs in IBD_set.items():

        # we could have multiple IBD segments with same start and end point,
        # but since we're tracking indvs by haplotype, shouldn't be a problem
        groups = []

        # first go through all pairs
        for pair in indv_pairs:
            joined = []
            joined_indicies = []
            for i,group in enumerate(groups):
                if pair[0] in group or pair[1] in group:
                    group.add(pair[0])
                    group.add(pair[1])
                    joined.append(group)
                    joined_indicies.append(i)

            # if we have joined 2 or more groups with this pair
            if len(joined) > 1:
                new_group = set().union(*joined)
                groups.append(new_group)
                joined_indicies = reversed(joined_indicies)

                # remove old groups
                for i in joined_indicies:
                    groups.pop(i)

            # if no groups joined, create a new one from this pair
            elif len(joined) == 0:
                groups.append({pair[0], pair[1]})

        # create an IBD instance for each group
        for group in groups:
            # this last argument is consistent with json file
            ibd = IBD(skeleton.chromosome, skeleton.start, skeleton.end, \
                skeleton.SNPs, skeleton.genetic_distance, str(min(group)))

            for indv_long in group:
                hap = int(indv_long[-1]) + 1 # gets a 0 or 1, we want 1 or 2
                indv = indv_long[:-2] # str

                # leave one out approach for testing
                if indv not in left_out:
                    # if both haplotypes are present, represent with a 0
                    if indv in ibd.get_indvs():
                        if ibd.get_hap(indv) != 0: # prob not necessary
                            if ibd.get_hap(indv) != hap:
                                ibd.set_hap(indv, 0)
                    else:
                        ibd.set_hap(indv, hap)
            IBDs.append(ibd)
    return IBDs

def ibd_to_indvs(ibds, ped):
    """Assign individuals IBDs (not currently accounting for homozygous)"""

    # give individuals ibds
    for ibd in ibds:
        # in the case of homoz, we can give IBD to both parents as well
        #to_add = set() SM: removing this right now

        # go through all individuals that have this IBD
        for indv, hap in ibd.get_indvs().items():
            if indv in ped.indvs: # adding check (for pedigree subsets)
                indv = ped.indvs[indv]

                # homozygous
                if str(hap) == "0":
                    indv.add_ibd("11", ibd)
                    indv.add_ibd("21", ibd)

                    # in this case we know the parents have it too
                    #indv.p.add_ibd("01", ibd) # SM: removing for now
                    #indv.m.add_ibd("01", ibd)
                    #to_add.add(indv.p.id)
                    #to_add.add(indv.m.id)

                # heterozygous
                else:
                    indv.add_ibd(str(hap) + "1", ibd)

                # set genotyped
                indv.genotyped = True
                ped.genotyped.add(indv)

        # special case (-1) when we don't know hap but we know parent had IBD
        #for parent_id in to_add: #SM: removing for now
        #    ibd.set_hap(parent_id, -1)

    # separate ibds (in cases when parents are genotyped), not super necessary
    separate_ibds(ped)

def separate_ibds(ped):
    """Assign IBD lists to proper parents based on number of shared IBDs"""
    THRESHOLD = 2

    # go through each individual, one at a time
    for indv in ped.genotyped:
        if indv == "0":
            continue

        if indv.married_in or indv.founder: # we can't know parents' haplotypes
            indv.move_IBDs() # move 1/2 to 3/4
            continue

        m = indv.m
        p = indv.p
        # if parents are genotyped, check relatedness to parents
        if p.genotyped == True and m.genotyped == True:

            p_ibds = p.get_IBDs("01") | p.get_IBDs("11") | p.get_IBDs("21")
            shared_p_1 = len(p_ibds & indv.get_IBDs("11"))
            shared_p_2 = len(p_ibds & indv.get_IBDs("21"))

            m_ibds = m.get_IBDs("01") | m.get_IBDs("11") | m.get_IBDs("21")
            shared_m_1 = len(m_ibds & indv.get_IBDs("11"))
            shared_m_2 = len(m_ibds & indv.get_IBDs("21"))

            # ibds["11"] are paternal and ibds["21"] are maternal, do nothing
            if shared_p_1 > shared_p_2 and shared_m_1 < shared_m_2:
                continue

            # ibds["21"] are paternal & ibds["11"] are maternal, swap ibd lists
            elif shared_p_1 < shared_p_2 and shared_m_1 > shared_m_2:
                indv.swap_IBDs() # swap 1&2
                continue

            # haplotype parents inconclusive
            else:

                p_dif = shared_p_1 - shared_p_2
                m_dif = shared_m_1 - shared_m_2

                if p_dif > THRESHOLD * m_dif:
                    continue
                elif m_dif > THRESHOLD * p_dif:
                    indv.swap_IBDs() # swap 1&2
                    continue
                else:
                    # still inconclusive
                    indv.move_IBDs() # 1/2 to 3/4
