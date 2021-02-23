"""
IBD Group object: wrapper for a list of IBDs, used in grouping algorithms
Authors: Kelly Finke, Michael Kourakos, Sara Mathieson
Date: 10/27/19
"""
# DONE
# python import
from uuid import uuid4

class IBD_Group:
    """Container object for a group of IBDs with a common sequence"""

    def __init__(self, ibd, ibd_alleles, color):
        """
        Constructor takes in a single IBD as the initial template of the group
        """
        self.ibds = {ibd}
        self.template = {}

        # initialize template with values from starting ibd
        for base, allele in ibd_alleles.items():
            self.template[base] = allele

        self.start = ibd.start
        self.end = ibd.end
        self.color = color
        self.secondary_colors = []
        self.homoz_stretches = [] # parts of the group that are homozygous

        # generate an unique id for each group and keep track of the 
        # group's latest version 
        self.id = uuid4()
        self.template_last_modified = 0
    # GETTERS
    def get_allele(self, base):
        """If base in template, return allele; else return ?"""
        if str(base) in self.template: # keys are str (something about json...)
            return self.template[str(base)]
        else:
            return "?"

    # OVERRIDE
    def __len__(self):
        return self.end - self.start

    def __str__(self):
        """Create string representation"""
        return ""+str(len(self.ibds))+" ibds, start: "+str(self.start)+ \
            ", end: "+str(self.end)

    def __lt__(self, other):
        if not isinstance(other, IBD_Group):
            raise TypeError("Cannot compare IBD object to anything but an IBD")
        return self.id < other.id

    def __eq__(self, other):
        if isinstance(other, IBD_Group):
            return (self.start == other.start) and (self.end == other.end) and \
                (self.template == other.template)
        return False

    # SETTERS
    def add(self, ibd, ibd_alleles):
        """Add an IBD to this group"""
        self.ibds.add(ibd)
        # TODO not currently checking that this does/doesn't conflict with
        # existing base
        for base, allele in ibd_alleles.items():
            self.template[base] = allele
       # as template is modified
        self.template_last_modified +=1

        self.start = min(self.start, ibd.start)
        self.end = max(self.end, ibd.end)

    # GETTERS
    def get_group_template_id(self):
        return self.id, self.template_last_modified

    def combine_with(self, group, override=False):
        """
        Combine two groups (need to do lots of checks to make sure this is okay)
        Returns True if they could be combined, False otherwise
        """
        # make sure we are not combining identical groups!
        if self == group:
            return False

        # if they share a homozygous group and we're not overriding, NO
        if self.share_homoz(group) and override == False:
            return False

        # now combine
        for ibd in group.ibds:
            self.ibds.add(ibd)
        for base, allele in group.template.items():
            # TODO should we be checking for conflicts somehow?
            self.template[base] = allele

        # after merging this group with the other groups template
        self.template_last_modified +=1

        # modify start and end, homoz_stretches, and colors
        self.start = min(self.start, group.start)
        self.end = max(self.end, group.end)
        for s in group.homoz_stretches:
            self.homoz_stretches.append([s[0], s[1]])
        self.secondary_colors.append(group.color)
        self.secondary_colors.extend(group.secondary_colors)
        return True

    def share_homoz(self, group):
        """Check if these two groups share a homozygous stretch"""
        for s1 in self.homoz_stretches:
            for s2 in group.homoz_stretches:
                if s1[0] == s2[0] and s1[1] == s2[1]:
                    return True
        return False

    def duplicate(self, ibd, ibd_alleles, color):
        """Duplicate so we have two of each homozygous group"""
        #ibd = self.ibds[0] # TODO can we do this?
        new = IBD_Group(ibd, ibd_alleles, color)
        for ibd in self.ibds:
            new.ibds.add(ibd)
        for base, allele in self.template.items():
            new.template[base] = allele
        new.start = self.start
        new.end = self.end

        new.homoz_stretches = []
        for s in self.homoz_stretches:
            new.homoz_stretches.append([s[0], s[1]]) # deep copy

        new.seconary_colors = []
        for c in self.secondary_colors: # deep copy
            new.seconary_colors.append(c)
        return new

class Chrom:
    """Container object for a chromosome and its IBDs (allows conflicts)"""

    def __init__(self, ibd, ibd_alleles):
        """
        Constructor takes in a single IBD as the initial template of the chrom
        """
        self.ibds = {ibd} # set
        self.template = {} # dict

        # initialize template with values from starting ibd
        for base, allele in ibd_alleles.items():
            self.template[base] = [allele] # allow multiple alleles

        self.start = ibd.start
        self.end = ibd.end

    def add(self, ibd, ibd_alleles):
        self.ibds.add(ibd)
        #TODO make faster with lists
        for base, allele in ibd_alleles.items():
            #assert int(base) == base
            if base in self.template:
                if allele not in self.template[base]:
                    self.template[base].append(allele)
            else:
                self.template[base] = [allele]
        self.start = min(self.start, ibd.start)
        self.end = max(self.end, ibd.end)
