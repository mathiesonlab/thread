"""
Individual object: contains individual identifiers and IBD segments
Authors: Kelly Finke, Michael Kourakos, Sara Mathieson
Date: 8/13/19
"""
# DONE

class Individual:
    """
    Individual represents an individual member of the pedigree and contains
    links to immeditately related individuals and a set of IBD objects.
    """

    def __init__(self, id, father, mother, sex): # all strs
        self.id = id # unique individual id (str)
        self.m_id = mother # str id of mother
        self.p_id = father # str id of father
        self.p = None # Indivdual obj containing father
        self.m = None # Individual obj containing mother
        self.sex = int(sex) # 1 or 2 (changing to int to be consistent with hap)

        # Couple obj with both parents (None if Indv is founder or married-in)
        self.parents = None
        # True if Indv is a founder (indv and spouse have no parents in the ped)
        self.founder = False
        # True if Indv is married in (has no parents in the pedigree)
        self.married_in = False

        # relationships
        self.couples = [] # Couple object(s) containing indivdual and spouse(s)
        self.children = []
        self.checked = False # used in find_relations method (PedigreeTree)

        # probability individual has an IBD (vector of length 3: [p0, p1, p2])
        # Note: this gets overriden many times as we try out this individual
        # having many IBDs should not be used later on
        self._pibd = None

        # sequence information
        self.genotyped = False # True if we have real genome data from Indv
        self._IBDs = {"00": {}, # unknown haplotype, uncertain
                    "01": set(),# unknown haplotype, certain
                    "10": {},   # paternal haplotype, uncertain
                    "11": set(),# paternal haplotype, certain
                    "20": {},   # maternal haplotype, uncertain
                    "21": set(),# maternal haplotype, certain
                    "30": {},   # undetermined parent haplotype 1, uncertain
                    "31": set(),# undetermined parent haplotype 1, certain
                    "40": {},   # undetermined parent haplotype 2, uncertain
                    "41": set() # undetermined parent haplotype 2, certain
                }

    # OVERRIDE
    def __str__(self):
        return str("id: %s, p: %s, m: %s" %(self.id, self.p_id, self.m_id))

    def __eq__(self, other):
        if isinstance(other, Individual):
            return self.id == other.id and self.m_id == other.m_id and \
                self.p_id == other.p_id
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.id)

    # GETTERS
    def get_pibd(self): return self._pibd

    def get_IBDs(self, key=None):
        """Return IBDs of a particular type (or all if key is None)"""
        if key == None:
            # only returns "certain" IBDs
            return self._IBDs['01'] | self._IBDs['11'] | self._IBDs['21'] | \
                self._IBDs['31'] | self._IBDs['41']
        else:
            return self._IBDs[key]

    def descendants(self):
        """Return a set of descendants of this individual (recursive)"""
        descend = set(self.children)
        for c in self.children:
            descend |= c.descendants()
        return descend

    # SETTERS
    def set_pibd(self, pibd):
        assert len(pibd) == 3
        assert abs(sum(pibd) - 1) < 1e-6
        self._pibd = pibd

    def reset_pibd(self):
        self._pibd = None

    def add_ibd(self, key, ibd):
        """Add certain IBDs only for right now"""
        assert key in self._IBDs.keys()
        self._IBDs[key].add(ibd)

    def move_IBDs(self):
        """This moves 1/2 IBDs to 3/4 (since we can't resolve parents)"""
        self._IBDs["31"] = self._IBDs["11"]
        self._IBDs["41"] =  self._IBDs["21"]
        self._IBDs["11"] = set()
        self._IBDs["21"] = set()

    def swap_IBDs(self):
        """Swap 1&2 IBDs since we can resolve parental haplotypes"""
        temp_11 = self._IBDs["11"]
        self._IBDs["11"] = self._IBDs["21"]
        self._IBDs["21"] = temp_11

    def remove_ibd(self, ibd):
        """If IBD is incorrectly placed, we can use this method to remove it"""
        removed = False
        for key in self._IBDs:
            if ibd in self._IBDs[key]:
                self._IBDs[key].remove(ibd)
                removed = True
        if not removed:
            print('HUGE PROBLEM! IBD NOT IN INDIVIDUAL!', self.id, str(ibd))
        return removed
