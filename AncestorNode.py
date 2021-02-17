"""
Used to find sources - keeps track of cohort of individuals that share the IBD.
Authors: Kelly Finke, Michael Kourakos, Sara Mathieson
Date: 10/27/19
"""
# DONE

class AncestorNode:

    def __init__(self, indv, cohort, children = False):

        self.indv = indv

        if type(cohort) == list:
            self.cohort = list(cohort)
        else:
            self.cohort = [cohort]

        if not children:
            self.children = []
        elif type(children) == list:
            self.children = children
        else:
            self.children = [children]

        # these are for a specific IBD (we compute them later)
        self.num_paths = None
        self.prob = None

    def add_child(self, child):
        self.children.append(child)
        self.cohort.extend(child.cohort)

    def has_child(self, child_id):
        for cur_child in self.children:
            if cur_child.indv.id == child_id:
                return True
        return False

    # OVERRIDE
    def __str__(self):
        children_str = ""
        cohort_str = ""
        for child in self.children:
            children_str += child.indv.id + ", "
        for indv in self.cohort:
            cohort_str += indv.id + ", "
        return str("id: %s, children: %s cohort: %s" % (self.indv.id, \
            children_str, cohort_str))

    def __eq__(self, other):
        if isinstance(other, AncestorNode):
            return self.indv == other.indv
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.indv)
