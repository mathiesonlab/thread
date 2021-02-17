"""
Pedigree Tree Object: holds pedigree tree data structure constructed from file
of family relationships.

    file format:

    ID FATHER MOTHER SEX            (<- header line (ignored in code))
    id father's_id mother's_id 1    (<-sex: 1 = male, 2 = female)
    ...                             (<- continue for all individuals)

    NOTE: any id listed as a father or mother must have its own entry as well
    NOTE: if father and/or mother does not exit in pedigree, input zeroes
    for the id of the missing parent. DO NOT use "0" as the id of any indivudal
    as it is reserved for parents of married-in indivudals.

Authors: Kelly Finke, Michael Kourakos, Sara Mathieson
Date: 10/27/19
"""
# DONE

# our imports
from Individual import Individual
from Couple import Couple
from AncestorNode import AncestorNode
import IBD

class PedigreeTree:
    """
    PedigreeTree contains an undirected graph that represents a family tree
    """

    def __init__(self, ped_filename):
        ped_file = open(ped_filename, 'r')
        ped_data = ped_file.readlines()
        ped_file.close()

        self.founders = []
        self.married_in = []
        self.genotyped = set()
        self.indvs = self.construct_individuals(ped_data[1:])
        self.find_relations()

    def construct_individuals(self, ped_data):
        """
        Reads pedigree file and creates individuals
        note: the "0" id is reserved for parents of married-in individuals
        """
        indvs = {} # str id -> Individual
        for indv in ped_data:
            indv = indv.strip().split()
            # ID, FATHER, MOTHER, SEX (all strings)
            indvs[indv[0]] = Individual(indv[0], indv[1], indv[2], indv[3])

        # represents parents of any founder/married-in
        indvs["0"] = "0"
        return indvs

    def get_married_in(self):
        """
        Return a list of all married-in individuals, identified by having "0"'s
        in place of parents.
        """
        lst = []
        for id_num in self.indvs:
            if id_num == "0":
                continue
            if self.indvs[id_num].married_in:
                  lst.append(id_num)
        return lst

    def find_relations(self):
        """
        Given a list of the most recent generation of a pedigree (all
        individuals which have no children), finds relations between all
        members of that pedigree.
        """
        # find parents, children, etc
        for indv in self.indvs.values():
            self.recursive_find_relations(indv)

    def recursive_find_relations(self, indv):
        if indv == "0":
            return

        # assign parents for each individual
        if indv.p == None or indv.m == None:
            indv.p = self.indvs[indv.p_id]
            indv.m = self.indvs[indv.m_id]

        # base case, don't need to recurse up (no more ancestors in tree)
        if (indv.p == "0" or indv.m == "0"):
            if not (indv.married_in or indv.founder) and len(indv.couples) > 0:
                # if we have spouse information, determine whether indv
                # is married in or founder
                for couple in indv.couples:
                    if couple.p.p == "0" and couple.m.p == "0":
                        self.founders.append(couple)
                        couple.p.founder = True
                        couple.m.founder = True
                if not indv.founder:
                    indv.married_in = True
                    self.married_in.append(indv)
            return

        # base case (only nones that have already been checked have a spouse)
        if indv.checked:
            return

        # recursive case
        indv.checked = True
        indv.p.children.append(indv)
        indv.m.children.append(indv)

        # add Couple information to both parents
        new_couple = Couple(indv.p, indv.m)
        if not new_couple in indv.p.couples:
            indv.p.couples.append(new_couple)
        if not new_couple in indv.m.couples:
            indv.m.couples.append(new_couple)
        indv.parents = new_couple

        # recurse
        self.recursive_find_relations(indv.p)
        self.recursive_find_relations(indv.m)

    def trim_redundant_ancestors(self,ancestor_tree, cohort, verbose=False):
        """
        helper method for find_collective_ca
        takes dictionary of ancestor tree nodes
        returns list of ancestors that are non-redundant sources
        (redundant source: an individual that is an ancestor of another sources
        without any unique paths)
        """
        sources = {}
        for ancestor, node in ancestor_tree.items():
            is_source = True
            for indv in cohort:
                if indv in self.indvs.keys():
                    if not self.indvs[indv] in node.cohort:
                        is_source = False

            if is_source:
                children_paths = [len(c.cohort) for c in node.children]
                #if verbose:
                #    print('is source!')
                #    print(children_paths, len(node.cohort))
                if len(children_paths) > 0: # check for pedigree subsets
                    # is non-redundant
                    if len(node.cohort) > max(children_paths):
                        sources[ancestor] = node
                        # if verbose:
                        #    print('is non-redundant!')

        #if verbose:
        #    print(sources)
        return sources

    def combine_couples(self, ancestor_tree):
        """
        helper method for find_collective_ca
        takes dictionary of ancestor tree nodes
        returns list of ancestors with married couples stored as single node
        """
        sources = {}
        for ancestor, node in ancestor_tree.items():
            #if verbose:
            #    print([str(couple) for couple in node.indv.couples])
            # no remarriage
            if len(node.indv.couples) == 1:
                couple = node.indv.couples[0]
                if couple.id in sources:
                    continue
                if couple.p.id in ancestor_tree and couple.m.id in \
                    ancestor_tree:
                    sources[couple.id] = AncestorNode(couple, node.cohort, \
                        node.children)

                # this could happen if ancestor is genotyped and has IBD
                else:
                    sources[ancestor] = node
            else:
                sources[ancestor] = node

        return sources

    def propogate_cohort(self, tree, indv, cohort):
        """
        Recursively propogate any cohort additions upwards
        that is, in case grandparents have already been created
        and child has been added to ancestor
        (only the case in cross-generational marriages)
        """
        if type(indv.p) == Individual and indv.p.id in tree.keys() and tree[indv.p.id].has_child(indv.id)::
            tree[indv.p.id].cohort.extend(cohort)
            self.propogate_cohort(tree, indv.p, cohort)
        if type(indv.m) == Individual and indv.m.id in tree.keys() and tree[indv.m.id].has_child(indv.id):
            tree[indv.m.id].cohort.extend(cohort)
            self.propogate_cohort(tree, indv.m, cohort)

    def find_collective_ca(self, indvs, verbose=False):
        """
        Takes list of individual ids (strings) and returns the common ancestor
        of the individuals if one exists, otherwise returns all common
        ancestors if multiple exist
        ASSUMPTION: no married-ins are closely related
        """
        # TODO use haplotype values from indvs to be more specific?
        ancestor_tree = {} # ancestor id -> AncestorNode
        queue = [] # next individuals to track

        # first queue all cohort indvs
        for indv in indvs:
            if indv in self.indvs: # check for pedigree subsets
                indv = self.indvs[indv]
                queue.append(indv)
                ancestor_tree[indv.id] = AncestorNode(indv,indv)

        #if verbose:
        #    print(ancestor_tree)

        # build ancestor tree
        while len(queue) != 0:
            indv = queue.pop(0)
            assert(type(indv) == Individual)
            curr_node = ancestor_tree[indv.id]

            # check if married in
            if indv.parents == None or str(indv.parents) == "0":
                continue

            # add to parent nodes
            if indv.p.id not in ancestor_tree.keys():
                # parent node has not been created, initialize parent node
                ancestor_tree[indv.p.id] = AncestorNode(indv.p, \
                    curr_node.cohort, curr_node)
                # add parent to queue
                if indv.p not in queue: # TODO too slow?
                    queue.append(indv.p)
            else:
                # parent node has been created, increment parent's children +
                # cohort
                ancestor_tree[indv.p.id].add_child(curr_node)
                self.propogate_cohort(ancestor_tree, indv.p, curr_node.cohort)

            if indv.m.id not in ancestor_tree.keys():
                # parent node has not been created, initialize parent node
                ancestor_tree[indv.m.id] = AncestorNode(indv.m, \
                    curr_node.cohort, curr_node)
                # add parent to queue
                if indv.m not in queue: # TODO too slow?
                    queue.append(indv.m)
            else:
                # parent node has been created, increment parent's children +
                # cohort
                ancestor_tree[indv.m.id].add_child(curr_node)
                self.propogate_cohort(ancestor_tree, indv.m, curr_node.cohort)

            # TODO break off early if all remaining queue members are sources?

        #if verbose:
        #    print(ancestor_tree)

        # get only most recent ancestor on each path
        sources = self.trim_redundant_ancestors(ancestor_tree, indvs, verbose)
        # if a married pair is a source, treat them as a unit
        sources = self.combine_couples(sources)
        #if verbose:
        #    print(sources)

        return sources

    def descendence_paths(self, sources, cohort):
        """
        Returns list of sources separated into all possible descendence paths.
        """
        if len(sources) == 0:
            print("no sources, cannot find paths")
            return

        # create & fill a dictionary of paths for each individual in the cohort
        cohort_paths = {}
        for indv in cohort:
            cohort_paths[indv] = {}
            for source in sources.values():
                self.get_all_paths(source, indv, source.indv.sex, \
                    cohort_paths[indv])

        # synthesize paths
        source_paths = {}
        for source in sources.values():
            # initialize the paths with the first indv in the cohort
            syn_paths = False
            for indv in cohort:
                if cohort_paths[indv][source] != False:
                    syn_paths = cohort_paths[indv][source]
                if syn_paths != False:
                    break

            if syn_paths == False:
                print("HUGE PROBLEM! no paths for source?")

            # initialize the paths with the first indv in the cohort then add
            # paths to remaining cohort individuals
            for indv in cohort:
                if cohort_paths[indv][source] == False:
                    continue
                path = cohort_paths[indv][source][0]

                # TODO make union/intersect faster using sets for paths
                temp_syn_paths = [path | s_path for s_path in syn_paths]
                for path in cohort_paths[indv][source][1:]:
                    for s_path in syn_paths:
                        if len(path & s_path) < len(path):
                            temp_syn_paths.append(path | s_path)

                syn_paths = temp_syn_paths
            source_paths[source] = syn_paths
            # TODO stop when we've hit max paths

        # each source adds some new paths, but we may still have paths shared
        # by more than one source. TODO get rid of redundant paths
        """
        for source, paths in source_paths:
            for source2, paths2 in source_paths:
                for path in paths:
                    for path2 in paths2:
                        if len(union(path, path2))+1:
                            # TODO get rid of redundant paths
        """

        return source_paths

    def get_all_paths(self, node, target, hap, paths):
        """
        Recursively finds all possible paths from node to target
        node: starting node for path
        target: should be a cohort individual
        hap: haplotype in node(gender of node's ancestor or 0 if node is source)
        paths: a dictionary of paths specific to this one individual
                paths: node -> paths from node to target
                (each path has the form: [(node, hap), (node, hap), ...])
        """
        #-----------------------------------------------------------------
        # Base Cases:

        # path from node has already been undetermined
        if node in paths.keys():
            if paths[node]:
                final = []
                for path in paths[node]:
                    new_path = path.copy()
                    new_path.add((node, hap))
                    final.append(new_path)
                return final
            else:
                return False

        # node is target
        if node.indv.id == target: # TODO string?
            paths[node] = [set()]
            final = [set()]
            final[0].add((node,hap))
            return final

        # node is leaf but not target (node has no children)
        if len(node.children) == 0 and node.indv.id != target:
            paths[node] = False
            return False

        # node does not have target in its cohort
        if target in self.indvs: # check for pedigree subsets
            if not self.indvs[target] in node.cohort:
                paths[node] = False
                return False
        #------------------------------------------------------------------
        # Recursive Case:

        output_paths = []
        for child in node.children:
            child_paths = self.get_all_paths(child, target, node.indv.sex, \
                paths)
            if child_paths:
                for path in child_paths:
                    if path:
                        output_paths.append(path)

        # check if dead end
        if len(output_paths) == 0:
            paths[node] = False
            return False

        paths[node] = output_paths
        final = []
        for path in output_paths:
            new_path = path.copy()
            new_path.add((node, hap))
            final.append(new_path)
        return final

    # TODO make method to create entire descendence array once instead of
    # multiple times

# TODO move the functions below somewhere else (different file or inside class)

def find_num_desc_paths(cohort):
    """Use idea of multi-sets to find number of paths per source."""
    counts = {}
    for c in cohort:
        if c in counts:
            counts[c] += 1
        else:
            counts[c] = 1
    num_paths = 1
    for v in counts.values():
        num_paths *= v
    return num_paths

'''def print_p(paths):
    for path in paths:
        strg = "      "
        for item in path:
            strg += item[0].indv.id + "(" + str(item[1])+  "), "
        print(strg)

def print_path_dict(big_dict):
    print("----")
    for source, paths in big_dict.items():
        print(source)
        if paths:
            for path in paths:
                strg = "      "
                for item in path:
                    strg += item[0].indv.id + "(" + str(item[1])+  "), "
                print(strg)
        else:
            print("     ", paths)'''
