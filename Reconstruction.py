"""
Class to hold our reconstruction. Holds IBDs, pedigree, and current groups for
each individual, which are updated throughout the algorithm.
Authors: Kelly Finke, Michael Kourakos, Sara Mathieson
Date: 10/30/19
"""

# python imports
import random
from tqdm import tqdm

# our imports
import IBD_content as IBDc
import PedigreeTree as pt
import probabilities
import recombination_points as reco

# TODO keep track of uncertain nodes

# globals
MAX_CONFLICTS = 5
PATH_THRESH = 100000 # if over this number of paths, ignore source for now

class Reconstruction:

    def __init__(self, IBDs, ibd_dict, ped, reconstruct_ids, ancestral_ids, \
        chrom, SNP_lst):
        self.IBDs = IBDs
        self.ibd_dict = ibd_dict
        self.ped = ped
        self.chrom = chrom # int

        # parse chrom info
        start = SNP_lst[0]
        end = SNP_lst[-1]
        num_snps = len(SNP_lst)
        self.chrom_len = [start, end, num_snps]
        self.SNP_lst = SNP_lst

        # ids are just keys of the dictionaries, but useful for initial steps
        self.reconstruct_ids = reconstruct_ids
        self.ancestral_ids = ancestral_ids
        self.reconstruct_groups = {} # key: individual ID, value: list of groups
        self.ancestral_groups = {}   # key: individual ID, value: list of groups

        self.all_conflicts = {}      # key: individual ID, value: conflicts

    def initialize_genotyped(self):
        """
        For the genotyped individuals, initialize their groups (only once).
        """
        num_bad_recon = 0
        for id in self.reconstruct_ids: # tqdm
            indv = self.ped.indvs[id]

            # for genotyped individual, use all certain IBDs (shouldn't be any
            # uncertain right now)
            ibds_certain = indv.get_IBDs()
            if len(ibds_certain) > 0:
                # testing the grouping algorithm on genotyped individuals
                #groups, conflicts = IBDc.group_IBDs(ibds_certain, id, \
                #    self.ibd_dict)
                # using the haplotypes for genotyped individuals!
                groups, conflicts = IBDc.initial_groups(ibds_certain, id, \
                    self.ibd_dict)

                # check we have two groups
                if len(groups) != 2:
                    print('HUGE PROBLEM: genotyped individual', id, 'with', \
                        len(groups), 'groups')
                    #reco.group_conflicts_diagram(groups, conflicts, id, \
                    #    self.chrom, self.chrom_len, "init")
                    num_bad_recon += 1
            else:
                print('HUGE PROBLEM: genotyped individual', id, 'with no IBDs')

            self.reconstruct_groups[id] = groups
            self.all_conflicts[id] = conflicts

        print("num_bad_geno", num_bad_recon)

    def group_all(self, create_images):
        """
        This is a core function we run during each iteration. At a high-level
        we try to add IBDs to already constructed individuals, and try to
        group IBDs from non-reconstructed individuals.
        """
        all_ids = self.reconstruct_ids + self.ancestral_ids

        for id in all_ids: # tqdm
            indv = self.ped.indvs[id]
            ibds_test = indv.get_IBDs()

            if id == 'b':
                print("IBDs placed in b:", len(ibds_test))
                for ibd in ibds_test:
                    print("\t", ibd)

            if len(ibds_test) > 0:

                # very important distinction! if reconstructed, try to add IBDs
                if id in self.reconstruct_ids:
                    groups = self.reconstruct_groups[id]
                    conflicts = self.all_conflicts[id]
                    IBDc.add_IBDs(id, groups, ibds_test, self.ibd_dict, \
                        conflicts)

                # if not reconstructed, group from scratch
                else:
                    groups, conflicts = IBDc.group_IBDs(ibds_test, id, \
                        self.ibd_dict)
                    self.ancestral_groups[id] = groups
                    self.all_conflicts[id] = conflicts

                # save images (this is time-consuming)
                if create_images:
                    print("\n-------------------------------------")
                    print("testing:", id)
                    print("-------------------------------------")
                    if indv.genotyped:
                        flag = "geno"
                    elif id in self.reconstruct_ids:
                        flag = "good"
                    else:
                        flag = "bad"

                    # flip a coin so we don't have so many images
                    #if random.random() < 0.1:
                    #    reco.group_conflicts_diagram(groups, conflicts, id, \
                    #        self.chrom, self.chrom_len, flag)

            # we want this so ancestral IDs matches ancestral groups
            else:
                self.ancestral_groups[id] = [] # no groups yet

    def move_to_constructed(self, good_indvs):
        """
        After each iteration, update reconstruction by adding newly
        reconstructed indvs to "good" and removing them from ancestral
        """
        self.reconstruct_ids.extend(good_indvs)
        for good in good_indvs:
            self.ancestral_ids.remove(good)
            self.reconstruct_groups[good] = self.ancestral_groups[good]
            del self.ancestral_groups[good]

    def no_conflict(self, id, ibd):
        """
        If IBD assigned to individual conflicts with existing groups, we will
        eventually reject the corresponding source.
        """
        if id in self.reconstruct_groups:
            ibd_alleles = self.ibd_dict[str(ibd)]
            groups = self.reconstruct_groups[id]
            num_conflicts = 0

            # look through all groups
            for g in groups:
                overlap,conflicts = IBDc.pairwise_conflicts(g, ibd, ibd_alleles)
                if len(conflicts) > MAX_CONFLICTS:
                    num_conflicts += 1
                else:
                    return True # matches one group

            # if it conflicts with every group, reject
            if num_conflicts == len(groups) and len(groups) >= 2:
                return False
        return True

    def test_source(self, ibd, path_list):
        """
        For all certain nodes from this source (i.e. on all paths), see if
        placing an IBD will cause a conflict. If yes, we can DEFINITELY reject
        this source.
        """

        # first find nodes on all paths
        certain_nodes = path_list[0]
        for path in path_list:
            certain_nodes = certain_nodes & path

        # look through each certain node (note that this doesn't include the
        # source itself (often a couple) and does not include individuals
        # that could get the IBD from either parent (usually a sequenced indv)
        for a in certain_nodes:
            # make sure placing this ibd does not conflict with certain results
            if not self.no_conflict(a[0].indv.id, ibd):
                return False, []
        return True, certain_nodes

    def identify_sources(self):
        """
        Finds the possible sources of the given IBDs using common ancestors of
        individuals possessing IBDs. Then we add the IBD to individuals on all
        paths from all sources (if there are no "bad" sources with many paths).
        """
        perfect_count = 0
        no_sources_count = 0

        source_counts = []
        num_path_counts = []

        # we could sort on a variety of features, but right now this step is
        # independent for each IBD, so no sorting needed
        print("Identifying sources...")
        for ibd in self.IBDs: # tqdm

            # find all potential sources
            # TODO if we already know haplotype, only search that parent
            verbose = False
            ancestors = self.ped.find_collective_ca(list( \
                ibd.get_indvs().keys()))
            source_counts.append(len(ancestors))

            if len(ancestors.keys()) == 0:
                print("HUGE PROBLEM: zero sources!", ibd, ibd.get_indvs())
                no_sources_count += 1
                continue

            if len(ancestors.keys()) == 1:
                perfect_count += 1

            # find an estimate of the number of paths using cohort counts
            # Note: this might be higher or lower than the actual number of
            # paths, depending on path redundancy or multiple haps/path
            num_bad = 0
            for key in ancestors:
                num_paths = pt.find_num_desc_paths(ancestors[key].cohort)
                ancestors[key].num_paths = num_paths # set path estimate
                if num_paths > PATH_THRESH:
                    num_bad += 1
                if num_paths <= 1000:
                    num_path_counts.append(num_paths)

            # assign ancestors to this IBD (we will prune/sort later)
            ibd.set_sources(ancestors)

            # if no ancestors have reasonable numbers of paths, record this
            if len(ancestors.keys()) == num_bad:
                no_sources_count += 1

            # this finds "super certain" placements: individuals on all paths
            # from all sources
            if num_bad == 0:
                paths = self.ped.descendence_paths(ancestors, \
                    list(ibd.get_indvs().keys()))

                # find nodes on all paths from all sources
                k = 0
                for source, path_list in paths.items():
                    # this ensures certain_nodes are in all paths
                    if k == 0:
                        certain_nodes = path_list[0]
                    for path in path_list:
                        certain_nodes = certain_nodes & path
                    k += 1

                # after finding certain nodes, check for conflicts, then place
                for a in certain_nodes:
                    # make sure placing this ibd does not conflict with results
                    if not self.no_conflict(a[0].indv.id, ibd):
                        print('HUGE ISSUE: super certain IBD conflicts')
                        input('enter')

                # assign individuals this IBD
                for a in certain_nodes:
                    cindv = a[0].indv
                    # if -1, don't know hap, so switch to 0 at this point only!
                    if a[1] == -1:
                        chap = "0"
                    else:
                        chap  = str(a[1])
                    #print("indv, hap", cindv.id, a[1])
                    assert cindv.id in self.ped.indvs.keys()
                    cindv.add_ibd(chap+"1", ibd) # add IBD to indv ("1" certain)
                    ibd.set_hap(cindv.id, a[1])  # add indv to IBD as well

        print('perfect', perfect_count)
        print('no_sources_count', no_sources_count)
        print('total', len(self.IBDs))

    def root_all_IBDs(self, IBD_list): # don't use self.IBDs (changes over time)
        """
        Based on the sources for each IBD, choose a source (usually based on
        small number of paths, but could be based on probabilities). Then give
        the IBD to all individuals on all paths from this chosen source.
        """
        not_rooted = 0 # IBDs where we couldn't find a good source

        # keep track of times when the most probable source is also the source
        # with fewest paths
        path_prob_same = 0
        path_prob_denom = 0

        print("Rooting IBDs...")
        for ibd in IBD_list: # tqdm
            sources = ibd.get_sources()

            # given each source, what is the probability of observing cohort?
            #print("\nIBD", ibd.get_indvs())
            for source in sources:
                prob = probabilities.psource(source, ibd, self)
                sources[source].prob = prob # set prob estimate
                #print(source, prob, sources[source].num_paths, sep='\t')
            #print('done one set of sources')

            rooted = False  # found a source that is not obviously bad
            run_out = False # ran out of sources with small numbers of paths
            while not rooted and len(sources) > 0 and not run_out:

                # take the source with the fewest paths
                min_path = min(sources.keys(), key=lambda s: \
                    sources[s].num_paths)
                '''max_path = max(sources.keys(), key=lambda s: \
                    sources[s].num_paths)'''
                # take the source with the highest probability
                max_prob = max(sources.keys(), key=lambda s: sources[s].prob)

                # check if they are the same (for diagnostics)
                if min_path == max_prob:
                    path_prob_same += 1
                path_prob_denom += 1

                # here is where we can make choices! right now: fewest paths
                best_source = min_path
                #best_source = max_path

                # make sure we don't have too many paths
                if sources[best_source].num_paths > PATH_THRESH:
                    run_out = True
                    #del sources[best_source]

                else:
                    # compute paths for just this source
                    mini_dict = {best_source: sources[best_source]}
                    paths = self.ped.descendence_paths(mini_dict, \
                        list(ibd.get_indvs().keys()))
                    path_list = paths[sources[best_source]]

                    # TODO here is where we could randomly choose a path

                    # make sure source doesn't conflict with existing results
                    result, certain_nodes = self.test_source(ibd, path_list)

                    # if no conflict, add to certain nodes
                    if result:

                        # here (and elsewhere) I want to add indv to IBDs so
                        # that we can remove it later if necessary
                        #certain_ids = []
                        for a in certain_nodes:
                            cindv = a[0].indv
                            #certain_ids.append(cindv.id)
                            # if -1, don't know hap, so switch to 0
                            if a[1] == -1:
                                chap = "0"
                            else:
                                chap  = str(a[1])
                            #print("indv, hap", cindv.id, a[1])
                            assert cindv.id in self.ped.indvs.keys()
                            # add IBD to indv ("1" certain)
                            cindv.add_ibd(chap+"1", ibd)
                            # add indv to IBD as well
                            ibd.set_hap(cindv.id, a[1])

                        # when only one source (i.e. not couple), add IBD to it
                        if "&" not in best_source:
                            sindv = self.ped.indvs[best_source]
                            sindv.add_ibd("01", ibd) # unknown hap, certain
                            ibd.set_hap(best_source, -1) # unknown hap

                        rooted = True
                        ibd.change_source(best_source)

                    # if there was a conflict, reject this source and move on
                    else:
                        del sources[best_source]

            if not rooted:
                not_rooted += 1

        print("not rooted", not_rooted)
        if path_prob_denom > 0:
            print("frac path&prob same:", path_prob_same, path_prob_denom, \
                path_prob_same/path_prob_denom)
        else:
            print("frac path&prob same:", path_prob_same, path_prob_denom)

    def find_bad_sources(self):
        """
        Based on previously grouped IBDs, identify which groupings are strong
        (i.e. 2 long groups with lots of support), and which groups are weak
        (i.e. multiple groups or two strong with some other poorly placed IBDs)
        """

        # keep track of good indvs (basically reconstructed)
        good_indvs = []
        mixed_indvs = []
        bad_ibds = set() # these IBDs definitely have bad sources

        for id in self.ancestral_ids:
            indv = self.ped.indvs[id] # Individual object
            groups = self.ancestral_groups[id]
            lengths = [len(g) for g in groups]
            num_ibds = [len(g.ibds) for g in groups]

            # 3 main cases: "perfect: two strong groups", "mixed: two strong
            # groups + some small incorrectness", "mess: many groups"
            result, i, j = self.analyze_groups(groups)

            if result == "perfect":
                # flag individual as "good"
                good_indvs.append(id)
                # give support to sources that placed IBDs here
                for g in groups:
                    for ibd in g.ibds:
                        ibd.add_vote_plus() # vote for current source

            # here we also need to remove ibds that were in the "bad" groups!
            elif result == "mixed":
                # flag individual as "good"
                good_indvs.append(id)
                mixed_indvs.append(id)
                # give support to 2 strong groups and reject others
                for g in [groups[i], groups[j]]: # should have been sorted
                    #assert good_group(g)
                    for ibd in g.ibds:
                        ibd.add_vote_plus() # vote for current source

                # later we will change their sources and remove from individuals
                for k in range(len(groups)):
                    if k != i and k != j:
                        g = groups[k]
                        for ibd in g.ibds:
                            bad_ibds.add(ibd) # set
                            # need to also remove from indv, etc
                            if not indv.genotyped:
                                # TODO why do we sometimes have IBDs not in indv
                                result = indv.remove_ibd(ibd)
                                if result:
                                    ibd.remove_indv(id)

                # omit other groups
                self.ancestral_groups[id] = [groups[i], groups[j]]
                # TODO also remove conflicts?

            elif result == "mess":
                # reduce support for all IBDs, but some could still be good
                for g in groups:
                    for ibd in g.ibds:
                        ibd.add_vote_minus() # vote against current source

            else:
                # do nothing if we have only one group
                # TODO think about what to do if we have one strong group
                pass

        print("mixed", len(mixed_indvs), mixed_indvs)
        return good_indvs, bad_ibds

    def analyze_groups(self, groups):
        """See which situation we're in: 2 strong groups, mixed, or a mess"""
        # TODO could sort on number of IBDs or other criteria
        groups = sorted(groups, key=lambda g: len(g), reverse=True)

        if len(groups) == 2:
            if self.good_group2(groups[0]) and self.good_group2(groups[1]):
                return "perfect", 0, 1
            else:
                # TODO maybe later on, if one is good, save that one
                return "mess", None, None

        elif len(groups) > 2:
            result, i, j = self.two_strong(groups)
            if result:
                return "mixed", i, j
            else:
                return "mess", None, None

        # only one group, but could be good (TODO)
        else:
            return "one"

    def good_group2(self, group):
        """Thresholds for a good group when there are exactly 2 groups"""

        #chr_len = self.chrom_len[1] - self.chrom_len[0]
        chr_snps = self.chrom_len[2]
        group_snps = len(group.template)
        # few IBDs, but long and good coverage
        if group_snps >= chr_snps*0.9 and len(group.ibds) >= 1:
            return True
        # more IBDs, bit less strict on length
        if group_snps >= chr_snps*0.7 and len(group.ibds) >= 3:
            return True
        # not strict on length, but many IBDs
        return group_snps >= chr_snps*0.5 and len(group.ibds) >= 10

    # no longer using: using #SNPs as a more robust measure of what fraction
    # of the chrom the IBD takes up
    '''def good_coverage(self, group, thresh):
        """Estimate the coverage based on the number of SNPs"""
        num_snps = len(group.template)
        ibd_len = len(group)
        #print(ibd_len/num_snps)
        # for sequenced, we have about 474 SNPs/Mb
        return (num_snps/ibd_len)*1000000 >= thresh'''

    def two_strong(self, groups):
        """
        Return true (along with group indices) if we found two strong groups,
        relative to other groups
        """
        chr_len = self.chrom_len[1] - self.chrom_len[0]
        groups = sorted(groups, key=lambda g: len(g.ibds), reverse=True)

        # second and third groups difference: either more IBDs or longer
        if len(groups[1].ibds) >= 2*len(groups[2].ibds) or len(groups[1]) >= \
            2*len(groups[2]):

            # first two both good
            if self.good_group2(groups[0]) and self.good_group2(groups[1]):
                    return True, 0, 1
        return False, None, None

    def remove_bad_sources(self, bad_ibds):
        """Remove IBDs with bad sources""" # TODO also remove conflicts?

        # go through bad IBDs, removing them from ancestral individuals and
        # their groups
        for ibd in bad_ibds:
            # individuals that have this IBD
            indv_ids = list(ibd.get_indvs().keys())
            for id in indv_ids:

                # remove ibd from indv and indv from ibd
                if id in self.ancestral_groups: # this check is very important!
                    indv = self.ped.indvs[id]
                    indv.remove_ibd(ibd)
                    ibd.remove_indv(id)
                    #del ibd.indvs[id]

                    # also remove ibd from any groups of the individual
                    groups = self.ancestral_groups[id]
                    for g in groups:
                        if ibd in g.ibds:
                            g.ibds.remove(ibd)

            # remove here, later on we will reassign source for this IBD
            ibd.remove_source(ibd.get_curr_source())

    def bad_from_votes(self, IBDs):
        """
        Check if many individuals downvoted this IBD. Currently we don't use
        this information.
        """
        bad_ibds = []
        for ibd in IBDs:
            # if we have more than twice as many votes for bad as good
            if ibd.get_vote_minus() > ibd.get_vote_plus()*2 and \
                ibd.get_vote_minus() >= 10:
                #print(ibd, ibd.vote_plus, ibd.vote_minus)
                bad_ibds.append(ibd)

        return bad_ibds

    def write_vcf(self, out_filename, genotyped_only=False):
        """
        For all non-genotyped individuals, write out the parts we could
        reconstruct to a VCF-like file
        """
        all_ids = self.reconstruct_ids + self.ancestral_ids

        # loop through all SNPs
        out_file = open(out_filename,'w')
        for snp in self.SNP_lst:
            #print(snp)
            allele_lst = []
            for id in all_ids:
                indv = self.ped.indvs[id]

                # either write out all genotyped or no genotyped
                if (not indv.genotyped and not genotyped_only) or (indv.genotyped and genotyped_only):

                    # reconstructed: add snp if part of groups
                    if id in self.reconstruct_ids:
                        groups = self.reconstruct_groups[id]
                        hap1 = groups[0].get_allele(snp)
                        hap2 = groups[1].get_allele(snp)
                        allele_lst.extend([hap1,hap2])

                    # unreconstructed: add missing data
                    else:
                        allele_lst.extend(["?","?"])

            out_file.write(str(snp) + " " + " ".join(allele_lst) + "\n")

        out_file.close()

    def write_ped(self, out_filename, genotyped_only=False, left_out=None):
        """
        For all non-genotyped individuals, write out the parts we could
        reconstruct to a PED-like file
        """
        all_ids = self.reconstruct_ids + self.ancestral_ids

        # loop through all indvs
        out_file = open(out_filename,'w')
        #print("writing file\n\n")

        for id in all_ids:
            indv = self.ped.indvs[id]

            # either write out all genotyped or no genotyped or write just one
            if (not indv.genotyped and not genotyped_only and left_out == None) or (indv.genotyped \
                and genotyped_only and left_out == None) or (left_out == id):

                indv_lst = self.write_one(id, indv)
                #print("indv_lst", indv_lst)
                out_file.write(" ".join(indv_lst) + "\n")

        out_file.close()

    def write_one(self, id, indv):
        # for ped format
        indv_lst = ["1", id, indv.p_id, indv.m_id, str(indv.sex), "-9"]

        for snp in self.SNP_lst:

            # reconstructed: add snp if part of groups
            if id in self.reconstruct_ids:
                groups = self.reconstruct_groups[id]
                hap1 = groups[0].get_allele(snp)
                if len(groups) > 1:
                    hap2 = groups[1].get_allele(snp)
                indv_lst.extend([hap1,hap2])

            # unreconstructed: add missing data
            else:
                indv_lst.extend(["?","?"])

        #print("indv_lst2", indv_lst)
        return indv_lst
