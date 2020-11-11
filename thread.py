"""
Main file to run entire reconstruction pipeline on a given pedigree and sequence
data.
Authors: Kelly Finke, Michael Kourakos, Sara Mathieson
Date: 10/30/19
"""

# python imports
import json
import optparse
import sys

# our imports
import IBD
from PedigreeTree import PedigreeTree
from Reconstruction import Reconstruction
import util

def parse_args(description):

    parser = optparse.OptionParser(description=description)
    parser.add_option("-g", "--germ_filename", type="string", \
        help="input GERMLINE .match file")
    parser.add_option("-s", "--struct_filename", type="string", \
        help="input .txt file with pedigree structure")
    parser.add_option("-m", "--map_filename", type="string", \
        help="input .map file")
    parser.add_option("-j", "--json_filename", type="string", \
        help="input json file (dictionary) of IBDs")
    parser.add_option("-p", "--ped_filename", type="string", \
        help="output .ped file of reconstructed individuals")
    parser.add_option("-x", "--max_prob", action="store_true",dest="max_prob", \
        help="include -x flag to change source assignment to max probability")

    mandatories = ["germ_filename", "struct_filename", "map_filename", \
        "json_filename"]
    (opts, args) = parser.parse_args()
    for m in mandatories:
        if not opts.__dict__[m]:
            parser.print_help()
            sys.exit()
    return opts

def main():
    """Run entire pipeline"""

    # retrieve arguments
    args = parse_args("pedigree args")

    # construct pedigree data structure
    ped = PedigreeTree(args.struct_filename)

    # leave one out approach (sequenced individual we wish to reconstruct)
    left_out = []

    # construct IBDs data structures (type: list)
    IBDs = IBD.get_IBDs(args.germ_filename, left_out) # toggle to leave out
    #IBDs = IBDs[:100] # testing for speed

    # assign IBDs to individuals in pedigree
    IBD.ibd_to_indvs(IBDs, ped)

    # cohort size info
    cohort_dist = [len(ibd.get_indvs()) for ibd in IBDs]
    #print("cohort min size", min(cohort_dist))
    #print("cohort max size", max(cohort_dist))

    ############################################################################
    # Pre-process IDs

    # find chromosome and SNP info
    chrom, SNP_lst = util.parse_chrom(args.map_filename)
    print("CHROM", chrom)

    # get IBD sequence information to use throughout
    with open(args.json_filename) as json_data:
        ibd_dict = json.load(json_data)

    # find indvs with no genotyped descendants (won't be able to reconstruct)
    no_seqd = util.no_seq_descendants(ped.indvs)

    # separate into genotyped and not initially (later on these will be
    # reconstructed and "yet to be constructed"
    # the idea is that hopefully we move most the ids to reconstruct at the end!
    reconstruct_ids = []
    ancestral_ids = []
    for id,indv in ped.indvs.items():
        if id == "0":
            continue
        if indv.genotyped:
            reconstruct_ids.append(id)
        elif id not in no_seqd:
            ancestral_ids.append(id)

    print("----------------")
    num_geno = len(reconstruct_ids)
    print("Genotyped:", num_geno)
    print("No genotyped descendants:", len(no_seqd))
    print("To reconstruct:", len(ancestral_ids))
    print("----------------")

    # finally, create a Reconstruction object that we will modify throughout
    reconstruction = Reconstruction(IBDs, ibd_dict, ped, reconstruct_ids, \
        ancestral_ids, chrom, SNP_lst, source_max_prob=args.max_prob)

    ############################################################################
    # FIRST: things we do only once

    # reconstruction first step: group IBDs from genotyped individuals
    print("Grouping IBDs from genotyped individuals...")
    reconstruction.initialize_genotyped()

    # identify sources and add them to IBDs
    print("Identify sources for IBDs...")
    reconstruction.identify_sources()

    ############################################################################
    # NOW ENTER ITERATIVE LOOP

    # now reconstruct!
    print("Reconstructing...")
    building = True
    num_recon = len(reconstruction.reconstruct_ids)

    # continue looping until we are not building new individuals
    iter = 0
    while building:
        print("\nITERATION", iter)

        # start from source with smallest number of paths and work from there
        # the first time, we use IBDs, then we need to use bad_ibds
        if iter == 0:
            reconstruction.root_all_IBDs(IBDs)
        else:
            reconstruction.root_all_IBDs(bad_ibds)

        # group unsequenced individuals
        print("Grouping IBDs from remaining individuals...")
        reconstruction.group_all(False)

        print("Analyzing groups...")
        good_indvs, bad_ibds = reconstruction.find_bad_sources()
        print('num good indvs', len(good_indvs))
        print(good_indvs)
        print('num bad ibds', len(bad_ibds))

        # TODO if we have conflicts with "super certain" IBDs that is bad -
        # shouldn't remove them

        # after noting these conflicts, update!
        # 1) move good_indvs out of ancestral and into reconstructed, groups too
        reconstruction.move_to_constructed(good_indvs)

        # 2) remove bad ibds from all individuals, and remove associated source
        # TODO double check we don't remove from good reconstructions
        reconstruction.remove_bad_sources(bad_ibds)

        # 3) go through votes and determine bad IBDs
        # TODO we could remove them from the groups as well, but not certain
        other_bad_ibds = reconstruction.bad_from_votes(IBDs)
        print("num bad ibds from votes", len(other_bad_ibds))

        # now repeat rooting with bad_ibds (no leave one out)
        print("Non-geno Reconstructed", len(reconstruction.reconstruct_ids) - \
            num_geno)

        # check if we are finished
        if len(reconstruction.reconstruct_ids) == num_recon:
            building = False
        else:
            num_recon = len(reconstruction.reconstruct_ids)
        iter += 1

    # after, optionally save images TODO shouldn't regroup here for runtime
    reconstruction.group_all(False)

    # write out .ped file (optional)
    if args.ped_filename != None:
        reconstruction.write_ped(args.ped_filename)

if __name__ == "__main__":
    main()
