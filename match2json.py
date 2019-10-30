"""
Convert .match file (from GERMLINE) to .json for ease of use later.
Authors: Kelly Finke, Michael Kourakos, Sara Mathieson
Date: 10/28/19
"""
# DONE

# python imports
import json
import optparse
import sys

# our imports
import IBD
from PedigreeTree import PedigreeTree

def parse_args(description):

    parser = optparse.OptionParser(description=description)
    parser.add_option("-g", "--germ_filename", type="string", \
        help="input GERMLINE .match file")
    parser.add_option("-s", "--struct_filename", type="string", \
        help="input .txt file with pedigree structure")
    parser.add_option("-m", "--map_filename", type="string", \
        help="input .map file")
    parser.add_option("-p", "--ped_filename", type="string", \
        help="input .ped file")
    parser.add_option("-j", "--json_filename", type="string", \
        help="output json file (dictionary) of IBDs")

    mandatories = ["germ_filename", "struct_filename", "map_filename", \
        "ped_filename", "json_filename"]
    (opts, args) = parser.parse_args()
    for m in mandatories:
        if not opts.__dict__[m]:
            parser.print_help()
            sys.exit()
    return opts

def store_sequences(IBD_lst, map_filename, ped_filename, pedigree_tree, \
    germ_filename, json_filename):
    """Store IBDs in a json format for ease of use later on"""
    IBDs = IBD_lst[:] # deep copy

    m_data = []
    with open(map_filename, 'r') as m_file:
        for line in m_file:
            tokens = line.strip().split()
            m_data.append(int(tokens[3]))

    big_dict = {}
    j = 0
    with open(json_filename, "w+") as j_file:
        with open(ped_filename, 'r') as p_file:
            for line in p_file:
                if len(IBDs) == 0:
                    break
                indv_id = line.strip().split()[1]
                #print(j, indv_id)
                j += 1
                try:
                    indv = pedigree_tree.indvs[indv_id]
                except:
                    continue
                if indv == "0":
                    continue

                indv_ibds = list(indv.get_IBDs())
                genotype = line.strip().split()[6:]

                for ibd in indv_ibds:
                    if ibd in IBDs:
                        big_dict[ibd.id] = {}
                        IBDs.remove(ibd)
                        haplotype = int(ibd.get_indvs()[indv_id]) - 1
                        if haplotype == -1:
                            haplotype = 0 #if the indv is homozygous for the IBD
                        start_index = binary_search(m_data, int(ibd.start))
                        end_index = binary_search(m_data, int(ibd.end))

                        for i in range(start_index, end_index + 1):
                            # tried changing to int, but json writes str
                            base = m_data[i]
                            allele = genotype[i * 2 + haplotype]
                            big_dict[ibd.id][base] = allele

        json.dump(big_dict, j_file, indent=4)

def binary_search(lst, value):
    """Helper for storing IBDs"""
    low = 0
    high = len(lst) - 1
    while low <= high:
        mid = (low + high) // 2
        if lst[mid] == value:
            return mid
        elif lst[mid] > value:
            high = mid - 1
        else:
            low = mid + 1
    return low

def main():

    # retrieve arguments
    args = parse_args("convert match file to json")

    # construct pedigree data structure
    ped = PedigreeTree(args.struct_filename)

    # construct IBDs data structures (type: list)
    left_out = []
    IBDs = IBD.get_IBDs(args.germ_filename, left_out) # toggle to leave out

    # assign IBDs to individuals in pedigree
    IBD.ibd_to_indvs(IBDs, ped)

    # create json file
    store_sequences(IBDs, args.map_filename, args.ped_filename, ped, \
        args.germ_filename, args.json_filename)

main()
