"""
Utility: computes global results and aids in file parsing.
Authors: Kelly Finke, Michael Kourakos, Sara Mathieson
Date: 10/16/19
"""
# DONE

def parse_chrom(chrom_filename):
    """
    Parse genetic-map-like file to obtain locations of genotyped SNPs. We will
    output these locations later in the final VCF file.
    """
    chrom = None
    SNP_lst = []

    chrom_file = open(chrom_filename,'r')
    for line in chrom_file:
        # chrom, SNP name, reco rate, SNP location
        tokens = line.strip().split()
        chrom = int(tokens[0])
        SNP_lst.append(int(tokens[-1]))
    chrom_file.close()

    return chrom, SNP_lst

def no_seq_descendants(indvs):
    """Find individuals who have no sequenced descendants (can't reconstruct)"""
    no_seq_all = []
    seq_geno = [] # genotyped, with genotyped descendants

    for id,indv in indvs.items():
        if id == "0":
            continue # "parent" of founders and married-in
        descend = indv.descendants()
        no_seq = True
        num_seq = 0
        for d in descend:
            if d.genotyped:
                no_seq = False
                num_seq += 1
        if no_seq and not indv.genotyped:
            no_seq_all.append(id)
        if not no_seq and indv.genotyped:
            seq_geno.append(id + ": " + str(num_seq))

    print("seq geno", "\n".join(seq_geno))
    return no_seq_all
