"""
Computes probability IBD came from a specific source.
Authors: Kelly Finke, Michael Kourakos, Sara Mathieson
Date: 8/13/19
"""
# DONE

# TODO maybe move to Reconstruction?
import math

def ptransmit(Fp, Mp, preco):
    """
    Compute the probability children have IBD, given the parental probabilities
    of having the IBD (takes into account recombination rate).
    preco = prob of recombination within IBD segment (thus breaking it up)
    """
    p0F = Fp[0]
    p1F = Fp[1]
    p2F = Fp[2]
    p0M = Mp[0]
    p1M = Mp[1]
    p2M = Mp[2]

    # no recombination breaking up IBD
    pnoreco = 1-preco

    # 0, 1, or 2 copies of the IBD
    newP0 = (p0F + 0.5*p1F*(1+preco)) * (p0M + 0.5*p1M*(1+preco))
    newP1 = (p0F + 0.5*p1F*(1+preco)) * (0.5*p1M*pnoreco + p2M) + \
        (0.5*p1F*pnoreco + p2F)*(p0M + 0.5*p1M*(1+preco))
    newP2 = (0.5*p1F*pnoreco + p2F) * (0.5*p1M*pnoreco + p2M)
    return [newP0, newP1, newP2]

def preco_cM(cM):
    """
    Given genetic distance in cM, compute probability of recombination.
    See https://en.wikipedia.org/wiki/Centimorgan
    """
    return (1-math.exp(-2*cM/100))/2

def psource(source, ibd, reconstruction):
    """
    Estimate probability that indv (type: Individual) is the source of this IBD
    """
    # first: figure out probability of recombination
    preco = preco_cM(ibd.genetic_distance)

    # keep track of individuals we've changed pibd for (reset to None after)
    to_reset = set()

    # couple
    if '&' in source:
        split = source.index('&')
        id1 = source[:split]
        indv = reconstruction.ped.indvs[id1]
        id2 = source[split+1:]
        if id1 in ibd.get_indvs(): # source has IBD
            return 1 # p=1
        if id2 in ibd.get_indvs(): # source's spouse has IBD
            return 1 # p=1

    # single
    else:
        indv = reconstruction.ped.indvs[source]
        if source in ibd.get_indvs():
            return 1

    total_P = {}
    num_have_ibd = 0

    # don't use directly, lists are mutable!
    queue = [ch.id for ch in indv.children]

    # assume 1 copy of the IBD
    indv.set_pibd([0, 1, 0])
    to_reset.add(indv)
    while len(queue) != 0:
        # if we have already seen this indv (i.e. loop), then update
        # automatically
        nextid = queue.pop(0)
        next = reconstruction.ped.indvs[nextid]

        # checks on parents (F/M arbitrary)
        F = next.parents.p
        M = next.parents.m
        if F.get_pibd() == None and M.get_pibd() != None:
            F.set_pibd([1, 0, 0])
            to_reset.add(F)
        elif M.get_pibd() == None and F.get_pibd() != None:
            M.set_pibd([1, 0, 0])
            to_reset.add(M)
        elif F.get_pibd() == None and M.get_pibd() == None:
            # TODO maybe this is okay if we've considered many children?
            print('HUGE PROBLEM! both None', F.get_pibd(), M.get_pibd())

        # calculate and set probability of having IBD
        next.set_pibd(ptransmit(F.get_pibd(), M.get_pibd(), preco))
        to_reset.add(next)
        if next.genotyped:

            if next.id in ibd.get_indvs():
                if next.id not in total_P:
                    num_have_ibd += 1
                if ibd.get_hap(next.id) == 0: # int
                    # take prob of two copies
                    total_P[next.id] = next.get_pibd()[2]
                else:
                    # take prob of one copy
                    total_P[next.id] = next.get_pibd()[1]

            # TODO is this the right thing? not considering individuals that
            # don't have the IBD when normalizing probabilities
            else:
                pass

        # this enforces that we stop when we've hit sequenced individual
        else:
            queue.extend([ch.id for ch in next.children])

    # reset pibds we've changed
    for individual in to_reset:
        individual.reset_pibd()

    probs = total_P.values()
    if len(probs) == 0:
        return 0

    # normalize by number of genotyped individuals (TODO: is this a good thing?)
    return sum(probs)/len(probs)
