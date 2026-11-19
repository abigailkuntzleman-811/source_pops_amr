"""IBD and admixture proportions from the simulated admixture models.

Shared by the Figure 6 scripts.

IBD from tskit: ts.ibd_segments reports a pair's shared ancestry in pieces
that break wherever the local tree changes, even inside one long tract. A
length cutoff (min_span) applied to those raw pieces throws away real IBD, so
here touching pieces are joined first and the cutoff is applied to whole
tracts, which approximates what hap-ibd detects in the real data. Only shared
ancestors younger than MAX_TIME generations count as IBD; older ones almost
never leave a tract of MIN_CM (mean tract length is 50/t cM).
"""

REC_RATE = 1e-8                    # recombination rate in the simulations
BP_PER_CM = 1 / (REC_RATE * 100)   # 1 Mb = 1 cM
MIN_CM = 2.0                       # hap-ibd minimum output length
MAX_TIME = 500                     # generations

# Models drawn in Figure 6: label -> (folder under SIM_DIR, tract file pattern)
MODELS = {
    "Standard":                  ("no_bottleneck_donor",     "{pop}_7tracts_rep{i}_stitched.csv.gz"),
    "Single 10% European donor": ("eur_donor",               "{pop}_10tracts_rep{i}_stitched.csv.gz"),
    "Two 10% European donors":   ("mxl_double_donor_recent", "{pop}_merged_89tracts_rep{i}.csv.gz"),
    "Six 10% European donors":   ("mxl_donor6",              "{pop}_totalIBSmerged_rep{i}.csv.gz"),
}


def _pop_id(ts, pop):
    return next(p.id for p in ts.populations() if p.metadata.get("name") == pop)


def admix_props(tracts, ts, pop):
    """European admixture proportion per individual: tract length / sequence length."""
    pop_id = _pop_id(ts, pop)
    totals = tracts.groupby("hap_id")["length"].sum()
    return {ind.id: totals.get(ind.id, 0) / ts.sequence_length
            for ind in ts.individuals()
            if ts.node(ind.nodes[0]).population == pop_id}


def ibd_per_individual(ts, pop, min_cm=MIN_CM, max_time=MAX_TIME):
    """Total IBD (cM) each individual in `pop` shares with everyone else in `pop`."""
    pop_id = _pop_id(ts, pop)
    inds = [ind.id for ind in ts.individuals()
            if ts.node(ind.nodes[0]).population == pop_id]
    node_to_ind = {n: i for i in inds for n in ts.individual(i).nodes}

    segments = ts.ibd_segments(within=list(node_to_ind), max_time=max_time,
                               store_segments=True)
    ibd = {i: 0.0 for i in inds}
    for (n1, n2), seg_list in segments.items():
        i1, i2 = node_to_ind[n1], node_to_ind[n2]
        if i1 == i2:                         # an individual's own two copies
            continue
        tracts = []                          # join touching pieces into tracts
        for s in sorted(seg_list, key=lambda s: s.left):
            if tracts and tracts[-1][1] == s.left:
                tracts[-1][1] = s.right
            else:
                tracts.append([s.left, s.right])
        cm = sum((r - l) / BP_PER_CM for l, r in tracts
                 if (r - l) / BP_PER_CM >= min_cm)
        ibd[i1] += cm
        ibd[i2] += cm
    return ibd


def ibd_per_individual_original(ts, pop):
    """The original get_ibd_dict, for comparison: min_span on raw pieces, cM."""
    pop_id = _pop_id(ts, pop)
    inds = [ind.id for ind in ts.individuals()
            if ts.node(ind.nodes[0]).population == pop_id]
    node_to_ind = {n: i for i in inds for n in ts.individual(i).nodes}
    segments = ts.ibd_segments(within=list(node_to_ind), min_span=1_000_000,
                               store_segments=True)
    ibd = {i: 0.0 for i in inds}
    for (n1, n2), seg_list in segments.items():
        bp = sum(s.right - s.left for s in seg_list)
        ibd[node_to_ind[n1]] += bp / BP_PER_CM
        ibd[node_to_ind[n2]] += bp / BP_PER_CM
    return ibd


def ibd_pairs(ts, pop, min_cm=MIN_CM, max_time=MAX_TIME):
    """Total IBD (cM) for every pair of individuals in `pop`, 0 if they share none.

    Same tract definition as ibd_per_individual. With haploid samples, each
    individual is one haplotype, so these are haplotype pairs.
    """
    from itertools import combinations

    pop_id = _pop_id(ts, pop)
    inds = [ind.id for ind in ts.individuals()
            if ts.node(ind.nodes[0]).population == pop_id]
    node_to_ind = {n: i for i in inds for n in ts.individual(i).nodes}

    segments = ts.ibd_segments(within=list(node_to_ind), max_time=max_time,
                               store_segments=True)
    ibd = {pair: 0.0 for pair in combinations(sorted(inds), 2)}
    for (n1, n2), seg_list in segments.items():
        i1, i2 = node_to_ind[n1], node_to_ind[n2]
        if i1 == i2:
            continue
        tracts = []
        for s in sorted(seg_list, key=lambda s: s.left):
            if tracts and tracts[-1][1] == s.left:
                tracts[-1][1] = s.right
            else:
                tracts.append([s.left, s.right])
        ibd[(min(i1, i2), max(i1, i2))] += sum(
            (r - l) / BP_PER_CM for l, r in tracts if (r - l) / BP_PER_CM >= min_cm)
    return list(ibd.values())
