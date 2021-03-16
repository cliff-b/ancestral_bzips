# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from Bio import Phylo
from Bio.Phylo import BaseTree
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import itertools


def read_tree(treefile, lower=True, silent=True):
    """Takes a phylogenetic tree file
    reads it in and returns the entire tree as well as
    nodes and tips separately
    Takes:
        treefile: filename of a rooted tree in .newick format
        silent: Boolean, set to True to get the number of nodes and tips in the tree
    Returns:
        tree: a tree object
        nodes: a list of clades that are the nodes
        tips: a list of clades that are the tips"""
    tree = Phylo.read(treefile, 'newick', rooted=True)
    nodes = Phylo.BaseTree.TreeMixin.get_nonterminals(tree)
    tips = Phylo.BaseTree.TreeMixin.get_terminals(tree)
    if not silent:
        print('Reading in tree..')
        print(treefile + ' contains ' + str(len(nodes)) + ' nodes')
        print(treefile + ' contains ' + str(len(tips)) + ' taxa')
    if lower is True:
        for t in tree.clade:
            if type(t.name) is str:
                t.name = t.name.lower()
        for t in tips:
            t.name = t.name.lower()
    return tree, tips, nodes


def load_experiment(filename):
    """Loads a file in .csv format with interaction scores in column
    twelve. Makes a dictionary of sets of the peptides to link to the 
    interaction score. Theoretically used to fix aliasing present both in names and hidden from names'
    Takes:
        filename: string, of a file with peptides in cols 1, 2 and interactions in 12
    returns:
         data_dic:   a dictionary of frozen sets of strings, key =frozenset({'pep1', 'pep2'}) value = interaction score"""
    f = open(filename, 'r')
    lines = f.readlines()
    data_dic = {}
    for i in lines[1:]:
        line = i.split()[0].split(',')
        if len(line[0].split('+')) > 1 or len(line[1].split('+')) > 1:
            for j in line[0].split('+'):
                for k in line[1].split('+'):
                    try:
                        data_dic[frozenset([j.lower(), k.lower()])] = np.average([data_dic[frozenset([j.lower(), k.lower()])], float(line[12])])
                    except KeyError:
                        data_dic[frozenset([j.lower(), k.lower()])] = float(line[12])
        else:
            try:
                data_dic[frozenset([line[0].lower(), line[1].lower()])] = np.average([data_dic[frozenset([line[0].lower(), line[1].lower()])], float(line[12])])
            except KeyError:
                data_dic[frozenset([line[0].lower(), line[1].lower()])] = float(line[12])
    return data_dic


def read_alignment(ali_file, seq_range):
    """Opens alignment in phylip format (easier than fasta)
        trims the sequence to just the AAs in seq_range
        self note: #90:131 is bzip domain  #65:131 is what I use for ASR, with 8cat
        Takes: 
            ali_file: file in phylip format
            seq_range: a numeric two item list 
        Returns: 
            dictionary of key = species, values = sequences"""
    seq_dict = {}
    f = open(ali_file, 'r')
    lines = f.readlines()
    f.close()
    for i in lines[1:]:  # Ignore first line in Phylip format
        seq_dict[i.split()[0]] = i.split()[1][seq_range[0]:seq_range[1]]
        if 'X' in i.split()[1][seq_range[0]:seq_range[1]]:
            print(i.split()[0] + ' has an X amino acid somewhere')
    return seq_dict    


def find_duplication_nodes(tree, duplication_file):
    """function to read in a list of extant proteins with a duplicated ancestor
    list is of form 'pep1,pep2', function locates the last common ancestor on the tree
    and returns a list of those clades
    Takes:
        tree: a tree object
        duplication_file: filename of a file containing proteins that were duplicated on the tree object
    Returns:
       duplication_nodes: a list of clade objects
    """
    # currently the duplication nodes returned are 173 269 268 310 317 323 248 239 248 174 175 200 201 195
    # 190 296 294 293 320 219 211 203
    duplication_nodes = []
    f = open(duplication_file, 'r')
    lines = f.readlines()
    f.close()
    pairs = lines[0].split()
    print(pairs)
    for i in pairs:
        i = i.lower()
        duplication_nodes.append(Phylo.BaseTree.TreeMixin.common_ancestor(tree, [i.split(',')[0], i.split(',')[1]]))
    return duplication_nodes


def find_orthologs(tree, data_dic, duplication_list):
    """Takes a tree and interaction data and gives back a list of orthologs
        importantly these orthologs are only ones where the homodimers of the two
        interact and the ancestral proteins homodimers interact
    Takes:
        data_dic: a dictionary frozensets to floats of interaction data in the form key = frozenset({'173', '174'}), value = 1.0
        tip_dic: dictionary of strings to strings to convert names from lowercase in the data_dic to mixedcase in the tree
        duplication_list: a list of clade objects that were points of gene duplication
    Returns:
        ortho_dic: a dictionary of key = frozenset and value = list of two numbers, distance between orthologs and whether they interact or not """
    print('Finding orthologs')
    ortho_dic = {}
    yes, no, missingpoints = 0, 0, 0
    for i in data_dic.keys():
        if len(list(i)) > 1:
            protein1, protein2 = list(i)
            ancestor = Phylo.BaseTree.TreeMixin.common_ancestor(tree, protein1, protein2)
            try:
                # check that the homodimers for protein1, protein2 and ancestor all homodimerize
                if ancestor not in duplication_list and data_dic[frozenset([ancestor])] == 1 and data_dic[frozenset([protein1])] == 1 and data_dic[frozenset([protein2])] == 1:
                    if data_dic[i] == 1:
                        ortho_dic[i] = [Phylo.BaseTree.TreeMixin.distance(tree,  protein1, protein2), 1]
                        yes += 1
                    elif data_dic[i] == 0:
                        ortho_dic[i] = [Phylo.BaseTree.TreeMixin.distance(tree, protein1, protein2), 0]
                        no += 1
            except KeyError:
                missingpoints += 1
    print('Done finding orthologs', '\nNumber of noninteractors: ', no, '\nNumber of interactors:', yes)
    print("Missing points: ", missingpoints)
    return ortho_dic


def pick_nonoverlap(tree, pairs):
    """finds a set of proteins from the same tree that do not share any branches on the tree and records their interaction scores and branch length
    Takes:
        tree: a tree object
        pairs: a dictionary of frozensets of protein pairs on the tree to lists of two numerics, in the form of key = frozenset({'173', '174'}),
                values = [total distance between the two proteins in the key, whether they interact]
    Returns:
        values:
            a list of numerics of whether the proteins interacted or not
        distances:
            a list of numerics of how far apart the proteins were"""
    for i in pairs.keys():
        pairs[i].append(Phylo.BaseTree.TreeMixin.trace(tree, list(i)[0], list(i)[1]) + [list(i)[0]])
    values = []
    distances = []
    keys = list(pairs.keys())
    while len(keys) > 0:
        key = np.random.choice(keys)
        values.append(pairs[key][1])
        distances.append(pairs[key][0])
        # Here we test for overlap in path. if i has any shared nodes with key then we delete it from the dicitonary
        for i in keys:
            try:
                if len(intersection(pairs[i][2], pairs[key][2])) > 0 and i != key:
                    del pairs[i]
            except KeyError:
                pass
        del pairs[key]
        keys = list(pairs.keys())
    return values, distances    


def convert_data_dic_to_species_names(data_dic, nodes, tips):
    """The data piped in from R is not in clade objects or anywhere close to it. This function finds where the things was on the tree and gets them in tree form
    Takes:
        data_dic: a dictonary of interaction scores, data_dic[frozenset([str(protein1), str(protein2)])] to values = floats, 0.0, 0.5, 1.0 depending on the interaction
        nodes: a list of nodes from the tree we're using
        tips: a list of tips from the tree we're using
    Returns:
        new_dic: a dictionary similar to data_dic but where all the names are lower case and represented by clade objects so new_dic[frozenset([clade1, clade2])] to floats
        """
    new_dic = {}
    for t in tips:
        t.name = t.name.lower()
    for i in data_dic.keys():
        val = data_dic[i]
        newkey = []
        for j in list(i):
            if j.isdigit():
                k = nodes[int(j)-172]
                newkey.append(k)
            else:
                try:
                    newkey.append([x for x in tips if x.name == j][0])
                except KeyError:
                    pass
        new_dic[frozenset(newkey)] = val
    return new_dic


def neutral_evo(tree, data_dic, duplication_nodes, matched_dic):
    bins = [0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2, 2.5, 3]
    filtered_orthologs = find_orthologs(tree, data_dic, duplication_nodes)
    # distances_tot = [dist[0] for dist in filtered_orthologs.values()]
    # values_tot = [score[1] for score in filtered_orthologs.values()]
    # pairs = list(filtered_orthologs.keys())
    plt.figure()
    distances_tot, values_tot, pairs = [], [], []
    for i in filtered_orthologs.keys():
        dist, val = filtered_orthologs[i]
        distances_tot.append(dist)
        values_tot.append(val)
        pairs.append(i)

    bin_contents = [[] for _ in range(len(bins))]
    bin_values = [[] for _ in range(len(bins))]

    indices = np.digitize(distances_tot, bins)
    for c, i in enumerate(indices):
        bin_contents[i-1].append(pairs[c])
        bin_values[i-1].append(values_tot[c])
    average_val, average_dist, average_no = [], [], []
    tot_yerr_pos, tot_yerr_neg = [], []
    samples = 30
    
    bin_average = [np.mean(i) for i in bin_values]
    print("the bin average was: ", bin_average)
    plt.plot(bins, bin_average, 'o', color='grey', alpha=0.4)
    plt.figure()
    with open("../200209_orthologs_ints.csv", 'w') as f:
        f.write("bin_num,sample,values,distances,num_samps,error1,error2\n")
        for s in range(samples):
            random_val, random_dist, random_no = [], [], []
            yerr_pos, yerr_neg = [], []
            for i in bin_contents:
                #new_dic = dict(filtered_orthologs)
                new_dic = {}
                for j in i:
                    new_dic[j] = filtered_orthologs[j]
                values, distances = pick_nonoverlap(tree, new_dic)
                random_val.append(1 - np.mean(values))
                random_dist.append(np.mean(distances))
                random_no.append(len(values))
                errors = clopper_pearson(np.sum(values), random_no[-1], alpha=0.32)
                yerr_pos.append(1 - errors[1])
                yerr_neg.append(1 - errors[0])
                line = str(bin_contents.index(i)) + "," + str(s) + "," + str(np.mean(values)) + "," + str(np.mean(distances)) + "," + str(len(values)) + "," + str(errors[0]) + "," + str(errors[1]) +'\n'
                f.write(line)
        f.close()


        plt.plot(bins, random_val, 'o', color='grey', alpha=0.01)

        tot_yerr_pos.append(yerr_pos)
        tot_yerr_neg.append(yerr_neg)  
        average_val.append(random_val)
        average_dist.append(random_dist)
        average_no.append(random_no)

    distances_tot, values_tot, pairs = [], [], []
    pardict = {}
    for i in matched_dic.keys():
        if len(matched_dic[i]) > 1:
            for j in matched_dic[i]:
                try:
                    anc = Phylo.BaseTree.TreeMixin.common_ancestor(tree, list(j)[0], list(j)[1])
                    if data_dic[frozenset([anc])] == 1 and data_dic[frozenset([list(j)[0]])] and data_dic[frozenset([list(j)[1]])] == 1:
                        val = data_dic[j]
                        if val == 0.0 or val == 1.0:
                            dist = Phylo.BaseTree.TreeMixin.distance(tree, list(j)[0], list(j)[1])
                            distances_tot.append(dist)
                            values_tot.append(val)
                            pairs.append(j)
                            pardict[j] = [dist, val]
                except:
                    pass
        # distances_tot.append(dist)
        # values_tot.append(val)
        # pairs.append(i)

    bin_contents = [[] for _ in range(len(bins))]
    bin_values = [[] for _ in range(len(bins))]

    indices = np.digitize(distances_tot, bins)
    for c, i in enumerate(indices):
        bin_contents[i - 1].append(pairs[c])
        bin_values[i - 1].append(values_tot[c])
    average_val, average_dist, average_no = [], [], []
    tot_yerr_pos, tot_yerr_neg = [], []
    samples = 30

    bin_average = [np.mean(i) for i in bin_values]
    print("the bin average was: ", bin_average)
    plt.plot(bins, bin_average, 'o', color='grey', alpha=0.4)
    plt.figure()
    with open("../200209_paralogs_ints.csv", 'w') as f:
        f.write("bin_num,sample,values,distances,num_samps,error1,error2\n")
        for s in range(samples):
            random_val, random_dist, random_no = [], [], []
            yerr_pos, yerr_neg = [], []
            for i in bin_contents:
                # new_dic = dict(filtered_orthologs)
                new_dic = {}
                for j in i:
                    new_dic[j] = pardict[j]
                values, distances = pick_nonoverlap(tree, new_dic)
                random_val.append(1 - np.mean(values))
                random_dist.append(np.mean(distances))
                random_no.append(len(values))
                errors = clopper_pearson(np.sum(values), random_no[-1], alpha=0.32)
                yerr_pos.append(1 - errors[1])
                yerr_neg.append(1 - errors[0])
                line = str(bin_contents.index(i)) + "," + str(s) + "," + str(np.mean(values)) + "," + str(
                    np.mean(distances)) + "," + str(len(values)) + "," + str(errors[0]) + "," + str(errors[1]) + '\n'
                f.write(line)
        f.close()
    xs = []
    for c, i in enumerate(bins):
        xs.append(i + (bins[c] - i)/2.0)
    result = scipy.optimize.leastsq(rev_firs_order_error, [0.2, 0.2], args=(xs[:8], np.mean(average_val, axis=0)[1:9])) #, np.mean(average_no, axis=0)[1:9]))
    print("the result was: ", result)
    print("the sum was: ", sum(rev_firs_order_error(result[0], xs[:8], np.mean(average_val, axis=0)[1:9])))  # , np.mean(average_no, axis=0)[1:9])))
    plt.plot(xs, rev_firs_order(np.array(xs), result[0][0])) #, result[0][1]))
    plt.errorbar(xs, np.mean(average_val, axis=0), yerr=[np.mean(average_val, axis=0) - np.nanmean(tot_yerr_neg, axis=0), np.nanmean(tot_yerr_pos, axis=0) - np.mean(average_val, axis=0)], marker='o', linestyle='None')
    plt.xlabel('Branch length')
    plt.ylabel("Percent of interactions lost")
    plt.figure()
    plt.errorbar(np.mean(average_dist, axis=0), np.mean(average_no, axis=0), yerr=np.std(average_val, axis=0))
    plt.show()

def paralog_ancestors(matched_dic, gene_tree):
    """Finds the ancestors of proteins that exist in the same organism
    Takes:
        matched_dic: a dictionary of key = species name, values = list of frozen sets of protein pairs that existed in that species
        gene_tree: phylo.basetree object of the gene relation tree
    Returns:
        paralog_anc_list: a list containing a frozenset of proteins and their LCA"""
    paralog_anc_list = []
    for i in matched_dic:
        for j in matched_dic[i]:
            if len(j) > 1:
                lca = Phylo.BaseTree.TreeMixin.common_ancestor(gene_tree, [list(j)[0], list(j)[1]])
                paralog_anc = [frozenset([list(j)[0], list(j)[1]]), lca]
                if paralog_anc not in paralog_anc_list:
                    paralog_anc_list.append(paralog_anc)
    return paralog_anc_list


def find_matched_ancestors(species_dic, species_tree, tree, duplication_nodes):
    """ Finds all pairs of time matched ancestors and pairs of extant paralogs within the same species.
    Takes:
        species_dic: a dictionary of species and their proteins in the form key = species, value = list of proteins from that species
        species_tree: a tree object containing the species relationships
        tree: a tree object containing the protein interacitons
        duplication nodes: a list containing clades where duplications happened
    Returns:
        matched_ancestors: a list of frozensets containing time matched ancestral nodes
        matched_dic: a dictionary of species to a list of proteins pairs in those species."""
    speciation_dic, matched_dic = {}, {}  # speciation_dic is a dictionary where the key is the node on the species tree and the items are all nodes on my gene tree that correspond to it.
    matched_ancestors, LCAs = [], []
    total = 0

    # Build the speciation_dic by looping through combinations of proteins from different combinations of species and
    # finding their ancestor, then finding the corresponding ancestral species and add that protein to their ancestral list
    species_combs = list(itertools.combinations(species_dic.keys(), 2))
    for i in species_combs:
        combs = [(x, y) for x in species_dic[i[0]] for y in species_dic[i[1]]]  # pick two proteins from our different species
        for k in combs:
            ancestor = Phylo.BaseTree.TreeMixin.common_ancestor(tree, [k[0], k[1]])
            if ancestor not in LCAs and ancestor not in duplication_nodes:  # if the ancestor isn't a duplication node and we haven't seen it before
                speciationnode = Phylo.BaseTree.TreeMixin.common_ancestor(species_tree, [i[0], i[1]])  # get the node on the species tree
                try:
                    speciation_dic[speciationnode].append(ancestor)  # You have a speciation node, associate with it the particular ancestral protein that species had
                except KeyError:
                    speciation_dic[speciationnode] = [ancestor]
                LCAs.append(ancestor)

    # Populate matched_dic with extant protein pairs from the same species, then from pairs from species in speciation_dic
    for i in species_dic.keys():
        if len(species_dic[i]) == 1:
            try:
                clade = list(tree.find_clades(species_dic[i][0]))[0]
                matched_dic[i] = [frozenset([clade])]
            except KeyError:
                pass
        else:
            clades = [list(tree.find_clades(x))[0] for x in species_dic[i]]
            matched_dic[i] = [frozenset(x) for x in list(itertools.combinations(clades, 2))]
    for i in speciation_dic.keys():
        if len(speciation_dic[i]) == 1:
            matched_dic[i] = [frozenset(speciation_dic[i])]
        else:
            pairwise_combs = list(itertools.combinations(speciation_dic[i], 2))
            matched_dic[i] = [frozenset(x) for x in list(itertools.combinations(speciation_dic[i], 2))]
            for x in pairwise_combs:
                matched_ancestors.append(frozenset(x))
            total += len(pairwise_combs)

    print('There are ' + str(total) + ' time-matched pairs of ancestors on this tree')
    return matched_ancestors, matched_dic




def rev_firs_order_error(params, t, values):
    kf, kb = params
    weights = [1] * len(t)
    model = np.array([rev_firs_order(i, kf) for i in t])
    error = []
    for c, i in enumerate(weights):
        error.append((values[c]-model[c])**2 * i)
    return error   


def rev_firs_order(t, kf):
    """reversible first order model"""
    a = 1 - np.exp(-1 * kf * t)
    return a

def get_species(tips):
    """
    Figures out what species we have in our tree. Should possibly convert to tip_dic
    Takes:
        tips: list of tips from clades
    Returns:
         species_dic: a dictionary of species names of the form key = species name, values = list of proteins from that species
    """
    species_dic = {}
    clades = ['VBP', 'HLF_insects', 'HLF_', 'HLF', 'TEF', 'E4BP4', 'DBP', 'Par1', 'Par2', 'PAR', 'New', 'Weird']
    clades = [c.lower() for c in clades]
    for i in tips:
        i = str(i)
        found = False
        numeral = ''
        if i[-1].isdigit():
            numeral = i[-1]
            i = i[:-1]
        for j in clades:
            if j in i:
                try:
                    species_dic[i[len(j):]].append(i + numeral)
                    found = True
                except KeyError:
                    species_dic[i[len(j):]] = [i + numeral]
                    found = True
                break
        if not found:
            try:
                species_dic[i].append(i + numeral)
            except KeyError:
                species_dic[i] = [i + numeral]
    return species_dic

def clopper_pearson(k, n, alpha=0.32):
    """
    http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    alpha confidence intervals for a binomial distribution of k expected successes on n trials
    Clopper Pearson intervals are a conservative estimate.
    """
    lo = scipy.stats.beta.ppf(alpha/2, k, n-k+1)
    hi = scipy.stats.beta.ppf(1 - alpha/2, k+1, n-k)
    return [lo, hi]


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


if __name__ == "__main__":
    mytree, mytips, mynodes = read_tree('../190918_EA_tree1only.txt')
    mydata_dic = load_experiment('../200120_most_aliases.csv')
    newdata_dic = convert_data_dic_to_species_names(mydata_dic, mynodes, mytips)
    myduplication_nodes = find_duplication_nodes(mytree, '../Gene_duplications.txt')
    myspecies = get_species(mytips)
    myspecies_tree, species_tips, species_nodes = read_tree('../speciestree_node_names.newick.nwk')
    mymatched_ancestors, mymatched_dic = find_matched_ancestors(myspecies, myspecies_tree, mytree, myduplication_nodes)
    # orthos =find_orthologs(tree, data_dic, tip_dic, duplication_nodes)
    # nonovers = pick_nonoverlap(tree,orthos,tip_dic)
    neutral_evo(mytree, newdata_dic, myduplication_nodes, mymatched_dic)
    # print(orthos)
    # print(duplication_nodes)
    
    # print(data_dic)
    # print(tips)
    # for i in range(10):
    #     print("Heres some data",list(data_dic)[i], " = ", data_dic[list(data_dic)[i]])
    # print("length of dict", len(data_dic))
    # aligns = read_alignment('../reformatted_alignment.phy', [25,66])
    # print(aligns.keys())
#    ancestral_seqs,altall_seqs,pp=read_ancestral_seqs('../MLtree', silent = False)  
#    print("ancestral seqs: ", ancestral_seqs)
#    print("altall_seqs: ", altall_seqs)
#    print("pp: ", pp)
#    alias = get_alias_dic()
#    print("\n alias is: ", alias