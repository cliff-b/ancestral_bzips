# -*- coding: utf-8 -*-
import Bio.Phylo as Phylo
from Bio.Phylo import BaseTree
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize
import glob
import itertools
from matplotlib import cm
import math
import copy
from scipy.interpolate import CubicSpline

def read_tree(treefile, silent=True):
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
    return tree, tips, nodes


def load_experiment(filename):
    """Loads a file in .csv format with interaction scores in column
    twelve. Makes a dictionary of sets of the peptides to link to the
    interaction score. Theoretically used to fix aliasing present both in names and hidden from names'
    Takes:
        filename: string, of a file with peptides in cols 1, 2 and interactions in 12
    returns:
         data_dic:   a dictionary of frozen sets of strings, key =frozenset({'pep1', 'pep2'}) value = interaction score
         """
    f = open(filename, 'r')
    # alias_dic= {}#get_alias_dic()
    lines = f.readlines()  # [0].split('\r')
    data_dic = {}
    #    print("so many lines", lines[1:10])
    for i in lines[1:]:

        line = i.split()[0].split(',')

        if len(line[0].split('+')) > 1 or len(line[1].split('+')) > 1:
            for j in line[0].split('+'):
                for k in line[1].split('+'):
                    #
                    #                        if j.lower() in  alias_dic.keys() and k.lower() in  alias_dic.keys():
                    #
                    #                            for l in alias_dic[j.lower()]:
                    #                                for m in alias_dic[k.lower()]:
                    #                                    try:
                    #                                        data_dic[frozenset([l,m])]= np.average(
                    #                                        [data_dic[frozenset([l,m])],float(line[12])])
                    #                                    except KeyError:
                    #                                        data_dic[frozenset([l,m])]=float(line[12])
                    #                        elif j.lower() in  alias_dic.keys():
                    #                            for l in alias_dic[j.lower()]:
                    #                                try:
                    #                                    data_dic[frozenset([l,k.lower()])]= np.average(
                    #                                    [data_dic[frozenset([l,line[1].lower()])],float(line[12])])
                    #                                except KeyError:
                    #                                    data_dic[frozenset([l,k.lower()])]=float(line[12])
                    #                        elif line[1].lower() in  alias_dic.keys():
                    #                            for l in alias_dic[line[1].lower()]:
                    #                                try:
                    #                                    data_dic[frozenset([j.lower(),l])]= np.average(
                    #                                    [data_dic[frozenset([j.lower(),l])],float(line[12])])
                    #                                except KeyError:
                    #                                    data_dic[frozenset([j.lower(),l])]=float(line[12])

                    try:
                        data_dic[frozenset([j.lower(), k.lower()])] = np.average(
                            [data_dic[frozenset([j.lower(), k.lower()])], float(line[12])])
                    except KeyError:
                        data_dic[frozenset([j.lower(), k.lower()])] = float(line[12])
        else:
            #            if False:
            #
            #
            #            if line[1].lower() in  alias_dic.keys() and line[0].lower() in  alias_dic.keys():
            #
            #                for j in alias_dic[line[0].lower()]:
            #                    for k in alias_dic[line[1].lower()]:
            #                        try:
            #                            data_dic[frozenset([j,k])]= np.average(
            #                               [data_dic[frozenset([j,k])],float(line[12])])
            #                        except KeyError:
            #                            data_dic[frozenset([j,k])]=float(line[12])
            #            elif line[0].lower() in  alias_dic.keys():
            #                for j in alias_dic[line[0].lower()]:
            #                    try:
            #                        data_dic[frozenset([j,line[1].lower()])]= np.average(
            #                           [data_dic[frozenset([j,line[1].lower()])],float(line[12])])
            #                    except KeyError:
            #                        data_dic[frozenset([j,line[1].lower()])]=float(line[12])
            #            elif line[1].lower() in  alias_dic.keys():
            #                for j in alias_dic[line[1].lower()]:
            #                    try:
            #                        data_dic[frozenset([line[0].lower(),j])]= np.average(
            #                           [data_dic[frozenset([line[0].lower(),j])],float(line[12])])
            #                    except KeyError:
            #                        data_dic[frozenset([line[0].lower(),j])]=float(line[12])

            try:
                data_dic[frozenset([line[0].lower(), line[1].lower()])] = np.average(
                    [data_dic[frozenset([line[0].lower(), line[1].lower()])], float(line[12])])
            except KeyError:
                data_dic[frozenset([line[0].lower(), line[1].lower()])] = float(line[12])
    #    homodimers=[]
    #    for i in data_dic.keys():
    #        if  len(list(i))==1:
    #            homodimers.append(data_dic[i])
    #    plt.figure()
    #    plt.hist(data_dic.values(), bins=50)
    #    plt.figure()
    #    plt.hist(homodimers, bins=20)
    #    plt.show()
    return data_dic


# def get_alias_dic():
#    ancestral_seqs,altall_seqs,pp=read_ancestral_seqs('../MLtree')     
#    #seq_dict=read_alignment('../Combined.phy',[90,131])
#    seq_dict=read_alignment('../reformatted_alignment.phy',[25,66])
#
#    alias_dic={}        
#    for i in  ancestral_seqs.keys():
#        for j in  seq_dict.keys():
#                    if ancestral_seqs[i] ==seq_dict[j]:
#                        try:    
#                            alias_dic[i].append(j.lower())
#                        except KeyError:
#    
#                             alias_dic[i.lower()]=[j.lower()]
#                        try:    
#                            alias_dic[j.lower()].append(i.lower())
#                        except KeyError:
#    
#                             alias_dic[j.lower()]=[i.lower()]
#                             
#    for i in alias_dic.keys():
#        print (i, alias_dic[i])
#    return alias_dic


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
        duplication_nodes.append(Phylo.BaseTree.TreeMixin.common_ancestor(tree, [i.split(',')[0], i.split(',')[1]]))
    return duplication_nodes


def find_orthologs(tree, data_dic, tip_dic, duplication_list):
    """Takes a tree and interaction data and gives back a list of orthologs
        importantly these orthologs are only ones where the homodimers of the two
        interact and the ancestral proteins homodimers interact
    Takes:
        data_dic: a dictionary frozensets to floats of interaction data in the form key = frozenset({'173', '174'}),
                    value = 1.0
        tip_dic: dictionary of strings to strings to convert names from lowercase in the data_dic to
                    mixedcase in the tree
        duplication_list: a list of clade objects that were points of gene duplication
    Returns:
        ortho_dic: a dictionary of key = frozenset and value = list of two numbers, distance between
                    orthologs and whether they interact or not """
    print('Finding orthologs')
    ortho_dic = {}
    yes, no, missingpoints = 0, 0, 0
    for i in data_dic.keys():
        if len(list(i)) > 1:
            protein1, protein2 = list(i)
            ancestor = Phylo.BaseTree.TreeMixin.common_ancestor(tree, tip_dic[protein1], tip_dic[protein2])
            # This is a gross way to get from tree to tip_dic name
            anc_in_data_dic = frozenset([list(tip_dic.keys())[list(tip_dic.values()).index(ancestor)]])
            try:
                # check that the homodimers for protein1, protein2 and ancestor all homodimerize
                if ancestor not in duplication_list and data_dic[anc_in_data_dic] == 1 and data_dic[frozenset([protein1])] == 1 and data_dic[frozenset([protein2])] == 1:
                    if data_dic[i] == 1:
                        ortho_dic[i] = [Phylo.BaseTree.TreeMixin.distance(tree, tip_dic[protein1], tip_dic[protein2]), 1]
                        yes += 1
                    elif data_dic[i] == 0:
                        ortho_dic[i] = [Phylo.BaseTree.TreeMixin.distance(tree, tip_dic[protein1], tip_dic[protein2]), 0]
                        no += 1
            except KeyError:
                missingpoints += 1
    print('Done finding orthologs', '\nNumber of noninteractors: ', no, '\nNumber of interactors:', yes)
    print("Missing points: ", missingpoints)
    return ortho_dic


def pick_nonoverlap(pairs):
    """finds a set of proteins from the same tree that do not share any branches on the tree and records their
    interaction scores and branch length the algorithm takes a random pair, records values, and removes those that
     overlap with it, then takes a new pair
    Takes:
        pairs: a dictionary of frozensets of protein pairs on the tree to lists of two numerics and a nested list of
        numerics, in the form of key = frozenset({'173', '174'}), values = [total distance between the two proteins in
        the key,  whether they interact, [trace of clade objects between the two proteins on a tree]]
        ie = sample[
    Returns:
        values:
            a list of numerics of whether the proteins interacted or not
        distances:
            a list of numerics of how far apart the proteins were
    """
    values, distances = [], []
    keys = list(pairs.keys())

    while len(keys) > 0:
        key = np.random.choice(keys)
        values.append(pairs[key][1])
        distances.append(pairs[key][0])
        # Here we test for overlap in path. if i has any shared nodes with key then we delete it from the dictionary
        for i in keys:
            try:
                if len(intersection(pairs[i][2], pairs[key][2])) > 0 and i != key:
                    del pairs[i]
            except KeyError:
                pass
        del pairs[key]
        keys = list(pairs.keys())
    return values, distances


def tip_dict(tip_list, node_list):
    """this is necessary to deal with the annoying formatting of the names
    takes in a list of clade objects of the tips of a tree, and either converts the name to
    lower case if it's a tip or a number if its the node
    Takes:
        tip_list: a list of clade objects that are tips  Clade(branch_length=0.041194, name='HLFCryptotermes_secundus')
        node_list: a list of clade objects that are not tips ie Clade(branch_length=0.025239)
    Returns:
        name_dic: a dictionary that converts tips between key: uncapitalized and value: capitalized
        and nodes between their key: number and value: clade object
        NB: why do we have a dictionary with values that are two different objects?"""
    name_dic = {}
    for i in tip_list:
        name_dic[str(i).lower()] = str(i)
    for c, i in enumerate(node_list):
        name_dic[str(c + 172)] = i
    return name_dic


def neutral_evo(tree, data_dic, tip_dic, duplication_nodes, samples=20):
    """Major function that creates three plots for how orthologs continue interacting
    over time. Takes a tree and finds the orthologs that homodimerize, divides them into
    different bins based on their cumulative branch length and then measures the average of each bin
    Then uses a sampling technique to get error-bars on these averages, by sampling orthologs with non-overlapping
    paths. Takes the odds of interacting from each of these 20 times and then averages them. It then produces a
    single term exponential model fit with least squares to the averages.
    Takes:
        tree: a tree object
        data_dic:a dictionary frozensets to floats of interaction data in the form key = frozenset({'173', '174'}),
                    value = 1.0
        tip_dic: dictionary of strings to strings to convert names from lowercase in the data_dic to
                mixedcase in the tree
        duplication_nodes: a list of clade objects that were points of gene duplication
        samples: Number of times to run the sampling. 20 seems to be fairly constant
    Returns:
        Nothing. It does make three graphs however
        Graph 1: Average percent of orthologs still interacting
        Graph 2: Average percent of interactions lost, with clopper pearson errorbars and fit with a simple exp model
        Graph 3: Number of orthologs per bin per sample
    """
    bins = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2, 2.5, 3]
    distances_tot, values_tot, pairs = [], [], []
    bin_contents, bin_values = [], []
    average_val, average_dist, average_no = [], [], []
    tot_yerr_pos, tot_yerr_neg = [], []

    filtered_orthologs = find_orthologs(tree, data_dic, tip_dic, duplication_nodes)
    for i in filtered_orthologs.keys():
        dist, val = filtered_orthologs[i]
        distances_tot.append(dist)
        values_tot.append(val)
        pairs.append(i)
        a = tip_dic[list(i)[0]]
        b = tip_dic[list(i)[1]]
        filtered_orthologs[i].append(Phylo.BaseTree.TreeMixin.trace(tree, a, b) + [a])

    for i in range(len(bins)):
        bin_contents.append([])
        bin_values.append([])
    indices = np.digitize(distances_tot, bins) - 1
    for c, i in enumerate(indices):
        bin_contents[i].append(pairs[c])
        bin_values[i].append(values_tot[c])

    bin_average = [np.mean(i) for i in bin_values]
    print("the bin average was: ", bin_average)
    plt.plot(bins, bin_average, 'o', color='grey')

    for s in range(samples):
        random_val, random_dist, random_no = [], [], []
        yerr_pos, yerr_neg = [], []
        for i in bin_contents:
            new_dic = {}
            for j in i:
                new_dic[j] = filtered_orthologs[j]
            values, distances = pick_nonoverlap(new_dic)
            random_val.append(1 - np.mean(values))
            random_dist.append(np.mean(distances))
            random_no.append(len(values))
            errors = clopper_pearson(np.sum(values), random_no[-1], alpha=0.55)
            yerr_pos.append(1 - errors[1])
            yerr_neg.append(1 - errors[0])
        tot_yerr_pos.append(yerr_pos)
        tot_yerr_neg.append(yerr_neg)
        average_val.append(random_val)
        average_dist.append(random_dist)
        average_no.append(random_no)

    xs = []
    for i in range(len(bins[:-1])):
        xs.append((bins[i] + bins[i + 1]) / 2.0)
    result = scipy.optimize.leastsq(rev_firs_order_error, np.array([0.2]),
                                    args=(xs[:8], np.mean(average_val, axis=0)[:8]))
    print("the result was: ", result)
    print("the sum was: ", sum(
        rev_firs_order_error(result[0], xs[:8], np.mean(average_val, axis=0)[:8])))
    plt.figure()
    plt.plot(xs, rev_firs_order(np.array(xs), result[0][0]))
    plt.errorbar(bins, np.mean(average_val, axis=0),
                 yerr=[np.mean(average_val, axis=0) - np.nanmean(tot_yerr_neg, axis=0),
                       np.nanmean(tot_yerr_pos, axis=0) - np.mean(average_val, axis=0)], marker='o', linestyle='None')
    plt.xlabel('Branch length')
    plt.ylabel("Percent of interactions lost")
    plt.figure()
    plt.errorbar(np.mean(average_dist, axis=0), np.mean(average_no, axis=0), yerr=np.std(average_val, axis=0))

    plt.show()


def rev_firs_order_error(kf, t, values):
    """Computes the least squared error of an array of values which
    it fits to a first order model
    Takes:
        kf: a numeric for the exponential term of the model
        t: an array of n time points for the model's x-axis
        values: true values for the y axis
    Returns:
        error: an 1d array of n values of the least squared error between the model and the data
    """

    model = np.array([rev_firs_order(i, kf) for i in t])
    error = []
    for c in range(len(values)):
        error.append((values[c] - model[c]) ** 2)
    error = np.ravel(error)
    return error


def rev_firs_order(t, kf):
    """Simple function takes time t and exponent kf and plugs them into
     a function that is a reversible first order model
     Takes:
        t: a time or general x-axis point
        kf: the exponential term
    Returns:
        a: the value of y in the function"""
    a = 1 - np.exp(-1 * (kf * t))
    return a


def clopper_pearson(k, n, alpha=0.32):
    """
    http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    alpha confidence intervals for a binomial distribution of k expected successes on n trials
    Clopper Pearson intervals are a conservative estimate.
    Takes:
        k = integer, number of k successes
        n = integer, number of n trials
        alpha = two sided confidence
    Returns:
        a two element list with the form [lower estimate, higher estimate]
    """
    lo = scipy.stats.beta.ppf(alpha / 2, k, n - k + 1)
    hi = scipy.stats.beta.ppf(1 - alpha / 2, k + 1, n - k)
    return [lo, hi]


def intersection(lst1, lst2):
    """Returns elements that appear in both lists
    Takes:
        lst1, lst2: lists of arbitrary objects
    Returns:
        lst3: list of arbitrary objects"""
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def sort_list(list1, list2):
    zipped_pairs = zip(list2, list1)

    z = [x for _, x in sorted(zipped_pairs)]

    return z

# def read_ancestral_seqs(treename, silent=True):
#    '''Legacy code as far as I can tell
#        Takes:
#            treename: a folder containing various tree information
#    '''
#    ancestral_seqs={}
#    altall_seqs={}
#    nodefiles=glob.glob(treename+'/tree1/*.dat') #list of Lazarus node files for that tree
#    pp={}
#    pp_persite=[]
#    ASRtree,tips, nodes=read_tree(treename+'/tree1/tree2') #get tree from Lazarus with reoptimized BLs, currently I have to manually make this file because lazarus formats it stupidly
#    ASRtree=ASRtree.as_phyloxml() # convert to format in which I can edit node labels
#    #print (ASRtree.format('newick'))
#    for i in nodefiles:
#        f=open(i, 'r')
#        lines=f.readlines()
#        nodenumber=i.split('\\')[-1].split('.')[0][4:]
#        ancestral_seqs[nodenumber]=''
#        altall_seqs[nodenumber]=''
#        pp_perfile=[]
#
#        for j in lines:
#            if int(j.split()[0])>25 and int(j.split()[0])<67: #only use bzip domain
#                ancestral_seqs[nodenumber]+=j.split()[1]
#                pp_perfile.append(float(j.split()[2]))
#                try:
#                    if float(j.split()[4])>=0.2:
#                        altall_seqs[nodenumber]+=j.split()[3]
#                    else:
#                        altall_seqs[nodenumber]+=j.split()[1]
#                except IndexError:
#                    altall_seqs[nodenumber]+=j.split()[1]
#        if ancestral_seqs[nodenumber]==altall_seqs[nodenumber]:
#            del altall_seqs[nodenumber]
#        pp[nodenumber]=np.average(pp_perfile)
#        pp_persite.append(pp_perfile)
#
#    if silent==False:
#        print ('There are '+str(len(set(ancestral_seqs.values())))+' unique ML sequences and '+ str(len(altall_seqs))+' altall sequences on '+ treename)
#        plt.figure()
#        plt.imshow(np.array([x for _,x in sorted(zip(pp,pp_persite))]))
#        plt.show()
#
#
#    f=open(treename+'/tree1/Ancestral_seqs.fasta','w')
#    for i in ancestral_seqs.keys():
#        f.writelines('>'+i+'\n')
#        f.writelines(ancestral_seqs[i]+'\n')
#    f.close()
#
#    return [ancestral_seqs,altall_seqs,pp]


def get_species(tips):
    """
    Figures out what species we have in our tree. Should possibly convert to tip_dic
    Takes:
        tips: list of tips from clades
    Returns:
         species_dic: a dictionary of species names of the form key = species name, values = list of proteins from that species
    """
    species_dic = {}
    clades = ['VBP', 'HLF_insects', 'HLF_', 'HLF', 'TEF', 'E4BP4', 'DBP', 'PAR', 'Par1', 'Par2', 'New', 'Weird']
    # A = str(tips[1])
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


def find_matched_extant(species_dic):
    matched_extant = []
    for i in species_dic.keys():
        combs = list(itertools.combinations(species_dic[i], 2))
        for j in combs:
            matched_extant.append(j)
    print('There are ' + str(len(matched_extant)) + ' extant paralog pairs on this tree')
    return matched_extant


def find_partial_orthologs(data_dic, matched_ancestors, matched_extant, tree, nodes):
    print('Finding partial orthologs')
    # This is currently very stupid - it only uses extant pairs and their ancestors
    partial_ortho_dic = {}
    yes = 0
    no = 0
    filterd_matched = []  # This is to make the code faster. Only keep pairs that don't interact
    for i in matched_extant:
        print(i)
        pair = frozenset([i[0].lower(), i[1].lower()])
        if pair in data_dic.keys() and data_dic[pair] == 0:
            filterd_matched.append(i)

    pair_combs = itertools.combinations(filterd_matched, 2)

    c = 0
    f = open('Partial_orthologs.txt', 'w')
    for i in pair_combs:
        c += 1
        print(c)
        try:  # this is to deal with the find_ancestral_pair function not giving output when the pair is dumb
            ancestral_pair, mismatched_pairs, mismatched_distance, ancestors_descendants = find_ancestral_pair(tree, i[0], i[1], nodes)

            if ancestral_pair in data_dic.keys() and data_dic[ancestral_pair] == 0 and len(list(
                    ancestral_pair)) == 2:  # If you have the ancestral pair in the data, if it didn't interact, and if it was actually a pair, not just a dimer. This last one we could skip
                for d, j in enumerate(mismatched_pairs):
                    if data_dic[j] == 1:
                        partial_ortho_dic[j] = [mismatched_distance[d], 1, ancestors_descendants[d]]
                    elif data_dic[j] == 0:
                        partial_ortho_dic[j] = [mismatched_distance[d], 0, ancestors_descendants[d]]

        except TypeError:
            pass
    print('Done finding partial orthologs')
    print('Writing partial orthologs file')
    for i in partial_ortho_dic.keys():
        p1, p2 = list(i)
        distance, value, ancestors_descendants = partial_ortho_dic[i]
        f.writelines(p1 + ' ' + p2 + ' ' + str(distance) + ' ' + str(value) + ' ' + ancestors_descendants[0][0] + ',' +
                     ancestors_descendants[0][1] + ' ' + ancestors_descendants[1][0] + ',' + ancestors_descendants[1][
                         1] + '\n')
    print('Done writing partial orthologs file')
    return partial_ortho_dic


def read_partial_ortholog_file(tree, nodes, data_dic):
    """If partaial ortholog data is already generated, just read in the file"""
    f = open('Partial_orthologs.txt', 'r')
    lines = f.readlines()
    partial_ortholo_dic = {}
    for i in lines:
        line = i.split()
        trace1 = trace = Phylo.BaseTree.TreeMixin.trace(tree, nodes[(int(line[4].split(',')[1]) - 172)],
                                                        line[4].split(',')[
                                                            0])  # I'm not adding the start node here to the trace for a reason here, because it can be in multiple comparisons without a problem. We care about branches, not nodes, so this is fine
        trace2 = trace = Phylo.BaseTree.TreeMixin.trace(tree, nodes[(int(line[5].split(',')[1]) - 172)],
                                                        line[5].split(',')[
                                                            0])  # The try statments below add two interactions per pair that we get by matching a protein to its mismatched ancestor. I will have to explain this in a figure......
        try:
            distance = Phylo.BaseTree.TreeMixin.distance(tree, nodes[(int(line[4].split(',')[1]) - 172)],
                                                         line[4].split(',')[0])
            partial_ortholo_dic[frozenset([line[5].split(',')[1], line[4].split(',')[0]])] = [distance, data_dic[
                frozenset([line[5].split(',')[1].lower(), line[4].split(',')[0].lower()])], trace1]
        except KeyError:
            print('Dont have it')
        try:
            distance = Phylo.BaseTree.TreeMixin.distance(tree, nodes[(int(line[5].split(',')[1]) - 172)],
                                                         line[5].split(',')[0])
            partial_ortholo_dic[frozenset([line[4].split(',')[1], line[5].split(',')[0]])] = [distance, data_dic[
                frozenset([line[4].split(',')[1].lower(), line[5].split(',')[0].lower()])], trace2]
        except KeyError:
            print('Dont have it')
        partial_ortholo_dic[frozenset([line[0], line[1]])] = [float(line[2]), float(line[3]),
                                                              trace1 + trace2]  # This is the interactionbetween e.g shark E4BP4 with human HLF.
    print(len(partial_ortholo_dic))
    f.close()
    bins = (0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2, 2.5, 3, 3.5)

    distances_tot = []
    values_tot = []
    pairs = []

    for i in partial_ortholo_dic.keys():
        dist, val, trace = partial_ortholo_dic[i]
        distances_tot.append(dist)
        values_tot.append(val)
        pairs.append(i)

    bin_contents = []
    bin_values = []
    for i in range(len(bins) + 1):
        bin_contents.append([])
        bin_values.append([])
    indices = np.digitize(distances_tot, bins)
    for c, i in enumerate(indices):
        bin_contents[i].append(pairs[c])
        bin_values[i].append(values_tot[c])
    average_val = []
    average_dist = []
    average_no = []
    tot_yerr_pos = []
    tot_yerr_neg = []

    samples = 100
    for s in range(samples):
        print('sample', s)
        random_val = []
        random_dist = []
        random_no = []
        yerr_pos = []
        yerr_neg = []
        for i in bin_contents:
            new_dic = {}
            for j in i:
                new_dic[j] = partial_ortholo_dic[j]
            # distances, values = randomly_sample_partial_orthologs(new_dic)
            values, distances = pick_nonoverlap(new_dic)
            random_val.append(np.mean(values))
            random_dist.append(np.mean(distances))
            random_no.append(len(values))
            errors = clopper_pearson(np.sum(values), random_no[-1], alpha=0.32)
            yerr_pos.append(errors[1])
            yerr_neg.append(errors[0])

        average_val.append(random_val)
        average_dist.append(random_dist)
        average_no.append(random_no)
        tot_yerr_pos.append(yerr_pos)
        tot_yerr_neg.append(yerr_neg)
    xs = []
    for c, i in enumerate(bins[:-1]):
        xs.append(i + (bins[c + 1] - i) / 2.0)
    xs.append(bins[-1])

    plt.figure()
    plt.errorbar(xs, np.mean(average_val, axis=0)[1:],
                 yerr=[np.mean(average_val, axis=0)[1:] - np.nanmean(tot_yerr_neg, axis=0)[1:],
                       np.nanmean(tot_yerr_pos, axis=0)[1:] - np.mean(average_val, axis=0)[1:]], marker='o',
                 linestyle='None')
    plt.xlabel('Branch length')
    plt.ylabel('Proportion of pairs that interact again')
    plt.savefig('partial_orthologs.png')
    plt.figure()
    plt.errorbar(np.mean(average_dist, axis=0), np.mean(average_no, axis=0), yerr=np.std(average_val, axis=0))

    plt.show()
    return partial_ortholo_dic


def find_ancestral_pair(tree, pair1, pair2, nodes):
    '''Tries to find the ancetral pair offtwo time matched pairs tkes the tree, two pairs and the nodes of the tree. It returns the ancestral par'''
    protein1, protein2 = list(pair1)
    protein3, protein4 = list(pair2)

    if protein1.isdigit():
        protein1 = (nodes.index(protein1) - 172)
    if protein2.isdigit():
        protein2 = (nodes.index(protein1) - 172)
    if protein3.isdigit():
        protein3 = (nodes.index(protein1) - 172)
    if protein4.isdigit():
        protein4 = (nodes.index(protein1) - 172)

    # Make all pairwise ancestors
    ancestor1 = Phylo.BaseTree.TreeMixin.common_ancestor(tree, [protein1, protein3])
    ancestor2 = Phylo.BaseTree.TreeMixin.common_ancestor(tree, [protein1, protein4])
    ancestor3 = Phylo.BaseTree.TreeMixin.common_ancestor(tree, [protein2, protein3])
    ancestor4 = Phylo.BaseTree.TreeMixin.common_ancestor(tree, [protein2, protein4])

    # Calculate the distance between protein1 and the ancestor with its two potential orthologs.
    distance1 = Phylo.BaseTree.TreeMixin.trace(tree, protein1, ancestor1)
    distance2 = Phylo.BaseTree.TreeMixin.trace(tree, protein1, ancestor2)
    distance3 = Phylo.BaseTree.TreeMixin.trace(tree, protein2, ancestor3)
    distance4 = Phylo.BaseTree.TreeMixin.trace(tree, protein2, ancestor4)
    # the
    if len(distance1) < len(distance2) and len(distance1) > 1 and len(distance1) > 1 and len(distance3) > 1 and len(
            distance4) > 1:
        ancestral_pair = frozenset([str(nodes.index(ancestor1) + 172), str(nodes.index(ancestor4) + 172)])
        mismatched_pairs = [frozenset([protein1.lower(), protein4.lower()]),
                            frozenset([protein2.lower(), protein3.lower()])]
        mismatched_distance = [
            Phylo.BaseTree.TreeMixin.distance(tree, protein1, ancestor1) + Phylo.BaseTree.TreeMixin.distance(tree,
                                                                                                             protein4,
                                                                                                             ancestor4),
            Phylo.BaseTree.TreeMixin.distance(tree, protein2, ancestor4) + Phylo.BaseTree.TreeMixin.distance(tree,
                                                                                                             protein3,
                                                                                                             ancestor1)]
        ancestors_descendants = [
            [[protein1, str(nodes.index(ancestor1) + 172)], [protein4, str(nodes.index(ancestor4) + 172)]],
            [[protein2, str(nodes.index(ancestor4) + 172)], [protein3, str(nodes.index(ancestor1) + 172)]]]
        return ancestral_pair, mismatched_pairs, mismatched_distance, ancestors_descendants

    elif len(distance2) < len(distance1) and len(distance1) > 1 and len(distance1) > 1 and len(distance3) > 1 and len(
            distance4) > 1:
        ancestral_pair = frozenset([str(nodes.index(ancestor2) + 172), str(nodes.index(ancestor3) + 172)])
        mismatched_pairs = [frozenset([protein1.lower(), protein3.lower()]),
                            frozenset([protein2.lower(), protein4.lower()])]
        mismatched_distance = [
            Phylo.BaseTree.TreeMixin.distance(tree, protein1, ancestor2) + Phylo.BaseTree.TreeMixin.distance(tree,
                                                                                                             protein3,
                                                                                                             ancestor3),
            Phylo.BaseTree.TreeMixin.distance(tree, protein2, ancestor3) + Phylo.BaseTree.TreeMixin.distance(tree,
                                                                                                             protein4,
                                                                                                             ancestor2)]
        ancestors_descendants = [
            [[protein1, str(nodes.index(ancestor2) + 172)], [protein3, str(nodes.index(ancestor3) + 172)]],
            [[protein2, str(nodes.index(ancestor3) + 172)], [protein4, str(nodes.index(ancestor2) + 172)]]]
        return ancestral_pair, mismatched_pairs, mismatched_distance, ancestors_descendants


def run_partial_orthologs(tree, tips, data_dic, duplication_nodes, nodes, find_partials=False):
    # duplication_nodes = find_duplication_nodes(tree, '../Gene_duplications.txt') # georg's most recent version of this takes more parameters
    # ancestral_seqs, altall_seqs, pp = read_ancestral_seqs('../MLtree')
    species_dic = get_species(tips)  # functional although perhaps should use tip_dic like everything else
    species_tree, species_tips, species_nodes = read_tree('../speciestree_node_names.newick.nwk')
    matched_ancestors, matched_dic = find_matched_ancestors(species_dic, species_nodes,
                                                                                     species_tree, mytree,
                                                                                     myduplication_nodes, mynodes)
    matched_extant = find_matched_extant(species_dic)
    if find_partials == True:  # Only do this if you have to, the calc is slow
        find_partial_orthologs(data_dic, matched_ancestors,
                               matched_extant)  # find partial orthologs and write partial ortholog file

    read_partial_ortholog_file(tree, nodes, data_dic)


##################################################################################################################
#
# Code specific for paralog tree here
#
##################################################################################################################
def find_xy_max(coordinates):
    """finds the dimensions of the tree for plotting circles"""
    xs = [coordinates[i][0] for i in coordinates.keys()]
    ys = [coordinates[i][1] for i in coordinates.keys()]
    return [max(xs), abs(min(ys))]


def make_resolve_dic():
    tree1, tips1, nodes1 = read_tree('../MLtree/tree1/tree2')
    tree, tips, nodes = read_tree('../gene_tree_names.txt')
    resolve_dic = {}
    for c, i in enumerate(tips1):
        resolve_dic[str(tips[c])] = str(i)
    return resolve_dic


def make_inverse_resolve_dic():
    tree1, tips1, nodes1 = read_tree('../MLtree/tree1/tree2')
    tree, tips, nodes = read_tree('../gene_tree_names.txt')
    inverse_resolve_dic = {}
    for c, i in enumerate(tips1):
        inverse_resolve_dic[str(i)] = str(tips[c])
    return inverse_resolve_dic


def convert_data_dic_to_species_names(data_dic, nodes):
    new_dic = {}
    for i in data_dic.keys():
        val = data_dic[i]
        newkey = []
        for j in list(i):
            if j.isdigit():
                k = nodes[int(j)-172]
                newkey.append(k)
            else:
                j = j.lower()
                newkey.append(j)
        new_dic[frozenset(newkey)] = val
    return new_dic


def convert_single_species_node_to_data_dic(node_name, tree_nodes):
    if str(node_name) == 'Clade':
        ret = tree_nodes.index(node_name)+172
    else:
        ret = node_name[:12]
    return ret


def parabola(point1, point2):
    """"Makes a parabola on the x-axis from two points, and returns it.
    Takes:
        point1: a pair of coordinates [x, y] where both are floats
        point2: a pair of coordinates [x, y] where both are floats
    Returns:
        y: a list of the x coordinates of the parabola
        x: a list of the y coordinates of the parabola
    """

    point3 = [max([point1[0], point2[0]]) + (abs((point1[1] - point2[1])/2)**2), np.average([point1[1], point2[1]])]

    x = np.linspace(min(point1[1], point2[1], point3[1]), max(point1[1], point2[1], point3[1]), 20)
    xs = sort_list([point1[0], point2[0], point3[0]], [point1[1], point2[1], point3[1]])
    ys = [point1[1], point2[1], point3[1]]
    ys.sort()
    cs = CubicSpline(ys, xs, axis=1)
    y = cs(x)
    return y, x


def match_nodes(resolve_dic):
    tree1, tips1, nodes1 = read_tree('../MLtree/tree1/tree2')
    tree, tips, nodes = read_tree('../gene_tree_names_resolve.txt.nwk')
    tiplen = len(tips) + 1
    node_dict = {}
    realtips = set([])
    for i in tips1:
        realtips.add(str(i))

    for i in nodes1:
        tips1 = Phylo.BaseTree.TreeMixin.get_terminals(i)
        set1 = set([])
        for j in tips1:
            set1.add(str(j))
        for j in nodes:
            tips2 = Phylo.BaseTree.TreeMixin.get_terminals(j)
            set2 = set([])
            for k in tips2:
                try:
                    set2.add(resolve_dic[str(k)])
                except KeyError:
                    set2.add(str(k))

            if set1 & set2 == set1 and (set2 - set1) & realtips == set([]):
                if str(nodes1.index(i) + 172) not in node_dict.values():
                    node_dict[str(nodes.index(j) + tiplen)] = str(nodes1.index(i) + 172)
                    # forgive us our sins lord
                elif Phylo.BaseTree.TreeMixin.count_terminals(j) < Phylo.BaseTree.TreeMixin.count_terminals(
                        nodes[int(list(node_dict.keys())[
                                      list(node_dict.values()).index(str(nodes1.index(i) + 172))]) - tiplen]):
                    del node_dict[list(node_dict.keys())[list(node_dict.values()).index(str(nodes1.index(i) + 172))]]
                    node_dict[str(nodes.index(j) + tiplen)] = str(nodes1.index(i) + 172)
    print(len(nodes), len(nodes1), len(
        node_dict))  # This is odd. I have 230 nodes in the resolution tree and 170 on the normal one. How can the dictionary have 202 entries? Some nodes get matched to two nodes. This fucks up if the earliest branching node is a loss node!
    return node_dict


def get_species_tip_coordinates(tree, nodes, tips):
    """ Takes a tree and it's corresponding tips and nodes to get the cartesian coordinates of each clade.
    Builds a dictionary of clades to coordinates by enumerating the tips and then looping through the nodes and locating
    the nodes in relation to those tips.
    Takes:
        tree: a tree object
        nodes: a list of clade objects that are nodes on the tree
        tips: a list of clade objects that are tips on the tree
    Returns:
        corr: a dictionary of strings of clades to lists of numerics. key = str(clade object), value = [x-coord, y-coord]
    """
    corr = {}
    for c, i in enumerate(tips):
        corr[str(i)] = [Phylo.BaseTree.TreeMixin.distance(tree, nodes[0], i), c * -1]
    for nod in Phylo.BaseTree.TreeMixin.get_nonterminals(tree, order='postorder'):
        xval = min(corr[str(nod[0])][0], corr[str(nod[1])][0]) - 1  # go one to the left of the current left most
        yval = np.mean([corr[str(nod[0])][1], corr[str(nod[1])][1]])  # average between child nodes
        corr[str(nod)] = [xval, yval]
    return corr


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
                    speciation_dic[str(speciationnode)].append(ancestor)  # You have a speciation node, associate with it the particular ancestral protein that species had
                except KeyError:
                    speciation_dic[str(speciationnode)] = [ancestor]
                LCAs.append(ancestor)

    # Populate matched_dic with extant protein pairs from the same species, then from pairs from species in speciation_dic
    for i in species_dic.keys():
        matched_dic[i] = [frozenset(i) for i in list(itertools.combinations(species_dic[i], 2))]
    for i in speciation_dic.keys():
        if len(speciation_dic[i]) == 1:
            matched_dic[str(i)] = [frozenset(speciation_dic[i])]
        else:
            pairwise_combs = list(itertools.combinations(speciation_dic[i], 2))
            matched_dic[str(i)] = [frozenset(x) for x in pairwise_combs]
            for x in pairwise_combs:
                matched_ancestors.append(frozenset(x))
            total += len(pairwise_combs)

    print('There are ' + str(total) + ' time-matched pairs of ancestors on this tree')
    return matched_ancestors, matched_dic


def sort_out_resolved_nodes_and_species(matched_dic):
    """This sorts out the species names from the tree resolution and assocaites names from the resolved tree with nodes from the unresolved tree. Maybe."""
    resolve_dic = make_resolve_dic()
    node_dict = match_nodes(resolve_dic)
    mapping = {}
    for i in matched_dic.keys():
        for c, j in enumerate(matched_dic[i]):
            listi = list(j)
            for d, k in enumerate(listi):
                if str(k) in resolve_dic.keys():
                    listi[d] = resolve_dic[str(k)]
                    mapping[str(resolve_dic[k])] = k
                elif k in node_dict.keys():
                    listi[d] = node_dict[k]
                    mapping[node_dict[k]] = k
                elif str(k).isdigit():  # This i have to do in order to correct for nodes absent from the nomral tree but low enough in number that they are on the other tree.
                    listi[d] = listi[d] + 'n'
            matched_dic[i][c] = frozenset(listi)
    print("here's where we allow break points to let us think about the choices we made")
    return matched_dic, mapping


def to_gene_rec_form(clade):
    types = ['VBP', 'HLF', 'HLF_', 'HLF_insects',  'TEF', 'E4BP4', 'DBP', 'PAR', 'Par1', 'Par2', 'New', 'Weird']  # order of the list matters because I'm a bad programmer. Big
    lastnum = clade[-1:]
    islastnum = False
    thistype = ""
    if lastnum.isdigit():
        islastnum = True
    for i in types:
        if i in clade:
            thistype = i

    if islastnum:
        newname = clade[len(thistype):-1] + thistype + lastnum
    elif thistype == 'HLF_':
        thistype = 'HLF'
        newname = clade[(len(thistype) + 1):] + thistype  # I have no idea why this is necessary but the name of daphnia pulex between the systems needs it
    else:
        newname = clade[len(thistype):] + thistype
    return newname


def order_lookup(list_of_proteins_to_plot, order_list):
    tree, tips, nodes = read_tree('../190918_EA_tree1only.txt')
    tree1, tips1, nodes1 = read_tree('../gene_tree_names_resolve.txt.nwk')
    new_loptp = []
    node_set = set()
    for i in list_of_proteins_to_plot:
        for j in list(i):
            if 'LOST' in str(j):
                j = j.replace('LOST', '')
            node_set.add(j)
    for i in node_set:
        # position = -1
        zpos = 1
        terms, ints = [], []
        node_loc = i
        if type(i) is str and i[0] == 'r':
            pass
        else:
            if not isinstance(i, Phylo.BaseTree.Clade) and i[0] != 'r':
                node_loc = list(tree.find_clades(name = str(i)))[0]
            if node_loc.is_terminal():
                node_loc = to_gene_rec_form(str(node_loc))
                terms.append(node_loc)
            else:
                termspretree = node_loc.get_terminals()
                for clade in termspretree:
                    clade = to_gene_rec_form(str(clade))
                    terms.append(clade)
        if len(terms) > 1:
            new_clade = Phylo.BaseTree.TreeMixin.common_ancestor(tree1,[str(i) for i in terms])
        elif len(terms) == 1:
            new_clade = terms[0]
        elif len(terms) == 0 and i[0] == 'r':
            new_clade = i
        else:
            print("shit, we got no clades")

        trace = Phylo.BaseTree.TreeMixin.trace(tree1, nodes1[0], str(new_clade))
        trace = [str(j) for j in trace]
        ints = intersection(trace, order_list)
        # if len(ints) == 1:
        #     position = order_list.index(ints[0])
        #     for j in trace[trace.index(ints[0]):trace.index(str(new_clade))+1]:
        #         loc_clade = list(tree1.find_clades(name = j))[0]
        #         if len(loc_clade.clades) > 0 and str(loc_clade.clades[1]) in trace:
        #             zpos +=1
        #     new_loptp.append([i, position, zpos])

        if len(ints) == 1:
            position = order_list.index(ints[0])
            for j in new_loptp:
                if j[1] == position:
                    zpos += 1
            new_loptp.append([i, position, zpos])
    return new_loptp


def add_loss_nodes(matched_dic):
    """Takes lost genes and adds them to the parent that would have contained them in matched_dic. The dictionary of
     parents must be built first as biopython doesn't seem to have that function. Data manipulation in circle needs
     "LOST" to be a part of lost proteins names.
    Takes:
        matched_dic: a dictionary of species to the genes that they hold. key = species name/node name, value = list frozensets of pairs of proteins in the key species
    Returns:
         matched_dic: a dictionary of species to the genes that they hold. key = species name/node name, value = list frozensets of pairs of proteins in the key species
         now with lost proteins
    """
    tree, tips, nodes = read_tree('../speciestree_node_names.newick.nwk')
    tree1, tips1, nodes1 = read_tree('../gene_tree_names_resolve.txt.nwk')
    parent_dic = {}
    for i in tree.find_clades():
        trace = Phylo.BaseTree.TreeMixin.trace(tree, nodes[0], i)
        if len(trace) > 2:
            parent_dic[str(i)] = trace[-2]
    # take the resolved gene tree and searches it for lost tips. It adds the lost proteins to the parent of the species that lost them.
    for i in tips1:
        if str(i)[-4:] == 'LOST':
            trace = Phylo.BaseTree.TreeMixin.trace(tree1, nodes1[0], i)
            try:
                matched_dic[str(parent_dic[str(i).split('*')[0]])].append(frozenset([str(trace[-2]) + "LOST"]))
            except KeyError:
                matched_dic[str(parent_dic[str(i).split('*')[0]])] = [frozenset([str(trace[-2]) +"LOST"])]
    return matched_dic


def circle(node, matched_dic, coor, size, ax, order_list, data_dic, treenodes):
    stretch = 0.4
    ext_x, ext_y = size
    scale = 0.010
    x, y = coor
    node_set = set([])
    # for i in matched_dic[node]:
    #     for j in list(i):
    #         if 'LOST' in j and j.split('*')[0] != node:  #
    #             node_set.add(j)
    #         elif 'LOST' not in j:
    #             node_set.add(j)
    # # node_list=list(node_set)
    for i in matched_dic[node]:
        for j in list(i):
            if type(j) is str:
                j = j.lower()
            node_set.add(j)


    node_list = list(node_set)
    for i in order_list:
        if type(i[0]) is str:
            i[0] = i[0].lower()

    if len(node_list) > 0:
        #angle = 360.0 / len(node_list)
        # xs = [np.cos(math.radians(i * angle)) * scale * ext_x + x for i in range(len(node_list))]
        #ys = [np.sin(math.radians(i * angle)) * scale * ext_y * stretch + y for i in range(len(node_list))]
        angle = 360.0 / 5
        xs = [np.cos(math.radians(i * angle)) * scale * ext_x + x for i in range(5)]
        ys = [np.sin(math.radians(i * angle)) * scale * ext_y * stretch + y for i in range(5)]

        #make hexagon outlines
        if order_list:
            numhexs = max(order_list, key=lambda x: x[2])
            xlevs, ylevs = [], []
            for j in range(numhexs[2]):
                xlev1, ylev1 = [], []
                for i in range(len(xs)):
                    xlev1.append(xs[i] + (xs[i] - x)*(j+1))
                    ylev1.append(ys[i] + (ys[i] - y)*(j+1))
                xlev1.append(xs[0] + (xs[0] - x)*(j+1))
                ylev1.append(ys[0] + (ys[0] - y)*(j+1))
                xlevs.append(xlev1)
                ylevs.append(ylev1)
            for j in range(len(xlevs)):
                ax.plot(xlevs[j], ylevs[j], color = 'grey', alpha =0.4, zorder = 2)

        cols = []
        colnams = []
        xpoints, ypoints = [], []
        for i in order_list:
            xpoints.append(xs[i[1]] + (xs[i[1]] - x )* i[2])
            ypoints.append(ys[i[1]] + (ys[i[1]] - y )* i[2])
            nam = convert_single_species_node_to_data_dic(i[0], treenodes)
            colnams.append([str(nam), xs[i[1]] + (xs[i[1]] - x) * i[2], ys[i[1]] + (ys[i[1]] - y ) * i[2]])
            try:
                if data_dic[frozenset([i[0]])] == 1:
                    cols.append((255/256, 49/256, 49/256, 1))#hexcodes for colors. This is red.
                elif data_dic[frozenset([i[0]])] == 0:
                    cols.append((0,145/256,244/256,1))
            except KeyError:
                if i[0][:2] == 'r7': #yes this should be more clear that its a LOST node but i deleted that earlier and they do all start with r7
                    cols.append((1,0.75,0,1))
                else:
                    cols.append(cm.Greys(0.5))
        ax.scatter(xpoints, ypoints, color=cols, s=8, zorder=10, edgecolors='black')
        # for j in range(len(colnams)):
        #     if str(colnams[j][0]).isdigit() and int(colnams[j][0]) > 170 and int(colnams[j][0]) < 177:
        #         ax.text(colnams[j][1]+0.2, colnams[j][2] -0.33, s = colnams[j][0], fontsize=10, color = "red")
        #     else:
        #         #print("text being printed ", colnams[j], " num cols = ", len(colnams))
        #         ax.text(colnams[j][1]+0.2, colnams[j][2] -0.33, s = colnams[j][0], fontsize=5)
        # Draw lines between the nodes all the nodes colored by whether there is an interaction or not
        for c, d in itertools.combinations(range(len(order_list)), 2):
            co_x = [xs[order_list[c][1]] + (xs[order_list[c][1]] - x) * order_list[c][2], xs[order_list[d][1]] + (xs[order_list[d][1]] - x) * order_list[d][2]]
            co_y = [ys[order_list[c][1]] + (ys[order_list[c][1]] - y) * order_list[c][2], ys[order_list[d][1]] + (ys[order_list[d][1]] - y) * order_list[d][2]]
            try:
                if not ((order_list[c][1] == order_list[d][1]) and abs(order_list[c][2] - order_list[d][2]) > 1):
                    # forgive me god
                    if data_dic[frozenset([order_list[c][0], order_list[d][0]])] == 1.0:
                        ax.plot(co_x, co_y, color=(255/256, 49/256, 49/256, 1), zorder=3)
                    elif data_dic[frozenset([order_list[c][0], order_list[d][0]])] == 0.5:
                        ax.plot(co_x, co_y, color=(196/256, 0, 237/256,1), zorder=3)
                    elif data_dic[frozenset([order_list[c][0], order_list[d][0]])] == 0.0:
                        ax.plot(co_x, co_y, color=(0,145/256,244/256,1), zorder=3)
                else:
                    try:
                        xarc, yarc = parabola([co_x[0],co_y[0]], [co_x[1], co_y[1]])
                        if data_dic[frozenset([order_list[c][0], order_list[d][0]])] == 1.0:
                            arccol = (255/256, 49/256, 49/256, 1)
                        elif data_dic[frozenset([order_list[c][0], order_list[d][0]])] == 0.5:
                            arccol = (196/256, 0, 237/256,1)
                        elif data_dic[frozenset([order_list[c][0], order_list[d][0]])] == 0.0:
                            arccol = (0,145/256,244/256,1)
                        else:
                            arccol = (0,1,0,1)
                        ax.plot(xarc, yarc, color=arccol, alpha=1, zorder=4)
                    except:
                        print("something in the parabolas went poorly")

            except KeyError:
                pass
        xleft, xright = ax.get_xlim()
        ybottom, ytop = ax.get_ylim()
        ratio = 2.5
        ax.set_aspect(abs((xright - xleft) / (ybottom - ytop)) * ratio)
        ax.set_axis_off()
    else:
        try:
            val = data_dic[frozenset([str(node).lower()])]
            if val == 1.0:
                ax.scatter(coor[0], coor[1],  color=(255/256, 49/256, 49/256, 1),zorder = 10, s=8, edgecolors='black')
            elif val == 0.0:
                ax.scatter(coor[0], coor[1], color=(0, 145 / 256, 244 / 256, 1), zorder=10, s=8, edgecolors='black')
        except KeyError:
            pass


def plot_species_interactions(species_coordinates, matched_dic, ax, species_tree, tree, tree_nodes, data_dic):
    ext_x, ext_y = find_xy_max(species_coordinates)
    inverse_resolve_dic = make_inverse_resolve_dic()

    # duplication_nodes = find_duplication_nodes(tree, '../Gene_duplications.txt', 0, True, inverse_resolve_dic)
    # duplication_nodes = find_duplication_nodes(tree, '../Gene_duplications.txt', 0, True, inverse_resolve_dic)
    duplication_nodes = find_duplication_nodes(tree, '../Gene_duplications.txt')
    order_dic = {}  # this is used to keep track of the orientation of the dots across node in the interactogram
    # # First let's find the LCA of metazoa
    deepest = Phylo.BaseTree.TreeMixin.common_ancestor(species_tree, ['Drosophila_melanogaster', 'Actinia_tenebrosa'])
    #deepest = str(species_nodes.index(ancestor) + 53)
    #deepest = ancestor
    deepest_set = set()
    for i in matched_dic[str(deepest)]:
        deepest_set.update(list(i))
    order_list = list(deepest_set)  # Now we have a fixed order at the deppest node.
    descendants = Phylo.BaseTree.TreeMixin.get_terminals(deepest)
    done = set([])
    # for i in descendants:
    #     trace = Phylo.BaseTree.TreeMixin.trace(species_tree, deepest, i)
    #     for c, j in enumerate(trace[1:]):
    #         if j not in done:
    #             try:
    #                 order = match_order(order_dic[str(trace[c])], matched_dic[str(j)], tree, tips, nodes, duplication_nodes, mapping)
    #             except ValueError:
    #                 try:
    #                     order = match_order(order_dic[str(trace[c])], matched_dic[str(j)], tree, tips, nodes, duplication_nodes, mapping)
    #                 except:
    #                     print("there is so much failyre", "order_dic ", str(species_nodes.index(trace[c]) + 53), "\nmatched_dic presumably unprintable", "\nc / j: ", c, j)
    #             except KeyError:
    #                 print("Well that shit didn't work\n", "order_dic ", str(trace[c]), "\nmatched_dic ", str(j), "\nc / j: ", c, j)
    #
    #             done.add(j)
    #             try:
    #                 order_dic[str(j)] = order
    #             except ValueError:
    #                 order_dic[str(j)] = order
    #             except:
    #                 print("we're adding nothing to the dic because there isn't anything in matched at 101n")
    #
    # stretch = 1
    tree, tips, nodes = read_tree('../speciestree_node_names.newick.nwk')
    tree1, tips1, nodes1 = read_tree('../gene_tree_names_resolve.txt.nwk')
    cleanorder_list = []
    for j in order_list:
        if type(j) == str:
            j =  j.replace('LOST', '')
            cleanorder_list.append(str(j))
    for j in order_list:
        if not type(j) == str:
            terms = j.get_terminals()
            tree1terms = []
            for k in terms:
                clade = to_gene_rec_form(str(k))
                tree1terms.append(clade)
            j = Phylo.BaseTree.TreeMixin.common_ancestor(tree1, tree1terms)
            cleanorder_list.append(str(j))
    cleanorder_list = ['r7443','r7425','n5243','r7408','n5062',]
    for i in matched_dic.keys():
        #if len(matched_dic[i]) > 0:
        neworder = order_lookup(matched_dic[i], cleanorder_list)
        # for x in neworder:
        #     [x[1] for x in neworder]
        # if len(set([x[2] for [x[1] for x in neworder] in neworder])) != len([x[2] for x in neworder]):
        #     print("fuck")
        coor = species_coordinates[str(i)]
        try:
            circle(i, matched_dic, coor, [ext_x + 1, ext_y + 1], ax, neworder, data_dic, tree_nodes)
        except KeyError:
            circle(i, matched_dic, coor, [ext_x + 1, ext_y + 1], ax, [], data_dic, tree_nodes)
    #plt.savefig('species_interactions.eps')
    return done


def draw_species_tree(input_tree, nodes, tips, coordinates, duplication_nodes, subfig):
    """Draws the tree from a list of coordinates and labels the tips
     Finds each tip and draw a line between them and the parent. Skips if the line has ben already drawn
     The coordinate moving is brain damaging, but seems correct
     Takes:
        input_tree: The tree we're drawing
        nodes: a list of nodes on the tree
        tips: a list of tips on the tree
        coordinates: a dictionary of species names to list of xaxis, yaxis points. key = Species name/ancestral species number, value = [xaxis point, yaxis point]
        duplication nodes: PASS
        subfig: a matplotlib subplot to add bits too
    Returns:
        Nothing
     """
    exclude = set([])
    col = 'black'
    for i in tips:
        coordinate1 = coordinates['n102']  # this should be the most ancestral node
        trace = Phylo.BaseTree.TreeMixin.trace(input_tree, nodes[0], i)
        # Follow the trace from the most ancestral node to tip i and draw each coordinate pair
        for j in trace:
            if j in tips:
                coordinate2 = coordinates[str(j)]
                genus = str(j)[0]
                species = str(j).split("_")[1]
                subfig.text(coordinate2[0] + 1.5, coordinate2[1] - 0.33, genus + ". " + species, fontsize=5)
            else:
                coordinate2 = coordinates[str(j)]
                #subfig.text(coordinate2[0] + .25, coordinate2[1], str(j), fontsize=5)
            if frozenset([frozenset(coordinate1), frozenset(coordinate2)]) not in exclude:
                subfig.plot([coordinate1[0], coordinate1[0]], [coordinate1[1], coordinate2[1]], color=col, zorder=1)
                subfig.plot([coordinate1[0], coordinate2[0]], [coordinate2[1], coordinate2[1]], color=col, zorder=1)
            exclude.add(frozenset([frozenset(coordinate1), frozenset(coordinate2)]))
            if j in duplication_nodes:
                col = 'red'
            else:
                col = 'grey'
            coordinate1 = coordinate2
    return


# def make_interactogram():
#     ancestral_seqs,altall_seqs,pp=read_ancestral_seqs('../MLtree')
#     tree, tips, nodes=read_tree('../gene_tree_names_resolve.txt.nwk')
#     inverse_resolve_dic=make_inverse_resolve_dic()
#     #duplication_nodes=find_duplication_nodes(tree,'../Gene_duplications.txt',0, True, inverse_resolve_dic)
#     duplication_nodes = find_duplication_nodes(tree, '../Gene_duplications.txt')
#     done=set([])
#     #the for loop here renames the LOST taxa to be so they don't have duplicate names
#     for i in tips:
#         if str(i) in done:
#             count=1
#             while str(i)+str(count) in done:
#                 count+=1
#             done.add(str(i)+str(count))
#             i.__setattr__('name',str(i)+str(count))
#         else:
#             done.add(str(i))
#     species_dic=get_species_resolve(tips)
#     matched_ancestors, matched_extant,matched_dic,age_dic=find_matched_ancestors (species_dic, tree, duplication_nodes, nodes,ancestral_seqs, altall_seqs) # This retunrs node numbers in the numbering of the resolution tree.
#     fig=plt.figure()
#     ax = fig.add_subplot(111)
#     species_tree, species_tips, species_nodes=read_tree('../speciestree_node_names.newick.nwk')
#     species_coordinates=get_species_tip_coordinates(species_tree,species_nodes,species_tips)
#     species_tree, species_tips, species_nodes=read_tree('../speciestree_node_names.newick.nwk')
#     draw_species_tree(species_tree,species_nodes,species_tips,species_coordinates,[],ax)
#     matched_dic,mapping=sort_out_res1olved_nodes_and_species(matched_dic)
#     matched_dic=add_loss_nodes(matched_dic)
#     species_tree, species_tips, species_nodes=read_tree('../speciestree_node_names.newick.nwk')
#     plot_species_interactions(species_coordinates,matched_dic,ax,species_tree,species_nodes, tree,tips, nodes,duplication_nodes, mapping)
#     plt.show()


if __name__ == "__main__":
    mytree, mytips, mynodes = read_tree('../190918_EA_tree1only.txt')
    mytip_dic = tip_dict(mytips, mynodes)
    # print(tips)
    # print(nodes)
    # someday we'll get there
    mydata_dic = load_experiment('../200120_most_aliases.csv')
    myduplication_nodes = find_duplication_nodes(mytree, '../Gene_duplications.txt')
    print("dupe nodes", myduplication_nodes)
    species = get_species(mytips)
    print("species are: ", list(species.keys()))
    species_tree, species_tips, species_nodes = read_tree('../speciestree_node_names.newick.nwk')  # contains species relationships
    # print("immediate species tree: ", species_tree)
    # print("immediate species tips: ", species_tips)
    matched_ancestors, matched_dic = find_matched_ancestors(species, species_tree, mytree, myduplication_nodes)
    exts = find_matched_extant(species)
    #e print("post find_matched species tree: ",species_tree)
    # print("post find_matched species tips: ",species_tips)
    # run_partial_orthologs(mytree, mytips,mydata_dic, myduplication_nodes,mynodes)
    # find_partial_orthologs(mydata_dic, matched_ancestors, matched_extant, mytree, mynodes)
    # orthos =find_orthologs(tree, data_dic, tip_dic, duplication_nodes)
    # nonovers = pick_nonoverlap(tree,orthos,tip_dic)
    # neutral_evo(mytree, mydata_dic, mytip_dic, myduplication_nodes)
    species_coordinates = get_species_tip_coordinates(species_tree, species_nodes, species_tips)
    fig = plt.figure(figsize=(10,25))
    ax = fig.add_subplot(111)
    #matched_dic, mapping = sort_out_resolved_nodes_and_species(matched_dic)
    matched_dic = add_loss_nodes(matched_dic)
    newdict = convert_data_dic_to_species_names(mydata_dic, mynodes)#print("late species tree: ", species_tree)
    #

    draw_species_tree(species_tree, species_nodes, species_tips, species_coordinates, myduplication_nodes, ax)
    done = plot_species_interactions(species_coordinates, matched_dic, ax, species_tree, mytree, mynodes, newdict)


    figure = plt.gcf()
    figure.savefig('../200121_full_paralogsclean.pdf', figsize = (20, 50))
    plt.show()
    # for myi in range(10):
    #     print("Here's some data", list(mydata_dic)[myi], " = ", mydata_dic[list(mydata_dic)[myi]])
    # print("length of dict", len(mydata_dic))
    # aligns = read_alignment('../reformatted_alignment.phy', [25, 66])
    # print(aligns.keys())
#    ancestral_seqs,altall_seqs,pp=read_ancestral_seqs('../MLtree', silent = False)  
#    print("ancestral seqs: ", ancestral_seqs)
#    print("altall_seqs: ", altall_seqs)
#    print("pp: ", pp)
#    alias = get_alias_dic()
#    print("\n alias is: ", alias
