# -*- coding: utf-8 -*-
import Bio.Phylo as Phylo
from Bio.Phylo import BaseTree
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize
import glob
import itertools

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
                if ancestor not in duplication_list and data_dic[anc_in_data_dic] == 1 and data_dic[
                        frozenset([protein1])] == 1 and data_dic[frozenset([protein2])] == 1:
                    if data_dic[i] == 1:
                        ortho_dic[i] = [Phylo.BaseTree.TreeMixin.distance(tree, tip_dic[protein1], tip_dic[protein2]),
                                        1]
                        yes += 1
                    elif data_dic[i] == 0:
                        ortho_dic[i] = [Phylo.BaseTree.TreeMixin.distance(tree, tip_dic[protein1], tip_dic[protein2]),
                                        0]
                        no += 1
            except KeyError:
                missingpoints += 1
    print('Done finding orthologs', '\nNumber of noninteractors: ', no, '\nNumber of interactors:', yes)
    print("Missing points: ", missingpoints)
    return ortho_dic


def pick_nonoverlap(tree, pairs, tip_dic):
    """finds a set of proteins from the same tree that do not share any branches on the tree and records their
    interaction scores and branch length the algorithm takes a random pair, records values, and removes those that
     overlap with it, then takes a new pair
    Takes:
        tree: a tree object
        pairs: a dictionary of frozensets of protein pairs on the tree to lists of two numerics,
               in the form of key = frozenset({'173', '174'}), values = [total distance between the two proteins in
               the key,  whether they interact]
        tip_dic: tip_dic: dictionary of strings to strings to convert names from lowercase in the data_dic to
                mixedcase in the tree
    Returns:
        values:
            a list of numerics of whether the proteins interacted or not
        distances:
            a list of numerics of how far apart the proteins were
    """
    for i in pairs.keys():
        a = tip_dic[list(i)[0]]
        b = tip_dic[list(i)[1]]
        pairs[i].append(Phylo.BaseTree.TreeMixin.trace(tree, a, b) + [a])
    values = []
    distances = []
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
            values, distances = pick_nonoverlap(tree, new_dic, tip_dic)
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
    a = 1 - np.exp(-1*(kf*t))
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

def read_ancestral_seqs(treename, silent=True):
   '''Legacy code as far as I can tell
       Takes:
           treename: a folder containing various tree information
   '''
   ancestral_seqs={}
   altall_seqs={}
   nodefiles=glob.glob(treename+'/tree1/*.dat') #list of Lazarus node files for that tree
   pp={}
   pp_persite=[]
   ASRtree,tips, nodes=read_tree(treename+'/tree1/tree2') #get tree from Lazarus with reoptimized BLs, currently I have to manually make this file because lazarus formats it stupidly
   ASRtree=ASRtree.as_phyloxml() # convert to format in which I can edit node labels
   #print (ASRtree.format('newick'))
   for i in nodefiles:
       f=open(i, 'r')
       lines=f.readlines()
       nodenumber=i.split('\\')[-1].split('.')[0][4:]
       ancestral_seqs[nodenumber]=''
       altall_seqs[nodenumber]=''
       pp_perfile=[]

       for j in lines:
           if int(j.split()[0])>25 and int(j.split()[0])<67: #only use bzip domain
               ancestral_seqs[nodenumber]+=j.split()[1]
               pp_perfile.append(float(j.split()[2]))
               try:
                   if float(j.split()[4])>=0.2:
                       altall_seqs[nodenumber]+=j.split()[3]
                   else:
                       altall_seqs[nodenumber]+=j.split()[1]
               except IndexError:
                   altall_seqs[nodenumber]+=j.split()[1]
       if ancestral_seqs[nodenumber]==altall_seqs[nodenumber]:
           del altall_seqs[nodenumber]
       pp[nodenumber]=np.average(pp_perfile)
       pp_persite.append(pp_perfile)

   if silent==False:
       print ('There are '+str(len(set(ancestral_seqs.values())))+' unique ML sequences and '+ str(len(altall_seqs))+' altall sequences on '+ treename)
       plt.figure()
       plt.imshow(np.array([x for _,x in sorted(zip(pp,pp_persite))]))
       plt.show()


   f=open(treename+'/tree1/Ancestral_seqs.fasta','w')
   for i in ancestral_seqs.keys():
       f.writelines('>'+i+'\n')
       f.writelines(ancestral_seqs[i]+'\n')
   f.close()

   return [ancestral_seqs,altall_seqs,pp]


def get_species(tips):
    species_dic = {}
    clades = ['VBP', 'HLF_insects', 'HLF_', 'HLF', 'TEF', 'E4BP4', 'DBP', 'PAR', 'Par1', 'Par2', 'New', 'Weird']
    #A = str(tips[1])
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
        if found == False:
            try:
                species_dic[i].append(i + numeral)

            except KeyError:
                species_dic[i] = [i + numeral]

    return species_dic


def find_matched_ancestors(species_dic, tree, duplication_nodes, nodes):
    ''' Finds all pairs of time matched ancestors and pairs of extant paralogs within the same species. '''

    species_tree, species_tips, species_nodes = read_tree(
        '../speciestree_node_names.newick.nwk')  # contains species relationships
    speciation_dic = {}  # this is a dictionary where the key is the node on the species tree and the items are all nodes on my gene tree that correspond to it.
    LCAs = []
    total = 0
    matched_ancestors = []
    tiplen = Phylo.BaseTree.TreeMixin.count_terminals(tree) + 1
    # pick one species
    species_combs = list(itertools.combinations(species_dic.keys(), 2))
    ASRtree1 = tree.as_phyloxml()  # We l eventually have a tree in which each node is labelled '1' if it's part of a time matched comparison and '0' if it's not.
    for i in ASRtree1.find_clades():  # Set all the node labels of our new tree to 0
        if i.is_terminal() == False:
            i._set_confidence(0)
    matched_dic = {}
    matched_extant = []
    for i in species_dic.keys():
        matched_dic[i] = []
        combs = list(itertools.combinations(species_dic[i], 2))
        for j in combs:
            matched_extant.append(j)
            matched_dic[i].append(frozenset([j[0], j[1]]))
    for i in species_combs:

        combs = [(x, y) for x in species_dic[i[0]] for y in species_dic[i[1]]]  # pick two species

        for k in combs:  # go throug the combinations of their extant sequnces (one from each of the two)
            ancestor = Phylo.BaseTree.TreeMixin.common_ancestor(tree, [k[0], k[1]])  # get the ancestor

            if ancestor not in duplication_nodes and ancestor not in LCAs:

                try:
                    speciationnode = Phylo.BaseTree.TreeMixin.common_ancestor(species_tree,
                                                                              [i[0], i[1]])  # get the speciation node

                    try:
                        speciation_dic[str(species_nodes.index(speciationnode))].append(
                            k)  # You have a speciation node, associate with it the particular ancestor
                    except KeyError:
                        speciation_dic[str(species_nodes.index(speciationnode))] = [
                            k]  # You have a speciation node, associate with it the particular ancestor
                    LCAs.append(ancestor)
                except ValueError:
                    print('one of the species is not on the tree', i)
    age_dic = {}
    for i in speciation_dic.keys():  # Now go through the speciation dictm and  for each entry, make all pairwise combinations (each part of a pair is two sequences that define a LCA)

        pairwise_combs = list(itertools.combinations(speciation_dic[i], 2))
        matched_dic[str(int(i) + 53)] = []
        age_dic[i] = Phylo.BaseTree.TreeMixin.distance(species_tree, species_nodes[int(i)], species_nodes[0])
        if len(speciation_dic[i]) == 1:
            node1 = str(nodes.index(Phylo.BaseTree.TreeMixin.common_ancestor(tree, [speciation_dic[i][0][0],
                                                                                    speciation_dic[i][0][
                                                                                        1]])) + tiplen)  # Now I fucking get the problem. For nodes I don't have I don't correct the numbers. If they are small enough they stay
            matched_dic[str(int(i) + 53)].append(frozenset([node1]))
        for j in pairwise_combs:
            Phylo.BaseTree.TreeMixin.common_ancestor(ASRtree1, [j[0][0], j[0][1]])._set_confidence(
                1)  # Node labels are set to 1 if this node is part of any time-macthed comparison
            Phylo.BaseTree.TreeMixin.common_ancestor(ASRtree1, [j[1][0], j[1][1]])._set_confidence(1)
            node1 = str(nodes.index(Phylo.BaseTree.TreeMixin.common_ancestor(tree, [j[0][0], j[0][
                1]])) + tiplen)  # Get the node correspnding to the first pair of sequences
            node2 = str(nodes.index(Phylo.BaseTree.TreeMixin.common_ancestor(tree, [j[1][0], j[1][
                1]])) + tiplen)  # Get the node correspnding to the second pair of sequences
            matched_ancestors.append(frozenset([node1, node2]))
            matched_dic[str(int(i) + 53)].append(frozenset([node1, node2]))
        total += len(pairwise_combs)

    print('There are ' + str(total) + ' time-matched pairs of ancestors on this tree')
    print('There are ' + str(len(matched_extant)) + ' extant paralog pairs on this tree')

    f = open('../MLtree/tree1/time-matched_nodes.newick', 'w')
    f.writelines(
        ASRtree1.format('newick').replace('):', ')'))  # stupid trick so the tree is formatted right for fig tree
    f.close()
    return matched_ancestors, matched_extant, matched_dic, age_dic


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
            a = 1
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
    ''' If partaial ortholog data is already generated, just read in the file'''
    f = open('Partial_orthologs.txt', 'r')
    lines = f.readlines()
    partial_ortholo_dic = {}
    yes = []
    no = []
    distances = []
    values = []
    done = set([])
    for i in lines:
        line = i.split()
        trace1 = trace = Phylo.BaseTree.TreeMixin.trace(tree, nodes[(int(line[4].split(',')[1]) - 172)],
                                                        line[4].split(',')[
                                                            0])  # I'm not adding the start node here to the trace for a reason here, because it can be in multiple comparisons without a problem. We care about branches, not nodes, so this is fine
        trace2 = trace = Phylo.BaseTree.TreeMixin.trace(tree, nodes[(int(line[5].split(',')[1]) - 172)],
                                                        line[5].split(',')[0])
        # The try statments below add two interactions per pair that we get by matching a protein to its mismatched ancestor. I will have to explain this in a figure......
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

    distances = []
    values = []
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
            distances, values = randomly_sample_partial_orthologs(new_dic)
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

def randomly_sample_partial_orthologs(partial_ortholo_dic):
    keys=partial_ortholo_dic.keys()
    distances=[]
    values=[]
    done=set([])
    while len(keys)>0:
        key=np.random.choice(keys)
        distances.append(partial_ortholo_dic[key][0])
        values.append(partial_ortholo_dic[key][1])
        for j in keys:
            if j!=key :
                if  j in partial_ortholo_dic.keys()and  key in partial_ortholo_dic.keys():
                    if  len(intersection(partial_ortholo_dic[key][2],partial_ortholo_dic[j][2]))>0:
                        del partial_ortholo_dic[j] # Delete all tratces that overlap with it
        del partial_ortholo_dic[key]
        keys=partial_ortholo_dic.keys()
    return distances, values


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
    #duplication_nodes = find_duplication_nodes(tree, '../Gene_duplications.txt') # georg's most recent version of this takes more parameters
    #ancestral_seqs, altall_seqs, pp = read_ancestral_seqs('../MLtree')
    species_dic = get_species(tips)
    matched_ancestors, matched_extant, matched_dic, age_dic = find_matched_ancestors(species_dic, tree,
                                                                                     duplication_nodes, nodes)
    if find_partials == True:  # Only do this if you have to, the calc is slow
        find_partial_orthologs(data_dic, matched_ancestors,
                               matched_extant)  # find partial orthologs and write partial ortholog file

    read_partial_ortholog_file()




if __name__ == "__main__":
    mytree, mytips, mynodes = read_tree('../190918_EA_tree1only.txt')
    mytip_dic = tip_dict(mytips, mynodes)
    # print(tips)
    # print(nodes)
    # someday we'll get there
    mydata_dic = load_experiment('../190919_medianEA66k-1.csv')
    myduplication_nodes = find_duplication_nodes(mytree, '../Gene_duplications.txt')
    print("dupe nodes", myduplication_nodes)
    species = get_species(mytips)
    print("species are: ", species)
    # orthos =find_orthologs(tree, data_dic, tip_dic, duplication_nodes)
    # nonovers = pick_nonoverlap(tree,orthos,tip_dic)
    # neutral_evo(mytree, mydata_dic, mytip_dic, myduplication_nodes)

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
