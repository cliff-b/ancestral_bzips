# -*- coding: utf-8 -*-
import Bio.Phylo as Phylo
from Bio.Phylo import BaseTree
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize


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
#    ASRtree,tips, nodes=read_tree(treename+'/tree1/tree2') #get tree from Lazarus with reoptimized BLs,
#    currently I have to manually make this file because lazarus formats it stupidly
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
#        print ('There are '+str(len(set(ancestral_seqs.values())))+' unique ML sequences and '+ str(len(altall_seqs))+'
#        altall sequences on '+ treename)
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


if __name__ == "__main__":
    mytree, mytips, mynodes = read_tree('../190918_EA_tree1only.txt')
    mytip_dic = tip_dict(mytips, mynodes)
    # print(tips)
    # print(nodes)
    # someday we'll get there
    mydata_dic = load_experiment('../190919_medianEA66k-1.csv')
    myduplication_nodes = find_duplication_nodes(mytree, '../Gene_duplications.txt')
    # orthos =find_orthologs(tree, data_dic, tip_dic, duplication_nodes)
    # nonovers = pick_nonoverlap(tree,orthos,tip_dic)
    neutral_evo(mytree, mydata_dic, mytip_dic, myduplication_nodes)

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
