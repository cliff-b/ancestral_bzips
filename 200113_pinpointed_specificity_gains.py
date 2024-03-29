import Bio.Phylo as Phylo
from Bio.Phylo import BaseTree
import numpy as np
from scipy import stats
import itertools
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import math
from scipy.interpolate import CubicSpline
from matplotlib import cm
import random


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
        filename: string, of a csv file with peptides in cols 1, 2 and interactions in 12
    returns:
         data_dic:   a dictionary of frozen sets of strings, key =frozenset({'pep1', 'pep2'}) value = interaction score
         """
    f = open(filename, 'r')
    lines = f.readlines()
    data_dic, raw_dic = {}, {}
    for i in lines[1:]:
        line = i.split()[0].split(',')
        if len(line[0].split('+')) > 1 or len(line[1].split('+')) > 1:
            for j in line[0].split('+'):
                for k in line[1].split('+'):
                    try:
                        data_dic[frozenset([j.lower(), k.lower()])] = np.average([data_dic[frozenset([j.lower(), k.lower()])], float(line[12])])
                        raw_dic[frozenset([j.lower(), k.lower()])] = np.average([raw_dic[frozenset([j.lower(), k.lower()])], float(line[4])])
                    except KeyError:
                        data_dic[frozenset([j.lower(), k.lower()])] = float(line[12])
                        raw_dic[frozenset([j.lower(), k.lower()])] = float(line[4])
        else:
            try:
                data_dic[frozenset([line[0].lower(), line[1].lower()])] = np.average([data_dic[frozenset([line[0].lower(), line[1].lower()])], float(line[12])])
                raw_dic[frozenset([line[0].lower(), line[1].lower()])] = np.average([raw_dic[frozenset([line[0].lower(), line[1].lower()])], float(line[4])])
            except KeyError:
                data_dic[frozenset([line[0].lower(), line[1].lower()])] = float(line[12])
                raw_dic[frozenset([line[0].lower(), line[1].lower()])] = float(line[4])
    return data_dic, raw_dic


def clade_object_to_readable(clade, tree):
    """Takes a clade object and turns it into a human readable number or leaves it as a string
    Takes:
        clade: a clade object
        tree: the tree it came from, presumably the gene tree we're using here.
    Returns:
        a string of the clade
        """
    if type(clade.name) is str:
        return clade.name
    else:
        return str(tree.get_nonterminals().index(clade) + 172)


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


def get_descendants(tree, node, includeself=False):
    """Gets all the descendants of a node and returns them as a list
    Takes:
        tree: a tree object that the node is from
        node: a clade object that is the node that we want the descendants of
        includeself: boolean, if true the node will be included in the list of descendants
    Returns:
        descendants: a list of clade objects that are the descentants of a node"""
    descendants = []
    nodetips = node.get_terminals()
    for i in nodetips:
        trace = Phylo.BaseTree.TreeMixin.trace(tree, node, i)
        for j in trace[1:]:
            if j not in descendants:
                descendants.append(j)
    if includeself:
        descendants.append(node)
    return descendants


def get_tree_coordinates(tree, nodes, tips, cladogram=False):
    """ Takes a tree and it's corresponding tips and nodes to get the cartesian coordinates of each clade.
    Builds a dictionary of clades to coordinates by enumerating the tips and then looping through the nodes and locating
    the nodes in relation to those tips.
    Takes:
        tree: a tree object
        nodes: a list of clade objects that are nodes on the tree
        tips: a list of clade objects that are tips on the tree
        cladogram: boolean for whether the tree's branchlengths should be used or if it should be converted to a cladogram
    Returns:
        corr: a dictionary of strings of clades to lists of numerics. key = str(clade object), value = [x-coord, y-coord]
    """
    corr = {}
    for c, i in enumerate(tips):
        corr[i] = [Phylo.BaseTree.TreeMixin.distance(tree, nodes[0], i), c * -1]
    for nod in Phylo.BaseTree.TreeMixin.get_nonterminals(tree, order='postorder'):
        if cladogram:
            xval = min(corr[nod[0]][0], corr[nod[1]][0]) - 1  # go one to the left of the current left most
        else:
            try:
                xval = min(corr[nod[0]][0], corr[nod[1]][0]) - min(nod[0].branch_length, nod[1].branch_length)
            except TypeError:
                xval = min(corr[nod[0]][0], corr[nod[1]][0])-1
        yval = np.mean([corr[nod[0]][1], corr[nod[1]][1]])  # average between child nodes
        corr[nod] = [xval, yval]
    return corr


def sort_list(list1, list2):
    """Helper function that takes two lists and combines them to be sorted by the first list"""
    zipped_pairs = zip(list2, list1)
    z = [x for _, x in sorted(zipped_pairs)]
    return z


def clade_labels(clade, labtype, node, input_tree):
    """
    Takes information on the tree and turns it into human readable information about what things are what on the tree
    Takes:
        clade: the clade object for the individual point on the tree
        labtype: a string either 'species' or gene to allow the formatting for either species or gene trees
        node: Boolean, is this a node or not.
        input_tree: only necessary for nodes on gene trees, as this gets lookedup for what number it is
    Returns:
          labs: a string representing what clade was at a certain point
    """
    cladetypes = ['VBP', 'HLF_insects', 'HLF_', 'HLF', 'TEF', 'E4BP4', 'DBP', 'PAR', 'Par1', 'Par2', 'New', 'Weird']  # this order matters because I am a bad programmer
    if labtype == "species" and not node:
        genus = str(clade)[0]
        species = str(clade).split("_")[1]
        labs = genus + ". " + species
    elif labtype == "gene" and not node:
        ctype = ""
        for i in cladetypes:
            if i in str(clade):
                ctype = i
        genus = str(clade)[len(ctype):][0]
        species = str(clade)[len(ctype):].split("_")[1]
        labs = ctype + ' ' + genus + '. ' + species
    elif labtype == "gene" and node:
        labs = Phylo.BaseTree.TreeMixin.get_nonterminals(input_tree).index(clade) + 172
    else:
        raise ValueError("Labtype was not gene or species or node was not false")
    return labs


def draw_tree(input_tree, nodes, tips, coordinates, subfig):
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
        coordinate1 = coordinates[input_tree.root]  # this should be the most ancestral node. Might cause failures with weird trees
        trace = Phylo.BaseTree.TreeMixin.trace(input_tree, nodes[0], i)
        # Follow the trace from the most ancestral node to tip i and draw each coordinate pair
        for j in trace:
            if j in tips:
                coordinate2 = coordinates[j]
            else:
                coordinate2 = coordinates[j]
            if frozenset([frozenset(coordinate1), frozenset(coordinate2)]) not in exclude:
                subfig.plot([coordinate1[0], coordinate1[0]], [coordinate1[1], coordinate2[1]], color=col, zorder=1)
                subfig.plot([coordinate1[0], coordinate2[0]], [coordinate2[1], coordinate2[1]], color=col, zorder=1)
            exclude.add(frozenset([frozenset(coordinate1), frozenset(coordinate2)]))
            coordinate1 = coordinate2
    return


def parabola(point1, point2, power):
    """"Makes a parabola on the x-axis from two points, and returns it. Importantly this fails if both points are at the same Y point
    Takes:
        point1: a pair of coordinates [x, y] where both are floats
        point2: a pair of coordinates [x, y] where both are floats
    Returns:
        y: a list of the x coordinates of the parabola
        x: a list of the y coordinates of the parabola
    """
    point3 = [max([point1[0], point2[0]]) + (abs((point1[1] - point2[1])/2)**power), np.average([point1[1], point2[1]])]
    x = np.linspace(min(point1[1], point2[1], point3[1]), max(point1[1], point2[1], point3[1]), 20)
    xs = sort_list([point1[0], point2[0], point3[0]], [point1[1], point2[1], point3[1]])
    ys = [point1[1], point2[1], point3[1]]
    ys.sort()
    cs = CubicSpline(ys, xs, axis=1)
    y = cs(x)
    return y, x


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
    matched_ancestors = []
    LCAs = set()
    total = 0

    # Build the speciation_dic by looping through combinations of proteins from different combinations of species and
    # finding their ancestor, then finding the corresponding ancestral species and add that protein to their ancestral list
    species_combs = list(itertools.combinations(species_dic, 2))
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
                LCAs.add(ancestor)

    # Populate matched_dic with extant protein pairs from the same species, then from pairs from species in speciation_dic
    for i in species_dic:
        if len(species_dic[i]) == 1:
            clade = list(tree.find_clades(species_dic[i][0]))[0]
            matched_dic[i] = [frozenset([clade])]
        else:
            clades = [list(tree.find_clades(x))[0] for x in species_dic[i]]  # Get the clade objects for use as our output
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


def paralog_ancestors(matched_dic, gene_tree):
    """Finds the ancestors of proteins that exist in the same organism
    Takes:
        matched_dic: a dictionary of key = species name, values = list of frozen sets of protein pairs that existed in that species
        gene_tree: phylo.basetree object of the gene relation tree
    Returns:
        paralog_anc_list: a list containing a frozenset of proteins and their LCA"""
    paralog_anc_dic = {}
    for i in matched_dic:
        for j in matched_dic[i]:
            if len(j) > 1:
                lca = Phylo.BaseTree.TreeMixin.common_ancestor(gene_tree, [list(j)[0], list(j)[1]])
                paralog_anc_dic[j] = {'ancestor': frozenset([lca])}
    return paralog_anc_dic


def paralog_timing(paralog_anc_list, data_dic, tree):
    """Finds the interaction score for paralogs with homodimerizing ancestor and the distance between them
    Takes:
        paralog_anc_list: a list of [frozenset([protein paralogs]), ancestral protein]
        data_dic: a dictionary of keys = frozenset([protein pairs]) values = interaction score. Currently pairs are in the tree format
        tree: the gene tree
    Returns:
        intscorelist: a list that has built on paralog_anc_list to be [frozenset([protein paralogs]), ancestral protein, interaction score, distance between paralogs]"""
    ignoretiming = False
    if not all([data_dic[k] in (1.0, 0.5, 0.0, -1.0) for k in data_dic]):
        ignoretiming = True

    for i in paralog_anc_list:
        protein1, protein2 = list(i)
        try:
            if data_dic[paralog_anc_list[i]['ancestor']] == 1 or ignoretiming:  # and data_dic[frozenset([protein1])] == 1 and data_dic[frozenset([protein2])] == 1:
                distance = Phylo.BaseTree.TreeMixin.distance(tree, protein1, protein2)
                paralog_anc_list[i]['int_score'] = data_dic[i]
                paralog_anc_list[i]['distance'] = distance
        except KeyError:
            print("Key ", clade_object_to_readable(list(paralog_anc_list[i]['ancestor'])[0], tree), " parent of ", list(i), " was not in data dic.")
    print("Number of things in paralog_anc_list ", len(paralog_anc_list))
    return paralog_anc_list


def plot_specificity(intscores, coors, ax):
    """
    Draws interactions as parabolas between points or as points on a tree to in the color of their interaction
    Takes:
        intscores: a dumb list of lists situation. [[frozenset([clade_object1, clade_object2]), ancestral_clade], interaction_score, distance]
        tree: a tree object this should be plotted on
        coors, a coordinate dictionary of clades to a list of values, coor[clade] = [x-coor, y-coor]
        ax: a plot to plot things on to.
    Returns:
        Nothing
    """
    if all([intscores[i]['int_score'] in (1.0, 0.5, 0.0, -1.0) for i in intscores if 'int_score' in intscores[i]]):
        arccol = {1.0: (255 / 256, 49 / 256, 49 / 256, 1), 0.5: (196 / 256, 0, 237 / 256, 1), 0.0: (0, 145 / 256, 244 / 256, 1), -1.0: (0.5, 0.5, 0.5, 1)}
    else:
        onlyscores = [intscores[i]['int_score'] for i in intscores if 'int_score' in intscores[i]]
        scoremax = max(onlyscores)
        scoremin = min(onlyscores)
        arccol = {-1.0: (0.5, 0.5, 0.5, 1)}
        for score in onlyscores:
            if score != -1:
                normscore = (score - scoremin) / (scoremax - scoremin)
                if normscore <= 0.5:
                    arccol[score] = (196/256 * normscore*2, 145/256 - 145/256 * normscore*2, 244/256 - 7/256*normscore*2, 1)
                elif 0.5 < normscore <= 1.0:
                    arccol[score] = (196/256 + 59/256*(normscore-0.5)*2, 0 + 49/256*(normscore-0.5)*2, 237/256 - 188/256*(normscore-0.5)*2, 1)
                else:
                    print('whoopsies')
    for i in intscores:
        if len(i) > 1:
            try:
                cox1, coy1 = coors[list(i)[0]]  # I'm really sorry but I need the clade object.
                cox2, coy2 = coors[list(i)[1]]
                xarc, yarc = parabola([cox1, coy1], [cox2, coy2], .25)
                ax.plot(xarc, yarc, color=arccol[intscores[i]['int_score']], alpha=0.4, zorder=4)
            except KeyError:
                print("couldn't find", list(i)[0], "in coordic")
        else:
            try:
                cox, coy = coors[list(i)[0]]
                ax.scatter(cox, coy, color=arccol[intscores[i]['int_score']], s=12, zorder=5, edgecolor='black')
            except KeyError:
                pass
    return


def plot_specificity_rects(tree, coors, ax, intscores):
    """
    A function to draw rectangles across a tree between two nodes that gain specificity on a tree
    Takes:
        intscores: a dumb list of lists situation. [[frozenset([clade_object1, clade_object2]), ancestral_clade], interaction_score, distance]
        tree: a tree object this should be plotted on
        coors, a coordinate dictionary of clades to a list of values, coor[clade] = [x-coor, y-coor]
        ax: a plot to plot things on to.
    Returns:
        Nothing
    """

    coldots = [cm.jet(x / len(intscores)) for x in range(len(intscores))]
    random.shuffle(coldots)
    for c, i in enumerate(intscores):
        coordinate1 = coors[list(i)[0]]
        ax.scatter(coordinate1[0], coordinate1[1] + random.randrange(0, 10)/10, color=coldots[c], s=28, zorder=10, edgecolors='black')
        passed_lca = False
        for j in Phylo.BaseTree.TreeMixin.trace(tree, list(i)[0], list(i)[1]):
            if passed_lca:
                coordinate2 = coors[j]
                ax.plot([coordinate1[0], coordinate2[0]], [coordinate2[1], coordinate2[1]], alpha=0.4, color=coldots[c], zorder=1, linewidth=6)
                ax.plot([coordinate1[0], coordinate1[0]], [coordinate1[1], coordinate2[1]], alpha=0.4, color=coldots[c], zorder=1, linewidth=6)
                coordinate1 = coordinate2
            else:
                coordinate2 = coors[j]
                ax.plot([coordinate2[0], coordinate2[0]], [coordinate1[1], coordinate2[1]], alpha=0.4, color=coldots[c], zorder=1, linewidth=6)
                ax.plot([coordinate1[0], coordinate2[0]], [coordinate1[1], coordinate1[1]], alpha=0.4, color=coldots[c], zorder=1, linewidth=6)
                coordinate1 = coordinate2
            if Phylo.BaseTree.TreeMixin.common_ancestor(tree, list(i)[0], list(i)[1]) == j:
                passed_lca = True
        ax.scatter(coordinate1[0], coordinate1[1], color=coldots[c], s=28, zorder=10, edgecolors='black')
    return


def all_paralog_plot_wrapper(tree, nodes, tips, coors, intscores, filename, descendscores=[], data_dic=[], graphtype=''):
    """Wrapper function to plot an interaction tree
    Takes:
        tree: a tree object this should be plotted on
        nodes: a list of clade objects that are nodes from the tree
        tips: a list of clade object that are tips from the tree
        coors, a coordinate dictionary of clades to a list of values, coor[clade] = [x-coor, y-coor]
        intscores: a dumb list of lists situation. [[frozenset([clade_object1, clade_object2]), ancestral_clade], interaction_score, distance]
        filename: The name to save the file to. Needs to have an extention figure.savefig can use.
    Returns:
        Nothing
    """
    fig = plt.figure(figsize=(18, 18))
    ax = fig.add_subplot(111)
    draw_tree(tree, nodes, tips, coors, ax)
    if graphtype == 'rects':
        plot_specificity_rects(tree, coors, ax, intscores)
    elif graphtype == 'circles':
        plot_species_circles(coors, intscores, ax, tree, descendscores, data_dic)  # this is a bad hack, descendscores is actually the tips
    else:
        plot_specificity(intscores, coors, ax)
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ratio = 2.5
    ax.set_aspect(abs((xright - xleft) / (ybottom - ytop)) * ratio)
    ax.set_axis_off()
    figure = plt.gcf()
    figure.savefig(filename, figsize=(100, 30))
    plt.close()
    # plt.show()
    return


def basal_paralog_int_loss(intscores, tree, score_for_loss=0):
    """Takes the intscores list and finds the two paralogs that are most closely related to the ancestor node for each ancestral node and trims the intscores list to just those
    Takes:
        intscores: a dumb list of lists situation. [[frozenset([clade_object1, clade_object2]), ancestral_clade], interaction_score, distance]
        tree: the tree object that this occurred on
    Returns:
        trimmedscores: a dumb list of lists situation. [[frozenset([clade_object1, clade_object2]), ancestral_clade], interaction_score, distance]
    """
    trimmedscores = {}
    ancestor_dic = {}
    for i in intscores:
        if 'int_score' in intscores[i]:
            if intscores[i]['int_score'] == score_for_loss:
                if intscores[i]['ancestor'] not in ancestor_dic:
                    ancestor_dic[intscores[i]['ancestor']] = [i]
                else:
                    ancestor_dic[intscores[i]['ancestor']].append(i)
    for i in ancestor_dic:
        mindist = 100
        minpoint = "This should not occur"
        for j in ancestor_dic[i]:
            if Phylo.BaseTree.TreeMixin.distance(tree, list(j)[0], list(j)[1]) < mindist:
                minpoint = j
                mindist = Phylo.BaseTree.TreeMixin.distance(tree, list(j)[0], list(j)[1])
        trimmedscores[minpoint] = intscores[minpoint]
        trimmedscores[minpoint]['mindist'] = mindist
    print("Found ", len(trimmedscores), " paralog pairs that gain specificity")
    return trimmedscores


def write_odds_of_regain(tree, trim_score, matched_dic, data_dict):
    score_dic = {1.0: "Y", 0.5: "W", 0.0: "N"}
    alls = []
    matchedvals = [item for sublist in list(matched_dic.values()) for item in sublist]
    for i in trim_score:
        n1 = list(i[0][0])[0]
        n2 = list(i[0][0])[1]
        n1desc = get_descendants(tree, n1)
        n2desc = get_descendants(tree, n2)
        perms = [frozenset([x, y]) for x in n1desc for y in n2desc]
        for j in perms:
            try:
                score = data_dict[j]
                all_dist = Phylo.BaseTree.TreeMixin.distance(tree, list(j)[0], list(j)[1])
                alls.append([clade_object_to_readable(i[0][1], tree), clade_object_to_readable(list(j)[0], tree), clade_object_to_readable(list(j)[1], tree), score, score_dic[score], all_dist, "all_desc"])
                if j in matchedvals:
                    alls.append([clade_object_to_readable(i[0][1], tree), clade_object_to_readable(list(j)[0], tree), clade_object_to_readable(list(j)[1], tree), score, score_dic[score], all_dist, "match_desc"])
            except KeyError:
                pass
        trace = Phylo.BaseTree.TreeMixin.trace(tree, n1, n2)
        trace.insert(0, n1)
        tracecombsn1 = [frozenset([x, y]) for x in trace for y in n1desc]
        for j in tracecombsn1:
            try:
                score = data_dict[j]
                all_dist = Phylo.BaseTree.TreeMixin.distance(tree, list(j)[0], list(j)[1])
                alls.append([clade_object_to_readable(i[0][1], tree), clade_object_to_readable(list(j)[0], tree), clade_object_to_readable(list(j)[1], tree), score, score_dic[score], all_dist, "trace_w_n1desc"])
            except KeyError:
                pass
        tracecombsn2 = [frozenset([x, y]) for x in trace for y in n2desc]
        for j in tracecombsn2:
            try:
                score = data_dict[j]
                all_dist = Phylo.BaseTree.TreeMixin.distance(tree, list(j)[0], list(j)[1])
                alls.append([clade_object_to_readable(i[0][1], tree), clade_object_to_readable(list(j)[0], tree), clade_object_to_readable(list(j)[1], tree), score, score_dic[score], all_dist, "trace_w_n2desc"])
            except KeyError:
                pass
    with open('../200129_Regains.csv', 'w') as f:
        f.write("Dupe_node,X_peptide,Y_peptide,intscore,interaction,distance,int_subdivision\n")
        for line in alls:
            f.write(str(line).replace("'", "")[1:-1] + '\n')
    f.close()
    return


def write_csv(scores, keylist, filename, input_tree):
    """writes scores and the proteins that caused them and their ancestor to a csv file
        Takes:
            scores: a dumb list of lists situation. [[frozenset([clade_object1, clade_object2]), ancestral_clade], interaction_score, distance]
             filename: the name of the csv file
             input_tree: the tree object this is coming from
        Returns:
            Nothing"""
    with open(filename, 'w') as f:
        f.write('xpep,ypep,')
        f.write(str(keylist))
        f.write('\n')
        for i in scores:
            if len(list(i)) == 2:
                line = clade_object_to_readable(list(i)[0], input_tree) + ',' + clade_object_to_readable(list(i)[1], input_tree) + ','
            elif len(list(i)) == 1:
                line = clade_object_to_readable(list(i)[0], input_tree) + ',' + clade_object_to_readable(list(i)[0],  input_tree) + ','
            for key in keylist:
                if type(scores[i][key]) is Phylo.Newick.Clade:
                    line = line + clade_object_to_readable(scores[i][key], input_tree) + ','
                elif type(scores[i][key]) is frozenset:
                    for j in list(scores[i][key]):
                        line = line + clade_object_to_readable(j, input_tree) + ','
                elif type(scores[i][key]) is tuple:
                    for tup in scores[i][key]:
                        line = line + str(tup) + '+'
                    line = line[:-1] + ','

                else:
                    line = line + str(scores[i][key]) + ','
            line = line[:-1] + '\n'
            f.write(line)
    f.close()
    return


def trace_specificity_gain(intscore, tree, data_dic):
    protein1, protein2 = list(intscore)
    spectrace, tracedescend = {}, {}
    trace = Phylo.BaseTree.TreeMixin.trace(tree, protein1, protein2)
    trace.insert(0, protein1)
    for pair in list(itertools.combinations(trace, 2)):
        try:
            interact = data_dic[frozenset([pair[0], pair[1]])]
            spectrace[frozenset([pair[0], pair[1]])] = {'int_score': interact}
        except KeyError:
            pass
    for i in trace:
        try:
            interact = data_dic[frozenset([i])]
            spectrace[frozenset([i])] = {'int_score': interact}
        except KeyError:
            try:
                interact = data_dic[frozenset([str(i).lower()])]
                print('odd')
                spectrace[frozenset([i])] = {'int_score': interact}
            except KeyError:
                spectrace[frozenset([i])] = {'int_score': -1}
    return spectrace


def plot_species_circles(species_coordinates, matched_dic, ax, species_tree,  tips, data_dic):
    """produces a species tree with proteins cluster around the tips, with heterodimers as lines between homodimeric points
    includes an ordering so that the prefixes that appear in the clades function are the order around the circle that
    proteins appear, and there's an outside circle that colors them"""
    clades = ['VBP', 'HLF', 'TEF',  'DBP', 'Par', 'New', 'Weird', 'E4BP4', '']
    clades = [c.lower() for c in clades]
    arccol = {1.0: (255 / 256, 49 / 256, 49 / 256, 1), 0.5: (196 / 256, 0, 237 / 256, 1), 0.0: (0, 145 / 256, 244 / 256, 1), -1.0: (0.5, 0.5, 0.5, 1)}
    scale, stretch = 0.022, 8
    tab10 = cm.get_cmap('tab10')
    ext_x, ext_y = [max(list(species_coordinates.values())[0]) + 1, abs(min(list(species_coordinates.values())[1])) + 1]
    for specie in matched_dic:
        proteins = list(set(x for l in matched_dic[specie] for x in l))
        protein2 = []
        for p in proteins:
            for cnt, c in enumerate(clades):
                if str(p).startswith(c):
                    protein2.append((p, cnt))
                    break
        proteins = sorted(protein2, key=lambda x: x[1])
        justprots = [prot[0] for prot in proteins]
        if proteins[0][0] in tips:
            for clade in Phylo.BaseTree.TreeMixin.find_clades(species_tree, name=specie):
                x, y = species_coordinates[clade]
            numspecs = len(proteins)
            angle = 360.0 / numspecs
            xs = [np.cos(math.radians(i * angle)) * scale * ext_x + x for i in range(numspecs)]
            ys = [np.sin(math.radians(i * angle)) * scale * ext_y * stretch + y for i in range(numspecs)]
            for cnt, gene in enumerate(proteins):
                try:
                    ax.scatter(xs[cnt], ys[cnt], color=arccol[data_dic[frozenset([gene[0]])]], s=20, zorder=10, edgecolors='black')
                except KeyError:
                    ax.scatter(xs[cnt], ys[cnt], color=arccol[-1.0], s=20, zorder=10, edgecolors='black')
                except IndexError:
                    print('shit bad')
                wed = Wedge((x, y), .4, (cnt * angle - angle/2), ((cnt + 1) * angle - angle/2), width=0.1, color=tab10.colors[gene[1]])
                ax.add_patch(wed)
            for pair in matched_dic[specie]:
                if len(pair) == 2:
                    try:
                        ax.plot([xs[justprots.index(list(pair)[0])], xs[justprots.index(list(pair)[1])]], [ys[justprots.index(list(pair)[0])], ys[justprots.index(list(pair)[1])]], color=arccol[data_dic[pair]], zorder=3)
                    except KeyError:
                        pass
    for cnt, c in enumerate(clades):
        wed = Wedge((10, -cnt * 1.3), .4, 60, 120, width=0.1, color=tab10.colors[cnt])
        ax.add_patch(wed)
        ax.text(10.4, -cnt*1.3, c)

    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ratio = 2.5
    ax.set_aspect(abs((xright - xleft) / (ybottom - ytop)) * ratio)
    return


if __name__ == "__main__":
    #
    # Processing of data into usable paralogs and trees
    #
    mydata_dic, myraws = load_experiment('../200120_most_aliases.csv')
    mytree, mytips, mynodes = read_tree('../190918_EA_tree1only.txt')
    myduplication_nodes = find_duplication_nodes(mytree, '../Gene_duplications.txt')
    myspecies = get_species(mytips)
    myspecies_tree, species_tips, species_nodes = read_tree('../speciestree_node_names.newick.nwk')
    mymatched_ancestors, mymatched_dic = find_matched_ancestors(myspecies, myspecies_tree, mytree, myduplication_nodes)
    mypars = paralog_ancestors(mymatched_dic, mytree)
    mynewdict = convert_data_dic_to_species_names(mydata_dic, mynodes, mytips)
    mynewraw = convert_data_dic_to_species_names(myraws, mynodes, mytips)
    myscores = paralog_timing(mypars, mynewdict, mytree)
    mycoors = get_tree_coordinates(mytree, mynodes, mytips)
    myspeccoors = get_tree_coordinates(myspecies_tree, species_nodes, species_tips)
    mytrimmedscores = basal_paralog_int_loss(myscores, mytree)

    figstomake = ['0A']

    if '0A' in figstomake:
        all_paralog_plot_wrapper(myspecies_tree, species_nodes, species_tips, myspeccoors, mymatched_dic, '../figures/201106_extant_species_ints.pdf', mytips, mynewdict, 'circles')

    #
    # Figure 2A ....... currently
    # All paralogs plotted in spagettigram fashion
    #
    if '2A' in figstomake:
        all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, myscores, '../figures/201006_all_paralogs_spaghetti.pdf', [], [], 'spaget')

    #
    # Figure 2B ...... currently
    # rectangle paths showing the 8 traces that specificity evolved over between paralogs
    #
    if '2B' in figstomake:
        all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, mytrimmedscores, '../figures/200505_traced_specificity_gains.pdf', [], [], 'rects')
    #
    # Figure S1 A-H ....... currently
    # spaghettigrams for each path that specificity was gained over for paralogs with all interactions on that path
    #
    if 'S1' in figstomake:
        for mycount, score in enumerate(mytrimmedscores):
            specs = trace_specificity_gain(score, mytree, mynewdict)
            all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, specs, str('../figures/200506_traced_specgain' + clade_object_to_readable(list(score)[0], mytree) + clade_object_to_readable(list(score)[1], mytree) + '.pdf'), [], [], 'spaget')
    #
    # Figure S2 A-H ........ currently
    # spagettigrams for each path specificty was gained over for paralogs with only paralogous interactions
    #
    if 'S2' in figstomake:
        for mycount, score in enumerate(mytrimmedscores):
            specs = trace_specificity_gain(score, mytree, mynewdict)
            parspecs = {}
            for key in specs:
                if len(key) == 1 or key in mypars:
                    parspecs[key] = specs[key]
            all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, parspecs, str('../figures/200506_spaghetti_paralogs' + clade_object_to_readable(list(score)[0], mytree) + clade_object_to_readable(list(score)[1], mytree) + '.pdf'), [], [], 'spaget')
    #
    # Figure 2C distance from root that specificity occurred
    #
    #
    if 'S3' in figstomake:
        myparscopy = paralog_ancestors(mymatched_dic, mytree)
        myparsraws = paralog_timing(myparscopy, mynewraw, mytree)
        all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, myparsraws, '../figures/201006_all_paralogs_rawscore_spaghetti.pdf', [], [], 'spaget')

    ############ compute number of paralogs that occurred before specificity was gained
    for key in mytrimmedscores.keys():
        protein1, protein2 = list(key)
        mytrimmedscores[key]['proteinpair'] = (protein1, protein2)
        numpars = 0
        pars = []
        trace = Phylo.BaseTree.TreeMixin.trace(mytree, protein1, protein2)
        trace.insert(0, protein1)
        for pair in list(itertools.combinations(trace, 2)):
            if frozenset({pair[0], pair[1]}) in mypars.keys():
                try:
                    interact = mynewdict[frozenset([pair[0], pair[1]])]
                    numpars += 1
                    pars.append((frozenset({pair}), interact))
                    mytrimmedscores[key]['numpars'] = numpars
                    mytrimmedscores[key]['parnames'] = pars
                except KeyError:
                    numpars += 1
                    pars.append((frozenset({pair[0], pair[1]}), -1))
                    mytrimmedscores[key]['numpars'] = numpars
                    mytrimmedscores[key]['parnames'] = pars
        # add computation to figure out which side changes occured on
        interact1, interact2 = -1, -1
        try:
            interact1 = mynewdict[mytrimmedscores[key]['ancestor'].union(frozenset({protein1}))]  # it's dumb but the ancestor is already a frozenset
        except KeyError:
            pass
        try:
            interact2 = mynewdict[mytrimmedscores[key]['ancestor'].union(frozenset({protein2}))]
        except KeyError:
            pass
        mytrimmedscores[key]['anc_to_diverge'] = (interact1, interact2)
    #
    #   Figure 2E
    #
    write_csv(mytrimmedscores, ('ancestor','distance','numpars'), '../data_for_R/200929_num_paralogs_till_spec.csv', mytree)
    write_csv(mytrimmedscores, ('distance','anc_to_diverge'), '../data_for_R/200722_anc_to_diverge_ints.csv', mytree)

    ######## calculate number of interacting paralogs
    for key in mypars.keys():
        protein1, protein2 = list(key)
        mypars[key]['proteinpair'] = (protein1, protein2)
        try:
            mypars[key]['int_score'] = mynewdict[key]
        except KeyError:
            mypars[key]['int_score'] = -1
        try:
            mypars[key]['distance'] = Phylo.BaseTree.TreeMixin.distance(mytree, list(key)[0], list(key)[1])
        except:
            myscores[key]['distance'] = -1
    write_csv(myscores, ('int_score','distance'), '../data_for_R/200722_allpar_score_and_dist.csv', mytree)

    homodimers = {}
    for node in mynodes:
        homodimers[frozenset([node])] = {}
        try:
            homodimers[frozenset([node])]['intscore'] = mynewdict[frozenset([node])]
            homodimers[frozenset([node])]['extant'] = 'No'
        except KeyError:
            homodimers[frozenset([node])]['intscore'] = -1
            homodimers[frozenset([node])]['extant'] = 'No'
    for tip in mytips:
        homodimers[frozenset([tip])] = {}
        try:
            homodimers[frozenset([tip])]['intscore'] = mynewdict[frozenset([tip])]
            homodimers[frozenset([tip])]['extant'] = 'Yes'
        except KeyError:
            homodimers[frozenset([tip])]['intscore'] = -1
            homodimers[frozenset([tip])]['extant'] = 'Yes'

    write_csv(homodimers, ('intscore','extant'), '../data_for_R/200930_homodimer_counting.csv', mytree)

    if '3D' in figstomake:
        alltreebits =  mynodes + mytips
        distcors = {}
        for k1 in alltreebits:
            for k2 in alltreebits:
                distcors[frozenset([k1, k2])] = {}
                distcors[frozenset([k1, k2])]['pearson'] = ''
                distcors[frozenset([k1, k2])]['para'] = 'ortholog'
                distcors[frozenset([k1, k2])]['spearman'] =''
                distcors[frozenset([k1, k2])]['pearsonlog'] = ''
                distcors[frozenset([k1, k2])]['parentscore'] = ''
                distcors[frozenset([k1, k2])]['paraparents'] = ''
                distcors[frozenset([k1, k2])]['ancestor'] = Phylo.BaseTree.TreeMixin.common_ancestor(mytree, [k1, k2])
                distcors[frozenset([k1, k2])]['isparentchild'] = 'No'
                distcors[frozenset([k1, k2])]['k1mean'] = ''
                distcors[frozenset([k1, k2])]['k2mean'] = ''
                if frozenset([k1, k2]) in mypars:
                    distcors[frozenset([k1, k2])]['para'] = 'paralogs'
                    try:
                        partrace = Phylo.BaseTree.TreeMixin.trace(mytree, k1, k2)
                        if len(partrace) > 1:
                            distcors[frozenset([k1, k2])]['parentscore'] = mynewdict[frozenset([partrace[0], partrace[-2]])]
                            distcors[frozenset([k1, k2])]['paraparents'] = (clade_object_to_readable(partrace[0], mytree), clade_object_to_readable(partrace[-2], mytree))
                        else:
                            print('this trace was missing ', [clade_object_to_readable(x, mytree) for x in partrace])
                    except:
                        print('whooopsies')
                elif k1 == k2:
                    distcors[frozenset([k1, k2])]['para'] = 'homodimer'
                try:
                    parenttrace = Phylo.BaseTree.TreeMixin.trace(mytree, k1, k2)
                    if len(parenttrace) == 1:
                        distcors[frozenset([k1, k2])]['isparentchild'] = 'Yes'
                except:
                    pass
                try:
                    distcors[frozenset([k1, k2])]['distance'] = Phylo.BaseTree.TreeMixin.distance(mytree, k1, k2)
                    distcors[frozenset([k1, k2])]['raw'] = mynewraw[frozenset([k1, k2])]
                    distcors[frozenset([k1, k2])]['score'] = mynewdict[frozenset([k1, k2])]
                except:
                    distcors[frozenset([k1, k2])]['raw'] = -1
                    distcors[frozenset([k1, k2])]['score'] = -1
                if k1 != k2:
                    k1ints, k2ints = [], []
                    for k3 in alltreebits:
                        try:
                            k1val = mynewraw[frozenset([k1, k3])]
                            k2val = mynewraw[frozenset([k2, k3])]
                            k1ints.append(k1val)
                            k2ints.append(k2val)
                        except KeyError:
                            pass
                    if k1ints:
                        cor = np.array([k1ints, k2ints])
                        distcors[frozenset([k1, k2])]['pearson'] = np.corrcoef(cor)[1][0]
                        logs = np.log(cor.clip(min=0.001))
                        distcors[frozenset([k1, k2])]['pearsonlog'] = np.corrcoef(logs)[1][0]
                        rho, p = stats.spearmanr(cor, axis = 1)
                        distcors[frozenset([k1, k2])]['spearman'] = rho
                        distcors[frozenset([k1, k2])]['k1mean'], distcors[frozenset([k1, k2])]['k2mean'] = np.average(cor,1)


        write_csv(distcors, ('raw', 'distance', 'pearson', 'para', 'pearsonlog','spearman','score','parentscore', 'paraparents', 'ancestor','isparentchild'), '../data_for_R/200930_interactionprofiles_cors.csv', mytree)

    badbranches = {}
    for i in distcors:
        if distcors[i]['isparentchild'] == 'Yes' and distcors[i]['pearsonlog'] != '' and float(distcors[i]['pearsonlog']) < 0.25:
           badbranches[i] = 'Yes'
    all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, badbranches, '../figures/201102_branches_where_things_change.pdf', [], [],'rects')

    # tracedcoeffs = {}
    # for i in mytrimmedscores:
    #     mytrace = Phylo.BaseTree.TreeMixin.trace(mytree, list(i)[.,0], list(i)[1])

    sharedbranches = {}
    for i in badbranches:
        if i in mytrimmedscores:
            sharedbranches[i] = 'Yes'
    all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, sharedbranches, '../figures/201103_branches_where_things_change_and_gain_spec.pdf', [], [], 'rects')


    parentchild = {}
    for i in distcors:
        if distcors[i]['isparentchild'] == 'Yes':
            parentchild[i] = {}
            parentchild[i]['intscore'] = distcors[i]['score']
    all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, parentchild,
                             '../figures/201103_parent_child_spaghetti.pdf', [],
                             [], 'spaget')


