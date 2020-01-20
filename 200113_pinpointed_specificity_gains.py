import Bio.Phylo as Phylo
from Bio.Phylo import BaseTree
import numpy as np
import itertools
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import matplotlib.colors as cm
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


def draw_tree(input_tree, nodes, tips, coordinates, subfig, labtype):
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
                labs = clade_labels(j, labtype, False, input_tree)
                subfig.text(coordinate2[0] + 1.5, coordinate2[1] - 0.33, labs, fontsize=5)
            else:
                coordinate2 = coordinates[j]
                labs = clade_labels(j, labtype, True, input_tree)
                subfig.text(coordinate2[0] + .25, coordinate2[1], labs, fontsize=5)
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


def paralog_timing(paralog_anc_list, data_dic, tree):
    """Finds the interaction score for paralogs with homodimerizing ancestor and the distance between them
    Takes:
        paralog_anc_list: a list of [frozenset([protein paralogs]), ancestral protein]
        data_dic: a dictionary of keys = frozenset([protein pairs]) values = interaction score. Currently pairs are in the tree format
        tree: the gene tree
    Returns:
        intscorelist: a list that has built on paralog_anc_list to be [frozenset([protein paralogs]), ancestral protein, interaction score, distance between paralogs]"""
    intscorelist = []
    for i in paralog_anc_list:
        protein1, protein2 = list(i[0])[0], list(i[0])[1]
        try:
            if data_dic[frozenset([i[1]])] == 1:  # and data_dic[frozenset([protein1])] == 1 and data_dic[frozenset([protein2])] == 1:
                distance = Phylo.BaseTree.TreeMixin.distance(tree, protein1, protein2)
                if type(protein1) is str:
                    protein1 = protein1.lower()
                if type(protein2) is str:
                    protein2 = protein2.lower()
                intscore = data_dic[frozenset([protein1, protein2])]
                intscorelist.append([i, intscore, distance])
        except KeyError:
            print("Key ", i[1], " parent of ", list(i[0]), " was not in data dic.")
    print(len(intscorelist), " intscorelength and paralog_anc_listlength ", len(paralog_anc_list))
    return intscorelist


def info_on_specificity_gains(intscores):
    """Helper function to give statistics on what interactions gained specificity
    Takes:
        intscores: a retarded list of list of list situation. [[frozenset([protein1, protein2]), ancestor_clade], interaction_score, distance]"""
    numspec_gains = 0
    distances = []
    for i in intscores:
        if i[1] == 0.0:
            numspec_gains += 1
            distances.append(i[2])
    distances.sort()
    print("short distances ", distances[1:10])
    for j in distances[1:10]:
        for i in intscores:
            if i[2] == j:
                print("this is a nearby guy ", i)
    print("Of the ", len(distances), " distances, the shortest was ", min(distances), " the longest was ", max(distances), " and the average was ", np.mean(distances))


def plot_specificity(intscores, tree, coors, ax):
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
    arccol = {1.0: (255 / 256, 49 / 256, 49 / 256, 1), 0.5: (196 / 256, 0, 237 / 256, 1), 0.0: (0, 145 / 256, 244 / 256, 1)}
    for i in intscores:
        if len(i[0][0]) > 1:
            try:
                cox1, coy1 = coors[list(tree.find_clades(list(i[0][0])[0]))[0]]  # I'm really sorry but I need the clade object.
                cox2, coy2 = coors[list(tree.find_clades(list(i[0][0])[1]))[0]]
                xarc, yarc = parabola([cox1, coy1], [cox2, coy2], .25)
                ax.plot(xarc, yarc, color=arccol[i[1]], alpha=0.4, zorder=4)
            except KeyError:
                print("couldn't find", i[0][0], "in coordic")
        else:
            try:
                cox, coy = coors[list(i[0][0])[0]]
                ax.scatter(cox, coy, color=arccol[i[1]], s=12, zorder=10, edgecolor='black')
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
        coordinate1 = coors[list(tree.find_clades(list(i[0][0])[0]))[0]]
        ax.scatter(coordinate1[0], coordinate1[1] + random.randrange(0, 10)/10, color=coldots[c], s=28, zorder=10, edgecolors='black')
        passed_lca = False
        for j in Phylo.BaseTree.TreeMixin.trace(tree, list(tree.find_clades(list(i[0][0])[0]))[0], list(tree.find_clades(list(i[0][0])[1]))[0]):  # I'm so sorry, but I don't have a better way of getting clade objects
            if passed_lca:
                coordinate2 = coors[j]
                ax.plot([coordinate1[0], coordinate2[0]], [coordinate2[1], coordinate2[1]], alpha=0.4, color=(35/256, 138/256, 141/256, 1), zorder=1, linewidth=6)
                ax.plot([coordinate1[0], coordinate1[0]], [coordinate1[1], coordinate2[1]], alpha=0.4, color=(35/256, 138/256, 141/256, 1), zorder=1, linewidth=6)
                coordinate1 = coordinate2
            else:
                coordinate2 = coors[j]
                ax.plot([coordinate2[0], coordinate2[0]], [coordinate1[1], coordinate2[1]], alpha=0.4, color=(35/256, 138/256, 141/256, 1), zorder=1, linewidth=6)
                ax.plot([coordinate1[0], coordinate2[0]], [coordinate1[1], coordinate1[1]], alpha=0.4, color=(35/256, 138/256, 141/256, 1), zorder=1, linewidth=6)
                coordinate1 = coordinate2
            if Phylo.BaseTree.TreeMixin.common_ancestor(tree, list(tree.find_clades(list(i[0][0])[0]))[0], list(tree.find_clades(list(i[0][0])[1]))[0]) == j:
                passed_lca = True
        ax.scatter(coordinate1[0], coordinate1[1] + random.randrange(0, 10)/10, color=coldots[c], s=28, zorder=10, edgecolors='black')
    return


def all_paralog_plot_wrapper(tree, nodes, tips, coors, intscores, filename):
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
    fig = plt.figure(figsize=(40, 18))
    ax = fig.add_subplot(111)
    draw_tree(tree, nodes, tips, coors, ax, "gene")
    plot_specificity(intscores, tree, coors, ax)
    # plot_specificity_rects(tree, coors, ax, intscores)
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ratio = 2.5
    ax.set_aspect(abs((xright - xleft) / (ybottom - ytop)) * ratio)
    ax.set_axis_off()
    figure = plt.gcf()
    figure.savefig(filename, figsize=(100, 30))
    # plt.show()
    return


def basal_paralog_int_loss(intscores, tree):
    """Takes the intscores list and finds the two paralogs that are most closely related to the ancestor node for each andestral node and trims the intscores list to just those
    Takes:
        intscores: a dumb list of lists situation. [[frozenset([clade_object1, clade_object2]), ancestral_clade], interaction_score, distance]
        tree: the tree object that this occurred on
    Returns:
        trimmedscores: a dumb list of lists situation. [[frozenset([clade_object1, clade_object2]), ancestral_clade], interaction_score, distance]
    """
    trimmedscores = []
    ancestor_dic = {}
    for i in intscores:
        if i[1] == 0.0:
            if i[0][1] not in ancestor_dic:
                ancestor_dic[i[0][1]] = [i[0][0]]
            else:
                ancestor_dic[i[0][1]].append(i[0][0])
    for i in ancestor_dic:
        mindist = 100
        minpoint = "This should not occur"
        for j in ancestor_dic[i]:
            if Phylo.BaseTree.TreeMixin.distance(tree, list(j)[0], list(j)[1]) < mindist:
                minpoint = j
                mindist = Phylo.BaseTree.TreeMixin.distance(tree, list(j)[0], list(j)[1])
        trimmedscores.append([[minpoint, i], 0.0, mindist])
    print("Found ", len(trimmedscores), " paralog pairs that gain specificity")
    print(trimmedscores)
    return trimmedscores


def write_scores(scores, filename, input_tree):
    """writes scores and the proteins that caused them and their ancestor to a csv file
        Takes:
            scores: a dumb list of lists situation. [[frozenset([clade_object1, clade_object2]), ancestral_clade], interaction_score, distance]
             filename: the name of the csv file
             input_tree: the tree object this is coming from
        Returns:
            Nothing"""
    with open(filename, 'w') as f:
        f.write("X_peptide,Y_peptide,ancestor,interaction,distance\n")
        for i in scores:
            if type(list(i[0][0])[0]) is Phylo.Newick.Clade:
                xpep = Phylo.BaseTree.TreeMixin.get_nonterminals(input_tree).index(list(i[0][0])[0]) + 172
            else:
                xpep = str(list(i[0][0])[0]).lower()
            if type(list(i[0][0])[1]) is Phylo.Newick.Clade:
                ypep = Phylo.BaseTree.TreeMixin.get_nonterminals(input_tree).index(list(i[0][0])[1]) + 172
            else:
                ypep = str(list(i[0][0])[1]).lower()
            ancestor = Phylo.BaseTree.TreeMixin.get_nonterminals(input_tree).index(i[0][1]) + 172
            line = str(xpep) + ',' + str(ypep) + "," + str(ancestor) + "," + str(i[1]) + "," + str(i[2]) + '\n'
            f.write(line)
    f.close()
    return


def trace_specificity_gain(intscore, tree, nodes, tips, coors, data_dic, count):
    filename = '../200120_specificity_gains' + str(count) + '.pdf'
    protein1, protein2 = list(intscore[0][0])[0], list(intscore[0][0])[1]
    spectrace = []
    trace = Phylo.BaseTree.TreeMixin.trace(tree, protein1, protein2)
    trace.append(list(tree.find_clades(protein1))[0])
    for pair in list(itertools.combinations(trace, 2)):
        try:
            interact = data_dic[frozenset([pair[0], pair[1]])]
            spectrace.append([[frozenset([pair[0], pair[1]]), 0], interact, 0])
        except KeyError:
            pass
    for i in trace:
        try:
            interact = data_dic[frozenset([i])]
            spectrace.append([[frozenset([i]), 0], interact, 0])
        except KeyError:
            try:
                interact = data_dic[frozenset([str(i).lower()])]
                spectrace.append([[frozenset([i]), 0], interact, 0])
            except KeyError:
                pass
    all_paralog_plot_wrapper(tree, nodes, tips, coors, spectrace, filename)


if __name__ == "__main__":
    mydata_dic = load_experiment('../190919_medianEA66k-1.csv')
    mytree, mytips, mynodes = read_tree('../190918_EA_tree1only.txt')
    myduplication_nodes = find_duplication_nodes(mytree, '../Gene_duplications.txt')
    myspecies = get_species(mytips)
    myspecies_tree, species_tips, species_nodes = read_tree('../speciestree_node_names.newick.nwk')
    mymatched_ancestors, mymatched_dic = find_matched_ancestors(myspecies, myspecies_tree, mytree, myduplication_nodes)
    mypars = paralog_ancestors(mymatched_dic, mytree)
    mynewdict = convert_data_dic_to_species_names(mydata_dic, mynodes, mytips)
    myscores = paralog_timing(mypars, mynewdict, mytree)
    mycoors = get_tree_coordinates(mytree, mynodes, mytips)

    info_on_specificity_gains(myscores)
    mytrimmedscores = basal_paralog_int_loss(myscores, mytree)
    # write_scores(trimmedscores, "../200113_basal_paralogs_with_specgain.csv", mytree)
    # all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, trimmedscores, '../200114_basal_specificity_gains.pdf')

    for mycount, score in enumerate(mytrimmedscores):
        trace_specificity_gain(score, mytree, mynodes, mytips, mycoors, mynewdict, mycount)
