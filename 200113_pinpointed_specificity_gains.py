import Bio.Phylo as Phylo
from Bio.Phylo import BaseTree
import numpy as np
import itertools
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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


def specificity_after_duplication(tree, data_dic, duplication_nodes):
    """"""
    duplengths = {}
    for i in duplication_nodes:
        try:
            mindist = 100
            minpath = ''
            if data_dic[frozenset([i])] == 1.0:
                descendants = get_descendants(tree, i)
                for j in list(itertools.combinations(descendants, 2)):
                    if not any(x in Phylo.BaseTree.TreeMixin.trace(tree, i, j[1])[1:] for x in Phylo.BaseTree.TreeMixin.trace(tree, i, j[0])[1:]):  # make sure they're not part of the same lineage
                        try:
                            if data_dic[frozenset([j[0], j[1]])] == 0.0:
                                curdist = Phylo.BaseTree.TreeMixin.distance(tree, j[0], j[1])
                                if curdist < mindist:
                                    mindist = curdist
                                    minpath = Phylo.BaseTree.TreeMixin.trace(tree, j[0], j[1])
                                    minpath.insert(0, j[0])
                        except KeyError:
                            pass
                            #  print("bad pair", clade_object_to_readable(j[0], tree),"  ", clade_object_to_readable(j[1], tree))
            else:
                print(clade_object_to_readable(i, tree), " was not found a homodimer")
            duplengths[i] = [mindist, minpath]
        except KeyError:
            print(i, " was not in data_dic")
    print("pause momento")
    intscores = []
    # [[frozenset([clade_object1, clade_object2]), ancestral_clade], interaction_score, distance]
    for i in duplengths:
        if duplengths[i][0] != 100:
            intscores.append([[frozenset([duplengths[i][1][0], duplengths[i][1][-1]]), i], 0.0, duplengths[i][0]])
    return duplengths, intscores


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
                # labs = clade_labels(j, labtype, False, input_tree)
                # subfig.text(coordinate2[0] + 1.5, coordinate2[1] - 0.33, labs, fontsize=5)
            else:
                coordinate2 = coordinates[j]
                # labs = clade_labels(j, labtype, True, input_tree)
                # subfig.text(coordinate2[0] + .25, coordinate2[1], labs, fontsize=5)
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
    for i in paralog_anc_list:
        protein1, protein2 = list(i)
        try:
            if data_dic[paralog_anc_list[i]['ancestor']] == 1:  # and data_dic[frozenset([protein1])] == 1 and data_dic[frozenset([protein2])] == 1:
                distance = Phylo.BaseTree.TreeMixin.distance(tree, protein1, protein2)
                paralog_anc_list[i]['int_score'] = data_dic[i]
                paralog_anc_list[i]['distance'] = distance
        except KeyError:
            print("Key ", clade_object_to_readable(list(paralog_anc_list[i]['ancestor'])[0], tree), " parent of ", list(i), " was not in data dic.")
    print("Number of things in paralog_anc_list ", len(paralog_anc_list))
    return paralog_anc_list


def just_get_homos(data_dic, tree):
    intscorelist = []
    for i in tree.get_terminals():  # this is really dumb but there's not an obvious way to interate through the entire tree.
        try:
            intscorelist.append([[frozenset([i]), 0], data_dic[frozenset([i])], 0])
        except KeyError:
            intscorelist.append([[frozenset([i]), 0], -1, 0])
    for i in tree.get_nonterminals():
        try:
            intscorelist.append([[frozenset([i]), 0], data_dic[frozenset([i])], 0])
        except KeyError:
            intscorelist.append([[frozenset([i]), 0], -1, 0])
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


def plot_specificity(intscores, tree, coors, ax, descendscores, data_dic):
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
    arccol = {1.0: (255 / 256, 49 / 256, 49 / 256, 1), 0.5: (196 / 256, 0, 237 / 256, 1), 0.0: (0, 145 / 256, 244 / 256, 1), -1.0: (0.5, 0.5, 0.5, 1)}
    for i in intscores:
        if len(i) > 1:
            try:
                cox1, coy1 = coors[list(i)[0]]  # I'm really sorry but I need the clade object.
                cox2, coy2 = coors[list(i)[1]]
                xarc, yarc = parabola([cox1, coy1], [cox2, coy2], .25)
                ax.plot(xarc, yarc, color=arccol[intscores[i]['int_score']], alpha=0.4, zorder=4)
            except:
                print("couldn't find", i[0][0], "in coordic")
        else:
            try:
                cox, coy = coors[list(i)[0]]
                ax.scatter(cox, coy, color=arccol[intscores[i]['int_score']], s=12, zorder=5, edgecolor='black')
            except:
                pass
            # try:
            #     cox, coy = coors[list(i[0][0])[1]]
            #     ax.scatter(cox, coy, color=arccol[data_dic[frozenset([list(i[0][0])[1]])]], s=12, zorder=3, edgecolor='black')
            # except KeyError:
            #     pass
    # for i in intscores:
    #     if len(i[0][0]) == 1:
    #         circ1 =  patches.Ellipse(coors[list(i[0][0])[0]], width = 0.25, height= 3, fill = False, edgecolor = 'Blue', linewidth = 2, zorder = 1)
    #         ax.add_patch(circ1)
    #         break
    # circ2 =  patches.Ellipse(coors[list(intscores[-1][0][0])[0]], width = 0.25, height= 3, fill = False, edgecolor = 'Blue', linewidth = 2, zorder = 1)
    # ax.add_patch(circ2)
    # for i in descendscores:
    #     if len(i[0][0]) > 1:
    #         try:
    #             cox1, coy1 = coors[list(tree.find_clades(list(i[0][0])[0]))[0]]  # I'm really sorry but I need the clade object.
    #             cox2, coy2 = coors[list(tree.find_clades(list(i[0][0])[1]))[0]]
    #             xarc, yarc = parabola([cox1, coy1], [cox2, coy2], .25)
    #             ax.plot(xarc, yarc, color=arccol[i[1]], alpha=0.4, zorder=4)
    #         except KeyError:
    #             print("couldn't find", i[0][0], "in coordic")
    #     else:
    #         try:
    #             cox, coy = coors[list(i[0][0])[0]]
    #             ax.scatter(cox, coy, color=arccol[i[1]], s=12, zorder=10, edgecolor='black')
    #         except KeyError:
    #             pass
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
    linecol = (35/256, 138/256, 141/256, 1)
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


def all_paralog_plot_wrapper(tree, nodes, tips, coors, intscores, filename, descendscores=[], data_dic=[], rects=True):
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
    draw_tree(tree, nodes, tips, coors, ax, "gene")
    #plot_specificity(intscores, tree, coors, ax, descendscores, data_dic)
    if rects:
        plot_specificity_rects(tree, coors, ax, intscores)
    else:
        plot_specificity(intscores, tree, coors, ax, descendscores, data_dic)
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ratio = 2.5
    #ax.set_aspect(abs((xright - xleft) / (ybottom - ytop)) * ratio)
    ax.set_axis_off()
    figure = plt.gcf()
    figure.savefig(filename, figsize=(100, 30))
    plt.close()
    # plt.show()
    return

def dupes_intscores(node, tree, data_dic):
    nlist = get_descendants(tree, node, True)
    nodescore = []
    for i in list(itertools.combinations(nlist, 2)):
        try:
            score =  data_dic[frozenset([i[0], i[1]])]
            dist = Phylo.BaseTree.TreeMixin.distance(tree, i[0], i[1])
            intsc = [[frozenset([i[0], i[1]]), node], score, dist]
            nodescore.append(intsc)
        except KeyError:
            pass
    for i in nlist:
        try:
            score =  data_dic[frozenset([i])]
            intsc = [[frozenset([i]), node], score, 0.0]
            nodescore.append(intsc)
        except KeyError:
            pass
    return nodescore

def calc_avg_dist(clade, tree):
    dists = []
    desc = Phylo.BaseTree.TreeMixin.get_terminals(clade)
    for term in desc:
        d = Phylo.BaseTree.TreeMixin.distance(tree, clade, term)
        dists.append(d)
    return np.mean(dists)





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

def write_possible_paralogs(tree, dupe_nodes, trim_score, matched_dic):
    with open('../200129.2_specgains.csv', 'w') as f:
        f.write("Dupe_node,distance,num_pars,specific,time,parlist\n")
        for i in dupe_nodes:
            dists = []
            myspec = "No"
            spectime  = 0
            leftarm = get_descendants(tree, i.clades[0], includeself=True)
            rightarm = get_descendants(tree, i.clades[1], includeself=True)
            perms = [frozenset([x, y]) for x in leftarm for y in rightarm]
            pars = []
            matchedvals = [item for sublist in list(matched_dic.values()) for item in sublist]
            for j in perms:
                if j in matchedvals:
                    dist = Phylo.BaseTree.TreeMixin.distance(tree, list(j)[0], list(j)[1])
                    dists.append(dist)
                    pars.append(j)
            if i in [x[0][1] for x in trim_score]:
                myspec = "Yes"
                spectime = [x[2] for x in trim_score if x[0][1] == i][0] #wow i'm bad at programming
            mydist = np.mean(dists)
            line = clade_object_to_readable(i, tree) + ',' + str(mydist) + ',' + str(len(dists)) + ',' + myspec + ',' + str(spectime) + "," + str([(clade_object_to_readable(list(x)[0], tree), clade_object_to_readable(list(x)[1], tree)) for x in pars]) + '\n'
            f.write(line)
        f.close()


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
            f.write(str(line).replace("'","")[1:-1] + '\n')
    f.close()
    return


def pars_and_ancs(parscores, i, tree, data_dic):
    parsandanc_onspec =[]
    trace = Phylo.BaseTree.TreeMixin.trace(tree, list(i[0][0])[0], list(i[0][0])[1])
    trace.append(list(i[0][0])[0])
    for j in list(itertools.combinations(trace, 2)):
        for k in parscores:
            if frozenset([j[0], j[1]]) == k[0]:
                intscore = -1
                try:
                    intscore = data_dic[k[0]]
                except:
                    pass
                #parsandanc_onspec.append([[k[0], i[0][1]], intscore, 0.0])
    ancscore1, ancscore2 = -1, -1
    try:
        ancscore1 = data_dic[frozenset([list(i[0][0])[0], i[0][1]])]
    except:
        pass
    try:
        ancscore2 = data_dic[frozenset([list(i[0][0])[1], i[0][1]])]
    except:
        pass
    parsandanc_onspec.append([[frozenset([list(i[0][0])[0], i[0][1]]), i[0][1]], ancscore1, 0.0])
    parsandanc_onspec.append([[frozenset([list(i[0][0])[1], i[0][1]]), i[0][1]], ancscore2, 0.0])
    return parsandanc_onspec


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
            if list(i[0][0])[0] not in input_tree.get_terminals():
                xpep = Phylo.BaseTree.TreeMixin.get_nonterminals(input_tree).index(list(i[0][0])[0]) + 172
            else:
                xpep = str(list(i[0][0])[0]).lower()
            if list(i[0][0])[1] not in input_tree.get_terminals():
                ypep = Phylo.BaseTree.TreeMixin.get_nonterminals(input_tree).index(list(i[0][0])[1]) + 172
            else:
                ypep = str(list(i[0][0])[1]).lower()
            ancestor = Phylo.BaseTree.TreeMixin.get_nonterminals(input_tree).index(i[0][1]) + 172
            line = str(xpep) + ',' + str(ypep) + "," + str(ancestor) + "," + str(i[1]) + "," + str(i[2]) + '\n'
            f.write(line)
    f.close()
    return


def trace_specificity_gain(intscore, tree, nodes, tips, coors, data_dic, count, descendants = False):
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
    # if descendants:
    #     para1terms = trace[0].get_terminals()
    #     para1desc = []
    #     for term1 in para1terms:
    #         paratrace = Phylo.BaseTree.TreeMixin.trace(tree, trace[0], term1)
    #         for clade in paratrace:
    #             if clade not in para1desc:
    #                 para1desc.append(clade)
    #     para2terms = trace[-1].get_terminals()
    #     para2desc = []
    #     for term2 in para2terms:
    #         paratrace2 = Phylo.BaseTree.TreeMixin.trace(tree, trace[-1], term2)
    #         for clade in paratrace2:
    #             if clade not in para2desc:
    #                 para2desc.append(clade)
    #     for pair in [[x, y] for x in trace for y in para1desc]:
    #         try:
    #             interact = data_dic[frozenset([pair[0], pair[1]])]
    #             tracedescend.append([[frozenset([pair[0], pair[1]]), 0], interact, 0])
    #         except KeyError:
    #             pass
    #     for pair in [[x, y] for x in trace for y in para2desc]:
    #         try:
    #             interact = data_dic[frozenset([pair[0], pair[1]])]
    #             tracedescend.append([[frozenset([pair[0], pair[1]]), 0], interact, 0])
    #         except KeyError:
    #             pass
    #     for i in para1desc:
    #         try:
    #             interact = data_dic[frozenset([i])]
    #             tracedescend.append([[frozenset([i]), 0], interact, 0])
    #         except KeyError:
    #             try:
    #                 interact = data_dic[frozenset([str(i).lower()])]
    #                 tracedescend.append([[frozenset([i]), 0], interact, 0])
    #             except KeyError:
    #                 tracedescend.append([[frozenset([i]), 0], -1, 0])
    #     for i in para2desc:
    #         try:
    #             interact = data_dic[frozenset([i])]
    #             tracedescend.append([[frozenset([i]), 0], interact, 0])
    #         except KeyError:
    #             try:
    #                 interact = data_dic[frozenset([str(i).lower()])]
    #                 tracedescend.append([[frozenset([i]), 0], interact, 0])
    #             except KeyError:
    #                 tracedescend.append([[frozenset([i]), 0], -1, 0])
    return spectrace
    #all_paralog_plot_wrapper(tree, nodes, tips, coors, spectrace, filename, tracedescend, [], False)


if __name__ == "__main__":
    mydata_dic = load_experiment('../200120_most_aliases.csv')
    mytree, mytips, mynodes = read_tree('../190918_EA_tree1only.txt')
    myduplication_nodes = find_duplication_nodes(mytree, '../Gene_duplications.txt')
    myspecies = get_species(mytips)
    myspecies_tree, species_tips, species_nodes = read_tree('../speciestree_node_names.newick.nwk')
    mymatched_ancestors, mymatched_dic = find_matched_ancestors(myspecies, myspecies_tree, mytree, myduplication_nodes)
    mypars = paralog_ancestors(mymatched_dic, mytree)
    mynewdict = convert_data_dic_to_species_names(mydata_dic, mynodes, mytips)
    myscores = paralog_timing(mypars, mynewdict, mytree)
    mycoors = get_tree_coordinates(mytree, mynodes, mytips)

    #info_on_specificity_gains(myscores)
    mytrimmedscores = basal_paralog_int_loss(myscores, mytree)

    print("pauser")
    #rectangle paths showing where specificity occurred
    all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, mytrimmedscores, '../figures/200505_traced_specificity_gains.pdf')
    #spaghettigrams for each path
    for mycount, score in enumerate(mytrimmedscores):
       specs = trace_specificity_gain(score, mytree, mynodes, mytips, mycoors, mynewdict, mycount, descendants=False)
       all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, specs, str('../figures/200506_traced_specgain' + clade_object_to_readable(list(score)[0], mytree)+ clade_object_to_readable(list(score)[1], mytree) + '.pdf'), [], [], False)

    for mycount, score in enumerate(mytrimmedscores):
       specs = trace_specificity_gain(score, mytree, mynodes, mytips, mycoors, mynewdict, mycount, descendants=False)
       parspecs = {}
       for key in specs:
           if len(key) == 1 or key in mypars:
               parspecs[key] = specs[key]
       all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, parspecs, str('../figures/200506_spaghetti_paralogs' + clade_object_to_readable(list(score)[0], mytree)+ clade_object_to_readable(list(score)[1], mytree) + '.pdf'), [], [], False)

    # plotscores = pars_and_ancs(mypars, mytrimmedscores, mytree, mynewdict)
    # allints = []
    # for i in mynewdict:
    #     allints.append([[i, 'anc'], mynewdict[i], 'dist'])
    # plotscores = []
    # for c, i in enumerate(mytrimmedscores):
    #     plotscores.extend(pars_and_ancs(mypars, i, mytree, mynewdict))
    #     #all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, plotscores, '../200220_pars_withspecgain'+str(c)+'.pdf', [], [])
    #
    # write_scores(plotscores, "200220_anctoregains.csv", mytree)

    #mytrimmedscores = dupes_intscores(myduplication_nodes, mytree, mynewdict)
    #write_possible_paralogs(mytree, myduplication_nodes, mytrimmedscores, mymatched_dic)
    #write_odds_of_regain(mytree, mytrimmedscores, mymatched_dic, mynewdict)
    # altdict = {}
    # with open('../Altalls.fasta') as f:
    #     counter = 0
    #     for line in f:
    #         name = line.rstrip()
    #         seq = f.readline().rstrip()
    #         name = name[1:-6]
    #         for k in name.split(sep="+"):
    #             altdict[frozenset([k])] = seq
    # altdictoclade = convert_data_dic_to_species_names(altdict, mynodes, mytips)
    #
    # aliasdic = {}
    # with open("../sequence_aliases.txt") as f:
    #     for line in f:
    #         name, alias = line.rstrip().split()
    #         for k in alias.split(sep ="+"):
    #             aliasdic[frozenset([k])] = name
    #
    # nodedict = {}
    # with open('../Ancestors.fasta') as f:
    #     counter = 0
    #     for line in f:
    #         name = line.rstrip()
    #         seq = f.readline().rstrip()
    #         name = name[1:]
    #         for k in name.split(sep="+"):
    #             nodedict[frozenset([k])] = seq
    #             if frozenset([k]) in aliasdic.keys():
    #                 nodedict[frozenset([aliasdic[k]])] = seq
    # nodedictoclade = convert_data_dic_to_species_names(nodedict, mynodes, mytips)
    # #print("altdic", altdict)
    # alttop = {}
    # with open('../alt_top.fasta') as f:
    #     for line in f:
    #         name = line.rstrip()
    #         seq = f.readline().rstrip()
    #         name = name[1:-7]
    #         for k in name.split(sep="+"):
    #             alttop[frozenset([k])] = seq
    # alttopclade = convert_data_dic_to_species_names(alttop, mynodes, mytips)
    #
    # altpairs = {}
    # with open('../Sequence_pairs_cleaned.txt') as f:
    #     for line in f:
    #         if 'altphy' in line:
    #             xpair, ypair = line.rstrip().split(sep=',')
    #             xpair = xpair[:-6]
    #             ypair = ypair[:-6]
    #             altpairs[frozenset([xpair, ypair])] = "A"
    # alttoppairs = convert_data_dic_to_species_names(altpairs, mynodes, mytips)
    #
    #
    # fig = plt.figure(figsize=(18, 18))
    # ax = fig.add_subplot(111)
    # draw_tree(mytree, mynodes, mytips, mycoors, ax, "gene")
    # # for i in altdictoclade.keys():
    # #     co = mycoors[list(i)[0]]
    # #     ax.scatter(co[0], co[1] +.3, s = 20,color= 'red', zorder = 7, edgecolor="black")
    # # for i in nodedictoclade.keys():
    # #     co = mycoors[list(i)[0]]
    # #     ax.scatter(co[0], co[1] - .3, s = 20,color= 'blue', zorder = 7, edgecolor="black")
    # for i in alttopclade.keys():
    #     co = mycoors[list(i)[0]]
    #     ax.scatter(co[0], co[1]+.3, s = 35, color = 'gold', zorder = 7, edgecolors="black")
    # for i in alttoppairs.keys():
    #     if len(list(i)) > 1:
    #         cox1, coy1 = mycoors[list(i)[0]]
    #         cox2, coy2 = mycoors[list(i)[1]]
    #         xarc, yarc = parabola([cox1, coy1], [cox2, coy2], .25)
    #         ax.plot(xarc, yarc, color='teal', alpha=0.4, zorder=4)
    #     else:
    #         co = mycoors[list(i)[0]]
    #         ax.scatter(co[0], co[1]-.3, s=35, color='teal', zorder=7, edgecolors="black")
    # ax.plot([1.5, 2], [-160, -160], color = 'black')
    # ax.set_axis_off()
    # figure = plt.gcf()
    # #plt.show()
    # figure.savefig('../200303_alttoppology.pdf', figsize=(100, 30))
    # plt.close()



    #myduplication_nodes.reverse()
    # for mycount, mydupnode in enumerate(myduplication_nodes):
    #    mytrimmedscores = dupes_intscores(mydupnode, mytree, mynewdict)
    #    filename = "../figures/200128_duplicationnode_trace" + clade_object_to_readable(mydupnode, mytree) + '.pdf'
    #    all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, mytrimmedscores, filename, [], False)


    #alldups, allscorestrim = specificity_after_duplication(mytree, mynewdict, myduplication_nodes)
    #print(alldups)
    #write_scores(myscores, "../200205_allparalogs.csv", mytree)

    #myhomos = just_get_homos(mynewdict, mytree)
    #all_paralog_plot_wrapper(mytree, mynodes, mytips, mycoors, myhomos, '../200120_homodimers.pdf')
