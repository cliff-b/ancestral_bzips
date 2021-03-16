# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from Bio import Phylo
import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy.stats
import random

def read_tree(treefile, silent=True):
    '''Takes a phylogenetic tree file
    reads it in and returns the entire tree as well as
    nodes and tips separately
    Takes: 
        treefile: filename of a rooted tree in .newick format
        silent: Boolean, set to True to get the number of nodes and tips in the tree
        
    Returns:
        tree: a tree object
        nodes: a list of clades that are the nodes
        tips: a list of clades that are the tips
    '''
    tree=Phylo.read(treefile, 'newick', rooted=True)
    nodes=Phylo.BaseTree.TreeMixin.get_nonterminals(tree)
    tips=Phylo.BaseTree.TreeMixin.get_terminals(tree)
    if silent==False:
        print('Reading in tree..')
        print(treefile+' contains '+str(len(nodes))+' nodes')
        print(treefile+' contains '+str(len(tips))+' taxa')
    return tree, tips, nodes

def load_experiment(filename):
    '''Loads a file in .csv format with interaction scores in column
    twelve. Makes a dictionary of sets of the peptides to link to the 
    interaction score. Theoretically used to fix aliasing present both in names and hidden from names'
    Takes:
        filename: string, of a file with peptides in cols 1, 2 and interactions in 12
    returns:
         data_dic:   a dictionary of frozen sets of strings, key =frozenset({'pep1', 'pep2'}) value = interaction score'''
    f=open(filename,'r')
    #alias_dic= {}#get_alias_dic()
    lines=f.readlines() #[0].split('\r')
    data_dic={}
#    print("so many lines", lines[1:10])
    for i in lines[1:]:
    
        line=i.split()[0].split(',') 
            
        if len(line[0].split('+'))>1 or len(line[1].split('+'))>1:
                for j in line[0].split('+'):
                    for k in line[1].split('+'):
#                        
#                        if j.lower() in  alias_dic.keys() and k.lower() in  alias_dic.keys():
#            
#                            for l in alias_dic[j.lower()]:
#                                for m in alias_dic[k.lower()]:
#                                    try:
#                                        data_dic[frozenset([l,m])]= np.average([data_dic[frozenset([l,m])],float(line[12])])
#                                    except KeyError:
#                                        data_dic[frozenset([l,m])]=float(line[12])                                      
#                        elif j.lower() in  alias_dic.keys():
#                            for l in alias_dic[j.lower()]:
#                                try:
#                                    data_dic[frozenset([l,k.lower()])]= np.average([data_dic[frozenset([l,line[1].lower()])],float(line[12])])
#                                except KeyError:
#                                    data_dic[frozenset([l,k.lower()])]=float(line[12])
#                        elif line[1].lower() in  alias_dic.keys():
#                            for l in alias_dic[line[1].lower()]:
#                                try:
#                                    data_dic[frozenset([j.lower(),l])]= np.average([data_dic[frozenset([j.lower(),l])],float(line[12])])
#                                except KeyError:
#                                    data_dic[frozenset([j.lower(),l])]=float(line[12])                

                        try:
                            data_dic[frozenset([j.lower(),k.lower()])]= np.average([data_dic[frozenset([j.lower(),k.lower()])],float(line[12])])
                        except KeyError:
                            data_dic[frozenset([j.lower(),k.lower()])]=float(line[12])
        else:
#            if False:
#        
#
##            if line[1].lower() in  alias_dic.keys() and line[0].lower() in  alias_dic.keys():
##
##                for j in alias_dic[line[0].lower()]:
##                    for k in alias_dic[line[1].lower()]:
##                        try:
##                            data_dic[frozenset([j,k])]= np.average([data_dic[frozenset([j,k])],float(line[12])])
##                        except KeyError:
##                            data_dic[frozenset([j,k])]=float(line[12])                                      
##            elif line[0].lower() in  alias_dic.keys():
##                for j in alias_dic[line[0].lower()]:
##                    try:
##                        data_dic[frozenset([j,line[1].lower()])]= np.average([data_dic[frozenset([j,line[1].lower()])],float(line[12])])
##                    except KeyError:
##                        data_dic[frozenset([j,line[1].lower()])]=float(line[12])
##            elif line[1].lower() in  alias_dic.keys():
##                for j in alias_dic[line[1].lower()]:
#                    try:
#                        data_dic[frozenset([line[0].lower(),j])]= np.average([data_dic[frozenset([line[0].lower(),j])],float(line[12])])
#                    except KeyError:
#                        data_dic[frozenset([line[0].lower(),j])]=float(line[12])
          

            try:
                data_dic[frozenset([line[0].lower(),line[1].lower()])]= np.average([data_dic[frozenset([line[0].lower(),line[1].lower()])],float(line[12])])
            except KeyError:
                data_dic[frozenset([line[0].lower(),line[1].lower()])]=float(line[12])
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

#def get_alias_dic():   
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


#def read_ancestral_seqs(treename, silent=True):
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


def read_alignment(ali_file,seq_range):
    '''Opens alignment in phylip format (easier than fasta)
        trims the sequence to just the AAs in seq_range
        self note: #90:131 is bzip domain  #65:131 is what I use for ASR, with 8cat
        Takes: 
            ali_file: file in phylip format
            seq_range: a numeric two item list 
        Returns: 
            dictionary of key = species, values = sequences'''
    seq_dict={}
    f=open(ali_file,'r')
    lines=f.readlines()
    f.close()
    for i in lines[1:]: #Ignore first line in Phylip format
        seq_dict[i.split()[0]]=i.split()[1][seq_range[0]:seq_range[1]]
        if 'X' in i.split()[1][seq_range[0]:seq_range[1]]:
            print( i.split()[0]+' has an X amino acid somewhere')
    
    return seq_dict    

def find_duplication_nodes(tree, duplication_file):
    '''function to read in a list of extant proteins with a duplicated ancestor
    list is of form 'pep1,pep2', funciton locates the last common ancestor on the tree
    and returns a list of those clades
    Takes:
        tree: a tree object
        duplication_file: filename of a file containing proteins that were duplicated on the tree object
    Returns:
       duplication_nodes: a list of clade objects
    '''
    #currently the duplication nodes returned are 173 269 268 310 317 323 248 239 248 174 175 200 201 195 190 296 294 293 320 219 211 203
    duplication_nodes=[]
    
    f=open(duplication_file, 'r')
    lines=f.readlines()
    f.close()
    pairs=lines[0].split()
    print(pairs)
    for i in pairs:
        duplication_nodes.append(Phylo.BaseTree.TreeMixin.common_ancestor(tree,[i.split(',')[0], i.split(',')[1]]))
    return duplication_nodes    


def find_orthologs(tree, data_dic, tip_dic, duplication_list):
    '''Takes a tree and interaction data and gives back a list of orthologs
        importantly these orthologs are only ones where the homodimers of the two
        interact and the ancestral proteins homodimers interact
    Takes:
        data_dic: a dictionary frozensets to floats of interaction data in the form key = frozenset({'173', '174'}), value = 1.0
        tip_dic: dictionary of strings to strings to convert names from lowercase in the data_dic to mixedcase in the tree
        duplication_list: a list of clade objects that were points of gene duplication
    Returns:
        ortho_dic: a dictionary of key = frozenset and value = list of two numbers, distance between orthologs and whether they interact or not '''
    print ('Finding orthologs')
    ortho_dic={}
    yes, no, missingpoints =0, 0 ,0
    for i in data_dic.keys():
        if len(list(i))>1:
           protein1,protein2=list(i)
           ancestor=Phylo.BaseTree.TreeMixin.common_ancestor(tree, tip_dic[protein1],tip_dic[protein2])
           #This is a gross way to get from tree to tip_dic name
           anc_in_data_dic = frozenset([list(tip_dic.keys())[list(tip_dic.values()).index(ancestor)]])
           try: 
               # check that the homodimers for protein1, protein2 and ancestor all homodimerize
                if ancestor not in duplication_list and data_dic[anc_in_data_dic]==1 and data_dic[frozenset([protein1])]==1 and data_dic[frozenset([protein2])]==1:
                    if data_dic[i]==1:
                        ortho_dic[i]=[Phylo.BaseTree.TreeMixin.distance(tree,  tip_dic[protein1],tip_dic[protein2]), 1]
                        yes+=1
                    elif data_dic[i]==0:
                        ortho_dic[i]=[Phylo.BaseTree.TreeMixin.distance(tree,  tip_dic[protein1],tip_dic[protein2]),0]
                        no+=1
           except KeyError:
                missingpoints += 1
    print ('Done finding orthologs','\nNumber of noninteractors: ', no,'\nNumber of interactors:', yes)
    print("Missing points: ", missingpoints)
    return ortho_dic

def pick_nonoverlap(tree, pairs, tip_dic):
    '''finds a set of proteins from the same tree that do not share any branches on the tree
    Takes:
        tree: a tree object
        pairs: a dictionary of protein pairs on the tree, key = frozenset({'173', '174'})
    '''
    #print("the pairs",pairs)
    #tip_dic=tip_dict(tips)
    #i don't think probs_dic does things
    #probs_dic={}
    #keys=[]
    print("its britney bitch")
    for i in pairs.keys():
        #print("current i",i)
        A= tip_dic[list(i)[0]]
        B= tip_dic[list(i)[1]]
        #I'm pretty sure this is an error and we're adding B twice to the list and rather we need to have +[A]
        #pairs[i].append(Phylo.BaseTree.TreeMixin.trace(tree,A,B)+[B])
        pairs[i].append(Phylo.BaseTree.TreeMixin.trace(tree,A,B)+[A])
        #keys.append(i)
        #probs_dic[i]=pairs[i][0]
#    print(probs_dic)
    values=[]
    distances=[]
    keys=list(pairs.keys())
#    tot=sum(probs_dic.values())
#    probs=np.array([probs_dic[i]/float(tot) for i in keys])
#    probs=probs/sum(probs)
    #keys is a list of what accesses dictionary pairs
    while len(keys)>0:
        print("length keys: ", len(keys))#, print(keys))
        #print("alll keys",keys)
        #choose a random key from the list of keys that access the dictionary
        key=np.random.choice(keys)#,p=probs)
        #key = random.choice(keys)
        print("chosen key" ,key)
        #add the interaction score for the value to the list of values
        values.append(pairs[key][1])
        #add the distance between the two proteins of the pair to the list of distances
        distances.append(pairs[key][0])
        #what the fuck are we doing here
        #values.append(1)
        #distances.append(0)
        #keys
        #Here we test for overlap if i in keys does not have more than one shared node with key then we don't delete it from the dicitonary
        for i in keys:
            try:
                #if len(intersection(pairs[i][2],pairs[key][2]))>1 and i!=key:
                if len(intersection(pairs[i][2],pairs[key][2]))>0 and i!=key:
                    del pairs[i]
                    #probs_dic[i]
            except KeyError:
                    1
        del pairs[key]
#        del probs_dic[key]
        keys=list(pairs.keys())
#        tot=sum(probs_dic.values())
        #probs=np.array([probs_dic[i]/float(tot) for i in keys])**2
        #probs=probs/sum(probs)
    print("distances: ", distances, "values: ", values)
    return values, distances    


def tip_dict(tip_list, node_list):
    '''this is necessary to deal with the annoying formatting of the names
    takes in a list of clade objects of the tips of a tree, and either converts the name to
    lower case if it's a tip or a number if its the node
    Takes:
        tip_list: a list of clade objects that are tips ie Clade(branch_length=0.041194, name='HLFCryptotermes_secundus')
        node_list: a list of clade objects that are not tips ie Clade(branch_length=0.025239)
    Returns:
        name_dic: a dictionary that converts tips between key: uncapitalized and value: capitalized 
        and nodes between their key: number and value: clade object
        NB: why do we have a dictionary with values that are two different objects?'''
    name_dic={}
    name_dic2={}
    for i in tip_list:
       name_dic[str(i).lower()]=str(i)
       name_dic2[str(i)] = str(i).lower()
    for c,i in enumerate(node_list):
        name_dic[str(c+172)]=i
        name_dic2[i] = str(c+172)
    return name_dic, name_dic2

def neutral_evo(tree, data_dic, tip_dic, duplication_nodes):
    bins=[0.0,0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.5,2,2.5,3]
    filtered_orthologs=find_orthologs(tree, data_dic, tip_dic, duplication_nodes)
    distances_tot=[]
    values_tot=[]
    pairs=[]
    plt.figure()

    for i in filtered_orthologs.keys():
        dist,val=filtered_orthologs[i]
        distances_tot.append(dist)
        values_tot.append(val)
        #distances_tot.append(0)
        #values_tot.append(1)
        pairs.append(i)
        #pairs.append(i)

    bin_contents=[]
    bin_values=[]
    for i in range(len(bins)+1):
        bin_contents.append([])    
        bin_values.append([]) 
    indices=np.digitize(distances_tot,bins)
    print("the most indicies was :", max(indices))
    for c,i in enumerate(indices):
        bin_contents[i].append(pairs[c])
        bin_values[i].append(values_tot[c])
    average_val=[]
    average_dist=[]
    average_no=[]
    tot_yerr_pos=[]
    tot_yerr_neg=[]
    samples=20
    
    bin_average=[np.mean(i) for i in bin_values[1:]]
    print("the bin average was: ", bin_average)
    plt.plot( bins, bin_average, 'o', color='grey', alpha=0.1)
    plt.figure()
    for s in range(samples):
        random_val=[]
        random_dist=[]
        random_no=[]
        yerr_pos=[]
        yerr_neg=[]
        for i in bin_contents:
            new_dic={}
            for j in i:
                new_dic[j]=filtered_orthologs[j]
            values,distances= pick_nonoverlap(tree,new_dic, tip_dic)
            random_val.append(1-np.mean(values))
            random_dist.append(np.mean(distances))
            random_no.append(len(values))
            errors=clopper_pearson(np.sum(values),random_no[-1],alpha=0.32)
            yerr_pos.append(1-errors[1])
            yerr_neg.append(1-errors[0])
        #plt.plot( bins, random_val[1:], 'o', color='grey', alpha=0.01)

        tot_yerr_pos.append(yerr_pos)
        tot_yerr_neg.append(yerr_neg)  
        average_val.append(random_val)
        average_dist.append(random_dist)
        average_no.append(random_no)
    #t=[]
    val=[]
    xs=[]
    for c,i in enumerate(bins[:-1]):
        xs.append(i+(bins[c+1]-i)/2.0)
    xs.append(bins[-1])
    result=scipy.optimize.leastsq(rev_firs_order_error, [0.2, 0.2], args=(xs[:8], np.mean(average_val, axis=0)[1:9], np.mean(average_no, axis=0)[1:9]))
    print("the result was: ", result)
    print("the sum was: ", sum(rev_firs_order_error(result[0],xs[:8], np.mean(average_val, axis=0)[1:9], np.mean(average_no, axis=0)[1:9])))
    #plt.plot(xs,rev_firs_order(np.array(xs),result[0][0],result[0][1]))
    plt.errorbar( xs, np.mean(average_val, axis=0)[1:],yerr=[np.mean(average_val, axis=0)[1:]-np.nanmean(tot_yerr_neg, axis=0)[1:], np.nanmean(tot_yerr_pos, axis=0)[1:]-np.mean(average_val, axis=0)[1:]], marker='o', linestyle='None')
    plt.xlabel('Branch length')
    plt.ylabel("Percent of interactions lost")
    plt.figure()
    plt.errorbar(np.mean(average_dist, axis=0), np.mean(average_no, axis=0),yerr=np.std(average_val, axis=0))
  

    plt.show()


def rev_firs_order_error(params,t, values, weights):
    kf,kb=params
    weights=[ 1 for i in t]
    model=np.array([rev_firs_order(i,kf,kb) for i in t])
    error=[]
    for c, i in enumerate(weights):
        error.append((values[c]-model[c])**2 *i) 
    return error   

def rev_firs_order(t,kf,kb):
    '''reversible first order model'''    
    A0=1
    kb=0
    Aeq=(kb*A0)/(kf+kb)
    A=1-(A0-(Aeq))*np.exp(-1*(kf+kb)*t)-(Aeq)
    return A

def clopper_pearson(k,n,alpha=0.32):
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

if __name__ =="__main__":
    #print("in main")
    tree, tips, nodes=read_tree('../190918_EA_tree1only.txt')
    tip_dic,tips2 = tip_dict(tips, nodes)
    #print(tips)
    #print(nodes)
    #someday we'll get there
    data_dic = load_experiment('../190919_medianEA66k-1.csv') 
    duplication_nodes=find_duplication_nodes(tree,'../Gene_duplications.txt')
    #orthos =find_orthologs(tree, data_dic, tip_dic, duplication_nodes)
    #nonovers = pick_nonoverlap(tree,orthos,tip_dic)
    neutral_evo(tree, data_dic, tip_dic, duplication_nodes)
    #print(orthos)
    #print(duplication_nodes)
    
    #print(data_dic)
    #print(tips)
    for i in range(10):
        print("Heres some data",list(data_dic)[i], " = ", data_dic[list(data_dic)[i]])
    print("length of dict", len(data_dic))
    aligns=read_alignment('../reformatted_alignment.phy',[25,66])
    #print(aligns.keys())
#    ancestral_seqs,altall_seqs,pp=read_ancestral_seqs('../MLtree', silent = False)  
#    print("ancestral seqs: ", ancestral_seqs)
#    print("altall_seqs: ", altall_seqs)
#    print("pp: ", pp)
#    alias = get_alias_dic()
#    print("\n alias is: ", alias