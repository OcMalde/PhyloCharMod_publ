#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Load and interpret a seadog .output reconciliation file and build a module presence map during gene evolution
"""

import argparse
import os
import random
import time
import requests
import shutil
import concurrent.futures
import matplotlib as plt
from pathlib import Path
from natsort import natsorted
from ete3 import Tree
from itertools import chain
from dna_features_viewer import GraphicFeature, GraphicRecord

#==============================================================================
# Conserved module class 
#==============================================================================

class c_module:
    """
    Conserved module class
    Represents a module and build his itol string
    
    Attributes
    ----------
    name : str
        Name of the module
    start : str
        Start position of the module
    end : str
        End position of the module
    gene : str
        Name of the gene where the module segment come from
    module_node_name : str
        Name of the gene node corresponding
    shape : str
        Type of shape associated to the module in itol visualisation
    color : str
        Color associated to the module in itol visualisation
    
    Methods
    -------
    get_shapeColor 
        Associate a `shape` and a `color` to the module
    itol_str
        Generate an itol string formated in itol domain dataset format (unique shape/color)
    itol_str_gain
        Generate an itol string formated in itol domain dataset format (green box)
    itol_str_lost
        Generate an itol string formated in itol domain dataset format (red box)
        
    """

    def __init__(self, name, start, end, gene, module_node_name):
        # Known
        self.name = name
        self.start = start
        self.end = end
        self.gene = gene
        self.module_node_name = module_node_name
        # To define
        self.shape = None
        self.color = None

    def __repr__(self):
        return f"{self.module_node_name}"

    def __str__(self):
        return f"{self.module_node_name}"

    def get_shapeColor(self, dict_moduleShapeColor):
        self.shape = dict_moduleShapeColor[self.name][0]
        self.color = dict_moduleShapeColor[self.name][1]

    def itol_str(self) -> str:
        if self.start == None: 
            return ""
        else: 
            return f"{self.shape}|{round(int(self.start),2)}|{round(int(self.end),2)}|{self.color}|{self.name}"
    
    def itol_str_gain(self) -> str:
        if self.start == None: 
            return ""
        else: 
            return f"RE|{round(int(self.start),2)}|{round(int(self.end),2)}|#00FF00|{self.name}"
        
    def itol_str_lost(self) -> str:
        if self.start == None: 
            return ""
        else: 
            return f"RE|{round(int(self.start),2)}|{round(int(self.end),2)}|#8B0000|{self.name}"

#==============================================================================
# Make dict module : (shape,color)
#==============================================================================

def make_dict_Module_ShapeColor(module_list) -> dict:
    """
    From a module list, make an attribution of a shape and color for each module
    
    Parameters
    ----------
    module_list : list
        A list of modules
    """
    shape_list = ["RE","HH","HV",
                    "EL","DI","TR",
                    "TL","PL","PR",
                    "PU","PD","OC",
                    "GP"]
    return {module : (random.choice(shape_list), random_hexaColor())
                for module in module_list
            }

def make_dict_Domains_ShapeColor(domain_list) -> dict:
    """
    From a domain list, make attribution of simple defined grey shape of each domain
    
    Parameters
    ----------
    domain_list : list
        A list of domains
        
    Returns
    -------
    dict_domain_shape_color : dict
        A dictionary with domain name as key and shape/color (tuple of str) as value
    """
    dict_domain_shape_color = {}
    for domain in domain_list:
        # All are ligh grey
        color = "#D3D3D3"
        # TSP are elispes
        if domain == "PS50092" or domain == "PF00090" or domain == "PF19030": shape = "EL" 
        # Plac are diamond
        elif domain == "PS50900": shape = "DI"
        # CUB are triangle
        elif domain == "PS01180": shape = "TR"
        # Ig are pentagram
        elif domain == "PS50835": shape = "PU"
        # GON domain are left pointing triangle
        elif domain == "PF08685": shape = "TL"
        # Kunitz and associated domains are gap
        elif domain == "PF00014" or domain == "PF16626": shape = "GP"
        # NTR (like 5) are gap
        elif domain == "PF01759": shape = "GP"
        # If nothing specified, shape is rectangle
        else: shape = "RE"
        dict_domain_shape_color[domain] = (shape, color)
    return dict_domain_shape_color

def random_hexaColor() -> str:
    """
    Generate an string with a random hexadecimal color code
    
    Returns
    -------
    col : str
        A hexadecimal string corresponding a random color
    """
    col = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
    return col

#==============================================================================
# Reconciliation mapping class
#==============================================================================

class c_mapping:
    """
    Reconciliation mapping class
    Represents any DGS mapping
    DG or GS
    
    Attributes
    ----------
    entity_1 : str
        First element described in the mapping (i.e., a tree node)
    entity_2 : str
        Second element described in the mapping (i.e., a node of another tree)
    event : str
        Type of event associated to the mapping
    
    Methods
    -------
    get_geneNode 
        Get the node in the gene tree, corresponding a gene entity (could be `entity_1` or `entity_2`)
    get_modNode 
        Get the node in the module tree, corresponding a module entity (could be `entity_1` or `entity_2`)
        
    """

    def __init__(self, seadog_mapping_line):
        # Split and stock the elements of the raw seadog mapping line
        # Stock : entity1 / entity2 / kind of event
        self.entity_1 = seadog_mapping_line.split(":")[0]
        self.entity_2 = seadog_mapping_line.split("->")[1].split(" ")[1]
        self.event = seadog_mapping_line.split(":")[1].split(",")[0][1:]
        # If the event is a transfer, we add a recipient
        if "transfer" in self.event:
            self.recipient = seadog_mapping_line.split("Recipient ->")[1].split(" ")[1]
    
    def __repr__(self):
        return f"{self.entity_1} <== {self.event} ==> {self.entity_2}"

    def __str__(self):
        return f"{self.entity_1} <== {self.event} ==> {self.entity_2}"

    def get_geneNode(self, gene_tree):
        """
        Check if one of the entity is a gene, if its the case, return his node
        """
        for node in gene_tree.traverse("postorder"):
            if self.entity_1 in node.name:
                return node
            elif self.entity_2 in node.name:
                return node

    def get_modNode(self, dic_module_modTree):
        """
        Check if one of the entity is a module, if its the case, return his node
        """
        for module, module_tree in dic_module_modTree.items():
            for node in module_tree.traverse():
                if self.entity_1 in node.name:
                    return node
                elif self.entity_2 in node.name:
                    return node

#==============================================================================
# Read Seadog output
#==============================================================================

def read_seadogO(file_name, gene_tree_file) -> tuple:
    """
    Read the seadog output file and build :
        - a gene Tree
        - a dict {gene : [domain_list]}
        
    Parameters
    ----------
    file_name : str
        Name of a txt file formated as an SEADOG-MD output file
    gene_tree_file : str
        Name of the gene tree file corresponding, in newick format
        
    Returns
    -------
    gene_tree : Tree
        Gene tree as a ete3.Tree object
    dict_module_mappingList : dict
        Dictionary with module name as key and list of mapping (c_mapping object) as value
    dict_module_modTree : dict
        Dictionary with module name as key and module tree (ete3.Tree) as value
    gs_mapping_list : list
        A list of mapping object from Gene-Species reconciliation
    """
    
    dict_gene_moduleList = {}
    dict_gene_spGeneEvent = {}
    dict_module_modTree = {}
    dict_module_mappingList = {}
    all_module_list = []
    list_module_gene_recipient = []
    reconc_check = (0,0)

    with open(file_name, "r") as so_file:
        for line in so_file:
            # Extract the differents modules trees
            #(B13|265|269|_NP001193502.1_0,B13|265|269|_NP891550.1_1)D2_0_1;
            if line.startswith("(") and "D" in line:
                module_tree = Tree(f"({line.replace(';','')});", format=1)
                module_name = module_tree.get_leaves()[0].name.split("|")[0]
                # Add it to the dict
                dict_module_modTree[module_name] = module_tree
            # Load the gene tree as a ete3 tree object
            elif line.startswith("(") and ")G" in line:
                reconc_gene_tree = Tree(f"({line.replace(';','')});", format=1)
                length_gene_tree = Tree(str(gene_tree_file), format=1)#8
                gene_tree = full_gene_tree(reconc_gene_tree, length_gene_tree)
            # Indicator to know in which reconciliation is the actual mapping line
            elif "Reconciliation between gene trees and species tree:" in line:
                reconc_check = (1,1)
                # Save last module mappings in their dict
                module = [mapp.entity_1 for mapp in mg_mapping_list if len(mapp.entity_1.split("|")) > 1][0].split("|")[0]
                dict_module_mappingList[module] = mg_mapping_list
                # Mapping list for gene_species
                gs_mapping_list = []
            elif "Reconciliation between domain trees and gene trees:" in line:
                reconc_check = (1,0)
                # Init first mapping list
                mg_mapping_list = []
            
            #====================================================================================
            # Gene-moduel reconciliations infos
            #====================================================================================
            # Extract mappings, per modules
            # We want to know when we change of module mappings
            if "Domain Tree:" in line and reconc_check == (1,0) and len(mg_mapping_list) > 0:
                # Get the module name
                module = [mapp.entity_1 for mapp in mg_mapping_list if len(mapp.entity_1.split("|")) > 1][0].split("|")[0]
                # We change of module mappings, so save the finished modules mapping in the dict
                dict_module_mappingList[module] = mg_mapping_list
                # And then reset the mapping list for the next module
                mg_mapping_list = []
            # For this, use the Mappings lines
            elif "Mapping" in line and reconc_check == (1,0):
                # Create an mapping object for the mapping line
                mapping_object = c_mapping(line)
                # Add it to the all mapping list
                mg_mapping_list.append(mapping_object)
            
            #====================================================================================
            # Gene Specie Reconciliations infos
            #====================================================================================
            # Take the mapping lines
            elif "Mapping" in line and reconc_check == (1,1):
                # Build a mapping object
                mapping_object = c_mapping(line)
                # Add it to the list
                gs_mapping_list.append(mapping_object)
    
    return gene_tree, dict_module_mappingList, dict_module_modTree, gs_mapping_list

#====================================================================================
# Infers modules genes compositions
#====================================================================================

def infers_modulesCompo(gene_tree, dict_module_mappingList, dict_module_modTree) -> tuple:
    """
    Use the modules trees topollogy, the gene tree topology, 
    and their corresponding mappings from the reconciliation
    to infers the modules composition of each ancestral gene (internal nodes)
    Return a dict {gene : [module list]}
    
    Parameters
    ----------
    gene_tree : Tree
        Gene tree as a ete3.Tree object
    dict_module_mappingList : dict
        Dictionary with module name as key and list of mapping (c_mapping object) as value
    dict_module_modTree : dict
        Dictionary with module name as key and module tree (ete3.Tree) as value
        
    Returns
    -------
    dict_gene_modulesList : dict
        Dictionary with gene (str) as key and list of modules (c_module object) as value
    dict_module_mappingList : dict
        Dictionary with module name as key and list of mapping (c_mapping object) as value (completed)
    """
    # Init the dict ; gene : [modules list]
    dict_gene_modulesList = {n.name : [] for n in gene_tree.traverse()}
    # Attribute the colors and shapes to the modules
    dict_Module_ShapeColor = make_dict_Module_ShapeColor(dict_module_modTree.keys())
    # Doing it module per module
    for module, modTree in dict_module_modTree.items():
        # Use only "parts" of the big dict
        mapping_list = dict_module_mappingList[module]
        sub_dic_module_modTree = {m : mt for m, mt in dict_module_modTree.items() if m == module}
        # Get the repartition of the actual module
        dict_gene_moduleCount, gene_node_list = module_gene_inference(module, modTree, gene_tree, mapping_list, sub_dic_module_modTree)
        # Complete the inferred mapping dict, we want to keep the infferred nodes mapping informations
        for inferred_mapping in gene_node_list:
            dict_module_mappingList[module].append(c_mapping(f"{module}: Inferred, Mapping -> {inferred_mapping.name} inferred by logics topology tree"))
        mapping_list = dict_module_mappingList[module]
        # Use the repartition of the actual module, to complete the gene : module list dict
        for gene, count in dict_gene_moduleCount.items():
            # Use the gene node object
            gene_node = gene_tree&gene
            # Get the corresponding mapping
            module_mapping = [mapping for mapping in mapping_list if mapping.get_geneNode(gene_tree) == gene_node][0]
            # General informations
            module_string = module_mapping.entity_1
            split_module_string = module_string.split("_")[0].split("|")
            module_node_name = module_mapping.entity_1
            # Leaf case
            if gene_node.is_leaf():
                module_start = split_module_string[1]
                module_end = split_module_string[2]
            else:
                module_start = None
                module_end = None
            # Build module object 
            module_object = c_module(module, module_start, module_end, gene, module_node_name)
            module_object.get_shapeColor(dict_Module_ShapeColor)
            # Uniq list / presence of modules (could have more, dont consider them for now)
            if module_object.name not in [m.name for m in dict_gene_modulesList[gene]]:
                dict_gene_modulesList[gene].append(module_object)
    # Return the dict ;  gene (str) : modules (module ojbect) list (list) / And the completed module mappingList dict (with inferred mapping)
    return dict_gene_modulesList, dict_module_mappingList

def module_gene_inference(module, module_tree, gene_tree, mapping_list, dic_module_modTree) -> tuple:
    """
    Use the module tree topology, the gene tree topology, 
    and their mapping gene <-> module
    to infers the presence of the modules to the ancestral genes
    Return a dict : {gene : module_count}
    
    Parameters
    ----------
    module : str
        Name of the module
    module_tree : Tree
        The module tree as a ete3.Tree object
    gene_tree : Tree
        Gene tree as a ete3.Tree object
    mapping_list : list
        A list of mapping (c_mapping object)
    dict_module_modTree : dict
        Dictionary with module name as key and module tree (ete3.Tree) as value
        
    Returns
    -------
    dict_gene_moduleCount : dict
        Dictionary with gene as key (str) and module count as value (int)
    infers_node_list : list
        A list of the gene nodes where the presence of the module have been infered
    """
    dict_gene_moduleCount = {}
    infers_node_list = []
    # Iterate on the module tree, from the root (preorder)
    for mod_node in module_tree.iter_descendants("preorder"):
        # Get the corrsponding mapping (1 mapping per module node)
        actualNode_mapping = [mapping for mapping in mapping_list if mapping.get_modNode(dic_module_modTree) == mod_node]
        if len(actualNode_mapping) > 0:
            actualNode_mapping = actualNode_mapping[0]
            # Get the corrsponding gene node
            gene_node = actualNode_mapping.get_geneNode(gene_tree)
            # Gene node list for which we have an module presence / number info
            gene_node_list = []
            # Leaf event
            if actualNode_mapping.event == "Leaf":
                # Add it to the actual gene
                gene_node_list.append(gene_node)
            # Co-divergence event
            elif actualNode_mapping.event == "Co-divergence":
                # Add it the actual gene
                gene_node_list.append(gene_node)
                # Add it the gene nodes beetween this modules events, in his child (=the whole descendance)
                infers_nodes = descendance_mapping(mod_node, module_tree, gene_tree, mapping_list, dic_module_modTree)
                gene_node_list.extend(infers_nodes)
                infers_node_list.extend(infers_nodes)
            # Duplication event
            elif actualNode_mapping.event == "Domain duplication":
                # Add it the actual gene (2 times)
                # In duplication, a domains, and his chil, are mapped to the same gene
                # WARNING : that means we cant and dont have to make a descendance mapping
                gene_node_list.append(gene_node)
                gene_node_list.append(gene_node)
                # Add it the gene nodes beetween this modules events, in his child (=the whole descendance)
                #gene_node_list.extend(descendance_mapping(mod_node, module_tree, gene_tree, mapping_list, dic_module_modTree))
            # If the corresponding mapping event is an intra gene transfer
            elif actualNode_mapping.event == "Intra-gene-tree domain transfer":
                gene_node_list.append(gene_node)
                recipient_gene_node = get_any_geneNode(actualNode_mapping.recipient, gene_tree)
                # Add it to the genes nodes beetween the recipient gene and his most ancien descendant with module informations
                #infers_nodes = gene_descendance(recipient_gene_node, gene_tree, module_tree, mapping_list, dic_module_modTree)
                # A transfer is a basic module tree separation in 2 : co-divergence and transfer
                # The first is the normal evolution to the gene descendant (co-divergence)
                # The second is the evolution from the recipient to his descendants (teleportation in gene tree then gene evolution)
                # !!! not really this : we want from child1 to gene_module_desc and child2 to gene_module desc, not mod_gene to gene_module_desc !!
                infers_nodes = transfer_descendance_mapping(mod_node, recipient_gene_node, module_tree, gene_tree, mapping_list, dic_module_modTree)
                gene_node_list.extend(infers_nodes)
                infers_node_list.extend(infers_nodes)
                if recipient_gene_node not in infers_node_list:
                    infers_node_list.append(recipient_gene_node)
                    gene_node_list.append(recipient_gene_node)
            # Add dict entry for all genes in the node list
            for g_node in gene_node_list:
                # Infers the presence of the module to this gene (= add it to the dict)
                if g_node.name not in dict_gene_moduleCount.keys():
                    dict_gene_moduleCount[g_node.name] = 1
                # Or add an entry 
                else:
                    dict_gene_moduleCount[g_node.name] += 1
        else: print(f"No mapping found for {mod_node.name} !")
    # Return the dict ; gene (str) : module count (int), And the gene nodes list (the nodes we infers module presence)
    return dict_gene_moduleCount, infers_node_list

def descendance_mapping(module_node, module_tree, gene_tree, mapping_list, dic_module_modTree) -> list:
    """
    Check if the module child of a module, is in a gene that is the direct child of the gene 
    Module 1 <-> Gene A
    Module 2 (child of M1 in module tree) <-> Gene B
    It's possible that Gene B is not the direct child of Gene,
    -> it could be a more distant descendants !
    If it's the case, we want the path beetween Gene A and Gene B in the gene tree
    => They are the descendants of the gene with the last modules events, and the ancestors of the modules childs of the next event
    (2 childs, so 2 paths, we return the list of all nodes in theses)
    
    Parameters
    ----------
    module_node : Tree
        Node of the module as a ete3.Tree object
    module_tree : Tree
        The module tree as a ete3.Tree object
    gene_tree : Tree
        Gene tree as a ete3.Tree object
    mapping_list : list
        A list of mapping (c_mapping object)
    dict_module_modTree : dict
        Dictionary with module name as key and module tree (ete3.Tree) as value
        
    Returns
    -------
    geneNode_list : list
        List of gene node (as ete3.Tree object)
    """
    geneNode_list = []
    # We first want the mapping of our module of interest, so we know his corrsponding gene node
    interest_module_node_mapping = [mapping for mapping in mapping_list if mapping.get_modNode(dic_module_modTree) == module_node][0]
    interest_gene_node = interest_module_node_mapping.get_geneNode(gene_tree)
    # We look at the modules childs of the module of interest (so in module tree)
    for module_child in module_node.children:
        # We want the mapping of the module children, to know his corresponding gene (gene node)
        child_module_node_mapping = [mapping for mapping in mapping_list if mapping.get_modNode(dic_module_modTree) == module_child]
        if len(child_module_node_mapping) > 0:
            child_module_node_mapping = child_module_node_mapping[0]
            child_gene_node = child_module_node_mapping.get_geneNode(gene_tree)
            # Get the path (in gene tree), beetween the gene_node of our module of interest, and the gene_node of his children we"re looking at
            path_gene_child_parent = ete3_nodesPath(gene_tree, interest_gene_node.name, child_gene_node.name)
            geneNode_list.extend(path_gene_child_parent)
        else: print(f"No mapping found for {module_child.name} !")
    # Return the gene node (ete3 node object) list
    return geneNode_list

def transfer_descendance_mapping(module_node, recipient_gene_node, module_tree, gene_tree, mapping_list, dic_module_modTree) -> list:
    """
    Children of a module that mapp a transfer intra gene, are 1) the recipient and 2) his normal evolution
    Module I <-> Gene I
    Child of I in module tree are : 1) recipient 2) gene descendant
    It's possible that next module event is not the direct child of the recipient or the common gene evol,
    -> it could be a more distant descendants !
    If it's the case, we want the path beetween Recipient and Next in the gene tree, and path beetween normal child and Next descendant with event
    => They are the descendants of the gene with the last modules events, and the ancestors of the modules childs of the next event
    (2 childs, so 2 paths, we return the list of all nodes in theses)
    
    Parameters
    ----------
    module_node : Tree
        Node of the module as a ete3.Tree object
    recipient_gene_node : Tree
        Node of the recipient gene as a ete3.Tree object
    module_tree : Tree
        The module tree as a ete3.Tree object
    gene_tree : Tree
        Gene tree as a ete3.Tree object
    mapping_list : list
        A list of mapping (c_mapping object)
    dict_module_modTree : dict
        Dictionary with module name as key and module tree (ete3.Tree) as value
        
    Returns
    -------
    geneNode_list : list
        List of gene node (as ete3.Tree object)
    """
    geneNode_list = []
    # We first want the mapping of our module of interest, so we know his corrsponding gene node
    interest_module_node_mapping = [mapping for mapping in mapping_list if mapping.get_modNode(dic_module_modTree) == module_node][0]
    interest_gene_node = interest_module_node_mapping.get_geneNode(gene_tree)
    # We look at the modules childs of the module of interest (so in module tree)
    for module_child in module_node.children:
        # We want the mapping of the module children, to know his corresponding gene (gene node)
        child_module_node_mapping = [mapping for mapping in mapping_list if mapping.get_modNode(dic_module_modTree) == module_child]
        if len(child_module_node_mapping) > 0:
            child_module_node_mapping = child_module_node_mapping[0]
            child_gene_node = child_module_node_mapping.get_geneNode(gene_tree)
            # The first is the evolution from the recipient to his descendants (teleportation in gene tree then gene evolution)
            # Recipient is an ancestor of child gene_node
            if child_gene_node in recipient_gene_node:
                # Get the path (in gene tree), beetween the recipient gene_node of our module of interest, and the gene_node of his children we"re looking at
                path_gene_child = ete3_nodesPath(gene_tree, recipient_gene_node.name, child_gene_node.name)
                geneNode_list.extend(path_gene_child)
            # If the recipient is the child himself
            elif child_gene_node.name == recipient_gene_node.name:
                # If it's a leaf, don't mind, leaf mappin will add it
                if child_gene_node.is_leaf():
                    continue
                else:
                    # Else, it's an internal node so we add it
                    geneNode_list.append(child_gene_node)
            # The first is the normal evolution to the gene descendant (co-divergence)
            else:
                # Get the path (in gene tree), beetween the gene_node of our module of interest, and the gene_node of his children we"re looking at
                if interest_gene_node.name == child_gene_node.name:
                    geneNode_list.append(child_gene_node)
                else:
                    path_gene_child_parent = ete3_nodesPath(gene_tree, interest_gene_node.name, child_gene_node.name)
                    geneNode_list.extend(path_gene_child_parent)
        else: print(f"No mapping found for {module_child.name} !")
    # Return the gene node (ete3 node object) list
    return geneNode_list

def gene_descendance(gene_node, gene_tree, module_tree, mapping_list, dic_module_modTree) -> list:
    """
    Get the nodes list of all the descendants a gene node
    Beetween the gene node, and his most ancestral descendant with module information (=closest)
    
    Parameters
    ----------
    gene_node : Tree
        Node of the gene as a ete3.Tree object
    gene_tree : Tree
        Gene tree as a ete3.Tree object
    module_tree : Tree
        The module tree as a ete3.Tree object
    mapping_list : list
        A list of mapping (c_mapping object)
    dict_module_modTree : dict
        Dictionary with module name as key and module tree (ete3.Tree) as value
        
    Returns
    -------
    geneNode_list : list
        List of gene node (as ete3.Tree object)
    """
    geneNode_list = []
    # Make an full module node list (= all module nodes)
    module_node_list = [node for node in module_tree.traverse()]
    # Mapped them to their corresponding gene node (= gene nodes list of gene with modules informations)
    module_node_mapped_list = [mapping for mapping in mapping_list if mapping.get_modNode(dic_module_modTree) in module_node_list]
    g_node_module_mapped_list = [mapping.get_geneNode(gene_tree) for mapping in module_node_mapped_list]
    # Iterate on the genes nodes descendants of our gene node of interest (per level)
    for g_node in gene_node.iter_descendants("levelorder"):
        # If this gene node have module informations (on the mapped gene nodes module list)
        if g_node in g_node_module_mapped_list:
            # If this node is a descendant of one that we already got (so not the most ancient), dont mind it
            if g_node in geneNode_list:
                continue
            # If it's the more ancient descendants with module events informations
            else:
                geneNode_list.extend(ete3_nodesPath(gene_tree, gene_node.name, g_node.name))
    return geneNode_list

#====================================================================================
# Some complementes / tools functions
#====================================================================================

def get_any_geneNode(gene_name, gene_tree) -> object:
    """
    Search the gene node that corresponds to the gene name
    
    Parameters
    ----------
    gene_name : str
        Name of a gene
    gene_tree : Tree
        Gene tree as a ete3.Tree object
        
    Returns
    -------
    node : Tree
        A gene node as ete3.Tree object
    """
    for node in gene_tree.traverse():
        if gene_name in node.name:
            return node

def ete3_nodesPath(tree, nodeAnc_name, nodeChild_name) -> list:
    """
    Get the nodes between 2 nodes, one child (recent), and on ancestral
    using the topology oe their tree
    
    Parameters
    ----------
    tree : Tree
        A tree as a ete3.Tree object
    nodeAnc_name : str
        Name corresponding a node of the tree
    nodeChild_name : str
        Name corresponding a node of the tree
        
    Returns
    -------
    path : list
        A list of nodes (as ete3.Tree object)
    """
    # Get the node corresponding the name
    child_node = tree&nodeChild_name
    anc_node = tree&nodeAnc_name
    # Init the path list as empty
    path = []
    # Iterate from the child
    node = child_node
    # While his parent is not the ancestral
    while node.up != anc_node:
        # The iterator node become the parent of the previus / init
        node = node.up
        # Add the actual node to the path
        path.append(node)
    return path

def complemente_dict(dict_k_vlist, key, value) -> dict:
    """
    If not key, add a it and init the value list
    else add the value to the value list
    
    Parameters
    ----------
    dict_k_vlist : dict
        A dictionary with k as key and vlist as value (a list)
    key
        A key of the `dict_k_vlist` dictionary
    value
        A value to append to the value list associated to the `key`
        
    Returns
    -------
    dict_k_vlist : dict
        A dictionary with k as key and vlist as value (a list), (completed) 
    """
    if key not in dict_k_vlist:
        dict_k_vlist[key] = [value]
    else:
        dict_k_vlist[key].append(value)
    return dict_k_vlist

def complemente_dict_nonAdditif(dict_k_vlist, key, value) -> dict:
    """
    If not key, add a it and init the value list
    else add the value to the value list
    Non additif version : if its a domain to add, and its already present in the list
    we dont add it (to avoid copy)
    In theory, not needed cause we trust our previous choices
    
    Parameters
    ----------
    dict_k_vlist : dict
        A dictionary with k as key and vlist as value (a list)
    key
        A key of the `dict_k_vlist` dictionary
    value
        A value to append to the value list associated to the `key`
        
    Returns
    -------
    dict_k_vlist : dict
        A dictionary with k as key and vlist as value (a list), (completed) 
    """
    if key not in dict_k_vlist:
        dict_k_vlist[key] = [value]
    else:
        # In case we wand to add a c_module object, is there is already one with same name
        # Dont add it, cause its the same module entity
        # Just another instance of the same module
        if type(value) is c_module:
            if value.name not in [i.name for i in dict_k_vlist[key]]:
                dict_k_vlist[key].append(value)
        else:
            dict_k_vlist[key].append(value)
    return dict_k_vlist

#==============================================================================
# Modules composition change 
#==============================================================================

def make_dict_module_change(dict_gene_domainList, gene_tree) -> dict:
    """
    Use the tree topology, and the nodes modules composition to make a dict:
    { gene : [ [news modules] , [lost modules] ] }
    A composition change, is the differences in module of a node, with his parent
    (so can have news or lost modules that his parent had)
    
    Parameters
    ----------
    dict_gene_domainList : dict
        A dictionary with gene name (str) as key and list of domains/modules as value (being a c_module object)
    gene_tree : Tree
        A tree as ete3.Tree object, including gene nodes with names of the  genes in `dict_gene_domainList`
        
    Returns
    -------
    dict_gene_moduleChange : dict
        A dictionary with gene name (str) as key and list of modules changes [list of modules gained and list of modules lost] as value 
    """
    # Init the dict
    dict_gene_moduleChange = {}
    # Iterate on the gene_tree
    for node in gene_tree.traverse("postorder"):
        # We look at all (leafs + internals)
        # Get the module list of the protein at the node
        interest_modules_list = [module.name for module in dict_gene_domainList[node.name]]
        # If this node is the root, he have no parent, its the first
        if node.is_root():
            parent_modules_list = []
        else:
            # Get the parent of this proteins (to compare the 2 modules list)
            parent_node = (gene_tree&node.name).up
            # Get the module list of that protein at parent node
            if parent_node.is_root() == False:
                parent_modules_list = [module.name for module in dict_gene_domainList[parent_node.name]]
            else:
                parent_modules_list= []
        # Init the count
        news_module, lost_module = [], []
        # Count the gain (in interest, but not in parent)
        for module in interest_modules_list:
            if module not in parent_modules_list: news_module.append(module)
        # Count the lost (in parent, not in interest)
        for module in parent_modules_list:
            if module not in interest_modules_list: lost_module.append(module)
        # Add it to the dict
        dict_gene_moduleChange[node.name] = [news_module, lost_module]
    # Return the dict
    return dict_gene_moduleChange

#==============================================================================
# Look the positions of modules on their leafs, for each nodes
#==============================================================================
 
def look_at_modules_on_leafs(dict_gene_moduleChange, dict_gene_moduleList, gene_tree) -> dict:
    """
    Here, we look at the modules change, and specifically at how are thses modules on the leafs of their changes
    Return a dict
    {gene : { dict leafs : modules of the module change }}}
    
    Parameters
    ----------
    dict_gene_moduleChange : dict
        A dictionary with gene name (str) as key and list of modules changes [list of modules gained and list of modules lost] as value 
    dict_gene_moduleList : dict
        A dictionary with gene name (str) as key and list of modules as value (being a c_module object)
    gene_tree : Tree
        A tree as ete3.Tree object
        
    Returns
    -------
    dict_gene_LeafsModulesThatChanged : dict
        A dictionary with gene name (str) as key and a dict of modules changes at leaves as value 
    """
    # Init the out dict
    dict_gene_LeafsModulesThatChanged = {}
    # For each change in the history (= a gene node)
    for gene, moduleChange in dict_gene_moduleChange.items():
        # Init the leaf dict for this gene
        leaf_dict = {}
        # Get the changes
        win_modules, lost_modules = moduleChange[0], moduleChange[1]
        # Get his leaves descendants (using tree topology)
        # For all leaves (to see lost)
        for leaf in gene_tree:
            # Get the leaf name, its an actual protein (with a sequence and all)
            actual_gene = leaf.name
            # Get his module decomposition
            actual_gene_modulesList = dict_gene_moduleList[actual_gene]
            # Search for the win modules at the nodes, in this actual protein modules list
            win_modules_corresponding = [module for module in actual_gene_modulesList if module.name in win_modules]
            # Search for the lost modules at the nodes, in this actual protien modules list
            lost_modules_corresponding = [module for module in actual_gene_modulesList if module.name in lost_modules]
            # Add it to the leaf dict
            # Descendant case
            if leaf in gene_tree&gene:
                leaf_dict[f"desc_{leaf.name}"] = [win_modules_corresponding, lost_modules_corresponding, actual_gene_modulesList]
            # Non descendant case
            else:
                leaf_dict[f"nonDesc_{leaf.name}"] = [win_modules_corresponding, lost_modules_corresponding, actual_gene_modulesList]
        # Add it to the gene : leaf change dict
        dict_gene_LeafsModulesThatChanged[gene] = leaf_dict
    # Return the dict
    return dict_gene_LeafsModulesThatChanged

def modules_on_leafs_output(dict_gene_LeafsModulesThatChanged, dict_gene_moduleList, dict_gene_domainsList, directory) -> None:
    """
    Make an output directory to have an idea of modules implicated on modules changes
    Module at an internal nodes are not consistent, so wee look at their correspondance on the real / true gene
    ( a file per gene / node of change)
    
    Parameters
    ----------
    dict_gene_LeafsModulesThatChanged : dict
        A dictionary with gene name (str) as key and a dict of modules changes at leaves as value 
    dict_gene_moduleList : dict
        A dictionary with gene name (str) as key and list of modules as value (being a c_module object)
    dict_gene_domainsList : dict
        A dictionary with gene name (str) as key and list of domains as value (being a c_module object)
    directory : str
        Name of the output directory
    """
    # Make the dict if its not existent
    if not os.path.exists(directory):
        os.makedirs(directory)
    # For each gene / node of modules change
    executor = concurrent.futures.ProcessPoolExecutor()
    # pdf version
    #futures = [executor.submit(gene_pdf, gene, leaves_change, dict_gene_domainsList, directory) for gene, leaves_change in dict_gene_LeafsModulesThatChanged.items()]
    # Itol domains version
    futures = [executor.submit(gene_itol, gene, leaves_change, dict_gene_domainsList, dict_gene_moduleList, directory) for gene, leaves_change in dict_gene_LeafsModulesThatChanged.items()]
    concurrent.futures.wait(futures)

def gene_itol(gene, leaves_change, dict_gene_domainsList, dict_gene_PresentmoduleList, directory) -> None:
    """
    Create itol modules changes / domains known
    
    Parameters
    ----------
    gene : str
        Name of a gene
    leaves_change : dict
        Dictionary with leaf node name (str) as key and list of module change [list of modules gained, list of modules lost] as value
    dict_gene_domainsList : dict
        A dictionary with gene name (str) as key and list of domains as value (being a c_module object)
    dict_gene_PresentmoduleList : dict
        A dictionary with gene name (str) as key and list of modules as value (being a c_module object)
    directory : str
        Name of the output directory
    """
    # Dict with only module gained / lost for each leaf
    dict_gene_moduleList = {}
    # For each leaf that is a descendant of this gene
    for leaf, modules_change  in leaves_change.items():
        leaf_name = "_".join(leaf.split("_")[1:])
        dict_gene_moduleList[leaf_name] = [[],[]]
        # For each modules (win)
        for mod in modules_change[0]:
            # Add the module gained to the list of module gained (for itol file generation)
            dict_gene_moduleList[leaf_name][0].append(mod)
        # For each modules (lost)
        for mod in modules_change[1]:
             # Add the module lost to the list of module lost (for itol file generation)
            dict_gene_moduleList[leaf_name][1].append(mod)
    # Itol output domains file
    write_itol_domain_gained_lost(f"{directory}/itolModulesThatChanged_{gene}_only_mod.txt", dict_gene_domainsList, dict_gene_moduleList, gene)
    write_itol_present(f"{directory}/itolModulesPresent_{gene}_only_mod.txt", dict_gene_domainsList, dict_gene_PresentmoduleList, gene)
    #write_itol_domain_and_modules_gained_lost(f"{directory}/itolModulesThatChanged_{gene}_dom_mod.txt", dict_gene_domainsList, dict_gene_moduleList, gene)
       
def gene_pdf(gene, leaves_change, dict_gene_domainsList, directory) -> None:
    """
    Create pdf(s) for 1 gene (so multiple pdf, 1 per leaves)
    
    Parameters
    ----------
    gene : str
        Name of a gene
    leaves_change : dict
        Dictionary with leaf node name (str) as key and list of module change [list of modules gained, list of modules lost] as value
    dict_gene_domainsList : dict
        A dictionary with gene name (str) as key and list of domains as value (being a c_module object)
    directory : str
        Name of the output directory
    """
    # Dict with only module gained / lost for each leaf
    dict_gene_moduleList = {}
    # For each leaf that is a descendant of this gene
    for leaf, modules_change  in leaves_change.items():
        # Init visu variables
        features = []
        min_pos, max_pos = 9999999999999999, 0
        # For all modules (paloma found)
        """
        for mod in modules_change[2]:
            features.append(GraphicFeature(start=int(mod.start), end=int(mod.end), strand=+1, color="#DCDCDC", label=mod.name))
            if int(mod.end) > max_pos: max_pos = int(mod.end)
            if int(mod.start) < min_pos: min_pos = int(mod.start)
        """
        leaf_name = "_".join(leaf.split("_")[1:])
        dict_gene_moduleList[leaf_name] = [[],[]]
        # For each modules (win)
        for mod in modules_change[0]:
            features.append(GraphicFeature(start=int(mod.start), end=int(mod.end), strand=+1, color="#32CD32", label=mod.name))
            # Add the module gained to the list of module gained (for itol file generation)
            dict_gene_moduleList[leaf_name][0].append(mod)
        # For each modules (lost)
        for mod in modules_change[1]:
            features.append(GraphicFeature(start=int(mod.start), end=int(mod.end), strand=+1, color="#B22222", label=mod.name))
             # Add the module lost to the list of module lost (for itol file generation)
            dict_gene_moduleList[leaf_name][1].append(mod)
        # For known domains
        for domain in dict_gene_domainsList[leaf_name]:
            features.append(GraphicFeature(start=int(domain.start), end=int(domain.end), strand=+1, color="#00BFFF", label=domain.name))
            if int(domain.end) > max_pos: max_pos = int(domain.end)
            if int(domain.start) < min_pos: min_pos = int(domain.start)
        # Build record and build figure
        record = GraphicRecord(sequence_length=max_pos, features=features)
        cropped_record = record.crop((min_pos, max_pos))
        ax, _ = cropped_record.plot(figure_width=20)
        # Save figure on file
        ax.figure.savefig(f"{directory}/modulesThatChanged_{gene}_{leaf}.pdf", bbox_inches="tight")
        plt.pyplot.close(ax.figure)
        
    # Itol output domains file
    #write_itol_domain_gained_lost(f"{directory}/itolModulesThatChanged_{gene}_all.txt", dict_gene_domainsList, dict_gene_moduleList, gene)
    write_itol_domain_and_modules_gained_lost(f"{directory}/itolModulesThatChanged_{gene}_dom_mod.txt", dict_gene_domainsList, dict_gene_moduleList, gene)

#==============================================================================
# Generation of "simpliest" things for itol files generation
#==============================================================================

def make_module_gene_recipient_list(dict_gene_moduleObjList, dict_module_mappingList) -> list:
    """
    Make a list : [module_object, gene, recipient]
    
    Parameters
    ----------
    dict_gene_moduleObjList : dict
        A dictionary with gene name (str) as key and list of modules as value (being a c_module object)
    dict_module_mappingList : dict
        Dictionary with module name as key and list of mapping (c_mapping object) as value
        
    Returns
    -------
    list_module_gene_recipient : list
        A list of list [c_module object, gene, recipient]
    """
    list_module_gene_recipient = []
    # Iterate on mapping events, to find transfer
    for module, mappingList in dict_module_mappingList.items():
        for mapping in mappingList:
            if mapping.event == "Intra-gene-tree domain transfer":
                # Take all info of the mapping transfer object
                gene = mapping.entity_2
                recipient = mapping.recipient
                # Search for the corresponding module object
                for o_mod in dict_gene_moduleObjList[gene]:
                    if o_mod.name == module:
                        o_module = o_mod
                        break
                list_module_gene_recipient.append([o_module, gene, recipient])
    return list_module_gene_recipient

#==============================================================================
# Load complementary annotations
#==============================================================================

def load_pastml_ancestral_states(pastml_tab, gene_tree) -> dict:
    """
    Read a pastML combined_ancestral_states.tab file
    https://github.com/evolbioinfo/pastml
    get a header, and first column is the node name, all the other are presence or absence of a carac
    (ex : GOA termes)
    Return a dict :
    { gene : [GO list] }
    
    Parameters
    ----------
    pastml_tab : str
        A name of a tab file , e.g., pastML combined_ancestral_states.tab file
    gene_tree : Tree
        A tree as ete3.Tree object
        
    Returns
    -------
    dict_nodeName_annotationsList : dict
        Dictionary with gene node name (str) as key and list of annotations (str) as value
    """
    # Init out dict
    dict_nodeName_annotationsList = {}
    # Init a header list (to know the order of the carac / columns)
    header_list = []
    # Open and read the tab formated file
    with open(pastml_tab, "r") as pastml_file:
        # Iterate on his line
        for line in pastml_file:
            # Split it
            splited_line = line.replace("\n","").split("\t")
            # If its the first, line, make the header list
            if len(header_list) == 0:
                header_list = splited_line
            # Else, its a gene entry, so we use it to complete our dict
            else:
                gene_name = ""
                # Check for name correspondance in out tree
                for gene in gene_tree.traverse("postorder"):
                    if splited_line[0] in gene.name:
                        gene_name = gene.name
                        annotations_list = [header_list[indice] for indice in range(len(splited_line)) if splited_line[indice] == "1"]
                        # Agglomerate all possible presence of the 2 scenario (if present 1 scenario, considers it as present)
                        if gene_name in dict_nodeName_annotationsList:
                            old_annotations_list = dict_nodeName_annotationsList[gene_name]
                            complete_annotationsList = []
                            for annot in header_list:
                                if annot in [*annotations_list, *old_annotations_list]:
                                    complete_annotationsList.append(annot)
                            # TEST 
                            #dict_nodeName_annotationsList[gene_name] = complete_annotationsList
                        else:
                            dict_nodeName_annotationsList[gene_name] = annotations_list
    # Add as empty list the node without infos
    for gene in gene_tree.traverse("postorder"):
        if gene.name not in dict_nodeName_annotationsList:
            dict_nodeName_annotationsList[gene.name] = []
    return dict_nodeName_annotationsList

def make_dict_annotation_change(dict_nodeName_annotationsList, gene_tree) -> dict:
    """
    Use the tree topology, and the nodes annotation to make a dict:
    { gene : [ [news annotation] , [lost annotation] ] }
    A composition change, is the differences in annotation of a node, with his parent
    (so can have news or lost annotations that his parent had)
    
    Parameters
    ----------
    dict_nodeName_annotationsList : dict
        Dictionary with gene node name (str) as key and list of annotations (str) as value
    gene_tree : Tree
        A tree as ete3.Tree object
        
    Returns
    -------
    dict_gene_annotationsChange : dict
        Dictionary with gene node name (str) as key and list of changes (list, [list of gained annotations, list of lost annotations]) as value
    """
    # Init the dict
    dict_gene_annotationsChange = {}
    # Iterate on the gene_tree
    for node in gene_tree.traverse("postorder"):
        # We look at all (leafs + internals)
        if node.name in dict_nodeName_annotationsList:
            # Get the annotation list of the protein at the node
            interest_annotations_list = dict_nodeName_annotationsList[node.name]
            # If this node is the root, its the first (no parent)
            if node.is_root():
                parent_annotations_list = []
            else:
                # Get the parent of this proteins (to compare the 2 annotations list)
                parent_node = (gene_tree&node.name).up
                parent_annotations_list = []
                # Get the annotation list of that protein at parent node
                if parent_node.is_root() == False and parent_node.name in dict_nodeName_annotationsList:
                    parent_annotations_list = dict_nodeName_annotationsList[parent_node.name]
                else:
                    parent_annotations_list = []
            # Init the count
            news_annotations, lost_annotations = [], []
            # Count the gain (in interest, but not in parent)
            for annot in interest_annotations_list:
                if annot not in parent_annotations_list: news_annotations.append(annot)
            # Count the lost (in parent, not in interest)
            for annot in parent_annotations_list:
                if annot not in interest_annotations_list: lost_annotations.append(annot)
            # Add it to the dict
            dict_gene_annotationsChange[node.name] = [news_annotations, lost_annotations]
    # Return the dict
    return dict_gene_annotationsChange

def write_function_module_change(dict_gene_moduleChange, dict_nodeName_annotationsChange, filename) -> None:
    """
    Write a csv file with nodes with annotations / functions change, 
    and the module composition change at these gene
    
    Parameters
    ----------
    dict_gene_moduleChange : dict
        A dictionary with gene name (str) as key and list of modules changes [list of modules gained and list of modules lost] as value 
    dict_nodeName_annotationsChange : dict
        Dictionary with gene node name (str) as key and list of changes (list, [list of gained annotations, list of lost annotations]) as value
    filename : str
        Name of the output file
    """
    # Make the file 
    with open(filename, "w+") as csv_file:
        # Header
        csv_file.write("gene,function(s)_win,function(s)_lost,module(s)_win,module(s)_lost\n")
        # For each gene
        for gene, annotations_change in dict_nodeName_annotationsChange.items():
            # Write a line (if there is annotations change)
            if len(annotations_change[0]) + len(annotations_change[1]) > 0:
                csv_file.write(f"{gene},{'|'.join(annotations_change[0])},{'|'.join(annotations_change[1])},{'|'.join(dict_gene_moduleChange[gene][0])},{'|'.join(dict_gene_moduleChange[gene][1])}\n")
                
def write_complete_function_module_change(dict_gene_moduleChange, dict_nodeName_annotationsChange, filename) -> None:
    """
    Write a csv file with nodes with all annotations / functions change, 
    and the module composition change at these gene
    (even if no annotation change)
    
    Parameters
    ----------
    dict_gene_moduleChange : dict
        A dictionary with gene name (str) as key and list of modules changes [list of modules gained and list of modules lost] as value 
    dict_nodeName_annotationsChange : dict
        Dictionary with gene node name (str) as key and list of changes (list, [list of gained annotations, list of lost annotations]) as value
    filename : str
        Name of the output file
    """
    # Make the file 
    with open(filename, "w+") as csv_file:
        # Header
        csv_file.write("gene,function(s)_win,function(s)_lost,module(s)_win,module(s)_lost\n")
        # For each gene
        for gene, annotations_change in dict_nodeName_annotationsChange.items():
            # Write all line
            csv_file.write(f"{gene},{'|'.join(annotations_change[0])},{'|'.join(annotations_change[1])},{'|'.join(dict_gene_moduleChange[gene][0])},{'|'.join(dict_gene_moduleChange[gene][1])}\n")

def write_function_module_change_expanded(gene_tree, dict_gene_moduleList, dict_gene_moduleChange, dict_nodeName_annotationsChange, filename) -> None:
    """
    Write a csv file with nodes with annotations / functions changes, 
    but expanded to the perspectives of all their descendants (ie leaves)
    
    Parameters
    ----------
    gene_tree : Tree
        A tree as a ete3.Tree object
    dict_gene_moduleList : dict
        A dictionary with gene name (str) as key and list of modules as value (being a c_module object)
    dict_gene_moduleChange : dict
        A dictionary with gene name (str) as key and list of modules changes [list of modules gained and list of modules lost] as value 
    dict_nodeName_annotationsChange : dict
        Dictionary with gene node name (str) as key and list of changes (list, [list of gained annotations, list of lost annotations]) as value
    filename : str
        Name of the output file
    """
    # Make the file
    with open(filename, "w+") as csv_file:
        # Header
        csv_file.write("gene,leaf,function(s)_win,function(s)_lost,module(s)_win,module(s)_lost\n")
        # For each gene
        for gene, annotations_change in dict_nodeName_annotationsChange.items():
            # Get his descendants (ie leaves)
            for node in gene_tree.traverse():
                if node.name == gene:
                    descendants_leaves = node.get_leaves()
            if "_" in gene:
                # For each of the desceandants of this gene
                for leaf in descendants_leaves:
                    anc_modules_gain = dict_gene_moduleChange[gene][0] 
                    anc_modules_lost = dict_gene_moduleChange[gene][1]
                    leaf_modules = [mod.name for mod in dict_gene_moduleList[leaf.name]]
                    leaf_modules_gain = []
                    leaf_modules_lost = []
                    for module in leaf_modules:
                        if module in anc_modules_gain:
                            leaf_modules_gain.append(module)
                        if module in anc_modules_lost:
                            leaf_modules_lost.append(module)
                    # write a line (all cases !)
                    csv_file.write(f"{gene},{leaf.name},{'|'.join(annotations_change[0])},{'|'.join(annotations_change[1])},{'|'.join(leaf_modules_gain)},{'|'.join(leaf_modules_lost)}\n")

#==============================================================================
# Gene tree function
#==============================================================================

def full_gene_tree(reconc_gene_tree, length_gene_tree) -> object:
    """
    Fusion the informations of the gene trees
    Being ; 
    1) the internal nodes names of the reconcilied gene tree (from dgs)
    2) the branch length from the phylogenetic tree (before dgs)
    
    Parameters
    ----------
    reconc_gene_tree : Tree
        A tree as a ete3.Tree, with internal nodes names
    length_gene_tree : Tree
        A tree as ete3.Tree, with branch length
        
    Returns
    -------
    full_gene_tree : Tree
        A tree as ete3.Tree, with internal nodes names and branch length
    """
    for node in length_gene_tree.traverse():
        for reconc_node in reconc_gene_tree.traverse():
            if [n.name.split('_')[0] for n in node.get_leaves()] == [n.name.split('_')[0] for n in reconc_node.get_leaves()]:
                node.name = reconc_node.name
    full_gene_tree = length_gene_tree
    t = Tree()
    t.add_child(full_gene_tree)
    full_gene_tree = t
    return full_gene_tree

#==============================================================================
# Wrtie csv file for the species-gene event from their reconciliation
#==============================================================================

def write_sp_gene_event(dgs_reconciliation_output_fn, gene_tree_fn) -> str:
    """
    Write csv directly from the seadog output and reconc gene tree
    (usefull for launch pastml before integrate all)
    
    Parameters
    ----------
    dgs_reconciliation_output_fn : str
        Name of a DGS-reconciliation file in SEADOG-MD output format
    gene_tree_fn : str
        Name of a tree file in newick format
        
    Returns
    -------
    reconc_gene_tree : str
        Name of a tree file in newick format
    filename : str
       Name of a Gene-Species events file in csv format
    """
    gene_tree, dict_module_mappingList, dict_module_modTree, gs_mapping_list = read_seadogO(dgs_reconciliation_output_fn, gene_tree_fn)
    dict_gene_spGeneEvent = {gs_mapping.entity_1 : gs_mapping.event for gs_mapping in gs_mapping_list}
    filename = Path(f"{dgs_reconciliation_output_fn.parents[0]}/{dgs_reconciliation_output_fn.stem}_sp_gene_event.csv").resolve()
    write_spGeneEvent_csv(dict_gene_spGeneEvent, filename)
    reconc_gene_tree = Path(f"{dgs_reconciliation_output_fn.parents[0]}/{dgs_reconciliation_output_fn.stem}_gene.tree").resolve()
    gene_tree.write(format=1, outfile=reconc_gene_tree,  format_root_node=True)
    return reconc_gene_tree, filename

def write_spGeneEvent_csv(dict_gene_spGeneEvent, filename) -> None:
    """
    Write a csv file
    gene (= node name of the the gene tree), event gene-specie (Speciation or Gene duplication, or just a leaf)
    
    Parameters
    ----------
    dict_gene_spGeneEvent : dict
        Dictionary with gene name (str) as key and type of Gene-Species event (str) as value
    filename : str
        Name of the output file, csv format
    """
    # Create the csv file
    with open(filename, "w+") as csv_file:
        # For each entry of the dict
        for gene, sp_gene_event in dict_gene_spGeneEvent.items():
            # Write a line of the csv (1 dict entry = 1 csv line)
            csv_file.write(f"{gene},{sp_gene_event}\n")

#==============================================================================
# Write csv files, for each nodes, modules composition, and modules composition change
#==============================================================================

def write_module_compo_csv(dict_gene_moduleList, filename) -> None:
    """
    Write a csv file
    node_name, module composition list (split by |)
    
    Parameters
    ----------
    dict_gene_moduleList : dict
        A dictionary with gene name (str) as key and list of modules as value (being a c_module object)
    filename : str
        Name of the output file, csv format
    """
    # Create the csv file
    with open(filename, "w+") as csv_file:
        # Init header
        csv_file.write("id,modules_composition\n")
        # Each entry of the dict is a csv line
        for gene, module_list in dict_gene_moduleList.items():
            csv_file.write(f"{gene},{'|'.join([module.name for module in module_list])}\n")

def write_module_change_csv(dict_gene_moduleChange, filename) -> None:
    """
    Write a csv file
    node_name, modules gained, modules lost (modules spliter being |)
    Use the gene moduleChange dict
    
    Parameters
    ----------
    dict_gene_moduleChange : dict
        A dictionary with gene name (str) as key and list of modules changes [list of modules gained and list of modules lost] as value 
    filename : str
        Name of the output file, csv format
    """
    # Create the csv file
    with open(filename, "w+") as csv_file:
        # Init header
        csv_file.write("id,news_modules,lost_modules\n")
        # Iterate on the gene, moduleChange dict (contains all we need)
        for gene, moduleChange in dict_gene_moduleChange.items():
            # A moduleChange is a list with in 0 the news, in 1 the lost modules
            news_module = moduleChange[0]
            lost_module = moduleChange[1]
            # Write this entry in the csv file
            csv_file.write(f"{gene},{'|'.join(news_module)},{'|'.join(lost_module)}\n")

#==============================================================================
# Read domains csv file
#==============================================================================

def load_domains_from_csv(gene_list, domains_file) -> dict:
    """
    Load the domains coordinates from a csv file
    csv file is : geneName (no node ID), domainsName, start, stop
    Return a dict
    { gene : [ domain_dict list ] }
    A domain of the domain_dict list is :
    { domain_name : [ start, end ]}
    
    Parameters
    ----------
    gene_list : list
        A list of gene node name (str)
    domains_file : str
        Name of domains file in csv format (geneName (no node ID), domainsName, start, stop)
        
    Returns
    -------
    dict_gene_domainsList : dict
        Dictionary with gene node name (str) as key and list of the associated domains as value
    """
    # Init the dict, with full gene_list (with node id)
    dict_gene_domainsList = {gene : [] for gene in gene_list}
    # Open and read line per line the input csv_file
    with open(domains_file, "r") as csv_file:
        # Iterate on csv file
        for line in csv_file:
            # Split the line
            splited_line = line.replace("\n", "").split(",")
            g_name, d_name, d_start, d_end = splited_line[0], splited_line[1], splited_line[2], splited_line[3]
            # Get the gene name with node id, from the raw gene name of the csv file (nod node id)
            for gene in gene_list:
                if g_name in gene: full_geneName = gene
            # Build the domains dict (1 csv line = 1 domain = 1 dict)
            dict_domain_posList = {f"{d_name}|{d_start}|{d_end}" : [d_start, d_end]}
            # Add it to the gene domains list dict
            dict_gene_domainsList[full_geneName].append(dict_domain_posList)
    # Return the dict
    return dict_gene_domainsList

def domains_as_modules(dict_gene_domainsList) -> dict:
    """
    Build a dict gene : domains list, 
    with domain being module object (for modules operations ...)
    
    Parameters
    ----------
    dict_gene_domainsList : dict
        Dictionary with gene node name (str) as key and list of the associated domains as value (a domain being a dictionary)
        
    Returns
    -------
    dict_gene_domainsOList : dict
        Dictionary with gene node name (str) as key and list of the associated domains as value (a domain being a c_module object)
    """
    dict_gene_domainsOList = {}
    domain_list_uniq = []
    for gene, domains_list in dict_gene_domainsList.items():
        domains_o_list = []
        for domain in domains_list:
            for n, pos in domain.items():
                name = n.split("|")[0]
                start, end = pos[0], pos[1]
            domain_object = c_module(name, start, end, gene, name)
            domains_o_list.append(domain_object)
            if name not in domain_list_uniq:
                domain_list_uniq.append(name)
        dict_gene_domainsOList[gene] = domains_o_list
    dict_Module_ShapeColor = make_dict_Domains_ShapeColor(domain_list_uniq)
    for gene, domains_o_list in dict_gene_domainsOList.items():
        for d_o in domains_o_list:
            d_o.get_shapeColor(dict_Module_ShapeColor)
    return dict_gene_domainsOList

def dict_domains_in_modules(dict_gene_moduleList, dict_gene_domainsList) -> dict:
    """
    Get the modules subdivisions of all the domains
    Here, we want to know the modules that are present in the domains, using proteics sequences positions
    Return a dict :
    { gene : { domain_modulesList_dict } }
    
    Parameters
    ----------
    dict_gene_moduleList : dict
        Dictionary with gene node name (str) as key and list of the associated modules as value (a module being a c_module object)
    dict_gene_domainsList : dict
        Dictionary with gene node name (str) as key and list of the associated domains as value (a domain being a dictionary)
        
    Returns
    -------
    dict_gene_domainsModules : dict
        Dictionary with gene node name (str) as key and a dictionary of domains and modules in it as value (domain being a key, list of module as a value)
    """
    # Init dict
    dict_gene_domainsModules = {gene : {} for gene in dict_gene_domainsList.keys()}
    # Iterate on gene (we want to consider each gene separately)
    for gene, domains_list in dict_gene_domainsList.items():
        # For each domains of this gene
        for domain in domains_list:
            # Get the domains info
            for d_name, d_pos in domain.items():
                # Get the module subdivisions of this domain, using seq positions
                module_subdiv_list = domain_in_modules(d_pos, dict_gene_moduleList[gene])
            # Add this module to the dict
            dict_gene_domainsModules[gene][d_name] = module_subdiv_list
    # Return the dict
    return dict_gene_domainsModules

def infers_domains_from_modules(dict_gene_moduleList, dict_gene_domainsModules) -> dict:
    """
    Use the known domains / modules compositions
    (1 domains is a sum of modules) 
    To search for non indentified domains
    Our gene is a sum of modules, we want to scan it for the presence of non known
    sum of modules, that can represents a domains not found yet
    Particulary useful for ancestral genes, that have non domains
    
    Parameters
    ----------
    dict_gene_moduleList : dict
        Dictionary with gene node name (str) as key and list of the associated modules as value (a module being a c_module object)
    dict_gene_domainsModules : dict
        Dictionary with gene node name (str) as key and a dictionary of domains and modules in it as value (domain being a key, list of module as a value)
        
    Returns
    -------
    dict_gene_infDomainsModules : dict
        Dictionary with gene node name (str) as key and a dictionary of infered domains (based on modules compositions)
    """
    # Init the new dict for inferred domaind modules (empty, so full inferred domains)
    dict_gene_infDomainsModules = {gene : {} for gene in dict_gene_moduleList.keys()}
    # Begin by making a list with all the known domains (as modules list)
    all_dom = []
    for gene, domains in dict_gene_domainsModules.items():
        for d_name, d_modules in domains.items():
            all_dom.append([d_name, d_modules])
    # For each gene
    for gene, module_list in dict_gene_moduleList.items():
        # For each domains of the all domains list (all dom knowns)
        for domain in all_dom:
            # Init a list of the genes module of the actual domain found in this gene
            d_modules_found = []
            # Look 1 module at a time
            for d_module in domain[1]:
                # For the actual module, look at all the modules of the gene to compare
                for g_module in module_list:
                    # It the actual domains module is in the gene module list (compa only name)
                    if d_module.name.split("|")[0] == g_module.name.split("|")[0]:
                        # Add the module to the list of the modules of the actual domain that are found in this gene
                        d_modules_found.append(g_module)
            # If all the modules of the actual domains are found in this gene
            if len(d_modules_found) == len(domain[1]) and len(domain[1]) > 0:
                # This domains exist in this gene (by his sum of modules definition)
                new_domain_name, new_domain_modList = f"{domain[0].split('|')[0]}|0|0", d_modules_found
                # Add it the new dict of inferred domains (if not already in, using "true name")
                if new_domain_name.split("|")[0] not in [d.split("|")[0] for d in dict_gene_infDomainsModules[gene].keys()]:
                    dict_gene_infDomainsModules[gene][new_domain_name] = new_domain_modList
    # Return the infered domains list (using sum of modules domain definition)
    return dict_gene_infDomainsModules

def domain_in_modules(domain_pos_list, moduleList) -> list:
    """
    Get the module subdivisions of the given domain
    (Protein specific)
    a domains is a start / stop on the proteic sequence
    for each module, we got his start / stop of the same proteic sequence
    Here, we want to know the modules that are present in the domains, using proteics sequences positions
    Return a list of the modules present in this domain
    
    Parameters
    ----------
    domain_pos_list : list
        A list of start and stop position of a domain
    moduleList : list
        A list of modules, a module being a c_module object
            
    Returns
    -------
    d_moduleList : list
        A list of selected modules, a module being a c_module object
    """
    # Init the module list of this domain
    d_moduleList = []
    # Get positions variable
    dom_start, dom_end = domain_pos_list[0], domain_pos_list[1]
    # Iterate on the known modules of this gene (protein sequence)
    for mod in moduleList:
        # Look for overlapp mod / domain
        if int(mod.end) >= int(dom_start) and int(dom_end) >= int(mod.start):
            # Add it to the module list of the domain
            d_moduleList.append(mod)
    # Return this list
    return d_moduleList

def modules_not_in_domains(dict_gene_moduleList, dict_gene_domainsModules) -> dict:
    """
    For each genes, get the modules that are not in a domains
    (= free modules)
    Return them as a dictionnary
    { gene : [free modules list] }
    
    Parameters
    ----------
    dict_gene_moduleList : dict
        Dictionary with gene node name (str) as key and list of the associated modules as value (a module being a c_module object)
    dict_gene_domainsModules : dict
        Dictionary with gene node name (str) as key and a dictionary of domains and modules in it as value (domain being a key, list of module as a value)
            
    Returns
    -------
    dict_gene_freeModules : dict
        Dictionary with gene node name (str) as key and a list of modules as value (a module being a c_module object)
    """
    # Init the dict
    dict_gene_freeModules = {gene : [] for gene in dict_gene_moduleList.keys()}
    # Iterate on the genes and modules list of the genes
    for gene, modules_list in dict_gene_moduleList.items():
        # Extract all the values list into 1 concat list (= all modules into domains for this gene)
        all_mod_in_dom = [l for j in dict_gene_domainsModules[gene].values() for l in j]
        # Iterate on all modules of the gene
        for module in modules_list:
            # If the actuel module is not in a domains
            if module not in all_mod_in_dom:
                # He's free, so we add it to the free module list of the gene
                dict_gene_freeModules[gene].append(module)
    # Return the dict
    return dict_gene_freeModules

def ghost_domains(dict_gene_domainsModules) -> dict:
    """
    For each gene, get the ghost domains
    They are domain without any modules (so lost / false informations)
    Return a dict
    { gene : [ ghost domains list ]}
    
    Parameters
    ----------
    dict_gene_domainsModules : dict
        Dictionary with gene node name (str) as key and a dictionary of domains and modules in it as value (domain being a key, list of module as a value)
            
    Returns
    -------
    dict_gene_ghostDomains : dict
        Dictionary with gene node name (str) as key and a list of domains as value
    """
    # Init the dict
    dict_gene_ghostDomains = {gene : [] for gene in dict_gene_domainsModules.keys()}
    # Iterate on genes
    for gene, domainsModule in dict_gene_domainsModules.items():
        # Iterate on all the domains of the actual gene
        for domain, modules_list in domainsModule.items():
            # If there is no modules correspinding this domains
            if len(modules_list) == 0:
                # Add it to the dict, as a ghost domain
                dict_gene_ghostDomains[gene].append(domain)
    # Return the dict
    return dict_gene_ghostDomains

#==============================================================================
# Get list of genes containing a module composition / any modules of the composition
#==============================================================================

def genes_containing_modulesComp_any(modules_composition, dict_gene_moduleList) -> list:
    """
    Search the genes containing at least 1 module of a module composition
    
    Parameters
    ----------
    modules_composition : list
        A list of modules (a module being a c_module object)
    dict_gene_moduleList : dict
        Dictionary with gene node name (str) as key and list of the associated modules as value (a module being a c_module object)
            
    Returns
    -------
    gene_list : list
        A list of gene node name (str)
    """
    # Init the gene list
    gene_list = []
    # For each gene
    for gene, module_list in dict_gene_moduleList.items():
        # For the actual module, look at all the modules of the gene to compare
        for g_module in module_list:
            # Look at all the modules of the module composition of interest
            for i_module in modules_composition:
                # It the actual domains module is in the gene module list (compa only name)
                if i_module.name.split("|")[0] == g_module.name.split("|")[0]:
                    # Add the gene to the list of the gene with this module
                    if gene not in gene_list: gene_list.append(gene)
    # Return the gene list
    return gene_list

def genes_containing_modulesComp_whole(modules_composition, dict_gene_moduleList) -> list:
    """
    Search the genes containing the whole module composition
    
    Parameters
    ----------
    modules_composition : list
        A list of modules (a module being a c_module object)
    dict_gene_moduleList : dict
        Dictionary with gene node name (str) as key and list of the associated modules as value (a module being a c_module object)
            
    Returns
    -------
    gene_list : list
        A list of gene node name (str)
    """
    # Init the gene list
    gene_list = []
    # For each gene
    for gene, module_list in dict_gene_moduleList.items():
        # Init the modules found count
        modules_found = 0
        # For the actual module, look at all the modules of the gene to compare
        for g_module in module_list:
            # Look at all the modules of the module composition of interest
            for i_module in modules_composition:
                # It the actual domains module is in the gene module list (compa only name)
                if i_module.name.split("|")[0] == g_module.name.split("|")[0]:
                    # Add the gene to the list of the gene with this module
                    if gene not in gene_list: 
                        gene_list.append(gene)
                        modules_found += 1
        if modules_found != len(modules_composition):
            if gene in gene_list: gene_list.remove(gene)
    # Return the gene list
    return gene_list

#==============================================================================
# Write output itol files
#==============================================================================

def write_itol_modulesComposChange(dict_gene_moduleChange, filename) -> None:
    """
    Make itol pie chart file with number of modules composotions change
    at internals nodes, 
    1 modules, gain or lost = 1 change
    using the gene tree topology
    (what change in the module composition at the node, in comparaison of the parents)
    
    Parameters
    ----------
    dict_gene_moduleChange : dict
        A dictionary with gene name (str) as key and list of modules changes [list of modules gained and list of modules lost] as value 
    filename : str
        Name of the output file (itol txt format)
    """
    # Init the itol PIE Chart file
    itol_str = "DATASET_PIECHART\n"
    itol_str += "SEPARATOR COMMA\nDATASET_LABEL,Pie : modules gains / losts\nCOLOR,#990000\n"
    # Init 2 color, 1 for gain, 1 for lost
    itol_str += "FIELD_COLORS,#00FF00,#990000\nFIELD_LABELS,gain,lost\n"
    # Dataset time (info for all internals nodes)
    itol_str += "DATA\n"
    # All the change info are on the dict gene moduleChange
    for gene, moduleChange in dict_gene_moduleChange.items():
        # A moduleChange is :[  [news modules] , [lost modules] ]
        news_module = moduleChange[0]
        lost_module = moduleChange[1]
        if len(gene) > 1:
            itol_str += f"{gene},0,{len(news_module + lost_module)},{len(news_module)},{len(lost_module)}\n"
    # Write the itol file
    with open(filename, "w+") as itol_file:
        itol_file.write(itol_str)

def write_itol_annotationsChange(dict_gene_annotationsChange, filename) -> None:
    """
    Make itol pie chart file with number of annotations change
    at internals nodes, 
    1 annotation, gain or lost = 1 change
    using the gene tree topology
    (what change in the annotations composition at the node, in comparaison of the parents)
    
    Parameters
    ----------
    dict_gene_annotationsChange : dict
        Dictionary with gene node name (str) as key and list of changes (list, [list of gained annotations, list of lost annotations]) as value
    filename : str
        Name of the output file (itol txt format)
    """
    itol_str = "DATASET_PIECHART\n"
    itol_str += "SEPARATOR COMMA\nDATASET_LABEL,Pie : annotations gains / losts\nCOLOR,#0000FF\n"
    itol_str += "FIELD_COLORS,#0000FF,#A9A9A9\nFIELD_LABELS,gain,lost\n"
    itol_str += "DATA\n"
    for gene, annotationsChange in dict_gene_annotationsChange.items():
        news_annotation = annotationsChange[0]
        lost_annotation = annotationsChange[1]
        if len(gene) > 1:
            itol_str += f"{gene},0.5,{len(news_annotation + lost_annotation)},{len(news_annotation)},{len(lost_annotation)}\n"
    with open(filename, "w+") as itol_file:
        itol_file.write(itol_str)

def write_itol_moduleNb(filename, dict_gene_moduleList) -> None:
    """
    Itol simple bar file, with bar corresponding the number of modules at the leafs
    
    Parameters
    ----------
    filename : str
        Name of the output file (itol txt format)
    dict_gene_moduleList : dict
        Dictionary with gene node name (str) as key and list of the associated modules as value (a module being a c_module object)
    """
    itol_str = "DATASET_SIMPLEBAR\n"
    itol_str += "SEPARATOR COMMA\nDATASET_LABEL,Modules number\nCOLOR,#00FF00\n"
    itol_str += "WIDTH,500\nSHOW_VALUE,1\n"
    itol_str += "DATA\n"
    for gene, moduleList in dict_gene_moduleList.items():
        if len(gene) > 1:
            itol_str += f"{gene},{len(moduleList)}\n"
    with open(filename, "w+") as itol_file:
        itol_file.write(itol_str)
        
def write_itol_mod_heatmap(filename, dict_gene_moduleList) -> None:
    """
    Itol heatmap file, with presence / absence of modules at the leafs
    
    Parameters
    ----------
    filename : str
        Name of the output file (itol txt format)
    dict_gene_moduleList : dict
        Dictionary with gene node name (str) as key and list of the associated modules as value (a module being a c_module object)
    """
    itol_str = "DATASET_HEATMAP\n"
    itol_str += "SEPARATOR COMMA\nDATASET_LABEL,Modules presence\nCOLOR,#2F4F4F\n"
    all_modules = []
    for gene, module_list in dict_gene_moduleList.items():
        for module in module_list:
            if module.name not in all_modules:
                all_modules.append(module.name)
    all_modules = natsorted(list(set(all_modules)))
    itol_str += f"FIELD_LABELS,{','.join(all_modules)}\n"
    itol_str += "STRIP_WIDTH,2\nSHOW_LABELS,0\n"
    itol_str += f"COLOR_MIN,#FFFFFF\n"
    itol_str += f"COLOR_MAX,#00ff00\n"
    itol_str += "DATA\n"
    for gene, module_list in dict_gene_moduleList.items():
        if len(gene) > 1:
            itol_str += f"{gene}"
            for module in all_modules:
                if module in [m.name for m in module_list]:
                    itol_str += ",1"
                else:
                    itol_str += ",0"
            itol_str += "\n"
    with open(filename, "w+") as itol_file:
        itol_file.write(itol_str)
        
def write_itol_annot_heatmap(filename, dict_gene_annotationsList) -> None:
    """
    Itol heatmap file, with presence / absence of modules at the leafs
    
    Parameters
    ----------
    filename : str
        Name of the output file (itol txt format)
    dict_gene_annotationsList : dict
        Dictionary with gene node name (str) as key and list of the associated annotations as value (an annotation being a str)
    """
    itol_str = "DATASET_HEATMAP\n"
    itol_str += "SEPARATOR COMMA\nDATASET_LABEL,Annotations presence\nCOLOR,#2F4F4F\n"
    all_annot = []
    for gene, annotations_list in dict_gene_annotationsList.items():
        for annot in annotations_list:
            if annot not in all_annot:
                all_annot.append(annot)
    itol_str += f"FIELD_LABELS,{','.join(all_annot)}\n"
    itol_str += "STRIP_WIDTH,2\nSHOW_LABELS,0\n"
    itol_str += f"COLOR_MIN,#FFFFFF\n"
    itol_str += f"COLOR_MAX,#2986cc\n"
    itol_str += "DATA\n"
    for gene, annot_list in dict_gene_annotationsList.items():
        if len(gene) > 1:
            itol_str += f"{gene}"
            for annot in all_annot:
                if annot in annot_list:
                    itol_str += ",1"
                else:
                    itol_str += ",0"
            itol_str += "\n"
    with open(filename, "w+") as itol_file:
        itol_file.write(itol_str)
        
def write_itol_presenceShape(filename, dict_gene_moduleList) -> None:
    """
    Itol shapes file, with presence / absence of modules at the leafs
    
    Parameters
    ----------
    filename : str
        Name of the output file (itol txt format)
    dict_gene_moduleList : dict
        Dictionary with gene node name (str) as key and list of the associated modules as value (a module being a c_module object)
    """
    itol_str = "DATASET_EXTERNALSHAPE\n"
    itol_str += "SEPARATOR COMMA\nDATASET_LABEL,Modules presence\nCOLOR,#2F4F4F\n"
    all_modules = []
    for gene, module_list in dict_gene_moduleList.items():
        for module in module_list:
            if module.name not in all_modules:
                all_modules.append(module.name)
    all_modules = natsorted(list(set(all_modules)))
    itol_str += f"FIELD_COLORS{',#2F4F4F'*len(all_modules)}\n"
    itol_str += f"FIELD_LABELS,{','.join(all_modules)}\n"
    #itol_str += f"LEGEND_SHAPES,{','.join('1'*len(all_modules))}\n"
    itol_str += "SHAPE_SPACING,0\n"
    itol_str += "SHAPE_TYPE,1\n"
    itol_str += "DATA\n"
    for gene, module_list in dict_gene_moduleList.items():
        if len(gene) > 1:
            itol_str += f"{gene}"
            for module in all_modules:
                if module in [m.name for m in module_list]:
                    itol_str += ",1"
                else:
                    itol_str += ",0"
            itol_str += "\n"
    with open(filename, "w+") as itol_file:
        itol_file.write(itol_str)

def write_itol_domain(filename, dict_gene_domainList, label="Module composition") -> None:
    with open(filename, "w+") as itol_file:
        itol_file.write("DATASET_DOMAINS\nSEPARATOR COMMA\n")
        itol_file.write(f"DATASET_LABEL,{label}\n")
        itol_file.write("COLOR,#ff0000\n")
        itol_file.write(f"SHOW_INTERNAL,0\nSHOW_DOMAIN_LABELS,0\nBORDER_WIDTH,1\n")
        itol_file.write("WIDTH,1000\nHEIGHT_FACTOR,0.6\n")
        itol_file.write("DATA\n")
        for gene, domainList in dict_gene_domainList.items():
            if gene.startswith("G"):
                continue
            else:
                max_pos = 0
                for domain in domainList:
                    if int(domain.end) > int(max_pos):
                        max_pos = domain.end
                #max_pos = int(max_pos) - 1
                if int(max_pos) > 1:
                    domainList_noEnd = [domain for domain in domainList if domain.name != "end" and int(domain.end) != int(domain.start)]
                    itol_file.write(f"{gene},{max_pos}")
                    if len(domainList_noEnd) > 0:
                        itol_file.write(f",{','.join([domain.itol_str() for domain in domainList_noEnd])}")
                    itol_file.write("\n")
        
def write_itol_present(filename, dict_gene_domainList, dict_gene_moduleList, internal_gene_name) -> None:
    with open(filename, "w+") as itol_file:
        itol_file.write("DATASET_DOMAINS\nSEPARATOR COMMA\nSHOW_INTERNAL,0\nSHOW_DOMAIN_LABELS,0\nBORDER_WIDTH,1\n")
        itol_file.write(f"DATASET_LABEL,{internal_gene_name} present\n")
        itol_file.write("WIDTH,1000\nHEIGHT_FACTOR,0.4\nBACKBONE_HEIGHT,0\n")
        itol_file.write("MARGIN,-1000\n")
        itol_file.write("COLOR,#ff0000\nDATA\n")
        color = f"#BC8F8F"
        for gene, domainList in dict_gene_domainList.items():
            if gene.startswith("G"):    
                continue
            else:
                max_pos = 0
                for domain in domainList:
                    if int(domain.end) > int(max_pos):
                        max_pos = domain.end
                if gene in dict_gene_moduleList and int(max_pos) > 1:
                    moduleList = [module for module in dict_gene_moduleList[gene] 
                                  if int(module.end) != int(module.start) 
                                  and module.name in [m.name for m in dict_gene_moduleList[internal_gene_name]]
                                  ]
                    itol_file.write(f"{gene},{max_pos}")
                    if len(moduleList) > 0:
                        itol_file.write(f",{','.join([itol_square_color(module.start, module.end, module.name, color) for module in moduleList])}")
                    itol_file.write("\n")        
                    
def itol_square_color(start, end, name, color) -> str:
    return f"RE|{round(int(start),2)}|{round(int(end),2)}|{color}|{name}"
        
def write_itol_domain_gained_lost(filename, dict_gene_domainList, dict_gene_moduleList, internal_gene_name) -> None:
    with open(filename, "w+") as itol_file:
        itol_file.write("DATASET_DOMAINS\nSEPARATOR COMMA\nSHOW_INTERNAL,0\nSHOW_DOMAIN_LABELS,0\nBORDER_WIDTH,1\n")
        itol_file.write(f"DATASET_LABEL,{internal_gene_name} gained\n")
        itol_file.write("WIDTH,1000\nHEIGHT_FACTOR,0.4\nBACKBONE_HEIGHT,0\n")
        itol_file.write("MARGIN,-1000\n")
        itol_file.write("COLOR,#ff0000\nDATA\n")
        for gene, domainList in dict_gene_domainList.items():
            if gene.startswith("G"):    
                continue
            else:
                max_pos = 0
                for domain in domainList:
                    if int(domain.end) > int(max_pos):
                        max_pos = domain.end
                if gene in dict_gene_moduleList and int(max_pos) > 1:
                    moduleList_gain = [module for module in dict_gene_moduleList[gene][0] if int(module.end) != int(module.start)]
                    moduleList_lost = [module for module in dict_gene_moduleList[gene][1] if int(module.end) != int(module.start)]
                    itol_file.write(f"{gene},{max_pos}")
                    if len(moduleList_gain) > 0:
                        itol_file.write(f",{','.join([module.itol_str_gain() for module in moduleList_gain])}")
                    #if len(moduleList_lost) > 0:
                    #    itol_file.write(f",{','.join([module.itol_str_lost() for module in moduleList_lost])}") 
                    itol_file.write("\n")
        
def write_itol_domain_and_modules_gained_lost(filename, dict_gene_domainList, dict_gene_moduleList, internal_gene_name) -> None:
    with open(filename, "w+") as itol_file:
        itol_file.write("DATASET_DOMAINS\nSEPARATOR COMMA\nSHOW_INTERNAL,0\nSHOW_DOMAIN_LABELS,0\nBORDER_WIDTH,1\n")
        itol_file.write(f"DATASET_LABEL,{internal_gene_name} gained lost\n")
        itol_file.write("WIDTH,1000\nHEIGHT_FACTOR,0.6\n")
        itol_file.write("COLOR,#ff0000\nDATA\n")
        for gene, domainList in dict_gene_domainList.items():
            if gene.startswith("G"):    
                continue
            else:
                max_pos = 0
                for domain in domainList:
                    if int(domain.end) > int(max_pos):
                        max_pos = domain.end
                if gene in dict_gene_moduleList and int(max_pos) > 1:
                    moduleList_gain = [module for module in dict_gene_moduleList[gene][0] if int(module.end) != int(module.start)]
                    moduleList_lost = [module for module in dict_gene_moduleList[gene][1] if int(module.end) != int(module.start)]
                    domainList_noEnd = [domain for domain in domainList if domain.name != "end"]
                    itol_file.write(f"{gene},{max_pos}")
                    if len(domainList_noEnd) > 0:
                        itol_file.write(f",{','.join([domain.itol_str() for domain in domainList_noEnd])}")
                    if len(moduleList_gain) > 0:
                        itol_file.write(f",{','.join([module.itol_str_gain() for module in moduleList_gain])}")
                    if len(moduleList_lost) > 0:
                        itol_file.write(f",{','.join([module.itol_str_lost() for module in moduleList_lost])}") 
                    itol_file.write("\n")

def write_itol_popup(filename, gene_tree, dict_gene_domainList, dict_gene_moduleChange, dict_nodeName_annotationsList, dict_nodeName_annotationsChange) -> None:
    # Init 
    ncbi_url = "https://www.ncbi.nlm.nih.gov/protein/"
    quickGO_url = "https://www.ebi.ac.uk/QuickGO/term/"
    uniprot_url = "https://www.uniprot.org/uniprot/"
    itol_string = "POPUP_INFO\nSEPARATOR COMMA\nDATA\n"
    # Iterate on all nodes
    for node in gene_tree.iter_descendants():
        pred_string = ""
        # Get module composition, and parent module compositions, to compare them
        if node.name in dict_gene_domainList:
            interest_modules_list = [module.name for module in dict_gene_domainList[node.name]]
        # Get the change
        news_module = dict_gene_moduleChange[node.name][0]
        lost_module = dict_gene_moduleChange[node.name][1]
        if node.is_leaf():
            gene_name = node.name.split("_")[0].replace("P", "P_")
        else:
            gene_name = node.name
        itol_string += f"{node.name}, Protein informations, <h1>{gene_name}</h1>"
        # Node type
        if node.is_leaf():
            itol_string += f"<p style='color:blue'> {gene_name} -> Actual protein <a target='_blank' href='{ncbi_url}{gene_name}'>Ncbi protein link</a></p>"
        else:
            itol_string += f"<p style='color:blue'> {gene_name} -> Internal protein (ancestor inferred by phylo) </p>"
        # GOA annotations change (the full list is far to long, so only go for the changements)
        itol_string += f"<B style='color:black'> Annotation(s) ({len(dict_nodeName_annotationsList[node.name])}) </B> <p style='color:black'>"
        for func in dict_nodeName_annotationsList[node.name]:
            if "GO" in func: 
                func_url = quickGO_url
                uniprot = func.replace("O","O:")
                name = func.replace("O","O:")
            else: 
                func_url = uniprot_url
                if "_" in func:
                    uniprot = func.split("_")[0]
                    name = func.split("_")[1]
                else:
                    uniprot = ""
                    name = func
            itol_string += f"<a target='_blank' href='{func_url}{uniprot}'>{name} </a>"
        itol_string += "</p>"
        itol_string += f"<B style='color:green'> Annotation(s) win ({len(dict_nodeName_annotationsChange[node.name][0])}) </B> <p style='color:green'>"
        for func in dict_nodeName_annotationsChange[node.name][0]:
            if "GO" in func: 
                func_url = quickGO_url
                uniprot = func.replace("O","O:")
                name = func.replace("O","O:")
            else: 
                func_url = uniprot_url
                if "_" in func:
                    uniprot = func.split("_")[0]
                    name = func.split("_")[1]
                else:
                    uniprot = ""
                    name = func
            itol_string += f"<a target='_blank' href='{func_url}{uniprot}'>{name} </a>"
        itol_string += "</p>"
        itol_string += f"<B style='color:red'> Annotation(s) lost ({len(dict_nodeName_annotationsChange[node.name][1])}) </B> <p style='color:red'>"
        for func in dict_nodeName_annotationsChange[node.name][1]:
            if "GO" in func: 
                func_url = quickGO_url
                uniprot = func.replace("O","O:")
                name = func.replace("O","O:")
            else: 
                func_url = uniprot_url
                if "_" in func:
                    uniprot = func.split("_")[0]
                    name = func.split("_")[1]
                else:
                    uniprot = ""
                    name = func
            itol_string += f"<a target='_blank' href='{func_url}{uniprot}'>{name} </a>"
        itol_string += "</p>"
        # Modules compositions and their change
        interest_modules_list = natsorted(list(set(interest_modules_list)))
        news_module =  natsorted(list(set(news_module)))
        lost_module = natsorted(list(set(lost_module)))
        itol_string += f"<B style='color:black'> Module composition ({len(interest_modules_list)}) </B> <p style='color:black'>{' '.join(interest_modules_list)}</p>"
        itol_string += f"<B style='color:green'> Module(s) win ({len(news_module)}) </B> <p style='color:green'>{' '.join(news_module)}</p>"
        itol_string += f"<B style='color:red'> Module(s) lost ({len(lost_module)}) </B> <p style='color:red'>{' '.join(lost_module)}</p>"
        itol_string += "\n"
    with open(filename, "w+") as itol_file:
        itol_file.write(itol_string)

def write_itol_annotation(filename, dict_nodeName_annotations) -> None:
    # Init a dict with ipp : color
    dic_annot_itolStr = {}
    # Complete it, so each uniq ipp got a specific color (using uniq value list)
    for annot in list(set(val for dic in dict_nodeName_annotations.values() for val in dic)):
        # Attribute a symbol and a color for each annot (no given position and label)
        dic_annot_itolStr[annot] = f"{random.choice(['1','2','3','4','5'])},1,{random_hexaColor()},1"
    # Init the itol file
    with open(filename, "w+") as itol_file:
        itol_file.write("DATASET_SYMBOL\nSEPARATOR COMMA\nDATASET_LABEL,Func Annotations\nCOLOR,#ffff00\nMAXIMUM_SIZE,5\nDATA\n")
        # Now write the itol nodes infos, based on the nodeName_annot dic
        for nodeName, annot_list in dict_nodeName_annotations.items():
            # Init a position in case of multiple annot at a node
            pos = 0
            for annot in annot_list:
                itol_file.write(f"{nodeName},{dic_annot_itolStr[annot]},{pos},{annot}\n")
                pos += 0.15

def write_itol_connectTransfer(filename, list_module_gene_recipient) -> None:
    """
    Connection itol file for module intra-gene transfer
    """
    # Init the itol str
    itol_str = "DATASET_CONNECTION\n"
    itol_str += "SEPARATOR COMMA\nDATASET_LABEL,Module transfer\nCOLOR,#4682B4\n"
    itol_str += "DRAW_ARROWS,1\nCURVE_ANGLE,-65\n"
    # Data time
    itol_str += "DATA\n"
    for transfer in list_module_gene_recipient:
        itol_str += f"{transfer[1]},{transfer[2]},0.1,#4682B4,normal,{transfer[0].name}\n"
    # Write file
    with open(filename, "w+") as itol_file:
        itol_file.write(itol_str)

def write_itol_reconcSpGene(filename, dict_gene_spGeneEvent) -> None:
    """
    Branch symbol itol file, containing info the specie / gene reconciliations
    Events could be : Leaf mapping (dont really care) / Speciation / Gene duplication
    """
    #Init itol file
    itol_str = "DATASET_SYMBOL\n"
    itol_str += "SEPARATOR COMMA\nDATASET_LABEL,Genes-Sp events\nCOLOR,#BA55D3\n"
    itol_str += "MAXIMUM_SIZE,15\n"
    # Data time
    itol_str += "DATA\n"
    for gene, event in dict_gene_spGeneEvent.items():
        # We just care about the speciation and gene duplication
        # And we want to discrimine the 2 of them
        if event == "Speciation":
            itol_str += f"{gene},2,1,#8A2BE2,1,0.5,{event}\n"
        elif event == "Gene duplication":
            itol_str += f"{gene},3,3,#FF00FF,1,0.5,{event}\n"
    # Write file
    with open(filename, "w+") as itol_file:
        itol_file.write(itol_str)

def write_itol_geneMapping(filename, gene_list) -> None:
    """
    Branch symbol itol file 
    Highlight the gene of this list 
    """
    # Init itol file
    itol_str = f"DATASET_SYMBOL\n"
    itol_str += f"SEPARATOR COMMA\nDATASET_LABEL,{filename}\nCOLOR,#BA55D3\n"
    itol_str += f"MAXIMUM_SIZE,15\n"
    # Data time
    itol_str += "DATA\n"
    for gene in gene_list:
        itol_str += f"{gene},2,1,#8b008b,1,0,{filename}\n"
    # Write file
    with open(filename, "w+") as itol_file:
        itol_file.write(itol_str)

#==============================================================================
# Batch itol upload
#==============================================================================

def itol_uploader(zip_directory, api_k, project_n) -> None:
    """
    Use bash uploader api to uploade :
    gene tree + all his annotations files (in a zip directory)
    """
    # My itol api access key
    api_keys = str(api_k)
    # Param request
    batch_uploader = "https://itol.embl.de/batch_uploader.cgi"
    project_name = str(project_n)
    # Set param
    with open(zip_directory, "rb") as f:
        files = {"zipFile" : f}
        param = {
            "APIkey" : api_keys,
            "projectName" : project_name
            }
        # Request
        time.sleep(1)
        resp = requests.post(batch_uploader, files=files, data=param)

#==============================================================================
# Argparse parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("seadog_output",
                        help = "Seadog output (.output) text file",
                        type=str)
    parser.add_argument("gene_tree",
                        help = "Gene tree file (with branch length)",
                        type=str)
    parser.add_argument("--itol",
                        help = "Prepare zip file that could be upload on a registered Itol account",
                        action="store_true")
    parser.add_argument("--pastml_tab",
                        help = "Add an annotation pastml tab file for the differents nodes of the gene tree",
                        type=str)
    parser.add_argument("--pastml_directory",
                        help = "Use the pastml directory (with marginal proba)",
                        type=str)
    parser.add_argument("--domains_csv",
                        help = "Add an domains csv file, at the format : gene_name, domain_name, start, end",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

if __name__ == "__main__":
    args = parser()  
    # Read and extract seadog reconciliation output infos
    seadog_output = Path(args.seadog_output).resolve()
    gene_tree_file = Path(args.gene_tree).resolve()
    gene_tree, dict_module_mappingList, dict_module_modTree, gs_mapping_list = read_seadogO(seadog_output, gene_tree_file)

    # Infers modules compositions
    dict_gene_moduleList, dict_module_mappingList = infers_modulesCompo(gene_tree, dict_module_mappingList, dict_module_modTree)    
    dict_gene_moduleChange = make_dict_module_change(dict_gene_moduleList, gene_tree)

    # Some preparations for itols files 
    list_module_gene_recipient = make_module_gene_recipient_list(dict_gene_moduleList, dict_module_mappingList)
    dict_gene_spGeneEvent = {gs_mapping.entity_1 : gs_mapping.event for gs_mapping in gs_mapping_list}
    
    # If there is an annotations (= GOA terms) pastml tab file, make it in a dict gene : annotations
    if args.pastml_tab:
        dict_nodeName_annotationsList = load_pastml_ancestral_states(args.pastml_tab, gene_tree)
    else:
        dict_nodeName_annotationsList = {gene : [] for gene in dict_gene_moduleList.keys()}
    # We also want the GO termes changes on the tree, specially the gain (same algo / different function as module change)
    dict_nodeName_annotationsChange = make_dict_annotation_change(dict_nodeName_annotationsList, gene_tree)

    # If there is an domains file (ex Prosite / pfam) csv file
    dict_gene_domainsList = {gene : [] for gene in dict_gene_moduleList}
    if args.domains_csv:
        dict_gene_domainsList = load_domains_from_csv(dict_gene_moduleList.keys(), args.domains_csv) 
        dict_gene_domainsOList = domains_as_modules(dict_gene_domainsList)
        # Experimentals tests (time consuming)
        #dict_gene_domainsModules = dict_domains_in_modules(dict_gene_moduleList, dict_gene_domainsList)
        #dict_gene_infDomainsModules = infers_domains_from_modules(dict_gene_moduleList, dict_gene_domainsModules)
        #dict_gene_ghostDomains = ghost_domains(dict_gene_domainsModules)
        #dict_gene_freeModules = modules_not_in_domains(dict_gene_moduleList, dict_gene_domainsModules)

    # Write itols files
    directory = f"{seadog_output.parents[0]}/visuReconc_{seadog_output.stem}"
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Look at the modules change
    dict_gene_LeafsModulesThatChanged = look_at_modules_on_leafs(dict_gene_moduleChange, dict_gene_moduleList, gene_tree)
    #leaf_itol_directory = f"{seadog_output.parents[0]}/modulesChange_{seadog_output.stem}"
    modules_on_leafs_output(dict_gene_LeafsModulesThatChanged, dict_gene_moduleList, dict_gene_domainsOList, directory)

    write_itol_domain(f"{directory}/itolModules_{seadog_output.stem}.txt", dict_gene_moduleList, "Module composition")
    write_itol_domain(f"{directory}/itolDomains_{seadog_output.stem}.txt", dict_gene_domainsOList, "Pfam, Prosite domains")
    write_itol_popup(f"{directory}/itolPopup_{seadog_output.stem}.txt", 
                        gene_tree, dict_gene_moduleList, 
                        dict_gene_moduleChange, 
                        dict_nodeName_annotationsList,
                        dict_nodeName_annotationsChange)
    write_itol_modulesComposChange(dict_gene_moduleChange, f"{directory}/itol_modules_PieGainsLost_{seadog_output.stem}.txt")
    write_itol_annotationsChange(dict_nodeName_annotationsChange, f"{directory}/itol_ppi_PieGainsLost_{seadog_output.stem}.txt")
    write_itol_moduleNb(f"{directory}/itolBarModulesNb_{seadog_output.stem}.txt", dict_gene_moduleList)
    write_itol_mod_heatmap(f"{directory}/itolModPresence_{seadog_output.stem}.txt", dict_gene_moduleList)
    write_itol_annot_heatmap(f"{directory}/itolAnnotPresence_{seadog_output.stem}.txt", dict_nodeName_annotationsList)
    write_itol_connectTransfer(f"{directory}/itolModTransfer_{seadog_output.stem}.txt", list_module_gene_recipient)
    write_itol_reconcSpGene(f"{directory}/itolSpGeneEvents_{seadog_output.stem}.txt", dict_gene_spGeneEvent)
    write_itol_annotation(f"{directory}/itolGOt_{seadog_output.stem}.txt", {n : aL for n, aL in dict_nodeName_annotationsList.items() if n in gene_tree})
    
    # Write the gene tree
    gene_tree.write(format=1, outfile=f"{directory}/geneReconc_{seadog_output.stem}.tree", format_root_node=True)
    
    # Write the gene-species event .csv file (event at each gene nodes)
    write_spGeneEvent_csv(dict_gene_spGeneEvent, f"specieGeneEvent_{seadog_output.stem}.csv")

    # Write the .csv module compos and modules change
    write_module_compo_csv(dict_gene_moduleList, f"modulesCompo_{seadog_output.stem}.csv")
    write_module_change_csv(dict_gene_moduleChange, f"modulesChange_{seadog_output.stem}.csv")

    # Write the .csv function apparition and corresponding modules change
    write_function_module_change(dict_gene_moduleChange, dict_nodeName_annotationsChange, f"functionChange_moduleChange_{seadog_output.stem}.csv")
    write_complete_function_module_change(dict_gene_moduleChange, dict_nodeName_annotationsChange, f"complete_functionChange_moduleChange_{seadog_output.stem}.csv")
    write_function_module_change_expanded(gene_tree, dict_gene_moduleList, dict_gene_moduleChange, dict_nodeName_annotationsChange, f"functionChange_moduleChange_expand_{seadog_output.stem}.csv")

    # Zip gene tree with all his annotation files, and upload it on my itol account
    if args.itol:
        zip_directory = directory
        shutil.make_archive(zip_directory, "zip", directory) 
        #itol_uploader(f"{zip_directory}.zip", api_k, project_n)
        



