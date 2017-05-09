#!/usr/bin/python

__author__ = "Zhaolong Yu"
__copyright__ = "Copyright 2017"
__credits__ = ["Zhaolong Yu"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Zhaolong Yu"
__email__ = "zhaolong.yu@yale.edu"

### Usage:      python3 centrality_calculation.py -i <input file> -a <annotation file> -n <dip ppi file>
### Example:    python3 centrality_calculation.py -i input.txt -a map_table.txt -n dip.txt


import argparse
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import networkx as nx
import mygene
import seaborn

#os.chdir("/Users/jerome/Projects/protein_network")

parser = argparse.ArgumentParser(description='Centrality Calculation V1.0')
parser.add_argument('-i', '--input', help='input mutation file', required=True)
parser.add_argument('-a', '--annotation', help='uniprot annotation file', required=True)
parser.add_argument('-n', '--network', help='dip protein-protein interaction network', required=True)
args = parser.parse_args()

## open the PPI file
fo_dip = open(args.network,"r")
fo_dip.readline()
edges = []
dip_uniprot = {}
uniprot_dip = {}
for rawline in fo_dip:
    line = rawline.strip().split("\t")
    list1 = line[0].split("|")
    list2 = line[1].split("|")
    if len(list1)>0 and len(list2)>0:
        node1 = line[0].split("|")[0]
        node2 = line[1].split("|")[0]
        edges.append([node1,node2])
    ref = list1[-1].split(":")
    if len(list1)>1 and ref[0] == "uniprotkb":
        if list1[0] not in dip_uniprot.keys():
            dip_uniprot[list1[0]] = ref[1]
            uniprot_dip[ref[1]] = list1[0]
    ref = list2[-1].split(":")
    if len(list2)>1 and ref[0] == "uniprotkb":
        if list2[0] not in dip_uniprot.keys():
            dip_uniprot[list2[0]] = ref[1]
            uniprot_dip[ref[1]] = list2[0]

        
    
fo_dip.close()

## genes in the PPI network
all_genes = []
#for item in edges:
#    if item[0] not in all_genes:
#        all_genes.append(item[0])
#    if item[1] not in all_genes:
#        all_genes.append(item[1])
for gene in dip_uniprot.items():
    all_genes.append(gene[1])
    
## output the gene list
fo_gene = open("all_genes.txt","w")
for gene in all_genes:
    fo_gene.write(gene + "\n")
fo_gene.close()
    
## read the snp file
#fo_snp = open("Z.3DStruct_annotation.txt","r")
fo_snp = open(args.input,"r")

snp_genes = []
for rawline in fo_snp:
    line = rawline.strip().split("\t")
    gene = line[2].split(";")[1]
    snp_genes.append(gene)
fo_snp.close()   

## get the proteins containing SNPs 
unique_snp_genes = np.unique(snp_genes)
  
## gene annotation file
#fo_map = open("map_table.txt","r")
fo_map = open(args.annotation,"r")
fo_map.readline()
gene_uniprot = {}
uniprot_gene = {}
for rawline in fo_map:
    line = rawline.strip().split("\t")
    gene_uniprot[line[1]] = line[0]
    uniprot_gene[line[0]] = line[1]
    
fo_map.close()

## get the gene list
gene_list = []
for gene in all_genes:
    if gene in uniprot_gene.keys():
        gene_list.append(uniprot_gene[gene])

## gene set partitioning
unique_nonsnp_genes = set(gene_list)-set(unique_snp_genes)       
unique_snp_genes = set(unique_snp_genes).intersection(set(gene_list))        
    

## building the network
g = nx.Graph() 
g.add_edges_from(edges)   

g.number_of_nodes()  

snp_degree_centrality = {}
nonsnp_degree_centrality = {}

snp_betweenness_centrality = {}
nonsnp_betweenness_centrality = {}

degree_dict = {}
degree_dict = nx.degree(g)

## output the nodes with top10 degrees 
f = open('degree_top10.csv','w')
f.write("uniprot_id" + "," + "degree" + "\n")   
sorted_keys = sorted(degree_dict, key=degree_dict.get, reverse=True)
for r in sorted_keys[:10]:
    f.write(dip_uniprot[r] + "," + str(degree_dict[r]) + "\n")
f.close()
    
## calculate the degree centrality
tmp_dc = nx.degree_centrality(g)

## result dict for degree centrality
for gene in unique_nonsnp_genes:
    nonsnp_degree_centrality[gene] = tmp_dc[uniprot_dip[gene_uniprot[gene]]]

for gene in unique_snp_genes:
    snp_degree_centrality[gene] = tmp_dc[uniprot_dip[gene_uniprot[gene]]]

## calculate the betweenness centrality
tmp_bc = nx.betweenness_centrality(g,k=4901,normalized=True)

## result dict for betweenness centrality    
for gene in unique_nonsnp_genes:
    nonsnp_betweenness_centrality[gene] = tmp_bc[uniprot_dip[gene_uniprot[gene]]]
     
for gene in unique_snp_genes:
    snp_betweenness_centrality[gene] = tmp_bc[uniprot_dip[gene_uniprot[gene]]]

## result list (for further analysis)    
snp_degree_centrality_list = []
nonsnp_degree_centrality_list = []

snp_betweenness_centrality_list = []
nonsnp_betweenness_centrality_list = []
      
for gene in snp_degree_centrality.keys():
    snp_degree_centrality_list.append(float(snp_degree_centrality[gene]))
    
for gene in nonsnp_degree_centrality.keys():
    nonsnp_degree_centrality_list.append(float(nonsnp_degree_centrality[gene]))
    
for gene in snp_betweenness_centrality.keys():
    snp_betweenness_centrality_list.append(float(snp_betweenness_centrality[gene]))
    
for gene in nonsnp_betweenness_centrality.keys():
    nonsnp_betweenness_centrality_list.append(float(nonsnp_betweenness_centrality[gene]))

## plot the value distribution 
plt.hist(snp_degree_centrality_list)
plt.title("Degree Centrality of Proteins with SNPs")    

plt.hist(nonsnp_degree_centrality_list)
plt.title("Degree Centrality of Proteins without SNPs")    

plt.hist(snp_betweenness_centrality_list)
plt.title("Betweenness Centrality of Proteins with SNPs")    

plt.hist(nonsnp_betweenness_centrality_list)
plt.title("Betweenness Centrality of Proteins without SNPs")    

## sorting results and output the top10 list
f = open('nonsnp_degree_centrality.csv','w')
f.write("uniprot_id" + "," + "value" + "\n")
sorted_keys = sorted(nonsnp_degree_centrality, key=nonsnp_degree_centrality.get, reverse=True)
for r in sorted_keys:
    f.write(gene_uniprot[r] + "," + str(nonsnp_degree_centrality[r]) + "\n")
f.close()

f = open('snp_degree_centrality.csv','w')
f.write("uniprot_id" + "," + "value" + "\n")    
sorted_keys = sorted(snp_degree_centrality, key=snp_degree_centrality.get, reverse=True)
for r in sorted_keys:
    f.write(gene_uniprot[r] + "," + str(snp_degree_centrality[r]) + "\n")
f.close()

f = open('nonsnp_betweenness_centrality.csv','w')
f.write("uniprot_id" + "," + "value" + "\n")    
sorted_keys = sorted(nonsnp_betweenness_centrality, key=nonsnp_betweenness_centrality.get, reverse=True)
for r in sorted_keys:
    f.write(gene_uniprot[r] + "," + str(nonsnp_betweenness_centrality[r]) + "\n")
f.close()

f = open('snp_betweenness_centrality.csv','w')
f.write("uniprot_id" + "," + "value" + "\n")        
sorted_keys = sorted(snp_betweenness_centrality, key=snp_betweenness_centrality.get, reverse=True)
for r in sorted_keys:
    f.write(gene_uniprot[r] + "," + str(snp_betweenness_centrality[r]) + "\n")
f.close()

## betweenness centrality table
f = open('betweenness_centrality_all.csv','w')
f.write("uniprot_id" + "," + "value" + "\n")    
sorted_keys = sorted(tmp_bc, key=tmp_bc.get, reverse=True)
for r in sorted_keys:
    if r in dip_uniprot.keys():
        f.write(dip_uniprot[r] + "," + str(tmp_bc[r]) + "\n")
f.close()  

## degree centrality table
f = open('degree_centrality_all.csv','w')
f.write("uniprot_id" + "," + "value" + "\n")    
sorted_keys = sorted(tmp_dc, key=tmp_dc.get, reverse=True)
for r in sorted_keys:
    if r in dip_uniprot.keys():
        f.write(dip_uniprot[r] + "," + str(tmp_dc[r]) + "\n")
f.close()  
