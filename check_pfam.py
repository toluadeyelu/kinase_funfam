import networkx as nx
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles

#import csv
#import cPickle
#import lzma
#import gzip

#HG = nx.read_gpickle("C:/Users/Adeyelu-Tolu/Dropbox/human_db.pkl")

pfam_dict = defaultdict(list)
uniprot_pfam = defaultdict(list)
all_kinase = set()
#for c, line in enumerate(open("C:/Users/Adeyelu-Tolu/Dropbox/pfam_map.tsv")):
for c, line in enumerate(open("/cath/homes2/ucbttad/Dropbox/pfam_map.tsv")):
    if c ==0:continue
    line = line.strip().split("\t")
    pfam = line[0]
    uniprot = line[1]
    all_kinase.add(uniprot)
    pfam_dict[pfam].append(uniprot)
    
#for c,line in enumerate(open("C:/Users/Adeyelu-Tolu/Dropbox/kinase_uniprot.csv")):
for c, line in enumerate(open("/cath/homes2/ucbttad/Dropbox/kinase_uniprot.csv")):
    if c==0:continue
    line = line.strip().split("\t")
    uniprot_id = line[1]
    for k, v in pfam_dict.items():
        for i in v:
            if i!=uniprot_id:continue
            #print i,k 

genes_disease_map = defaultdict(list)
for c, line in enumerate(open("/cath/homes2/ucbttad/Dropbox/disease_genes.txt")):
#for c,line in enumerate(open("C:/Users/Adeyelu-Tolu/Dropbox/disease_genes.txt")):
    line = line.strip().split("\t")
    disease = line[0]
    genes = line[2]
    genes_disease_map[genes].append(disease)
gene_uniprot ={}
for c, line in enumerate(open("/cath/homes2/ucbttad/Dropbox/kinase_gene_mapping.tsv")):
#for c, line in enumerate(open("C:/Users/Adeyelu-Tolu/Dropbox/kinase_gene_mapping.tsv")):
    if c==0:continue
    line = line.strip().split("\t")
    uniprot_id = line[0]
    kin_genes = line[2]
    gene_uniprot[kin_genes]=uniprot_id
    #print kin_genes, genes_disease_map.get(kin_genes)

intside_proteins = set()
'''Proteins involved with side effects downloaded from InTside'''
for c, line in  enumerate(open("/cath/homes2/ucbttad/Dropbox/intside_associated_proteins.tsv")):
    if c == 0:continue
    line = line.strip().split("\t")
    sd_prot = line[2]
    intside_proteins.add(sd_prot)


targeted_kinases = set()
for c , line in enumerate(open("/cath/homes2/ucbttad/Dropbox/drug-ff_overreptest_over.tsv")):
    if c ==0:continue
    line = line.strip().split("\t")
    target = line[0]
    targeted_kinases.add(target)

'''substrate to uniprot_ID'''
kinsub_map = {}
#for c, line in enumerate(open("C:/Users/Adeyelu-Tolu/Dropbox/kinsubstrate_map.tsv")):
for c, line in enumerate(open("/cath/homes2/ucbttad/Dropbox/kinsubstrate_map.tsv")):
    line = line.strip().split("\t")
    uniprot_id = line[0]
    gene_id = line[1]
    gene_id = gene_id[:-6]
    kinsub_map[gene_id]=uniprot_id

kinase_substrate = defaultdict(list)

kin_sub = nx.Graph()
subst_set = set()
kin_set = set()

#for c, line in enumerate(open("C:/Users/Adeyelu-Tolu/Dropbox/rawKSI.csv")):
for c, line in enumerate(open("/cath/homes2/ucbttad/Dropbox/rawKSI.csv")):
    if c==0:continue
    line = line.strip().split("\t")
    kin_gene = line[0]
    subst = line[1]
    kin_gene = gene_uniprot.get(kin_gene)
    if kin_gene is None:continue
    subst = kinsub_map.get(subst)
    if subst is None:continue
    kin_sub.add_edge(kin_gene, subst)
    #print kin_gene, subst
    subst_set.add(subst)
    kin_set.add(kin_gene)
    kinase_substrate[kin_gene].append(subst)
all_essentials = set()
#for c,line in enumerate(open("C:/Users/Adeyelu-Tolu/Dropbox/essential_genes.tsv")):
for c, line in enumerate(open("/cath/homes2/ucbttad/Dropbox/essential_genes.tsv")):
    if c ==0:continue
    line = line.strip().split("\t")
    org = line[7]
    if org !="Homo sapiens":continue
    essen = line[12]
    all_essentials.add(essen)

'''plotting venn diagram for the intersections of the set'''
s=(all_essentials, targeted_kinases, intside_proteins)
v = venn3(subsets=s, set_labels= ("Essential genes","Targeted Kinases","Side effect Proteins"))
#v.get_label_by_id()
#v.get_patch_by_id()

plt.show()

exit()


print len(all_essentials & targeted_kinases)
print len(intside_proteins & all_essentials)
print len(intside_proteins & targeted_kinases)
print len(targeted_kinases & intside_proteins & all_essentials)
print len(all_essentials)
print len(targeted_kinases)
print len(intside_proteins)