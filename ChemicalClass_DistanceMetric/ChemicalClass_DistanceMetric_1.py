#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 10:38:12 2018

@author: madeleineernst
"""
####################################################################################################
#                                                                                                  #
#  Create binary matrices of ClassyFire chemical classes at the direct parent and subclass level   #
#                                                                                                  #
####################################################################################################

import pandas as pd

cf = pd.read_csv('ClassyFire/ClassyFire_Output_forCytoscape.tsv',sep='\t') 
ft = pd.read_csv('FeatureTable_Rhamnaceae.csv',sep=',') 
ft = ft.rename(columns = {'row ID':'cluster.index'})

comb = pd.merge(cf, ft, on="cluster.index")

rem = ['cluster.index', 'CF_componentindex', 'CF_substituents',
       'CF_score', 'CF_kingdom_scores', 'CF_superclass_scores',
       'CF_class_scores', 'CF_subclass_scores', 'CF_directparent_scores',
       'CF_substituents_scores']

comb = comb.drop(rem, 1)

taxons = ['CF_kingdom','CF_superclass','CF_class','CF_subclass','CF_directparent']
hierarchical_df = comb.groupby(taxons).sum() #sum or whatever is most appropiate for your data

hierarchical_df.to_csv("hierarchical_ClassyFire.tsv",sep='\t',index=True)


########### Create binary matrix of ClassyFire chemical classes at the direct parent level
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('hierarchical_ClassyFire.tsv',sep='\t')
taxons = ['CF_kingdom','CF_superclass','CF_class','CF_subclass','CF_directparent']
df = df[taxons]
df['type'] = df[['CF_kingdom','CF_superclass','CF_class','CF_subclass']].values.tolist()

en = pd.get_dummies(df.type.apply(lambda x: pd.Series([i for i in x])))
en.index = df.CF_directparent
en = en[~en.index.duplicated(keep='first')]

en.to_csv("classlist_directparents.tsv",sep='\t',index=True)

########### Create binary matrix of ClassyFire chemical classes at the subclass level
import pandas as pd

df = pd.read_csv('hierarchical_ClassyFire.tsv',sep='\t')
taxons = ['CF_kingdom','CF_superclass','CF_class','CF_subclass']
df = df[taxons]
df['type'] = df[['CF_kingdom','CF_superclass','CF_class']].values.tolist()

en = pd.get_dummies(df.type.apply(lambda x: pd.Series([i for i in x])))
en.index = df.CF_subclass
en = en[~en.index.duplicated(keep='first')]

en.to_csv("classlist_subclass.tsv",sep='\t',index=True)

####################################################################################################
#                                                                                                  #
# Create feature tables of ClassyFire chemical classes at the direct parent and subclass level     #
#                                                                                                  #
####################################################################################################

########### create count matrix of chemical subclasses     
import pandas as pd

cf = pd.read_csv('ClassyFire/ClassyFire_Output_forCytoscape.tsv',sep='\t') 
ft = pd.read_csv('FeatureTable_Rhamnaceae.csv',sep=',') 

#introduce cut-off of min. intensity 1000
rowID = list(ft['row ID'])
ft[ft < 1000] = 0
ft['row ID'] = rowID
#introduce cut-off of min. intensity 1000

ft = ft.rename(columns = {'row ID':'cluster.index'})

comb = pd.merge(cf, ft, on="cluster.index")

subcl = comb.CF_subclass.unique()

subcl_df = []
for i in range(len(subcl)):
    sel = comb.loc[comb['CF_subclass'] == subcl[i]] 
    out = sel.astype(bool).sum(axis=0)
    out = out.to_dict()
    subcl_df.append(out)

df = pd.DataFrame(subcl_df)
df.insert(loc=0, column='id', value= list(subcl))
df = df.drop('cluster.index', 1)
df = df[df.columns.drop(list(df.filter(regex='CF')))] # this dataframe contains samples in columns and subclasses in rows, the numeric values describe number of molecules within the corresponding subclasses

df.to_csv("featuretable_subclasses_cutoff1000.tsv",sep='\t',index=False)


########### create count matrix of direct parents           
import pandas as pd

cf = pd.read_csv('ClassyFire/ClassyFire_Output_forCytoscape.tsv',sep='\t') 
ft = pd.read_csv('FeatureTable_Rhamnaceae.csv',sep=',') 

#introduce cut-off of min. intensity 1000
rowID = list(ft['row ID'])
ft[ft < 1000] = 0
ft['row ID'] = rowID
#introduce cut-off of min. intensity 1000

ft = ft.rename(columns = {'row ID':'cluster.index'})

comb = pd.merge(cf, ft, on="cluster.index")

subcl = comb.CF_directparent.unique()

subcl_df = []
for i in range(len(subcl)):
    sel = comb.loc[comb['CF_directparent'] == subcl[i]] 
    out = sel.astype(bool).sum(axis=0)
    out = out.to_dict()
    subcl_df.append(out)

df = pd.DataFrame(subcl_df)
df.insert(loc=0, column='id', value= list(subcl))
df = df.drop('cluster.index', 1)
df = df[df.columns.drop(list(df.filter(regex='CF')))] # this dataframe contains samples in columns and subclasses in rows, the numeric values describe number of molecules within the corresponding direct parents

df.to_csv("featuretable_directparent_cutoff1000.tsv",sep='\t',index=False)

# Use ChemicalClass_DistanceMetric_2.R to calculate chemically informed distance matrices based on ClassyFire chemical ontology

######################################
# Create mapping file                #
######################################

# read metadata file and distance matrix
import pandas as pd

dist = pd.read_csv('BioSynDist_cutoff1000_directparents.tsv',sep='\t') 
md = pd.read_csv('MetaData_Rhamnaceae.txt',sep='\t') 

md['filename'] = [w.replace(' ', '.') for w in md['filename']]
md['filename'] = [w.replace('(', '.') for w in md['filename']]
md['filename'] = [w.replace(')', '.') for w in md['filename']]

md = md.rename(columns = {'filename':'#SampleID'})
md.insert(loc=1, column='BarcodeSequence', value= 'GATACA')
md.insert(loc=2, column='LinkerPrimerSequence', value= 'GATACA')
md.insert(loc=3, column='Description', value= 'Metabolome')
                
md.to_csv("MappingFile_BioSynDist.txt", sep='\t',index=False)

######################################
# Create Emperor plot                #
######################################

# Execute in command line
# source activate qiime2-2018.2
# qiime tools import --input-path BioSynDist_cutoff1000_directparents.tsv --output-path BioSynDist_cutoff1000.qza --type 'DistanceMatrix'
# qiime diversity pcoa --i-distance-matrix BioSynDist_cutoff1000.qza --o-pcoa BioSynDist_pcoa_cutoff1000.qza
# qiime emperor plot --i-pcoa BioSynDist_pcoa_cutoff1000.qza --m-metadata-file MappingFile_BioSynDist.txt --o-visualization BioSynDist_visualization_cutoff1000.qzv

# qiime tools view BioSynDist_visualization_cutoff1000.qzv
