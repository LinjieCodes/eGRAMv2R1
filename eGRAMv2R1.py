'''
Several notes:
1.  There are two kinds of input files - (a) user data and (b) system data. The former has three files - (a) expression profile, (b) lncRNA target
    prediction, (c) TF target prediction. LncRNAs' targets are predicted using LongTarget (or the more rapid version Farsim); TFs' targets are
    predicted using any method, such as CellOracle. Since the primary goal of eGRAM is analyzing lncRNA function in transcriptional regulation, the
    TF target prediction file is optional. System data include human/mouse KEGG and wikipathway pathways downloaded from the KEGG and Wikipathways
    websites.
2.  All variables begin with their data type: df_, list_, dict_, dict2D_, set_.
3.  It is the user's responsibility to determine correct gene symbols and gene types, and keep gene symbols and types consistent in the input files.
    In the given examples, targets of lncRNAs and TFs are only lncRNA genes, protein-coding genes, and TF genes. The user can explore other gene
    types (but they are barely annotated in pathway databases).
4.  eGRAM performs transcription analysis based on both correlated gene expression and lncRNA/TF target prediction. For small sample size, Pearson
    correlation is recommended, for large sample size, Spearman correlation can also be chosen.
5.  eGRAM identifies "all-mutually-correlated" lncRNAs and "all-mutually-correlated" TFs from lncRNAs and TFs in the gene expression file. We call
    this method "COCO", which is akin to clustering genes into clusters allowing overlap. Clustering with overlapping is reasonable because
    a lncRNA may work with different partners in different situations to regulate transcription, so do TFs. Two command-line arguments,
    --lncCorr and --tfCorr, work with the option. Since there are many clustering algorithm, future versions will allow the use of other algorithms
    (e.g. FLAME).
6.  To make eGRAM able to handle lncRNAs and TFs of any names in any numbers, we organize lncRNAs and TFs into dictionaries, each lncRNA/TF is the
    key and its COCO set (a list) is the value.
7.  Several pending issues
    (1) The relations between regulators and targets may be (a) co-upregulation, (b) co-downregulation, (c) differential expression - regulator high
    and target low, (d) differential expression - regulator low and target high. --- we have not explore this???
    (2) Whether a cutoff is needed for controlling the smallest COCO set size. For example, a set of <=3 lncRNAs or a set of <=2 TFs may not indicate
    a reliable regulator set.
    (3) There are three kinds of modules - lncRNAs' target sets, TFs' target sets, and merged target sets (intersection and union). Pathway enrichment
    analysis can be applied to each of the sets to reveal whether lncRNAs' target sets and TFs' target sets are enriched in different pathways.
    (4) Examining the intersection or union of merged sets is also biologically different, indicating co-regulation coordinated regulation by lncRNAs
    and TFs, respectively.
8.  A gene with all 0s is removed, but otherwise kept. 0s in RNA-seq datasets are assumed true value, so no imputation is performed. If scRNA-seq data
    is analyzed, it is the user's responsibility to perform imputation using any method. We assume that imputation is finished before running eGRAM.

Main running steps:
    step 1 - read in user files and command line arguments
    step 2 - read in kegg files and wikipathway files
    step 3 - computing lncRNA-gene correlations, lncRNA-lncRNA correlations, and upon which, lncRNAsets
    step 4 - computing TF-gene correlations, TF-TF correlations, and upon which, TFsets.
    step 5 - generate shared lncRNA target sets upon df_lncRNA_gene_DBS and df_lncRNA_gene_corr
    step 6 - generate shared TF target sets upon df_TF_gene_DBS and df_TF_gene_corr
    step 7 - identify typically co-regulated modules and independently-regulated modules
    step 8 - perform pathway enrichment analysis for module genes in df_path_module
    step 9 - write results to files
    step 10 - generate files for cytoscape

An explanatory command line:
    eGRAMv2R1.py --exp 2024May-DEG-exp-A549-2513WT.csv --lncDBS 2024May-lncRNA-DBS-A549-2513.csv --tfDBS 2024May-TF-DBS-A549-2513.csv --species 1
                 --lncCutoff 100 --tfCutoff 8 --lncCorr 0.6 --tfCorr 0.6 --moduleSize 50 --corr Pearson --cluster COCO --fdr 0.01 --out outA5492513WT
The user can output intermediate results to a specific file by adding "> file_name" at the end of the command line in ther Terminal mode.
'''


#!/usr/bin/python
import numpy as np
import pandas as pd
import re
import os
import sys, getopt
from scipy.stats import hypergeom
#from sklearn.impute import SimpleImputer

########################################################################################################################
# step 2 - read in kegg files and wikipathway files
def read_kegg(geneFile, pathFile, linkFile):
    '''
    Read kegg genes, pathways and pathway-gene links.
    '''
    gene_symbols = {}  # debug info: gene_symbols={dict:0} {}, 0 indicating no element.
    kegg_pathway = {}  # debug info: gene_symbols={dict:0} {}.
    kegg_link = {}
    gene_id_set = set()  # debug info: gene_id_set={set:0} {}.

    with open(geneFile) as f:
        for line in f:  # debug: first = 'hsa:102466751\tmiRNA\t1:complement(17369..17436)\tMIR6859-1, hsa-mir-6859-1; microRNA 6859-1\n'
            cols = line.strip('\r\n').split('\t')  # debug: ['hsa:102466751', 'miRNA', '1:complement(17369..17436)', 'MIR6859-1, hsa-mir-6859-1; microRNA 6859-1']
            gene_id = cols[0]  # debug: gene_id = 'hsa:102466751', the 1st element in cols, i.e., cols[0].
            descri = cols[-1]  # descri is the last element of cols, e.g., 'MIR6859-1, hsa-mir-6859-1; microRNA 6859-1'.
            if (',' in descri) or (';' in descri):
                symbol_descri = descri.split(';')[0]
                symbols = symbol_descri.split(", ")
                if gene_id not in gene_symbols:  #
                    gene_symbols[gene_id] = set()
                for symbol in symbols:
                    gene_symbols[gene_id].add(symbol)
                gene_id_set.add(gene_id)

    with open(pathFile) as f:
        for line in f:
            pathID, pathwayName = line.strip('\r\n').split('\t')
            pathID = pathID.replace("path:", '')
            pathwayName = pathwayName[:pathwayName.rfind(' -')]
            kegg_pathway[pathID] = pathwayName

    with open(linkFile) as f:
        for line in f:
            pathID, gene_id = line.strip('\r\n').split('\t')
            pathID = pathID.replace("path:", '')
            if pathID in kegg_pathway:
                if gene_id in gene_id_set:
                    if pathID not in kegg_link:
                        kegg_link[pathID] = set()
                    kegg_link[pathID].add(gene_id)

    totalGeneNum = len(gene_id_set)  # totalGeneNum = 22101

    kegg_gene = {}
    for gene_id in gene_symbols:
        for gene_symbol in gene_symbols[gene_id]:
            kegg_gene[gene_symbol] = gene_id

    return kegg_gene, kegg_pathway, kegg_link, totalGeneNum
	
def read_wikip(wikiDir, kegg_gene):
    '''
    Read genes and pathways in Wikipathways.
    '''
    wiki_pathway = {}
    wiki_link = {}

    pathway_files = os.listdir(wikiDir)
    pathway_comp = re.compile(' Name="(.+?)" ')
    gene_comp = re.compile(' TextLabel="(.+?)" ')
    for pathway_file in pathway_files:
        pathway_id = pathway_file.split('_')[-2]
        with open(wikiDir + pathway_file, encoding='utf-8') as f:
            for line in f:
                pathway_re = re.search(pathway_comp, line)
                if pathway_re:
                    pathway = pathway_re.group(1)
                    wiki_pathway[pathway_id] = pathway_re.group(1)
                gene_re = re.search(gene_comp, line)
                if gene_re:
                    symbol = gene_re.group(1)
                    if symbol in kegg_gene:
                        gene_id = kegg_gene[symbol]
                        if pathway_id in wiki_pathway:
                            if pathway_id not in wiki_link:
                                 wiki_link[pathway_id] = set()
                            wiki_link[pathway_id].add(gene_id)

    return wiki_pathway, wiki_link


########################################################################################################################
# step 3 - computing lncRNA-gene correlations, lncRNA-lncRNA correlations, and upon which, lncRNAsets
def  generate_lncRNA_files(df_gene_sample_exp):
    df_lncRNA_sample_exp = df_gene_sample_exp.copy(deep=True)
    # When deep=True (default), a new object is created.
    df_lncRNA_sample_exp = df_lncRNA_sample_exp.drop(df_lncRNA_sample_exp[df_lncRNA_sample_exp['genetype'] != 'lncRNA'].index)
    # dropping non-lncRNA rows.
    df_lncRNA_sample_exp = df_lncRNA_sample_exp.reset_index()
    # "drop" operation does change the original index.
    # This .reset_index() re-sets index (also keeps the original index in df_gene_sample_exp).
    df_lncRNA_sample_exp.drop(["index"], axis=1, inplace=True)
    # drop the old index.
    list_lncRNA = df_lncRNA_sample_exp['genesymbol_GRCh38']

    df_lncRNA_sample_exp.drop(["genesymbol_GRCh38", "genetype"], axis=1, inplace=True)

    df_sample_lncRNA_exp = df_lncRNA_sample_exp.T
    # correlation is computed between columns.
    df_lncRNA_lncRNA_corr = df_sample_lncRNA_exp.corr(method='pearson')

    df_lncRNA_lncRNA_corr.columns = list_lncRNA
    # This replaces column names (but not column content) with lncRNA names in list_lncRNA (i.e., re-index columns).
    df_gene_sample_exp2 = df_gene_sample_exp.copy(deep=True)
    # df_gene_sample_exp2 is a temporary dataframe
    list_gene = df_gene_sample_exp2['genesymbol_GRCh38']
    list_type = df_gene_sample_exp2['genetype']

    df_gene_sample_exp2.drop(["genesymbol_GRCh38", "genetype"], axis=1, inplace=True)

    df_sample_gene_exp = df_gene_sample_exp2.T
    print("\nafter transferred\n", df_sample_gene_exp)

    if arg_Pearson == 'Spearman':
        df_gene_gene_corr = df_sample_gene_exp.corr(method='spearman')
    else:
        df_gene_gene_corr = df_sample_gene_exp.corr(method='pearson')
    print(df_gene_gene_corr)

    df_gene_gene_corr.columns = list_gene
    # Using gene names in list_gene to replace the 0-N column names.
    print("\nafter putting column names\n", df_gene_gene_corr)

    for j in range(len(list_gene)):
        if list_type[j] != 'lncRNA':
            df_gene_gene_corr.drop(columns=list_gene[j], inplace=True)
    print("\nafter deleting non-lncRNA columns\n", df_gene_gene_corr)

    df_lncRNA_gene_corr = df_gene_gene_corr.copy(deep=True)
    df_lncRNA_gene_corr.index = list_gene

    return df_lncRNA_lncRNA_corr, df_lncRNA_gene_corr, list_lncRNA

def generate_lncRNA_sets(list_lncRNA, df_lncRNA_lncRNA_corr):
    dict_lncRNAsets = {}
    # For each lncRNA in list_lncRNA, a dict is generated with the lncRNA as the key and a list as the value.
    # The list contains all-mutually-correlated lncRNAs of this lncRNA.
    # Process lncRNAs in df_lncRNA_lncRNA_corr row by row (lncRNAs in list_lncRNA and df_lncRNA_lncRNA_corr have the same order).
    for row in df_lncRNA_lncRNA_corr.index:
        deleted = []
        for col in range(0, len(list_lncRNA)):
            if df_lncRNA_lncRNA_corr.iat[row, col] < arg_LncRNACorrCutoff:
                df_lncRNA_lncRNA_corr.iat[row, col] = -100  # -100 is a number for marking the current column
                df_lncRNA_lncRNA_corr.iat[col, row] = -100  # -100 is a number for marking the corresponding row（diagonal elements are equal）
                deleted.append(list_lncRNA[col])            # Put this lncRNA to deleted[]

        # For lncRNAs with high pairwise corr with the current lncRNA (row), the two-for-loop check/ensure all-mutually-correlated lncRNAs.
        for i in range(0, len(list_lncRNA)):
            if i != row and list_lncRNA[i] not in deleted:
                for j in range(0, len(list_lncRNA)):
                    if list_lncRNA[j] not in deleted:
                        if df_lncRNA_lncRNA_corr.iat[i, j] < arg_LncRNACorrCutoff:
                            df_lncRNA_lncRNA_corr.iat[i, j] = -100  # -100 for marking the newly identified invalid lncRNA
                            deleted.append(list_lncRNA[j])          # Put the newly identified invalid lncRNA into deleted[]
                    else:
                        df_lncRNA_lncRNA_corr.iat[i, j] = -100      # If the col name is already in deleted[], we still set the value to -100.

        lncRNAset = set(list_lncRNA) - set(deleted)
        if list_lncRNA[row] not in deleted:
            lncRNAlist = list(lncRNAset)
        else:
            lncRNAlist = []
        dict_lncRNAsets.update({list_lncRNA[row]: lncRNAlist})
        print(row, {list_lncRNA[row]: lncRNAlist})

    return dict_lncRNAsets

########################################################################################################################
# step 4 - computing TF-gene correlations, TF-TF correlations, and upon which, TFsets.
def generate_TF_files(df_gene_sample_exp):
    df_TF_sample_exp = df_gene_sample_exp.copy(deep=True)
    # When deep=True (default), a new object is created with a copy of the calling object’s data and indices.
    df_TF_sample_exp = df_TF_sample_exp.drop(df_TF_sample_exp[df_TF_sample_exp['genetype'] != 'TF'].index)
    # drop non-TF rows.
    df_TF_sample_exp = df_TF_sample_exp.reset_index()
    # "drop" operation does change the original index.
    # This .reset_index() re-sets index and also keeps the original index in df_gene_sample_exp.
    df_TF_sample_exp.drop(["index"], axis=1, inplace=True)
    # drop the old index.
    list_TF = df_TF_sample_exp['genesymbol_GRCh38']
    # get all TFs in the first column that are gene names (this list also obtains the column name "genesymbol_GRCh38")
    df_TF_sample_exp.drop(["genesymbol_GRCh38", "genetype"], axis=1, inplace=True)
    # drop "genesymbol_GRCh38" and "genetype", because all fileds should be numeric to generate the corr dataframe.
    df_sample_TF_exp = df_TF_sample_exp.T
    # correlation is computed between columns.
    df_TF_TF_corr = df_sample_TF_exp.corr(method='pearson')

    df_TF_TF_corr.columns = list_TF
    # This replaces column names (but not column content) with TF names in list_TF (i.e., re-index columns).
    df_gene_sample_exp2 = df_gene_sample_exp.copy(deep=True)
    # df_gene_sample_exp2 is a temporary dataframe
    list_gene = df_gene_sample_exp2['genesymbol_GRCh38']
    list_type = df_gene_sample_exp2['genetype']

    df_gene_sample_exp2.drop(["genesymbol_GRCh38", "genetype"], axis=1, inplace=True)

    df_sample_gene_exp = df_gene_sample_exp2.T
    print("\nafter transferred\n", df_sample_gene_exp)

    if arg_Pearson == 'Spearman':
        df_gene_gene_corr = df_sample_gene_exp.corr(method='spearman')
    else:
        df_gene_gene_corr = df_sample_gene_exp.corr(method='pearson')
    print(df_gene_gene_corr)

    df_gene_gene_corr.columns = list_gene
    # Using gene names in list_gene to replace the 0-N column names.
    print("\nafter putting column names\n", df_gene_gene_corr)

    for j in range(len(list_gene)):
        if list_type[j] != 'TF':
            df_gene_gene_corr.drop(columns=list_gene[j], inplace=True)
    print("\nafter deleting non-TF columns\n", df_gene_gene_corr)

    df_TF_gene_corr = df_gene_gene_corr.copy(deep=True)
    df_TF_gene_corr.index = list_gene

    return df_TF_gene_corr, df_TF_TF_corr, list_TF

def generate_TF_sets(list_TF, df_TF_TF_corr):
    dict_TFsets = {}
    # For each TF in list_TF, a dict is generated with the TF as key and a list as the value.
    # The list contains all-mutually-correlated TFs of this TF.
    # Process TFs in df_TF_TF_corr row by row. Note TFs in list_TF and df_TF_TF_corr have the same order.
    for row in df_TF_TF_corr.index:  # This row-for-loop check all TFs (rows)
        deleted = []
        for col in range(0, len(list_TF)):
            if df_TF_TF_corr.iat[row, col] < arg_TFCorrCutoff:
                df_TF_TF_corr.iat[row, col] = -100
                df_TF_TF_corr.iat[col, row] = -100
                deleted.append(list_TF[col])

        # For TFs with high pairwise corr with the current TF (row), the two-for-loop check/ensure all-mutually-correlated TFs.
        for i in range(0, len(list_TF)):
            if i != row and list_TF[i] not in deleted:
                for j in range(0, len(list_TF)):
                    if list_TF[j] not in deleted:
                        if df_TF_TF_corr.iat[i, j] < arg_TFCorrCutoff:
                            df_TF_TF_corr.iat[i, j] = -100
                            deleted.append(list_TF[j])
                    else:
                        df_TF_TF_corr.iat[i, j] = -100

        TFset = set(list_TF) - set(deleted)
        if list_TF[row] not in deleted:
            TFlist = list(TFset)
        else:
            TFlist = []
        dict_TFsets.update({list_TF[row]: TFlist})
        print(row, {list_TF[row]: TFlist})

    return dict_TFsets

########################################################################################################################
# step 5 - generate shared lncRNA target sets upon df_lncRNA_gene_DBS and df_lncRNA_gene_corr
def generate_lncRNA_target(df_lncRNA_gene_DBS, df_lncRNA_gene_corr, dict_lncRNAsets):
    '''
    Upon df_lncRNA_gene_corr and df_lncRNA_gene_DBS, we get targets of every lncRNA set, which are the
    intersection of set_lncRNAtarget_DBS (DBS>cutoff1) and set_lncRNAtarget_CORR (|corr|>cutoff2).
    cutoff2 (lncRNA-gene-corr) = lncRNA-lncRNA-corr (arg_LncRNACorrCutoff). A separable cutoff can be considered.
    :
    traverse key in dict_lncRNAsets
       traverse lncRNA in lncRNAlist
          get targets of lncRNA from df_lncRNA_gene_DBS,  append targets to dict_lncRNAtarget_DBS
          get targets of lncRNA from df_lncRNA_gene_CORR, append targets to dict_lncRNAtarget_CORR
       get overlap of the two target sets, append the overlap to dict_lncRNAtarget_shared
    '''
    list_lncRNACORR = df_lncRNA_gene_corr.columns
    list_lncRNADBS = df_lncRNA_gene_DBS.columns

    dict_lncRNAtarget_DBS = {}
    dict_lncRNAtarget_CORR = {}
    dict_lncRNAtarget_shared = {}

    rows1 = df_lncRNA_gene_corr.shape[0]
    cols1 = df_lncRNA_gene_corr.shape[1]
    rows2 = df_lncRNA_gene_DBS.shape[0]
    cols2 = df_lncRNA_gene_DBS.shape[1]

    for key, value in dict_lncRNAsets.items():
        set_tpmCORR = set()
        set_tmpDBS = set()
        if key in list_lncRNACORR and value != []:  # if the COCO set of this lncRNA is non-empty
            for lncRNA in value:                    # then process all lncRNAs in the lncRNAlist
                if lncRNA in df_lncRNA_gene_corr.columns:
                    for target in df_lncRNA_gene_corr.index:
                        if abs(df_lncRNA_gene_corr[lncRNA][target]) > arg_LncRNACorrCutoff:
                            set_tpmCORR.add(target)
        dict_lncRNAtarget_CORR.update({key: set_tpmCORR})

        if key in list_lncRNADBS and value != []:
            for lncRNA in value:
                if lncRNA in df_lncRNA_gene_DBS.columns:
                    for target in df_lncRNA_gene_DBS.index:
                        if df_lncRNA_gene_DBS[lncRNA][target] > arg_LncRNADBSCutoff:
                            set_tmpDBS.add(target)
        dict_lncRNAtarget_DBS.update({key: set_tmpDBS})

    # dict_lncRNAtarget_CORR and dict_lncRNAtarget_DBS have the same keys and key order.
    for key, value1 in dict_lncRNAtarget_DBS.items():
        if key in dict_lncRNAtarget_CORR:
            value2 = dict_lncRNAtarget_CORR[key]
            shared = value1 & value2
            dict_lncRNAtarget_shared.update({key: shared})

    return dict_lncRNAtarget_shared

########################################################################################################################
# step 6 - generate shared TF target sets upon df_TF_gene_DBS and df_TF_gene_corr
def generate_TF_target(df_TF_gene_DBS, df_TF_gene_corr, dict_TFsets):
    list_TFCORR = df_TF_gene_corr.columns
    list_TFDBS  = df_TF_gene_DBS.columns

    dict_TFtarget_DBS = {}
    dict_TFtarget_CORR = {}
    dict_TFtarget_shared = {}

    rows1 = df_TF_gene_corr.shape[0]
    cols1 = df_TF_gene_corr.shape[1]
    rows2 = df_TF_gene_DBS.shape[0]
    cols2 = df_TF_gene_DBS.shape[1]

    for key, value in dict_TFsets.items():
        set_tmpCORR = set()
        set_tmpDBS = set()
        if key in list_TFCORR and value != []:
            for TF in value:
                if TF in df_TF_gene_corr.columns:
                    for target in df_TF_gene_corr.index:
                        if abs(df_TF_gene_corr[TF][target]) > arg_TFCorrCutoff:
                            set_tmpCORR.add(target)
        dict_TFtarget_CORR.update({key: set_tmpCORR})

        if key in list_TFDBS and value != []:
            for TF in value:
                if TF in df_TF_gene_DBS.columns:
                    for target in df_TF_gene_DBS.index:
                        if df_TF_gene_DBS[TF][target] > arg_TFDBSCutoff:
                            set_tmpDBS.add(target)
        dict_TFtarget_DBS.update({key: set_tmpDBS})

    # dict_TFtarget_CORR and dict_TFtarget_DBS have the same keys and key order.
    for key, value1 in dict_TFtarget_DBS.items():
        if key in dict_TFtarget_CORR:
            value2 = dict_TFtarget_CORR[key]
            shared = value1 & value2
            dict_TFtarget_shared.update({key: shared})

    return dict_TFtarget_shared

########################################################################################################################
# step 7 - identify modules and merge lncRNAs/TFs that regulate the same modules
def collect_modules(dict_lncRNAtarget_shared, dict_TFtarget_shared):
    '''
    This function collects all lncRNAs' modules (lncRNA target sets) and TFs' modules (TF target sets) into a dict,
    with target genes in the module (a tuple) as key and lncRNAs/TFs that co-regulate the module as values.
    '''
    dict_module_regulators = {} #

    # zhu
    # dict_lncRNAtarget_shared = {'AC090152.1': {'LRMDA', 'FGFR2', 'ABLIM1', 'TNFSF9'}, 'LINC02783': set(), 'Z99572.1': set(), 'AC010889.1': {'GABRA3', 'ABLIM1'}, 'ZNF790-AS1': {'COL18A1-AS1', 'MED14'}, 'AC125807.2': set()}
    # dict_TFtarget_shared     = {'ZNF790': {'GABRA3', 'FGFR2', 'CEACAM6', 'TNFSF9'}, 'SOX21': {'LINC02577', 'COL18A1-AS1', 'MED14'}, 'ALX1': {'DENND2A', 'TSPAN7'}, 'HNF4A': {'DENND2A', 'TSPAN7'}, 'NR5A2': set()}
    #dict_xx_shared = {'1': {'1', '2', '3', '4'}, '5': set(), '6': set(), '7': {'8', '9'}, '10': {'11', '12'}, '13': set()}
    #dict_zz_shared = {'a': {'b', 'd', 'e', 'f'}, 'g': {'h', 'i', 'j'}, 'k': {'l', 'm'}, 'n': {'o', 'p'}, 'q': set()}
    #allModules = [tuple(sorted(module)) for module in list(dict_xx_shared.values()) + list(dict_zz_shared.values()) if len(module) > 1]

    # We sort the order of target genes in modules so that we can conveniently compare and determine whether two modules are the same
    allModules = [tuple(sorted(module)) for module in list(dict_lncRNAtarget_shared.values()) + list(dict_TFtarget_shared.values()) if len(module) > arg_ModuleSize]
    for module in allModules:
        dict_module_regulators[module] = set()

    for lncRNA in dict_lncRNAtarget_shared:
        module = tuple(sorted(dict_lncRNAtarget_shared[lncRNA]))
        # dict_lncRNAtarget_shared has 102 elements, each has a key (a lncRNA) and a value (the set of lncRNAs corresponding to the key)
        # module is the key, now the key of lncRNA set is the value, better to use the value of lncRNA set as the value
        if module in dict_module_regulators:
            dict_module_regulators[module].add(lncRNA + '(*)')  # we use * to indicate lncRNA

    for TF in dict_TFtarget_shared:
        module = tuple(sorted(dict_TFtarget_shared[TF]))
        if module in dict_module_regulators:
            dict_module_regulators[module].add(TF + '(#)')  # we use # to indicate TF

    print("dict_module_regulators=", dict_module_regulators)
    return dict_module_regulators

########################################################################################################################
# step 8 - perform pathway enrichment analysis for module genes ### in df_path_module
def pathway_analysis(dict_module_regulators, kegg_gene, kegg_pathway, kegg_link, totalGeneNum, arg_fdrCutoff):
    dict_module_kegg = {}
    results = []
    pvalues = []
    
    for module in dict_module_regulators:
        module_geneIDs = {}  
        #a gene may have serveral symbols
        #we convert its symbol which user assign to KEGG geneID to obtain the interset between the module and a pathway
        for gene in module:
            if gene in kegg_gene:
                geneID = kegg_gene[gene]
                module_geneIDs[geneID] = gene
                
        moduleSize = len(module_geneIDs)
        
        for pathID in kegg_link:
            interSet = set(module_geneIDs.keys()) & kegg_link[pathID]
            if interSet:
                hitNum = len(interSet)
                thisPathGeneNum = len(kegg_link[pathID])
                
                # Indentify significantly enriched pathways using hypergeometric distribution test
                pVal = hypergeom.sf(hitNum-1, totalGeneNum, thisPathGeneNum, moduleSize)
                if pVal<0.05:
                    pathName = kegg_pathway[pathID]
                    hitGenes = set([module_geneIDs[geneID] for geneID in interSet])
                    results.append([module, pathID, pathName, hitGenes, pVal])
                    pvalues.append(pVal)
                
    #Benjamini-Hochberg method for p-value correction
    pValues = np.asfarray(pvalues)
    by_descend = pValues.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(pValues)) / np.arange(len(pValues), 0, -1)
    fdr_values = np.minimum(1, np.minimum.accumulate(steps * pValues[by_descend]))
    fdrs = []
    for i in range(len(fdr_values[by_orig])):
        fdrs.append(fdr_values[by_orig][i])
        
    for i in range(len(fdrs)):
        results[i].append(fdrs[i])
        
    for module, pathID, pathName, hitGenes, pVal, fdr in results:
        if fdr < arg_fdrCutoff:
            if module not in dict_module_kegg:
                dict_module_kegg[module] = [(fdr, pathID, pathName, hitGenes)]
            else:
                dict_module_kegg[module].append((fdr, pathID, pathName, hitGenes))

    return dict_module_kegg
    
########################################################################################################################
# step 9 - write eGRAM results to files in the directory "eGRAM_results/"
def write_modules(dict_module_regulators, dict_module_kegg, dict_module_wiki, outFile):
    with open('eGRAM_results/' + outFile, 'w') as f:
        f.write('\t'.join(['ModuleID',
                           'LncRNA(*)/TF(#)',
                           'Target gene',
                           'KEGG pathway, pathwayID, and FDR',
                           'Hit genes of top one significant KEGG pathway',
                           'WikiPathway, pathwayID, and FDR',
                           'Hit genes of top one significant Wikipathway'])+'\n')
        moduleID = 0
        for module in dict_module_regulators:
            moduleID += 1
            if module in dict_module_kegg:
                this_module_kegg = []
                this_module_kegg_hitGenes = set()
                for fdr, pathID, pathName, hitGenes in sorted(dict_module_kegg[module]):
                    this_module_kegg.append(pathName + ' (' + pathID + ', ' + str('{0:3f}'.format(fdr)) + ')')
                    this_module_kegg_hitGenes = this_module_kegg_hitGenes | hitGenes
                pathway_kegg = '; '.join(this_module_kegg)
                hitGenes_kegg = ', '.join(this_module_kegg_hitGenes)
            else:
                pathway_kegg = hitGenes_kegg = ''
               
            if module in dict_module_wiki:
                this_module_wiki = []
                this_module_wiki_hitGenes = set()
                for fdr, pathID, pathName, hitGenes in sorted(dict_module_wiki[module]):
                    #this_module_wiki.append(pathName + ' (' + pathID + ', ' + str(fdr) + ')')
                    this_module_wiki.append(pathName + ' (' + pathID + ', ' + str('{0:3f}'.format(fdr)) + ')')
                    this_module_wiki_hitGenes = this_module_wiki_hitGenes | hitGenes
                pathway_wiki = '; '.join(this_module_wiki)
                hitGenes_wiki =  ', '.join(this_module_wiki_hitGenes)
            else:
                pathway_wiki = hitGenes_wiki = ''
               
            regulators = ' | '.join(dict_module_regulators[module])
            targetGenes = ' | '.join(module)
            
            f.write('\t'.join(['ME'+str(moduleID),
                               regulators,
                               targetGenes,
                               pathway_kegg,
                               hitGenes_kegg,
                               pathway_wiki,
                               hitGenes_wiki])+'\n')

# Generate and write files for module visulization using cytoscape in the 'eGRAM_results/' directory.
def generate_cytosc_file(moduleFile):
    if '/' not in moduleFile:
        moduleFile = 'eGRAM_results/' + moduleFile
    fileName = moduleFile.split('/')[-1]
    f_re1 = open('eGRAM_results/' + fileName + '_edge', 'w')
    f_re1.write('\t'.join(['lncRNA', 'module', 'weight'])+'\n')
    f_re2 = open('eGRAM_results/' + fileName + '_group', 'w')
    f_re2.write('\t'.join(['module', 'Keggpathway', 'Wikipathway'])+'\n')
    allLncs = set()
    allTFs = set()
    with open(moduleFile) as f_mod:
        f_mod.readline()
        for line in f_mod:
            cols = line.strip('\n').split('\t')
            moduleID = cols[0]
            regulators = cols[1].split(' | ')
            
            kegg_results = cols[3].split('; ')
            kegg_pathways = []
            for kegg_result in kegg_results:
                kegg_pathways.append(kegg_result[:kegg_result.find(' (')])
            kegg_pathways_str = '; '.join(kegg_pathways)
            
            wiki_results = cols[5].split('; ')
            wiki_pathways = []
            for wiki_result in wiki_results:
                wiki_pathways.append(wiki_result[:wiki_result.find(' (')])
            wiki_pathways_str = '; '.join(wiki_pathways)
            
            for regulator in regulators:
                regulator_symbol = regulator.split('(')[0]
                f_re1 .write('\t'.join([regulator_symbol, moduleID, 'NA']) + '\n')
                if '(*)' in regulator:
                    allLncs.add(regulator_symbol)
                elif '(#)' in regulator:
                    allTFs.add(regulator_symbol)
            f_re2.write('\t'.join([moduleID, kegg_pathways_str, wiki_pathways_str]) + '\n')
        for lnc in allLncs:
            f_re2.write('\t'.join([lnc, 'LncRNA']) + '\n')
        for TF in allTFs:
            f_re2.write('\t'.join([TF, 'TF']) + '\n')
    f_re1.close()
    f_re2.close()

# Modify the resulting file name to prevent overwriting existing files
def check_output(arg_OutFile, resultDir):
    maxFileNum = -1
    for file in os.listdir(resultDir):
        if file==arg_OutFile and maxFileNum==-1:
            maxFileNum = 0
        elif arg_OutFile+'_' in file:
            this_fileNum = file.replace(arg_OutFile+'_', '')
            if this_fileNum.isdigit():
                this_fileNum = int(this_fileNum)
                if this_fileNum > maxFileNum:
                    maxFileNum = this_fileNum
    if maxFileNum == -1:
        return
    else:
        newFileNum = maxFileNum + 1
        newOutFile = arg_OutFile + '_' + str(newFileNum)
        return newOutFile


########################################################################################################################
def usage():
    print('Usage: python eGRAMv2R1.py --exp gene_exp_file --lncDBS lncRNA_DBS_file --tfDBS TF_DBS_file --species 1 --lncCutoff lncRNA-DBS-cutoff --tfCutoff TF-DBS-cutoff --lncCorr lncRNA-corr-cutoff --tfCorr TF-corr-cutoff --moduleSize modulesize --corr Pearson/Spearman --cluster FLAME/CoCo --fdr 0.01 --out output')
    print()
    print('Example: python eGRAMv2R1.py --exp 2024May-DEG-exp-A549-2513WT.csv --lncDBS 2024May-lncRNA-DBS-A549-2513.csv --tfDBS 2024May-TF-DBS-A549-2513.csv --species 1 --lncCutoff 100 --tfCutoff 8 --lncCorr 0.5 --tfCorr 0.5 --moduleSize 50 --corr Pearson --cluster COCO --fdr 0.01 --out myout')
    print()
    print("""Parameters:
        --exp gene_exp_file is required; if genetype does not have TF, all TF-related functions are not performed. 
        --lncDBS     lncRNA_DBS_file is required; integer or 0. 
        --tfDBS      TF_DBS_file is not a must; required if gene_exp_file contains TF.
        --species    1=human, 2=mouse.
        --lncCutoff  the DBS threshold for determining lncRNAs' targets (defult=100).
        --tfCutoff   the DBS threshold for determining TFs's targets (default=10).
        --lncCorr    the correlation threshold for lncRNA/target and lncRNA/lncRNA (default=0.5).
        --tfCorr     the correlation threshold for TF/target and TF/TF (default=0.5).
        --moduleSize the minimum gene number of a module (default=50).
        --corr       Pearson/Spearman (default=pearson).
        --cluster    clustering (default=COCO, see 2018Saelons).
        --fdr        the FDR value to determine significance of pathway enrichment (default 0.01).
        --out        output file name.
        """)
########################################################################################################################
########################################################################################################################
if __name__ == '__main__':
    param = sys.argv[1:]

    try:
        opts, args = getopt.getopt(param, '-h', ['exp=', 'lncDBS=', 'tfDBS=', 'species=', 'lncCutoff=', 'tfCutoff=', 'lncCorr=', 'tfCorr=', 'moduleSize=', 'corr=' , 'cluster=', 'fdr=', 'out='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    # step 1 - read in user files and command line arguments
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit(2)
        elif opt == '--exp':
            df_gene_sample_exp = pd.read_csv(arg)
        elif opt == '--lncDBS':
            df_lncRNA_gene_DBS = pd.read_csv(arg, index_col = [0])
        elif opt == '--tfDBS':
            df_TF_gene_DBS = pd.read_csv(arg, index_col = [0])
        elif opt == '--species':
            arg_SpeciesPath = int(arg)
        elif opt == '--lncCutoff':
            arg_LncRNADBSCutoff = int(arg)
        elif opt == '--tfCutoff':
            arg_TFDBSCutoff = int(arg)
        elif opt == '--lncCorr':
            arg_LncRNACorrCutoff = float(arg)
        elif opt == '--tfCorr':
            arg_TFCorrCutoff = float(arg)
        elif opt == '--moduleSize':
            arg_ModuleSize = int(arg)
        elif opt == '--corr':
            arg_Pearson = arg
        elif opt == '--cluster':
            arg_Clustering = arg
        elif opt == '--fdr':
            arg_fdrCutoff = float(arg)
        elif opt == '--out':
            arg_OutFile = arg

    # default arguments
    if 'arg_SpeciesPath' not in dir():
        arg_SpeciesPath = 1
    if 'arg_LncRNADBSCutoff' not in dir():
        arg_LncRNADBSCutoff = 100
    if 'arg_TFDBSCutoff' not in dir():
        arg_TFDBSCutoff = 10
    if 'arg_LncRNACorrCutoff' not in dir():
        arg_LncRNACorrCutoff = 0.5
    if 'arg_TFCorrCutoff' not in dir():
        arg_TFCorrCutoff = 0.5
    if 'arg_ModuleSize' not in dir():
        arg_ModuleSize = 50
    if 'arg_Pearson' not in dir():
        arg_Pearson = 'Pearson'
    if 'arg_Clustering' not in dir():
        arg_Clustering = 'COCO'
    if 'arg_fdrCutoff' not in dir():
        arg_fdrCutoff = 0.01
    if arg_SpeciesPath == 1:
        keggGeneFile = "KEGG/keggGene-hsa"
        keggPathFile = "KEGG/keggPathway-hsa"
        keggLinkFile = "KEGG/keggLink-hsa"
        wikiPath = 'WikiPathways/wikipathways-20240310-gpml-Homo_sapiens/'
    elif arg_SpeciesPath == 2:
        keggGeneFile = "KEGG/keggGene-mmu"
        keggPathFile = "KEGG/keggPathway-mmu"
        keggLinkFile = "KEGG/keggLink-mmu"
        wikiPath = 'WikiPathways/wikipathways-20240310-gpml-Mus_musculus/'
    else:
        sys.exit(2)

########################################################################################################################
    # step 2 - read in kegg files and wikipathway files
    kegg_gene, kegg_pathway, kegg_link, totalGeneNum = read_kegg(keggGeneFile, keggPathFile, keggLinkFile)

    wiki_pathway, wiki_link = read_wikip(wikiPath, kegg_gene)
    
    resultDir = 'eGRAM_results/'
    if not os.path.exists(resultDir):
        os.mkdir(resultDir)

    # step 3 - computing lncRNA-gene correlations, lncRNA-lncRNA correlations, and upon which, lncRNAsets
    # i.e., clustering lncRNAs into 'COCO' sets (this could be done by introduing an outside program such as FLAME)
    df_lncRNA_lncRNA_corr, df_lncRNA_gene_corr, list_lncRNA = generate_lncRNA_files(df_gene_sample_exp)
    dict_lncRNAsets = generate_lncRNA_sets(list_lncRNA, df_lncRNA_lncRNA_corr)

    # step 4 - computing TF-gene correlations, TF-TF correlations, and upon which, TFsets.
    # This is the counterpart of the step 3 for TFs.
    df_TF_gene_corr, df_TF_TF_corr, list_TF = generate_TF_files(df_gene_sample_exp)
    dict_TFsets = generate_TF_sets(list_TF, df_TF_TF_corr)

    # step 5 - generate shared lncRNA target sets upon df_lncRNA_gene_DBS and df_lncRNA_gene_corr
    dict_lncRNAtarget_shared = generate_lncRNA_target(df_lncRNA_gene_DBS, df_lncRNA_gene_corr, dict_lncRNAsets)

    # step 6 - generate shared TF target sets upon df_TF_gene_DBS and df_TF_gene_corr
    if 'df_TF_gene_DBS' in dir():
        dict_TFtarget_shared = generate_TF_target(df_TF_gene_DBS, df_TF_gene_corr, dict_TFsets)
    else:
        dict_TFtarget_shared = {}

    # step 7 - identify modules and merge lncRNAs/TFs that regulate the same modules
    dict_module_regulators = collect_modules(dict_lncRNAtarget_shared, dict_TFtarget_shared)
    
    # step 8 - perform pathway enrichment analysis for module genes in df_path_module. KEGG gene annotation is used.
    dict_module_kegg = pathway_analysis(dict_module_regulators, kegg_gene, kegg_pathway, kegg_link, totalGeneNum, arg_fdrCutoff)
    dict_module_wiki = pathway_analysis(dict_module_regulators, kegg_gene, wiki_pathway, wiki_link, totalGeneNum, arg_fdrCutoff)
    
    # step 9 - write results to files
    newOutFile = check_output(arg_OutFile, resultDir)
    if newOutFile:
        arg_OutFile = newOutFile
        print('A result file with the same name already exists, the newly generated file is named with '+ newOutFile + '.')
    write_modules(dict_module_regulators, dict_module_kegg, dict_module_wiki, arg_OutFile)
    
    generate_cytosc_file(arg_OutFile)