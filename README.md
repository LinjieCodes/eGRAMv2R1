# A tool for identifying gene modules and regulatory lncRNAs and TFs

The eGRAMv2R1 program identifies gene modules comprising co-expressed genes, their regulatory lncRNAs and their regulatory TFs based on lncRNA/DNA bindings, TF/DNA bindings and gene expression correlations.

There are two kinds of input files - (a) user data and (b) system data. The former has three files - (a) expression profile, (b) lncRNA target prediction, (c) TF target prediction. LncRNAs' targets are predicted using LongTarget (or the more rapid version Fasim); TFs' targets are predicted using any method, such as CellOracle. Since the primary goal of eGRAM is analyzing lncRNA function in transcriptional regulation, the TF target prediction file is optional. System data include human/mouse KEGG and wikipathway pathways downloaded from the KEGG and Wikipathways websites.

# Requirements
1. **Python**: >=3.7.0

2. **numpy**: >=1.21.6

3. **pandas**: >=1.3.5

4. **scipy**: >=1.7.3

5. **OS**: the eGRAMv2R1 code has been tested on Linux system.

# Data
1. **2024May-DEG-exp-A549-2513WT.csv**  --  the gene expression matrix.

2. **2024May-lncRNA-DBS-A549-2513.csv**  --  the DNA binding matrix (i.e. target genes) of lncRNAs.

3. **2024May-TF-DBS-A549-2513.csv**  --  the DNA binding matrix (i.e. target genes) of TFs.

4. **KEGG**  --  the KEGG pathway annotation of human and mouse.

5. **WikiPathways**  --  the WikiPathways annotation of human and mouse. The compressed package needs to be decompressed before running the program.

# Usage
Here is a command line to run the eGRAMv2R1 program:

```
'Example: python eGRAMv2R1.py --exp 2024May-DEG-exp-A549-2513WT.csv --lncDBS 2024May-lncRNA-DBS-A549-2513.csv --tfDBS 2024May-TF-DBS-A549-2513.csv --species 1 --lncCutoff 100 --tfCutoff 8 --lncCorr 0.5 --tfCorr 0.5 --moduleSize 50 --corr Pearson --cluster COCO --fdr 0.01 --out myout'
```

# Help information
Here is a brief explanation of the command line arguments:

```
Parameters      Functions
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
```

# Bug reports
Please send comments and bug reports to JL.linjie@outlook.com or zhuhao@smu.edu.cn.

# Related website
To obtain details about lncRNA/DNA binding prediction using LongTarget, please go to our website http://www.gaemons.net/LongTarget.
