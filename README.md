

## EasyMultiProfiler TuTorial

### Install

Easily install

```R
if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
remotes::install_github("liubingdong/EasyMultiProfier")
```

Recommend install

```R
setRepositories(addURLs = c(BioCsoft = "https://bioconductor.org/packages/3.18/bioc",
                  BioCann = "https://bioconductor.org/packages/3.18/data/annotation"))  
options(timeout = 600000000) 
install.packages("remotes")
install.packages("BiocManager")
install.packages("dendextend")
install.packages("AnnotationDbi")
remotes::install_github("liubingdong/EasyMultiProfier")
```

### Import data 

1. Complete import model

   ```R
   meta_data <- read.table('col.txt',header = T,row.names = 1)
   
   dfmap <-read.table('dfmap.txt',header = T) |> 
     dplyr::mutate(assay = factor(assay)) ## factor is necessary
   
   tax_data <- EMP_taxonomy_import('tax.txt')
   ko_data <- EMP_function_import('ko.txt',type = 'ko')
   
   ec_data <- EMP_function_import('ec.txt',type = 'ec')
   
   metbol_data <- EMP_normal_import('metabol.txt',
                                    dfmap = dfmap,assay  = 'untarget_metabol')
   
   geno_data <- EMP_normal_import('tran.txt',dfmap = dfmap,assay  = 'host_gene')
   
   #### Note that the naming here must be consistent with that of dfmap
   objlist <- list("taxonomy" = tax_data,
                   "geno_ko" = ko_data,
                   "geno_ec" = ec_data,
                   "untarget_metabol" = metbol_data,
                   "host_gene" = geno_data)
   
   MAE <- MultiAssayExperiment::MultiAssayExperiment(objlist, meta_data, dfmap)
   ```

2. Quickly import model ( Only for one single omics data)

   ```R
   meta_data <- read.table('col.txt',header = T,row.names = 1)
   tran_data <- EMP_easy_import('tran.txt',coldata = meta_data,type = 'normal')
   ```

### Prepare the demo data from EasyMultiProfiler

```R 
Data(MAE)
```

### EMP_assay_extract

1. Extract the assay data of one existed experiment from MultiAssaayExperiment

 ```R
 MAE |>
    EMP_assay_extract('taxonomy')
MAE |>
    EMP_assay_extract('geno_ko')
 ```


2. Search for specific features according to the rowdata

```R
MAE |>
    	EMP_assay_extract('geno_ec',pattern = '1.1.1.1',pattern_ref = 'feature',exact = T)
MAE |>
    	EMP_assay_extract('geno_ko',pattern = 'mtlD',pattern_ref = 'Name',exact = F)
```

### EMP_rowdata_extract

1. Extract the rowdata of one existed experiment from MultiAssaayExperiment

```R
MAE |>
  EMP_rowdata_extract('taxonomy')
MAE |>
  EMP_rowdata_extract('geno_ko')  
```

###### 2. Extract the rowdata of MultiAssaayExperiment

```R
MAE |>
  EMP_rowdata_extract(experiment = NULL) -> total_row_data
dim(total_row_data)
```

### EMP_coldata_extract

1. Extract the coldata/meta-data/sample-info/patient-info of one existed experiment from MultiAssaayExperiment

```R
MAE |>
  EMP_coldata_extract('taxonomy')
```

2. Extract all coldata

```R
MAE |>
  EMP_coldata_extract(action = 'get')  ### when action = get, output is a tibble.
MAE |>
  EMP_coldata_extract(action = 'add')  ### when action = add, output is a EMPT object 
```


### EMP_feature_convert

1. for gene ID convert

```R
MAE |>
  EMP_feature_convert(experiment = 'host_gene',
                      from = 'SYMBOL',to='ENTREZID',species = 'Human')
```

The built-in database only supports Human, Mouse, Pig, Zebrafish.
Other species could utilize OrgDb to convert.
More OrgDb is on the  https://bioconductor.org/packages/release/BiocViews.html#___OrgDb

```R
library(org.Hs.eg.db)
MAE |>
	EMP_feature_convert(experiment = 'host_gene',
                      from = 'SYMBOL',to='ENTREZID',OrgDb = org.Hs.eg.db)
```

2. for compound ID convert

```R
MAE |>
  EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
               estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |> 
  EMP_feature_convert(from = 'KEGG',to='HMDB')
```

### EMP_impute

1. For coldata

```R
MAE |>
  EMP_assay_extract('geno_ec') |>
  EMP_impute(assay=F,coldata=T,rowdata=F)
```

Support formula, such as only impute SAS and SDS

```R
MAE |>
  EMP_coldata_extract(action ='add') |>
  EMP_impute(.formula = SAS+SDS ~ .,action = 'get') 
```

2. For assay

```R
MAE |>
  EMP_coldata_extract(action ='add') |>
  EMP_impute(.formula = SAS+SDS ~ .,action = 'get') 
```

3. For rowdata (Not Not recommended)

```R
MAE |>
  EMP_assay_extract('geno_ec') |>
  EMP_impute(assay=F,coldata=F,rowdata=T)
```




### EMP_identify_assay

1. Consider the minnum relative abundance and min ratio specilally for microbial data. First, the abundance below this threshold will be converted to 0 according to the set minimum species relative abundance. Secondly, the core species will be required to meet the requirement that the occurrence rate of the species group in at least one group is higher than the preset threshold, and the remaining species will Then it will be judged as "rare species" and filtered.
   Note: If absolute abundance is provided as input, it will be automatically converted to relative abundance for filtering purposes during calculations. However, the output will remain in absolute abundance.

```R
MAE |>
  EMP_assay_extract('taxonomy') |>
  EMP_identify_assay(estimate_group = 'Group',method = 'default',
                     min=0.01,min_ratio = 0.7) ### consider the Group
MAE |>
  EMP_assay_extract('taxonomy') |>
  EMP_identify_assay(method = 'default') ### consider all samples belong to one group
```

2. Consider the minnum counts abundance and min ratio specially for expression data

```R
MAE |>
  EMP_assay_extract('geno_ec') |>
  EMP_identify_assay(method = 'edgeR',min = 10,min_ratio = 0.7,estimate_group = 'Group') 
```


### EMP_modify_assay

For some special cases, to modify the assay data and EMP_identify_assay works better in most cases.

1. Change the expression value which is 0 into 0.0001

  ```R
MAE |>
  EMP_assay_extract('geno_ec') |>
  EMP_modify_assay('==0',pseudocount=0.0001)
  ```

2. Change the counts which is below some threshold

  ```R
MAE |>
  EMP_assay_extract('taxonomy') |>
  EMP_modify_assay('<10',pseudocount=0) 
  ```


### EMP_rrarefy

Rarefythe data according to the lowest abundance.

```R
MAE |>
  EMP_rrarefy(experiment = 'taxonomy')
```

Only show the depth

```R
MAE |>
  EMP_rrarefy(experiment = 'taxonomy',only_show_depth=T) 
```

Set a specific threshold

```R
MAE |>
  EMP_rrarefy(experiment = 'taxonomy',raresize=1000) 
```


### EMP_decostand

Transfer data into relative format.

```R
MAE |>
  EMP_decostand(experiment = 'taxonomy',method = 'relative') 
```

Transfer data into centered log ratio ("clr") format.

```R
MAE |>
  EMP_decostand(experiment = 'geno_ko',method = 'clr',pseudocount=0.0001)
```

Transfer data into logformat.

```R
MAE |>
  EMP_decostand(experiment = 'geno_ec',method = 'log',logbase = 2) 
```

### EMP_collapse

merge assay data accoding to duplicate coldata.

```R
MAE |>
  EMP_collapse(experiment = 'untarget_metabol',collapse_by='col',
               estimate_group = 'Group',method = 'mean',collapse_sep = '+') 
```

merge assay data accoding to duplicate rowdata.

```R
MAE |> EMP_rowdata_extract('untarget_metabol')
MAE |>
  EMP_collapse(experiment = 'untarget_metabol',
               collapse_by='row',na_string = c("NA", "null", "","-"),
               estimate_group = 'MS2kegg',method = 'mean',collapse_sep = '+') 
```

combie two collapse method

```R
MAE |>
  EMP_collapse(experiment = 'untarget_metabol',
               collapse_by='row',na_string = c("NA", "null", "","-"),
               estimate_group = 'MS2kegg',method = 'mean',collapse_sep = '+') |>
  EMP_collapse(collapse_by='col',estimate_group = 'Group',
               method = 'mean',collapse_sep = '+')
```

### EMP_adjust_abudance

combat_seq method 

```R
MAE |>
  EMP_collapse(experiment = 'untarget_metabol',
               collapse_by='row',na_string = c("NA", "null", "","-"),
               estimate_group = 'MS2kegg',method = 'mean',collapse_sep = '+') |>
  EMP_collapse(collapse_by='col',estimate_group = 'Group',method = 'mean',collapse_sep = '+')
```

combat method 

```R
MAE |>
  EMP_assay_extract(experiment='geno_ko') |>
  EMP_adjust_abudance(.factor_unwanted = 'Region',.factor_of_interest = 'Group',
                      method = 'combat',action = 'add') 
```

limma_remove_batch_effect

```R
MAE |>
  EMP_assay_extract(experiment='geno_ko') |>
  EMP_adjust_abudance(.factor_unwanted = 'Region',.factor_of_interest = 'Group',
                      method = 'limma_remove_batch_effect') 
```



### EMP_alpha_analysis

```R
MAE |>
  EMP_assay_extract(experiment='taxonomy') |> 
  EMP_alpha_analysis()
MAE |>
  EMP_assay_extract(experiment='geno_ec') |> 
  EMP_alpha_analysis()
```

### EMP_cluster_analysis

1. Cluster the samples according to the assay data

  ```R
MAE |>
  EMP_assay_extract(experiment = 'geno_ec') |>
  EMP_cluster_analysis()
MAE |>
  EMP_assay_extract(experiment = 'geno_ec') |>
  EMP_cluster_analysis(h=0.15) |> ### identify the outlier samples
  EMP_filter(cluster != 1) ### filter away the outlier samples
  ```

2. Cluster the samples according to the coldata

  ```R
MAE |> 
  EMP_coldata_extract(action = 'add') |> ### transfer the coldata to asaay
  EMP_impute(assay = T) |> ### impute missing value
  EMP_cluster_analysis(method = 'ward.D2',distance='clark',h=0.2) 
  ```

3. Cluster the features according to the assay data

  ```R
MAE |> 
  EMP_assay_extract(experiment = 'geno_ec',
                    pattern='1.1.1.1',pattern_ref='feature') |>
  EMP_cluster_analysis(rowdata = T,h=0.8)
  ```

### EMP_cor_analysis

```R
k1 <- MAE |>
  EMP_assay_extract('taxonomy') |>
  EMP_collapse(estimate_group = 'Genus',collapse_by = 'row') |>
  EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
  EMP_filter(feature_condition = pvalue<0.05)
k2 <- MAE |>
  EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
               estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |>
  EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
  EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) > 1.5)

(k1 + k2) |> EMP_cor_analysis(method = 'spearman') |>
  EMP_heatmap_plot() ###### Visualization
```


### EMP_diff_analysis

t.test or wilcox.test

```R
MAE |>
  EMP_decostand(experiment = 'taxonomy',method = 'relative',pseudocount=0.0001) |>
  EMP_diff_analysis(method = 't.test',estimate_group = 'Group',p.adjust = 'fdr')
```

```R
MAE |>
  EMP_decostand(experiment = 'taxonomy',method = 'relative',pseudocount=0.0001) |>
  EMP_diff_analysis(method = 't.test',estimate_group = 'Group',p.adjust = 'fdr')
```

DESeq2

```R
MAE |>
  EMP_decostand(experiment = 'geno_ec',method = 'integer') |>
  EMP_diff_analysis(method='DESeq2',.formula = ~Group)  

MAE |>
  EMP_decostand(experiment = 'geno_ec',method = 'integer') |>
  EMP_diff_analysis(method='DESeq2',
                    .formula = ~Region+Group)  # Eliminate the batch_effect in DESeq2
```

edgeR_quasi_likelihood

```R
MAE |>
  EMP_assay_extract(experiment = 'geno_ec') |>
  EMP_diff_analysis(method='edgeR_quasi_likelihood',
                    .formula = ~0+Group,
                    estimate_group = c('Group_B','Group_A')) # Set the comparison order.
```

### EMP_dimension_analysis

PCA

```R
MAE |>
  EMP_dimension_analysis(experiment = 'taxonomy',method = 'pca') 
```

Pcoa

```R
MAE |>
  EMP_dimension_analysis(experiment = 'taxonomy',method = 'pcoa',distance = 'bray')
```

Pls

```R
MAE |>
  EMP_dimension_analysis(experiment = 'untarget_metabol',
                         method = 'pls',estimate_group = 'Group')
MAE |>
  EMP_collapse(experiment = 'untarget_metabol',
               na_string=c('NA','null','','-'),
               estimate_group = 'MS2kegg',
               method = 'sum',collapse_by = 'row') |> # Get the kegg compound
  EMP_dimension_analysis(method = 'pls',estimate_group = 'Group')
```


OPLS

```R
MAE |>
  EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
               estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |>
  EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
  EMP_filter(feature_condition = pvalue < 0.05) |>
  EMP_dimension_analysis(method = 'opls',estimate_group = 'Sex') |>
  EMP_scatterplot(estimate_group='Sex',show='p12html',ellipse=0.6) ###### Visualization
```

### EMP_enrich_analysis

Make the enrichment after EMP_diff_analysis

```R
MAE |>
  EMP_assay_extract(experiment = 'geno_ec') |>
  EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
  EMP_enrich_analysis(keyType ='ec',KEGG_Type = 'KEGG',
                      pvalue<0.05,pvalueCutoff=1,species = 'all') 
```

Make the enrichment after EMP_diff_analysis and Visualization

```R
MAE |>
  EMP_assay_extract(experiment = 'geno_ec') |>
  EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
  EMP_enrich_analysis(keyType ='ec',KEGG_Type = 'KEGG',
                      pvalue<0.05,pvalueCutoff=1,species = 'all') |>
  EMP_dotplot()

MAE |>
  EMP_assay_extract(experiment = 'geno_ec') |>
  EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
  EMP_enrich_analysis(keyType ='ec',KEGG_Type = 'KEGG',
                      pvalue<0.05,pvalueCutoff=1,species = 'all') |>
  EMP_netplot()
```

Transcriptomic data

```R
MAE |>
  EMP_assay_extract(experiment = 'host_gene') |>
  EMP_feature_convert(from = 'symbol',to='entrezid',species='Human') |>
  EMP_diff_analysis(method = 'DESeq2',
                    .formula = ~Group,p.adjust = 'fdr') |> 
  EMP_enrich_analysis(keyType ='entrezid',
                      KEGG_Type ='KEGG',pvalue<0.05,
                      pvalueCutoff=0.05,species = 'hsa') |>
  EMP_dotplot()
```


### EMP_marker_analysis

To estimate the improtance of feature by Boruta algorithm

```R
MAE |>
  EMP_marker_analysis(experiment = 'geno_ec',method = 'boruta',
                      estimate_group = 'Group') |>
  EMP_filter(feature_condition = Boruta_decision!= 'Rejected') ###### select the Confirmed and Tentative feature
```

regression or classify by randomforest

```R
MAE |>
  EMP_marker_analysis(experiment = 'geno_ec',method = 'randomforest',
                      estimate_group = 'Group') 
MAE |>
  EMP_marker_analysis(experiment = 'geno_ec',method = 'randomforest',
                      estimate_group = 'Education_Years') 
```

regression or classify by xgboost

1. regression

```R
MAE |>
  EMP_marker_analysis(experiment = 'geno_ec',method = 'xgboost',
                      xgboost_run = "regression",
                      estimate_group = 'Education_Years',objective = 'reg:linear')

```

2. For two categories classify

```R
MAE |>
  EMP_marker_analysis(experiment = 'geno_ec',method = 'xgboost',
                      xgboost_run = "classify",
                      estimate_group = 'Group',objective = 'binary:logistic')
```

3. For multible categories classify

```R
MAE |>
  EMP_marker_analysis(experiment = 'geno_ec',method = 'xgboost',
                      xgboost_run = "classify",
                      estimate_group = 'Status',objective = 'multi:softmax',
                      num_class=3) ###### num_class is necessary
```

Lasso regression

```R
MAE |>
  EMP_marker_analysis(experiment = 'geno_ko',
                      method = 'lasso',estimate_group = 'Education_Years') |>
  EMP_filter(feature_condition = lasso_coe >0) ### Select the imprortant feature
```

### EMP_GSEA_analysis

Based on cor analysis

```R
MAE |>
  EMP_GSEA_analysis(experiment = 'geno_ko',method='cor',
                    estimate_group = 'BMI',cor_method = 'spearman',
                    threshold_r = 0.3,threshold_p = 0.05, ###### filter by coe and pvalue
                    pvalueCutoff = 0.05,keyType = 'ko')
```

Based on diff analysis

```R
MAE |>
  EMP_diff_analysis(experiment = 'geno_ko',method='DESeq2',.formula = ~0+Group,
                    group_level=c('Group_A','Group_B')) |>
  EMP_GSEA_analysis(method='log2FC',pvalue<0.05,
                    keyType = 'ko',KEGG_Type = 'KEGG')
```

Based on signal2Noise

```R
MAE |>
  EMP_GSEA_analysis(experiment = 'geno_ko',method='signal2Noise',
                    estimate_group = 'Group',
                    pvalueCutoff = 0.05,keyType = 'ko')
```

Visualization

```R
MAE |>
  EMP_GSEA_analysis(experiment = 'geno_ko',method='signal2Noise',
                    estimate_group = 'Group',
                    pvalueCutoff = 0.05,keyType = 'ko') |>
  EMP_curveplot(geneSetID='map00680')
MAE |>
  EMP_GSEA_analysis(experiment = 'geno_ko',method='signal2Noise',
                    estimate_group = 'Group',
                    pvalueCutoff = 0.05,keyType = 'ko') |>
  EMP_dotplot(color='p.adjust',showCategory=10) 
MAE |>
  EMP_GSEA_analysis(experiment = 'geno_ko',method='signal2Noise',
                    estimate_group = 'Group',
                    pvalueCutoff = 0.05,keyType = 'ko') |>
  EMP_netplot(showCategory=5) 
```

### EMP_WGCNA_cluster_analysis

```R
MAE |>
  EMP_assay_extract('geno_ec') |>
  EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85)
MAE |>
  EMP_assay_extract('geno_ko') |>
  EMP_WGCNA_cluster_analysis(RsquaredCut = 0.8,mergeCutHeight=0.4)
```

### EMP_WGCNA_cor_analysis

From one experiment

```R
WGCNA_COR_result <-MAE |>
  EMP_assay_extract('geno_ec')  |> 
  EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
  EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)  |>
  EMP_WGCNA_cor_analysis(coldata_to_assay = c('BMI','PHQ9','GAD7','HAMD','SAS','SDS'),
                         method='spearman') 
```

Visualization

```R
MAE |>
  EMP_assay_extract('geno_ec')  |> 
  EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
  EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)  |>
  EMP_WGCNA_cor_analysis(coldata_to_assay = c('BMI','PHQ9','GAD7','HAMD','SAS','SDS'),
                         method='spearman') |>
  EMP_heatmap_plot(palette = 'Spectral')
```

Filter the interesting module and make the enrichment analysis

```R
MAE |>
  EMP_assay_extract('geno_ec')  |> 
  EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
  EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)  |>
  EMP_WGCNA_cor_analysis(coldata_to_assay = c('BMI','PHQ9','GAD7','HAMD','SAS','SDS'),
                         method='spearman') |>
  EMP_heatmap_plot(palette = 'Spectral') |>
  EMP_filter(feature_condition = WGCNA_color == 'brown' ) |> 
  EMP_diff_analysis(method = 'DESeq2',.formula = ~Group) |>
  EMP_enrich_analysis(keyType = 'ec',KEGG_Type = 'MKEGG') |>
  EMP_dotplot()
```

From two different experiments

```R
k1 <- MAE |>
  EMP_assay_extract('geno_ec')  |> 
  EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
  EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)

k2 <- MAE |>
  EMP_assay_extract('host_gene',pattern = c('A1BG','A1CF','A2MP1','AACS'),pattern_ref = 'feature')

(k1 + k2) |>
  EMP_WGCNA_cor_analysis(method='spearman') |>
  EMP_heatmap_plot(palette = 'Spectral') 
```

### EMP_multi_analysis

prepare the result

```R
k1 <- MAE |> 
  EMP_assay_extract('geno_ec') |> 
  EMP_diff_analysis(method = 'DESeq2',.formula = ~Group) |> 
  EMP_enrich_analysis(pvalue < 0.05,keyType = 'ec')
```

Use different differential analysis methods to mimic the statistical results of the same feature across different cohorts.

```R
k2 <- MAE |> 
  EMP_assay_extract('geno_ec') |> 
  EMP_diff_analysis(method = 'wilcox.test',estimate_group = 'Group') |>
  EMP_enrich_analysis(pvalue < 0.05,keyType = 'ec')

k3 <- MAE |> 
  EMP_assay_extract('geno_ko') |> 
  EMP_diff_analysis(method = 'DESeq2',.formula = ~Group) |> 
  EMP_enrich_analysis(pvalue < 0.05,keyType = 'ko')

k4 <- MAE |>
  EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
               estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |>
  EMP_diff_analysis(method = 'DESeq2',.formula = ~Group) |> 
  EMP_enrich_analysis(pvalue < 0.05,keyType = 'cpd')
```

Melt diff analysis of the same feature from two result

```R
(k1+k2) |> 
  EMP_multi_analysis(method = 'feature',combineMethod='edgington',p.adjust = 'BH') 
```

Melt diff analysis of the same feature from two result and make a combine enrichment

```R
(k1+k2) |> 
  EMP_multi_analysis(method = 'same_feature_enrich',keyType = 'ec',combineFun='ActivePathways') 
```

 Melt diff analysis of the different feature from multi-results and make a combine enrichment

```R
(k1+k4) |> 
  EMP_multi_analysis(method = 'diff_feature_enrich')
```

Visualization

```R
(k1+k2) |> 
  EMP_multi_analysis(method = 'same_feature_enrich',keyType = 'ec',combineFun='ActivePathways') |>
  EMP_dotplot()
```

### EMP_boxplot

For assay

```R
MAE |> 
  EMP_assay_extract('host_gene',pattern = 'A1BG',pattern_ref = 'feature') |>
  EMP_boxplot(method='t.test',estimate_group='Group')
```

For alpha analysis

```R
MAE |> 
  EMP_assay_extract('taxonomy') |> 
  EMP_alpha_analysis() |>
  EMP_boxplot(method='t.test',estimate_group='Group') 
```

### EMP_scatterplot

```R
MAE |> 
 EMP_assay_extract('taxonomy') |> 
 EMP_collapse(estimate_group = 'Species',collapse_by = 'row') |>
 EMP_dimension_analysis(method = 'pcoa',distance = 'bray',estimate_group = 'Group') |>
 EMP_scatterplot(show='p12html') ### eg. p12,p12html,p23,p23htm
```

### EMP_dotplot

```R
MAE |> 
 EMP_assay_extract('taxonomy') |> 
 EMP_collapse(estimate_group = 'Species',collapse_by = 'row') |>
 EMP_dimension_analysis(method = 'pcoa',distance = 'bray',estimate_group = 'Group') |>
 EMP_scatterplot(show='p12html') ### eg. p12,p12html,p23,p23htm
```

### EMP_netplot

```R
MAE |>
  EMP_GSEA_analysis(experiment = 'geno_ko',method='signal2Noise',
                    estimate_group = 'Group',
                    pvalueCutoff = 0.05,keyType = 'ko') |>
  EMP_netplot(showCategory=10) 
```

### EMP_curveplot

```R
MAE|>
  EMP_GSEA_analysis(experiment = 'geno_ko',method='signal2Noise',
                    estimate_group = 'Group',
                    pvalueCutoff = 0.05,keyType = 'ko') |>
  EMP_curveplot(geneSetID='map00680')
```

### EMP_heatmap_plot

For cor analysis

```R
k1 <- MAE |>
  EMP_assay_extract('taxonomy') |>
  EMP_collapse(estimate_group = 'Genus',collapse_by = 'row') |>
  EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
  EMP_filter(feature_condition = pvalue<0.05)

k2 <- MAE |>
  EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
               estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |>
  EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
  EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) > 1.5)

(k1 + k2) |> EMP_cor_analysis(method = 'spearman') |>
  EMP_heatmap_plot() ###### Visualization
```

For WGCNA

````R
MAE |>
  EMP_assay_extract('geno_ec')  |> 
  EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
  EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,mergeCutHeight=0.4)  |>
  EMP_WGCNA_cor_analysis(coldata_to_assay = c('BMI','PHQ9','GAD7','HAMD','SAS','SDS'),method='spearman') |>
  EMP_heatmap_plot(palette = 'Spectral') 
````

### EMP_volcanol_plot

```R
MAE |>
  EMP_assay_extract(experiment = 'geno_ec') |>
  EMP_diff_analysis(method='DESeq2',.formula = ~Group)  |>
  EMP_volcanol_plot(show='html')
```

### EMP_filter

For  MultiAssayExperiment

```R
MAE |>
  EMP_filter(sample_condition = BMI>20 & Sex == 'M') |>
  EMP_summary()
```

For  EMPT

```R
MAE |>
  EMP_assay_extract(experiment = 'host_gene') |>
  EMP_identify_assay(method = 'edgeR',estimate_group = 'Group') |>
  EMP_diff_analysis(method='DESeq2',.formula = ~Group) |>
  EMP_filter(feature_condition = fdr < 0.05)

MAE |>
  EMP_assay_extract(experiment = 'taxonomy') |>
  EMP_alpha_analysis() |>
  EMP_filter(sample_condition = shannon >3 | invsimpson >19)
```

Precise selection

```R
MAE |>
  EMP_assay_extract(experiment = 'taxonomy') |>
  EMP_alpha_analysis() |>
  EMP_filter(sample_condition = shannon >3 | invsimpson >19,
             filterSample = c('P11774','P31579'),
             action = 'kick') # Accurately kick samples after sample_condition
```

### EMP_history

###### from EMPT

```R
MAE |>
  EMP_assay_extract(experiment = 'taxonomy') |>
  EMP_alpha_analysis() |>
  EMP_history()
```

###### from EMP

```R
k1 <- MAE |>
  EMP_assay_extract('taxonomy') |>
  EMP_collapse(estimate_group = 'Genus',collapse_by = 'row') |>
  EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
  EMP_filter(feature_condition = pvalue<0.05)

k2 <- MAE |>
  EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
               estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') |>
  EMP_diff_analysis(method='DESeq2', .formula = ~Group) |>
  EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) > 1.5)

(k1 + k2) |> EMP_cor_analysis(method = 'spearman') |>
  EMP_heatmap_plot()  |>
  EMP_history()  ###### get all the history about what happened to the object
```


### EMP_result

Obtain the result from EMPT

```R
MAE |> 
  EMP_assay_extract('geno_ec') |>
  EMP_alpha_analysis()|> 
  EMP_diff_analysis(method = 'DESeq2',.formula = ~Group) |>
  EMP_enrich_analysis(pvalue<0.05,keyType='ec',pvalueCutoff=0.05) -> result

diff_re <- result |> EMP_result(info = 'EMP_diff_analysis')
alpha_re <- result |> EMP_result(info = 'EMP_alpha_analysis')
enrich_re <- result |> EMP_result(info = 'EMP_enrich_analysis')
```



Inject external result into EMPT

1. get a EMPT object

```R
MAE |> 
  EMP_assay_extract('taxonomy') |>
  EMP_collapse(estimate_group = 'Genus',method = 'sum',
               collapse_by = 'row',action = 'add') -> obj  
```

2. get the raw data from the EMPT

```R 
MAE |> 
  EMP_assay_extract('taxonomy') |>
  EMP_collapse(estimate_group = 'Genus',method = 'sum',
               collapse_by = 'row',action = 'get') -> assay_data  
```

3. caculate the result from other packages

```R
assay_data <- assay_data |> tibble::column_to_rownames('primary')
shannon_index <- vegan::diversity(assay_data,index = 'shannon') 
new_result <- tibble::tibble(primary=names(shannon_index),new_shannon=shannon_index)
```

4. inject the new result into EMPT object

```R
EMP_result(obj,
           value_name = 'new_alpha',
           affect_when_sample_changed=0,
           affect_when_feature_changed=1,
           attribute='primary',
           attribute2='normal',source='user_import') <- new_result

obj |> 
  EMP_filter(sample_condition  = new_shannon >2)
```

