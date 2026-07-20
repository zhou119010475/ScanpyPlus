# ScanpyPlus
A set of files to do single-cell analysis in python

To use any of these script collections, just run these two lines in your *python kernel / Jupyter notebook*:
```python
sys.path.append('/home/ubuntu/tools/ScanpyPlus')
import Scanpyplus
```
## Citation:
[He, Lim and Sun et al.](https://www.cell.com/cell/fulltext/S0092-8674(22)01415-5)
A human fetal lung cell atlas uncovers proximal-distal gradients of differentiation and key regulators of epithelial fates

<details>
  <summary><b>DeepTree feature selection</b></summary>

![image](https://user-images.githubusercontent.com/4110443/146441826-a4079e4c-c9de-4d93-9ebe-3e1c07227eb1.png)


Among the functions in *Scanpyplus*, there's also a function to do feature gene selection (*DeepTree* algorithm). It removes garbage among highly variable genes, mitigate batch effect if you remove garbage batch by batch, and increases signal-to-noise ratio of the top PCs to promote rare cell type discovery.

[Here](https://nbviewer.jupyter.org/github/Peng-He-Lab/ScanpyPlus/tree/master/DeepTree_algorithm_demo.ipynb) is a [notebook](https://github.com/Peng-He-Lab/ScanpyPlus/tree/master/DeepTree_algorithm_demo.ipynb) to use DeepTree algorithm to "de-noise" highly-variable genes and improve initial clustering. 

A *MATLAB* implementation can be found [here](https://github.com/brianpenghe/Matlab-genomics).

This algorithm can be potentially used to reduce batch effect when fearing overcorrection, especially comparing conditions or time points. Two notebooks are provided showing "soft integration" of [embryonic limb](Soft_integration_limb.ipynb) and [pancreas](Soft_integration_pancreas.ipynb) data.

More conventionally, `HVGbyBatch` calls `sc.pp.highly_variable_genes` per batch, stores per-batch boolean flags and an aggregate count `highly_variable_n`
```python
Scanpyplus.HVGbyBatch(adata, batch_key='sample', min_disp=0.5)
```

`HVG_cutoff(adata, range_int=10, cutoff=5000, HVG_var='highly_variable_n', fig_size=(8,6))` Plots HVG count vs. threshold and prints the minimal threshold meeting your target count. This can also be used to determine the n for DeepTree-derived genes across multiple samples.
```python
i = Scanpyplus.HVG_cutoff(adata, cutoff=5000)
```

`HVG_Venn_Upset` compares HVG sets across provided flags (e.g., per-batch).
```python
Scanpyplus.HVG_Venn_Upset(adata, ['highly_variable_A','highly_variable_B'])
```

`DeepTree` clusters HVG expression (via `snsCluster`), cuts the gene dendrogram, and flags genes that fall into sufficiently large clades (`Deep`). Returns `[bdata, test, test1, test2]` (filtered object plus Seaborn cluster maps).
```python
bdata, test, test1, test2 = Scanpyplus.DeepTree(
    adata, MouseC1ColorDict2={False:'#000000', True:'#00FFFF'},
    cell_type='leiden', gene_type='highly_variable', Cutoff=0.8, CladeSize=2)
```

`DeepTree_per_batch` runs DeepTree within each batch on that batch's HVGs and stores `Deep_<batch>` flags plus aggregate `Deep_n`.
```python
Scanpyplus.DeepTree_per_batch(adata, batch_key='batch', Cutoff=0.8, CladeSize=2)
```

`DeepTree2` is a lighter scipy-only variant that returns an AnnData with a `Deep` gene flag (no Seaborn plots).
```python
bdata = Scanpyplus.DeepTree2(adata, Cutoff=0.8, CladeSize=2)
```
</details>

<details>
  <summary><b>Doublet Cluster Labeling (DouCLing)</b></summary>

![unnamed](https://user-images.githubusercontent.com/4110443/146441371-e7b4bec2-9e87-4a9d-98ad-3f3401ce13ed.jpg)

There are 4 types of doublets:

![image](https://user-images.githubusercontent.com/4110443/146040113-1c1b27e6-453e-48fa-a4e8-786ff8c759ec.png)

Cross-sample doublets can usually be identified by hastags or genetic backgrounds. Theoretically, (n-1)/n of doublets can be identified as cross-sample doublets when n samples with different hashtags/genetics are pooled equally.

Heterotypic doublets can sometimes trick data scientists into thinking they are a new type of dual-feature cell type like NKT cells etc. 
Heterotypic doublets are usually identified by matching individual cells to synthetic doublets regardless of manually curated clusters. Algorithms like Scrublet can remove a substantial part of doublet cells but not all of them. The survivor doublets can still aggregate into tiny clusters picked up by the annotaters when doing subclustering. Doublets of rarer cell types are also often missed, which obscures the discoveries of new cell types and states.
To leverage the input from biologists' manual parsing and the increased sensitivity of cluster-average signatures, I introduce here an alternative approach to facilitate heterotypic doublet cluster identification. This approach scans through individual tiny clusters and look for its "Parent 2" that gives it a unique feature that's different from its sibling subclusters sharing the same "Parent 1". 
[A notebook using published PBMC data](https://nbviewer.jupyter.org/github/Peng-He-Lab/ScanpyPlus/tree/master/DOUblet_Cluster_Labeling.ipynb) is provided.

`DouCLing(adata, hi_type, lo_type, rm_genes=[], print_marker_genes=False, fraction_threshold=0.6)` scores high-resolution clusters against low-resolution compartments and returns a table with putative parents, hypergeometric p-values, and an `Is_doublet_cluster` flag.
```python
scores = Scanpyplus.DouCLing(adata, hi_type='leiden_R', lo_type='leiden',
                             rm_genes=['TYMS','MKI67'], fraction_threshold=0.6)
```
</details>

<details>
  <summary><b>Other functions in Scanpyplus:</b></summary>

### An alternative way to call doublet subclusters based on *Scrublet* and [the gastrulation paper](https://www.nature.com/articles/s41586-019-0933-9)
`Bertie(adata,Resln=1,batch_key='batch')` was written with the help from [K. Polanski](https://github.com/ktpolanski). This script aggregates *Scrublet* scores from subclusters and makes threshold cuts based on subcluster p-values. And this is done batch by batch.
```python
Scanpyplus.Bertie(adata, Resln=1, batch_key='sample')
```
A variant version `Bertie_preclustered` applies the same subcluster-enrichment logic using an existing cluster assignment (e.g., Leiden).
```python
Scanpyplus.Bertie_preclustered(adata, cluster_key='leiden', batch_key='donor')
```

### Color utilities:
You can extract the color dict of a variable from an anndata object using `ExtractColor(adata,obsKey='louvain',keytype=int)`, 

and manually edit the color dict, and then use it to update colors in `adata.uns` using `UpdateUnsColor`.

```python
Scanpyplus.UpdateUnsColor(adata, ColorDict, obsKey='cell_type')
```

You can also force a category to render white using `MakeWhite`
```python
Scanpyplus.MakeWhite(adata, obsKey='condition', whiteCat='background')
```

### Metadata / obs (observation) helpers:
You can plot sankey graph between two variables of an anndata object using `ScanpySankey`. 

```python
Scanpyplus.ScanpySankey(adata, 'louvain', 'cell_type')
```

`orderGroups(adata, groupby='leiden')` returns dendrogram-ordered category labels for a grouping key.
```python
order = Scanpyplus.orderGroups(adata, groupby='louvain')
```

`remove_barcode_suffix` removes the suffix after the '-' in the cell (barcode) name.
```python
Scanpyplus.remove_barcode_suffix(adata)
```

`CopyMeta` copies obs/var metadata columns from one object to another (optionally overwriting).
```python
Scanpyplus.CopyMeta(ref, query, overwrite=False)
```

`AddMeta` stores a dataframe of obs values per each cell into an object.
```python
Scanpyplus.AddMeta(adata, df)
```

`AddMetaBatch` reads a dataframe of obs values per batch into an object. This format of metadata (rows are batch names, columns are obs categories) is more common, compact and human readable that is usually stored in *Excel* spreadsheets.
```python
Scanpyplus.AddMetaBatch(adata, batch_table, batch_key='donor')
```

`ExtractMetaBatch` is the reverse of `AddMetaBatch`: returns the modal value of each obs column per batch.
```python
Scanpyplus.ExtractMetaBatch(adata, batch_key='batch')
```

`ConvertString` casts selected obs columns to string.
```python
Scanpyplus.ConvertString(adata, ['sample', 'donor'])
```

`MapCategories` renames categorical levels in an obs column using a mapping dict.
```python
Scanpyplus.MapCategories(adata, 'cell_type', {'old':'new'})
```

`dropmeta` drops specified columns from `adata.obs` if present.
```python
Scanpyplus.dropmeta(adata, ['tmp_score', 'unused'])
```

`SubclusterAll` Leiden-subclusters each parent cluster and writes hierarchical labels (e.g. `0_1`) to `result_key`.
```python
Scanpyplus.SubclusterAll(adata, parent_key='leiden', resolution=0.2, result_key='leiden_R')
```

### Gene / var (variable)  metadata:
`OrthoTranslate` translates mouse genes to human orthologs and filter out poorly conserved genes, based on ortholog table that can be derived from Biomart etc.

### Converting file types:
`file2gz` creates .gz files which is useful for creating artificial 10X files.
```python3
Scanpyplus.file2gz('matrix.mtx')
```

`Scanpy2MM` saves an *anndata* into *MatrixMarket* form (optionally writing embeddings into metadata via `write2Dobsm`).
```python3
Scanpyplus.Scanpy2MM(adata, prefix='mm_out/', write2Dobsm=['X_umap'])
```

`mtx2df` reads *MatrixMarket* files into a dataframe.
```python
Scanpyplus.mtx2df('matrix.mtx', 'features.tsv', 'barcodes.tsv')
```

### Matrix utilities:
`GetRaw` copy .raw (or layer) to .X 
```python
Scanpyplus.GetRaw(adata)
```

`CalculateRaw` → rebuild a raw count matrix from log-transformed counts based on `n_counts` 
```python
Scanpyplus.CalculateRaw(adata, scaling_factor=10000)
```

`CalculateRawAuto` reverses log1p without `n_counts` by treating each cell's minimum nonzero value as 1.
```python
Scanpyplus.CalculateRawAuto(adata)
```

`CheckGAPDH` Quickly inspects first few entries for a given gene, optionally from a sparse matrix.
```python
Scanpyplus.CheckGAPDH(adata, sparse=True, gene='GAPDH')
```

`FindSimilarGenes` Computes gene–gene correlations and returns a sorted similarity series.
```python
Scanpyplus.FindSimilarGenes(adata, genename='XIST')
```

For large matrices, cells can be `DownSample`d based on labels such as cell types.
```python
Scanpyplus.DownSample(adata, 'cell_type', downsampleTo=500)
```

Sometimes `PseudoBulk` profiles are also useful to generate, whether it's the mean, median or max.
```python
Scanpyplus.PseudoBulk(adata, group_key='sample', layer='counts')
```

### Embedding utilities:
`ShiftEmbedding` creates a platter that juxtaposes subsets of the data (batches, stages etc.) to visualize side by side.
```python
Scanpyplus.ShiftEmbedding(adata, 'batch', embedding='X_umap')
```

`CopyEmbedding` copies the embedding of one object to another.
```python
Scanpyplus.CopyEmbedding(ref, query, embedding='X_umap')
```

### Plotting stacked barplots of cell-type/condition proportions:
`celltype_per_stage_plot` and `stage_per_celltype_plot` plot horizontal and vertical bar plots respectively based on two metadata variables (cell type and stage, for example).

### UMAP plotting utilities:
`Plot3DimUMAP` generates a 3D plot (by *plotly*) of the UMAP after sc.tl.umap produces the 3D coordinates.

`plot_umap_with_labels` Plots UMAP with categorical labels drawn on top; auto-adjusts text to avoid overlap and saves PNG/PDF.
```python
Scanpyplus.plot_umap_with_labels(adata, 'cell_type', 'plot.png', 'plot.pdf')
```

`PlotCrosstab` plots a heatmap of a crosstab / confusion matrix between two obs keys.
```python
Scanpyplus.PlotCrosstab(adata, row_key='Predicted', col_key='cell_type',
                        outfile='confusion_matrix.png')
```

### Gene-level calculation and plotting:
`returnDEres` extracts a tidy DE dataframe (scores, logFC, p-values, optionally pts) from `adata.uns['rank_genes_groups']`, optionally removing mito/ribo genes.
```python
Scanpyplus.returnDEres(adata, column='0', key='rank_genes_groups')
```

`DEmarkers` calculates, filters and plots differentially expressed genes between two populations.

`GlobalMarkers` calculates marker genes for every cell cluster and filters them.


`GPT_annotation_genes` Ranks genes per cluster (Wilcoxon) and returns a formatted text prompt listing top markers per cluster for quick human/GPT annotation.  
```python
Scanpyplus.GPT_annotation_genes(adata, leiden_key='leiden', top_n=10)
```

`ClusterGenes` transposes a log-transformed *adata* object and performs clustering and dimension reduction to classify genes.

`Dotplot2D` plots the expression levels of a gene across two metadata categories (e.g. samples and cell types). It can be used to trace maternal contimation by plotting XIST and check a key gene's expression patterns against cell types and age etc.

### Seaborn utilities:
`snsSplitViolin` plots splitviolin plots for two populations.
```python
Scanpyplus.snsSplitViolin(adata, ['GATA3'], celltype='leiden', celltypelist=['0','1'])
```

`snsCluster` plots clustermaps using an *anndata object* as input. This has been helped by Bao Zhang from [Zhang lab](https://github.com/ZhangHongbo-Lab)

`markSeaborn` marks specific genes on a *Seaborn* plot.

`extractSeabornRows` extracts the rowlabels of a *Seaborn object* and saves into a *Series*.

### Plotting Venn / UpSet diagram:
`UpSetFromLists` (also available in *pandasPlus*) plots upset plots (bar plots of each category of intersections) from lists of lists. `HVG_Venn_Upset` is the AnnData-oriented variant for comparing HVG flags.

### Treemap:
`Treemap(adata, output='temp', branchlist=['project','batch'], width=1000, height=700, title='title')`
Treemap of obs counts by hierarchical categories; saves PDF and CSV.
```python
Scanpyplus.Treemap(adata, output='out/treemap', branchlist=['donor','cell_type'])
```

### Label transfer / logistic regression:
`LogisticRegressionCellType` can learn the defining features of a variable (such as cell type) of the reference object and predict the corresponding labels of a query object. 

The saved model files can also be re-used to predict a new query object in future by `LogisticPrediction`.

Helpers for saved models:
- `LoadLogitModel(model_addr)` / `LoadLogitGenes(genecsv)` load a pickled model and its gene list.
- `ExtractLogitScores(adata, model, CT_genes)` returns per-class probabilities and `lr_score`.
- `Model2Coeff(model_pkl, genelistcsv)` returns a coefficient dataframe (classes × genes).
```python
Scanpyplus.LogisticPrediction(adata, 'model.pkl', 'genes.csv', scores=True, compute_umap=False)
```

### QC and data loading:
`QC` annotates mt/ribo/hb genes (human or mouse) and runs `sc.pp.calculate_qc_metrics`.
```python
Scanpyplus.QC(adata, species='human')
```

`read_folder` loads 10X / AnnData files from sample subdirectories (optional concat). `read_numbered_folders` is a convenience wrapper for numbered directories.
```python
adatas = Scanpyplus.read_folder('/path/to/samples', data_pattern='filtered_feature_bc_matrix.h5',
                                batch_key='sample', concat=False)
adata = Scanpyplus.read_numbered_folders('/path/to/samples', min_dir=1, max_dir=12, concat=True)
```

`iRODS_stats_starsolo(samples)` downloads STARsolo Gene/cr3 matrices via `iget`, computes cell count and median n_counts per library, and returns a QC table.

`OverlapBarcode` reports shared barcode prefixes across a list of AnnData objects.
```python
Scanpyplus.OverlapBarcode([adata1, adata2], ['s1', 's2'])
```

`extract_barcode_prefix` returns the barcode string before `-`.

`import_souporcell` imports Souporcell cluster / status assignments into `adata.obs`.
```python
Scanpyplus.import_souporcell(adata, 'clusters.tsv', prefix='')
```

`run_celltypist` wraps CellTypist annotation (and optional UMAP plots).
```python
Scanpyplus.run_celltypist(adata, model_name='Immune_All_Low.pkl', majority_voting=True)
```

### Spatial / iStar inputs:
Helpers to export Visium-style AnnData into iStar input files:
- `iStar_extract_locs_raw` → `locs-raw.tsv`
- `iStar_extract_counts_raw` → `cnts.tsv`
- `iStar_extract_scalefactors` → `pixel-size-raw.txt` and `radius-raw.txt`
- `iStar_generate_all_inputs(h5ad_path, formatted_name, output_dir='.')` runs all of the above.
```python
Scanpyplus.iStar_generate_all_inputs('sample.h5ad', formatted_name='sample', output_dir='istar_in')
```

### Also mirrored from pandasPlus:
`show_graph_with_labels`, `DF2Ann`, `UpSetFromLists`, `zscore`, `Ginni`, `cellphonedb_p2adjMat`, `cellphonedb_n_interaction_Mat`, `cellphonedb_mat_per_interaction`, and `ListsOverlap` are available on `Scanpyplus` as well (see **Functions in pandasPlus** below).
</details>

<details>
  <summary><b>Functions in pandasPlus:</b></summary>

`DF2Ann` converts a dataframe into an *anndata* object.
```python
pandasPlus.DF2Ann(df)
```

`UpSetFromLists` plots an upset plot (barplot of Venn diagram intersections) based on lists of lists.
```python
pandasPlus.UpSetFromLists([list_a, list_b], labels=['A','B'])
```

`show_graph_with_labels` plots an interaction graph using edges to represent connection strength (edges shown when value ≥ 0.9).
```python
pandasPlus.show_graph_with_labels(adjacency_df)
```

Dataframe values can also be used to calculate `zscore` and `Ginni` coefficients.
```python
pandasPlus.zscore(df, dropna=True, axis=1)
pandasPlus.Ginni(df)  # appends a 'Ginni' row
```

`cellphonedb_p2adjMat`, `cellphonedb_n_interaction_Mat`, and `cellphonedb_mat_per_interaction` reformat CellPhoneDB `pvalues.txt` outputs into cell-type × cell-type matrices (significant pair strings, counts, or a single interaction's p-values).
```python
pandasPlus.cellphonedb_p2adjMat('pvalues.txt', pval=0.05)
pandasPlus.cellphonedb_n_interaction_Mat('pvalues.txt', pval=0.05)
pandasPlus.cellphonedb_mat_per_interaction('EGFR_TGFB1', 'pvalues.txt')
```

`ListsOverlap(A, B)` returns a dataframe of shared-element counts between each list in `A` and each list in `B` (lists of lists).
```python
pandasPlus.ListsOverlap([['a','b'], ['b','c']], [['b'], ['a','c']])
```
</details>

<details>
  <summary><b>Other scripts in this repository:</b></summary>

`Pythonplus.CheckSource(Function2Check)` prints the source code of a Python function (useful for inspecting Callables interactively).

`cellranger_summary.py` merges Cell Ranger `outs/summary.csv` files across sample subdirectories into one TSV:
```bash
python cellranger_summary.py --base_dir /path/to/cellranger_outs --output summary.tsv
```
</details>
