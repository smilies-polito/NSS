# Neuronal Spike Shapes (NSS)

This repository contains the code of the Neuronal Spike Shapes (NSS), a simple approach for analyzing the electrophysiological profiles of cells based on their Action Potential (AP) waveforms.

The NSS method aims to explore the heterogeneity of cell types and states by summarizing the AP waveform into a triangular representation complemented by a set of derived electrophysiological (EP) features.

## Release notes

NSS v1.0: 

* First release of NSS.

## How to cite

### NSS Primary publication

* Martini, L., Amprimo, G., Di Carlo, S., Olmo, G., Ferraris, C., Savino, A. and Bardini, R., Neuronal Spike Shapes (NSS): a simple approach to study electrophysiological data for heterogeneity investigation, 2023 (submitted to Elsevier Computers in Biology and Medicine).

NSS validation relies on two datasets of murine cortical interneurons.

### PatchSeqDataset
It leverages the patch-seq technique, that combines patch-clamp with single-cell RNA-Seq (scRNA-Seq) to collect EP and transcriptomic profiles of 4,200 mouse visual cortical GABAergic interneurons, reconstructing the morphological conformation of 517 of them. Primary publication:

* Gouwens, N. W., Sorensen, S. A., Baftizadeh, F., Budzillo, A., Lee, B. R., Jarsky, T., ... & Zeng, H. (2020). Integrated morphoelectric and transcriptomic classification of cortical GABAergic cells. Cell, 183(4), 935-953, https://doi.org/10.1016/j.cell.2020.09.057.


### PatchClampDataset
It provides EP data from 1,938 neurons of adult mouse visual cortical neurons and morphological reconstructions for 461 of them. Primary publication:

* Gouwens, N.W., Sorensen, S.A., Berg, J. et al. Classification of electrophysiological and morphological neuron types in the mouse visual cortex. Nat Neurosci 22, 1182–1195 (2019), https://doi.org/10.1038/s41593-019-0417-0.



## Experimental setup

Follow these steps to setup for reproducing the experiments provided in _Martini et al., 2023_.
1) Install `Singularity` from https://docs.sylabs.io/guides/3.0/user-guide/installation.html:
	* Install `Singularity` release 3.10.2, with `Go` version 1.18.4
	* Suggestion: follow instructions provided in _Download and install singularity from a release_ section after installing `Go`
	* Install dependencies from: https://docs.sylabs.io/guides/main/admin-guide/installation.html
2) Clone the NSS repository in your home folder
```
git clone https://github.com/smilies-polito/NSS.git
```
3) Move to the NSS source subfolder, and build the singularity container with 
```
mv NSS/source
sudo singularity build NSS.sif NSS.def
```
or
```
mv NSS/source
singularity build --fakeroot NSS.sif NSS.def
```

# Reproducing NSS analysis

In _Martini et al., 2023_, NSS validation relies on multimodal analysis over both datasets:

* **Hard Validation**: A EP and transcriptomic joint and multimodal analysis over _PatchSeqDataset_ 
* **Soft Validation**: A EP and cell-type-based joint analysis over _PatchClampDataset_ 

## Data required

In order to reproduce the analysis, it is necessary to gather and rename required data files, organizing them in the repository folders as provided below.

### Hard Validation: A EP and transcriptomic joint and multimodal analysis over _PatchSeqDataset_ 

Three files are required, and it is necessary to rename them and organize them in the repository folders as provided below.

1) The file reporting the EP features of all cells in the dataset and supporting EP analysis: 
```
NSS/data/PatchSeqDataset/PatchSeq_EP_features.csv
```
This file (provided in the repository) is generated with a custom script calling the DANDI (https://dandiarchive.org) APIs to access the DANDISET 000020 (https://doi.org/10.48324/dandi.000020/0.210913.1639). This script leverages the `dandi-cli` tool (10.5281/zenodo.3692138) to access raw EP data for each cell and compute EP features.

Credits:  Allen Institute for Brain Science (2020). Patch-seq recordings from mouse visual cortex. Available from https://dandiarchive.org/dandiset/000020/ and https://github.com/dandisets/000020.

Primary publication: Gouwens, N. W., Sorensen, S. A., Baftizadeh, F., Budzillo, A., Lee, B. R., Jarsky, T., ... & Zeng, H. (2020). Integrated morphoelectric and transcriptomic classification of cortical GABAergic cells. Cell, 183(4), 935-953, https://doi.org/10.1016/j.cell.2020.09.057.

2) The file reporting the metadata of the dataset, connecting EP and transcriptomic IDs of cells: 
```
NSS/data/PatchSeqDataset/PatchSeq_metadata.csv
```
This file (provided in the repository) can be downloaded at https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/0f/86/0f861fcb-36d5-4d3a-80e6-c3c04f34a8c7/20200711_patchseq_metadata_mouse.csv

Credits: Allen Institute for Brain Science (2020). Patch-seq recordings from mouse visual cortex. Available from https://portal.brain-map.org/explore/classes/multimodal-characterization/multimodal-characterization-mouse-visual-cortex.

Primary publication: Gouwens, N. W., Sorensen, S. A., Baftizadeh, F., Budzillo, A., Lee, B. R., Jarsky, T., ... & Zeng, H. (2020). Integrated morphoelectric and transcriptomic classification of cortical GABAergic cells. Cell, 183(4), 935-953, https://doi.org/10.1016/j.cell.2020.09.057.

3) The file reporting the scRNA-seq counts, measuring gene expression and supporting transcriptomic analysis: 
```
NSS/data/PatchSeqDataset/count.csv
```
This file (not provided in the repository) can be downloaded and extracted from https://data.nemoarchive.org/other/AIBS/AIBS_patchseq/transcriptome/scell/SMARTseq/processed/analysis/20200611/20200513_Mouse_PatchSeq_Release_count.v2.csv.tar.

Credits: Allen Institute for Brain Science (2020). Patch-seq recordings from mouse visual cortex. Available from https://data.nemoarchive.org/other/AIBS/AIBS_patchseq/transcriptome/scell/SMARTseq/processed/analysis/20200611/.

Primary publication: Gouwens, N. W., Sorensen, S. A., Baftizadeh, F., Budzillo, A., Lee, B. R., Jarsky, T., ... & Zeng, H. (2020). Integrated morphoelectric and transcriptomic classification of cortical GABAergic cells. Cell, 183(4), 935-953, https://doi.org/10.1016/j.cell.2020.09.057.

### Soft Validation: A EP and cell-type-based joint analysis over _PatchClampDataset_ 
Two files are required, and it is necessary to rename them and organize them in the repository folders as provided below.
1) The file reporting the EP features of all cells in the dataset and supporting EP analysis:
```
NSS/data/PatchClampDataset/PatchClamp_EP_features.csv
```
This file (provided in the repository) is generated with a custom script calling the Allen SDK APIs and using the IPFX library to access precumputed EP features for each cell in the Cell Types Database.

Credits:  Allen Institute for Brain Science (2023). Allen SDK. Available from https://github.com/alleninstitute/allensdk. Allen Institute for Brain Science (2023). IPFX. Available from https://github.com/alleninstitute/ipfx.

Primary publication: Gouwens, N.W., Sorensen, S.A., Berg, J. et al. Classification of electrophysiological and morphological neuron types in the mouse visual cortex. Nat Neurosci 22, 1182–1195 (2019), https://doi.org/10.1038/s41593-019-0417-0.

2) The file reporting the metadata of the dataset, including CRE cell lines-based cell type labels of cells: 
```
NSS/data/PatchClampDataset/PatchClamp_metadata.csv
```
This file (provided in the repository) can be downloaded at http://celltypes.brain-map.org/cell_types_specimen_details.csv. 

Credits:  Allen Institute for Brain Science (2023). Cell Types dataset. Available from http://celltypes.brain-map.org/data (DOWNLOAD CELL FEATURE DATA button).

Primary publication: Gouwens, N.W., Sorensen, S.A., Berg, J. et al. Classification of electrophysiological and morphological neuron types in the mouse visual cortex. Nat Neurosci 22, 1182–1195 (2019), https://doi.org/10.1038/s41593-019-0417-0.

## Reproducing the analysis interactively within the NSS Singularity container
To run analyses manually launch the NSS Singularity container, move to `/source`, and launch the scripts as follows.

First of all, launch the NSS Singularity container
```
cd source
singularity shell --no-home --bind  /local/path/to/NSS:/local/path/to/home/ NSS.sif
```
This will run a shell within the container, and the following prompt should appear:
```
Singularity>
```
Using this prompt, follow the steps below. 

### Hard Validation: EP and transcriptomic joint and multimodal analysis over _PatchSeqDataset_ 
1) Run the EP analysis indicating:
 
 * the type of dataset employed (`PatchSeqDataset`)
 * the desired clustering cardinality (`k=2` or `k=3`)
 
For example:
```
Singularity> python3 EP_analysis.py PatchSeqDataset 3
```
This script leverages EP files in `data/PatchSeqDataset/` to generate NSS clustering files and images of NSS-based embeddings labeled with NSS clustering labels in `output/PatchSeqDataset/`.

2) Run the transcriptomic analysis: 
```
Singularity> Rscript transcriptomic_analysis.r
```
This script leverages the count file in `data/PatchSeqDataset/` to generate transcriptomic clustering files and images of transcriptomic-based embeddings labeled with transcriptomics-based cell types labels in `output/PatchSeqDataset/`. It employs NSS labels from the NSS clustering files to generate transcriptomic-based embeddings labeled with NSS-based labels for both the entire dataset and the Sst subset. Then, performing Differential Expression (DE) analysis over NSS-based cell subsets it generates lists of DE genes for each NSS cluster.

3) Run the Gene Ontology (GO) Enrichment Analysis (EA)
```
Singularity> Rscript GO_enrichment_analysis.r
```

To reproduce the provided analysis in the paper, upload the list of DE genes for the EC0 NSS-based cluster (Pvalb-like EP profiles) from the whole dataset to the `ShinyGO` tool (http://bioinformatics.sdstate.edu/go/) and generate plots for the three ontologies employed, selecting them from the Pathway Database menu: `GO Biological Process`, `GO Cellular Component` and `GO Molecular Function`. Keep all the other parameters at default values. 

Note: since this step is not fully automatic, the repository provides pre-computed GO EA files to generate figures. 


### Soft Validation: EP and cell-type-based joint analysis over _PatchClampDataset_ 

1) Run the EP and cell types analysis indicating:
 
 * the type of dataset employed (`PatchClampDataset`)
 * the desired clustering cardinality (`k=2` or `k=3`)

For example:
```
Singularity> python3 EP_analysis.py PatchClampDataset 2
```
This script leverages files in `data/PatchClampDataset/` to generate NSS clustering files in and images of NSS-based embeddings labeled with NSS clustering labels. It generates cell type label-based groups in `output/PatchSeqDataset/Cre_lines_groups/`, and images of NSS-based embeddings labeled with Cre cell lines-based cell types labels in `output/PatchSeqDataset/`.

## Reproducing the analysis running the NSS Singularity container

To reproduce the analyses from _Martini et al., 2023_, run the `NSS.sif` container with the required commandline arguments: a keyword to indicate the type of analysis (`hardValidation` or `softValidation`), and the EP clustering cardinality (`k=2` or `k=3`).

### Hard Validation: EP and transcriptomic joint and multimodal analysis over _PatchSeqDataset_ 
```
singularity run --no-home --bind  /local/path/to/NSS:/local/path/to/home/ NSS.sif hardValidation
```

### Soft Validation: EP and cell-type-based joint analysis over _PatchClampDataset_ 
```
singularity run --no-home --bind  /local/path/to/NSS:/local/path/to/home/ NSS.sif softValidation
```

## Repository structure

```
|
├── data                                       // Data files
|    ├── PatchSeqDataset                       // Data files for the PatchSeqDataset
|    |    ├── PatchSeq_EP_features.csv         // EP features file for EP analysis of the PatchSeqDataset (Credits:  Allen Institute for Brain Science (2020). Patch-seq recordings from mouse visual cortex. Available from https://dandiarchive.org/dandiset/000020/ and https://github.com/dandisets/000020.)
|    |    ├── PatchSeq_metadata.csv            // metadata file for EP analysis of the PatchSeqDataset (Credits: Allen Institute for Brain Science (2020). Patch-seq recordings from mouse visual cortex. Available from https://portal.brain-map.org/explore/classes/multimodal-characterization/multimodal-characterization-mouse-visual-cortex.
|    |    └── count.csv                        // scRNA-seq data for transcriptomic analysis of the PatchSeqDataset (Credits: Allen Institute for Brain Science (2020). Patch-seq recordings from mouse visual cortex. Available from https://data.nemoarchive.org/other/AIBS/AIBS_patchseq/transcriptome/scell/SMARTseq/processed/analysis/20200611/.)
|    | 
|    └── PatchClampDataset                     // Data files for the PatchClampDataset
|         ├── PatchClamp_EP_features.csv       // EP features file for the PatchClampDataset (Credits:  Allen Institute for Brain Science (2023). Allen SDK. Available from https://github.com/alleninstitute/allensdk. Allen Institute for Brain Science (2023). IPFX. Available from https://github.com/alleninstitute/ipfx.)
|         └── PatchClamp_metadata.csv          // metadata file for the PatchClampDataset, including cell lines based cell type labels for cell types analysis of PatchClampDataset (Credits:  Allen Institute for Brain Science (2023). Cell Types dataset. Available from http://celltypes.brain-map.org/data - DOWNLOAD CELL FEATURE DATA button.)
|    
|   
├── source                                    // Scripts for EP, transcriptomic and cellTypes analyses
|    ├── EP_analysis.py                       // Python script for EP analysis
|    ├── cell_lines_types_analysis.py         // Python script for cell lines-based cell type labels analysis
|    ├── transcriptomic_analysis.R            // R script for transcriptomic analysis
|    ├── GO_enrichment_analysis.R             // R script for cell types analysis
|    └── NSS.def                              // NSS Singularity recipe

│
├── output                                     // Output of the NSS analysis
|    ├── PatchSeqDataset                       // Output files and images for the PatchSeqDataset
|    |    ├── NSS_clusters                     // Files linking cell IDs and NSS clustering labels
|    |    |    ├── NSS_clusters_k3.csv         // File linking cell IDs and NSS clustering labels for k=3 clustering
|    |    |    └── NSS_clusters_k2.csv         // File linking cell IDs and NSS clustering labels for k=2 clustering
|    |    |
|    |    ├──  GO_enrichment_analysis                   // Files supporting the gene ontology enrichment analysis
|    |    |    ├── NSS_cluster0_k3_DE_genes.csv         // list of DE genes in NSS cluster EC0 (k=3) from the whole dataset for gene-ontology enrichment analysis
|    |    |    ├── enrichment_cellular.csv              // data for gene-ontology enrichment analysis images (downloaded from http://bioinformatics.sdstate.edu/go/, GO Cellular Component)
|    |    |    ├── enrichment_molecular.csv             // data for gene-ontology enrichment analysis images (downloaded from http://bioinformatics.sdstate.edu/go/, GO Molecular Function)
|    |    |    └── enrichment_biological.csv            // data for gene-ontology enrichment analysis images (downloaded from http://bioinformatics.sdstate.edu/go/, GO Biological Process)
|    |    |
|    |    ├── MISC                                      // miscellaneous files supporting intermediate analysis steps
|    |    |    └── ...                                       
|    |    |
|    |    └── IMAGES                                                                     // Images generated by the analysis 
|    |         ├── histograms_feature_distributions_slope_deep.png                       // Slope_deep feature distribution histograms (Figure 4, top left)
|    |         ├── histograms_feature_distributions_dv_ratio.png                         // dv_ratio feature distribution histograms (Figure 4, top right)
|    |         ├── histograms_feature_distributions_up_down_ratio.png                    // UpDown_ratio feature distribution histograms (Figure 4, bottom left)
|    |         ├── histograms_feature_distributions_dv_deep.png                          // dv_deep feature distribution histograms (Figure 4, bottom right)
|    |         ├── silhouette.png                                                        // Silhouette plot for the identification of optimal number of clusters (Figure 5)
|    |         ├── NSS_embedding_NSS_labels_k3.html                                      // Interactive file showing NSS labels on NSS embedding for k=2 (Figure 6, left)
|    |         ├── NSS_embedding_NSS_labels_k2.html                                      // Interactive file showing NSS labels on NSS embedding for k=3 (Figure 6, right)
|    |         ├── correlation_matrix_k2.png                                             // Spearman’s correlation matrix of NSS features to clustering label for k=2 (Figure 7, top)
|    |         ├── correlation_matrix_k3.png                                             // Spearman’s correlation matrix of NSS features to clustering label for k=3 (Figure 7, bottom)
|    |         ├── KDE_plot_k2_slope_deep.png                                            // KDE plots of the distribution of features slope_deep for k=2 (Figure 8, bottom left)
|    |         ├── KDE_plot_k2_dv_ratio.png                                              // KDE plots of the distribution of features slope_deep and dv_ratio for k=2 (Figure 8, bottom right)
|    |         ├── KDE_plot_k3_slope_deep.png                                            // KDE plots of the distribution of features slope_deep for k=3 (Figure 8, top left)
|    |         ├── KDE_plot_k3_dv_ratio.png                                              // KDE plots of the distribution of features slope_deep and dv_ratio for k=3 (Figure 8, top right)
|    |         ├── bar_plots_mean_accuracy_decrease.png                                  //  Mean accuracy decrease in the cluster label prediction done by RF model when shuffling values of NSS feature (Figure 9)
|    |         ├── transcriptomic_embedding_marker_expression.pdf                        // Cell type markers expression levels on transcriptomic embedding (Figure 11)
|    |         ├── transcriptomic_embedding_transcriptomic_labels.pdf                    // Transcriptional labels provided by PatchSeqDataset metadata on transcriptomic embedding (Figure 12)
|    |         ├── transcriptomic_embedding_NSS_labels_k3.pdf                            // NSS labels on transcriptomic embedding for k=2 (Figure 13, left)
|    |         ├── transcriptomic_embedding_NSS_labels_k2.pdf                            // NSS labels on transcriptomic embedding for k=2 (Figure 13, right)
|    |         ├── transcriptomic_embedding_transcriptomic_labels_Sst_subset.pdf         // Transcriptional labels provided by PatchSeqDataset metadata on transcriptomic embedding of the Sst subset (Figure 14, left)
|    |         ├── transcriptomic_embedding_NSS_labels_k3_Sst_subset.pdf                 // NSS labels for k=3 on transcriptomic embedding of the Sst subset (Figure 14, right)
|    |         ├── GO_biological_plot.pdf                                                // Gene ontology enrichment analysis plot (downloaded from http://bioinformatics.sdstate.edu/go/, GO Biological Process) (Figure 15, top)
|    |         ├── GO_molecular_plot.pdf                                                 // Gene ontology enrichment analysis plot (downloaded from http://bioinformatics.sdstate.edu/go/, GO Molecular Function) (Figure 15, middle)
|    |         ├── GO_cellular_plot.pdf                                                  // Gene ontology enrichment analysis plot (downloaded from http://bioinformatics.sdstate.edu/go/, GO Cellular Component) (Figure 15, bottom)
|    |         ├── transcriptomic_Kcnc2_violin_plot.pdf                                  // Violin plots showing the Kcnc2 gene expression levels distributions across cells in each NSS cluster with k=3 (Figure 16, third row)
|    |         └── transcriptomic_Kcnn2_violin_plot.pdf                                  // Violin plots showing the Kcnn2 gene expression levels distributions across cells in each NSS cluster with k=3 (Figure 16, fourth row)
|    |
|    |    
|    └── PatchClampDataset                     // Output files and images for the PatchClampDataset
|         ├── NSS_clusters                     // Files linking cell IDs and NSS clustering labels
|         |    ├── NSS_clusters_k3.csv         // File linking cell IDs and NSS clustering labels for k=3 clustering
|         |    └── NSS_clusters_k2.csv         // File linking cell IDs and NSS clustering labels for k=2 clustering
|         |
|         ├── cre_lines_clusters            // Files linking cell IDs and Cre lines-based cell types labels
|         |    ├── lines_k2_cl0.csv         // File linking cell IDs and Cre lines-based cell types labels for NSS cluster EC0 with k=2
|         |    ├── lines_k2_cl1.csv         // File linking cell IDs and Cre lines-based cell types labels for NSS cluster EC1 with k=2
|         |    ├── lines_k3_cl0.csv         // File linking cell IDs and Cre lines-based cell types labels for NSS cluster EC0 with k=3
|         |    ├── lines_k3_cl1.csv         // File linking cell IDs and Cre lines-based cell types labels for NSS cluster EC1 with k=3
|         |    └── lines_k3_cl2.csv         // File linking cell IDs and Cre lines-based cell types labels for NSS cluster EC2 with k=3
|         |
|         ├── cre_lines_all              // Files linking cell IDs and Cre lines-based cell types labels for all cells
|         |    └── cre_lines.csv         // File linking cell IDs and Cre lines-based cell types labels for all cells
|         |
|         ├── MISC                               // miscellaneous files supporting intermediate analysis steps
|         |    ├── feature_analysis.csv          // File supporting feature analysis
|         |    ├── cell_lines_types.csv          // File supporting Cre lines-based cell types analysis
|         |    ├── lines_k2_cl0_type.csv         // File supporting Cre lines-based cell types analysis for NSS cluster EC0 with k=2
|         |    ├── lines_k2_cl1_type.csv         // File supporting Cre lines-based cell types analysis for NSS cluster EC1 with k=2
|         |    ├── lines_k3_cl0_type.csv         // File supporting Cre lines-based cell types analysis for NSS cluster EC0 with k=3
|         |    ├── lines_k3_cl1_type.csv         // File supporting Cre lines-based cell types analysis for NSS cluster EC1 with k=3
|         |    └── lines_k3_cl2_type.csv         // File supporting Cre lines-based cell types analysis for NSS cluster EC2 with k=3
|         |
|         └── IMAGES                                              // Images generated by the analysis 
|              ├── NSS_embedding_Cre_lines_labels_k2.html         // Interactive file showing Transcriptional labels provided by PatchSeqDataset metadata on NSS embedding (Figure 17, top center)
|              ├── NSS_embedding_Cre_lines_labels_k3.html         // Interactive file showing Transcriptional labels provided by PatchSeqDataset metadata on NSS embedding (Figure 17, top center)
|              ├── NSS_embedding_NSS_labels_k3.html               // Interactive file showing NSS labels on NSS embedding for k=3 (Figure 17, bottom left)
|              └── NSS_embedding_NSS_labels_k2.html               // Interactive file showing NSS labels on NSS embedding for k=2 (Figure 17, bottom right)
|
|    
└── README.md                                                     // This README file          
```
  
