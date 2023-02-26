# UW-Biostatistics-Capstone-Project

*Group Project*

Sponsored by NanoString Technologies Inc. 

Project: Benchmarking and Evaluation Cell Typing Methods for Spatial Transcriptomic Technologies

## Benchmarking supervised cell typing 

We identify the cell types in a kidney biopsy from a lupus nephritis patient using reference profiles from the Human Cell Atlas, which provide average expression profiles for each kidney cell type based on their single cell RNA-seq studies. Since there is no definitive gold standard cell type data for this sample, we cannot determine with certainty whether a cell type call is accurate. To evaluate the quality of our cell typing results, we rely on known kidney biology and use two metrics. 

The first metric evaluates the accuracy of glomerular cell assignments to glomeruli, as kidneys contain substructures called glomeruli that comprise distinct cell types. Our expectation is that 100% of podocytes and glomerular endothelial cells will be assigned to glomeruli, while 0% of proximal tubule cells will be assigned to glomeruli. 

Second, we record the intensity of a few "marker proteins" that are highly informative of cell type for each cell. We compare the marker protein expression differences between cell types that are either positive or negative for the markers, such as CD45 for immune cells versus other cell types, or PanCK for cytokeratin+ versus cytokeratin- cell types. These differences are scored using logarithmic transformation.  

To evaluate the performance of our supervised cell typing approach, we compare InsituType to other methods such as SingleR, scVI, CHETAH, SeuratV3, SingleCellNet, scPred, and SVM. Additionally, we conduct all analyses in a consistent computing environment and record computation times.

> Methods

- InSitutype: Developed by Nanostring, Insitutype utilizes a likelihood model to enable supervised cell typing from reference datasets via a Bayes classifier and unsupervised or semi-supervised cell typing via an Expectation Maximization algorithm. Insitutype is designed to address the challenges of sparse and multi-modal spatial transcriptomics data, and it employs an escalating subsampling scheme to handle large datasets efficiently. Compared to existing cell typing methods, Insitutype offers several advantages, including its ability to handle large datasets, its use of a likelihood model to weigh the evidence from every transcript in a cell, its incorporation of alternative data types such as images and spatial context, and its ability to identify new clusters alongside reference cell types. 

- SingleR: An unbiased cell typing method for scRNA-seq by leveraging reference transcriptomic datasets of pure cell types to infer the cell of origin of each single cell independently. 

- CHETAH (Characterization of Cell Types Aided by Hierarchical classification): A scRNA-seq classifier by hierarchical clustering of the reference data. The classification tree enables a step-wise, top-to-bottom classification.  

- SeuratV3: A single-cell transcriptomics classifier can anchor diverse datasets together, enabling us to integrate single-cell measurements not only across scRNA-seq, but also across 	different modalities (e.g. scATAC-seq).  

- SingleCellNet: A random forest classifier to learn cell type-specific gene pairs from cross-platform and cross-species datasets and thus quantitatively assesses cell identity at a single-cell resolution. 


