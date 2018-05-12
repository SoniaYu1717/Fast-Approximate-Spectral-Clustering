# Fast-Approximate-Spectral-Clustering
The paper ‘Fast Approximate Spectral Clustering’ proposed to use K-means and RP trees as a data preprocessing step before the actual spectral clustering step, to improve the overall runtime efficicency while upholding the clustering quality.

In this project, to extend the research of Donghui Yan et al., we applied and evaluated **Partitioning Around Medoids (PAM)** and **Entropy Weighted K-means (EWKM)** as the data preprocessing step.

To run the program, type the following command in R

	>source("mean.R")

Note:
	kernllab, MASS, gtools, cluster and wskm are requested libraries. To install these packages, use the R standard method install.packages()
