# Methods for differential abundance analyses of mass cytometry data

To repeat the real data analyses for the MEF reprogramming time course:

- Download the FCS files from Cytobank for accession number 43324, and put them in `refdata/Cytobank_43324`.
- Enter `refdata` and run `preprocess_43324.R` to transform and gate on the intensities, followed by `gen_objects_43324.R` to obtain the hypersphere counts.
- Enter `real/analysis` and run (in order) `detect_da.R`, to detect DA hyperspheres; `plot_da_initial.R`, to construct the _t_-SNE plots; `plot_da_extra.R`, to construct subpopulation- and time-point-specific plots; and `plot_dispersions.R`, to construct the dispersion plot. 
Note that the _t_-SNE plots are stochastic and may not be the same as in the paper.

Alternatively, you can look at `workflow.Rmd` in `real/workflow`, which provides an annotated workflow of the processing and analysis of the Oct4-GFP time course data.

To run the additional analyses for the time course data:

- Enter `real/neighbours` and run `get_neighbors.R` to construct plots of distances to nearest neighbours.
- Enter `real/robustness` and run `robustness.R` to obtain DA statistics with different parameter settings.
- Enter `real/clustering` and run `citrus_test.R` to obtain DA clusters from CITRUS.
Then run `map_onto_tsne_citrus.R` to map the CITRUS cluster centres onto the _t_-SNE plot from the hypersphere analysis.

To run the analyses for the BMMC data:

- Download the FCS files from Cytobank for accession number 44185, and put them in `refdata/Cytobank_44185`.
- Enter `refdata` and run `preprocess_44185.R` to transform and gate on the intensities, followed by `gen_objects_44185.R` to obtain the hypersphere counts.
- Enter `real/secondary` and run (in order)  `detect_da.R`, to detect DA hyperspheres; and `plot_da.R`, to construct the _t_-SNE plots.

To repeat the simulations, enter `simulations` and run:

- `edgeR_check.R`, to check type I error control.
- `FDR_check.R`, to check spatial FDR control.
- `cluster_sim.R`, to check performance of clustering.
- `radius_sim.R`, to check the effect of increasing the radius on hypersphere positions.
- `shift_check.R`, to check the type I error control on intensity-shifted data.

