# Methods for differential abundance analyses of mass cytometry data

To repeat the analyses:

- Download the FCS files from Cytobank for accession number 43324, and put them in `refdata/Cytobank_43324`.
- Enter `refdata` and run `preprocess.R` to transform and gate on the intensities, followed by `gen_objects.R` to obtain the hypersphere counts.
- Enter `simulations` and run `bsubber.h` to perform simulations for type I error control (in `edgeR_check.R`) and spatial FDR control (in `FDR_check.R`).
- Enter `real/analysis` and run (in order) `detect_da.R`, to detect DA hyperspheres; `plot_da_initial.R`, to construct the _t_-SNE plots; `plot_da_extra.R`, to construct subpopulation- and time-point-specific plots; and `plot_dispersions.R`, to construct the dispersion plot.
- Enter `real/neighbours` and run `get_neighbors.R` to construct plots of distances to nearest neighbours.
- Enter `real/robustness` and run `robustness.R` to obtain DA statistics with different parameter settings.

