all: short.pdf supp.pdf 

short.pdf: short.tex ref.bib Figure1.pdf Figure2.pdf 
	pdflatex $<
	bibtex short.aux
	pdflatex $<
	pdflatex $<

supp.pdf: supp.tex \
	../real/neighbours/Cytobank_43324_4FI.pdf ../real/neighbours/Cytobank_43324_NG.pdf ../real/neighbours/Cytobank_43324_NN.pdf ../real/neighbours/nvd_Cytobank_43324_4FI.pdf \
	../real/analysis/pics/dispersions.pdf ../simulations/plot_alpha_power.pdf \
	pics/Figure_FDRpartitions.pdf ../simulations/plot_FDR.pdf \
	../real/analysis/pics/Cytobank_43324_4FI_markers.png ../real/analysis/pics/Cytobank_43324_NG_logFC.png ../real/analysis/pics/Cytobank_43324_NG_markers.png ../real/analysis/pics/Cytobank_43324_NN_logFC.png ../real/analysis/pics/Cytobank_43324_NN_markers.png \
   	mef_timing.pdf ../real/analysis/pics/novel_timecourse.pdf \
	../simulations/cluster_setup.pdf ../simulations/plot_cluster.pdf ../real/clustering/Cytobank_43324_4FI_citrus.png \
    pics/Figure_RadiusSim.pdf ../simulations/radius_setup.pdf \
	../refdata/pics/ecdf_wnorm2.pdf ../refdata/pics/ecdf_raw2.pdf ../refdata/pics/ecdf_rnorm2.pdf \
	pics/Figure_BMMC.pdf ../real/secondary/pics/IL10_markers.png \
	../simulations/plot_signal.pdf
	pdflatex $< 
	bibtex supp.aux
	pdflatex $< 
	pdflatex $< 

pics/Figure_RadiusSim.pdf: pics/figure_radius.tex ../simulations/plot_radius.pdf
	cd pics; pdflatex `basename $<`
	pdfcrop pics/figure_radius.pdf $@

pics/Figure_FDRpartitions.pdf: pics/spatial_fdr_partitions.svg
	cd pics; inkscape -A `basename $@` -D `basename $<`

pics/Figure_BMMC.pdf: ../real/secondary/pics/labelled_results.svg
	inkscape -A Figure_BMMC_temp.pdf -D $<
	gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default \
		    -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages \
			-dCompressFonts=true -r150 -sOutputFile=$@ Figure_BMMC_temp.pdf
	rm Figure_BMMC_temp.pdf

Figure1.pdf: pics/unified_figure.svg
	inkscape -A $@ -D $<

Figure2.pdf: ../real/analysis/pics/labelled_results.svg
	inkscape -A Figure5_temp.pdf -D $<
	gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default \
		    -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages \
			-dCompressFonts=true -r300 -sOutputFile=$@ Figure5_temp.pdf
	rm Figure5_temp.pdf

mef_timing.pdf: ../real/analysis/pics/extra_results.svg
	inkscape -A Figure6_temp.pdf -D $<
	gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default \
		    -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages \
			-dCompressFonts=true -r150 -sOutputFile=$@ Figure6_temp.pdf
	rm Figure6_temp.pdf


