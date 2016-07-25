all: description.pdf supp.pdf 

description.pdf: description.tex ref.bib Figure1.pdf Figure2.pdf Figure3.pdf Figure4.pdf Figure5.png Figure6.pdf
	pdflatex $<
	bibtex description.aux
	pdflatex $<
	pdflatex $<

supp.pdf: supp.tex description.pdf
	pdflatex $< 
	pdflatex $< 

Figure1.pdf: pics/pipeline.svg
	inkscape -A $@ -D $<

Figure2.pdf: pics/figure_edgeR.tex ../simulations/plot_alpha.pdf ../simulations/plot_power.pdf 
	cd pics; pdflatex `basename $<`
	pdfcrop pics/figure_edgeR.pdf $@

Figure3.pdf: pics/spatial_fdr.svg
	inkscape -A $@ -D $<

Figure4.pdf: ../simulations/plot_FDR.pdf
	cp $< $@

Figure5.png: ../real/analysis/pics/labelled_results.svg
	inkscape -e $@ -d 300 -D $<

Figure6.pdf: pics/figure_real_extra.tex ../real/analysis/pics/Cytobank_43324_4FI_induction.png ../real/analysis/pics/Cytobank_43324_4FI_withdrawal.png ../real/analysis/pics/extra_bar.png
	cd pics; pdflatex `basename $<`
	pdfcrop pics/figure_real_extra.pdf $@

