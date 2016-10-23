bsub -R "rusage[mem=16000]" -n 2 -e "loge.err" R CMD BATCH --no-save edgeR_check.R
bsub -R "rusage[mem=16000]" -n 2 -e "logf.err" R CMD BATCH --no-save FDR_check.R 
bsub -R "rusage[mem=16000]" -n 4 -e "logc.err" R CMD BATCH --no-save citrus_comp.R 
