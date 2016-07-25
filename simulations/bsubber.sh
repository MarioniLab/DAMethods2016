bsub -R "rusage[mem=16000]" -n 2 -e "loge.err" Rdevel CMD BATCH --no-save edgeR_check.R
bsub -R "rusage[mem=16000]" -n 2 -e "logf.err" Rdevel CMD BATCH --no-save FDR_check.R 
