bsub -R "rusage[mem=16000]" -n 1 -e "loge.err" -J DA_edgeR R CMD BATCH --no-save edgeR_check.R
bsub -R "rusage[mem=16000]" -n 1 -e "logf.err" -J DA_FDR R CMD BATCH --no-save FDR_check.R 
bsub -R "rusage[mem=16000]" -n 4 -e "logci.err" -J DA_citrus R CMD BATCH --no-save citrus_comp.R 
bsub -R "rusage[mem=16000]" -n 4 -e "logcl.err" -J DA_cluster R CMD BATCH --no-save cluster_sim.R 
bsub -R "rusage[mem=16000]" -n 1 -e "logs.err" -J DA_shift R CMD BATCH --no-save shift_check.R
bsub -R "rusage[mem=16000]" -n 1 -e "logr.err" -J DA_radius R CMD BATCH --no-save radius_sim.R

