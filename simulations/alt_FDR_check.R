# This checks that the spatial FDR is being properly controlled, in the opposite scenario of FDR_check.R.
# Namely, we set up a situation where the non-DA subpopulations are very tight and the DA subpopulation is relatively diffuse (but still small).
# This could theoretically trigger loss of control when measuring the FDR with large partitions, as each subpopulation gets assigned into one partition, regardless of volume.

require(cydar)
require(edgeR)
source("functions.R")
ofile <- "results_FDR_alt.txt"
existing <- FALSE

############################################ 
# Setting up

nmarkers <- 30
for (it in seq_len(50)) {
    all.hypers <- rbind(matrix(rnorm(1000*nmarkers, mean=runif(1), sd=0.1), ncol=nmarkers), # Scrambling location to avoid edge effects.
                        matrix(runif(1)*0.5+1, nrow=1000, ncol=nmarkers),
                        matrix(runif(1)*0.5+2, nrow=1000, ncol=nmarkers),
                        matrix(runif(1)*0.5+3, nrow=1000, ncol=nmarkers),
                        matrix(runif(1)*0.5+4, nrow=1000, ncol=nmarkers))
    colnames(all.hypers) <- paste0("X.", seq_len(nmarkers))
    all.pvalues <- c(rbeta(1000, 1, 100),
                     rep(runif(1), 1000), # Identical p-values because they're all on the same spot.
                     rep(runif(1), 1000),
                     rep(runif(1), 1000),
                     rep(runif(1), 1000))
    is.DA <- rep(c(TRUE, FALSE), c(1000, 4000))

    # Controlling the FDR, spatially or naively.
    all.results <- list()
    for (con in c("naive", "spatial")) {
        if (con=="naive") {
            qval <- p.adjust(all.pvalues, method="BH")
        } else {
            qval <- spatialFDR(all.hypers, all.pvalues)
        }
        is.sig <- qval <= 0.05
        all.results[[paste0(con, "_detected")]] <- sum(is.sig & is.DA)

        # Assessing the FDR using partitions of varying size.
        for (width in c(0.2, 0.4, 0.6, 0.8, 1)) { 
            partitions <- apply(all.hypers[is.sig,,drop=FALSE], 1, function(x) { paste(floor(x/width), collapse=".") })
            all.tests <- table(partitions)
            false.pos <- table(partitions[!is.DA[is.sig]])
            m <- match(names(false.pos), names(all.tests))
            fdr <- sum(false.pos/all.tests[m])/length(all.tests)
            all.results[[paste0(con, "_", width)]] <- fdr
        }
    }

    write.table(file=ofile, data.frame(rbind(unlist(all.results))),
            quote=FALSE, sep="\t", col.names=!existing, append=existing, row.names=FALSE)
    existing <- TRUE
    gc()
}

############################################ 
# Checking the results

blah <- read.table("results_FDR_alt.txt", header=TRUE)
colMeans(blah, na.rm=TRUE)

# In practice, loss of FDR control across the partitions doesn't seem to occur.
# This seems to be because the p-values for the tight non-DA subpopulation are highly correlated (if not identical).
# This means that the effective extent of multiple testing (in terms of number of tests) is lower than what would be expected for independent tests.
# Furthermore, despite our best efforts, even a small amount of diffuseness gets picked up in separate partitions.
# These two factors reduce the frequency of rejection of the non-DA subpopulations, and mitigate the effect on the observed FDR when it does occur.
# We also see an improvement in power with the spatial FDR with respect to picking up the DA hyperspheres.
# (The opposite would occur in FDR_check.R, but the loss of control with the naive method is so severe in that case, it's not worth mentioning.)

############################################ 
# End.
