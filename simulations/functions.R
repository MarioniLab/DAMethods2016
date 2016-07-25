# A function for simulation.

resampleCells <- function(cd, setting=1L) {
    current.exprs <- list()
    samples <- attributes(cd)$samples
    pooled <- t(cd[,order(attributes(cd)$sample.id, attributes(cd)$cell.id),drop=FALSE])
    
    # Sampling values to recapitulate the original library sizes.
    cells.per.sample <- tabulate(attributes(cd)$sample.id + 1)
    for (i in seq_along(samples)) {
        to.sample <- cells.per.sample[i]
        
        if (setting==1L) {
            # Totally random (Poisson noise).
            chosen <- sample(nrow(pooled), to.sample, replace=TRUE)
        } else {
            # Adding some Gamma noise. If it's too variable, however, it'll break,
            # as you get "stacks" of rare cells that make the p-value very small.
            # (These have large dispersions, but the dispersion trend itself is 
            # too high to detect outliers at this point; saddlepoint failure?)
            weights <- rgamma(nrow(pooled), 1e-2, 1e-2)
            chosen <- sample(nrow(pooled), to.sample, replace=TRUE, prob=weights)
        }
        
        # Mocking up things to save to file.
        current.exprs[[samples[i]]] <- pooled[chosen,,drop=FALSE]
    }
    return(current.exprs)
}

