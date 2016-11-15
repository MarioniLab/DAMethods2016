library(Biobase)
library(flowCore)

# A function for resampling cells for simulation.
resampleCells <- function(cd, setting=1L) {
    current.exprs <- list()
    samples <- sampleNames(cd)
    ci <- cellIntensities(cd)
    pooled <- t(ci[,order(cellData(cd)$sample.id, cellData(cd)$cell.id),drop=FALSE])
    
    # Sampling values to recapitulate the original library sizes.
    cells.per.sample <- tabulate(cellData(cd)$sample.id)
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

# Adds a bunch of cells at a specified point location to each sample, proportional to the number of cells.
addPointDifference <- function(current.exprs, chosen, loc, prop.DA) {
    for (i in chosen) {
        to.sample <- nrow(current.exprs[[i]])
        extras <- matrix(loc, round(to.sample*prop.DA), ncol(current.exprs[[i]]))
        current.exprs[[i]] <- rbind(current.exprs[[i]], extras)
    }
    return(current.exprs)
}

# Function to write to FCS file.
dumpToFile <- function(output.dir, current.exprs) { 
    out.files <- list()
    for (f in seq_along(current.exprs)) {
        curexp <- current.exprs[[f]]
        p <- AnnotatedDataFrame(data.frame(name=colnames(curexp), 
                                           desc=colnames(curexp),
                                           range=apply(curexp, 2, function(x) { diff(range(x)) }),
                                           minRange=apply(curexp, 2, min),
                                           maxRange=apply(curexp, 2, max)))
        ff <- flowFrame(curexp, p)
        fname <- file.path(output.dir, paste0(f, ".fcs"))
        suppressWarnings(write.FCS(file=fname, ff))
        out.files[[f]] <- fname
    }
    return(unlist(out.files))
}

