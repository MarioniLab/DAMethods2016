# This performs logicle-transformation and gating on intensities for each data set. 

library(ncdfFlow)
host.dir <- "."
for (dataset in 1:3) {

    if (dataset<=3L) {
        ref.dir <- "Cytobank_43324"
        if (dataset==1L) {
            extra <- "4FI"
        } else if (dataset==2L) {
            extra <- "NG"
        } else if (dataset==3L) {
            extra <- "NN"
        }
        ref.out.dir <- paste0(ref.dir, "_", extra)
        pattern <- sprintf("^%s_[0-9][0-9].fcs", extra)
        
        # Reading in the files.
        fcs <- list.files(file.path(host.dir, ref.dir), full=TRUE, pattern=pattern)
        nfs <- read.ncdfFlowSet(fcs[1])
        toignore <- c(match(c("Cell_length", "Time", "beadDist", "barcode"), colnames(nfs)), grep("^BC-[0-9]", colnames(nfs)))

        # Gate specifications.
        gate.up <- match(c("DNA-1(Ir191)Dd", "DNA-2(Ir193)Dd"), colnames(nfs))
        gate.down <- NULL
    }

    ########################################################
    # Estimating the parameters of the logicle transformation, after pooling all samples.
    new.exprs <- list()
    new.nfs <- suppressWarnings(suppressMessages(read.ncdfFlowSet(fcs)))
    for (i in seq_along(fcs)) {
        new.exprs[[i]] <- exprs(new.nfs[[i]])
    }
    new.exprs <- do.call(rbind, new.exprs)
    replaced.ff <- nfs[[1]]
    exprs(replaced.ff) <- new.exprs
    parameters(replaced.ff)$minRange <- apply(new.exprs, 2, min)
    parameters(replaced.ff)$maxRange <- apply(new.exprs, 2, max)
    parameters(replaced.ff)$range <- parameters(replaced.ff)$maxRange - pmax(parameters(replaced.ff)$minRange, 0)
    lgcl <- estimateLogicle(replaced.ff, channels=colnames(replaced.ff)[-toignore], type='data')
    trans.ff <- transform(replaced.ff, lgcl)
    trans.nfs <- transform(new.nfs, lgcl)
   
    ########################################################
    # Need to gate on specified aspects (gate.up, gate.down).

    if (!is.null(gate.up) || !is.null(gate.down)) {
        all.gates <- c(gate.up, gate.down)
        threshold.above <- rep(c(TRUE, FALSE), c(length(gate.up), length(gate.down)))
        med.per.gate <- apply(exprs(trans.ff)[,all.gates], 2, median) # Using the median and MAD across _all_ events.
        mad.per.gate <- apply(exprs(trans.ff)[,all.gates], 2, mad)
        lower.bound <- med.per.gate - 3 * mad.per.gate
        lower.bound[!threshold.above] <- -Inf
        upper.bound <- med.per.gate + 3 * mad.per.gate
        upper.bound[threshold.above] <- Inf
        toignore <- c(toignore, all.gates)
    } else {
        all.gates <- NULL        
    }

    # Running through the files, subsetting by the specified gates, and saving the FCS files.
    # This allows us to use the exact same transformed intensities for tools that requires FCS input.
    dir.create(ref.out.dir, showWarning=FALSE)
    fcs.dir <- file.path(ref.out.dir, "fcs") 
    dir.create(fcs.dir)
    saved <- list()
    for (x in sampleNames(trans.nfs)) {
        cur.fcs <- trans.nfs[[x]]
        if (!is.null(all.gates)) {
            tex <- t(exprs(cur.fcs)[,all.gates,drop=FALSE])
            stats <- tex < lower.bound | tex > upper.bound
            saved[[x]] <- rowMeans(stats)
            keep <- colSums(stats)==0L
            cur.fcs <- cur.fcs[keep,]
        }
        description(cur.fcs)["transformation"] <- "logicle" # otherwise read.FCS doesn't like it.
        write.FCS(cur.fcs, file.path(fcs.dir, x))
    }
    new.fcs <- file.path(fcs.dir, basename(fcs))
    if (length(saved)) {
        saved <- do.call(rbind, saved)
        write.table(saved, file=file.path(ref.out.dir, "gate_saved.tsv"), sep="\t", quote=FALSE, col.names=NA)
    }

    ########################################################
    # Cleaning out the memory to avoid overruns.
    
    rm(trans.ff, replaced.ff, trans.nfs, new.nfs, new.exprs)
    gc()
} 


