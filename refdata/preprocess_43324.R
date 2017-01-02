# This performs logicle-transformation and gating on intensities for each data set. 

library(ncdfFlow)
library(cydar)

host.dir <- "."
for (dataset in 1:3) {
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
    gate.up <- c("DNA-1(Ir191)Dd", "DNA-2(Ir193)Dd")

    ########################################################
    # Estimating the parameters of the logicle transformation, after pooling all samples.
    new.exprs <- list()
    new.nfs <- suppressWarnings(suppressMessages(read.ncdfFlowSet(fcs)))
    replaced.ff <- poolCells(new.nfs, equalize=FALSE) # for legacy reasons, probably should use equalize=TRUE in future.
    lgcl <- estimateLogicle(replaced.ff, channels=colnames(replaced.ff)[-toignore], type='data')
    trans.ff <- transform(replaced.ff, lgcl)
    trans.nfs <- transform(new.nfs, lgcl)
   
    ########################################################
    # Need to gate on specified aspects (gate.up, gate.down).

    gate.191 <- outlierGate(trans.ff, gate.up[1], type="lower")
    gate.193 <- outlierGate(trans.ff, gate.up[2], type="lower")

    # Running through the files, subsetting by the specified gates, and saving the FCS files.
    # This allows us to use the exact same transformed intensities for tools that requires FCS input.
    dir.create(ref.out.dir, showWarning=FALSE)
    fcs.dir <- file.path(ref.out.dir, "fcs") 
    dir.create(fcs.dir)
    saved <- list()

    for (x in sampleNames(trans.nfs)) {
        cur.fcs <- trans.nfs[[x]]
        remaining <- list()
        remaining$Total <- nrow(cur.fcs)
        cur.fcs <- Subset(cur.fcs, gate.191)
        remaining$Gate.191 <- nrow(cur.fcs)
        cur.fcs <- Subset(cur.fcs, gate.193)
        remaining$Gate.193 <- nrow(cur.fcs)
        saved[[x]] <- c(Total=remaining$Total, -diff(unlist(remaining))/remaining$Total)     
        description(cur.fcs)["transformation"] <- "logicle" # otherwise read.FCS doesn't like it.
        write.FCS(cur.fcs, file.path(fcs.dir, x))
    }
    new.fcs <- file.path(fcs.dir, basename(fcs))
        
    saved <- data.frame(do.call(rbind, saved))
    write.table(saved, file=file.path(ref.out.dir, "gate_saved.tsv"), sep="\t", quote=FALSE, col.names=NA)

    ########################################################
    # Cleaning out the memory to avoid overruns.
    
    rm(trans.ff, replaced.ff, trans.nfs, new.nfs, new.exprs)
    gc()
} 


