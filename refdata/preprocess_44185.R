#This performs logicle-transformation and gating on intensities for each data set. 

library(ncdfFlow)
library(cydar)
host.dir <- "."

for (dataset in 1:5) {
    ref.dir <- "Cytobank_44185"
    extra <- paste0("H", dataset)

    ref.out.dir <- paste0(ref.dir, "_", extra)
    pattern <- sprintf("^%s_NoDrug_(Basal1|IL-[0-9]+)_.*", extra)
    
    # Reading in the files.
    fcs <- list.files(file.path(host.dir, ref.dir), full=TRUE, pattern=pattern)
    nfs <- read.ncdfFlowSet(fcs[1])
    descriptions <- as.character(parameters(nfs[[1]])$desc)
    toignore <- c(match(c("Time", "Cell_length", "PhenoGraph"), descriptions), grep("^BC[0-9]", descriptions))
    
    # Choosing which ones to gate on.
    gate.up <- match(c("DNA1", "DNA2"), descriptions)
    gate.down <- match(c("Viability"), descriptions)

    ########################################################
    # Estimating the parameters of the logicle transformation, after pooling all samples.
    
    new.exprs <- list()
    new.nfs <- suppressWarnings(suppressMessages(read.ncdfFlowSet(fcs)))
    replaced.ff <- poolCells(new.nfs, equalize=FALSE) # set FALSE for legacy purposes only.
    lgcl <- estimateLogicle(replaced.ff, channels=colnames(replaced.ff)[-toignore], type='data', m=5) # m=5 to avoid error.
    trans.ff <- transform(replaced.ff, lgcl)
    trans.nfs <- transform(new.nfs, lgcl)
   
    ########################################################
    # Need to gate on specified aspects (gate.up, gate.down).

    all.names <- colnames(trans.ff)
    gate.dna1 <- outlierGate(trans.ff, all.names[gate.up[1]], type="lower")
    gate.dna2 <- outlierGate(trans.ff, all.names[gate.up[2]], type="lower")
    gate.viab <- outlierGate(trans.ff, all.names[gate.down], type="upper")

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
        cur.fcs <- Subset(cur.fcs, gate.dna1)
        remaining$DNA1 <- nrow(cur.fcs)
        cur.fcs <- Subset(cur.fcs, gate.dna2)
        remaining$DNA2 <- nrow(cur.fcs)
        cur.fcs <- Subset(cur.fcs, gate.viab)
        remaining$Viability <- nrow(cur.fcs)
        saved[[x]] <- c(Total=remaining$Total, -diff(unlist(remaining))/remaining$Total)     
       
        description(cur.fcs)["transformation"] <- "logicle" # otherwise read.FCS doesn't like it.
        write.FCS(cur.fcs, file.path(fcs.dir, x))
    }
    new.fcs <- file.path(fcs.dir, basename(fcs))
        
    saved <- do.call(rbind, saved)
    write.table(saved, file=file.path(ref.out.dir, "gate_saved.tsv"), sep="\t", quote=FALSE, col.names=NA)

    ########################################################
    # Cleaning out the memory to avoid overruns.
    
    rm(trans.ff, replaced.ff, trans.nfs, new.nfs, new.exprs)
    gc()
} 


