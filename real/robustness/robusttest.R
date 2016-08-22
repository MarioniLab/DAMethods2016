# Loading data and counts.
require(cydar)
require(ncdfFlow)

output.all <- output.da <- output.propA <- output.propB <- list()
for (dataset in c("Cytobank_43324_4FI", "Cytobank_43324_NG", "Cytobank_43324_NN")) {
    cd <- readRDS(file.path("../../refdata", paste0(dataset, "_raw.rds")))

    res.scenarios <- list()
    all.collected <- list()
    for (scenario in 1:5) {
        tol <- 0.5
        nn <- 50
        if (scenario==2L) {
            tol <- 0.4
        } else if (scenario==3L) {
            tol <- 0.6
        } else if (scenario==4L) {
            nn <- 20
        } else if (scenario==5L) {
            nn <- 100
        }

        # counting cells and saving them for posterity.
        out <- countCells(cd, BPPARAM=SerialParam(), downsample=10, tol=tol)

        # Setting up the design matrix.
        timings <- as.integer(sub(".*_([0-9]+).fcs", "\\1", colnames(out$counts)))
        design <- model.matrix(~splines::ns(timings, 3))

        # Testing for differential proportions across time, using the spline.
        require(edgeR)
        y <- DGEList(out$counts, lib.size=out$total, genes=out$coordinates)
        keep <- aveLogCPM(y) >= aveLogCPM(5, mean(out$total))
        y <- y[keep,]

        y <- estimateDisp(y, design)
        fit <- glmQLFit(y, design, robust=TRUE)
        res <- glmQLFTest(fit, coef=2:ncol(design))
        qvals <- spatialFDR(y$genes, res$table$PValue, neighbors=nn)

        all.collected[[scenario]] <- rownames(out$coordinates)[keep]
        res.scenarios[[scenario]] <- rownames(out$coordinates)[keep][qvals <= 0.05]
    }

    output.all[[dataset]] <- lengths(all.collected)
    output.da[[dataset]] <- lengths(res.scenarios)
    output.propA[[dataset]] <- c(length(setdiff(res.scenarios[[2]], res.scenarios[[1]])),
                                 length(setdiff(res.scenarios[[3]], res.scenarios[[1]])),
                                 length(setdiff(res.scenarios[[4]], res.scenarios[[1]])),
                                 length(setdiff(res.scenarios[[5]], res.scenarios[[1]])))
    output.propB[[dataset]] <- c(length(setdiff(res.scenarios[[1]], res.scenarios[[2]])),
                                 length(setdiff(res.scenarios[[1]], res.scenarios[[3]])),
                                 length(setdiff(res.scenarios[[1]], res.scenarios[[4]])),
                                 length(setdiff(res.scenarios[[1]], res.scenarios[[5]])))
}

# Writing results to file, in Latex table format.
full.names <- c(Cytobank_43324_4FI="Oct4-GFP", 
                Cytobank_43324_NG="Nanog-GFP", 
                Cytobank_43324_NN="Nanog-Neo")

ofile <- "effects.txt"
con <- file(ofile, open="w")
for (dataset in names(output.all)) {
    x <- paste(c(paste(c(full.names[[dataset]], "Total", output.all[[dataset]]), collapse=" & "),
      paste(c("", "Significant", output.da[[dataset]]), collapse=" & "),
      paste(c("", "Gained", "-", output.propA[[dataset]]), collapse=" & "),
      paste(c("", "Lost", "-", output.propB[[dataset]]), collapse=" & ")), "\\\\")
    x <- c(x, "\\hline")
    writeLines(x, con=con)
}
close(con)

