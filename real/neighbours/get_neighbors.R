require(ncdfFlow)
require(FNN)
set.seed(100)

firstin <- FALSE
for (dataset in c("Cytobank_43324_4FI", "Cytobank_43324_NG", "Cytobank_43324_NN")) {
    ref.dir <- file.path("../../refdata", dataset, "fcs")
    x <- read.ncdfFlowSet(list.files(ref.dir, full=TRUE))
    toignore <- c(match(c("Cell_length", "Time", "beadDist", "barcode"), colnames(x)), grep("^BC-[0-9]", colnames(x)), grep("DNA", colnames(x)))
    x <- x[,-toignore]

    # Telling us what it should be.
    if (dataset=="Cytobank_43324_4FI") {
        transfect <- "Oct4-GFP"
        nn <- 10
    } else if (dataset=="Cytobank_43324_NG") {
        transfect <- "Nanog-GFP"
        nn <- 1
    } else if (dataset=="Cytobank_43324_NN") {
        transfect <- "Nanog-Neo"
        nn <- 10
    }

    # Computing the NN distances for each cell in each sample.
    current <- exprs(x[[1]])
    chosen <- nn # c(5, 10, 20, 40)
    stuff <- get.knn(current, k=max(chosen))
    stuff$nn.dist <- stuff$nn.dist[,chosen,drop=FALSE]
    stuff$nn.index <- stuff$nn.index[,chosen,drop=FALSE]

    output.loc <- output.err <- list()
    nmarkers <- c(9, 16, 25, 36)
    for (y in nmarkers) {
        if (y < ncol(current)) { 
            xcur <- current[,sample(ncol(current), y)]
            collected <- list()
            for (i in seq_len(nrow(current))) {
                collected[[i]] <- sqrt(colSums((xcur[i,] - t(xcur[stuff$nn.index[i,],,drop=FALSE]))^2))
            }
            collected <- do.call(rbind, collected)
        } else {
            collected <- stuff$nn.dist
        }

        i <- length(output.loc)+1L
        output.loc[[i]] <- apply(collected, 2, mean)
        output.err[[i]] <- apply(collected, 2, sd) 
        # Using standard deviation, not standard error.
        # This is because we want to know how it varies across cells.
    }

    output.loc <- do.call(rbind, output.loc)
    output.err <- do.call(rbind, output.err)

    pdf(paste0(dataset, ".pdf"))
    plot(0,0, ylab="Euclidean distance", xlab="Number of markers", 
         cex.axis=1.2, cex.lab=1.4, ylim=c(0, 4), type="n", 
         xlim=c(5, 40), xaxt="n", main=transfect,cex.main=1.4)
    axis(1, at=nmarkers, label=nmarkers, cex.axis=1.2)

    hot.col <- "grey" # rev(grey.colors(length(nmarkers)))
    offsets <- 0 # c(-0.75, -0.25, 0.25, 0.75)
    for (i in seq_len(ncol(output.loc))) {
        xpos <- nmarkers+offsets[i]
        ypos <- output.loc[,i]
        yplus <- ypos + output.err[,i]
        yneg <- ypos - output.err[,i]
        segments(xpos, yneg, xpos, yplus)
        segments(xpos-0.1, yplus, xpos+0.1, yplus)
        segments(xpos-0.1, yneg, xpos+0.1, yneg)
        points(xpos, ypos, bg=hot.col[i], pch=21)
    }

    curve(sqrt(x)*0.5, add=TRUE, col="red", lwd=2)

#    if (!firstin) {
#        mytext <- paste(chosen, "neighbours")
#        legend(28, 1, pch=21, legend=mytext, pt.bg=hot.col, cex=1.2)
#        firstin <- TRUE
#    }
    dev.off()
}

########################################################
# A similar plot.

library(ncdfFlow)
dataset <- "Cytobank_43324_4FI"
ref.dir <- file.path("../../refdata", dataset, "fcs")
x <- read.ncdfFlowSet(list.files(ref.dir, full=TRUE))
toignore <- c(match(c("Cell_length", "Time", "beadDist", "barcode"), colnames(x)), grep("^BC-[0-9]", colnames(x)), grep("DNA", colnames(x)))
x <- x[,-toignore]

library(cydar)
cd <- prepareCellData(x)
distances <- neighborDistances(cd)
distances <- cbind(0, distances) # adding itself as part of the count.

pdf("nvd_Cytobank_43324_4FI.pdf")
boxplot(distances, horizontal=TRUE, xlab=expression(r[0]*sqrt(2)), yaxt="n", ylim=c(0, 0.9), ylab="Cell count", cex.axis=1.2, cex.lab=1.4, outline=FALSE)
axis(side=2, at=pretty(c(0, ncol(distances))), cex=1.2)
abline(v=0.5, col="red", lwd=2, lty=2)
dev.off()

