pdf("pics/dispersions.pdf")
plot(0, 0, ylab="NB dispersion", xlab=expression(Log[2]*"-average count"), type="n", xlim=c(2, 10), ylim=c(0.2, 1.8), cex.axis=1.2, cex.lab=1.4) 
colours <- heat.colors(3)
datasets <- c("Cytobank_43324_4FI", "Cytobank_43324_NG", "Cytobank_43324_NN")

for (i in seq_along(datasets)) { 
    dataset <- datasets[i]
    data <- readRDS(file.path("../../refdata", paste0(dataset, "_counts.rds")))
    res <- readRDS(paste0(dataset, "_res.rds"))
    ab <- res$results$AveLogCPM + log2(mean(data$totals)/1e6)
    o <- order(ab)

    # Telling us what it should be.
    if (dataset=="Cytobank_43324_4FI") {
        transfect <- "Oct4-GFP"
    } else if (dataset=="Cytobank_43324_NG") {
        transfect <- "Nanog-GFP"
    } else if (dataset=="Cytobank_43324_NN") {
        transfect <- "Nanog-Neo"
    }

    lines(ab[o], res$results$Dispersion[o], col="black", lwd=5)
    lines(ab[o], res$results$Dispersion[o], col=colours[i], lwd=1.5)
    j <- which.max(ab)
    text(ab[j], res$results$Dispersion[j], transfect, pos=4, cex=1.2) 
}
dev.off()
