stuff <- readRDS("../Cytobank_43324_4FI_res.rds")
cd <- readRDS("../../../refdata/Cytobank_43324_4FI_raw.rds")

require(cydar)
nonredundant <- findFirstSphere(stuff$coords, stuff$results$PValue)
labels <- character(sum(nonredundant))
labels[1] <- "Mouse embryonic fibroblasts"
labels[2] <- "SC4-like cells, IdU+"
labels[3] <- "ESC-like cells"
app <- interpretSpheres(stuff$coords[nonredundant,], cd, run=TRUE, labels=labels)

# Hit 'next' a couple of times, zoom out and save a screenshot to 'shiny.png'.
