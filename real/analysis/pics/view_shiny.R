stuff <- readRDS("../Cytobank_43324_4FI_res.rds")
cd <- readRDS("../../../refdata/Cytobank_43324_4FI.rds")
coords <- read.table("../Cytobank_43324_4FI_coords.txt")

require(cydar)
nonredundant <- findFirstSphere(stuff$coords, stuff$results$PValue)
is.sig <- stuff$results$FDR <= 0.05

labels <- character(sum(is.sig))
labels[which(nonredundant)[1:3]] <- c("Mouse embryonic fibroblasts",
                                      "SC4-like cells, IdU+",
                                      "ESC-like cells")
app <- interpretSpheres(cd[stuff$kept[is.sig],], select=nonredundant[is.sig], 
                        red.coords=coords, run=TRUE, labels=labels)

# Hit 'next' a couple of times, zoom out and save a full-page webshot to 'shiny.png'.
