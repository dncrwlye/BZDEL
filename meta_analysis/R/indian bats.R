batphy_for_rotl_update <- read_csv("~/Desktop/BDEL/BZDEL/meta_analysis/data/batphy_for_rotl update.csv")

text2 <- stri_extract_all_regex(text[1], "\\([A-Z][a-z]+ [a-z]+\\)")
text2 <- as.data.frame(text2)
text2[,2] <- gsub("\\(","", text2[,1])
text2[,3] <- gsub("\\)","", text2[,2])

text2$V3 =plyr::revalue(text2$V3,c("Rousettus egyptiacus"="Rousettus aegyptiacus",
                                   "Rousettus leschenaulti"="Rousettus leschenaultii",
                                   "Myotis blythii"="Myotis oxygnathus",
                                   "Pipistrellus affinis"="Falsistrellus affinis",
                                   "Pipistrellus cadornae" = "Hypsugo savii",
                                   "Pipistrellus savii" = "Hypsugo savii",
                                   "Tadarida plicata" = 'Chaerephon plicatus'))
                                   
text2$V3 <- gsub(" ","_", text2$V3)
