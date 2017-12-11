#! /usr/bin/R

## implementation of the PSICQUIC vignette

library(PSICQUIC);

psicquic = PSICQUIC();
providers(psicquic);

tbl = interactions(psicquic, id=c("TP53", "MYC"), species="9606");
dim(tbl);
table(tbl$provider);

tbl[ , c("provider", "type", "detectionMethod")];

tbl[grep("affinity", tbl$detectionMethod),
	c("type", "publicationID", "firstAuthor", "confidenceScore", "provider")];

tbl.myc = interactions(psicquic, "MYC", species = "9606", publicationID = "21150319");
dim(tbl.myc);
table(tbl.myc$provider);
table(tbl.myc$confidenceScore);

idMapper = IDMapper("9606");
tbl.myc = addGeneInfo(idMapper,tbl.myc);
print(head(tbl.myc$A.name));
print(head(tbl.myc$B.name));

tbl.3 = interactions(psicquic, id=c("ALK", "JAK3", "SHC3"), species="9606", quiet=TRUE);
tbl.3g = addGeneInfo(idMapper, tbl.3);
tbl.3gd = with(tbl.3g, as.data.frame(table(detectionMethod, type, A.name, B.name, provider)));
print(tbl.3gd = subset(tbl.3gd, Freq > 0));

