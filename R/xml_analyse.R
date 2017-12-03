#! /usr/bin/R

#===============================================================================
# process XML entries
# Copyright (C) 2017 Jens Kleinjung
# Dependencies are the libraries listed below, for 'pathview' run on the shell:
# source("http://bioconductor.org/biocLite.R")
# biocLite("pathview")
#===============================================================================

library("xml2");
library("bio3d");
library("pathview");

#_______________________________________________________________________________
## path of XML files
inPath = "../ASD/ASD_Release_062015_XF";
## character vector of all XML files (database entries)
filenames = list.files(path = inPath, full.names = TRUE, pattern = 'xml$');
## read contents of all XML files
datalist = lapply(filenames, function(x){read_xml(x)});
## convert contents to lists, making it a list (entries) of lists (contents)
## takes a minute
datalist.l = lapply(datalist, as_list);

#_______________________________________________________________________________
## get PDB identifiers
pdb.l = lapply(datalist.l, function(x) {x$PDB_List$PDB$PDB_ID[[1]]});
## as character vector; empty identifiers adopt "NULL"
pdb.v = as.character(pdb.l);
## ASD models
pdb.mod.v = pdb.v[substring(pdb.v, 1, 3) == "ASD"];
pdb.mod.unique.v = unique(pdb.mod.v);
## PDB models
pdb.pdb.v = pdb.v[substring(pdb.v, 1, 3) != "ASD"];
pdb.pdb.unique.v = unique(pdb.pdb.v);

#_______________________________________________________________________________
## download PDB structures
for (i in 1:length(pdb.pdb.unique.v)) {
	get.pdb(pdb.pdb.unique.v[[i]], path = "../pdb");
}

#_______________________________________________________________________________
## get KEGG identifiers
kegg.l = lapply(datalist.l, function(x) {x$KEGG_ID[[1]]});
## as vector; empty identifiers adopt "NULL"
kegg.v = as.character(kegg.l);
kegg.unique.v = unique(kegg.v);

#_______________________________________________________________________________
## download KEGG structures


#_______________________________________________________________________________
save.image();

#_______________________________________________________________________________
#_______________________________________________________________________________

