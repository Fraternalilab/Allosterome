#! /usr/bin/R

#===============================================================================
# process XML entries
# Copytight (C) 2017 Jens Kleinjung
#===============================================================================

library("xml2");

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
## get PDB entries
pdb.l = lapply(datalist.l, function(x) {x$PDB_List$PDB$PDB_ID[[1]]});
## as vector (entry 1473 has no PDB identifier),
##   therefore 'pdb.v' is one element shorter than 'pdb.l'
pdb.v = unlist(pdb.l);
## ASD models
pdb.mod.v = pdb.v[substring(pdb.v, 1, 3) == "ASD"];
pdb.mod.unique.v = unique(pdb.mod.v);
## PDB models
pdb.pdb.v = pdb.v[substring(pdb.v, 1, 3) != "ASD"];
pdb.pdb.unique.v = unique(pdb.pdb.v);


#_______________________________________________________________________________
save.image();

#_______________________________________________________________________________
#_______________________________________________________________________________

