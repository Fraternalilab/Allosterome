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
library("PSICQUIC");
library("parallel");

## number of cores
nCore = detectCores();

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
## PDBs
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
##   only once
#for (i in 1:length(pdb.pdb.unique.v)) {
#	get.pdb(file = pdb.pdb.unique.v[[i]], path = "../pdb");
#}

## extract FASTA sequences
##   only once
for (i in 1:length(pdb.pdb.unique.v)) {
	print(paste(i, pdb.pdb.unique.v[i]));
	if (pdb.pdb.unique.v[i] != "NULL") {
		pdbIn = read.pdb(file = paste("../pdb/", pdb.pdb.unique.v[i], ".pdb", sep = ""));
		fastaOut = pdbseq(pdbIn);
		write.fasta(seqs = fastaOut, ids = pdb.pdb.unique.v[[i]],
					file = paste("../pdb/", pdb.pdb.unique.v[[i]], ".fasta", sep = ""));
	}
}

#_______________________________________________________________________________
## pairwise alignment of all FASTA sequences
pdb.fasta.unique.v = paste("../pdb/", pdb.pdb.unique.v, ".fasta", sep = "");

## aligns two FASTA sequences and returns sequence ID
fastAlign = function(fasta.v, fasta1, fasta2, ali.m) {
	fastaIn1 = read.fasta(file = fasta.v[fasta1]);
	fastaIn2 = read.fasta(file = fasta.v[fasta2]);
	## bind FASTA inputs to single object
	fastaIn12 = seqbind(fastaIn1, fastaIn2);
	## align
	ali12 = seqaln(aln = fastaIn12);
	## assign sequence identity to result matrix
	ali.m[fasta1, fasta2] = seqidentity(ali12)[1, 2];
}

seqpair = combn(length(pdb.fasta.unique.v), 2);
seqpair.l = apply(seqpair, 2, as.list);
seqpair.short.l = seqpair.l[1:100];
aliresult.m = matrix(0, nrow = length(pdb.fasta.unique.v), ncol = length(pdb.fasta.unique.v));

## initiate cluster for parallel computation 
clu = makeCluster(nCore);
## make parallel functions see predefined variables
clusterExport(clu, c("fastAlign", "read.fasta", "seqbind", "seqaln", "seqidentity", "pdb.fasta.unique.v", "seqpair.short.l", "aliresult.m"));

lapply(seqpair.short.l, function(x) fastAlign(pdb.fasta.unique.v,
							as.numeric(unlist(x[1])),
							as.numeric(unlist(x[2])),
							aliresult.m));
#parLapply(clu, seqpair.short.l, function(x) fastAlign(pdb.fasta.unique.v, x[1], x[2], aliresult.m));
## save results
saveRDS(ali.m, "ali.RDS");
## release memory
stopCluster(clu);


#_______________________________________________________________________________
## effectors
#_______________________________________________________________________________
## get KEGG identifiers
kegg.l = lapply(datalist.l, function(x) {x$KEGG_ID[[1]]});
## as vector; empty identifiers adopt "NULL"
kegg.v = as.character(kegg.l);
kegg.unique.v = unique(kegg.v);

#_______________________________________________________________________________
## download KEGG structures
## we use Bioconductor's PSICQUIC package as database interface
psicquic = PSICQUIC();



#_______________________________________________________________________________
save.image(file = "xml_analyse.RData");

#===============================================================================

