##########################################################################
#                  NECESSARY SETTINGS & DEPENDENCIES                     #
##########################################################################

options(scipen = 999)

library(ggvenn)
library(readxl)
library(gplots)



##########################################################################
#                     FUNCTION DEFINITIONS                               #
##########################################################################

source(file="library/ddz-pd-phospho-lib.r")



##########################################################################
#                               MAIN                                     #
##########################################################################

cwd <- "E:/Projekte/Phospho2.1-pd3.0/output_phospho2.1-workflow"

prot.dat.orig <- as.data.frame(read_excel("E:/Projekte/Phospho2.1-pd3.0/Phospho2.1_Quan_with_isoforms-(1)_proteins.xlsx"))
dat.orig <- as.data.frame(read_excel("E:/Projekte/Phospho2.1-pd3.0/Phospho2.1_Quan_with_isoforms-(1)_phosphopeptides.xlsx"))

p.orig <- nrow(dat.orig)

# remove peptides without protein identification
# --> relevant biological information not available
# (no phosphosite localization, kinase substrate, pathways, ...)
if(length(which(is.na(dat.orig$`Master Protein Accessions`))) > 0) dat.orig <- dat.orig[-which(is.na(dat.orig$`Master Protein Accessions`)),]



#-------|
# NOTE: |
#-------|
# --> ENSURE THAT PEPTIDES ARE ORDERED BY SEQUENCE
dat.orig <- dat.orig[order(dat.orig$`Annotated Sequence`),]



accession.list <- unique(dat.orig[,"Master Protein Accessions"])
row.IDs <- vector(mode="character", length=length(accession.list))
for(i in 1:length(accession.list)){
  #idx <- grep(accession.list[i], dat.orig[,"Master Protein Accessions"], fixed=TRUE)
  idx <- grep(paste0("^",accession.list[i],"$"), dat.orig[,"Master Protein Accessions"])
  n.peps <- length(idx)
  row.IDs[idx] <- paste0(accession.list[i], "_peptide", 1:n.peps)
}
rownames(dat.orig) <- row.IDs



#-------|
# NOTE: |
#-------|
# nrow(dat.orig) = 13228
# length(unmatched.idx) = 0
# length(contamin.idx) = 32
# length(intersect(unmatched.idx, contamin.idx)) = 0
# --> nrow(dat.orig[-c(contamin.idx, unmatched.idx),]) = 13196 (not less since contamin = 0)
unmatched.idx <- giveUnmatchedPeps(pep.names=rownames(dat.orig), prot.vec=prot.dat.orig[,"Accession"])
contamin.idx <- giveContaminPeps(pep.dat=dat.orig)
if(length(unmatched.idx) > 0 || length(contamin.idx) > 0){
    dat.orig <- dat.orig[-unique(c(contamin.idx, unmatched.idx)),]
    row.IDs <- row.IDs[-unique(c(contamin.idx, unmatched.idx))]
}



#col.abundance.reps <- grep("^Abundance: F\\d+:", colnames(dat.orig), value=TRUE)
col.abundance.reps <- grep("^Abundances \\(Normalized", colnames(dat.orig), value=TRUE)
dat.abundance.reps <- data.matrix(dat.orig[,col.abundance.reps])
col.abundance.gr <- grep("^Abundances \\(Grouped\\):", colnames(dat.orig), value=TRUE)
dat.abundance.gr <- data.matrix(dat.orig[,col.abundance.gr])
col.abundance <- grep("^Abundances \\(by Bio", colnames(dat.orig), value=TRUE)
dat.abundance <- data.matrix(dat.orig[,col.abundance])
col.ratio <- grep("Abundance Ratio\\:", colnames(dat.orig), value=TRUE)
dat.ratio <- data.matrix(dat.orig[,col.ratio])
col.p.values <- grep("Abundance Ratio P-Value", colnames(dat.orig), value=TRUE)
p.values <- data.matrix(dat.orig[,col.p.values])
col.p.values.adj <- grep("Abundance Ratio Adj. P-Value", colnames(dat.orig), value=TRUE)
p.values.adj <- data.matrix(dat.orig[,col.p.values.adj])
col.mods <- grep("^Modifications$", colnames(dat.orig), value=TRUE)
mods <- dat.orig[,col.mods]
col.seq <- grep("Sequence", colnames(dat.orig), value=TRUE)
seq <- dat.orig[,col.seq]
col.pos <- grep("^Positions in Master", colnames(dat.orig), value=TRUE)
pos <- dat.orig[,col.pos]
col.mods.master <- grep("^Modifications in Master Proteins$", colnames(dat.orig), value=TRUE)
mods.master <- dat.orig[,col.mods.master]
col.missed.cleav <- grep("# Missed Cleavages", colnames(dat.orig), value=TRUE)
missed.cleav <- dat.orig[,col.missed.cleav]
col.psm.number <- grep("# PSMs", colnames(dat.orig), value=TRUE)
psm.number <- dat.orig[,col.psm.number]



rownames(dat.abundance.reps) <- row.IDs
rownames(dat.abundance.gr) <- row.IDs
rownames(dat.abundance) <- row.IDs
rownames(dat.ratio) <- row.IDs
rownames(p.values) <- row.IDs
rownames(p.values.adj) <- row.IDs
names(mods) <- row.IDs
names(seq) <- row.IDs
names(pos) <- row.IDs
names(mods.master) <- row.IDs
names(missed.cleav) <- row.IDs
names(psm.number) <- row.IDs



colnames(dat.abundance) <- gsub("^Abundances \\(by Bio.+\\: (\\d+\\.*\\d*)\\,([A-Z1-9]+)", "\\2_\\1", colnames(dat.abundance))



#p <- nrow(dat.orig)
timepoint.lab <- c("1", "2.5", "5", "15", "30", "60")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Match information from protein-centric PD-results file to |
# the measured phosphopeptides from the peptide-centric PD- |
# results file.                                             |
#------------------------------------------------------------

prot.info <- prepareProteinInfo(pep.names=rownames(dat.orig), prot.info=prot.dat.orig)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Create lists of protein accessions from current dataset |
# e.g. for UniProt ID-mapping service.                    |
#----------------------------------------------------------

mapPepToProt(pep.ids=rownames(dat.orig), iso.rm=TRUE, group.rm=TRUE, output.path=cwd)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Filtering NAs to have continuous time course data.  |
#------------------------------------------------------

na.idx_noNA <- c()
for(i in 1:nrow(dat.abundance.gr)){
  if(sum(is.na(dat.abundance.gr[i,2:7])) > 0){
    na.idx_noNA <- c(na.idx_noNA,i)  
  }
}

dat.abundance.reps.used <- dat.abundance.reps[-(na.idx_noNA),]
dat.abundance.gr.used <- dat.abundance.gr[-(na.idx_noNA),]
dat.abundance.used <- dat.abundance[-(na.idx_noNA),]
dat.ratio.used <- dat.ratio[-(na.idx_noNA),]
p.values.used <- p.values[-(na.idx_noNA),]
p.values.adj.used <- p.values.adj[-(na.idx_noNA),]
mods.used <- mods[-(na.idx_noNA)]
seq.used <- seq[-(na.idx_noNA)]
pos.used <- pos[-(na.idx_noNA)]
mods.master.used <- mods.master[-(na.idx_noNA)]
missed.cleav.used <- missed.cleav[-(na.idx_noNA)]
psm.number.used <- psm.number[-(na.idx_noNA)]

timepnt.n <- ncol(dat.ratio.used)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Write txt files containing different pieces of the data. |
#----------------------------------------------------------

output.tmp <- cbind(rownames(dat.abundance), dat.abundance)
colnames(output.tmp) <- c("Peptide ID", colnames(dat.abundance))
write.table(x=output.tmp, file=paste0(cwd, "/dat.abundance.txt"), col.names=T, row.names=F, quote=F, sep="\t")

output.tmp <- cbind(rownames(dat.abundance), p.values, dat.ratio, dat.abundance)
colnames(output.tmp) <- c("Peptide ID", colnames(p.values), colnames(dat.ratio), colnames(dat.abundance))
write.table(x=output.tmp, file=paste0(cwd, "/dat.abundance.extended.txt"), col.names=T, row.names=F, quote=F, sep="\t")

output.tmp <- cbind(rownames(dat.abundance.used), dat.abundance.used)
colnames(output.tmp) <- c("Peptide ID", colnames(dat.abundance.used))
write.table(x=output.tmp, file=paste0(cwd, "/dat.abundance.used.txt"), col.names=T, row.names=F, quote=F, sep="\t")

output.tmp <- cbind(rownames(dat.abundance.used), p.values.used, dat.ratio.used, dat.abundance.used)
colnames(output.tmp) <- c("Peptide ID", colnames(p.values.used), colnames(dat.ratio.used), colnames(dat.abundance.used))
write.table(x=output.tmp, file=paste0(cwd, "/dat.abundance.extended.used.txt"), col.names=T, row.names=F, quote=F, sep="\t")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Write Excel file containing information on    |
# phosphopeptides to the specified output path. |
#------------------------------------------------

preparePeptideInfo(
  pep.dat=dat.ratio.used,
  pep.cols1=cbind(seq.used, mods.used),
  prot.info=prot.info,
  pep.cols2=cbind(p.values.used, p.values.adj.used),
  output.path=cwd
)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Prepare regulated phosphopeptides for further |
# analysis steps.                               |
#------------------------------------------------

reg.peps <- list()
reg.peps.adj <- list()
reg.peps.herwig <- list()
reg.peps.pos <- list()
reg.peps.neg <- list()
reg.peps.idx <- list()
reg.peps.idx.pos <- list()
reg.peps.idx.neg <- list()
reg.peps.adj.idx <- list()
reg.peps.n.vec <- vector(mode="numeric", length=timepnt.n)
reg.peps.pos.n.vec <- vector(mode="numeric", length=timepnt.n)
reg.peps.neg.n.vec <- vector(mode="numeric", length=timepnt.n)
reg.peps.adj.n.vec <- vector(mode="numeric", length=timepnt.n)
reg.peps.all <- c()
reg.peps.all.adj <- c()
reg.peps.all.herwig <- c()
reg.direction.matrix <- matrix(0, nrow=nrow(dat.ratio.used), ncol=timepnt.n)
rownames(reg.direction.matrix) <- rownames(dat.ratio.used)
colnames <- timepoint.lab

# same as above but including "inter-time point" NAs
reg.peps.na <- list()
reg.peps.pos.na <- list()
reg.peps.neg.na <- list()
reg.peps.idx.na <- list()
reg.peps.n.vec.na <- vector(mode="numeric", length=timepnt.n)
reg.peps.pos.n.vec.na <- vector(mode="numeric", length=timepnt.n)
reg.peps.neg.n.vec.na <- vector(mode="numeric", length=timepnt.n)
reg.peps.all.na <- c()
reg.peps.all.pos.na <- c()
reg.peps.all.neg.na <- c()

for(i in 1:timepnt.n){
  reg.peps[[i]] <- regulatedPeptides(ratios=dat.ratio.used[,i], p.values=p.values.used[,i], direction="both", output.type="row.names")
  reg.peps.pos[[i]] <- regulatedPeptides(ratios=dat.ratio.used[,i], p.values=p.values.used[,i], direction="positive", output.type="row.names")
  reg.peps.neg[[i]] <- regulatedPeptides(ratios=dat.ratio.used[,i], p.values=p.values.used[,i], direction="negative", output.type="row.names")
  reg.peps.adj[[i]] <- regulatedPeptides(ratios=dat.ratio.used[,i], p.values=p.values.adj.used[,i], direction="both", output.type="row.names")
  reg.peps.herwig[[i]] <- regulatedPeptides(ratios=dat.ratio.used[,i], p.values=p.values.adj.used[,i], direction="herwig", output.type="row.names")
  #next row: hartwig list for comparison with total proteome
  #reg.peps.herwig[[i]] <- regulatedPeptides(ratios=dat.ratio.used[,i], p.values=p.values.used[,i], direction="both", output.type="row.names")
  
  reg.peps.idx[[i]] <- regulatedPeptides(ratios=dat.ratio.used[,i], p.values=p.values.used[,i], direction="both", output.type="boolean")
  reg.peps.idx.pos[[i]] <- regulatedPeptides(ratios=dat.ratio.used[,i], p.values=p.values.used[,i], direction="positive", output.type="boolean")
  reg.peps.idx.neg[[i]] <- regulatedPeptides(ratios=dat.ratio.used[,i], p.values=p.values.used[,i], direction="negative", output.type="boolean")
  
  reg.peps.adj.idx[[i]] <- regulatedPeptides(ratios=dat.ratio.used[,i], p.values=p.values.adj.used[,i], direction="both", output.type="boolean")
  
  reg.direction.matrix[reg.peps.idx.pos[[i]],i] <- 1
  reg.direction.matrix[reg.peps.idx.neg[[i]],i] <- -1 
  
  reg.peps.n.vec[i] <- length(reg.peps[[i]])
  reg.peps.adj.n.vec[i] <- length(reg.peps.adj[[i]])
  reg.peps.pos.n.vec[i] <- length(which(reg.peps.idx.pos[[i]]))
  reg.peps.neg.n.vec[i] <- length(which(reg.peps.idx.neg[[i]]))
  reg.peps.all <- c(reg.peps.all, reg.peps[[i]])
  reg.peps.all.adj <- c(reg.peps.all.adj, reg.peps.adj[[i]])
  reg.peps.all.herwig <- c(reg.peps.all.herwig, reg.peps.herwig[[i]])
  
  
  
  rows.tmp <- names(na.omit(dat.ratio[,i]))
  reg.peps.na[[i]] <- regulatedPeptides(ratios=dat.ratio[rows.tmp,i], p.values=p.values[rows.tmp,i], direction="both", output.type="row.names")
  reg.peps.pos.na[[i]] <- regulatedPeptides(ratios=dat.ratio[rows.tmp,i], p.values=p.values[rows.tmp,i], direction="positive", output.type="row.names")
  reg.peps.neg.na[[i]] <- regulatedPeptides(ratios=dat.ratio[rows.tmp,i], p.values=p.values[rows.tmp,i], direction="negative", output.type="row.names")
  reg.peps.idx.na[[i]] <- regulatedPeptides(ratios=dat.ratio[rows.tmp,i], p.values=p.values[rows.tmp,i], direction="both", output.type="boolean")
  reg.peps.n.vec.na[i] <- length(reg.peps.na[[i]])
  reg.peps.pos.n.vec.na[i] <- length(reg.peps.pos.na[[i]])
  reg.peps.neg.n.vec.na[i] <- length(reg.peps.neg.na[[i]])
  reg.peps.all.na <- c(reg.peps.all.na, reg.peps.na[[i]])
  reg.peps.all.pos.na <- c(reg.peps.all.pos.na, reg.peps.pos.na[[i]])
  reg.peps.all.neg.na <- c(reg.peps.all.neg.na, reg.peps.neg.na[[i]])
}
barplot(reg.peps.adj.n.vec)
barplot(reg.peps.n.vec)
barplot(reg.peps.n.vec.na)
reg.peps.all <- unique(reg.peps.all)
reg.peps.all.adj <- unique(reg.peps.all.adj)
reg.peps.all.herwig <- unique(reg.peps.all.herwig)
print(paste0("Number of unique reg. peptides across all time points: ", length(reg.peps.all)))
print(paste0("Number of unique reg. peptides across all time points (adjusted): ", length(reg.peps.all.adj)))
print(paste0("Number of unique reg. peptides across all time points (herwig): ", length(reg.peps.all.herwig)))

reg.peps.all.na <- unique(reg.peps.all.na)
reg.peps.all.pos.na <- unique(reg.peps.all.pos.na)
reg.peps.all.neg.na <- unique(reg.peps.all.neg.na)
print(paste0("Number of unique reg. peptides across all time points (incl. inter-tp NAs): ", length(reg.peps.all.na)))
print(paste0("Number of unique upreg. peptides across all time points (incl. inter-tp NAs): ", length(reg.peps.all.pos.na)))
print(paste0("Number of unique downreg. peptides across all time points (incl. inter-tp NAs): ", length(reg.peps.all.neg.na)))
print(paste0("Control - unique(c(upreg,downreg): ", length(unique(c(reg.peps.all.pos.na, reg.peps.all.neg.na)))))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Count phosphosites on identified/quantified/localized |
# phosphopeptides.                                      |
#--------------------------------------------------------

# TODO: Write output to file!?

phossite.dat.QOT.LOF <- phosphoPepsToSites(dat=dat.abundance.gr,
                                         mods.master=mods.master,
                                         mods=mods, psm.number=psm.number,
                                         reg.peps=NULL,
                                         quantified.only=T,
                                         localized.only=F,
                                         output.path=paste0(cwd, "/phosphoPepsToSites"),
                                         file.name="phossite.dat.QOT.LOF")
nrow(phossite.dat.QOT.LOF)

#---

phossite.dat.QOT.LOT <- phosphoPepsToSites(dat=dat.abundance.gr,
                                         mods.master=mods.master,
                                         mods=mods, psm.number=psm.number,
                                         reg.peps=NULL,
                                         quantified.only=T,
                                         localized.only=T,
                                         output.path=paste0(cwd, "/phosphoPepsToSites"),
                                         file.name="phossite.dat.QOT.LOT")
nrow(phossite.dat.QOT.LOT)

#---

phossite.dat.QOT.LOT.reg.all <- phosphoPepsToSites(dat=dat.abundance.gr,
                                         mods.master=mods.master,
                                         mods=mods, psm.number=psm.number,
                                         reg.peps=reg.peps.all.na,
                                         quantified.only=T,
                                         localized.only=T,
                                         output.path=paste0(cwd, "/phosphoPepsToSites"),
                                         file.name="phossite.dat.QOT.LOT.reg.all")
nrow(phossite.dat.QOT.LOT.reg.all)

#-------------------------------------------------------------------------------

phossite.dat.QOT.LOF.1 <- phosphoPepsToSites(dat=dat.abundance.gr,
                                           mods.master=mods.master,
                                           mods=mods, psm.number=psm.number,
                                           reg.peps=reg.peps.na[[1]],
                                           quantified.only=T,
                                           localized.only=F,
                                           output.path=paste0(cwd, "/phosphoPepsToSites"),
                                           file.name="phossite.dat.QOT.LOF.1")
phossite.dat.QOT.LOF.2 <- phosphoPepsToSites(dat=dat.abundance.gr,
                                           mods.master=mods.master,
                                           mods=mods, psm.number=psm.number,
                                           reg.peps=reg.peps.na[[2]],
                                           quantified.only=T,
                                           localized.only=F,
                                           output.path=paste0(cwd, "/phosphoPepsToSites"),
                                           file.name="phossite.dat.QOT.LOF.2")
phossite.dat.QOT.LOF.3 <- phosphoPepsToSites(dat=dat.abundance.gr,
                                           mods.master=mods.master,
                                           mods=mods, psm.number=psm.number,
                                           reg.peps=reg.peps.na[[3]],
                                           quantified.only=T,
                                           localized.only=F,
                                           output.path=paste0(cwd, "/phosphoPepsToSites"),
                                           file.name="phossite.dat.QOT.LOF.3")
phossite.dat.QOT.LOF.4 <- phosphoPepsToSites(dat=dat.abundance.gr,
                                           mods.master=mods.master,
                                           mods=mods, psm.number=psm.number,
                                           reg.peps=reg.peps.na[[4]],
                                           quantified.only=T,
                                           localized.only=F,
                                           output.path=paste0(cwd, "/phosphoPepsToSites"),
                                           file.name="phossite.dat.QOT.LOF.4")
phossite.dat.QOT.LOF.5 <- phosphoPepsToSites(dat=dat.abundance.gr,
                                           mods.master=mods.master,
                                           mods=mods, psm.number=psm.number,
                                           reg.peps=reg.peps.na[[5]],
                                           quantified.only=T,
                                           localized.only=F,
                                           output.path=paste0(cwd, "/phosphoPepsToSites"),
                                           file.name="phossite.dat.QOT.LOF.5")
phossite.dat.QOT.LOF.6 <- phosphoPepsToSites(dat=dat.abundance.gr,
                                           mods.master=mods.master,
                                           mods=mods, psm.number=psm.number,
                                           reg.peps=reg.peps.na[[6]],
                                           quantified.only=T,
                                           localized.only=F,
                                           output.path=paste0(cwd, "/phosphoPepsToSites"),
                                           file.name="phossite.dat.QOT.LOF.6")

nrow(phossite.dat.QOT.LOF.1)
nrow(phossite.dat.QOT.LOF.2)
nrow(phossite.dat.QOT.LOF.3)
nrow(phossite.dat.QOT.LOF.4)
nrow(phossite.dat.QOT.LOF.5)
nrow(phossite.dat.QOT.LOF.6)
#---
length(intersect(rownames(phossite.dat.QOT.LOF.1), rownames(phossite.dat.QOT.LOF.2)))
length(intersect(rownames(phossite.dat.QOT.LOF.2), rownames(phossite.dat.QOT.LOF.3)))
length(intersect(rownames(phossite.dat.QOT.LOF.3), rownames(phossite.dat.QOT.LOF.4)))
length(intersect(rownames(phossite.dat.QOT.LOF.4), rownames(phossite.dat.QOT.LOF.5)))
length(intersect(rownames(phossite.dat.QOT.LOF.5), rownames(phossite.dat.QOT.LOF.6)))
#---
length(intersect(rownames(phossite.dat.QOT.LOF.1), rownames(phossite.dat.QOT.LOF.6)))
#---
phos.site.tmp <- unique(c(rownames(phossite.dat.QOT.LOF.2),
                          rownames(phossite.dat.QOT.LOF.3),
                          rownames(phossite.dat.QOT.LOF.4),
                          rownames(phossite.dat.QOT.LOF.5),
                          rownames(phossite.dat.QOT.LOF.6))) 
length(setdiff(rownames(phossite.dat.QOT.LOF.1), phos.site.tmp))
#---
phos.site.tmp <- unique(c(rownames(phossite.dat.QOT.LOF.1),
                          rownames(phossite.dat.QOT.LOF.3),
                          rownames(phossite.dat.QOT.LOF.4),
                          rownames(phossite.dat.QOT.LOF.5),
                          rownames(phossite.dat.QOT.LOF.6))) 
length(setdiff(rownames(phossite.dat.QOT.LOF.2), phos.site.tmp))
#---
phos.site.tmp <- unique(c(rownames(phossite.dat.QOT.LOF.1),
                          rownames(phossite.dat.QOT.LOF.2),
                          rownames(phossite.dat.QOT.LOF.4),
                          rownames(phossite.dat.QOT.LOF.5),
                          rownames(phossite.dat.QOT.LOF.6))) 
length(setdiff(rownames(phossite.dat.QOT.LOF.3), phos.site.tmp))
#---
phos.site.tmp <- unique(c(rownames(phossite.dat.QOT.LOF.1),
                          rownames(phossite.dat.QOT.LOF.2),
                          rownames(phossite.dat.QOT.LOF.3),
                          rownames(phossite.dat.QOT.LOF.5),
                          rownames(phossite.dat.QOT.LOF.6))) 
length(setdiff(rownames(phossite.dat.QOT.LOF.4), phos.site.tmp))
#---
phos.site.tmp <- unique(c(rownames(phossite.dat.QOT.LOF.1),
                          rownames(phossite.dat.QOT.LOF.2),
                          rownames(phossite.dat.QOT.LOF.3),
                          rownames(phossite.dat.QOT.LOF.4),
                          rownames(phossite.dat.QOT.LOF.6))) 
length(setdiff(rownames(phossite.dat.QOT.LOF.5), phos.site.tmp))
#---
phos.site.tmp <- unique(c(rownames(phossite.dat.QOT.LOF.1),
                          rownames(phossite.dat.QOT.LOF.2),
                          rownames(phossite.dat.QOT.LOF.3),
                          rownames(phossite.dat.QOT.LOF.4),
                          rownames(phossite.dat.QOT.LOF.5))) 
length(setdiff(rownames(phossite.dat.QOT.LOF.6), phos.site.tmp))

#-------------------------------------------------------------------------------

phossite.dat.QOT.LOT.1 <- phosphoPepsToSites(dat=dat.abundance.gr,
                                           mods.master=mods.master,
                                           mods=mods, psm.number=psm.number,
                                           reg.peps=reg.peps.na[[1]],
                                           quantified.only=T,
                                           localized.only=T,
                                           output.path=paste0(cwd, "/phosphoPepsToSites"),
                                           file.name="phossite.dat.QOT.LOT.1")
phossite.dat.QOT.LOT.2 <- phosphoPepsToSites(dat=dat.abundance.gr,
                                           mods.master=mods.master,
                                           mods=mods, psm.number=psm.number,
                                           reg.peps=reg.peps.na[[2]],
                                           quantified.only=T,
                                           localized.only=T,
                                           output.path=paste0(cwd, "/phosphoPepsToSites"),
                                           file.name="phossite.dat.QOT.LOT.2")
phossite.dat.QOT.LOT.3 <- phosphoPepsToSites(dat=dat.abundance.gr,
                                           mods.master=mods.master,
                                           mods=mods, psm.number=psm.number,
                                           reg.peps=reg.peps.na[[3]],
                                           quantified.only=T,
                                           localized.only=T,
                                           output.path=paste0(cwd, "/phosphoPepsToSites"),
                                           file.name="phossite.dat.QOT.LOT.3")
phossite.dat.QOT.LOT.4 <- phosphoPepsToSites(dat=dat.abundance.gr,
                                           mods.master=mods.master,
                                           mods=mods, psm.number=psm.number,
                                           reg.peps=reg.peps.na[[4]],
                                           quantified.only=T,
                                           localized.only=T,
                                           output.path=paste0(cwd, "/phosphoPepsToSites"),
                                           file.name="phossite.dat.QOT.LOT.4")
phossite.dat.QOT.LOT.5 <- phosphoPepsToSites(dat=dat.abundance.gr,
                                           mods.master=mods.master,
                                           mods=mods, psm.number=psm.number,
                                           reg.peps=reg.peps.na[[5]],
                                           quantified.only=T,
                                           localized.only=T,
                                           output.path=paste0(cwd, "/phosphoPepsToSites"),
                                           file.name="phossite.dat.QOT.LOT.5")
phossite.dat.QOT.LOT.6 <- phosphoPepsToSites(dat=dat.abundance.gr,
                                           mods.master=mods.master,
                                           mods=mods, psm.number=psm.number,
                                           reg.peps=reg.peps.na[[6]],
                                           quantified.only=T,
                                           localized.only=T,
                                           output.path=paste0(cwd, "/phosphoPepsToSites"),
                                           file.name="phossite.dat.QOT.LOT.6")
nrow(phossite.dat.QOT.LOT.1)
nrow(phossite.dat.QOT.LOT.2)
nrow(phossite.dat.QOT.LOT.3)
nrow(phossite.dat.QOT.LOT.4)
nrow(phossite.dat.QOT.LOT.5)
nrow(phossite.dat.QOT.LOT.6)
#---
length(intersect(rownames(phossite.dat.QOT.LOT.1), rownames(phossite.dat.QOT.LOT.2)))
length(intersect(rownames(phossite.dat.QOT.LOT.2), rownames(phossite.dat.QOT.LOT.3)))
length(intersect(rownames(phossite.dat.QOT.LOT.3), rownames(phossite.dat.QOT.LOT.4)))
length(intersect(rownames(phossite.dat.QOT.LOT.4), rownames(phossite.dat.QOT.LOT.5)))
length(intersect(rownames(phossite.dat.QOT.LOT.5), rownames(phossite.dat.QOT.LOT.6)))
#---
length(intersect(rownames(phossite.dat.QOT.LOT.1), rownames(phossite.dat.QOT.LOT.6)))
#---







phos.site.tmp <- unique(c(rownames(phossite.dat.QOT.LOT.2),
                          rownames(phossite.dat.QOT.LOT.3),
                          rownames(phossite.dat.QOT.LOT.4),
                          rownames(phossite.dat.QOT.LOT.5),
                          rownames(phossite.dat.QOT.LOT.6))) 
length(setdiff(rownames(phossite.dat.QOT.LOT.1), phos.site.tmp))
grep("Q96B36", phos.site.tmp, value=T)
#---
phos.site.tmp <- unique(c(rownames(phossite.dat.QOT.LOT.1),
                          rownames(phossite.dat.QOT.LOT.3),
                          rownames(phossite.dat.QOT.LOT.4),
                          rownames(phossite.dat.QOT.LOT.5),
                          rownames(phossite.dat.QOT.LOT.6))) 
length(setdiff(rownames(phossite.dat.QOT.LOT.2), phos.site.tmp))
grep("Q96B36", phos.site.tmp, value=T)
#---
phos.site.tmp <- unique(c(rownames(phossite.dat.QOT.LOT.1),
                          rownames(phossite.dat.QOT.LOT.2),
                          rownames(phossite.dat.QOT.LOT.4),
                          rownames(phossite.dat.QOT.LOT.5),
                          rownames(phossite.dat.QOT.LOT.6))) 
length(setdiff(rownames(phossite.dat.QOT.LOT.3), phos.site.tmp))
grep("Q96B36", phos.site.tmp, value=T)
#---
phos.site.tmp <- unique(c(rownames(phossite.dat.QOT.LOT.1),
                          rownames(phossite.dat.QOT.LOT.2),
                          rownames(phossite.dat.QOT.LOT.3),
                          rownames(phossite.dat.QOT.LOT.5),
                          rownames(phossite.dat.QOT.LOT.6))) 
length(setdiff(rownames(phossite.dat.QOT.LOT.4), phos.site.tmp))
grep("Q96B36", phos.site.tmp, value=T)
#---
phos.site.tmp <- unique(c(rownames(phossite.dat.QOT.LOT.1),
                          rownames(phossite.dat.QOT.LOT.2),
                          rownames(phossite.dat.QOT.LOT.3),
                          rownames(phossite.dat.QOT.LOT.4),
                          rownames(phossite.dat.QOT.LOT.6))) 
length(setdiff(rownames(phossite.dat.QOT.LOT.5), phos.site.tmp))
grep("Q96B36", phos.site.tmp, value=T)
#---
phos.site.tmp <- unique(c(rownames(phossite.dat.QOT.LOT.1),
                          rownames(phossite.dat.QOT.LOT.2),
                          rownames(phossite.dat.QOT.LOT.3),
                          rownames(phossite.dat.QOT.LOT.4),
                          rownames(phossite.dat.QOT.LOT.5))) 
length(setdiff(rownames(phossite.dat.QOT.LOT.6), phos.site.tmp))
grep("Q96B36", phos.site.tmp, value=T)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

plotProfiles(x=dat.abundance.used,
             info=dat.orig[rownames(dat.abundance.used),c("Modifications","Modifications in Master Proteins","# PSMs","# Missed Cleavages")],
             prot.info=prot.info,
             plot=TRUE,
             output.path=paste0(cwd, "/Profiles"))

#mono.id.vec <- grep("1xPhospho", dat.orig[rownames(dat.abundance.used),"Modifications in Master Proteins"], value=F)
#plotProfiles(x=dat.abundance.used[mono.id.vec,],
plotProfiles(x=dat.abundance.used[1:100,],
#plotProfiles(x=dat.abundance.used,
             info=dat.orig[rownames(dat.abundance.used),c("Modifications","Modifications in Master Proteins","# PSMs","# Missed Cleavages")],
             prot.info=prot.info,
             plot=FALSE,
             output.path=paste0(cwd, "/Profiles_test"))

#library(stringr)
#str_count(pattern="1xPhospho", string=dat.orig[1:100,"Modifications in Master Proteins"])
mono.id.vec <- grep("1xPhospho", dat.orig[rownames(dat.abundance.used),"Modifications in Master Proteins"], value=F)
plotProfiles(x=dat.abundance.used[mono.id.vec,],
             info=dat.orig[,c("Modifications","Modifications in Master Proteins","# PSMs","# Missed Cleavages")],
             prot.info=prot.info,
             plot=TRUE,
             output.path=paste0(cwd, "/Profiles_mono_only"))

SupplFig6.id.vec <- c("Q07352_peptide3","P31751; P31749; Q9Y243_peptide1","P62753_peptide4","P23588_peptide11")
plotProfiles(x=dat.abundance.used[SupplFig6.id.vec,],
             info=dat.orig[,c("Modifications","Modifications in Master Proteins","# PSMs","# Missed Cleavages")],
             prot.info=prot.info,
             plot=TRUE,
             anonymize.ids=TRUE,
             output.path=paste0(cwd, "/Profiles_SupplFig4_only"))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Perform identification of stable sub-   |
# network & correlation plot for NetCore. |
#------------------------------------------

# ----> TODO: stable (netcore) subnet function.

# ----> TODO: correlaltion plot function.



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Compile & save 'phospho_herwig list' of |
# phosphopeptides & proteins.             |                         |
#------------------------------------------

pep.first.cols <- c("Master Protein Accessions",
                    "Master Protein Descriptions",
                    "Modifications",
                    "Modifications in Master Proteins",
                    "Protein Accessions")
pep.last.cols <- c("# Proteins",
                   "Annotated Sequence",
                   "# PSMs")
phospep.herwig <- cbind(
                    reg.peps.all.herwig,
                    dat.orig[reg.peps.all.herwig, pep.first.cols],
                    dat.ratio.used[reg.peps.all.herwig,],
                    p.values.used[reg.peps.all.herwig,],
                    dat.orig[reg.peps.all.herwig, pep.last.cols])
colnames(phospep.herwig) <- c("Peptide ID",
                              pep.first.cols,
                              colnames(dat.ratio.used),
                              colnames(p.values.used),
                              pep.last.cols)
write.table(x=phospep.herwig,
            file=paste0(cwd, "/herwig_phosphopeptidegroups.txt"),
            row.names=F,
            col.names=T,
            sep="\t")
  
prot.herwig <- prepareProteinInfo(pep.names=reg.peps.all.herwig, prot.info=prot.dat.orig)
write.table(x=prot.herwig,
            file=paste0(cwd, "/herwig_proteins.txt"),
            row.names=F,
            col.names=T,
            sep="\t")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Check replicate reproducibility.     |
#---------------------------------------

donor.ids <- unique(gsub("_.+", "", colnames(dat.abundance.used)))
sample.ids <- c()
for(i in 1:length(donor.ids)) {
  sample.ids <- c(sample.ids, paste(c("Sample, 0","Sample, 1","Sample, 2.5","Sample, 5","Sample, 15","Sample, 30","Sample, 60"),
                                    donor.ids[i], sep=", "))
}
assessReproducibility(x=dat.abundance.reps.used, sample.ids=sample.ids, output.path=cwd)

# time-patterns -> 0, 1, 2.5, 5, 15, 30, 60
#sample.ids2 <- grep(", 60,", colnames(dat.abundance.reps.used), value=T)
# sample-patterns -> AR7, AT8, AT13, AT18, AT22
#sample.ids2 <- grep(", AT22", colnames(dat.abundance.reps.used), value=T)
sample.ids2 <- colnames(dat.abundance.reps.used)
# -> cex.lab=0.55 or 0.85
assessReproducibility2(x=dat.abundance.reps.used[,sample.ids2], sample.ids=sample.ids2, cex.lab=0.55, output.path=cwd)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Supl. Figure X (NAs & contaminants). |
#---------------------------------------

col00 <- vector(mode="numeric", length=timepnt.n)
col0 <- vector(mode="numeric", length=timepnt.n)
col1 <- vector(mode="numeric", length=timepnt.n)
for(i in 1:timepnt.n) col00[i] <- length(dat.ratio[,i]) #number of all phos.peps
for(i in 1:timepnt.n) col0[i] <- length(na.omit(dat.ratio[,i])) #number of phos.peps after time point-specific NA filter
for(i in 1:timepnt.n) col1[i] <- length(dat.ratio.used[,i]) #number of phos.peps after across time NA filter

na.tab <- rbind(col1, col0-col1, col00-col0)
colnames (na.tab) <-  c("1 min", "2.5 min", "5 min", "15 min", "30 min", "60 min") 
rownames (na.tab) <- c('#Phosphopeptides.filtered.across.time.points', "#NAs.across.time.points", "#NAs.within.time.point")

removed.n <- p.orig - nrow(dat.ratio.used)

png(filename=paste0(cwd, "/number_NAs.png"), height=2000, width=2000, res=300)
    #mar = c(bottom, left, top, right)
    par(mar = c(5.1, 5.1, 4.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
    xx <- barplot(
              na.tab[1:3,],
              #border="black",
              #space=0.04,
              col = c("grey30", "grey60", "grey75"),
              font.axis=2,
              font.lab=2,
              main = paste0("Removing ", removed.n, " phosphopeptides"),
              ylim=c(0,16500),
              cex.axis=1.1,
              cex.lab=1.5,
              cex.main=1.5,
              cex.names=1.2,
              las=1
    )
    title(ylab="#Phosphopeptides", mgp=c(4,1,0),font.lab=2,cex.lab=1.5)
    text(x=xx, y=na.tab[1,], label=na.tab[1,], pos=1, cex=1.15, font=2, col="white")
    text(x=xx, y=(na.tab[1,]+na.tab[2,]+na.tab[3,]), label=(na.tab[1,]+na.tab[2,]+na.tab[3,]), pos=3, cex=1.15, font=2, col="black")

    legend("top",
       legend = c("Intra-time\npoint NAs", "Inter-time\npoint NAs", "Kept\npeptides"),
       col = c("grey75", "grey60", "grey30"),
       pch = 15,
       x.intersp = 1.1,
       bty = "n",
       text.width = 1.25,
       text.font = 2,
       pt.cex = 2.5,
       cex = 1,
       ncol = 3
    )
dev.off()

png(filename=paste0(cwd, "/number_NAs_noTitle.png"), height=2000, width=2000, res=300)
  #mar = c(bottom, left, top, right)
  par(mar = c(5.1, 5.1, 4.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
  xx <- barplot(
    na.tab[1:3,],
    col = c("grey30", "grey60", "grey75"),
    font.axis=2,
    font.lab=2,
    ylim=c(0,16500),
    cex.axis=1.1,
    cex.lab=1.5,
    cex.main=1.5,
    cex.names=1.2,
    las=1
  )
  title(ylab="#Phosphopeptides", mgp=c(4,1,0),font.lab=2,cex.lab=1.5)
  text(x=xx, y=na.tab[1,], label=na.tab[1,], pos=1, cex=1.15, font=2, col="white")
  text(x=xx, y=(na.tab[1,]+na.tab[2,]+na.tab[3,]), label=(na.tab[1,]+na.tab[2,]+na.tab[3,]), pos=3, cex=1.15, font=2, col="black")

  legend("top",
       legend = c("Intra-time\npoint NAs", "Inter-time\npoint NAs", "Kept\npeptides"),
       col = c("grey75", "grey60", "grey30"),
       pch = 15,
       x.intersp = 1.1,
       bty = "n",
       text.width = 1.25,
       text.font = 2,
       pt.cex = 2.5,
       cex = 1,
       ncol = 3
  )
dev.off()



png(filename=paste0(cwd, "/number_NAs_zoom.png "), height=2000, width=2000, res=300)
    #mar = c(bottom, left, top, right)
    par(mar = c(5.1, 5.1, 4.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
    xx <- barplot(
              na.tab[1:3,],
              #border="black",
              #space=0.04,
              col = c("grey30", "grey60", "grey75"),
              font.axis=2,
              font.lab=2,
              main = paste0("Removing ", removed.n, " phosphopeptides"),
              ylim=c(10000,14000),
              cex.axis=1.1,
              cex.lab=1.5,
              cex.main=1.5,
              cex.names=1.2,
              las=1
              )
    title(ylab="#Phosphopeptides", mgp=c(4,1,0),font.lab=2,cex.lab=1.5)
    text(x=xx, y=na.tab[1,], label=na.tab[1,], pos=1, cex=1.15, font=2, col="white")
    text(x=xx, y=(na.tab[1,]+na.tab[2,]), label=na.tab[2,], pos=1, cex=1.3, font=2, col="white")
    text(x=xx, y=(na.tab[1,]+na.tab[2,]+na.tab[3,]), label=na.tab[3,], pos=1, cex=1.3, font=2, col="white")
    text(x=xx, y=(na.tab[1,]+na.tab[2,]+na.tab[3,]), label=(na.tab[1,]+na.tab[2,]+na.tab[3,]), pos=3, cex=1.15, font=2, col="black")
    

    legend("top",
           legend = c("Intra-time\npoint NAs", "Inter-time\npoint NAs", "Kept\npeptides"),
           col = c("grey75", "grey60", "grey30"),
           pch = 15,
           x.intersp = 1.1,
           bty = "n",
           text.width = 1.25,
           text.font = 2,
           pt.cex = 2.5,
           cex = 1,
           ncol = 3
    )
dev.off()

png(filename=paste0(cwd, "/number_NAs_zoom_noTitle.png "), height=2000, width=2000, res=300)
  #mar = c(bottom, left, top, right)
  par(mar = c(5.1, 5.1, 4.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
  xx <- barplot(
    na.tab[1:3,],
    col = c("grey30", "grey60", "grey75"),
    font.axis=2,
    font.lab=2,
    ylim=c(10000,14000),
    cex.axis=1.1,
    cex.lab=1.5,
    cex.main=1.5,
    cex.names=1.2,
    las=1
  )
  title(ylab="#Phosphopeptides", mgp=c(4,1,0),font.lab=2,cex.lab=1.5)
  text(x=xx, y=na.tab[1,], label=na.tab[1,], pos=1, cex=1.15, font=2, col="white")
  text(x=xx, y=(na.tab[1,]+na.tab[2,]), label=na.tab[2,], pos=1, cex=1.3, font=2, col="white")
  text(x=xx, y=(na.tab[1,]+na.tab[2,]+na.tab[3,]), label=na.tab[3,], pos=1, cex=1.3, font=2, col="white")
  text(x=xx, y=(na.tab[1,]+na.tab[2,]+na.tab[3,]), label=(na.tab[1,]+na.tab[2,]+na.tab[3,]), pos=3, cex=1.15, font=2, col="black")


  legend("top",
       legend = c("Intra-time\npoint NAs", "Inter-time\npoint NAs", "Kept\npeptides"),
       col = c("grey75", "grey60", "grey30"),
       pch = 15,
       x.intersp = 1.1,
       bty = "n",
       text.width = 1.25,
       text.font = 2,
       pt.cex = 2.5,
       cex = 1,
       ncol = 3
  )
dev.off()

#---------
#---------

na.tab2 <- rbind(col0, col00-col0)
colnames (na.tab2) <-  c("1 min", "2.5 min", "5 min", "15 min", "30 min", "60 min") 
rownames (na.tab2) <- c('#Phosphopeptides.filtered.within.time.points', "#NAs.within.time.point")

png(filename=paste0(cwd, "/number_NAs_na.png"), height=2000, width=2000, res=300)
  #mar = c(bottom, left, top, right)
  par(mar = c(5.1, 5.1, 4.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
  xx <- barplot(
    na.tab2[1:2,],
    #border="black",
    #space=0.04,
    col = c("grey30", "grey75"),
    font.axis=2,
    font.lab=2,
    main = paste0("Removing time point-specific NAs"),
    ylim=c(0,16000),
    cex.axis=1.1,
    cex.lab=1.5,
    cex.main=1.5,
    cex.names=1.2,
    las=1
  )
  title(ylab="#Phosphopeptides", mgp=c(4,1,0),font.lab=2,cex.lab=1.5)
  text(x=xx, y=na.tab2[1,], label=na.tab2[1,], pos=1, cex=1.15, font=2, col="white")
  text(x=xx, y=(na.tab2[1,]+na.tab2[2,]), label=(na.tab2[1,]+na.tab2[2,]), pos=3, cex=1.15, font=2, col="black")

  legend("top",
       legend = c("Time point-specific NAs", "Kept peptides"),
       col = c("grey75", "grey30"),
       pch = 15,
       x.intersp = 1.1,
       bty = "n",
       text.width = 2.7,
       text.font = 2,
       pt.cex = 2.5,
       cex = 1,
       horiz=T
  )
dev.off()

png(filename=paste0(cwd, "/Fig_S1a_number_NAs_noTitle_na.png"), height=2000, width=2000, res=300)
  #mar = c(bottom, left, top, right)
  par(mar = c(5.1, 5.1, 4.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
  xx <- barplot(
    na.tab2[1:2,],
    #border="black",
    #space=0.04,
    col = c("grey30", "grey75"),
    font.axis=2,
    font.lab=2,
    ylim=c(0,16000),
    cex.axis=1.1,
    cex.lab=1.5,
    cex.main=1.5,
    cex.names=1.2,
    las=1
  )
  title(ylab="#Phosphopeptides", mgp=c(4,1,0),font.lab=2,cex.lab=1.5)
  text(x=xx, y=na.tab2[1,], label=na.tab2[1,], pos=1, cex=1.15, font=2, col="white")
  text(x=xx, y=(na.tab2[1,]+na.tab2[2,]), label=(na.tab2[1,]+na.tab2[2,]), pos=3, cex=1.15, font=2, col="black")

  legend("top",
       legend = c("Time point-specific NAs", "Kept peptides"),
       col = c("grey75", "grey30"),
       pch = 15,
       x.intersp = 1.1,
       bty = "n",
       text.width = 2.7,
       text.font = 2,
       pt.cex = 2.5,
       cex = 1,
       horiz=T
  )
dev.off()

png(filename=paste0(cwd, "/number_NAs_zoom_na.png "), height=2000, width=2000, res=300)
  #mar = c(bottom, left, top, right)
  par(mar = c(5.1, 5.1, 4.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
  xx <- barplot(
    na.tab2[1:2,],
    #border="black",
    #space=0.04,
    col = c("grey30", "grey75"),
    font.axis=2,
    font.lab=2,
    main = paste0("Removing time point-specific NAs"),
    ylim=c(10000,14000),
    cex.axis=1.1,
    cex.lab=1.5,
    cex.main=1.5,
    cex.names=1.2,
    las=1
  )
  title(ylab="#Phosphopeptides", mgp=c(4,1,0),font.lab=2,cex.lab=1.5)
  text(x=xx, y=na.tab2[1,], label=na.tab2[1,], pos=1, cex=1.15, font=2, col="white")
  text(x=xx, y=(na.tab2[1,]+na.tab2[2,]), label=na.tab2[2,], pos=1, cex=1.3, font=2, col="white")
  text(x=xx, y=(na.tab2[1,]+na.tab2[2,]), label=(na.tab2[1,]+na.tab2[2,]), pos=3, cex=1.15, font=2, col="black")
  
  
  legend("top",
       legend = c("Time point-specific NAs", "Kept peptides"),
       col = c("grey75", "grey30"),
       pch = 15,
       x.intersp = 1.1,
       bty = "n",
       text.width = 2.7,
       text.font = 2,
       pt.cex = 2.5,
       cex = 1,
       horiz=T
  )
dev.off()

png(filename=paste0(cwd, "/Fig_S1b_number_NAs_zoom_noTitle_na.png "), height=2000, width=2000, res=300)
#mar = c(bottom, left, top, right)
par(mar = c(5.1, 5.1, 4.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
xx <- barplot(
  na.tab2[1:2,],
  #border="black",
  #space=0.04,
  col = c("grey30", "grey75"),
  font.axis=2,
  font.lab=2,
  ylim=c(10000,14000),
  cex.axis=1.1,
  cex.lab=1.5,
  cex.main=1.5,
  cex.names=1.2,
  las=1
)
title(ylab="#Phosphopeptides", mgp=c(4,1,0),font.lab=2,cex.lab=1.5)
text(x=xx, y=na.tab2[1,], label=na.tab2[1,], pos=1, cex=1.15, font=2, col="white")
text(x=xx, y=(na.tab2[1,]+na.tab2[2,]), label=na.tab2[2,], pos=1, cex=1.3, font=2, col="white")
text(x=xx, y=(na.tab2[1,]+na.tab2[2,]), label=(na.tab2[1,]+na.tab2[2,]), pos=3, cex=1.15, font=2, col="black")


legend("top",
       legend = c("Time point-specific NAs", "Kept peptides"),
       col = c("grey75", "grey30"),
       pch = 15,
       x.intersp = 1.1,
       bty = "n",
       text.width = 2.7,
       text.font = 2,
       pt.cex = 2.5,
       cex = 1,
       horiz=T
)
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Supl. Figure X (Pie chart reg. vs. unreg. peptides) |
#----------------------------------------------------

total.number.reg.peptides <- length(reg.peps.all)
total.number.peptides <- nrow(dat.ratio.used)
total.number.unreg.peptides <- total.number.peptides - total.number.reg.peptides
png (paste0(cwd, "/pie_sig_phosphopeptides.png"), res=200, height = 1200, width=1200)
  par(lwd = 3)
  pie(
    c(total.number.reg.peptides, total.number.unreg.peptides),
    labels = c(paste0("Reg. phospho-\npeptides: ", total.number.reg.peptides), paste0("Unreg. phospho-  \npeptides: ", total.number.unreg.peptides)),
    col = c("grey50", "grey75"),
    border = "black",
    cex = 1.0,
    font = 2
  )
dev.off()

#---------

total.number.reg.peptides <- length(reg.peps.all.adj)
total.number.peptides <- nrow(dat.ratio.used)
total.number.unreg.peptides <- total.number.peptides - total.number.reg.peptides
png (paste0(cwd, "/pie_sig_phosphopeptides_adj.png"), res=200, height = 1200, width=1200)
  par(lwd = 3)
  pie(
    c(total.number.reg.peptides, total.number.unreg.peptides),
    labels = c(paste0("Reg. phospho-\npeptides: ", total.number.reg.peptides), paste0("Unreg. phospho-  \npeptides: ", total.number.unreg.peptides)),
    col = c("grey50", "grey75"),
    border = "black",
    cex = 1.0,
    font = 2
  )
dev.off()

#---------
#---------

total.number.reg.peptides <- length(reg.peps.all.na)
total.number.peptides <- nrow(dat.ratio)
total.number.unreg.peptides <- total.number.peptides - total.number.reg.peptides
png (paste0(cwd, "/Fig_S2_pie_sig_phosphopeptides_na.png"), res=200, height = 1200, width=1200)
par(lwd = 3)
  pie(
    c(total.number.reg.peptides, total.number.unreg.peptides),
    labels = c(paste0("Reg. phospho-\npeptides: ", total.number.reg.peptides), paste0("Unreg. phospho-  \npeptides: ", total.number.unreg.peptides)),
    col = c("grey50", "grey75"),
    border = "black",
    cex = 1.0,
    font = 2
  )
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Supl. Figure X (Pie chart proteins with vs. |
# proteins without reg. phosphopeptides)      |
#----------------------------------------------

total.reg.prots <- unique(gsub("\\_peptide\\d+", "", reg.peps.all))
prot.groups.idx <- grep("\\; ", total.reg.prots)
prot.groups <- grep("\\; ", total.reg.prots, value=TRUE)
splitted.accessions <- unlist(strsplit(grep("\\; ", prot.groups, value=TRUE), "\\; "))
if(length(prot.groups.idx) > 0){
  total.reg.prots <- total.reg.prots[-prot.groups.idx]
  total.reg.prots <- unique(c(total.reg.prots, splitted.accessions))
}

total.prots <- unique(gsub("\\_peptide\\d+", "", rownames(dat.ratio.used)))
prot.groups.idx <- grep("\\; ", total.prots)
prot.groups <- grep("\\; ", total.prots, value=TRUE)
splitted.accessions <- unlist(strsplit(grep("\\; ", prot.groups, value=TRUE), "\\; "))
if(length(prot.groups.idx) > 0){
  total.prots <- total.prots[-prot.groups.idx]
  total.prots <- unique(c(total.prots, splitted.accessions))
}

total.number.reg.prots <- length(total.reg.prots)
total.number.prots <- length(total.prots)
total.number.unreg.prots <- total.number.prots - total.number.reg.prots

png (paste0(cwd, "/pie_sig_phosphoproteins.png"), res=200, height = 1200, width=1200)
  par(lwd = 3)
  pie(
    c(total.number.reg.prots, total.number.unreg.prots),
    labels = c(paste0("Proteins with reg.\n    phosphopeptides: ", total.number.reg.prots, "       "), paste0("Proteins without reg.\nphosphopeptides: ", total.number.unreg.prots)),
    col = c("grey50", "grey75"),
    border = "black",
    cex = 1.0,
    font = 2
  )
dev.off()

#-----------------------

total.reg.prots <- unique(gsub("\\_peptide\\d+", "", reg.peps.all.adj))
prot.groups.idx <- grep("\\; ", total.reg.prots)
prot.groups <- grep("\\; ", total.reg.prots, value=TRUE)
splitted.accessions <- unlist(strsplit(grep("\\; ", prot.groups, value=TRUE), "\\; "))
if(length(prot.groups.idx) > 0){
  total.reg.prots <- total.reg.prots[-prot.groups.idx]
  total.reg.prots <- unique(c(total.reg.prots, splitted.accessions))
}

total.prots <- unique(gsub("\\_peptide\\d+", "", rownames(dat.ratio.used)))
prot.groups.idx <- grep("\\; ", total.prots)
prot.groups <- grep("\\; ", total.prots, value=TRUE)
splitted.accessions <- unlist(strsplit(grep("\\; ", prot.groups, value=TRUE), "\\; "))
if(length(prot.groups.idx) > 0){
  total.prots <- total.prots[-prot.groups.idx]
  total.prots <- unique(c(total.prots, splitted.accessions))
}

total.number.reg.prots <- length(total.reg.prots)
total.number.prots <- length(total.prots)
total.number.unreg.prots <- total.number.prots - total.number.reg.prots

png (paste0(cwd, "/pie_sig_phosphoproteins_adj.png"), res=200, height = 1200, width=1200)
  par(lwd = 3)
  pie(
    c(total.number.reg.prots, total.number.unreg.prots),
    labels = c(paste0("Proteins with reg.\nphosphopeptides: ", total.number.reg.prots, "       "), paste0("Proteins without reg.\nphosphopeptides: ", total.number.unreg.prots)),
    col = c("grey50", "grey75"),
    border = "black",
    cex = 1.0,
    font = 2
  )
dev.off()

#-----------------------

total.reg.prots <- unique(gsub("\\_peptide\\d+", "", reg.peps.all.na))
prot.groups.idx <- grep("\\; ", total.reg.prots)
prot.groups <- grep("\\; ", total.reg.prots, value=TRUE)
splitted.accessions <- unlist(strsplit(grep("\\; ", prot.groups, value=TRUE), "\\; "))
if(length(prot.groups.idx) > 0){
  total.reg.prots <- total.reg.prots[-prot.groups.idx]
  total.reg.prots <- unique(c(total.reg.prots, splitted.accessions))
}

total.prots <- unique(gsub("\\_peptide\\d+", "", rownames(dat.ratio)))
prot.groups.idx <- grep("\\; ", total.prots)
prot.groups <- grep("\\; ", total.prots, value=TRUE)
splitted.accessions <- unlist(strsplit(grep("\\; ", prot.groups, value=TRUE), "\\; "))
if(length(prot.groups.idx) > 0){
  total.prots <- total.prots[-prot.groups.idx]
  total.prots <- unique(c(total.prots, splitted.accessions))
}

total.number.reg.prots <- length(total.reg.prots)
total.number.prots <- length(total.prots)
total.number.unreg.prots <- total.number.prots - total.number.reg.prots

png (paste0(cwd, "/Fig_S3_pie_sig_phosphoproteins_na.png"), res=200, height = 1200, width=1200)
  par(lwd = 3)
  pie(
    c(total.number.reg.prots, total.number.unreg.prots),
    labels = c(paste0("Proteins with reg.\nphosphopeptides: ", total.number.reg.prots, "       "), paste0("  Proteins without reg.\nphosphopeptides: ", total.number.unreg.prots)),
    col = c("grey50", "grey75"),
    border = "black",
    cex = 1.0,
    font = 2
  )
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Figure 1 B) (Time course: number of reg. phosphopeptides). |
#-------------------------------------------------------------

col00 <- vector(mode="numeric", length=timepnt.n)
col0 <- vector(mode="numeric", length=timepnt.n)
col1 <- vector(mode="numeric", length=timepnt.n)
col2 <- vector(mode="numeric", length=timepnt.n)
for(i in 1:timepnt.n) col00[i] <- length(dat.ratio[,i]) #number of all phos.peps
for(i in 1:timepnt.n) col0[i] <- length(na.omit(dat.ratio[,i])) #number of phos.peps after time point-specific NA filter
for(i in 1:timepnt.n) col1[i] <- length(dat.ratio.used[,i]) #number of phos.peps after across time NA filter
for(i in 1:timepnt.n) col2[i] <- reg.peps.n.vec[i] #number of reg. phos.peps
                                                   

phospep.tab <- rbind (col1,col1-col2,col2)
colnames (phospep.tab) <-  c("1 min", "2.5 min", "5 min", "15 min", "30 min", "60 min") 
rownames (phospep.tab) <- c('#Phosphopeptides', "#Unreg.phosphopeptides", "#Sig.reg.phosphopeptides")

png(filename=paste0(cwd, "/number_regulated_phosphopeptides.png "), height=2000, width=2000, res=300)
  #mar = c(bottom, left, top, right)
  par(mar = c(5.1, 5.1, 4.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
  xx <- barplot(phospep.tab[2:3,],
              border="black",
              col = c("grey30", "grey75"),
              #space=0.04,
              font.axis=2,
              font.lab=2,
              main ="Regulated phosphopeptides\nin time course of insulin treatment",
              ylim=c(0,14000),
              cex.axis=1.1,
              cex.lab=1.5,
              cex.main=1.5,
              cex.names=1.2,
              las=1
  )
  title(ylab="#Phosphopeptides", mgp=c(4,1,0),font.lab=2,cex.lab=1.5)
  text(x=xx, y=phospep.tab[1,], label=phospep.tab[3,], pos=1, cex=1.5, font=2, col="white")
  text(x=xx, y=phospep.tab[1,], label=phospep.tab[1,], pos=3, cex=0.75, font=2, col="black")
  text(x=xx, y=phospep.tab[2,], label=phospep.tab[2,], pos=1, cex=0.75, font=2, col="white")
  legend("top",
       horiz = TRUE,
       legend = c("Regulated", "Unregulated"),
       col = c("grey75", "grey30"),
       pch = 15,
       bty = "n",
       pt.cex = 3,
       cex=1.5
  )
dev.off()

png(filename=paste0(cwd, "/number_regulated_phosphopeptides_noTitle.png "), height=2000, width=2000, res=300)
  #mar = c(bottom, left, top, right)
  par(mar = c(5.1, 5.1, 4.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
  xx <- barplot(phospep.tab[2:3,],
              border="black",
              col = c("grey30", "grey75"),
              #col = c("grey60", "red"),
              font.axis=2,
              font.lab=2,
              ylim=c(0,14499),
              cex.axis=1.1,
              cex.lab=1.5,
              cex.names=1.2,
              las=1
  )
  title(ylab="#Phosphopeptides", mgp=c(4,1,0),font.lab=2,cex.lab=1.5)
  text(x=xx, y=phospep.tab[1,], label=phospep.tab[3,], pos=3, cex=1.5, font=2, col="black")
  #text(x=xx, y=phospep.tab[1,], label=phospep.tab[3,], pos=3, cex=1.5, font=2, col="red")
  legend("top",
       horiz = TRUE,
       legend = c("Regulated", "Unregulated"),
       col = c("grey75", "grey30"),
       #col = c("red", "grey60"),
       pch = 15,
       bty = "n",
       pt.cex = 3,
       cex=1.5
  )
dev.off()

#---------------------

col00 <- vector(mode="numeric", length=timepnt.n)
col0 <- vector(mode="numeric", length=timepnt.n)
col1 <- vector(mode="numeric", length=timepnt.n)
col2 <- vector(mode="numeric", length=timepnt.n)
for(i in 1:timepnt.n) col00[i] <- length(dat.ratio[,i]) #number of all phos.peps
for(i in 1:timepnt.n) col0[i] <- length(na.omit(dat.ratio[,i])) #number of phos.peps after time point-specific NA filter
for(i in 1:timepnt.n) col1[i] <- length(dat.ratio.used[,i]) #number of phos.peps after across time NA filter
for(i in 1:timepnt.n) col2[i] <- reg.peps.adj.n.vec[i] #number of reg. phos.peps (adj. p-value)


phospep.tab <- rbind (col1,col1-col2,col2)
colnames (phospep.tab) <-  c("1 min", "2.5 min", "5 min", "15 min", "30 min", "60 min") 
rownames (phospep.tab) <- c('#Phosphopeptides', "#Unreg.phosphopeptides", "#Sig.reg.phosphopeptides")

png(filename=paste0(cwd, "/number_regulated_phosphopeptides_adj.png "), height=2000, width=2000, res=300)
  #mar = c(bottom, left, top, right)
  par(mar = c(5.1, 5.1, 4.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
  xx <- barplot(phospep.tab[2:3,],
              border="black",
              col = c("grey30", "grey75"),
              #space=0.04,
              font.axis=2,
              font.lab=2,
              main ="Regulated phosphopeptides\nin time course of insulin treatment",
              ylim=c(0,14000),
              cex.axis=1.1,
              cex.lab=1.5,
              cex.main=1.5,
              cex.names=1.2,
              las=1
  )
  title(ylab="#Phosphopeptides", mgp=c(4,1,0),font.lab=2,cex.lab=1.5)
  text(x=xx, y=phospep.tab[1,], label=phospep.tab[3,], pos=1, cex=0.75, font=2, col="white")
  text(x=xx, y=phospep.tab[1,], label=phospep.tab[1,], pos=3, cex=0.75, font=2, col="black")
  text(x=xx, y=phospep.tab[2,], label=phospep.tab[2,], pos=1, cex=0.75, font=2, col="white")
  legend("top",
       horiz = TRUE,
       legend = c("Regulated", "Unregulated"),
       col = c("grey75", "grey30"),
       pch = 15,
       bty = "n",
       pt.cex = 3,
       cex=1.5
  )
dev.off()

png(filename=paste0(cwd, "/number_regulated_phosphopeptides_adj_noTitle.png "), height=2000, width=2000, res=300)
  #mar = c(bottom, left, top, right)
  par(mar = c(5.1, 5.1, 4.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
  xx <- barplot(phospep.tab[2:3,],
              border="black",
              col = c("grey30", "grey75"),
              #col = c("grey60", "red"),
              font.axis=2,
              font.lab=2,
              ylim=c(0,14499),
              cex.axis=1.1,
              cex.lab=1.5,
              cex.names=1.2,
              las=1
  )
  title(ylab="#Phosphopeptides", mgp=c(4,1,0),font.lab=2,cex.lab=1.5)
  text(x=xx, y=phospep.tab[1,], label=phospep.tab[3,], pos=3, cex=1.5, font=2, col="black")
  #text(x=xx, y=phospep.tab[1,], label=phospep.tab[3,], pos=3, cex=1.5, font=2, col="red")
  legend("top",
       horiz = TRUE,
       legend = c("Regulated", "Unregulated"),
       col = c("grey75", "grey30"),
       #col = c("red", "grey60"),
       pch = 15,
       bty = "n",
       pt.cex = 3,
       cex=1.5
  )
dev.off()

#---------------------
#---------------------

col00 <- vector(mode="numeric", length=timepnt.n)
col0 <- vector(mode="numeric", length=timepnt.n)
col1 <- vector(mode="numeric", length=timepnt.n)
col2 <- vector(mode="numeric", length=timepnt.n)
for(i in 1:timepnt.n) col00[i] <- length(dat.ratio[,i]) #number of all phos.peps
for(i in 1:timepnt.n) col0[i] <- length(na.omit(dat.ratio[,i])) #number of phos.peps after time point-specific NA filter
for(i in 1:timepnt.n) col1[i] <- length(dat.ratio.used[,i]) #number of phos.peps after across time NA filter
for(i in 1:timepnt.n) col2[i] <- reg.peps.n.vec.na[i] #number of reg. phos.peps


phospep.tab <- rbind (col0,col0-col2,col2)
colnames (phospep.tab) <-  c("1 min", "2.5 min", "5 min", "15 min", "30 min", "60 min") 
rownames (phospep.tab) <- c('#Phosphopeptides', "#Unreg.phosphopeptides", "#Sig.reg.phosphopeptides")

png(filename=paste0(cwd, "/number_regulated_phosphopeptides_na.png "), height=2000, width=2000, res=300)
  #mar = c(bottom, left, top, right)
  par(mar = c(5.1, 5.1, 4.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
  xx <- barplot(phospep.tab[2:3,],
              border="black",
              col = c("grey30", "grey75"),
              #space=0.04,
              font.axis=2,
              font.lab=2,
              main ="Regulated phosphopeptides\nin time course of insulin treatment",
              ylim=c(0,14999),
              cex.axis=1.1,
              cex.lab=1.5,
              cex.main=1.5,
              cex.names=1.2,
              las=1
  )
  title(ylab="#Phosphopeptides", mgp=c(4,1,0),font.lab=2,cex.lab=1.5)
  text(x=xx, y=phospep.tab[1,], label=phospep.tab[3,], pos=1, cex=1.5, font=2, col="white")
  text(x=xx, y=phospep.tab[1,], label=phospep.tab[1,], pos=3, cex=0.75, font=2, col="black")
  text(x=xx, y=phospep.tab[2,], label=phospep.tab[2,], pos=1, cex=0.75, font=2, col="white")
  legend("top",
       horiz = TRUE,
       legend = c("Regulated", "Unregulated"),
       col = c("grey75", "grey30"),
       pch = 15,
       bty = "n",
       pt.cex = 3,
       cex=1.5
  )
dev.off()

png(filename=paste0(cwd, "/Fig_1b_number_regulated_phosphopeptides_na_noTitle.png "), height=2000, width=2000, res=300)
  #mar = c(bottom, left, top, right)
  par(mar = c(5.1, 5.1, 4.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
  xx <- barplot(phospep.tab[2:3,],
              border="black",
              col = c("grey30", "grey75"),
              #col = c("grey60", "red"),
              font.axis=2,
              font.lab=2,
              ylim=c(0,15499),
              cex.axis=1.1,
              cex.lab=1.5,
              cex.names=1.2,
              las=1
  )
  title(ylab="#Phosphopeptides", mgp=c(4,1,0),font.lab=2,cex.lab=1.5)
  text(x=xx, y=phospep.tab[1,], label=phospep.tab[3,], pos=3, cex=1.5, font=2, col="black")
  #text(x=xx, y=phospep.tab[1,], label=phospep.tab[3,], pos=3, cex=1.5, font=2, col="red")
  legend("top",
       horiz = TRUE,
       legend = c("Regulated", "Unregulated"),
       col = c("grey75", "grey30"),
       #col = c("red", "grey60"),
       pch = 15,
       bty = "n",
       pt.cex = 3,
       cex=1.5
  )
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Figure 1 C) (Time course: numbers of reg. phosphopeptides    |
#                           observed only as mono-, di-  &     |
#                           tri-phosphopeptides or as multiple |
#                           phospho-peptidoforms).             |
#---------------------------------------------------------------

exclude.missed.cleav <- FALSE

pep.phos.counts <- matrix(0, nrow=7, ncol=timepnt.n)
colnames(pep.phos.counts) <- c("1 min", "2.5 min", "5 min", "15 min", "30 min", "60 min")
rownames(pep.phos.counts) <- c("One phosphosite",
                               "Two phosphosites",
                               "Three phosphosites",
                               "Four phosphosites",
                               "Five phosphosites",
                               "Six phosphosites",
                               "Multiple peptidoforms")
pep.phos.percents <- pep.phos.counts
if(exclude.missed.cleav){
  missed.cleav.idx <- which(missed.cleav.used > 0)
  informative.missed.cleav.idx <- giveInformativeMissedCleav(missed.cleav=missed.cleav.idx, mods.master=mods.master.used)
  del.missed.cleav.idx <- setdiff(missed.cleav.idx, informative.missed.cleav.idx)
}
for(i in 1:timepnt.n){
  #na.idx <- which(is.na(p.values.used[,i]))
  #reg.idx <- ((dat.ratio.used[,i] >= 1.5) | (dat.ratio.used[,i] <= 2/3) & p.values.used[,i] <= 0.05)
  #if(length(na.idx) > 0) reg.idx[na.idx] <- TRUE
  reg.idx <- reg.peps.idx[[i]]
  seq.used.tmp <- seq.used[reg.idx]
  mods.used.tmp <- mods.used[reg.idx]
  multi.idx <- givePeptidoforms(seq.used.tmp)
  seq.used.tmp.multi <- seq.used.tmp[multi.idx]
  seq.used.tmp.solo <- seq.used.tmp[!multi.idx]
  mods.used.tmp.multi <- mods.used.tmp[multi.idx]
  mods.used.tmp.solo <- mods.used.tmp[!multi.idx]
  
  x <- cbind(names(seq.used.tmp[multi.idx]), names(mods.used.tmp.multi), seq.used.tmp[multi.idx], mods.used.tmp.multi)
  colnames(x) <- c("Names seq.used.tmp", "Names mods.used.tmp.multi", "seq.used.tmp", "mods.used.tmp.multi")
  write.table(x=x, file=paste0(cwd, "/mods.used.tmp.multi_", i, ".txt"), row.names=F, col.names=T, sep="\t")

  
  if(exclude.missed.cleav == FALSE){
    pep.phos.counts[1,i] <- length(grep("1xPhospho", mods.used.tmp.solo))
    pep.phos.counts[2,i] <- length(grep("2xPhospho", mods.used.tmp.solo))
    pep.phos.counts[3,i] <- length(grep("3xPhospho", mods.used.tmp.solo))
    pep.phos.counts[4,i] <- length(grep("4xPhospho", mods.used.tmp.solo))
    pep.phos.counts[5,i] <- length(grep("5xPhospho", mods.used.tmp.solo))
    pep.phos.counts[6,i] <- length(grep("6xPhospho", mods.used.tmp.solo))
    pep.phos.counts[7,i] <- length(mods.used.tmp.multi)
  }else if(exclude.missed.cleav == TRUE){
    pep.phos.counts[1,i] <- length(setdiff(grep("1xPhospho", mods.used.tmp.solo), missed.cleav.idx))
    pep.phos.counts[2,i] <- length(setdiff(grep("2xPhospho", mods.used.tmp.solo), missed.cleav.idx))
    pep.phos.counts[3,i] <- length(setdiff(grep("3xPhospho", mods.used.tmp.solo), missed.cleav.idx))
    pep.phos.counts[4,i] <- length(setdiff(grep("4xPhospho", mods.used.tmp.solo), missed.cleav.idx))
    pep.phos.counts[5,i] <- length(setdiff(grep("5xPhospho", mods.used.tmp.solo), missed.cleav.idx))
    pep.phos.counts[6,i] <- length(setdiff(grep("6xPhospho", mods.used.tmp.solo), missed.cleav.idx))
    pep.phos.counts[7,i] <- length(setdiff(mods.used.tmp.multi, missed.cleav.idx))
  }
  
  #pep.phos.percents[1,i] <- round(length(grep("1xPhospho", mods.used.tmp.solo)) / length(seq.used.tmp) * 100, 2)
  #pep.phos.percents[2,i] <- round(length(grep("2xPhospho", mods.used.tmp.solo)) / length(seq.used.tmp) * 100, 2)
  #pep.phos.percents[3,i] <- round(length(grep("3xPhospho", mods.used.tmp.solo)) / length(seq.used.tmp) * 100, 2)
  pep.phos.percents[1,i] <- round(pep.phos.counts[1,i] / length(seq.used.tmp) * 100, 2)
  pep.phos.percents[2,i] <- round(pep.phos.counts[2,i] / length(seq.used.tmp) * 100, 2)
  pep.phos.percents[3,i] <- round(pep.phos.counts[3,i] / length(seq.used.tmp) * 100, 2)
  pep.phos.percents[4,i] <- round(pep.phos.counts[4,i] / length(seq.used.tmp) * 100, 2)
  pep.phos.percents[5,i] <- round(pep.phos.counts[5,i] / length(seq.used.tmp) * 100, 2)
  pep.phos.percents[6,i] <- round(pep.phos.counts[6,i] / length(seq.used.tmp) * 100, 2)
  pep.phos.percents[7,i] <- round(length(mods.used.tmp.multi) / length(seq.used.tmp) * 100, 2)
  
  #print(paste0("regulated: ", length(seq.used.tmp)))
  #print(paste0("multi total: ", length(seq.used.tmp.multi)))
  #print(paste0("solo total: ", length(seq.used.tmp.solo)))
  #print(paste0("one.phos.tmp: ", one.phos.tmp))
  #print(paste0("two.phos.tmp: ", two.phos.tmp))
  #print(paste0("three.phos.tmp: ", three.phos.tmp))
  #print(paste0("multi.phos.tmp: ", multi.phos.tmp))
  #print("++++++++++++++++++++++++++")
}
print(pep.phos.counts)
print(pep.phos.percents)

write.table(x=cbind(rownames(pep.phos.counts),pep.phos.counts), file=paste0(cwd, "/peptide_phosphosite_numbers.txt"), row.names=F, col.names=T, sep="\t")

pep.phos.counts <- pep.phos.counts[apply(pep.phos.counts!=0, 1, all),]
pep.phos.percents <- pep.phos.percents[apply(pep.phos.percents!=0, 1, all),]
pep.numbers.timepoints <- apply(pep.phos.counts, 2, sum)

png(filename=paste0(cwd, "/peptide_phosphosite_numbers.png"), height=2000, width=2000, res=300)
bp <- barplot(
  pep.phos.counts,
  ylim=c(0,1299),
  col = c("grey30", "grey60", "grey75", "grey90"),
  font.axis=2,
  font.lab=2,
  cex.axis=1.1,
  cex.lab=1.5,
  cex.main=1.5,
  cex.names=1.2,
  las=1,
  main=paste0("Number of regulated mono-, di- and tri-\nphospho peptides")
)
text(x=bp, y=pep.phos.counts[1,], label=pep.phos.counts[1,], pos=1, cex=1.3, font=2, col="white")
text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,], label=pep.phos.counts[2,], pos=1, cex=1.3, font=2, col="white")
text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,]+pep.phos.counts[3,], label=pep.phos.counts[3,], pos=1, offset=0.1, cex=0.4, font=2, col="black")
text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,]+pep.phos.counts[3,]+pep.phos.counts[4,], label=pep.phos.counts[4,], pos=1, offset=0.15, cex=0.8, font=2, col="black")
text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,]+pep.phos.counts[3,]+pep.phos.counts[4,], label=pep.numbers.timepoints, pos=3, cex=1.0, font=2, col="black")

legend("top",
       legend = c("Mono", "Di", "Tri", "Other"),
       col = c("grey30", "grey60", "grey75", "grey90"),
       pch = 15,
       x.intersp = 1.0,
       bty = "n",
       text.width = 1.0,
       text.font = 2,
       pt.cex = 3,
       cex = 1.25,
       #horiz = TRUE,
       ncol = 4
)
dev.off()

png(filename=paste0(cwd, "/peptide_phosphosite_percentages.png"), height=2000, width=2000, res=300)
bp <- barplot(
  pep.phos.percents,
  ylim=c(0,119),
  col = c("grey30", "grey60", "grey75", "grey90"),
  font.axis=2,
  font.lab=2,
  cex.axis=1.1,
  cex.lab=1.5,
  cex.main=1.5,
  cex.names=1.2,
  las=1,
  main=paste0("Shares of regulated mono-, di- and tri-\nphospho peptides")
)
text(x=bp, y=pep.phos.percents[1,], label=paste0(pep.phos.percents[1,],"%"), pos=1, cex=1.0, font=2, col="white")
text(x=bp, y=pep.phos.percents[1,]+pep.phos.percents[2,], label=paste0(pep.phos.percents[2,],"%"), pos=1, cex=1.0, font=2, col="white")
text(x=bp, y=pep.phos.percents[1,]+pep.phos.percents[2,]+pep.phos.percents[3,], label=paste0(pep.phos.percents[3,],"%"), pos=1, offset=0.1, cex=0.55, font=2, col="black")
text(x=bp, y=pep.phos.percents[1,]+pep.phos.percents[2,]+pep.phos.percents[3,]+pep.phos.percents[4,], label=paste0(pep.phos.percents[4,],"%"), pos=1, offset=0.25, cex=1.0, font=2, col="black")

legend("top",
       legend = c("Mono", "Di", "Tri", "Multi"),
       col = c("grey30", "grey60", "grey75", "grey90"),
       pch = 15,
       x.intersp = 1.0,
       bty = "n",
       text.width = 1.0,
       text.font = 2,
       pt.cex = 3,
       cex = 1.25,
       #horiz = TRUE,
       ncol = 4
)
dev.off()

png(filename=paste0(cwd, "/peptide_phosphosite_numbers_noTitle.png"), height=2000, width=2000, res=300)
  #mar = c(bottom, left, top, right)
  par(mar = c(5.1, 5.1, 3.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
  bp <- barplot(
    pep.phos.counts,
    ylim=c(0,1399),
    col = c("grey30", "grey60", "grey75", "grey90"),
    font.axis=2,
    font.lab=2,
    cex.axis=1.1,
    cex.lab=1.5,
    cex.main=1.5,
    cex.names=1.2,
    las=1
  )
  text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,]+pep.phos.counts[3,]+pep.phos.counts[4,], label=pep.numbers.timepoints, pos=3, cex=1.3, font=2, col="black")

  legend("top",
       legend = c("Mono", "Di", "Tri", "Other"),
       col = c("grey30", "grey60", "grey75", "grey90"),
       pch = 15,
       x.intersp = 1.0,
       bty = "n",
       text.width = 1.0,
       text.font = 2,
       pt.cex = 3,
       cex = 1.25,
       #horiz = TRUE,
       ncol = 4
  )
dev.off()

#-----------------------

pep.phos.counts <- matrix(0, nrow=7, ncol=timepnt.n)
colnames(pep.phos.counts) <- c("1 min", "2.5 min", "5 min", "15 min", "30 min", "60 min")
rownames(pep.phos.counts) <- c("One phosphosite",
                               "Two phosphosites",
                               "Three phosphosites",
                               "Four phosphosites",
                               "Five phosphosites",
                               "Six phosphosites",
                               "Multiple peptidoforms")
pep.phos.percents <- pep.phos.counts
if(exclude.missed.cleav){
  missed.cleav.idx <- which(missed.cleav.used > 0)
  informative.missed.cleav.idx <- giveInformativeMissedCleav(missed.cleav=missed.cleav.idx, mods.master=mods.master.used)
  del.missed.cleav.idx <- setdiff(missed.cleav.idx, informative.missed.cleav.idx)
}
for(i in 1:timepnt.n){
  #na.idx <- which(is.na(p.values.used[,i]))
  #reg.idx <- ((dat.ratio.used[,i] >= 1.5) | (dat.ratio.used[,i] <= 2/3) & p.values.used[,i] <= 0.05)
  #if(length(na.idx) > 0) reg.idx[na.idx] <- TRUE
  reg.idx <- reg.peps.adj.idx[[i]]
  seq.used.tmp <- seq.used[reg.idx]
  mods.used.tmp <- mods.used[reg.idx]
  multi.idx <- givePeptidoforms(seq.used.tmp)
  seq.used.tmp.multi <- seq.used.tmp[multi.idx]
  seq.used.tmp.solo <- seq.used.tmp[!multi.idx]
  mods.used.tmp.multi <- mods.used.tmp[multi.idx]
  mods.used.tmp.solo <- mods.used.tmp[!multi.idx]
  
  x <- cbind(names(seq.used.tmp[multi.idx]), names(mods.used.tmp.multi), seq.used.tmp[multi.idx], mods.used.tmp.multi)
  colnames(x) <- c("Names seq.used.tmp", "Names mods.used.tmp.multi", "seq.used.tmp", "mods.used.tmp.multi")
  write.table(x=x, file=paste0(cwd, "/mods.used.tmp.multi_adj_", i, ".txt"), row.names=F, col.names=T, sep="\t")
  
  
  if(exclude.missed.cleav == FALSE){
    pep.phos.counts[1,i] <- length(grep("1xPhospho", mods.used.tmp.solo))
    pep.phos.counts[2,i] <- length(grep("2xPhospho", mods.used.tmp.solo))
    pep.phos.counts[3,i] <- length(grep("3xPhospho", mods.used.tmp.solo))
    pep.phos.counts[4,i] <- length(grep("4xPhospho", mods.used.tmp.solo))
    pep.phos.counts[5,i] <- length(grep("5xPhospho", mods.used.tmp.solo))
    pep.phos.counts[6,i] <- length(grep("6xPhospho", mods.used.tmp.solo))
    pep.phos.counts[7,i] <- length(mods.used.tmp.multi)
  }else if(exclude.missed.cleav == TRUE){
    pep.phos.counts[1,i] <- length(setdiff(grep("1xPhospho", mods.used.tmp.solo), missed.cleav.idx))
    pep.phos.counts[2,i] <- length(setdiff(grep("2xPhospho", mods.used.tmp.solo), missed.cleav.idx))
    pep.phos.counts[3,i] <- length(setdiff(grep("3xPhospho", mods.used.tmp.solo), missed.cleav.idx))
    pep.phos.counts[4,i] <- length(setdiff(grep("4xPhospho", mods.used.tmp.solo), missed.cleav.idx))
    pep.phos.counts[5,i] <- length(setdiff(grep("5xPhospho", mods.used.tmp.solo), missed.cleav.idx))
    pep.phos.counts[6,i] <- length(setdiff(grep("6xPhospho", mods.used.tmp.solo), missed.cleav.idx))
    pep.phos.counts[7,i] <- length(setdiff(mods.used.tmp.multi, missed.cleav.idx))
  }
  
  #pep.phos.percents[1,i] <- round(length(grep("1xPhospho", mods.used.tmp.solo)) / length(seq.used.tmp) * 100, 2)
  #pep.phos.percents[2,i] <- round(length(grep("2xPhospho", mods.used.tmp.solo)) / length(seq.used.tmp) * 100, 2)
  #pep.phos.percents[3,i] <- round(length(grep("3xPhospho", mods.used.tmp.solo)) / length(seq.used.tmp) * 100, 2)
  pep.phos.percents[1,i] <- round(pep.phos.counts[1,i] / length(seq.used.tmp) * 100, 2)
  pep.phos.percents[2,i] <- round(pep.phos.counts[2,i] / length(seq.used.tmp) * 100, 2)
  pep.phos.percents[3,i] <- round(pep.phos.counts[3,i] / length(seq.used.tmp) * 100, 2)
  pep.phos.percents[4,i] <- round(pep.phos.counts[4,i] / length(seq.used.tmp) * 100, 2)
  pep.phos.percents[5,i] <- round(pep.phos.counts[5,i] / length(seq.used.tmp) * 100, 2)
  pep.phos.percents[6,i] <- round(pep.phos.counts[6,i] / length(seq.used.tmp) * 100, 2)
  pep.phos.percents[7,i] <- round(length(mods.used.tmp.multi) / length(seq.used.tmp) * 100, 2)
  
  #print(paste0("regulated: ", length(seq.used.tmp)))
  #print(paste0("multi total: ", length(seq.used.tmp.multi)))
  #print(paste0("solo total: ", length(seq.used.tmp.solo)))
  #print(paste0("one.phos.tmp: ", one.phos.tmp))
  #print(paste0("two.phos.tmp: ", two.phos.tmp))
  #print(paste0("three.phos.tmp: ", three.phos.tmp))
  #print(paste0("multi.phos.tmp: ", multi.phos.tmp))
  #print("++++++++++++++++++++++++++")
}
print(pep.phos.counts)
print(pep.phos.percents)

write.table(x=cbind(rownames(pep.phos.counts),pep.phos.counts), file=paste0(cwd, "/peptide_phosphosite_numbers_adj.txt"), row.names=F, col.names=T, sep="\t")

pep.phos.counts <- pep.phos.counts[apply(pep.phos.counts!=0, 1, all),]
pep.phos.percents <- pep.phos.percents[apply(pep.phos.percents!=0, 1, all),]
pep.numbers.timepoints <- apply(pep.phos.counts, 2, sum)

png(filename=paste0(cwd, "/peptide_phosphosite_numbers_adj.png"), height=2000, width=2000, res=300)
bp <- barplot(
  pep.phos.counts,
  ylim=c(0,799),
  col = c("grey30", "grey60", "grey75", "grey90"),
  font.axis=2,
  font.lab=2,
  cex.axis=1.1,
  cex.lab=1.5,
  cex.main=1.5,
  cex.names=1.2,
  las=1,
  main=paste0("Number of regulated mono-, di- and tri-\nphospho peptides")
)
text(x=bp, y=pep.phos.counts[1,], label=pep.phos.counts[1,], pos=1, cex=1.3, font=2, col="white")
text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,], label=pep.phos.counts[2,], pos=1, cex=1.3, font=2, col="white")
text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,]+pep.phos.counts[3,], label=pep.phos.counts[3,], pos=1, offset=0.1, cex=0.4, font=2, col="black")
text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,]+pep.phos.counts[3,]+pep.phos.counts[4,], label=pep.phos.counts[4,], pos=1, offset=0.15, cex=0.8, font=2, col="black")
text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,]+pep.phos.counts[3,]+pep.phos.counts[4,], label=pep.numbers.timepoints, pos=3, cex=1.0, font=2, col="black")

legend("top",
       legend = c("Mono", "Di", "Tri", "Other"),
       col = c("grey30", "grey60", "grey75", "grey90"),
       pch = 15,
       x.intersp = 1.0,
       bty = "n",
       text.width = 1.0,
       text.font = 2,
       pt.cex = 3,
       cex = 1.25,
       #horiz = TRUE,
       ncol = 4
)
dev.off()

png(filename=paste0(cwd, "/peptide_phosphosite_percentages_adj.png"), height=2000, width=2000, res=300)
bp <- barplot(
  pep.phos.percents,
  ylim=c(0,119),
  col = c("grey30", "grey60", "grey75", "grey90"),
  font.axis=2,
  font.lab=2,
  cex.axis=1.1,
  cex.lab=1.5,
  cex.main=1.5,
  cex.names=1.2,
  las=1,
  main=paste0("Shares of regulated mono-, di- and tri-\nphospho peptides")
)
text(x=bp, y=pep.phos.percents[1,], label=paste0(pep.phos.percents[1,],"%"), pos=1, cex=1.0, font=2, col="white")
text(x=bp, y=pep.phos.percents[1,]+pep.phos.percents[2,], label=paste0(pep.phos.percents[2,],"%"), pos=1, cex=1.0, font=2, col="white")
text(x=bp, y=pep.phos.percents[1,]+pep.phos.percents[2,]+pep.phos.percents[3,], label=paste0(pep.phos.percents[3,],"%"), pos=1, offset=0.1, cex=0.55, font=2, col="black")
text(x=bp, y=pep.phos.percents[1,]+pep.phos.percents[2,]+pep.phos.percents[3,]+pep.phos.percents[4,], label=paste0(pep.phos.percents[4,],"%"), pos=1, offset=0.25, cex=1.0, font=2, col="black")

legend("top",
       legend = c("Mono", "Di", "Tri", "Multi"),
       col = c("grey30", "grey60", "grey75", "grey90"),
       pch = 15,
       x.intersp = 1.0,
       bty = "n",
       text.width = 1.0,
       text.font = 2,
       pt.cex = 3,
       cex = 1.25,
       #horiz = TRUE,
       ncol = 4
)
dev.off()

png(filename=paste0(cwd, "/peptide_phosphosite_numbers_adj_noTitle.png"), height=2000, width=2000, res=300)
#mar = c(bottom, left, top, right)
par(mar = c(5.1, 5.1, 3.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
bp <- barplot(
  pep.phos.counts,
  ylim=c(0,799),
  col = c("grey30", "grey60", "grey75", "grey90"),
  font.axis=2,
  font.lab=2,
  cex.axis=1.1,
  cex.lab=1.5,
  cex.main=1.5,
  cex.names=1.2,
  las=1
)
text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,]+pep.phos.counts[3,]+pep.phos.counts[4,], label=pep.numbers.timepoints, pos=3, cex=1.3, font=2, col="black")

legend("top",
       legend = c("Mono", "Di", "Tri", "Other"),
       col = c("grey30", "grey60", "grey75", "grey90"),
       pch = 15,
       x.intersp = 1.0,
       bty = "n",
       text.width = 1.0,
       text.font = 2,
       pt.cex = 3,
       cex = 1.25,
       #horiz = TRUE,
       ncol = 4
)
dev.off()

#-----------------------
#-----------------------

pep.phos.counts <- matrix(0, nrow=7, ncol=timepnt.n)
colnames(pep.phos.counts) <- c("1 min", "2.5 min", "5 min", "15 min", "30 min", "60 min")
rownames(pep.phos.counts) <- c("One phosphosite",
                               "Two phosphosites",
                               "Three phosphosites",
                               "Four phosphosites",
                               "Five phosphosites",
                               "Six phosphosites",
                               "Multiple peptidoforms")
pep.phos.percents <- pep.phos.counts
if(exclude.missed.cleav){
  missed.cleav.idx <- which(missed.cleav > 0)
  informative.missed.cleav.idx <- giveInformativeMissedCleav(missed.cleav=missed.cleav.idx, mods.master=mods.master)
  del.missed.cleav.idx <- setdiff(missed.cleav.idx, informative.missed.cleav.idx)
}
for(i in 1:timepnt.n){
  reg.idx <- which(rownames(dat.ratio) %in% reg.peps.na[[i]])
  seq.tmp <- seq[reg.idx]
  mods.tmp <- mods[reg.idx]
  multi.idx <- givePeptidoforms(seq.tmp)
  seq.tmp.multi <- seq.tmp[multi.idx]
  seq.tmp.solo <- seq.tmp[!multi.idx]
  mods.tmp.multi <- mods.tmp[multi.idx]
  mods.tmp.solo <- mods.tmp[!multi.idx]
  
  x <- cbind(names(seq.tmp[multi.idx]), names(mods.tmp.multi), seq.tmp[multi.idx], mods.tmp.multi)
  colnames(x) <- c("Names seq.tmp", "Names mods.tmp.multi", "seq.tmp", "mods.tmp.multi")
  write.table(x=x, file=paste0(cwd, "/mods.tmp.multi_", i, ".txt"), row.names=F, col.names=T, sep="\t")
  
  
  if(exclude.missed.cleav == FALSE){
    pep.phos.counts[1,i] <- length(grep("1xPhospho", mods.tmp.solo))
    pep.phos.counts[2,i] <- length(grep("2xPhospho", mods.tmp.solo))
    pep.phos.counts[3,i] <- length(grep("3xPhospho", mods.tmp.solo))
    pep.phos.counts[4,i] <- length(grep("4xPhospho", mods.tmp.solo))
    pep.phos.counts[5,i] <- length(grep("5xPhospho", mods.tmp.solo))
    pep.phos.counts[6,i] <- length(grep("6xPhospho", mods.tmp.solo))
    pep.phos.counts[7,i] <- length(mods.tmp.multi)
  }else if(exclude.missed.cleav == TRUE){
    pep.phos.counts[1,i] <- length(setdiff(grep("1xPhospho", mods.tmp.solo), missed.cleav.idx))
    pep.phos.counts[2,i] <- length(setdiff(grep("2xPhospho", mods.tmp.solo), missed.cleav.idx))
    pep.phos.counts[3,i] <- length(setdiff(grep("3xPhospho", mods.tmp.solo), missed.cleav.idx))
    pep.phos.counts[4,i] <- length(setdiff(grep("4xPhospho", mods.tmp.solo), missed.cleav.idx))
    pep.phos.counts[5,i] <- length(setdiff(grep("5xPhospho", mods.tmp.solo), missed.cleav.idx))
    pep.phos.counts[6,i] <- length(setdiff(grep("6xPhospho", mods.tmp.solo), missed.cleav.idx))
    pep.phos.counts[7,i] <- length(setdiff(mods.tmp.multi, missed.cleav.idx))
  }
  
  #pep.phos.percents[1,i] <- round(length(grep("1xPhospho", mods.tmp.solo)) / length(seq.tmp) * 100, 2)
  #pep.phos.percents[2,i] <- round(length(grep("2xPhospho", mods.tmp.solo)) / length(seq.tmp) * 100, 2)
  #pep.phos.percents[3,i] <- round(length(grep("3xPhospho", mods.tmp.solo)) / length(seq.tmp) * 100, 2)
  pep.phos.percents[1,i] <- round(pep.phos.counts[1,i] / length(seq.tmp) * 100, 2)
  pep.phos.percents[2,i] <- round(pep.phos.counts[2,i] / length(seq.tmp) * 100, 2)
  pep.phos.percents[3,i] <- round(pep.phos.counts[3,i] / length(seq.tmp) * 100, 2)
  pep.phos.percents[4,i] <- round(pep.phos.counts[4,i] / length(seq.tmp) * 100, 2)
  pep.phos.percents[5,i] <- round(pep.phos.counts[5,i] / length(seq.tmp) * 100, 2)
  pep.phos.percents[6,i] <- round(pep.phos.counts[6,i] / length(seq.tmp) * 100, 2)
  pep.phos.percents[7,i] <- round(length(mods.tmp.multi) / length(seq.tmp) * 100, 2)
  
  #print(paste0("regulated: ", length(seq.tmp)))
  #print(paste0("multi total: ", length(seq.tmp.multi)))
  #print(paste0("solo total: ", length(seq.tmp.solo)))
  #print(paste0("one.phos.tmp: ", one.phos.tmp))
  #print(paste0("two.phos.tmp: ", two.phos.tmp))
  #print(paste0("three.phos.tmp: ", three.phos.tmp))
  #print(paste0("multi.phos.tmp: ", multi.phos.tmp))
  #print("++++++++++++++++++++++++++")
}
print(pep.phos.counts)
print(pep.phos.percents)

write.table(x=cbind(rownames(pep.phos.counts),pep.phos.counts), file=paste0(cwd, "/peptide_phosphosite_numbers_na.txt"), row.names=F, col.names=T, sep="\t")

pep.phos.counts <- pep.phos.counts[apply(pep.phos.counts!=0, 1, all),]
pep.phos.percents <- pep.phos.percents[apply(pep.phos.percents!=0, 1, all),]
pep.numbers.timepoints <- apply(pep.phos.counts, 2, sum)

png(filename=paste0(cwd, "/peptide_phosphosite_numbers_na.png"), height=2000, width=2000, res=300)
  bp <- barplot(
    pep.phos.counts,
    ylim=c(0,2099),
    col = c("grey30", "grey60", "grey75", "grey90"),
    font.axis=2,
    font.lab=2,
    cex.axis=1.1,
    cex.lab=1.5,
    cex.main=1.5,
    cex.names=1.2,
    las=1,
    main=paste0("Number of regulated mono-, di- and tri-\nphospho peptides")
  )
  text(x=bp, y=pep.phos.counts[1,], label=pep.phos.counts[1,], pos=1, cex=1.3, font=2, col="white")
  text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,], label=pep.phos.counts[2,], pos=1, cex=1.3, font=2, col="white")
  text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,]+pep.phos.counts[3,], label=pep.phos.counts[3,], pos=1, offset=0.1, cex=0.4, font=2, col="black")
  text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,]+pep.phos.counts[3,]+pep.phos.counts[4,], label=pep.phos.counts[4,], pos=1, offset=0.15, cex=0.8, font=2, col="black")
  text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,]+pep.phos.counts[3,]+pep.phos.counts[4,], label=pep.numbers.timepoints, pos=3, cex=1.0, font=2, col="black")

  legend("top",
       legend = c("Mono", "Di", "Tri", "Other"),
       col = c("grey30", "grey60", "grey75", "grey90"),
       pch = 15,
       x.intersp = 1.0,
       bty = "n",
       text.width = 1.0,
       text.font = 2,
       pt.cex = 3,
       cex = 1.25,
       #horiz = TRUE,
       ncol = 4
  )
dev.off()

png(filename=paste0(cwd, "/peptide_phosphosite_percentages_na.png"), height=2000, width=2000, res=300)
  bp <- barplot(
    pep.phos.percents,
    ylim=c(0,119),
    col = c("grey30", "grey60", "grey75", "grey90"),
    font.axis=2,
    font.lab=2,
    cex.axis=1.1,
    cex.lab=1.5,
    cex.main=1.5,
    cex.names=1.2,
    las=1,
    main=paste0("Shares of regulated mono-, di- and tri-\nphospho peptides")
  )
  text(x=bp, y=pep.phos.percents[1,], label=paste0(pep.phos.percents[1,],"%"), pos=1, cex=1.0, font=2, col="white")
  text(x=bp, y=pep.phos.percents[1,]+pep.phos.percents[2,], label=paste0(pep.phos.percents[2,],"%"), pos=1, cex=1.0, font=2, col="white")
  text(x=bp, y=pep.phos.percents[1,]+pep.phos.percents[2,]+pep.phos.percents[3,], label=paste0(pep.phos.percents[3,],"%"), pos=1, offset=0.1, cex=0.55, font=2, col="black")
  text(x=bp, y=pep.phos.percents[1,]+pep.phos.percents[2,]+pep.phos.percents[3,]+pep.phos.percents[4,], label=paste0(pep.phos.percents[4,],"%"), pos=1, offset=0.25, cex=1.0, font=2, col="black")

  legend("top",
       legend = c("Mono", "Di", "Tri", "Multi"),
       col = c("grey30", "grey60", "grey75", "grey90"),
       pch = 15,
       x.intersp = 1.0,
       bty = "n",
       text.width = 1.0,
       text.font = 2,
       pt.cex = 3,
       cex = 1.25,
       #horiz = TRUE,
       ncol = 4
  )
dev.off()

png(filename=paste0(cwd, "/Fig_1c_peptide_phosphosite_numbers_noTitle_na.png"), height=2000, width=2000, res=300)
  #mar = c(bottom, left, top, right)
  par(mar = c(5.1, 5.1, 3.1, 1.1)) #default: mar = c(5.1, 4.1, 4.1, 2.1)
  bp <- barplot(
    pep.phos.counts,
    ylim=c(0,2099),
    col = c("grey30", "grey60", "grey75", "grey90"),
    font.axis=2,
    font.lab=2,
    cex.axis=1.1,
    cex.lab=1.5,
    cex.main=1.5,
    cex.names=1.2,
    las=1
  )
  text(x=bp, y=pep.phos.counts[1,]+pep.phos.counts[2,]+pep.phos.counts[3,]+pep.phos.counts[4,], label=pep.numbers.timepoints, pos=3, cex=1.3, font=2, col="black")

  legend("top",
       legend = c("Mono", "Di", "Tri", "Other"),
       col = c("grey30", "grey60", "grey75", "grey90"),
       pch = 15,
       x.intersp = 1.0,
       bty = "n",
       text.width = 1.0,
       text.font = 2,
       pt.cex = 3,
       cex = 1.25,
       #horiz = TRUE,
       ncol = 4
  )
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Supl Figure X (Time course: numbers of master proteins of    |
#                            protein groups of reg. phospho-   |
#                            peptides with one, two, three     |
#                            and multiple phosphosites         |
#                            reliably localized by PD).        |
#---------------------------------------------------------------

# col.seq <- grep("Sequence", colnames(dat.orig), value=TRUE)
# col.mods <- grep("^Modifications$", colnames(dat.orig), value=TRUE)
# col.pos <- grep("^Positions in Master", colnames(dat.orig), value=TRUE)
# col.mods.master <- grep("^Modifications in Master", colnames(dat.orig), value=TRUE)
# prot.phos.counts <- matrix(0, nrow=4, ncol=timepnt.n)
# colnames(prot.phos.counts) <- c("1 min", "2.5 min", "5 min", "15 min", "30 min", "60 min")
# rownames(prot.phos.counts) <- c("One phosphosite", "Two phosphosites", "Three phosphosites", "Multiple phosphosites")
# prot.phos.stat <- list()
# protein.list <- list()
# for(i in 1:timepnt.n){
#     #na.idx <- which(is.na(p.values.used[,i]))
#     #reg.idx <- ((dat.ratio.used[,i] >= 1.5) | (dat.ratio.used[,i] <= 2/3) & p.values.used[,i] <= 0.05)
#     #if(length(na.idx) > 0) reg.idx[na.idx] <- TRUE
#     reg.idx <- reg.peps.idx[[i]]
#     reg.rownames <- rownames(dat.ratio.used)[reg.idx]
#     phos.prot.tab <- dat.orig[reg.rownames, c(col.seq, col.mods, col.pos, col.mods.master)]
#     prot.phos.stat[[i]] <- giveProtPhosStatus(phos.prot.tab)
#     print(prot.phos.stat[[i]][[1]])
#     protein.list[[i]] <- names(prot.phos.stat[[i]])
#     for(j in 1:length(prot.phos.stat[[i]])){
#     #  if(prot.phos.stat[[i]][[j]] == 1) prot.phos.counts[1,i] <- prot.phos.counts[1,i] +1
#     #  else if(prot.phos.stat[[i]][[j]] == 2) prot.phos.counts[2,i] <- prot.phos.counts[2,i] +1
#     #  else if(prot.phos.stat[[i]][[j]] == 3) prot.phos.counts[3,i] <- prot.phos.counts[3,i] +1
#     #  else if(prot.phos.stat[[i]][[j]] > 3) prot.phos.counts[4,i] <- prot.phos.counts[4,i] +1
#       
#       if(length(prot.phos.stat[[i]][[j]]) == 1) prot.phos.counts[1,i] <- prot.phos.counts[1,i] +1
#       else if(length(prot.phos.stat[[i]][[j]]) == 2) prot.phos.counts[2,i] <- prot.phos.counts[2,i] +1
#       else if(length(prot.phos.stat[[i]][[j]]) == 3) prot.phos.counts[3,i] <- prot.phos.counts[3,i] +1
#       else if(length(prot.phos.stat[[i]][[j]]) > 3) prot.phos.counts[4,i] <- prot.phos.counts[4,i] +1
#     }
# }
# 
# assoc.protein.numbers.timepoints <- apply(prot.phos.counts, 2, sum)
# assoc.protein.number.total <- length(unique(c(protein.list[[1]], protein.list[[2]], protein.list[[3]], protein.list[[4]], protein.list[[5]], protein.list[[6]])))
# 
# prot.phos.percents <- prot.phos.counts
# for(i in 1:timepnt.n){
#   prot.phos.percents[,i] <- round(prot.phos.percents[,i] / sum(prot.phos.percents[,i]) * 100, 2)
# }
# 
# png(filename=paste0(cwd, "/protein_phosphosite_numbers.png "), height=2000, width=2000, res=300)
#       bp <- barplot(
#           prot.phos.counts,
#           ylim=c(0,1099),
#           col = c("grey30", "grey60", "grey75", "grey90"),
#           font.axis=2,
#           font.lab=2,
#           cex.axis=1.1,
#           cex.lab=1.5,
#           cex.main=1.5,
#           cex.names=1.2,
#           las=1,
#           main=paste0("Phosphosite numbers of ", assoc.protein.number.total, " proteins\nlinked to regulated peptides")
#       )
#       text(x=bp, y=prot.phos.counts[1,], label=prot.phos.counts[1,], pos=1, cex=1.3, font=2, col="white")
#       text(x=bp, y=prot.phos.counts[1,]+prot.phos.counts[2,], label=prot.phos.counts[2,], pos=1, cex=1.3, font=2, col="white")
#       text(x=bp, y=prot.phos.counts[1,]+prot.phos.counts[2,]+prot.phos.counts[3,], label=prot.phos.counts[3,], pos=1, cex=1.3, font=2, col="black")
#       text(x=bp, y=prot.phos.counts[1,]+prot.phos.counts[2,]+prot.phos.counts[3,]+prot.phos.counts[4,], label=prot.phos.counts[4,], pos=1, cex=1.3, font=2, col="black")
#       text(x=bp, y=prot.phos.counts[1,]+prot.phos.counts[2,]+prot.phos.counts[3,]+prot.phos.counts[4,], label=assoc.protein.numbers.timepoints, pos=3, cex=1.3, font=2, col="black")
# 
#       legend("top",
#              legend = c("1 site", "2 sites", "3 sites", "> 3 sites"),
#              col = c("grey30", "grey60", "grey75", "grey90"),
#              pch = 15,
#              x.intersp = 1.0,
#              bty = "n",
#              text.width = 1.15,
#              text.font = 2,
#              pt.cex = 3,
#              cex = 1.25,
#              #horiz = TRUE,
#              ncol = 4
#       )
# dev.off()
# 
# 
# 
# png(filename=paste0(cwd, "/protein_phosphosite_percentages.png "), height=2000, width=2000, res=300)
#   bp <- barplot(
#         prot.phos.percents,
#         ylim=c(0,119),
#         col = c("grey30", "grey60", "grey75", "grey90"),
#         font.axis=2,
#         font.lab=2,
#         cex.axis=1.1,
#         cex.lab=1.5,
#         cex.main=1.5,
#         cex.names=1.2,
#         las=1,
#         main=paste0("Shares of phosphosite numbers of\n", assoc.protein.number.total, " proteins linked to regulated peptides")
#       )
#       text(x=bp, y=prot.phos.percents[1,], label=paste0(prot.phos.percents[1,],"%"), pos=1, cex=1.0, font=2, col="white")
#       text(x=bp, y=prot.phos.percents[1,]+prot.phos.percents[2,], label=paste0(prot.phos.percents[2,],"%"), pos=1, cex=1.0, font=2, col="white")
#       text(x=bp, y=prot.phos.percents[1,]+prot.phos.percents[2,]+prot.phos.percents[3,], label=paste0(prot.phos.percents[3,],"%"), pos=1, cex=1.0, font=2, col="black")
#       text(x=bp, y=prot.phos.percents[1,]+prot.phos.percents[2,]+prot.phos.percents[3,]+prot.phos.percents[4,], label=paste0(prot.phos.percents[4,],"%"), pos=1, cex=1.0, font=2, col="black")
# 
#       legend("top",
#              legend = c("1 site", "2 sites", "3 sites", "> 3 sites"),
#              col = c("grey30", "grey60", "grey75", "grey90"),
#              pch = 15,
#              x.intersp = 1.0,
#              bty = "n",
#              text.width = 1.15,
#              text.font = 2,
#              pt.cex = 3,
#              cex = 1.25,
#              #horiz = TRUE,
#              ncol = 4
#       )      
# dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Figure 1 F) (Venn diagram of numbers of reg.      |
#              phosphopeptides within the time      |
#              periods Early, Intermediate & Late.) |
#----------------------------------------------------

# --> ACHTUNG: Folgendes nicht notwendig, weil alle 24 NA-Peptide von Zeitpunkt1
# in den anderen beiden Zeitpunkten von Early vorhanden!
#
#na.idx1 <- which(is.na(p.values.used[,1]))
#na.idx2 <- which(is.na(p.values.used[,2]))
#na.idx3 <- which(is.na(p.values.used[,3]))
#reg1 <- ((dat.ratio.used[,1] >= 1.5) | (dat.ratio.used[,1] <= 2/3) & p.values.used[,1] <= 0.05)
#reg2 <- ((dat.ratio.used[,2] >= 1.5) | (dat.ratio.used[,2] <= 2/3) & p.values.used[,2] <= 0.05)
#reg3 <- ((dat.ratio.used[,3] >= 1.5) | (dat.ratio.used[,3] <= 2/3) & p.values.used[,3] <= 0.05)
#if(length(na.idx1) > 0) reg1[na.idx1] <- TRUE
#if(length(na.idx2) > 0) reg2[na.idx2] <- TRUE
#if(length(na.idx3) > 0) reg3[na.idx3] <- TRUE
#Early2 <- unique(c(rownames(dat.ratio.used[reg1,]),
#                  rownames(dat.ratio.used[reg2,]),
#                  rownames(dat.ratio.used[reg3,])
#))

#Early <- unique(c(na.omit(rownames(dat.ratio.used[(dat.ratio.used[,1] >= 1.5) | (dat.ratio.used[,1] <= 2/3) & p.values.used[,1] <= 0.05,])),
#                  rownames(dat.ratio.used[(dat.ratio.used[,2] >= 1.5) | (dat.ratio.used[,2] <= 2/3) & p.values.used[,2] <= 0.05,]),
#))
#
#Intermediate <- unique(c(rownames(dat.ratio.used[(dat.ratio.used[,3] >= 1.5) | (dat.ratio.used[,3] <= 2/3) & p.values.used[,3] <= 0.05,]),
#                         rownames(dat.ratio.used[(dat.ratio.used[,4] >= 1.5) | (dat.ratio.used[,4] <= 2/3) & p.values.used[,4] <= 0.05,])
#))
#
#Late <- unique(c(rownames(dat.ratio.used[(dat.ratio.used[,5] >= 1.5) | (dat.ratio.used[,5] <= 2/3) & p.values.used[,5] <= 0.05,]),
#                 rownames(dat.ratio.used[(dat.ratio.used[,6] >= 1.5) | (dat.ratio.used[,6] <= 2/3) & p.values.used[,6] <= 0.05,])
#))

Early <- unique(c(reg.peps[[1]], reg.peps[[2]]))
Intermediate <- unique(c(reg.peps[[3]], reg.peps[[4]]))
Late <- unique(c(reg.peps[[5]], reg.peps[[6]]))

png (paste0(cwd, "/#Sig.Reg.PhosphoPeptides_venn_diagram_percentage.png"), res=200, height = 1200, width=1200)
    ggvenn(
      list(Early=Early, Intermediate=Intermediate, Late=Late), 
      #fill_color = c("#E69F00", "#56B4E9", "#009E73"),
      fill_color = c("darkred", "navy", "darkgreen"),
      fill_alpha = 0.5,
      stroke_size=0.5,
      set_name_size=5,
      text_size=3.5
    )
dev.off()

library("VennDiagram")
venn.diagram(
  x = list(Early=Early, Intermediate=Intermediate, Late=Late),
  category.names = c(paste0("Early (", length(Early), ")"),
                     paste0("Intermediate (", length(Intermediate), ") "),
                     paste0("Late (", length(Late), ")")),
  filename = paste0(cwd, "/Fig_1f_#Sig.Reg.PhosphoPeptides_venn_diagram.png"),
  output = TRUE,
  imagetype="png",
  height = 600, 
  width = 600, 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("darkred", "navy", "darkgreen"),
  fill = c(alpha("darkred",0.5), alpha("navy",0.5), alpha("darkgreen",0.5)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.45,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("darkred", "navy", "darkgreen"),
  #cat.col = c("#440154ff", '#21908dff', 'indianred2'),
  rotation = 1
)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Figure 1 D) (Volcanos) |
#-------------------------

up.5min <- c("Q9UJ70 2xPhospho [S74(99.3); S76(100)]", # NAGK-S74;S76
             "P31751 1xPhospho [T309(100)]; P31749 1xPhospho [T308(100)]; Q9Y243 1xPhospho [T305(100)]", #AKT2/1/3
             "P55196 1xPhospho [S1799(100)]", # AFDN-S1799
             #"P31751 2xPhospho [S474(100); S478(100)]", # AKT2
             "Q86SQ0 2xPhospho [S200(100); S]", # PHLDB2-S200;S
             "P53985 2xPhospho [T466(100); S467(100)]", # SLC16A1-T466;S467
             "Q9UQC2 1xPhospho [S264(99.4)]", # GAB2-S264
             "O43524 1xPhospho [S253(100)]") # FOXO3-S253
down.5min <- c("P29590-2 2xPhospho [S512(100); S/T]; P29590-4 2xPhospho [S512(100); S/T]; P29590 2xPhospho [S512(100); S/T]", # PML-S512
             "Q92974 2xPhospho [S151(99.3); S163(100)]", # ARHGEF2-S151;S163
             "O14974 2xPhospho [S507(100); S509(100)]", # PPP1R12A-S507;S509
             "Q9HB90 1xMet-loss [N-Term]; 1xPhospho [S15(100)]", # RRAGC-S15
             "Q01433 1xPhospho [S76(100)]", # AMPD2-S76
             "Q9H7C4 1xPhospho [S26(100)]", # SYNC-S26
             "Q13573 1xPhospho [S224(100)]") # SNW1-S224
#mark5 <- names(mods.master.used)[which(mods.master.used %in% c(up.5min,down.5min))]
mark5 <- names(mods.master)[which(mods.master %in% c(up.5min,down.5min))]
mark5 <- intersect(reg.peps.na[[3]], mark5)
#mark5 <- mark5[-grep("Q9Y2W1_peptide8",mark5)]
mark.txt.5min <- gsub("_peptide\\d+", "", mark5)
mark.txt2.5min <- returnUniProtGeneNames(uniprot.ids=mark.txt.5min,
                       table=FALSE,
                       iso.rm=TRUE,
                       source="E:/Projekte/Phospho2.1/ID-Mapping/20220702_Phospho2.1_UniProt_to_GeneName.tsv",
                       file.name=NULL,
                       output.path=NULL)
#mark.txt2.5min <- gsub("; ", "\n", mark.txt2.5min)
mark.txt2.5min <- gsub("^AKT2; AKT1; AKT3$", "AKT2", mark.txt2.5min)

cbind(mark.txt2.5min,rep(0,14))
move.txt.x <- rep(0,length(mark.txt2.5min))
move.txt.x[c(4)] <- 1.1 #"ARHGEF2" +
move.txt.x[c(7)] <- -0.5 #"FOXO3" +
move.txt.x[c(10)] <- -1.1 #"PPP1R12A" +
move.txt.x[c(14)] <- -0.6 #"PML" +
move.txt.y <- rep(0,length(mark.txt2.5min))
move.txt.y[c(4)] <- -0.4 #"ARHGEF2" +
move.txt.y[c(10)] <- -0.4 #"PPP1R12A" +
move.txt.y[c(14)] <- 0.5 #"PML" +

#plotVolcano(fc=dat.ratio.used[,3],
#            p.val=p.values.adj.used[,3],
#            reg=reg.peps.adj[[3]],
#            mark=mark5,
#            mark.txt=mark.txt2.5min,
#            move.txt.x=move.txt.x,
#            move.txt.y=move.txt.y,
#            col.idx=9,
#            main="5 min.",
#            filename.tag="_5min_used",
#            output.path=cwd)

plotVolcano(fc=dat.ratio[,3],
            p.val=p.values[,3],
            reg=reg.peps.na[[3]],
            mark=mark5,
            mark.txt=mark.txt2.5min,
            move.txt.x=move.txt.x,
            move.txt.y=move.txt.y,
            col.idx=2,
            main="5 min",
            filename.tag="_5min",
            output.path=cwd)

# upreg:
mark.txt2.5min <- gsub("^NAGK$", "NAGK-S74;S76", mark.txt2.5min)
mark.txt2.5min <- gsub("^AKT2$", "AKT2-T309", mark.txt2.5min)
mark.txt2.5min <- gsub("^AFDN$", "AFDN-S1799", mark.txt2.5min)
mark.txt2.5min <- gsub("^PHLDB2$", "PHLDB2-S200;S", mark.txt2.5min)
mark.txt2.5min <- gsub("^SLC16A1$", "SLC16A1-T466;S467", mark.txt2.5min)
mark.txt2.5min <- gsub("^GAB2$", "GAB2-S264", mark.txt2.5min)
mark.txt2.5min <- gsub("^FOXO3$", "FOXO3-S253", mark.txt2.5min)
# downreg:
mark.txt2.5min <- gsub("^PML$", "PML-S512", mark.txt2.5min)
mark.txt2.5min <- gsub("^ARHGEF2$", "ARHGEF2-S151;S163", mark.txt2.5min)
mark.txt2.5min <- gsub("^PPP1R12A$", "PPP1R12A-S507;S509", mark.txt2.5min)
mark.txt2.5min <- gsub("^RRAGC$", "RRAGC-S15", mark.txt2.5min)
mark.txt2.5min <- gsub("^AMPD2$", "AMPD2-S76", mark.txt2.5min)
mark.txt2.5min <- gsub("^SYNC$", "SYNC-S26", mark.txt2.5min)
mark.txt2.5min <- gsub("^SNW1$", "SNW1-S224", mark.txt2.5min)

plotVolcano(fc=dat.ratio[,3],
            p.val=p.values[,3],
            reg=reg.peps.na[[3]],
            mark=mark5,
            mark.txt=mark.txt2.5min,
            move.txt.x=move.txt.x,
            move.txt.y=move.txt.y,
            col.idx=2,
            main="5 min",
            filename.tag="_5min_withSites",
            output.path=cwd)

#----------------

up.60min <- c("E9PAV3 1xPhospho [S1112(100)]", # NACA-S1112
             #"E9PAV3 2xPhospho [S585(99.1); S/T]", # NACA-S585;S/T,
             "Q53EL6 1xPhospho [S76(100)]", # PDCD4-S76
             "Q9UHB6-4 1xPhospho [S225(100)]; Q9UHB6 1xPhospho [S225(100)]", # LIMA1-S225
             "Q9NRS6 1xPhospho [T148(100)]", # SNX15-T148
             "Q14152 1xPhospho [S584(100)]", # EIF3A-S584
             "Q96B36 1xPhospho [T246(100)]", # AKT1S1-T246
             "P02545-2 1xPhospho [S404(100)]; P02545 1xPhospho [S404(100)]") # LMNA-S404
down.60min <- c("Q9Y4H2 2xPhospho [S346(99.4); T350(100)]", # IRS2-T350
               "O14974 2xPhospho [S507(100); S509(100)]", # PPP1R12A-S507;S509
               "Q92609 2xPhospho [S554(100); S557(99.2)]", # TBC1D5-S554;S557
               "Q9BTU6 2xPhospho [S462(98.4); S468(100)]", # PI4K2A-S462;S468
               "Q9Y3L3 1xPhospho [S262(100)]", # SH3BP1-S262
               "Q8IV50 1xPhospho [S24(100)]", # LYSMD2-S24
               "P35269 1xPhospho [S433(100)]", # GTF2F1-S433
               "Q15424 1xPhospho [S234(98.3)]; Q14151 1xPhospho [S233(98.3)]") # SAFB-S234
#mark60 <- names(mods.master.used)[which(mods.master.used %in% c(up.60min,down.60min))]
mark60 <- names(mods.master)[which(mods.master %in% c(up.60min,down.60min))]
#mark60 <- intersect(reg.peps.adj[[6]], mark60)
mark60 <- intersect(reg.peps.na[[6]], mark60)
mark60 <- mark60[-grep("Q53EL6_peptide5",mark60)] # same phospep (Q53EL6_peptide3/5) 
                                                # occurs 2x due to missed cleav.
mark60 <- mark60[-grep("Q96B36_peptide4",mark60)] # same phospep (Q96B36_peptide3/4) 
# occurs 2x due to missed cleav.
mark.txt.60min <- gsub("_peptide\\d+", "", mark60)
mark.txt2.60min <- returnUniProtGeneNames(uniprot.ids=mark.txt.60min,
                                          table=FALSE,
                                          iso.rm=TRUE,
                                          source="E:/Projekte/Phospho2.1/ID-Mapping/20220702_Phospho2.1_UniProt_to_GeneName.tsv",
                                          file.name=NULL,
                                          output.path=NULL)
## Why for-loop???
#mark.txt2.60min <- vector(mode="character", length=length(mark.txt.60min))
#for(i in 1:length(mark.txt.60min)){
#    mark.txt2.60min[i] <- returnUniProtGeneNames(uniprot.ids=mark.txt.60min[i],
#                                    table=FALSE,
#                                    iso.rm=TRUE,
#                                    source="E:/Projekte/Phospho2.1/ID-Mapping/20220702_Phospho2.1_UniProt_to_GeneName.tsv",
#                                    file.name=NULL,
#                                    output.path=NULL)
#}
#mark.txt2.60min <- gsub("; ", "\n", mark.txt2.60min)
mark.txt2.60min <- gsub("SAFB; SAFB2", "SAFB", mark.txt2.60min)

cbind(mark.txt2.60min,rep(0,15))
move.txt.x <- rep(0,length(mark.txt2.60min))
move.txt.x[c(3)] <- 0.8 #NACA
move.txt.x[c(8)] <- -1 #PDCD4 +
move.txt.x[c(11)] <- -1 #PPP1R12A +
move.txt.x[c(15)] <- 1 #PI4K2A +
move.txt.y <- rep(0,length(mark.txt2.60min))
move.txt.y[c(1)] <- -1 #NACA +
move.txt.y[c(3)] <- -0.4 #NACA
move.txt.y[c(9)] <- -1 #TBC1D5 +
move.txt.y[c(11)] <- -0.4 #PPP1R12A +

#plotVolcano(fc=dat.ratio.used[,6],
#            p.val=p.values.adj.used[,6],
#            reg=reg.peps.adj[[6]],
#            mark=mark60,
#            mark.txt=mark.txt2.60min,
#            move.txt.x=move.txt.x,
#            move.txt.y=move.txt.y,
#            col.idx=11,
#            main="60 min.",
#            filename.tag="_60min",
#            output.path=cwd)

plotVolcano(fc=dat.ratio[,6],
            p.val=p.values[,6],
            reg=reg.peps.na[[6]],
            mark=mark60,
            mark.txt=mark.txt2.60min,
            move.txt.x=move.txt.x,
            move.txt.y=move.txt.y,
            col.idx=2,
            main="60 min",
            filename.tag="_60min",
            output.path=cwd)

# upreg:
mark.txt2.60min <- gsub("^NACA$", "NACA-S1112", mark.txt2.60min)
mark.txt2.60min <- gsub("^PDCD4$", "PDCD4-S76", mark.txt2.60min)
mark.txt2.60min <- gsub("^LIMA1$", "LIMA1-S225", mark.txt2.60min)
mark.txt2.60min <- gsub("^SNX15$", "SNX15-T148", mark.txt2.60min)
mark.txt2.60min <- gsub("^EIF3A$", "EIF3A-S584", mark.txt2.60min)
mark.txt2.60min <- gsub("^AKT1S1$", "AKT1S1-T246", mark.txt2.60min)
mark.txt2.60min <- gsub("^LMNA$", "LMNA-S404", mark.txt2.60min)
# downreg:
mark.txt2.60min <- gsub("^IRS2$", "IRS2-T350", mark.txt2.60min)
mark.txt2.60min <- gsub("^PPP1R12A$", "PPP1R12A-S507;S509", mark.txt2.60min)
mark.txt2.60min <- gsub("^TBC1D5$", "TBC1D5-S554;S557", mark.txt2.60min)
mark.txt2.60min <- gsub("^PI4K2A$", "PI4K2A-S462;S468", mark.txt2.60min)
mark.txt2.60min <- gsub("^SH3BP1$", "SH3BP1-S262", mark.txt2.60min)
mark.txt2.60min <- gsub("^LYSMD2$", "LYSMD2-S24", mark.txt2.60min)
mark.txt2.60min <- gsub("^GTF2F1$", "GTF2F1-S433", mark.txt2.60min)
mark.txt2.60min <- gsub("^SAFB$", "SAFB-S234", mark.txt2.60min)

plotVolcano(fc=dat.ratio[,6],
            p.val=p.values[,6],
            reg=reg.peps.na[[6]],
            mark=mark60,
            mark.txt=mark.txt2.60min,
            move.txt.x=move.txt.x,
            move.txt.y=move.txt.y,
            col.idx=2,
            main="60 min",
            filename.tag="_60min_withSites",
            output.path=cwd)

#----------------

#mark.small <- c(mark5[5], mark5[7], mark5[9], mark60[7])
mark.small <- c(mark5[5], mark5[7], mark60[7], mark5[14], mark60[15])
mark.txt.small.tmp <- c(mark.txt2.5min, mark.txt2.60min)
mark.txt.small <- c(grep("AKT2", mark.txt.small.tmp, value=T),
                    grep("FOXO3", mark.txt.small.tmp, value=T),
                    #grep("GAB2", mark.txt.small.tmp, value=T),
                    grep("PDCD4", mark.txt.small.tmp, value=T),
                    grep("PML", mark.txt.small.tmp, value=T),
                    grep("IRS2", mark.txt.small.tmp, value=T))


#----------------

#plotVolcano(fc=dat.ratio.used[,1],
#            p.val=p.values.adj.used[,1],
#            reg=reg.peps.adj[[1]],
#            mark=mark5,
#            mark.txt=NULL,
#            move.txt.x=NULL,
#            move.txt.y=NULL,
#            col.idx=9,
#            mark2=mark60,
#            col.idx2=11,
#            main="1 min.",
#            filename.tag="_1min",
#            output.path=cwd)

plotVolcano(fc=dat.ratio[,1],
            p.val=p.values[,1],
            reg=reg.peps.na[[1]],
            mark=mark.small,
            mark.txt=mark.txt.small,
            col.idx=2,
            main="1 min",
            filename.tag="_1min",
            output.path=cwd)

#----------------

#plotVolcano(fc=dat.ratio.used[,2],
#            p.val=p.values.adj.used[,2],
#            reg=reg.peps.adj[[2]],
#            mark=mark5,
#            mark.txt=NULL,
#            move.txt.x=NULL,
#            move.txt.y=NULL,
#            col.idx=9,
#            mark2=mark60,
#            col.idx2=11,
#            main="2.5 min.",
#            filename.tag="_2.5min",
#            output.path=cwd)

plotVolcano(fc=dat.ratio[,2],
            p.val=p.values[,2],
            reg=reg.peps.na[[2]],
            mark=mark.small,
            mark.txt=mark.txt.small,
            col.idx=2,
            main="2.5 min",
            filename.tag="_2.5min",
            output.path=cwd)

#----------------

#plotVolcano(fc=dat.ratio.used[,4],
#            p.val=p.values.adj.used[,4],
#            reg=reg.peps.adj[[4]],
#            mark=mark5,
#            mark.txt=NULL,
#            move.txt.x=NULL,
#            move.txt.y=NULL,
#            col.idx=9,
#            mark2=mark60,
#            col.idx2=11,
#            main="15 min.",
#            filename.tag="_15min",
#            output.path=cwd)

plotVolcano(fc=dat.ratio[,4],
            p.val=p.values[,4],
            reg=reg.peps.na[[4]],
            mark=mark.small,
            mark.txt=mark.txt.small,
            col.idx=2,
            main="15 min",
            filename.tag="_15min",
            output.path=cwd)

#----------------

#plotVolcano(fc=dat.ratio.used[,5],
#            p.val=p.values.adj.used[,5],
#            reg=reg.peps.adj[[5]],
#            mark=mark5,
#            mark.txt=NULL,
#            move.txt.x=NULL,
#            move.txt.y=NULL,
#            col.idx=9,
#            mark2=mark60,
#            col.idx2=11,
#            main="30 min.",
#            filename.tag="_30min",
#            output.path=cwd)

plotVolcano(fc=dat.ratio[,5],
            p.val=p.values[,5],
            reg=reg.peps.na[[5]],
            mark=mark.small,
            mark.txt=mark.txt.small, 
            col.idx=2,
            main="30 min",
            filename.tag="_30min",
            output.path=cwd)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Figure 1 D) (Volcano crowd plot) |
#--------------------------------------------------------

plotVolcanoCrowd(fc=dat.ratio.used, p.val=p.values.used, filename.prefix="Fig_1d_", output.path=cwd)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Figure 1 E) (PCA biplot of sig. reg. phosphopeptides) |
#--------------------------------------------------------

# --> best cluster separation?!
#fc.thresh1 <- 6
##fc.thresh2 <- 1/6
#p.thresh <- 0.01
#direction <- "positive"

#fc.thresh1 <- 15/10
#fc.thresh2 <- 10/15
#p.thresh <- 0.05
#direction <- "both"

fc.thresh1 <- 30/10
fc.thresh2 <- 10/30
p.thresh <- 0.01
direction <- "both"

sig.reg.pept.pca <- unique(c(
                      na.omit(regulatedPeptides(ratios=dat.ratio.used[,1], p.values=p.values.used[,1], r.thresh.pos=fc.thresh1, r.thresh.neg=fc.thresh2, p.thresh=p.thresh, direction=direction, output.type="row.names")),
                      regulatedPeptides(ratios=dat.ratio.used[,2], p.values=p.values.used[,2], r.thresh.pos=fc.thresh1, r.thresh.neg=fc.thresh2, p.thresh=p.thresh, direction=direction, output.type="row.names"),
                      regulatedPeptides(ratios=dat.ratio.used[,3], p.values=p.values.used[,3], r.thresh.pos=fc.thresh1, r.thresh.neg=fc.thresh2, p.thresh=p.thresh, direction=direction, output.type="row.names"),
                      regulatedPeptides(ratios=dat.ratio.used[,4], p.values=p.values.used[,4], r.thresh.pos=fc.thresh1, r.thresh.neg=fc.thresh2, p.thresh=p.thresh, direction=direction, output.type="row.names"),
                      regulatedPeptides(ratios=dat.ratio.used[,5], p.values=p.values.used[,5], r.thresh.pos=fc.thresh1, r.thresh.neg=fc.thresh2, p.thresh=p.thresh, direction=direction, output.type="row.names"),
                      regulatedPeptides(ratios=dat.ratio.used[,6], p.values=p.values.used[,6], r.thresh.pos=fc.thresh1, r.thresh.neg=fc.thresh2, p.thresh=p.thresh, direction=direction, output.type="row.names")))

plotPCA(x=dat.abundance.used[sig.reg.pept.pca,],
        group.legend.pos="bottomleft",
        idv.legend.pos="bottomright",
        filename.prefix="23-08-29_Fig_2a_",
        output.path=cwd)

plotPCA(x=dat.abundance.used,
        group.legend.pos="bottomleft",
        idv.legend.pos="bottomright",
        filename.prefix="23-08-29_Supl_Fig_3_",
        output.path=cwd)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Figure 5     (Sankey plot for B1 showing                |
#               dynamics of phosphopeptide regulation.)   |
#----------------------------------------------------------

#install.packages("remotes")
#remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(ggplot2)
library(dplyr)

col.tmp1 <- rep("Unreg_1", nrow(dat.ratio.used))
col.tmp2 <- rep("Unreg_2.5", nrow(dat.ratio.used))
col.tmp3 <- rep("Unreg_5", nrow(dat.ratio.used))
col.tmp4 <- rep("Unreg_15", nrow(dat.ratio.used))
col.tmp5 <- rep("Unreg_30", nrow(dat.ratio.used))
col.tmp6 <- rep("Unreg_60", nrow(dat.ratio.used))
reg.peps.dat <- data.frame(cbind(col.tmp1, col.tmp2, col.tmp3, col.tmp4, col.tmp5, col.tmp6))
names(reg.peps.dat) <- c("1 min", "2.5 min", "5 min", "15 min", "30 min", "60 min")

for(i in 1:timepnt.n){
  #reg.idx <- reg.peps.idx.pos[[i]]
  #reg.idx <- reg.peps.idx.neg[[i]]
  reg.idx <- reg.peps.idx[[i]]
  
  if(i > 1){
    for(j in 1:nrow(dat.ratio.used)){
      if(reg.idx[j] == TRUE){
        if(length(grep("Reg", reg.peps.dat[j,i-1]))>0 && reg.direction.matrix[j,i]==reg.direction.matrix[j,i-1]){
          #reg.memory <- as.numeric(gsub(".*\\_(\\d+)\\_.*", "\\1", reg.peps.dat[j,i-1]))
          reg.memory <- as.numeric(gsub("Reg_(\\d+)_.+", "\\1", reg.peps.dat[j,i-1]))
          #print(paste0("reg.memory: ", reg.memory))
          reg.peps.dat[j,i] <- paste0("Reg_", as.character(reg.memory+1),"_", timepoint.lab[i]) 
        }else{
          reg.peps.dat[j,i] <- paste0("Reg_1_", timepoint.lab[i])  
        }
      }  
    }
  }else{
    reg.peps.dat[reg.idx,i] <- paste0("Reg_1_", timepoint.lab[i])
  }
  
}

del.idx.sankey <- c()
for(i in 1:nrow(reg.peps.dat)){
  if(length(grep("Unreg_", reg.peps.dat[i,])) > 5) del.idx.sankey <- c(del.idx.sankey, i)
}
reg.peps.dat.out <- reg.peps.dat    # -----------------> NEU!
reg.peps.dat <- reg.peps.dat[-del.idx.sankey,]

# -----------------> NEU!
rownames(reg.peps.dat.out) <- rownames(dat.ratio.used)
uniprot.tmp <- gsub("_peptide\\d+$", "", rownames(reg.peps.dat.out))
uniprot.tmp <- gsub("-\\d+", "", uniprot.tmp)
gene.name.tmp <- returnUniProtGeneNames(uniprot.ids=uniprot.tmp,
                                        iso.rm=FALSE,
                                        source="E:/Projekte/Phospho2.1/ID-Mapping/20220702_Phospho2.1_UniProt_to_GeneName.tsv")
pval.min.tmp <- apply(p.values.used,1,min)
fc.max.tmp <- apply(abs(log2(dat.ratio.used)),1,max)
output.tmp <- cbind(rownames(reg.peps.dat.out), 
                    uniprot.tmp,
                    gene.name.tmp,
                    pval.min.tmp,
                    fc.max.tmp,
                    reg.peps.dat.out)
colnames(output.tmp)[1] <- "Phosphopeptide"
colnames(output.tmp)[2] <- "UniProt-ID"
colnames(output.tmp)[3] <- "Gene"
colnames(output.tmp)[4] <- "Min. P-Val"
colnames(output.tmp)[5] <- "Max. Log2-Fold Change"
write.table(x=output.tmp,
            #file=paste0(cwd, "/sankey.reg.peps.dat.both.txt"),
            #file=paste0(cwd, "/sankey.reg.peps.dat.upregulated.txt"),
            file=paste0(cwd, "/sankey.reg.peps.dat.downregulated.txt"),
            col.names=T,
            row.names=F,
            sep="\t")
x.tmp <- unique(pep2prot(output.tmp[grep("Reg_6_60", output.tmp[,8]), 2]))
write.table(x=x.tmp, file=paste0(cwd, "/sankey_reg_6_60.txt"), col.names=F, row.names=F)
x.tmp <- unique(pep2prot(output.tmp[grep("Reg_5_60", output.tmp[,8]), 2]))
write.table(x=x.tmp, file=paste0(cwd, "/sankey_reg_5_60.txt"), col.names=F, row.names=F)
x.tmp <- unique(pep2prot(output.tmp[grep("Reg_4_60", output.tmp[,8]), 2]))
write.table(x=x.tmp, file=paste0(cwd, "/sankey_reg_4_60.txt"), col.names=F, row.names=F)
x.tmp <- unique(pep2prot(output.tmp[grep("Reg_3_60", output.tmp[,8]), 2]))
write.table(x=x.tmp, file=paste0(cwd, "/sankey_reg_3_60.txt"), col.names=F, row.names=F)
x.tmp <- unique(pep2prot(output.tmp[grep("Reg_2_60", output.tmp[,8]), 2]))
write.table(x=x.tmp, file=paste0(cwd, "/sankey_reg_2_60.txt"), col.names=F, row.names=F)
x.tmp <- unique(pep2prot(output.tmp[grep("Reg_1_60", output.tmp[,8]), 2]))
write.table(x=x.tmp, file=paste0(cwd, "/sankey_reg_1_60.txt"), col.names=F, row.names=F)

# ontologies: "BP", "MF" & "CC"
timeCourseDotPlot(x=output.tmp,
                  time.course=6,
                  ontology="MF",
                  nodeSize=5,
                  fdr.thresh=0.05,
                  #min.anno=5,
                  #max.anno=700,
                  min.anno=5,
                  max.anno=700,
                  term.max=20,
                  cwd=cwd)

#timeCourseDotPlot2(x=output.tmp,
#                  time.course=3,
#                  ontology="CC",
#                  nodeSize=5,
#                  fdr.thresh=0.05,
#                  #min.anno=5,
#                  #max.anno=700,
#                  min.anno=5,
#                  max.anno=700,
#                  term.max=20,
#                  cwd=cwd)
# -----------------> NEU!

sankey.reg.peps.dat <- reg.peps.dat %>% make_long(names(reg.peps.dat))

dagg <- sankey.reg.peps.dat %>%
  dplyr::group_by(node) %>%
  tally()

sankey.reg.peps.dat2 <- merge(sankey.reg.peps.dat, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)
label <- sankey.reg.peps.dat2$node
sankey.reg.peps.dat2 <- cbind(sankey.reg.peps.dat2, label)

sankey.reg.peps.dat2$label <- gsub("Unreg_\\d+\\.*\\d*", "Unregulated", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_1_1$", "1 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_2_2.5$", "1 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_3_5$", "1 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_4_15$", "1 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_5_30$", "1 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_6_60$", "1 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_1_2.5$", "2.5 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_2_5$", "2.5 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_3_15$", "2.5 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_4_30$", "2.5 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_5_60$", "2.5 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_1_5$", "5 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_2_15$", "5 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_3_30$", "5 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_4_60$", "5 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_1_15$", "15 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_2_30$", "15 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_3_60$", "15 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_1_30$", "30 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_2_60$", "30 min", sankey.reg.peps.dat2$label)
sankey.reg.peps.dat2$label <- gsub("^Reg_1_60$", "60 min", sankey.reg.peps.dat2$label)



pl <- ggplot(sankey.reg.peps.dat2, aes(x = x
                                       , next_x = next_x
                                       , node = node
                                       , next_node = next_node
                                       , fill = factor(node)
                                       , label = paste0(label,"\nn=", n))
)
pl <- pl + geom_sankey(flow.alpha = 0.6
                       , node.color = "grey30"
                       ,show.legend = FALSE)
pl <- pl + geom_sankey_label(size = 3, color = "white", fill= "grey40")
pl <- pl + theme_bw()
pl <- pl + theme(legend.position = "none")
pl <- pl + theme(axis.title = element_blank()
                 , axis.text.y = element_blank()
                 , axis.ticks = element_blank()  
                 , panel.grid = element_blank())
pl <- pl + scale_fill_viridis_d()
pl <- pl + theme_sankey(base_size = 18)
pl <- pl + labs(title = "Time course of phosphopeptide regulation")
#pl <- pl + labs(subtitle = "using  David Sjoberg's ggsankey package")
#pl <- pl + labs(caption = "@techanswers88")
pl <- pl + labs(x = NULL)
pl <- pl + theme(legend.position = "none",
                 plot.title = element_text(hjust = .5))
pl2 <- pl + scale_fill_manual(values = c('Reg_1_1'    = "darkred",
                                        'Reg_2_2.5'  = "darkred",
                                        "Reg_3_5" = "darkred",
                                        "Reg_4_15" = "darkred",
                                        "Reg_5_30" = "darkred",
                                        
                                        'Reg_1_2.5'  = "darkorange3",
                                        "Reg_2_5" = "darkorange3",
                                        "Reg_3_15" = "darkorange3",
                                        "Reg_4_30" = "darkorange3",
                                        
                                        "Reg_1_5" = "darkgreen",
                                        "Reg_2_15" = "darkgreen",
                                        "Reg_3_30" = "darkgreen",
                                        
                                        "Reg_1_15" = "deeppink",
                                        "Reg_2_30" = "deeppink",
                                        
                                        "Reg_1_30" = "blue"
) )
pl2
png(filename=paste0(cwd, "/20230926_sankey_phosphosite_regulation.png "), height=2000, width=2000, res=300)
  pl2
dev.off()



pl2 <- pl + scale_fill_manual(values = c('Reg_1_1'    = "darkred",
                                         'Reg_2_2.5'  = "darkred",
                                         "Reg_3_5" = "darkred",
                                         "Reg_4_15" = "darkred",
                                         "Reg_5_30" = "darkred",
                                         
                                         'Reg_1_2.5'  = "grey40",
                                         "Reg_2_5" = "grey40",
                                         "Reg_3_15" = "grey40",
                                         "Reg_4_30" = "grey40",
                                         
                                         "Reg_1_5" = "grey40",
                                         "Reg_2_15" = "grey40",
                                         "Reg_3_30" = "grey40",
                                         
                                         "Reg_1_15" = "grey40",
                                         "Reg_2_30" = "grey40",
                                         
                                         "Reg_1_30" = "grey40"
) )
pl2
png(filename=paste0(cwd, "/20230926_sankey_phosphosite_regulation_RED_ONLY.png "), height=2000, width=2000, res=300)
  pl2
dev.off()



pl2 <- pl + scale_fill_manual(values = c('Reg_1_1'    = "grey40",
                                         'Reg_2_2.5'  = "grey40",
                                         "Reg_3_5" = "grey40",
                                         "Reg_4_15" = "grey40",
                                         "Reg_5_30" = "grey40",
                                         
                                         'Reg_1_2.5'  = "darkorange3",
                                         "Reg_2_5" = "darkorange3",
                                         "Reg_3_15" = "darkorange3",
                                         "Reg_4_30" = "darkorange3",
                                         
                                         "Reg_1_5" = "grey40",
                                         "Reg_2_15" = "grey40",
                                         "Reg_3_30" = "grey40",
                                         
                                         "Reg_1_15" = "grey40",
                                         "Reg_2_30" = "grey40",
                                         
                                         "Reg_1_30" = "grey40"
) )
pl2
png(filename=paste0(cwd, "/20230926_sankey_phosphosite_regulation_ORANGE_ONLY.png "), height=2000, width=2000, res=300)
  pl2
dev.off()



pl2 <- pl + scale_fill_manual(values = c('Reg_1_1'    = "grey40",
                                         'Reg_2_2.5'  = "grey40",
                                         "Reg_3_5" = "grey40",
                                         "Reg_4_15" = "grey40",
                                         "Reg_5_30" = "grey40",
                                         
                                         'Reg_1_2.5'  = "grey40",
                                         "Reg_2_5" = "grey40",
                                         "Reg_3_15" = "grey40",
                                         "Reg_4_30" = "grey40",
                                         
                                         "Reg_1_5" = "darkgreen",
                                         "Reg_2_15" = "darkgreen",
                                         "Reg_3_30" = "darkgreen",
                                         
                                         "Reg_1_15" = "grey40",
                                         "Reg_2_30" = "grey40",
                                         
                                         "Reg_1_30" = "grey40"
) )
pl2
png(filename=paste0(cwd, "/20230926_sankey_phosphosite_regulation_GREEN_ONLY.png "), height=2000, width=2000, res=300)
  pl2
dev.off()



pl2 <- pl + scale_fill_manual(values = c('Reg_1_1'    = "grey40",
                                         'Reg_2_2.5'  = "grey40",
                                         "Reg_3_5" = "grey40",
                                         "Reg_4_15" = "grey40",
                                         "Reg_5_30" = "grey40",
                                         
                                         'Reg_1_2.5'  = "grey40",
                                         "Reg_2_5" = "grey40",
                                         "Reg_3_15" = "grey40",
                                         "Reg_4_30" = "grey40",
                                         
                                         "Reg_1_5" = "grey40",
                                         "Reg_2_15" = "grey40",
                                         "Reg_3_30" = "grey40",
                                         
                                         "Reg_1_15" = "deeppink",
                                         "Reg_2_30" = "deeppink",
                                         
                                         "Reg_1_30" = "grey40"
) )
pl2
png(filename=paste0(cwd, "/20230926_sankey_phosphosite_regulation_PINK_ONLY.png "), height=2000, width=2000, res=300)
  pl2
dev.off()



pl2 <- pl + scale_fill_manual(values = c('Reg_1_1'    = "grey40",
                                         'Reg_2_2.5'  = "grey40",
                                         "Reg_3_5" = "grey40",
                                         "Reg_4_15" = "grey40",
                                         "Reg_5_30" = "grey40",
                                         
                                         'Reg_1_2.5'  = "grey40",
                                         "Reg_2_5" = "grey40",
                                         "Reg_3_15" = "grey40",
                                         "Reg_4_30" = "grey40",
                                         
                                         "Reg_1_5" = "grey40",
                                         "Reg_2_15" = "grey40",
                                         "Reg_3_30" = "grey40",
                                         
                                         "Reg_1_15" = "grey40",
                                         "Reg_2_30" = "grey40",
                                         
                                         "Reg_1_30" = "blue"
) )
pl2
png(filename=paste0(cwd, "/20230926_sankey_phosphosite_regulation_BLUE_ONLY.png "), height=2000, width=2000, res=300)
  pl2
dev.off()



pl2 <- pl + scale_fill_manual(values = c('Reg_1_1'    = "grey40",
                                         'Reg_2_2.5'  = "grey40",
                                         "Reg_3_5" = "grey40",
                                         "Reg_4_15" = "grey40",
                                         "Reg_5_30" = "grey40",
                                         
                                         'Reg_1_2.5'  = "grey40",
                                         "Reg_2_5" = "grey40",
                                         "Reg_3_15" = "grey40",
                                         "Reg_4_30" = "grey40",
                                         
                                         "Reg_1_5" = "grey40",
                                         "Reg_2_15" = "grey40",
                                         "Reg_3_30" = "grey40",
                                         
                                         "Reg_1_15" = "grey40",
                                         "Reg_2_30" = "grey40",
                                         
                                         "Reg_1_30" = "grey40"
) )
pl2
png(filename=paste0(cwd, "/20230926_sankey_phosphosite_regulation_GREY_ONLY.png "), height=2000, width=2000, res=300)
  pl2
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

library("ReactomePA")
#library("org.Hs.eg.db")
#library("biomaRt")

#n.paths <- 60 #both
#n.paths <- 282 #pos
#n.paths <- 98 #neg
n.paths <- 500 #pos/neg for pep/prot count barplots

#Early <- unique(c(reg.peps[[1]], reg.peps[[2]]))
#Intermediate <- unique(c(reg.peps[[3]], reg.peps[[4]]))
#Late <- unique(c(reg.peps[[5]], reg.peps[[6]]))
#write.table(x=Early, file=paste0(cwd, "/Reactome/reg_peps_early.txt"), row.names=F, col.names=F, sep="\t")
#write.table(x=Intermediate, file=paste0(cwd, "/Reactome/reg_peps_intermediate.txt"), row.names=F, col.names=F, sep="\t")
#write.table(x=Late, file=paste0(cwd, "/Reactome/reg_peps_late.txt"), row.names=F, col.names=F, sep="\t")

Early <- unique(c(reg.peps.pos[[1]], reg.peps.pos[[2]]))
Intermediate <- unique(c(reg.peps.pos[[3]], reg.peps.pos[[4]]))
Late <- unique(c(reg.peps.pos[[5]], reg.peps.pos[[6]]))
write.table(x=Early, file=paste0(cwd, "/Reactome/reg_peps_early_pos.txt"), row.names=F, col.names=F, sep="\t")
write.table(x=Intermediate, file=paste0(cwd, "/Reactome/reg_peps_intermediate_pos.txt"), row.names=F, col.names=F, sep="\t")
write.table(x=Late, file=paste0(cwd, "/Reactome/reg_peps_late_pos.txt"), row.names=F, col.names=F, sep="\t")

#Early <- unique(c(reg.peps.neg[[1]], reg.peps.neg[[2]]))
#Intermediate <- unique(c(reg.peps.neg[[3]], reg.peps.neg[[4]]))
#Late <- unique(c(reg.peps.neg[[5]], reg.peps.neg[[6]]))
#write.table(x=Early, file=paste0(cwd, "/Reactome/reg_peps_early_neg.txt"), row.names=F, col.names=F, sep="\t")
#write.table(x=Intermediate, file=paste0(cwd, "/Reactome/reg_peps_intermediate_neg.txt"), row.names=F, col.names=F, sep="\t")
#write.table(x=Late, file=paste0(cwd, "/Reactome/reg_peps_late_neg.txt"), row.names=F, col.names=F, sep="\t")


reg.prot.info.early <- prepareProteinInfo(pep.names=Early, prot.info=prot.dat.orig)
reg.prot.info.intermediate <- prepareProteinInfo(pep.names=Intermediate, prot.info=prot.dat.orig)
reg.prot.info.late <- prepareProteinInfo(pep.names=Late, prot.info=prot.dat.orig)
prot.accessions <- list(
                        Early=reg.prot.info.early$`UniProt ID`,
                        Intermediate=reg.prot.info.intermediate$`UniProt ID`,
                        Late=reg.prot.info.late$`UniProt ID`
                    )
entrez.gene.ids.all <- returnUniProtEntrezGeneIDs(uniprot.ids=unique(c(prot.accessions[[1]],prot.accessions[[2]],prot.accessions[[3]])),
                           table=FALSE,
                           source="E:/Projekte/Phospho2.1/ID-Mapping/20220702_Phospho2.1_UniProt_to_GeneID.tsv",
                           na.rm=TRUE,
                           file.name="entrez.genes.id_all",
                           output.path=paste0(cwd, "/Reactome"))
#overlap <- calculate.overlap(list(Early=Early, Intermediate=Intermediate, Late=Late))
reactome.res <- list()
entrez.gene.ids <- list()
#mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
for(i in 1:length(prot.accessions)){
    #entrez.gene.ids <- getBM(
    #                      attributes=c('entrezgene_id'),
    #                      filters='uniprotswissprot',
    #                      values=unique(gsub("\\-\\d+", "", x=prot.accessions[[i]])),
    #                      mart=mart
    #                    )
    #entrez.gene.ids[[i]] <- returnUniProtEntrezGeneIDs(uniprot.ids=prot.accessions[[i]], file.name=paste0("entrez.genes.id_", names(prot.accessions)[i]), output.path=paste0(cwd, "/Reactome"))
    entrez.gene.ids[[i]] <- returnUniProtEntrezGeneIDs(uniprot.ids=prot.accessions[[i]],
                               table=FALSE,
                               source="E:/Projekte/Phospho2.1/ID-Mapping/20220702_Phospho2.1_UniProt_to_GeneID.tsv",
                               na.rm=TRUE,
                               file.name=paste0("entrez.genes.id_", names(prot.accessions)[i]), output.path=paste0(cwd, "/Reactome"))
    reactome.res[[i]] <- enrichPathway(
                             entrez.gene.ids[[i]],
                             #entrez.gene.ids$entrezgene_id,
                             organism="human",
                             #pvalueCutoff=0.05,
                             pvalueCutoff=1.0,
                             pAdjustMethod = "BH",
                             qvalueCutoff = 1.0,
                             minGSSize = 10,
                             maxGSSize = 707
                         )
    print(entrez.gene.ids[[i]])
}

unique.term.number <- length(unique(c(reactome.res[[1]][1:n.paths,"ID"], reactome.res[[2]][1:n.paths,"ID"], reactome.res[[3]][1:n.paths,"ID"])))

write.table(x=reactome.res[[1]][1:n.paths,], file=paste0(cwd, "/Reactome/reactome_ORA_early_top", n.paths, ".txt"), sep="\t", row.names=FALSE, col.names=TRUE)
write.table(x=reactome.res[[2]][1:n.paths,], file=paste0(cwd, "/Reactome/reactome_ORA_intermediate_top", n.paths, ".txt"), sep="\t", row.names=FALSE, col.names=TRUE)
write.table(x=reactome.res[[3]][1:n.paths,], file=paste0(cwd, "/Reactome/reactome_ORA_late_top", n.paths, ".txt"), sep="\t", row.names=FALSE, col.names=TRUE)

#------

term.counts <- matrix(nrow=unique.term.number, ncol=19)
term.ids <- unique(c(reactome.res[[1]][1:n.paths,"ID"],
                            reactome.res[[2]][1:n.paths,"ID"],
                            reactome.res[[3]][1:n.paths,"ID"]))
rownames(term.counts) <- term.ids
colnames(term.counts) <- c("ID",
                           "Description",
                           "pvalueEarly",
                           "pvalueIntermediate",
                           "pvalueLate",
                           "p.adjustEarly",
                           "p.adjustIntermediate",
                           "p.adjustLate",
                           "BackgroundCount",
                           "TermCount",
                           "RegGeneCountEarly",
                           "RegGeneCountIntermediate",
                           "RegGeneCountLate",
                           "CountEarly",
                           "CountIntermediate",
                           "CountLate",
                           "RatioEarly",
                           "RatioIntermediate",
                           "RatioLate")
backgr.count <- gsub("\\d+\\/(\\d+)", "\\1", reactome.res[[1]][1,"BgRatio"])
reg.genes.early <- as.numeric(gsub("\\d+\\/(\\d+)", "\\1", reactome.res[[1]][1,"GeneRatio"]))
reg.genes.intermediate <- as.numeric(gsub("\\d+\\/(\\d+)", "\\1", reactome.res[[2]][1,"GeneRatio"]))
reg.genes.late <- as.numeric(gsub("\\d+\\/(\\d+)", "\\1", reactome.res[[3]][1,"GeneRatio"]))

for(i in 1:unique.term.number){
  if(term.ids[i] %in% reactome.res[[1]][,"ID"]){
    term.counts[term.ids[i],"ID"] <- reactome.res[[1]][term.ids[i],"ID"]
  }else if(term.ids[i] %in% reactome.res[[2]][,"ID"]){
    term.counts[term.ids[i],"ID"] <- reactome.res[[2]][term.ids[i],"ID"]  
  }else if(term.ids[i] %in% reactome.res[[3]][,"ID"]){
    term.counts[term.ids[i],"ID"] <- reactome.res[[3]][term.ids[i],"ID"]  
  }
  
  
  if(term.ids[i] %in% reactome.res[[1]][,"ID"]){
    term.counts[term.ids[i],"Description"] <- reactome.res[[1]][term.ids[i],"Description"]
  }else if(term.ids[i] %in% reactome.res[[2]][,"ID"]){
    term.counts[term.ids[i],"Description"] <- reactome.res[[2]][term.ids[i],"Description"]  
  }else if(term.ids[i] %in% reactome.res[[3]][,"ID"]){
    term.counts[term.ids[i],"Description"] <- reactome.res[[3]][term.ids[i],"Description"]  
  }
  
  

  if(term.ids[i] %in% reactome.res[[1]][,"ID"]){
    term.counts[term.ids[i],"pvalueEarly"] <- reactome.res[[1]][term.ids[i],"pvalue"]
  }else{
    term.counts[term.ids[i],"pvalueEarly"] <- "1"
  }
  if(term.ids[i] %in% reactome.res[[2]][,"ID"]){
    term.counts[term.ids[i],"pvalueIntermediate"] <- reactome.res[[2]][term.ids[i],"pvalue"]  
  }else{
    term.counts[term.ids[i],"pvalueIntermediate"] <- "1"
  }
  if(term.ids[i] %in% reactome.res[[3]][,"ID"]){
    term.counts[term.ids[i],"pvalueLate"] <- reactome.res[[3]][term.ids[i],"pvalue"]  
  }else{
    term.counts[term.ids[i],"pvalueLate"] <- "1"
  }
  
  
  
  if(term.ids[i] %in% reactome.res[[1]][,"ID"]){
    term.counts[term.ids[i],"p.adjustEarly"] <- reactome.res[[1]][term.ids[i],"p.adjust"]
  }else{
    term.counts[term.ids[i],"p.adjustEarly"] <- "1"
  }
  if(term.ids[i] %in% reactome.res[[2]][,"ID"]){
    term.counts[term.ids[i],"p.adjustIntermediate"] <- reactome.res[[2]][term.ids[i],"p.adjust"]  
  }else{
    term.counts[term.ids[i],"p.adjustIntermediate"] <- "1"
  }
  if(term.ids[i] %in% reactome.res[[3]][,"ID"]){
    term.counts[term.ids[i],"p.adjustLate"] <- reactome.res[[3]][term.ids[i],"p.adjust"]  
  }else{
    term.counts[term.ids[i],"p.adjustLate"] <- "1"
  }
  
  
  
  term.counts[term.ids[i],"BackgroundCount"] <- backgr.count
  
  
  if(term.ids[i] %in% reactome.res[[1]][,"ID"]){
    term.counts[term.ids[i],"TermCount"] <- gsub("(\\d+)\\/\\d+", "\\1", reactome.res[[1]][term.ids[i],"BgRatio"])
  }else if(term.ids[i] %in% reactome.res[[2]][,"ID"]){
    term.counts[term.ids[i],"TermCount"] <- gsub("(\\d+)\\/\\d+", "\\1", reactome.res[[2]][term.ids[i],"BgRatio"])  
  }else if(term.ids[i] %in% reactome.res[[3]][,"ID"]){
    term.counts[term.ids[i],"TermCount"] <- gsub("(\\d+)\\/\\d+", "\\1", reactome.res[[3]][term.ids[i],"BgRatio"])  
  }
  
  
  if(term.ids[i] %in% reactome.res[[1]][,"ID"]){
    term.counts[term.ids[i],"RegGeneCountEarly"] <- gsub("\\d+\\/(\\d+)", "\\1", reactome.res[[1]][term.ids[i],"GeneRatio"])
  }else{
    term.counts[term.ids[i],"RegGeneCountEarly"] <- reg.genes.early
  }
  if(term.ids[i] %in% reactome.res[[2]][,"ID"]){
    term.counts[term.ids[i],"RegGeneCountIntermediate"] <- gsub("\\d+\\/(\\d+)", "\\1", reactome.res[[2]][term.ids[i],"GeneRatio"])
  }else{
    term.counts[term.ids[i],"RegGeneCountIntermediate"] <- reg.genes.intermediate
  }
  if(term.ids[i] %in% reactome.res[[3]][,"ID"]){
    term.counts[term.ids[i],"RegGeneCountLate"] <- gsub("\\d+\\/(\\d+)", "\\1", reactome.res[[3]][term.ids[i],"GeneRatio"])
  }else{
    term.counts[term.ids[i],"RegGeneCountLate"] <- reg.genes.late
  }
  
  
  if(term.ids[i] %in% reactome.res[[1]][,"ID"]){
    term.counts[term.ids[i],"CountEarly"] <- reactome.res[[1]][term.ids[i],"Count"]
  }else{
    term.counts[term.ids[i],"CountEarly"] <- "0"
  }
  if(term.ids[i] %in% reactome.res[[2]][,"ID"]){
    term.counts[term.ids[i],"CountIntermediate"] <- reactome.res[[2]][term.ids[i],"Count"]
  }else{
    term.counts[term.ids[i],"CountIntermediate"] <- "0"
  }
  if(term.ids[i] %in% reactome.res[[3]][,"ID"]){
    term.counts[term.ids[i],"CountLate"] <- reactome.res[[3]][term.ids[i],"Count"]
  }else{
    term.counts[term.ids[i],"CountLate"] <- "0"
  }
  

  term.counts[term.ids[i],"RatioEarly"] <- as.numeric(term.counts[term.ids[i],"CountEarly"]) / as.numeric(term.counts[term.ids[i],"RegGeneCountEarly"])
  term.counts[term.ids[i],"RatioIntermediate"] <- as.numeric(term.counts[term.ids[i],"CountIntermediate"]) / as.numeric(term.counts[term.ids[i],"RegGeneCountIntermediate"])
  term.counts[term.ids[i],"RatioLate"] <- as.numeric(term.counts[term.ids[i],"CountLate"]) / as.numeric(term.counts[term.ids[i],"RegGeneCountLate"])
}
write.table(x=term.counts, file=paste0(cwd, "/Reactome/term.counts.txt"), sep="\t", row.names=FALSE, col.names=TRUE)

#----------
# UPDATED BEGIN !!!
#----------
 timeperiod.lab <- c("Early", "Intermediate", "Late")
 for(i in 1:3){
   output.tmp <- matrix(ncol=9)
   colnames(output.tmp) <- c("Reactome-ID",
                             "Description",
                             "GeneRatio",
                             "BgRatio",
                             "p-value",
                             "p.adjust",
                             "Entrez-ID",
                             "UniProt-ID",
                             "GeneName")
   
   uniprot.ids.updated <- updateUniProtIDs(uniprot.ids=prot.accessions[[i]],
                                           source="E:/Projekte/Phospho2.1/ID-Mapping/20220927_Phospho2.1_UniProt_to_UniProt.tsv")
   
   #for(j in 1:2){
   for(j in 1:n.paths){
     entrez.vec <- unlist(strsplit(reactome.res[[i]]$geneID[j], "/"))
     uniprot.vec <- returnEntrezUniProtIDs(entrez.ids=entrez.vec,
                                           source="E:/Projekte/Phospho2.1/ID-Mapping/20220702_Phospho2.1_UniProt_to_GeneID.tsv")
     tmp.overlap <- vecsets::vintersect(c(uniprot.ids.updated, uniprot.vec[duplicated(uniprot.vec)]), uniprot.vec, multiple=TRUE)
     if(length(tmp.overlap) == length(entrez.vec)){
       print("OK!")
     }else{
       print("PROBLEM!")
       print(paste0("entrez.vec: ", length(entrez.vec)))
       print(paste0("tmp.overlap: ", length(tmp.overlap)))
       break
     }
     gene.name.tab <- returnUniProtGeneNames(uniprot.ids=tmp.overlap,
                                             table=TRUE,
                                             source="E:/Projekte/Phospho2.1/ID-Mapping/20220702_Phospho2.1_UniProt_to_GeneName.tsv")
     rownames(gene.name.tab) <- gene.name.tab[,"From"]
     gene.name.vec <- c()
     for(k in 1:length(entrez.vec)){
         gene.name.vec[k] <- gene.name.tab[tmp.overlap[k], "To"]   
     }
     for(k in 1:length(entrez.vec)){
         output.tmp <- rbind(output.tmp, c(
                         reactome.res[[i]]$ID[j],
                         reactome.res[[i]]$Description[j],
                         reactome.res[[i]]$GeneRatio[j],
                         reactome.res[[i]]$BgRatio[j],
                         reactome.res[[i]]$pvalue[j],
                         reactome.res[[i]]$p.adjust[j],
                         entrez.vec[k],
                         tmp.overlap[k],
                         gene.name.vec[k]
                       ))
     }
   }
   if(sum(is.na(output.tmp[1,])) == ncol(output.tmp)) output.tmp <- output.tmp[-1,]
   #write.table(x=output.tmp, file=paste0(cwd, "/Reactome/reactome_ORA_", timeperiod.lab[i], "_complete_top", n.paths, ".txt"), sep="\t", row.names=FALSE, col.names=TRUE)
   write.table(x=output.tmp, file=paste0(cwd, "/Reactome/reactome_ORA_", timeperiod.lab[i], "_complete_top", n.paths, "_pos.txt"), sep="\t", row.names=FALSE, col.names=TRUE)
   #write.table(x=output.tmp, file=paste0(cwd, "/Reactome/reactome_ORA_", timeperiod.lab[i], "_complete_top", n.paths, "_neg.txt"), sep="\t", row.names=FALSE, col.names=TRUE)
 }
#----------
# UPDATED END !!!
#----------


del.idx1 <- c(
              grep("GTPase", x=reactome.res[[1]][1:n.paths,"Description"]),
              grep("Transport of Mature mRNA", x=reactome.res[[1]][1:n.paths,"Description"]),
              grep("Apoptotic", x=reactome.res[[1]][1:n.paths,"Description"]),
              grep("Programmed Cell Death", x=reactome.res[[1]][1:n.paths,"Description"]),
              grep("Signaling by BRAF and RAF fusions", x=reactome.res[[1]][1:n.paths,"Description"])
            )
del.idx2 <- c(
              grep("GTPase", x=reactome.res[[2]][1:n.paths,"Description"]),
              grep("Transport of Mature mRNA", x=reactome.res[[2]][1:n.paths,"Description"]),
              grep("Apoptotic", x=reactome.res[[2]][1:n.paths,"Description"]),
              grep("Programmed Cell Death", x=reactome.res[[2]][1:n.paths,"Description"]),
              grep("Signaling by BRAF and RAF fusions", x=reactome.res[[2]][1:n.paths,"Description"])
            )
del.idx3 <- c(
              grep("GTPase", x=reactome.res[[3]][1:n.paths,"Description"]),
              grep("Transport of Mature mRNA", x=reactome.res[[3]][1:n.paths,"Description"]),
              grep("Apoptotic", x=reactome.res[[3]][1:n.paths,"Description"]),
              grep("Programmed Cell Death", x=reactome.res[[3]][1:n.paths,"Description"]),
              grep("Signaling by BRAF and RAF fusions", x=reactome.res[[2]][1:n.paths,"Description"])
            )
keep1 <- c(
              grep("^RAC1 GTPase cycle$", x=reactome.res[[1]][1:n.paths,"Description"]),
              grep("^CDC42 GTPase cycle$", x=reactome.res[[1]][1:n.paths,"Description"]),
              grep("^Signaling by Rho GTPases$", x=reactome.res[[1]][1:n.paths,"Description"]),
              grep("Transport of Mature Transcript to Cytoplasm", x=reactome.res[[1]][1:n.paths,"Description"])
         )
keep2 <- c(
              grep("^RAC1 GTPase cycle$", x=reactome.res[[2]][1:n.paths,"Description"]),
              grep("^CDC42 GTPase cycle$", x=reactome.res[[2]][1:n.paths,"Description"]),
              grep("^Signaling by Rho GTPases$", x=reactome.res[[2]][1:n.paths,"Description"]),
              grep("Transport of Mature Transcript to Cytoplasm", x=reactome.res[[2]][1:n.paths,"Description"])
         )
keep3 <- c(
              grep("^RAC1 GTPase cycle$", x=reactome.res[[3]][1:n.paths,"Description"]),
              grep("^CDC42 GTPase cycle$", x=reactome.res[[3]][1:n.paths,"Description"]),
              grep("^Signaling by Rho GTPases$", x=reactome.res[[3]][1:n.paths,"Description"]),
              grep("Transport of Mature Transcript to Cytoplasm", x=reactome.res[[3]][1:n.paths,"Description"])
         )
del.idx1 <- setdiff(del.idx1, keep1)
del.idx2 <- setdiff(del.idx2, keep2)
del.idx3 <- setdiff(del.idx3, keep3)
terms.idx1 <- setdiff(1:n.paths, del.idx1)
terms.idx2 <- setdiff(1:n.paths, del.idx2)
terms.idx3 <- setdiff(1:n.paths, del.idx3)


terms <- unique(c(
              reactome.res[[1]]$Description[terms.idx1],
              reactome.res[[2]]$Description[terms.idx2],
              reactome.res[[3]]$Description[terms.idx3]
            ))
terms <- sort(terms)
pathway.enrichment.dat <- matrix(nrow=length(terms), ncol=length(prot.accessions))
rownames(pathway.enrichment.dat) <- terms
colnames(pathway.enrichment.dat) <- c("Early", "Intermediate", "Late")
for(i in 1:length(terms)){
  if(length(intersect(reactome.res[[1]]$Description, terms[i])) > 0){
      ped.tmp1 <- min(reactome.res[[1]][reactome.res[[1]]$Description %in% terms[i],"p.adjust"]) #--> min-workaround for redundand pathways (R-HSA-9694631, R-HSA-9683610)    
  }else{
      ped.tmp1 <- NA  
  }
  
  if(length(intersect(reactome.res[[2]]$Description, terms[i])) > 0){
    ped.tmp2 <- min(reactome.res[[2]][reactome.res[[2]]$Description %in% terms[i],"p.adjust"]) #--> min-workaround for redundand pathways (R-HSA-9694631, R-HSA-9683610)   
  }else{
    ped.tmp2 <- NA  
  }
  
  if(length(intersect(reactome.res[[3]]$Description, terms[i])) > 0){
    ped.tmp3 <- min(reactome.res[[3]][reactome.res[[3]]$Description %in% terms[i],"p.adjust"]) #--> min-workaround for redundand pathways (R-HSA-9694631, R-HSA-9683610)   
  }else{
    ped.tmp3 <- NA  
  }
  
  #print(intersect(reactome.res[[3]]$Description, terms[i]))
  #print(ped.tmp1)
  #print(ped.tmp2)
  #print(ped.tmp3)
  pathway.enrichment.dat[i,] <- c(ped.tmp1, ped.tmp2, ped.tmp3) 
  #print("----------------------------------")
}

#rownames(pathway.enrichment.dat)[grep("Apoptosis", rownames(pathway.enrichment.dat))] <- gsub("Apoptosis", "Apoptosis*", rownames(pathway.enrichment.dat)[grep("Apoptosis", rownames(pathway.enrichment.dat))])
#rownames(pathway.enrichment.dat)[grep("Transport of Mature Transcript to Cytoplasm", rownames(pathway.enrichment.dat))] <- gsub("Transport of Mature Transcript to Cytoplasm", "Transport of Mature Transcript to Cytoplasm*", rownames(pathway.enrichment.dat)[grep("Transport of Mature Transcript to Cytoplasm", rownames(pathway.enrichment.dat))])
#rownames(pathway.enrichment.dat)[grep("Signaling by Rho GTPases", rownames(pathway.enrichment.dat))] <- gsub("Signaling by Rho GTPases", "Signaling by Rho GTPases*", rownames(pathway.enrichment.dat)[grep("Signaling by Rho GTPases", rownames(pathway.enrichment.dat))])
pathway.enrichment.dat[is.na(pathway.enrichment.dat)] <- 1
png(paste0(cwd, "/Reactome/20230907_heatmap_reactome_enriched_pathways.png"),  width=4000, height=4400, res=300)
    #par(mar=c(0,10.1,1.1,10.1), oma=c(10,1,1,1))
    par(cex.main=2.25, cex.lab=1.5, cex.axis=1.5)
    heatmap.2(-log10(pathway.enrichment.dat), 
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          hclustfun = function(x) hclust(x, method="average"),
          col = colorpanel(100,"blue","white","red"),
          margins = c(10, 35),
          trace = "none", 
          lhei = c(2, 10),
          scale = c("row"),
          na.color="white",
          cexRow = 1.1,
          cexCol = 2,
          main = "Enriched pathways", 
          Colv = FALSE,
          Rowv = TRUE,
          sepwidth=c(0.001,0.001),
          sepcolor="grey",
          colsep=0:ncol(pathway.enrichment.dat),
          rowsep=0:nrow(pathway.enrichment.dat),
          key=T,
          keysize=1.0,
          key.xlab = "Z-score of -log10(FDR)",
          key.title = NA,
          key.ylab = NA
    )
dev.off()

heatmap.dat.out <- cbind(rownames(pathway.enrichment.dat),-log10(pathway.enrichment.dat))
write.table(x=heatmap.dat.out, file=paste0(cwd, "/Reactome/reactome_ORA_heatmap_data", n.paths, ".txt"), sep="\t", row.names=FALSE, col.names=TRUE)



#pathway.clust <- read.table(file=paste0(cwd, "/Reactome/20220211_reactome-heatmap-manual_sorting.txt"), header=TRUE, sep="\t", quote = "\"")
pathway.clust <- read.table(file=paste0(cwd, "/Reactome/20220810_Reactome_manual_sorting_manual_selection.txt"), header=TRUE, sep="\t", quote = "\"")
manual.sort <- pathway.clust[,1]
rownames(pathway.enrichment.dat) <- gsub("^Translocation of SLC2A4 \\(GLUT4\\) to the plasma membrane$", "Translocation of SLC2A4 \\(GLUT4\\) to plasma membrane", rownames(pathway.enrichment.dat))
manual.sort <- gsub("^Translocation of SLC2A4 \\(GLUT4\\) to the plasma membrane$", "Translocation of SLC2A4 \\(GLUT4\\) to plasma membrane", manual.sort)
png(paste0(cwd, "/Reactome/20230907_heatmap_reactome_enriched_pathways_ManualSorting.png"),  width=4000, height=4400, res=300)
    par(font=2, font.axis=2)
    heatmap.2(-log10(pathway.enrichment.dat[manual.sort,]), 
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          hclustfun = function(x) hclust(x, method="average"),
          dendrogram = "none",
          col = colorpanel(100,"blue","white","red"),
          #margins = c(33, 65), # for 3 arguments 
          margins = c(26, 55), # for 2 arguments 
          trace = "none", 
          lhei = c(0.5, 10),
          #lwid = c(1, 6, 1), # 3 arguments if color-code side bar
          lwid = c(0.5, 6), # col width: 2 arguments if color-code side bar
          scale = c("row"),
          na.color="white",
          cexRow = 2.3,
          cexCol = 3.5,
          adjCol = c(NA,0.5),
          #main = "Enriched pathways",
          colRow = pathway.clust$Color,
          Colv = FALSE,
          Rowv = FALSE,
          #RowSideColors=pathway.clust$Color, # for side bar showing color-code for term categories
          sepwidth=c(0.001,0.001),
          sepcolor="grey",
          colsep=0:ncol(pathway.enrichment.dat),
          rowsep=0:nrow(pathway.enrichment.dat),
          key=F,
          keysize=1.0,
          key.xlab = "Z-score of -log10(FDR)",
          key.title = NA,
          key.ylab = NA
    )
    
    legend("bottomright",
           inset = c(0.125,0), # x,y-shift from named position
           legend = unique(pathway.clust$Pathway.cluster),
           col = unique(pathway.clust$Color),
           pch = 15,
           bty = "n",
           pt.cex = 4,
           cex=1.8
    )
dev.off()

#-----------

png(paste0(cwd, "/Reactome/20230907_heatmap_reactome_enriched_pathways_ManualSorting_withKey.png"),  width=4000, height=4400, res=300)
par(font=2, font.axis=2)
heatmap.2(-log10(pathway.enrichment.dat[manual.sort,]), 
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          hclustfun = function(x) hclust(x, method="average"),
          dendrogram = "none",
          col = colorpanel(100,"blue","white","red"),
          #margins = c(33, 65), # for 3 arguments 
          #margins = c(26, 55), # for 2 arguments 
          #margins = c(1, 1), # for 2 arguments 
          trace = "none", 
          lhei = c(1.2, 10),
          #lwid = c(1, 1, 1), # 3 arguments if color-code side bar
          lwid = c(1, 1.2), # col width: 2 arguments if color-code side bar
          scale = c("row"),
          na.color="white",
          cexRow = 2.3,
          cexCol = 3.5,
          adjCol = c(NA,0.5),
          #main = "Enriched pathways",
          colRow = pathway.clust$Color,
          Colv = FALSE,
          Rowv = FALSE,
          #RowSideColors=pathway.clust$Color, # for side bar showing color-code for term categories
          sepwidth=c(0.001,0.001),
          sepcolor="grey",
          colsep=0:ncol(pathway.enrichment.dat),
          rowsep=0:nrow(pathway.enrichment.dat),
          key=T,
          keysize=1,
          density.info="none",
          key.par = list("cex" = 1.25),
          key.title = "Z-score of -log10(FDR)",
          key.xlab = NA,
          #key.xlab = "Z-score of -log10(FDR)",
          #key.title = NA,
          key.ylab = NA
)

legend("bottomright",
       inset = c(0.125,0), # x,y-shift from named position
       legend = unique(pathway.clust$Pathway.cluster),
       col = unique(pathway.clust$Color),
       pch = 15,
       bty = "n",
       pt.cex = 4,
       cex=1.8
)
dev.off()

#-----------

revigo.results <- read.table(file=paste0(cwd,"/PANTHER/Revigo_CC_OnScreenTable.tsv"),header=T,sep="\t")
revigo.term.ids <- revigo.results$TermID

giveTermClassi(file.path=paste0(cwd, "/PANTHER/pantherGeneList_CC_Early.txt"),
               output.path=paste0(cwd, "/PANTHER/term_share_CC_Early.txt"))
giveTermClassi(file.path=paste0(cwd, "/PANTHER/pantherGeneList_CC_Intermediate.txt"),
               output.path=paste0(cwd, "/PANTHER/term_share_CC_Intermediate.txt"))
giveTermClassi(file.path=paste0(cwd, "/PANTHER/pantherGeneList_CC_Late.txt"),
               output.path=paste0(cwd, "/PANTHER/term_share_CC_Late.txt"))

giveTermClassi(file.path=paste0(cwd, "/PANTHER/pantherGeneList_CC_all.txt"),
               term.filter=revigo.term.ids,
               output.path=paste0(cwd, "/PANTHER/term_share_CC_all_revigo.txt"))

#-----------

for(i in 1:nrow(term.counts)){
  
  if(length(intersect(manual.sort, term.counts[i,"Description"])) > 0){
    
    #print(term.counts[i,"Description"])
    png(paste0(cwd, "/Reactome/CountBarplots/barplot_", i, ".png"),  width=4000, height=4400, res=300)
      par(cex.main=3)
      barplot(as.numeric(term.counts[i,c("CountEarly","CountIntermediate","CountLate")]),
        main=paste0(term.counts[i,"Description"], " (", term.counts[i,"CountEarly"], ", ", term.counts[i,"CountIntermediate"], ", ", term.counts[i,"CountLate"], ")"),
        names.arg=c("Early", "Intermediate", "Late"),
        cex.axis=3,
        cex.names=3)
      dev.off()
    
  }
  
}

#-------------------------------------------------------------------------------

#sig.reg.pept <- unique(c(reg.peps.pos.na[[1]], 
#                         reg.peps.pos.na[[2]],
#                         reg.peps.pos.na[[3]],
#                         reg.peps.pos.na[[4]],
#                         reg.peps.pos.na[[5]],
#                         reg.peps.pos.na[[6]]))

sig.reg.pept <- reg.peps.all.na

#clueR.dat <- prepareClueRDat(sig.reg.pept=sig.reg.pept,
#                             dat.ratio.used=dat.ratio.used,
#                             mods.master.used=mods.master.used,
#                             psm.number.used=psm.number.used,
#                             organism="human",
#                             log=TRUE,
#                             ratios=TRUE,
#                             output.path=cwd,
#                             file.name="clueR.dat2")

#clueR.dat <- prepareClueRDat(sig.reg.pept=sig.reg.pept,
#                             dat.ratio.used=dat.ratio,
#                             mods.master.used=mods.master,
#                             psm.number.used=psm.number,
#                             organism="human",
#                             log=TRUE,
#                             ratios=TRUE,
#                             output.path=cwd,
#                             file.name="clueR.dat5")
#clueR.dat[is.na(clueR.dat)] <- 0

#a <- rownames(dat.abundance.gr.used)[!rowSums(!is.finite(dat.abundance.gr.used))]
#b <- intersect(a, sig.reg.pept)
#clueR.dat <- prepareClueRDat(sig.reg.pept=b,
#                             dat.ratio.used=dat.abundance.gr.used,
#                             mods.master.used=mods.master.used,
#                             psm.number.used=psm.number.used,
#                             organism="human",
#                             log=TRUE,
#                             ratios=FALSE,
#                             output.path=cwd,
#                             file.name="clueR.dat3")


clueR.dat <- prepareClueRDat(sig.reg.pept=sig.reg.pept,
                             dat.ratio.used=dat.abundance.gr,
                             mods.master.used=mods.master,
                             psm.number.used=psm.number,
                             organism="human",
                             log=TRUE,
                             ratios=FALSE,
                             output.path=cwd,
                             file.name="clueR.dat4")
clueR.dat[is.na(clueR.dat)] <- 0

number.PSPkinases <- length(PhosphoSite.human)
number.PSPannotations <- sum(lengths(PhosphoSite.human))
#number.PSPannotations <- length(unlist(PhosphoSite.human)) # --> same as above
number.PSPphossites <- length(unique(unlist(PhosphoSite.human)))
number.phossites <- length(unique(rownames(clueR.dat)))
number.PSPphossites.overlap <- length(intersect(unique(unlist(PhosphoSite.human)), unique(rownames(clueR.dat))))

#best <- performClueR(clueR.dat=clueR.dat,
#                     organism="human",
#                     rep=20,
#                     kRange=2:10,
#                     effectiveSize=c(5,100),
#                     pvalueCutoff=0.07,
#                     output.path=cwd)

best <- performClueR(clueR.dat=clueR.dat,
                     organism="human",
                     rep=20,
                     kRange=2:10,
                     effectiveSize=c(5,100),
                     pvalueCutoff=0.07,
                     output.path=cwd)

plotClusterClueR(clueR.dat=clueR.dat,
                 best=best,
                 scale=TRUE,
                 unlog=FALSE,
                 output.path=cwd)

#plotClusterClueR(clueR.dat=clueR.dat,
#                 best=best,
#                 scale=FALSE,
#                 unlog=FALSE,
#                 output.path=cwd)

#plotClusterClueR(clueR.dat=clueR.dat,
#                 best=best,
#                 scale=FALSE,
#                 unlog=TRUE,
#                 output.path=cwd)

#-----------------------------------------------

#devtools::install_github("PengyiYang/ClueR")
library(ClueR)
library(PhosR)
data(PhosphoSitePlus)
load(paste0(cwd, "/ClueR/ClueR-K/clueObj.RData"))

reg.phospeps <- vector(mode="list", length=6)
reg.prot.info.used <- vector(mode="list", length=6)
ksea.results <- vector(mode="list", length=6)
for(i in 1:6){
  #na.idx <- which(is.na(p.values.used[,i]))
  #reg <- reg.peps.idx.na[[i]]
  #if(length(na.idx) > 0) reg[na.idx] <- TRUE
  #reg.phospeps[[i]] <- dat.ratio.used[reg,]
  #reg.phospeps[[i]] <- dat.ratio[reg.peps.na[[i]],]
  reg.phospeps[[i]] <- dat.abundance.gr[reg.peps.na[[i]],]
  reg.prot.info.used[[i]] <- prepareProteinInfo(pep.names=rownames(reg.phospeps[[i]]), prot.info=prot.dat.orig)
  write.table(x=reg.prot.info.used[[i]], file=paste0(cwd, "/reg.prot.info.used", i, ".txt"), sep="\t")
  #clueR.dat.tmp <- prepareClueRDat(sig.reg.pept=rownames(reg.phospeps[[i]]),
  #                                 dat.ratio.used=dat.ratio,
  #                                 mods.master.used=mods.master.used,
  #                                 psm.number.used=psm.number.used,
  #                                 organism="human",
  #                                 ratios=TRUE)
  clueR.dat.tmp <- prepareClueRDat(sig.reg.pept=rownames(reg.phospeps[[i]]),
                                   dat.ratio.used=dat.abundance.gr,
                                   mods.master.used=mods.master,
                                   psm.number.used=psm.number,
                                   organism="human",
                                   log=TRUE,
                                   ratios=FALSE)
  ksea.results[[i]] <- enrichmentTest(clust=rownames(clueR.dat.tmp),
                                      annotation=clueObj$annotation,
                                      universe=rownames(clueR.dat),
                                      alter = "greater")
                                      #alter="two.sided")
                                      #alter = "less")
  write.table(x=ksea.results[i][[1]], file=paste0(cwd, "/ClueR/ClueR-KSEA/ClueR-KSEA_timepoint_", i, ".txt"), row.names=FALSE, sep="\t")
}
ksea.results[5][[1]][1:10,]



top.kinases <- c()
for(i in 1:6){
  #top.kinases <- c(top.kinases, ksea.results[i][[1]][,1]) # --> All kinases!!! 
  #top.kinases <- c(top.kinases, ksea.results[i][[1]][1:10,1]) 
  top.kinases <- c(top.kinases, ksea.results[i][[1]][as.numeric(ksea.results[i][[1]][,"pvalue"]) < 0.15,1])
  #top.kinases <- c(top.kinases, ksea.results[i][[1]][as.numeric(ksea.results[i][[1]][,"# of substrates"]) > 10,1])
}
top.kinases <- unique(top.kinases)



kinase.heatmap.dat <- matrix(nrow=length(top.kinases), ncol=6)
rownames(kinase.heatmap.dat) <- top.kinases
colnames(kinase.heatmap.dat) <- c("1 min", "2.5 min", "5 min", "15 min", "30 min", "60 min")
for(i in 1:6){
  tmp <- ksea.results[i][[1]][,1:3]
  rownames(tmp) <- ksea.results[i][[1]][,1]
  #kinase.heatmap.dat[,i] <- as.numeric(ksea.results[i][[1]][ksea.results[i][[1]][,1] %in% top.kinases, "pvalue"])
  kinase.heatmap.dat[,i] <- as.numeric(tmp[top.kinases, "pvalue"])
}



manual.sort2 <- c(
"PLK1",
"MARK2",
"CDK6",
"PDK3",
"LATS1",
"ATM",
"IKBKE",
"OXSR1",
"MYLK",
"CIT",
"AURKB",
"PDPK1",
"CDK1",
"CDK5",
"GSK3A",
"MAP2K1",
"TBK1", 
"MAPKAPK2",
"MTOR",
"MAPK14",
"RPS6KB2", 
"MAPK3",
"PRKCB",
"SGK1",
"AKT1",
"MAPK1",
"RPS6KA1",
"PRKCA",
"RPS6KB1",
"PRKCD",
"MAPK11",
"SGK3",
"PAK1",
"PRKG2",
"GSK3B"
)

manual.sort2.lab <- c(
  expression(bold(paste("PLK1"^"\u{2731}"))), #"PLK1*",
  expression(bold(paste("MARK2"^"(\u{2731})"))), #"MARK2(*)",
  expression(bold(paste("CDK6"^"(\u{2731})"))), #CDK6(*)",
  "PDK3",
  "LATS1",
  "ATM",
  "IKBKE",
  expression(bold(paste("OXSR1"^"(\u{2731})"))), #"OXSR1(*)",
  "MYLK",
  "CIT",
  expression(bold(paste("AURKB"^"(\u{2731})"))), #"AURKB(*)",
  expression(bold(paste("PDPK1"^"\u{2731}\u{2731}"))), #"PDPK1**",
  expression(bold(paste("CDK1"^"(\u{2731})"))), #"CDK1(*)",
  expression(bold(paste("CDK5"^"\u{2731}"))), #"CDK5*",
  "GSK3A",
  "MAP2K1",
  expression(bold(paste("TBK1"^"(\u{2731})"))), #"TBK1(*)", 
  expression(bold(paste("MAPKAPK2"^"\u{2731}\u{2731}"))), #"MAPKAPK2**",
  expression(bold(paste("MTOR"^"(\u{2731})"))), #"MTOR(*)",
  expression(bold(paste("MAPK14"^"\u{2731}\u{2731}"))), #"MAPK14**",
  "RPS6KB2", 
  expression(bold(paste("MAPK3"^"\u{2731}"))), #"MAPK3*",
  expression(bold(paste("PRKCB"^"\u{2731}"))), #"PRKCB*",
  expression(bold(paste("SGK1"^"\u{2731}\u{2731}"))), #"SGK1**",
  expression(bold(paste("AKT1"^"\u{2731}\u{2731}\u{2731}"))), #"AKT1***",
  expression(bold(paste("MAPK1"^"\u{2731}"))), #"MAPK1*",
  expression(bold(paste("RPS6KA1"^"\u{2731}\u{2731}"))), #"RPS6KA1**",
  "PRKCA",
  expression(bold(paste("RPS6KB1"^"\u{2731}\u{2731}\u{2731}"))), #"RPS6KB1***",
  expression(bold(paste("PRKCD"^"\u{2731}"))), #"PRKCD*",
  "MAPK11",
  "SGK3",
  "PAK1",
  "PRKG2",
  "GSK3B"
) 


png(filename=paste0(cwd, "/ClueR/ClueR-KSEA/20230907_ClueR-KSEA.png"), height=2000, width=2000, res=300)
  par(font.axis=2)  
  heatmap.2(x=-log10(kinase.heatmap.dat[manual.sort2,]),
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          hclustfun = function(x) hclust(x, method="average"),
          Colv = FALSE,
          Rowv = FALSE,
          dendrogram="none",
          margins = c(7.5, 15), # for 2 arguments
          lwid=c(6,15), #make column of dendrogram and key very small and other colum very big 
          lhei=c(0.2,5), #make row of key and other dendrogram very small and other row big
          cexRow = 0.85,
          cexCol = 1.5,
          labRow = manual.sort2.lab,
          adjRow = c(0,0.3),
          adjCol = c(NA,0.5),
          trace="none",
          scale="row",
          col = bluered(100),
          na.color="white",
          key=FALSE,
          #colRow = gsub("grey", "grey25", pathway.clust$Color),
          #RowSideColors=pathway.clust$Color,
          sepwidth=c(0.001,0.001),
          sepcolor="grey",
          colsep=0:ncol(kinase.heatmap.dat),
          rowsep=0:nrow(kinase.heatmap.dat[manual.sort2,])
    )
dev.off()

#---

png(filename=paste0(cwd, "/ClueR/ClueR-KSEA/20230907_ClueR-KSEA_withKey.png"), height=2000, width=2000, res=300)
  par(font.axis=2)  
  heatmap.2(x=-log10(kinase.heatmap.dat[manual.sort2,]),
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          hclustfun = function(x) hclust(x, method="average"),
          Colv = FALSE,
          Rowv = FALSE,
          dendrogram="none",
          #margins = c(7.5, 15), # for 2 arguments
          lwid=c(1,0.5), #make column of dendrogram and key very small and other colum very big 
          lhei=c(1.1,5), #make row of key and other dendrogram very small and other row big
          cexRow = 0.85,
          cexCol = 1.5,
          labRow = manual.sort2.lab,
          adjRow = c(0,0.3),
          adjCol = c(NA,0.5),
          trace="none",
          scale="row",
          col = bluered(100),
          na.color="white",
          key=TRUE,
          keysize = 1.0,
          density.info="none",
          key.title = "Z-score of -log10(p-value)",
          key.xlab = NA,
          key.par = list("cex"=1.1),
          #colRow = gsub("grey", "grey25", pathway.clust$Color),
          #RowSideColors=pathway.clust$Color,
          sepwidth=c(0.001,0.001),
          sepcolor="grey",
          colsep=0:ncol(kinase.heatmap.dat),
          rowsep=0:nrow(kinase.heatmap.dat[manual.sort2,])
  )
dev.off()

#---

kinase.barplot.dat <- matrix(nrow=length(top.kinases), ncol=6)
rownames(kinase.barplot.dat) <- top.kinases
colnames(kinase.barplot.dat) <- c("1 min", "2.5 min", "5 min", "15 min", "30 min", "60 min")
for(i in 1:6){
  tmp <- ksea.results[i][[1]][,1:3]
  rownames(tmp) <- ksea.results[i][[1]][,1]
  kinase.barplot.dat[,i] <- as.numeric(tmp[top.kinases,3])  
}

write.table(
  x=top.kinases,
  file=paste0(cwd, "/ClueR/ClueR-KSEA/top_kinases.txt"),
  sep="\t",
  row.names=FALSE,
  col.names=FALSE
)

for(i in 1:length(top.kinases)){
  png(filename=paste0(cwd, "/ClueR/ks-barplots/ks-barplot_", top.kinases[i], ".png"), height=2000, width=2000, res=300)
    barplot.obj <- barplot(kinase.barplot.dat[top.kinases[i],], main=paste0("Number of ", top.kinases[i], " substrates with\nat least one regulated phosphosite"))
    text(x=barplot.obj, y=kinase.barplot.dat[top.kinases[i],], label=kinase.barplot.dat[top.kinases[i],], pos=1, cex=1, col="red")
  dev.off()
}

# #------
# # --> same as above for 'early', 'intermediate' & 'late'
# ksea.results.EIL <- kseaTP2EIL(ksea.results=ksea.results)
# 
# manual.sort3 <- c(
#   "RPS6KB2",
#   "LIMK1",
#   "LIMK2",
#   "TESK1",
#   "OXSR1",
#   "MARK3",
#   "ROCK2",
#   "MAP2K6",
#   "MAP3K5",
#   "PDPK1",
#   "PRKAA1", 
#   "CAMK2A",
#   "CDK4",
#   "PRKG1",
#   "PRKD1",
#   "SGK1",
#   "MTOR",
#   "MAP2K1",
#   "MAPKAPK2",
#   "AKT1",
#   "PRKCB",
#   "MAPK3",
#   "RPS6KA1",
#   "RPS6KA3",
#   "MAPK1",
#   "RPS6KB1",
#   "MAPK8",
#   "MAPK9",
#   "PRKCD",
#   "PIM1",
#   "CDK1",
#   "PRKACA"
# )
# 
# kinase.heatmap.dat.EIL <- matrix(nrow=length(manual.sort3), ncol=3)
# rownames(kinase.heatmap.dat.EIL) <- manual.sort3
# colnames(kinase.heatmap.dat.EIL) <- c("Early", "Intermediate", "Late")
# for(i in 1:3){
#   tmp <- ksea.results.EIL[i][[1]][,1:3]
#   rownames(tmp) <- ksea.results.EIL[i][[1]][,1]
#   kinase.heatmap.dat.EIL[,i] <- as.numeric(tmp[manual.sort3, "pvalue"])
# }
# 
# png(filename=paste0(cwd, "/ClueR/ClueR-KSEA/ClueR-KSEA-EIL.png"), height=2000, width=2000, res=300)
#     heatmap.2(x=-log10(kinase.heatmap.dat.EIL[manual.sort3,]),
#           distfun=function(x) as.dist((1-cor(t(x)))/2),
#           hclustfun = function(x) hclust(x, method="average"),
#           Colv = FALSE,
#           Rowv = FALSE,
#           dendrogram="none",
#           margins = c(10, 10),
#           lwid=c(1,1.25), #make column of dendrogram and key very small and other colum very big 
#           lhei=c(0.2,5), #make row of key and other dendrogram very small and other row big
#           cexRow = 1.1,
#           cexCol = 1.5,
#           adjCol = c(NA,0.5),
#           trace="none",
#           scale="row",
#           col = bluered(100),
#           na.color="white",
#           key=FALSE,
#           #colRow = gsub("grey", "grey25", pathway.clust$Color),
#           #RowSideColors=pathway.clust$Color,
#           sepwidth=c(0.001,0.001),
#           sepcolor="grey",
#           colsep=0:ncol(kinase.heatmap.dat.EIL),
#           rowsep=0:nrow(kinase.heatmap.dat.EIL[manual.sort3,])
#     )
# dev.off()
# #------

# library
library(tidyverse)

circos.kinases <- c(
      "AKT1",
      "MAPK1",
      "CDK1",
      "MAPK3",
      "RPS6KA1",
      "RPS6KB1",
      "MAPK14",
      "MTOR",
      "MAPKAPK2",
      "SGK1")
circos.dat <- data.frame(
  individual=as.vector(t(kinase.barplot.dat[circos.kinases,]), mode="character"),
  group=rep(x=rownames(kinase.barplot.dat[circos.kinases,]), each=ncol(kinase.barplot.dat)),
  value=c(t(kinase.barplot.dat[circos.kinases,]))
)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
for (i in 1:length(circos.dat$group)) {
  circos.dat$group[i] <- gsub("(.*)", paste0(formatC(ceiling(i/ncol(kinase.barplot.dat)), width=2, format ="d", flag ="0"),"_\\1"), circos.dat$group[i])
}
to_add <- data.frame( matrix(NA, empty_bar*nlevels(as.factor(circos.dat$group)), ncol(circos.dat)) )
colnames(to_add) <- colnames(circos.dat)
to_add$group <- rep(levels(as.factor(circos.dat$group)), each=empty_bar)
circos.dat <- rbind(circos.dat, to_add)
circos.dat <- circos.dat %>% arrange(group)
circos.dat$id <- seq(1, nrow(circos.dat))
#for (i in 1:length(circos.dat$group)) {
#  circos.dat$group[i] <- gsub("^\\d+_", "", circos.dat$group[i])
#}

# Get the name and the y position of each label
label_data <- circos.dat
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- circos.dat %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

#
circos.dat$group <- gsub("^\\d+_", "", circos.dat$group)
label_data$group <- gsub("^\\d+_", "", label_data$group)
base_data$group <- gsub("^\\d+_", "", base_data$group)
grid_data$group <- gsub("^\\d+_", "", grid_data$group)

# Make the plot
p <- ggplot(circos.dat, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 6, xend = start, yend = 6), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 12, xend = start, yend = 12), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 18, xend = start, yend = 18), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 24, xend = start, yend = 24), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(circos.dat$id),4), y = c(6, 12, 18, 24), label = c("6", "12", "18", "24") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-30,40) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+2, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) + 

  # Add base line information
  geom_segment(data=base_data, mapping=aes(x=start, y=-1, xend=end, yend=-1), colour="black", alpha=1.0, size=0.6 , inherit.aes=FALSE ) +  
  geom_text(data=base_data,
            aes(x=title, y=-4, label=group),
            hjust=c(0.5,0.75,0.8,0.85,0.45,0.5,0.2,0.25,0.15,0.3),
            vjust=c(0.5,0.5,0.5,0.5,-0.15,0.25,0.5,0.5,0.5,0.5),
            colour="black",
            alpha=0.8,
            size=2.5,
            fontface="bold",
            inherit.aes = FALSE)

p

png(filename=paste0(cwd, "/ClueR/ClueR-KSEA/circular_KS-barplots.png"), height=3200, width=3200, res=600)
  p
dev.off()

#------

#-----------------------------------------------

# #devtools::install_github("PengyiYang/ClueR")
# library(ClueR)
# library(PhosR)
# #data(PhosphoSite) #PhosphoSitePlus-ClueR-Version: 1) overlap phos-sites with MS-data: 233 2) kinases: 206 3) substrates: 9830
# data(PhosphoSitePlus) #PhosphoSitePlus-PhosR-Version: 1) overlap phos-sites with MS-data: 275 2) kinases: 379 3) substrates: 11998
# clue.pepforms <- TRUE
# 
# sig.reg.pept <- unique(c(Early, Intermediate, Late))
# dat.ratio.used.sig <- dat.ratio.used[sig.reg.pept,]
# mods.master.used.sig <- mods.master.used[sig.reg.pept]
# psm.number.used.sig <- psm.number.used[sig.reg.pept]
# 
# clueR.dat.tmp <- matrix(ncol=ncol(dat.ratio.used.sig))
# psm.number.used.sig.clueR.tmp <- c()
# #test.idx <- c(85:100, 205:215)
# #test.idx <- 1100:1200
# #for(i in test.idx){
# for(i in 1:length(mods.master.used.sig)){
#   if(clue.pepforms == TRUE && (length(grep("]; ", mods.master.used.sig[i])) > 0)){
#       splitted.mods <- unlist(strsplit(mods.master.used.sig[i], "]; "))
#       splitted.mods[1:(length(splitted.mods)-1)] <- paste0(splitted.mods[1:(length(splitted.mods)-1)], "]")
#       
#       last.uniprot.id <- ""
#       for(j in 1:length(splitted.mods)){
#           if(length(grep("^\\d+x", splitted.mods[j]) > 0)){
#             splitted.mods[j] <- paste0(last.uniprot.id, " ", splitted.mods[j])
#           }else{
#             last.uniprot.id <- gsub("^([0-9A-Z]+) \\d+x.+", "\\1", splitted.mods[j])
#           }  
#       }
#       
#       for(j in 1:length(splitted.mods)){
#         clueR.dat.tmp <- rbind(clueR.dat.tmp, dat.ratio.used.sig[i,]) 
#         rownames(clueR.dat.tmp)[nrow(clueR.dat.tmp)] <- splitted.mods[j]
#         psm.number.used.sig.clueR.tmp <- c(psm.number.used.sig.clueR.tmp, psm.number.used.sig[i])
#         names(psm.number.used.sig.clueR.tmp)[length(psm.number.used.sig.clueR.tmp)] <- splitted.mods[j]
#       }
#   }else if(length(grep("]; ", mods.master.used.sig[i])) == 0){
#     clueR.dat.tmp <- rbind(clueR.dat.tmp, dat.ratio.used.sig[i,]) 
#     rownames(clueR.dat.tmp)[nrow(clueR.dat.tmp)] <- mods.master.used.sig[i]
#     psm.number.used.sig.clueR.tmp <- c(psm.number.used.sig.clueR.tmp, psm.number.used.sig[i])
#     names(psm.number.used.sig.clueR.tmp)[length(psm.number.used.sig.clueR.tmp)] <- mods.master.used.sig[i]
#   }
# }
# if(sum(is.na(clueR.dat.tmp[1,])) == ncol(clueR.dat.tmp)) clueR.dat.tmp <- clueR.dat.tmp[-1,] # <---------- ????
# 
# clueR.dat <- matrix(ncol=ncol(clueR.dat.tmp))
# psm.number.used.sig.clueR <- c()
# for(i in 1:nrow(clueR.dat.tmp)){
#   print(paste0("original: ", rownames(clueR.dat.tmp)[i]))
#   if(clue.pepforms == TRUE && (length(grep("); ", rownames(clueR.dat.tmp)[i])) > 0)){
#     uniprot.id <- gsub("^([0-9A-Z]+)-*\\d* .*", "\\1", rownames(clueR.dat.tmp)[i], perl=TRUE)
#     tmp <- gsub("^.* \\d+xPhospho ", "", rownames(clueR.dat.tmp)[i], perl=TRUE)
#     tmp <- gsub("\\(\\d+\\.*\\d*\\)", "", tmp)
#     #tmp <- gsub("(\\[|\\]|[A-Z])", "", tmp)
#     tmp <- gsub("(\\[|\\])", "", tmp)
#     splitted.sites <- unlist(strsplit(tmp, "; "))
#     for(j in 1:length(splitted.sites)){
#         print(paste0(uniprot.id, ":", splitted.sites[j]))
#         clueR.dat <- rbind(clueR.dat, clueR.dat.tmp[i,])
#         rownames(clueR.dat)[nrow(clueR.dat)] <- paste0(uniprot.id, ":", splitted.sites[j]) 
#         print(paste0("PSM number: ", psm.number.used.sig.clueR.tmp[i]))
#         psm.number.used.sig.clueR <- c(psm.number.used.sig.clueR, psm.number.used.sig.clueR.tmp[i])
#         names(psm.number.used.sig.clueR)[length(psm.number.used.sig.clueR)] <- paste0(uniprot.id, ":", splitted.sites[j]) 
#     }
#   }else if(is.na(rownames(clueR.dat.tmp)[i])){
#     print(paste0("NA rowname found for row: ", i)) 
#     print(rownames(clueR.dat.tmp)[i])
#     print(paste0("PSM number: ", psm.number.used.sig.clueR.tmp[i]))
#   }else if(length(grep("); ", rownames(clueR.dat.tmp)[i])) == 0){
#     uniprot.id <- gsub("^([0-9A-Z]+)-*\\d* .*", "\\1", rownames(clueR.dat.tmp)[i], perl=TRUE)
#     #site <- gsub(".*\\[[A-Z]{1}(\\d+)\\(\\d+\\.*\\d*\\)\\]$", "\\1", rownames(clueR.dat.tmp)[i], perl=TRUE)
#     site <- gsub(".*\\[([A-Z]{1}\\d+)\\(\\d+\\.*\\d*\\)\\]$", "\\1", rownames(clueR.dat.tmp)[i], perl=TRUE)
#     print(paste0(uniprot.id, ":", site))
#     clueR.dat <- rbind(clueR.dat, clueR.dat.tmp[i,])
#     rownames(clueR.dat)[nrow(clueR.dat)] <- paste0(uniprot.id, ":", site)
#     print(paste0("PSM number: ", psm.number.used.sig.clueR.tmp[i]))
#     psm.number.used.sig.clueR <- c(psm.number.used.sig.clueR, psm.number.used.sig.clueR.tmp[i])
#     names(psm.number.used.sig.clueR)[length(psm.number.used.sig.clueR)] <- paste0(uniprot.id, ":", site) 
#   }
#   print("------------------------------------------------------")
# }
# if(sum(is.na(clueR.dat[1,])) == ncol(clueR.dat)) clueR.dat <- clueR.dat[-1,]
# if(length(grep("\\[",rownames(clueR.dat))) > 0) clueR.dat <- clueR.dat[-grep("\\[",rownames(clueR.dat)), ]
# if(length(grep("/",rownames(clueR.dat))) > 0) clueR.dat <- clueR.dat[-grep("/",rownames(clueR.dat)), ]
# clueR.dat <- clueR.dat[-which(is.na(clueR.dat[,1])), ]
# 
# 
# del.idx <- c()
# #---
# unique.rows <- unique(rownames(clueR.dat))
# for(i in 1:length(unique.rows)){
#   if(length(grep(unique.rows[i], rownames(clueR.dat))) > 1){
#     tmp.idx.vec <- grep(unique.rows[i], rownames(clueR.dat))
#     print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#     print("tmp.idx.vec:")
#     print(tmp.idx.vec)
#     print("PSM-Numbers:")
#     print(psm.number.used.sig.clueR[tmp.idx.vec])
#     keep.idx <- which.max(psm.number.used.sig.clueR[tmp.idx.vec])
#     print("keep.idx:")
#     print(keep.idx)
#     tmp.idx.vec <- tmp.idx.vec[-keep.idx]
#     print("tmp.idx.vec:")
#     print(tmp.idx.vec)
#     del.idx <- c(del.idx, tmp.idx.vec)
#   }
# }
# #---
# #for(i in 2:(nrow(clueR.dat)-1)){
# #    if(length(del.idx) > 0) {
# #      cond1 <- length(which(rownames(clueR.dat[-del.idx,]) == rownames(clueR.dat)[i])) > 1
# #    }else{
# #      cond1 <- length(which(rownames(clueR.dat) == rownames(clueR.dat)[i])) > 1  
# #    }
# #    cond2 <- clueR.dat[i,] == clueR.dat[i+1,] 
# #    cond3 <- clueR.dat[i,] == clueR.dat[i-1,] 
# #    if(cond1 && cond2){
# #      del.idx <- c(del.idx,i+1)
# #    }else if(cond1 && cond3){
# #      del.idx <- c(del.idx,i)  
# #    }
# #}
# #---
# clueR.dat <- clueR.dat[-del.idx,]
# clueR.dat <- clueR.dat[!duplicated(clueR.dat),]
# 
# 
# 
# #id.vec <- gsub("\\:\\d+$", "", rownames(clueR.dat))
# id.vec <- gsub("\\:[A-Z]{1}\\d+$", "", rownames(clueR.dat))
# site.vec <- gsub("^[0-9A-Z]+\\:", "", rownames(clueR.dat))
# id.mapping <- returnUniProtGeneNames(uniprot.ids=id.vec, table=TRUE)
# for(i in 1:nrow(id.mapping)){
#   id.vec[which(id.vec == id.mapping[i,"From"])] <- id.mapping[i,"To"] 
# }
# rownames(clueR.dat) <- paste0(id.vec,";",site.vec,";")
# 
# clueR.dat <- clueR.dat[!duplicated(rownames(clueR.dat)),]
# 
# baseline.col <- rep(1, nrow(clueR.dat))
# clueR.dat <- cbind("Abundance Ratio: (0) / (0)"=baseline.col, clueR.dat)
# clueR.dat <- log2(clueR.dat)
# colnames(clueR.dat) <- c("0min", "1min", "2.5min", "5min", "15min", "30min", "60min")
# 
# write.table(cbind(rownames(clueR.dat),clueR.dat), file=paste0(cwd, "/ClueR/clueR.dat2.txt"), sep="\t", col.names = TRUE, row.names = FALSE)
# save(clueR.dat, file=paste0(cwd, "/ClueR/clueR.dat2.RData"), compress = "xz", compression_level = 9)


#-----------------------------------------------

#m <- as.matrix(dat.orig[,c(18:23)]) # 13196 rows
#rownames(m) <- dat.orig$`Master Protein Accessions`
####################################
#m <- m[-(na.idx_0NA_oneNA_60NA),]
####################################
#m <- m[rowSums(is.na(m)) != ncol(m), ] #12764 phosphopeptides
#print(paste0("Number of phosphopeptides: ", nrow(m)))
#h <- rownames(m)
#h <- strsplit(h, split = ";", fixed = T)
#h <- unlist (lapply (1:length (h), function (x) h[[x]][1]))
#h <- unlist (lapply (1:length (h), function (x) strsplit(h, split = "-", fixed = T)[[x]][1]))
#y <- lapply (1:length (h), function (x) UniProtToGeneNames$To[UniProtToGeneNames$From == h[x] ])

#y[lengths(y) == 0] <- NaN 
#length (unique (unlist (y)))
#print(paste0("Number of proteins: ", length (unique (unlist (y)))))

#-----------------------------------------------

#-------------------------------------------------------------------------------
#Comparison total number NA's
#tmp.dat.ratio <- dat.ratio
#tmp.pd.duplicate.ratios <- pd.duplicate.ratios
#tmp.dat.own.ratio <- dat.own.ratio

#tmp1 <- sum(c(tmp.dat.ratio) == 0.01, na.rm=TRUE)
#tmp2 <- sum(c(tmp.dat.ratio) == 100, na.rm=TRUE)
#tmp3 <- sum(is.na(tmp.dat.ratio))
#tmp4 <- length(c(dat.ratio)) - (tmp1+tmp2+tmp3) 

#tmp5 <- sum(c(tmp.pd.duplicate.ratios) == 0.01, na.rm=TRUE)
#tmp6 <- sum(c(tmp.pd.duplicate.ratios) == 100, na.rm=TRUE)
#tmp7 <- sum(is.na(tmp.pd.duplicate.ratios))
#tmp8 <- length(c(pd.duplicate.ratios)) - (tmp5+tmp6+tmp7)

#tmp9 <- sum(c(tmp.dat.own.ratio) == 0.01, na.rm=TRUE)
#tmp10 <- sum(c(tmp.dat.own.ratio) == 100, na.rm=TRUE) 
#tmp11 <- sum(is.na(tmp.dat.own.ratio))
#tmp12 <- length(c(dat.own.ratio)) - (tmp9+tmp10+tmp11)

#tmp13 <- cbind(c(tmp1,tmp2,tmp3,tmp4), c(tmp5,tmp6,tmp7,tmp8), c(tmp9,tmp10,tmp11,tmp12))
#colnames(tmp13) <- c("PD ratios", "PD ratios\n(self-calc)", "Abundance ratios\n(self-calc)")

#tmp.sum1 <- sum(c(tmp1,tmp2,tmp3))
#tmp.sum2 <- sum(c(tmp5,tmp6,tmp7))
#tmp.sum3 <- sum(c(tmp9,tmp10,tmp11))
##sum(c(tmp9,tmp10,tmp11)) / length(c(dat.ratio))*100

#png(filename=paste0(cwd, "/comparison_ratios1.png "), height=2000, width=2000, res=300)
#barplotobj <- barplot(tmp13, main="Comparison of ratios", las=1)
#text(x=barplot.obj, y=c(length(c(dat.ratio)), length(c(dat.ratio)), length(c(dat.ratio))), label=c(length(c(dat.ratio)), length(c(dat.ratio)), length(c(dat.ratio))), pos=1, col="black")
#dev.off()
##barplot(tmp13, ylim=c(0,10000), las=2)
#png(filename=paste0(cwd, "/comparison_ratios2.png "), height=2000, width=2000, res=300)
#barplot.obj <- barplot(tmp13, ylim=c(0,15000), main="Comparison of ratios\n(red: 0.01s, blue: 100s, dark green: NAs)", las=1)
#text(x=barplot.obj, y=c(tmp.sum1, tmp.sum2, tmp.sum3), label=c(tmp.sum1, tmp.sum2, tmp.sum3), pos=3, col="black")
#text(x=barplot.obj, y=c(tmp1, tmp5, tmp9), label=c(tmp1, tmp5, tmp9), pos=1, col="red")
#text(x=barplot.obj, y=c(tmp1+tmp2, tmp5+tmp6, tmp9+tmp10), label=c(tmp2, tmp6, tmp10), pos=1, col="blue")
#text(x=barplot.obj, y=c(tmp1+tmp2+tmp3, tmp5+tmp6+tmp7, tmp9+tmp10+tmp11), label=c(tmp3, tmp7, tmp11), pos=1, col="darkgreen")
#dev.off()

#-------------------------------------------------------------------------------

#output <- cbind(p.values[,1], p.values.own[,1], dat.ratio[,1], abs(p.values.own[,1] - dat.ratio[,1]))
#write.table(x=output, file=paste0(cwd, "/p-val_compare.txt"), sep="\t")
