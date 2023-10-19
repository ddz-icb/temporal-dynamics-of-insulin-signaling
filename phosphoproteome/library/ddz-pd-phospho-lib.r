#++++++++++++++++++++++ assessReproducibility(...) +++++++++++++++++++++++++++++

assessReproducibility <- function(x=NULL, sample.ids=NULL, output.path=NULL){
  
  if(!is.null(output.path)) {
    output.path.tmp <- paste0(output.path, "/Reproducibility")
  }
  if(!is.null(output.path) && !dir.exists(output.path.tmp)) {
    dir.create(output.path.tmp)  
  }
  
  cor.vec <- vector(mode="numeric", length=length(sample.ids))
  cv.vec <- vector(mode="numeric", length=length(sample.ids))
  for(i in 1:length(sample.ids)){
    print(paste0("i: ", i))
    replicates.id.vec <- grep(sample.ids[i], colnames(x))
    cor.vec.vec <- vector(mode="numeric", length=choose(length(replicates.id.vec),2))
    cv.vec.vec <- vector(mode="numeric", length=choose(length(replicates.id.vec),2))
    combi <- 1
    for(j in 1:(length(replicates.id.vec)-1)){
      for(k in (j+1):length(replicates.id.vec)){
        replicate1 <- x[,replicates.id.vec[j]]
        replicate2 <- x[,replicates.id.vec[k]]
        
        if(any(is.na(replicate1)) || any(is.na(replicate2))){
          na.ids1 <- which(is.na(replicate1))
          na.ids2 <- which(is.na(replicate2))
          na.ids <- unique(c(na.ids1, na.ids2))
          
          replicate1 <- replicate1[-na.ids]
          replicate2 <- replicate2[-na.ids]
        }
        
        cor.vec.vec[combi] <- cor(log2(replicate1), log2(replicate2))
        #cor.vec.vec[combi] <- cor(replicate1, replicate2)
        print(paste0("cor.vec.vec[", combi, "]: ", cor.vec.vec[combi]))
        
        mean.vec.tmp <- apply(cbind(replicate1, replicate2), 1, mean, na.rm=T)
        sd.vec.tmp <- apply(cbind(replicate1, replicate2), 1, sd, na.rm=T)
        cv.vec.tmp <- sd.vec.tmp/mean.vec.tmp*100
        
        cv.vec.vec[combi] <- median(cv.vec.tmp)
        print(paste0("cv.vec.vec[", combi, "]: ", cv.vec.vec[combi]))
        
        combi <- combi+1
      }
    }
    cor.vec[i] <- median(cor.vec.vec)
    cv.vec[i] <- median(cv.vec.vec)
  }
  
  print("--- FINAL RESULTS ---")
  print("cor.vec:")
  print(cor.vec)
  print("cv.vec:")
  print(cv.vec)
  
  output1 <- cbind(sample.ids,cor.vec)
  colnames(output1) <- c("Sample ID","Correlation of techn. replicates")
  print(output1)
  write.table(x=output1,
              file=paste0(output.path, "/Reproducibility/tech_repl_cor.txt"),
              row.names=F,
              col.names=T,
              sep="\t")
  
  output2 <- cbind(sample.ids,cv.vec)
  colnames(output2) <- c("Sample ID","CVs of techn. replicates")
  print(output2)
  write.table(x=output2,
              file=paste0(output.path, "/Reproducibility/tech_repl_cv.txt"),
              row.names=F,
              col.names=T,
              sep="\t")
  
  print(median(cor.vec))
  print(median(cv.vec))
  
  
  png(paste(output.path, "/Reproducibility/tech_repl_cv_boxplot.png", sep=""), width=3600, height=3600, pointsize=15, res=600)
    boxplot(cv.vec,
          border="white",
          main=paste0("Average CVs of phosphopeptides abundances \nbetween technical replicates (median: ",
                      round(median(cv.vec),2),
                      "%)")
          )
    boxplot(cv.vec, outline=F, add=T)
    stripchart(cv.vec,
             method="jitter",
             jitter=0.05,
             pch=19,
             cex=2.0,
             col=alpha("red",alpha=0.4),
             vertical=T,
             add=T
             )
  dev.off()
  
  
  png(paste(output.path, "/Reproducibility/tech_repl_cor_boxplot.png", sep=""), width=3600, height=3600, pointsize=15, res=600)
    boxplot(cor.vec,
          border="white",
          main=paste0("Pearson's r between all technical replicates\n(median: ",
                      round(median(cor.vec),2),
                      ")")
          )
    boxplot(cor.vec, outline=F, add=T)
    stripchart(cor.vec,
             method="jitter",
             jitter=0.05,
             pch=19,
             cex=2.0,
             col=alpha("red",alpha=0.4),
             vertical=T,
             add=T
    )
  dev.off()
}

#++++++++++++++++++++++ assessReproducibility(...) +++++++++++++++++++++++++++++

#++++++++++++++++++++++ assessReproducibility2(...) +++++++++++++++++++++++++++++

assessReproducibility2 <- function(x=NULL,
                                   sample.ids=NULL,
                                   cex.lab=NULL,
                                   modify.ids=TRUE,
                                   anonymize.ids=TRUE,
                                   mode="cor",
                                   output.path=NULL){
  
  library(wesanderson)
  
  if(!is.null(output.path)) {
    output.path.tmp <- paste0(output.path, "/Reproducibility")
  }
  if(!is.null(output.path) && !dir.exists(output.path.tmp)) {
    dir.create(output.path.tmp)  
  }
  
  #pal <- colorpanel(n=100,low="white", high="red")
  pal <- colorpanel(n=100,low="blue", mid="white", high="red")
  #pal <- as.vector(wes_palette(name="Zissou1", n=20, type="continuous"))
  #pal <- wes_palette(name="BottleRocket1", n=6, type="continuous")
  #pal <- wes_palette(name="BottleRocket2", n=6, type="continuous")
  #pal <- wes_palette(name="Rushmore", n=6, type="continuous")
  #pal <- wes_palette(name="Royal1", n=6, type="continuous")
  #pal <- wes_palette(name="Royal2", n=6, type="continuous")
  #pal <- wes_palette(name="Darjeeling1", n=6, type="continuous")
  #pal <- wes_palette(name="Darjeeling2", n=6, type="continuous")
  #pal <- wes_palette(name="Chevalier1", n=6, type="continuous")
  #pal <- wes_palette(name="FantasticFox1", n=6, type="continuous")
  #pal <- wes_palette(name="GrandBudapest1", n=6, type="continuous")
  #pal <- wes_palette(name="GrandBudapest2", n=6, type="continuous")
  #pal <- wes_palette(name="Zissou1", n=6, type="continuous")
  #pal <- wes_palette(name="Cavalcanti1", n=6, type="continuous")
  #pal <- wes_palette(name="Moonrise1", n=6, type="continuous")
  #pal <- wes_palette(name="Moonrise2", n=6, type="continuous")
  #pal <- wes_palette(name="Moonrise3", n=6, type="continuous")
  
  output <- matrix(nrow=length(sample.ids), ncol=length(sample.ids))
  if(modify.ids==TRUE){
    name.vec1 <- gsub("^.* Sample, (\\d+\\.*\\d*), [A-Z]+\\d+$","\\1",sample.ids) # extract time point vector
    print(name.vec1)
    name.vec2 <- gsub("^.*, ([A-Z]+\\d+)$","\\1",sample.ids) # extract sample id vector
    if(anonymize.ids==TRUE){
      
      # ananoymization of human time course data:
      # "AR7" ---> "A"
      # "AT13" ---> "B"
      # "AT18" ---> "C"
      # "AT22" ---> "D"
      # "AT8" ---> "E"
      
      id.tmp1 <- sort(unique(name.vec2))
      id.tmp2 <- name.vec2
      if(length(id.tmp1) <= 26){
        new.ids <- LETTERS
      }else{
        new.ids <- paste0("S", 1:length(id.tmp1))
      }
      for(i in 1:length(id.tmp1)){
        id.tmp2 <- gsub(id.tmp1[i], new.ids[i], id.tmp2)  
      }
      name.vec2 <- id.tmp2
    }
    print(name.vec2)
    name.vec3 <- rep(1:2, length(sample.ids)/2) # create replicate id vector
    print(name.vec3)
    name.vec <- vector(mode="character", length=length(sample.ids))
    for(i in 1:length(sample.ids)) name.vec[i] <- paste0(name.vec2[i],
                                                       "-",
                                                       name.vec1[i],
                                                       "-",
                                                       name.vec3[i],
                                                       collapse="-")
    print(cbind(name.vec2, name.vec1, name.vec3, sample.ids))
    #print(name.vec)
    #name.vec <- gsub("^Abundance: ", "", sample.ids)
    #name.vec <- gsub("Sample, \\d+\\.*\\d*, ", "", name.vec)
    colnames(output) <- name.vec
    rownames(output) <- name.vec
  }else{
    colnames(output) <- sample.ids
    rownames(output) <- sample.ids  
  }
  for(i in 1:length(sample.ids)){
    print(paste0("i: ", i))
    for(j in 1:length(sample.ids)){
        #print(paste0("j: ", j))
        replicate1 <- x[,i]
        replicate2 <- x[,j]
        
        if(any(is.na(replicate1)) || any(is.na(replicate2))){
          na.ids1 <- which(is.na(replicate1))
          na.ids2 <- which(is.na(replicate2))
          na.ids <- unique(c(na.ids1, na.ids2))
          
          replicate1 <- replicate1[-na.ids]
          replicate2 <- replicate2[-na.ids]
        }
        
        if(mode=="cor"){
          #output[i,j] <- cor(log2(replicate1), log2(replicate2))
          output[i,j] <- cor(replicate1, replicate2)
          #print(paste0("cor.vec.vec[", combi, "]: ", cor.vec.vec[combi]))
        }else if(mode=="cv"){
          mean.vec.tmp <- apply(cbind(replicate1, replicate2), 1, mean, na.rm=T)
          sd.vec.tmp <- apply(cbind(replicate1, replicate2), 1, sd, na.rm=T)
          cv.vec.tmp <- sd.vec.tmp/mean.vec.tmp*100
          output[i,j] <- median(cv.vec.tmp)
          #cv.vec.vec[combi] <- median(cv.vec.tmp)
          #print(paste0("cv.vec.vec[", combi, "]: ", cv.vec.vec[combi]))
        }
    }
  }
  print(output)
  
  #if(mode=="cor"){
  #  heatmap.breaks <- 0:100/100
  #}else if(mode=="cv"){
  #  heatmap.breaks <- 0:100
  #}
  
  heatmap.2(output, trace="none", col=pal)
  #heatmap.2(output, trace="none", breaks=heatmap.breaks, col=pal)
  
  png(paste0(output.path, "/Reproducibility/tech_repl_", mode, "_heatmap.png"), width=3600, height=3600, pointsize=15, res=300)
    heatmap.2(output,
              trace="none",
              #breaks=heatmap.breaks,
              col=pal,
              cexRow=cex.lab,
              cexCol=cex.lab,
              adjRow=c(0,NA),
              adjCol=c(NA,0.5),
              density.info="none")
  dev.off()
  
  write.table(x=cbind(rownames(output),output),
              file=paste0(output.path, "/Reproducibility/tech_repl_", mode, "_heatmap.txt"),
              col.names=T,
              row.names=F,
              sep="\t")
  
  return(output)
}

#++++++++++++++++++++++ assessReproducibility2(...) +++++++++++++++++++++++++++++

#++++++++++++++++++++++ giveContaminPeps(...) ++++++++++++++++++++++++++++++++++++

giveContaminPeps <- function(pep.dat=NULL){
  
  output <- grep("TRUE", pep.dat[,"Contaminant"])
  return(output)
  
}

#++++++++++++++++++++++ giveContaminPeps(...) ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ giveInformativeMissedCleav(...) ++++++++++++++++++++++++++++++++++++

giveInformativeMissedCleav <- function(missed.cleav=NULL, mods.master=NULL){
  
  library(stringr)
  
  missed.uniprot.idx <- gsub("_peptide\\d+$", "", names(mods.master)[missed.cleav])
  output <- c()
  
  for(i in 1:length(missed.cleav)){
    mods.curr <- mods.master[missed.cleav[i]]
    mods.all <- mods.master[grep(missed.uniprot.idx[i], names(mods.master), value=TRUE)]
    mods.all <- mods.all[names(mods.all) != names(missed.cleav[i])]
    
    #print(mods.curr)
    mods.curr2 <- stringr::str_extract_all(mods.curr, "\\[[A-Z0-9\\(\\)\\.\\;\\/\\s]+\\]", simplify=TRUE)
    #if(length(mods.curr2)==0) mods.curr2 <- matrix(data=NA, nrow=1, ncol=1)
    mods.curr3 <- gsub("(\\[|\\]|\\(\\d+\\.*\\d*\\))", "", mods.curr2[1,])
    mods.curr4 <- c()
    for(j in 1:length(mods.curr3)){
      if(length(grep("\\; ", mods.curr3[j])) > 0){
        mods.curr4 <- c(mods.curr4, unique(unlist(strsplit(mods.curr3[j], "\\; "))))
      }else{
        mods.curr4 <- c(mods.curr4, mods.curr3[j])
      }
    }
    mods.curr4 <- sort(unique(mods.curr4))
    mods.curr4 <- mods.curr4[mods.curr4 != ""]
    #print(mods.curr4)
    
    #print("++++++++++++++++++++++++++++")
    
    #print(mods.all)
    if(length(mods.all) > 0){
      mods.all2 <- stringr::str_extract_all(mods.all, "\\[[A-Z0-9\\(\\)\\.\\;\\/\\s]+\\]", simplify=TRUE)
      mods.all3 <- c()
      for(j in 1:nrow(mods.all2)) mods.all3 <- c(mods.all3, gsub("(\\[|\\]|\\(\\d+\\.*\\d*\\))", "", mods.all2[j,]))
      mods.all4 <- c()
      for(j in 1:length(mods.all3)){
        if(length(grep("\\; ", mods.all3[j])) > 0){
          mods.all4 <- c(mods.all4, unique(unlist(strsplit(mods.all3[j], "\\; "))))
        }else{
          mods.all4 <- c(mods.all4, mods.all3[j])
        }
      }
      mods.all4 <- sort(unique(mods.all4))
      mods.all4 <- mods.all4[mods.all4 != ""]
      #print(mods.all4)
      
      
      if(length(setdiff(mods.curr4, mods.all4)) > 0){
        #print("INFORMATIVE:")
        #print(setdiff(mods.curr4, mods.all4))
        output <- c(output, missed.cleav[i])
      }
      
    }else{
      output <- c(output, missed.cleav[i])
    }
    #print("-------------------------")
    #print("")
    #print("")
  }
  return(output)
}

#++++++++++++++++++++++ giveInformativeMissedCleav(...) ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ givePeptidoforms(...) ++++++++++++++++++++++++++++++++++++

givePeptidoforms <- function(seq=NULL){
  
  # Note: seq must be sorted, i.e. duplicates must occur consecutively
  multi.idx <- duplicated(seq)  #gives only duplicates without first occurrence
  for(i in 1:(length(multi.idx)-1)){
    if(multi.idx[i+1] == TRUE) multi.idx[i] <- TRUE 
  }
  return(multi.idx)
}

#++++++++++++++++++++++ givePeptidoforms(...) ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ giveProtPhosStatus(...) ++++++++++++++++++++++++++++++++++++

giveProtPhosStatus <- function(x){
  
  library(stringr)
  
  uniprot.ids <- gsub("_peptide\\d+$", "", rownames(x))
  #uniprot.ids <- gsub("\\-\\d+$", "", x=uniprot.ids)
  
  output <- list()
  for(i in 1:nrow(x)){
    #print(x[i, "Modifications in Master Proteins"])  
    
    mods <- stringr::str_extract_all(x[i, "Modifications in Master Proteins"], "\\[[A-Z0-9\\(\\)\\.\\;\\/\\s]+\\]", simplify=TRUE)
    #print(mods)
    mods2 <- gsub("(\\[|\\]|\\(\\d+\\.*\\d*\\))", "", mods[1,])
    #print(mods2)
    mods3 <- c()
    for(j in 1:length(mods2)){
      if(length(grep("\\; ", mods2[j])) > 0){
        mods3 <- c(mods3, unique(unlist(strsplit(mods2[j], "\\; "))))
      }else{
        mods3 <- c(mods3, mods2[j])
      }
    }
    mods3 <- sort(unique(mods3))
    mods3 <- mods3[mods3 != ""]
    #print(mods3)
    
    #print(uniprot.ids[i])
    if(length(grep("\\; ", uniprot.ids[i])) > 0){
      #print("giveProtPhosStatus2: HERE!!!")
      splitted.accessions <- unique(unlist(strsplit(grep("\\; ", uniprot.ids[i], value=TRUE), "\\; ")))
      for(j in 1:length(splitted.accessions)){
        if(exists(splitted.accessions[j], where = output)){
          output[[splitted.accessions[j]]] <- sort(unique(c(output[[splitted.accessions[j]]], mods3)))
        }else{
          output[[splitted.accessions[j]]] <- mods3
        }
        #print(splitted.accessions[j])
        #print(output[[splitted.accessions[j]]])
      }
    }else{
      if(exists(uniprot.ids[i], where = output)){
        output[[uniprot.ids[i]]] <- sort(unique(c(output[[uniprot.ids[i]]], mods3)))
      }else{
        output[[uniprot.ids[i]]] <- mods3
      }
      #print(uniprot.ids[i])
      #print(output[[uniprot.ids[i]]])
    }
    
    #print("-------------------------------")
    #print("")
    #print("")
  }
  
  return(output)
  
}

#++++++++++++++++++++++ giveProtPhosStatus(...) ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ giveTermClassi(...) ++++++++++++++++++++++++++++++++++++

giveTermClassi <- function(file.path=NULL, term.filter=NULL, output.path=NULL){
  x <- read.table(file=file.path, header=F, skip=1, sep="\t", quote = "\"")
  colnames(x) <- c("Gene ID",
                   "Mapped IDs",
                   "Gene Name/Gene Symbol/Persistent ID/Orthologs",
                   "PANTHER Family/Subfamily",
                   "PANTHER Protein Class",
                   "Species",
                   "Terms")
  
  terms <- c()
  for(i in 1:nrow(x)) terms <- c(terms, unlist(strsplit(x[i,"Terms"],"\\;")))
  print(paste0("number of annotations: ", length(terms)))
  
  unique.terms <- unique(grep("GO:\\d{7}", terms, value=T))
  unique.term.ids <- unique(gsub("^.+\\((GO:\\d{7})\\)$", "\\1", terms))
  
  print(paste0("number of unique terms: ", length(unique.terms)))
  
  term.counts <- rep(0, length(unique.terms))
  names(term.counts) <- unique.terms
  term.percent <- rep(0, length(unique.terms))
  names(term.percent) <- unique.terms
  for(i in 1:length(unique.terms)){
    term.counts[i] <- length(grep(unique.term.ids[i], terms))
    term.percent[i] <- round((term.counts[i] / length(terms))*100, 2) 
  }
  count.order <- order(term.counts, decreasing=T)
  
  output <- cbind(names(term.counts), unique.term.ids, term.counts, term.percent)[count.order,]
  colnames(output) <- c("Term", "GO ID", "Count", "Percent")
  
  if(!is.null(term.filter)){
    filter.bool <- unique.term.ids %in% term.filter
    output <- output[filter.bool,]
  }
  
  if(!is.null(output.path)){
    write.table(x=output, file=output.path, row.names=F, col.names=T, sep="\t")
  }
  print(output[1:25,])
  #sum(term.counts)
}

#++++++++++++++++++++++ giveTermClassi(...) ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ giveUnmatchedPeps(...) ++++++++++++++++++++++++++++++++++++

giveUnmatchedPeps <- function(pep.names=NULL, prot.vec=NULL){
  
  prot.accessions.list <- as.list(gsub("\\_peptide\\d+", "", pep.names))
  
  output <- c()
  for(i in 1:length(pep.names)){
    if(length(grep("\\; ", prot.accessions.list[[i]])>0)){
      tmp <- unlist(strsplit(prot.accessions.list[[i]], "\\; "))
      if(length(intersect(tmp, prot.vec))==0){
        output <- c(output, i)  
      }
    }else{
      if(!is.element(prot.accessions.list[[i]], prot.vec)){
        output <- c(output, i)
      }
    }
  }
  
  return(output)
  
}

#++++++++++++++++++++++ giveUnmatchedPeps(...) ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ kseaTP2EIL(...) ++++++++++++++++++++++++++++++++++++++++

# Merging kinase-substrate enrichment analysis results that are time point-specific
# (here: 6 time points) to results for "Early", "Intermediate" and "Late" time
# period-specific results. Current definitions:
# - Early: time points 1 & 2
# - Intermediate: time points 3 & 4
# - Late: time points 5 & 6

kseaTP2EIL <- function(ksea.results=NULL){

  ksea.results.EIL <- vector(mode="list", length=3)
  
  ksea.tmp1 <- ksea.results[1][[1]]
  ksea.tmp2 <- ksea.results[2][[1]]
  ksea.tmp3 <- ksea.results[3][[1]]
  ksea.tmp4 <- ksea.results[4][[1]]
  ksea.tmp5 <- ksea.results[5][[1]]
  ksea.tmp6 <- ksea.results[6][[1]]
  
  rownames(ksea.tmp1) <- ksea.tmp1[,"kinase"]
  rownames(ksea.tmp2) <- ksea.tmp2[,"kinase"]
  rownames(ksea.tmp3) <- ksea.tmp3[,"kinase"]
  rownames(ksea.tmp4) <- ksea.tmp4[,"kinase"]
  rownames(ksea.tmp5) <- ksea.tmp5[,"kinase"]
  rownames(ksea.tmp6) <- ksea.tmp6[,"kinase"]
  
  kinases <- sort(rownames(ksea.tmp1))
  
  for(i in 1:length(kinases)){
    if(as.numeric(ksea.tmp1[kinases[i],"pvalue"]) <= as.numeric(ksea.tmp2[kinases[i],"pvalue"])){
      ksea.results.EIL[1][[1]] <- rbind(ksea.results.EIL[1][[1]], ksea.tmp1[kinases[i],])  
    }else{
      ksea.results.EIL[1][[1]] <- rbind(ksea.results.EIL[1][[1]], ksea.tmp2[kinases[i],]) 
    }
    
    if(as.numeric(ksea.tmp3[kinases[i],"pvalue"]) <= as.numeric(ksea.tmp4[kinases[i],"pvalue"])){
      ksea.results.EIL[2][[1]] <- rbind(ksea.results.EIL[2][[1]], ksea.tmp3[kinases[i],])  
    }else{
      ksea.results.EIL[2][[1]] <- rbind(ksea.results.EIL[2][[1]], ksea.tmp4[kinases[i],]) 
    }
    
    if(as.numeric(ksea.tmp5[kinases[i],"pvalue"]) <= as.numeric(ksea.tmp6[kinases[i],"pvalue"])){
      ksea.results.EIL[3][[1]] <- rbind(ksea.results.EIL[3][[1]], ksea.tmp5[kinases[i],])  
    }else{
      ksea.results.EIL[3][[1]] <- rbind(ksea.results.EIL[3][[1]], ksea.tmp6[kinases[i],])
    }
  }
  
  return(ksea.results.EIL)
  
}

#++++++++++++++++++++++ kseaTP2EIL(...) ++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++ mapPepToProt(...) +++++++++++++++++++++++++++++++++++++++

mapPepToProt <- function(pep.ids=NULL, iso.rm=FALSE, group.rm=FALSE, output.path=NULL){
  
  prot.ids <- gsub("\\_peptide\\d+", "", pep.ids)
  if(iso.rm==TRUE){
    prot.ids <- unique(gsub("\\-\\d+", "", prot.ids))  
  }
  if(group.rm==TRUE){
    prot.groups.idx <- grep("\\; ", prot.ids)
    prot.groups <- grep("\\; ", prot.ids, value=TRUE)
    splitted.groups <- unlist(strsplit(grep("\\; ", prot.groups, value=TRUE), "\\; "))
    if(length(prot.groups.idx) > 0){
      prot.ids <- prot.ids[-prot.groups.idx]
      prot.ids <- unique(c(prot.ids, splitted.groups))
    }
  }
  if(!is.null(output.path)){
    write.table(x=prot.ids,
                file=paste0(output.path, "/protein.ids.txt"),
                row.names=FALSE,
                col.names=FALSE)
  }
  return(prot.ids)
  
}

#+++++++++++++++++++++ mapPepToProt(...) +++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ parseDEPOD(...) ++++++++++++++++++++++++++++++++++++++++

parseDEPOD <- function(file.path=NULL){

  DEPOD <- read.table(
    #file=paste0(cwd, "/ClueR/ClueR-PSEA/20220302_PPase_protSubtrates_201903_CORRECTED_BY_MT.txt"),
    file=file.path,
    header=TRUE,
    sep="\t",
    na.strings=c("NA", "N/A", "N/A ", "NULL"))

  colnames(DEPOD) <- c(
      "Phosphatase.UniProt",                                
      "Phosphatase.Entry",                                      
      "Substrate.UniProt",                                  
      "Substrate.Entry",                                        
      "Dephosphosites",                                               
      "AA.window",
      "Assay",                           
      "Reference",                                          
      "Reliability" 
  )
  
  DEPOD$Dephosphosites <- gsub("(\\.+) $", "\\1", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub("multiple", NA, DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub(",(S|T|H)", ", \\1", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub("Tyr-15 \\[Thr-14, Tyr-15\\]", "Thr-14, Tyr-15", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub("Ser-11 \\(Ser-10 in ref.\\)", "Ser-10, Ser-11", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub(";", ",", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub(" in [0-9A-Z]{6}", "", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub(" in isoform \\d+", "", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub(" in Isoform Tau-F \\(Tau-4\\)", "", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub("five serines of ", "", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub("and ", "", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub(" in isoform C", "", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub(" in Isoform TrkA-I", "", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub(" in ref.", "", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub(" \\(PRKCB2\\)", "", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub(" in Isoform Beta-II", "", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub(" in Isoform p52Shc", "", DEPOD$Dephosphosites)
  
  DEPOD$Dephosphosites <- gsub("His-", "H", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub("Ser-", "S", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub("Thr-", "T", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub("Tyr-", "Y", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub("His", "H", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub("Ser", "S", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub("Thr", "T", DEPOD$Dephosphosites)
  DEPOD$Dephosphosites <- gsub("Tyr", "Y", DEPOD$Dephosphosites)
  DEPOD <- DEPOD[order(DEPOD$Substrate.Entry),]
  DEPOD <- DEPOD[order(DEPOD$Phosphatase.Entry),]
  
  phosphatases <- unique(DEPOD$Phosphatase.Entry)
  DEPOD.PS <- vector(mode="list")
  
  for(i in 1:length(phosphatases)){
    DEPOD.tmp <- DEPOD[DEPOD$Phosphatase.Entry==phosphatases[i],]
    DEPOD.PS.tmp <- c()
    #print(DEPOD.tmp)
    
    for(j in 1:nrow(DEPOD.tmp)){
      #print(paste0("j: ", j))
      
      if(length(grep(",", DEPOD.tmp$Dephosphosites[j])) > 0){
        splitted.sites.tmp <- unlist(strsplit(DEPOD.tmp$Dephosphosites[j], ","))
        for(k in 1:length(splitted.sites.tmp)){
          substrate.tmp <- gsub(" ", "", DEPOD.tmp$Substrate.Entry[j])
          dephosphosite.tmp <- gsub(" ", "", splitted.sites.tmp[k])
          #print(paste0(substrate.tmp,";",dephosphosite.tmp,";"))
          DEPOD.PS.tmp <- c(DEPOD.PS.tmp, paste0(substrate.tmp,";",dephosphosite.tmp,";"))
        }
      }else{
        if(!is.na(DEPOD.tmp$Dephosphosites[j])){
          substrate.tmp <- gsub(" ", "", DEPOD.tmp$Substrate.Entry[j])
          dephosphosite.tmp <- gsub(" ", "", DEPOD.tmp$Dephosphosites[j])
          #print(paste0(substrate.tmp,";",dephosphosite.tmp,";"))
          DEPOD.PS.tmp <- c(DEPOD.PS.tmp, paste0(substrate.tmp,";",dephosphosite.tmp,";"))
        }
      }
    
    }
    #print(DEPOD.PS.tmp)
    if(!is.null(DEPOD.PS.tmp)){
      DEPOD.PS[[phosphatases[i]]] <- DEPOD.PS.tmp
    }
    #print(">>>>>>>>>>>>>>>>>>>>>>")
  }
  
  return(DEPOD.PS)
  
}

#++++++++++++++++++++++ parseDEPOD(...) ++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ pep2prot(...) ++++++++++++++++++++++++++++++++++++++++++

pep2prot <- function(pep=NULL){
  
  prot.tmp <- gsub("\\_peptide\\d+", "", pep)
  prot.tmp <- gsub("-\\d+", "", prot.tmp)
  prot <- c()
  for(i in 1:length(prot.tmp)){
    if(length(grep("\\; ", prot.tmp[i])>0)){
      tmp <- unlist(strsplit(prot.tmp[i], "\\; "))
      for(j in 1:length(tmp)){
          prot <- c(prot,tmp[j])  
      }
    }else{
      prot <- c(prot, prot.tmp[i])
    }
  }
  
  return(prot)
}

#++++++++++++++++++++++ pep2prot(...) ++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ performClueR(...) ++++++++++++++++++++++++++++++++++++

performClueR <- function(clueR.dat=NULL, organism=NULL, mode="kinase", KSorPSanno=NULL, rep=20, kRange=2:10, effectiveSize=c(5,100), pvalueCutoff=0.07, output.path=NULL){

  library(ClueR)
  
  
  if(mode=="kinase"){
    #library(PhosR)
    #data(PhosphoSite) #PhosphoSitePlus-ClueR-Version: 1) overlap phos-sites with MS-data: 233 2) kinases: 206 3) substrates: 9830
    #data(PhosphoSitePlus) #PhosphoSitePlus-PhosR-Version: 1) overlap phos-sites with MS-data: 283 2) kinases: 379 3) substrates: 11998
    PhosphoSite.human <- preparePhosSitePlus(organism=organism) #own PhosphoSitePlus-Version: 1) XXX 2) kinases: 403 3) substrates: 13664
    #PhosphoSite.mouse <- preparePhosSitePlus(organism="mouse")
    mode.folder <- "/ClueR/ClueR-K" 
  }else if(mode=="phosphatase"){
    PSanno <- KSorPSanno
    mode.folder <- "/ClueR/ClueR-P" 
  }
  
  
  rownames(clueR.dat) <- toupper(rownames(clueR.dat))
  
  
  if(organism == "human" && mode == "kinase") {
    clueObj <- ClueR::runClue(
                        Tc=clueR.dat,
                        annotation=PhosphoSite.human,
                        rep=rep,
                        kRange=kRange,
                        effectiveSize=effectiveSize,
                        pvalueCutoff=pvalueCutoff)
  }else if(organism == "human" && mode == "phosphatase"){
    clueObj <- ClueR::runClue(
                        Tc=clueR.dat,
                        annotation=PSanno,
                        rep=rep,
                        kRange=kRange,
                        effectiveSize=effectiveSize,
                        pvalueCutoff=pvalueCutoff)  
  }else if(organism == "mouse" && mode == "kinase"){
    clueObj <- ClueR::runClue(
                        Tc=clueR.dat,
                        annotation=PhosphoSite.mouse,
                        rep=rep,
                        kRange=kRange,
                        effectiveSize=effectiveSize,
                        pvalueCutoff=pvalueCutoff)  
  }else{
    stop("performClueR(...): Error! Not supported organism and/or analysis mode!")  
  }
      

  if(!is.null(output.path)){  
      save(clueObj,
          file=paste0(output.path, mode.folder, "/clueObj.RData"),
          compress = "xz",
          compression_level = 9)
  }
  
  xl <- "Number of clusters"
  yl <- "Enrichment score"
  png(filename=paste0(output.path, mode.folder, "/cluster_number_boxplots.png"), height=2000, width=2000, res=300)
    par(font.lab=2, font.axis=2, font.main=2)
    boxplot(
      clueObj$evlMat,
      col=rainbow(ncol(clueObj$evlMat)),
      las=2,
      cex.lab=1.25,
      cex.main=1.5,
      xlab=xl,
      ylab=yl,
      main="CLUE: Comparing diff. numbers of clusters"
    )
  dev.off()
  
  print(paste0("maxK: ", clueObj$maxK))
  
  png(filename=paste0(output.path, mode.folder, "/cluster_plots_maxK_best.png"), height=2000, width=2000, res=300)
      best <- ClueR::clustOptimal(clueObj, rep=rep/2, mfrow=c(2, 3))
  dev.off()
  for(i in 1:length(best$enrichList)){
    write.table(
      x=best$enrichList[[i]],
      file=paste0(output.path, mode.folder, "/best_enrichList_", names(best$enrichList)[i], ".txt"),
      row.names=FALSE,
      sep="\t"
    )
  }
  
  #user.maxK.vec <- c(3,4,5,6,7,8,9)
  user.maxK.vec <- kRange
  for(i in 1:length(user.maxK.vec)){
    png(filename=paste0(output.path, mode.folder, "/cluster_plots_", user.maxK.vec[i], ".png"), height=2000, width=2000, res=300)
      if(user.maxK.vec[i] == 2) clustOptimal(clueObj, rep=5, user.maxK=user.maxK.vec[i], mfrow=c(1, 2))
      else if(user.maxK.vec[i] <= 4) clustOptimal(clueObj, rep=5, user.maxK=user.maxK.vec[i], mfrow=c(2, 2))
      else if(user.maxK.vec[i] <= 6) clustOptimal(clueObj, rep=5, user.maxK=user.maxK.vec[i], mfrow=c(3, 2))
      else if(user.maxK.vec[i] <= 9) clustOptimal(clueObj, rep=5, user.maxK=user.maxK.vec[i], mfrow=c(3, 3))
      else if(user.maxK.vec[i] <= 12) clustOptimal(clueObj, rep=5, user.maxK=user.maxK.vec[i], mfrow=c(4, 3))
    dev.off()  
  }
  
  return(best)
}

#++++++++++++++++++++++ performClueR(...) ++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ phosphoPepsToSites(...) ++++++++++++++++++++++++++++++++

phosphoPepsToSites <- function(dat=NULL,
                             mods.master=NULL,
                             mods=NULL,
                             psm.number=NULL,
                             reg.peps=NULL,
                             quantified.only=TRUE,
                             localized.only=TRUE,
                             log=TRUE,
                             output.path=NULL,
                             file.name="phossite.dat"){
  
  
  if(is.null(reg.peps)){
    dat.reg <- dat
    mods.master.reg <- mods.master
    psm.number.reg <- psm.number
  }else{
    dat.reg <- dat[reg.peps,]
    mods.master.reg <- mods.master[reg.peps]
    psm.number.reg <- psm.number[reg.peps]
  }
  
  #-----------------------------------------------------------------------------
  # BEGIN: Split site/peptide group!
  # E.g., for "Q8IYB3; Q9UPU5_peptide1" the phosphosite
  # "Q8IYB3 1xPhospho [S653(100)]; Q9UPU5 1xPhospho [S1141(100)]" was found.
  # This is splitted into:
  # - "Q8IYB3 1xPhospho [S653(100)]" 
  # - "Q8IYB3 1xPhospho [S653(100)]"
  #-----------------------------------------------------------------------------
  
  phossite.dat.tmp <- matrix(ncol=ncol(dat.reg))
  psm.number.reg.clueR.tmp <- c()
  #test.idx <- c(85:100, 205:215)
  #test.idx <- 1100:1200
  #for(i in test.idx){
  for(i in 1:length(mods.master.reg)){
    print(mods.master.reg[i])
    if(length(grep("]; ", mods.master.reg[i])) > 0){
      splitted.mods <- unlist(strsplit(mods.master.reg[i], "]; "))
      splitted.mods[1:(length(splitted.mods)-1)] <- paste0(splitted.mods[1:(length(splitted.mods)-1)], "]")
      
      last.uniprot.id <- ""
      for(j in 1:length(splitted.mods)){
        if(length(grep("^\\d+x", splitted.mods[j]) > 0)){
          splitted.mods[j] <- paste0(last.uniprot.id, " ", splitted.mods[j])
        }else{
          last.uniprot.id <- gsub("^([0-9A-Z]+-*\\d*) \\d+x.+", "\\1", splitted.mods[j])
        }  
      }
      
      for(j in 1:length(splitted.mods)){
        phossite.dat.tmp <- rbind(phossite.dat.tmp, dat.reg[i,]) 
        rownames(phossite.dat.tmp)[nrow(phossite.dat.tmp)] <- splitted.mods[j]
        print(rownames(phossite.dat.tmp)[nrow(phossite.dat.tmp)])
        psm.number.reg.clueR.tmp <- c(psm.number.reg.clueR.tmp, psm.number.reg[i])
        names(psm.number.reg.clueR.tmp)[length(psm.number.reg.clueR.tmp)] <- splitted.mods[j]
      }
    }else if(length(grep("]; ", mods.master.reg[i])) == 0 && !is.na(mods.master.reg[i])){
      phossite.dat.tmp <- rbind(phossite.dat.tmp, dat.reg[i,]) 
      rownames(phossite.dat.tmp)[nrow(phossite.dat.tmp)] <- mods.master.reg[i]
      print(rownames(phossite.dat.tmp)[nrow(phossite.dat.tmp)])
      psm.number.reg.clueR.tmp <- c(psm.number.reg.clueR.tmp, psm.number.reg[i])
      names(psm.number.reg.clueR.tmp)[length(psm.number.reg.clueR.tmp)] <- mods.master.reg[i]
    }else if(is.na(mods.master.reg[i]) && localized.only == F){
      mod.tmp <- gsub(".*\\dxPhospho \\[(.+)\\].*$", "\\1", mods[i])
      pep.idx.tmp <- gsub(".+_peptide(\\d+)$", "\\1", names(mods[i]))
      uniprot.id.tmp <- gsub("^(.+)_peptide\\d+$", "\\1", names(mods[i]), perl=TRUE)
      if(length(grep("\\;", uniprot.id.tmp)) > 0){
        uniprot.id.tmp.vec <- unlist(strsplit(uniprot.id.tmp, "; ")) 
        for(j in 1:length(uniprot.id.tmp.vec)){
          phossite.dat.tmp <- rbind(phossite.dat.tmp, dat.reg[i,]) 
          rownames(phossite.dat.tmp)[nrow(phossite.dat.tmp)] <- paste0("unloc_", uniprot.id.tmp.vec[j], "_", mod.tmp, "_", pep.idx.tmp)
          print(rownames(phossite.dat.tmp)[nrow(phossite.dat.tmp)])
          psm.number.reg.clueR.tmp <- c(psm.number.reg.clueR.tmp, psm.number.reg[i])
          names(psm.number.reg.clueR.tmp)[length(psm.number.reg.clueR.tmp)] <- mods.master.reg[i]
        }
        #readline(prompt="Here-3!!! Press [enter] to continue!")
      }else{
        phossite.dat.tmp <- rbind(phossite.dat.tmp, dat.reg[i,]) 
        rownames(phossite.dat.tmp)[nrow(phossite.dat.tmp)] <- paste0("unloc_", uniprot.id.tmp, "_", mod.tmp, "_", pep.idx.tmp)
        print(rownames(phossite.dat.tmp)[nrow(phossite.dat.tmp)])
        psm.number.reg.clueR.tmp <- c(psm.number.reg.clueR.tmp, psm.number.reg[i])
        names(psm.number.reg.clueR.tmp)[length(psm.number.reg.clueR.tmp)] <- mods.master.reg[i]
      }
    }
    print("****************************************************************")
  }
  if(sum(is.na(phossite.dat.tmp[1,])) == ncol(phossite.dat.tmp)) phossite.dat.tmp <- phossite.dat.tmp[-1,] # <---------- ????
  #print(phossite.dat.tmp)
  #readline(prompt="Press [enter] to continue!")
  
  #-----------------------------------------------------------------------------
  # END: Split site/peptide group!
  #-----------------------------------------------------------------------------
  
  phossite.dat <- matrix(ncol=ncol(phossite.dat.tmp))
  psm.number.reg.clueR <- c()
  for(i in 1:nrow(phossite.dat.tmp)){
    print(paste0("original: ", rownames(phossite.dat.tmp)[i]))
    if(length(grep("); ", rownames(phossite.dat.tmp)[i])) > 0){
      uniprot.id <- gsub("^([0-9A-Z]+-*\\d*) .*", "\\1", rownames(phossite.dat.tmp)[i], perl=TRUE)
      tmp <- gsub("^.* \\d+xPhospho ", "", rownames(phossite.dat.tmp)[i], perl=TRUE)
      tmp <- gsub("\\(\\d+\\.*\\d*\\)", "", tmp)
      tmp <- gsub("(\\[|\\])", "", tmp)
      splitted.sites <- unlist(strsplit(tmp, "; "))
      for(j in 1:length(splitted.sites)){
        print(paste0(uniprot.id, ":", splitted.sites[j]))
        phossite.dat <- rbind(phossite.dat, phossite.dat.tmp[i,])
        rownames(phossite.dat)[nrow(phossite.dat)] <- paste0(uniprot.id, ":", splitted.sites[j]) 
        print(paste0("PSM number: ", psm.number.reg.clueR.tmp[i]))
        psm.number.reg.clueR <- c(psm.number.reg.clueR, psm.number.reg.clueR.tmp[i])
        names(psm.number.reg.clueR)[length(psm.number.reg.clueR)] <- paste0(uniprot.id, ":", splitted.sites[j]) 
      }
    }else if(length(grep("^unloc_", rownames(phossite.dat.tmp)[i])) > 0){
      #uniprot.id <- gsub("^unloc_([0-9A-Z]+-*\\d*)_.+", "\\1", rownames(phossite.dat.tmp)[i], perl=TRUE)
      uniprot.id <- gsub("^unloc_(.+)_.+_\\d+$", "\\1", rownames(phossite.dat.tmp)[i], perl=TRUE)
      site <- gsub(".+_(.+_\\d+)$", "\\1", rownames(phossite.dat.tmp)[i], perl=TRUE)
      
      print(paste0(uniprot.id, ":", site))
      phossite.dat <- rbind(phossite.dat, phossite.dat.tmp[i,])
      rownames(phossite.dat)[nrow(phossite.dat)] <- paste0(uniprot.id, ":", site)
      print(paste0("PSM number: ", psm.number.reg.clueR.tmp[i]))
      psm.number.reg.clueR <- c(psm.number.reg.clueR, psm.number.reg.clueR.tmp[i])
      names(psm.number.reg.clueR)[length(psm.number.reg.clueR)] <- paste0(uniprot.id, ":", site)
      
      #if(length(grep("\\;", uniprot.id)) > 0) readline(prompt="Here-2!!! Press [enter] to continue!")
    }else if(length(grep("); ", rownames(phossite.dat.tmp)[i])) == 0){
      uniprot.id <- gsub("^([0-9A-Z]+-*\\d*) .*", "\\1", rownames(phossite.dat.tmp)[i], perl=TRUE)
      site <- gsub(".*\\[((S|T|Y){1}\\d+)\\(\\d+\\.*\\d*\\)\\]$", "\\1", rownames(phossite.dat.tmp)[i], perl=TRUE)
      
      if(length(grep("xPhospho", site)) > 0){
        print(paste0(uniprot.id, ":", site))
        phos.count <- as.numeric(gsub(".+ (\\d+)xPhospho .+", "\\1", site, perl=T))
        actual.site <- gsub(".+\\d+xPhospho \\[(.+)\\]$", "\\1", site, perl=T)
        print(paste0("phos.count: ", phos.count))
        print(paste0("actual.site: ", actual.site))
        for(j in 1:phos.count){
          print(paste0(uniprot.id, ":", actual.site, "-", j))
          phossite.dat <- rbind(phossite.dat, phossite.dat.tmp[i,])  
          rownames(phossite.dat)[nrow(phossite.dat)] <- paste0(uniprot.id, ":", actual.site, "-", j)
          print(paste0("PSM number: ", psm.number.reg.clueR.tmp[i]))
          psm.number.reg.clueR <- c(psm.number.reg.clueR, psm.number.reg.clueR.tmp[i])
          names(psm.number.reg.clueR)[length(psm.number.reg.clueR)] <- paste0(uniprot.id, ":", actual.site, "-", j)
        }
        #readline(prompt="HERE-1. Press [enter] to continue!")  
      }else{
        print(paste0(uniprot.id, ":", site))
        phossite.dat <- rbind(phossite.dat, phossite.dat.tmp[i,])
        rownames(phossite.dat)[nrow(phossite.dat)] <- paste0(uniprot.id, ":", site)
        print(paste0("PSM number: ", psm.number.reg.clueR.tmp[i]))
        psm.number.reg.clueR <- c(psm.number.reg.clueR, psm.number.reg.clueR.tmp[i])
        names(psm.number.reg.clueR)[length(psm.number.reg.clueR)] <- paste0(uniprot.id, ":", site)
      }
      
    }
    print("------------------------------------------------------")
  }
  #print(grep("(S|T|Y)(-|_)\\d+$",rownames(phossite.dat),value=T))
  
  # VVV 1st row of phossite.dat always a NA-row --> ALWAYS NEEDED! VVV
  if(sum(is.na(phossite.dat[1,])) == ncol(phossite.dat)) {
    del.idx1 <- c(1)
    phossite.dat <- phossite.dat[-del.idx1,]
    #print(nrow(phossite.dat))
    #print(length(psm.number.reg.clueR)) #now should be = nrow(phossite.dat)
  }
  # ^^^ 1st row of phossite.dat always a NA-row --> ALWAYS NEEDED! ^^^
  
  # VVV Discard cases like 'P20042:P20042 1xAcetyl [N-Term]' VVV
  if(length(grep("\\[N-Term\\]",rownames(phossite.dat))) > 0) {
    del.idx2 <- grep("\\[N-Term\\]",rownames(phossite.dat))
    #print(rownames(phossite.dat)[del.idx2])
    #print(names(psm.number.reg.clueR)[del.idx2])
    phossite.dat <- phossite.dat[-del.idx2,]
    psm.number.reg.clueR <- psm.number.reg.clueR[-del.idx2]
  }
  # ^^^ VVV Discard cases like '' and '' ^^^
  
  # VVV Discard cases like Q5SW79:S/T/Y_2' & 'Q9UQ35:S/Y/T-1' VVV
  if(localized.only == T && length(grep("(S|T|Y)(-|_)\\d+$",rownames(phossite.dat))) > 0){
    del.idx3 <- grep("(S|T|Y)(-|_)\\d+$",rownames(phossite.dat))
    print(rownames(phossite.dat)[del.idx3])
    phossite.dat <- phossite.dat[-del.idx3,]
    psm.number.reg.clueR <- psm.number.reg.clueR[-del.idx3]
  }
  #readline(prompt="#NOT.QUANT! Press [enter] to continue!")
  # ^^^ VVV Discard cases like '' and '' ^^^
  
  # VVV Discard cases like Q5SW79:S/T/Y' & 'Q9UQ35:S' VVV
  if(localized.only == T && length(grep("\\:(S|T|Y)$",rownames(phossite.dat))) > 0){
    print("Here - 2")
    del.idx3 <- grep("(S|T|Y)$",rownames(phossite.dat))
    print(rownames(phossite.dat)[del.idx3])
    phossite.dat <- phossite.dat[-del.idx3,]
    psm.number.reg.clueR <- psm.number.reg.clueR[-del.idx3]
  }
  #readline(prompt="#NOT.QUANT! Press [enter] to continue!")
  # ^^^ VVV Discard cases like '' and '' ^^^
  
  if(quantified.only == T){
    del.idx4 <- c()
    for(i in 1:nrow(phossite.dat)){
      if(sum(is.na(phossite.dat[i,])) == ncol(phossite.dat)) del.idx4 <- c(del.idx4, i) 
    }
    print(del.idx4)
    if(!is.null(del.idx4)){
      phossite.dat <- phossite.dat[-del.idx4,]
      psm.number.reg.clueR <- psm.number.reg.clueR[-del.idx4]
    }
    #readline(prompt="#NOT.QUANT! Press [enter] to continue!")
  }
  
  del.idx <- c()
  unique.rows <- unique(rownames(phossite.dat))
  for(i in 1:length(unique.rows)){
    unique.row.re <- paste0("^", unique.rows[i],"$")
    if(length(grep(unique.row.re, rownames(phossite.dat))) > 1){
      tmp.idx.vec <- grep(unique.row.re, rownames(phossite.dat))
      tmp.idx.vec2 <- grep(unique.row.re, rownames(psm.number.reg.clueR))
      print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
      print("unique.rows:")
      print(unique.rows[i])
      print("tmp.idx.vec:")
      print(tmp.idx.vec)
      print("rownames(phossite.dat)[tmp.idx.vec]:")
      print(rownames(phossite.dat)[tmp.idx.vec])
      print("PSM-Numbers:")
      print(psm.number.reg.clueR[tmp.idx.vec])
      keep.idx <- which.max(psm.number.reg.clueR[tmp.idx.vec])
      print("keep.idx:")
      print(keep.idx)
      tmp.idx.vec <- tmp.idx.vec[-keep.idx]
      print("tmp.idx.vec:")
      print(tmp.idx.vec)
      del.idx <- c(del.idx, tmp.idx.vec)
    }
  }
  phossite.dat <- phossite.dat[-del.idx,]
  
  if(log==TRUE) phossite.dat <- log2(phossite.dat)
  colnames(phossite.dat) <- c("0min", "1min", "2.5min", "5min", "15min", "30min", "60min")
  
  if(is.null(output.path) == FALSE){
    write.table(cbind(rownames(phossite.dat),phossite.dat), file=paste0(output.path, "/", file.name, ".txt"), sep="\t", col.names = TRUE, row.names = FALSE)
  }
  
  return(phossite.dat)
}

#++++++++++++++++++++++ phosphoPepsToSites(...) ++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ plotClusterClueR(...) ++++++++++++++++++++++++++++++++++

plotClusterClueR <- function(clueR.dat=NULL, best=NULL, scale=TRUE, unlog=FALSE, mode="kinase", output.path=NULL){
  
  if(!is.null(output.path) && mode=="kinase") {
    output.path.tmp <- paste0(output.path, "/ClueR/ClueR-K/own_cluster_plots") 
    
  }else if(!is.null(output.path) && mode=="phosphatase"){
    output.path.tmp <- paste0(output.path, "/ClueR/ClueR-P/own_cluster_plots")
  }
  if(!is.null(output.path) && !dir.exists(output.path.tmp)) {
    dir.create(output.path.tmp)  
  }
  clust.n <- length(best$clustObj$size)
  
  for(i in 1:clust.n){
    clust.id <- paste0("^", i, "$")
    idx.tmp <- grep(clust.id, best$clustObj$cluster)
    dat.tmp <- clueR.dat[idx.tmp,]
    if(unlog == TRUE){
      dat.tmp <- 2^dat.tmp
    }
    
    if(scale == TRUE) {
      y.min <- min(t(apply(dat.tmp, 1, scale)))-0.5
      y.max <- max(t(apply(dat.tmp, 1, scale)))+0.5
      png(filename=paste0(output.path.tmp, "/own_cluster_plot_", i, ".png"), height=2000, width=2000, res=300)
        plot(
          scale(dat.tmp[1,]),
          type="l",
          ylim=c(y.min,y.max),
          main=paste0("Cluster ", i, "; size = ", nrow(dat.tmp)),
          xlab="Time course",
          ylab="Standardized profile",
          col=alpha("red",alpha=0.05),
          cex=3,
          cex.main=2,
          cex.lab=1.5,
          cex.axis=1.5,
          font.main=2,
          font.lab=2,
          font.axis=2,
          lwd=3)
        for(j in 2:nrow(dat.tmp)) {
          lines(scale(dat.tmp[j,]), type="l", col=alpha("red", alpha=0.05), lwd=3)
        }
      dev.off()
      
      output.tmp <- cbind(rownames(dat.tmp), dat.tmp)
      colnames(output.tmp)[1] <- "Site"
      write.table(
        x=output.tmp,
        file=paste0(output.path.tmp, "/profiles_cluster_", i, ".txt"),
        row.names=FALSE,
        sep="\t")
    }else{
      y.min <- min(dat.tmp)-0.5
      y.max <- max(dat.tmp)+0.5
      png(filename=paste0(output.path.tmp, "/own_cluster_plot_", i, ".png"), height=2000, width=2000, res=300)
      plot(
        dat.tmp[1,],
        type="l",
        main=paste0("Cluster ", i, "; size = ", nrow(dat.tmp)),
        xlab="Time course",
        ylab="Standardized profile",
        ylim=c(y.min,y.max),
        col=alpha("red",alpha=0.05),
        cex=3,
        cex.main=2,
        cex.lab=1.5,
        cex.axis=1.5,
        font.main=2,
        font.lab=2,
        font.axis=2,
        lwd=3)
      for(j in 2:nrow(dat.tmp)) {
        lines(dat.tmp[j,], type="l", col=alpha("red", alpha=0.05), lwd=3)
      }
      dev.off()
      
      output.tmp <- cbind(rownames(dat.tmp), dat.tmp)
      colnames(output.tmp)[1] <- "Site"
      write.table(
        x=output.tmp,
        file=paste0(output.path.tmp, "/profiles_cluster_", i, ".txt"),
        row.names=FALSE,
        sep="\t")
    }
    
  }
  
}

#++++++++++++++++++++++ plotClusterClueR(...) ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ plotPCA(...) +++++++++++++++++++++++++++++++++++++++++++

plotPCA <- function(x=NULL, group.legend.pos=NULL, idv.legend.pos=NULL, filename.prefix=NULL, output.path=NULL){
  
  pca.dat <- x
  if(!is.null(output.path) && !dir.exists(paste0(output.path, "/PCA"))) {
    dir.create(paste0(output.path, "/PCA"))
  }
  output.path <- paste0(output.path, "/PCA")
  
  print(paste0("original feature number: ", nrow(pca.dat)))
  pca.dat <- na.omit(pca.dat)
  print(paste0("feature number after removing features containing NAs: ", nrow(pca.dat)))
  
  idv.vec <- c("AR7", "AT13", "AT18", "AT22", "AT8")
  idv.vec2 <- c("A", "B", "C", "D", "E")
  pch.vec <- c(21,22,23,24,25)
  group0 <- grep("_0$", colnames(pca.dat), value=TRUE)
  group1 <- grep("_1$", colnames(pca.dat), value=TRUE)
  group2 <- grep("_2.5$", colnames(pca.dat), value=TRUE)
  group3 <- grep("_5$", colnames(pca.dat), value=TRUE)
  group4 <- grep("_15$", colnames(pca.dat), value=TRUE)
  group5 <- grep("_30$", colnames(pca.dat), value=TRUE)
  group6 <- grep("_60$", colnames(pca.dat), value=TRUE)
  n0 <- length(group0)
  n1 <- length(group1)
  n2 <- length(group2)
  n3 <- length(group3)
  n4 <- length(group4)
  n5 <- length(group5) 
  n6 <- length(group6) 
  
  group.vec <- c(
    "0 min",
    "1 min",
    "2.5 min",
    "5 min",
    "15 min",
    "30 min",
    "60 min"
  )
  
  col.vec <- c(
    adjustcolor("black", alpha=0.3),
    adjustcolor("red", alpha=0.3),
    adjustcolor("darkred", alpha=0.3),
    adjustcolor("deepskyblue4", alpha=0.3),
    adjustcolor("navy", alpha=0.3),
    adjustcolor("chartreuse3", alpha=0.3),
    adjustcolor("darkgreen", alpha=0.3),
    adjustcolor("magenta", alpha=0.3),
    adjustcolor("darkorange", alpha=0.3),
    adjustcolor("darkorchid", alpha=0.3),
    adjustcolor("darkolivegreen4", alpha=0.3),
    adjustcolor("gold", alpha=0.3),
    adjustcolor("plum", alpha=0.3),
    adjustcolor("yellow", alpha=0.3),
    adjustcolor("khaki", alpha=0.3),
    adjustcolor("seagreen", alpha=0.3),
    adjustcolor("saddlebrown", alpha=0.3),
    adjustcolor("lavender", alpha=0.3)
  )
  col.vec2 <- c(
    adjustcolor("black", alpha=0.8),
    adjustcolor("red", alpha=0.8),
    adjustcolor("darkred", alpha=0.8),
    adjustcolor("deepskyblue4", alpha=0.8),
    adjustcolor("navy", alpha=0.8),
    adjustcolor("chartreuse3", alpha=0.8),
    adjustcolor("darkgreen", alpha=0.8)
  )
  
  
  pcdat <- prcomp(t(pca.dat), center=TRUE, scale=TRUE)
  scores <- pcdat$x
  for (i in 1:2){
    for (j in i:2){
      if (i<j){
        XLIM <- c(-max(abs(scores[,i])), max(abs(scores[,i])))
        XLIM <- XLIM+(XLIM*0.1)
        YLIM <- c(-max(abs(scores[,j])), max(abs(scores[,j])))
        YLIM <- YLIM+(YLIM*0.1)
        #print(summary(pcdat)$importance)
        Xlab <- paste0("PC", i, " (", round(summary(pcdat)$importance[2,i]*100, 2), "%)")
        Ylab <- paste0("PC", j, " (", round(summary(pcdat)$importance[2,j]*100, 2), "%)")
        png(paste(output.path, "/", filename.prefix, "01_pca_", i, "_", j, "_filteredFeatures.png", sep=""), width=3600, height=3600, pointsize=15, res=600)
          plot(scores[group0,i], scores[group0,j], type="n", xlab=Xlab, ylab=Ylab, xlim=XLIM, ylim=YLIM, pch=15, col=col.vec[1], cex=2)
          
          for(k in 1:length(idv.vec)){
            points(scores[grep(idv.vec[k], group0, value=T),i], scores[grep(idv.vec[k], group0, value=T),j], pch=pch.vec[k], col=col.vec[1], bg=col.vec[1], cex=2)  
            points(scores[grep(idv.vec[k], group1, value=T),i], scores[grep(idv.vec[k], group1, value=T),j], pch=pch.vec[k], col=col.vec[2], bg=col.vec[2], cex=2)
            points(scores[grep(idv.vec[k], group2, value=T),i], scores[grep(idv.vec[k], group2, value=T),j], pch=pch.vec[k], col=col.vec[3], bg=col.vec[3], cex=2)
            points(scores[grep(idv.vec[k], group3, value=T),i], scores[grep(idv.vec[k], group3, value=T),j], pch=pch.vec[k], col=col.vec[4], bg=col.vec[4], cex=2)
            points(scores[grep(idv.vec[k], group4, value=T),i], scores[grep(idv.vec[k], group4, value=T),j], pch=pch.vec[k], col=col.vec[5], bg=col.vec[5], cex=2)
            points(scores[grep(idv.vec[k], group5, value=T),i], scores[grep(idv.vec[k], group5, value=T),j], pch=pch.vec[k], col=col.vec[6], bg=col.vec[6], cex=2)
            points(scores[grep(idv.vec[k], group6, value=T),i], scores[grep(idv.vec[k], group6, value=T),j], pch=pch.vec[k], col=col.vec[7], bg=col.vec[7], cex=2)
          }
          
          legend(group.legend.pos, legend=group.vec[1:7], col=col.vec[1:7], pch=16, title="Time points", cex=0.75, pt.cex=1.5, bg="transparent", bty="n")
          legend(idv.legend.pos, legend=idv.vec2[1:5], col="black", pch=pch.vec[1:5], title="Donors", cex=0.75, pt.cex=1.0, pt.lwd=1.25, pt.bg="grey90", bg="transparent", bty="n")
        dev.off()
      }
    }
  }
  
  for (i in 1:2){
    for (j in i:2){
      if (i<j){
        XLIM <- c(-max(abs(scores[,i])), max(abs(scores[,i])))
        XLIM <- XLIM+(XLIM*0.1)
        YLIM <- c(-max(abs(scores[,j])), max(abs(scores[,j])))
        YLIM <- YLIM+(YLIM*0.1)
        #print(summary(pcdat)$importance)
        Xlab <- paste0("PC", i, " (", round(summary(pcdat)$importance[2,i]*100, 2), "%)")
        Ylab <- paste0("PC", j, " (", round(summary(pcdat)$importance[2,j]*100, 2), "%)")
        png(paste(output.path, "/", filename.prefix, "02_pca_", i, "_", j, "_filteredFeatures.png", sep=""), width=3600, height=3600, pointsize=15, res=600)
          plot(scores[group0,i], scores[group0,j], type="n", xlab=Xlab, ylab=Ylab, xlim=XLIM, ylim=YLIM, pch=16, col=col.vec[1], cex=2)
          
          for(k in 1:length(idv.vec)){
            points(scores[grep(idv.vec[k], group0, value=T),i], scores[grep(idv.vec[k], group0, value=T),j], pch=pch.vec[k], col=col.vec[1], bg=col.vec[1], cex=2)  
            points(scores[grep(idv.vec[k], group1, value=T),i], scores[grep(idv.vec[k], group1, value=T),j], pch=pch.vec[k], col=col.vec[2], bg=col.vec[2], cex=2)
            points(scores[grep(idv.vec[k], group2, value=T),i], scores[grep(idv.vec[k], group2, value=T),j], pch=pch.vec[k], col=col.vec[3], bg=col.vec[3], cex=2)
            points(scores[grep(idv.vec[k], group3, value=T),i], scores[grep(idv.vec[k], group3, value=T),j], pch=pch.vec[k], col=col.vec[4], bg=col.vec[4], cex=2)
            points(scores[grep(idv.vec[k], group4, value=T),i], scores[grep(idv.vec[k], group4, value=T),j], pch=pch.vec[k], col=col.vec[5], bg=col.vec[5], cex=2)
            points(scores[grep(idv.vec[k], group5, value=T),i], scores[grep(idv.vec[k], group5, value=T),j], pch=pch.vec[k], col=col.vec[6], bg=col.vec[6], cex=2)
            points(scores[grep(idv.vec[k], group6, value=T),i], scores[grep(idv.vec[k], group6, value=T),j], pch=pch.vec[k], col=col.vec[7], bg=col.vec[7], cex=2)
            text(scores[grep(idv.vec[k], group0, value=T),i], scores[grep(idv.vec[k], group0, value=T),j], labels=group0[k], col=col.vec2[1], cex=0.4)
            text(scores[grep(idv.vec[k], group1, value=T),i], scores[grep(idv.vec[k], group1, value=T),j], labels=group1[k], col=col.vec2[2], cex=0.4)
            text(scores[grep(idv.vec[k], group2, value=T),i], scores[grep(idv.vec[k], group2, value=T),j], labels=group2[k], col=col.vec2[3], cex=0.4)
            text(scores[grep(idv.vec[k], group3, value=T),i], scores[grep(idv.vec[k], group3, value=T),j], labels=group3[k], col=col.vec2[4], cex=0.4)
            text(scores[grep(idv.vec[k], group4, value=T),i], scores[grep(idv.vec[k], group4, value=T),j], labels=group4[k], col=col.vec2[5], cex=0.4)
            text(scores[grep(idv.vec[k], group5, value=T),i], scores[grep(idv.vec[k], group5, value=T),j], labels=group5[k], col=col.vec2[6], cex=0.4)
            text(scores[grep(idv.vec[k], group6, value=T),i], scores[grep(idv.vec[k], group6, value=T),j], labels=group6[k], col=col.vec2[7], cex=0.4)
          }
          
          legend(group.legend.pos, legend=group.vec[1:7], col=col.vec[1:7], pch=16, title="Time points", cex=0.75, pt.cex=1.5, bg="transparent", bty="n")
          legend(idv.legend.pos, legend=idv.vec2[1:5], col="black", pch=pch.vec[1:5], title="Donors", cex=0.75, pt.cex=1.0, pt.lwd=1.25, pt.bg="grey90", bg="transparent", bty="n")
        dev.off()
      }
    }
  }  
  
}

#++++++++++++++++++++++ plotPCA(...) +++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ plotPCA2(...) +++++++++++++++++++++++++++++++++++++++++++

plotPCA2 <- function(x=NULL, group.legend.pos=NULL, idv.legend.pos=NULL, filename.prefix=NULL, output.path=NULL){
  
  pca.dat <- x
  if(!is.null(output.path) && !dir.exists(paste0(output.path, "/PCA"))) {
    dir.create(paste0(output.path, "/PCA"))
  }
  output.path <- paste0(output.path, "/PCA")
  
  print(paste0("original feature number: ", nrow(pca.dat)))
  #pca.dat[is.na(pca.dat)] <- 1
  pca.dat <- na.omit(pca.dat)
  print(paste0("feature number after removing features containing NAs: ", nrow(pca.dat)))
  pca.dat <- log2(pca.dat)
  
  #idv.vec <- c("AR7", "AT13", "AT18", "AT22", "AT8")
  #idv.vec2 <- c("A", "B", "C", "D", "E")
  #pch.vec <- c(21,22,23,24,25)
  group1 <- grep("LEAN\\d+_0$", colnames(pca.dat), value=TRUE)
  group2 <- grep("LEAN\\d+_1$", colnames(pca.dat), value=TRUE)
  group3 <- grep("ADIP\\d+_0$", colnames(pca.dat), value=TRUE)
  group4 <- grep("ADIP\\d+_1$", colnames(pca.dat), value=TRUE)
  group5 <- grep("T2D\\d+_0$", colnames(pca.dat), value=TRUE)
  group6 <- grep("T2D\\d+_1$", colnames(pca.dat), value=TRUE)
  n1 <- length(group1)
  n2 <- length(group2)
  n3 <- length(group3)
  n4 <- length(group4)
  n5 <- length(group5) 
  n6 <- length(group6) 
  
  group.vec <- c(
    "lean-basal",
    "lean-clamp",
    "adipose-basal",
    "adipose-clmap",
    "diabetic-basal",
    "diabetic-clamp"
  )
  
  col.vec <- c(
    adjustcolor("black", alpha=0.3),
    adjustcolor("red", alpha=0.3),
    adjustcolor("darkred", alpha=0.3),
    adjustcolor("deepskyblue4", alpha=0.3),
    adjustcolor("navy", alpha=0.3),
    adjustcolor("chartreuse3", alpha=0.3),
    adjustcolor("darkgreen", alpha=0.3),
    adjustcolor("magenta", alpha=0.3),
    adjustcolor("darkorange", alpha=0.3),
    adjustcolor("darkorchid", alpha=0.3),
    adjustcolor("darkolivegreen4", alpha=0.3),
    adjustcolor("gold", alpha=0.3),
    adjustcolor("plum", alpha=0.3),
    adjustcolor("yellow", alpha=0.3),
    adjustcolor("khaki", alpha=0.3),
    adjustcolor("seagreen", alpha=0.3),
    adjustcolor("saddlebrown", alpha=0.3),
    adjustcolor("lavender", alpha=0.3)
  )
  col.vec2 <- c(
    adjustcolor("black", alpha=0.8),
    adjustcolor("red", alpha=0.8),
    adjustcolor("darkred", alpha=0.8),
    adjustcolor("deepskyblue4", alpha=0.8),
    adjustcolor("navy", alpha=0.8),
    adjustcolor("chartreuse3", alpha=0.8),
    adjustcolor("darkgreen", alpha=0.8)
  )
  
  
  pcdat <- prcomp(t(pca.dat), center=TRUE, scale=TRUE)
  scores <- pcdat$x
  
  for (i in 1:2){
    for (j in i:2){
      if (i<j){
        XLIM <- c(-max(abs(scores[,i])), max(abs(scores[,i])))
        XLIM <- XLIM+(XLIM*0.1)
        YLIM <- c(-max(abs(scores[,j])), max(abs(scores[,j])))
        YLIM <- YLIM+(YLIM*0.1)
        #print(summary(pcdat)$importance)
        Xlab <- paste0("PC", i, " (", round(summary(pcdat)$importance[2,i]*100, 2), "%)")
        Ylab <- paste0("PC", j, " (", round(summary(pcdat)$importance[2,j]*100, 2), "%)")
        png(paste(output.path, "/", filename.prefix, "02_pca_", i, "_", j, "_filteredFeatures.png", sep=""), width=3600, height=3600, pointsize=15, res=600)
        plot(scores[group1,i], scores[group1,j], type="n", xlab=Xlab, ylab=Ylab, xlim=XLIM, ylim=YLIM, pch=16, col=col.vec[1], cex=2)
        
        points(scores[group1,i], scores[group1,j], pch=21, col=col.vec[2], bg=col.vec[2], cex=2)
        points(scores[group2,i], scores[group2,j], pch=24, col=col.vec[3], bg=col.vec[3], cex=2)
        points(scores[group3,i], scores[group3,j], pch=21, col=col.vec[4], bg=col.vec[4], cex=2)
        points(scores[group4,i], scores[group4,j], pch=24, col=col.vec[5], bg=col.vec[5], cex=2)
        points(scores[group5,i], scores[group5,j], pch=21, col=col.vec[6], bg=col.vec[6], cex=2)
        points(scores[group6,i], scores[group6,j], pch=24, col=col.vec[7], bg=col.vec[7], cex=2)

        text(scores[group1,i], scores[group1,j], labels=group1, col=col.vec2[2], cex=0.4)
        text(scores[group2,i], scores[group2,j], labels=group2, col=col.vec2[3], cex=0.4)
        text(scores[group3,i], scores[group3,j], labels=group3, col=col.vec2[4], cex=0.4)
        text(scores[group4,i], scores[group4,j], labels=group4, col=col.vec2[5], cex=0.4)
        text(scores[group5,i], scores[group5,j], labels=group5, col=col.vec2[6], cex=0.4)
        text(scores[group6,i], scores[group6,j], labels=group6, col=col.vec2[7], cex=0.4)
        
        #legend(group.legend.pos, legend=group.vec[1:7], col=col.vec[1:7], pch=16, title="Time points", cex=0.75, pt.cex=1.5, bg="transparent", bty="n")
        ##legend(idv.legend.pos, legend=idv.vec2[1:5], col="black", pch=pch.vec[1:5], title="Donors", cex=0.75, pt.cex=1.0, pt.lwd=1.25, pt.bg="grey90", bg="transparent", bty="n")
        dev.off()
      }
    }
  }  
  
  for (i in 1:2){
    for (j in i:2){
      if (i<j){
        XLIM <- c(-max(abs(scores[,i])), max(abs(scores[,i])))
        XLIM <- XLIM+(XLIM*0.1)
        YLIM <- c(-max(abs(scores[,j])), max(abs(scores[,j])))
        YLIM <- YLIM+(YLIM*0.1)
        #print(summary(pcdat)$importance)
        Xlab <- paste0("PC", i, " (", round(summary(pcdat)$importance[2,i]*100, 2), "%)")
        Ylab <- paste0("PC", j, " (", round(summary(pcdat)$importance[2,j]*100, 2), "%)")
        png(paste(output.path, "/", filename.prefix, "02_pca_", i, "_", j, "_filteredFeatures2.png", sep=""), width=3600, height=3600, pointsize=15, res=600)
        plot(scores[group1,i], scores[group1,j], type="n", xlab=Xlab, ylab=Ylab, xlim=XLIM, ylim=YLIM, pch=16, col=col.vec[1], cex=2)
        
        points(scores[group1,i], scores[group1,j], pch=21, col=col.vec[2], bg=col.vec[2], cex=2)
        points(scores[group2,i], scores[group2,j], pch=24, col=col.vec[3], bg=col.vec[3], cex=2)
        #points(scores[group3,i], scores[group3,j], pch=21, col=col.vec[4], bg=col.vec[4], cex=2)
        #points(scores[group4,i], scores[group4,j], pch=24, col=col.vec[5], bg=col.vec[5], cex=2)
        points(scores[group5,i], scores[group5,j], pch=21, col=col.vec[6], bg=col.vec[6], cex=2)
        points(scores[group6,i], scores[group6,j], pch=24, col=col.vec[7], bg=col.vec[7], cex=2)
        
        text(scores[group1,i], scores[group1,j], labels=group1, col=col.vec2[2], cex=0.4)
        text(scores[group2,i], scores[group2,j], labels=group2, col=col.vec2[3], cex=0.4)
        #text(scores[group3,i], scores[group3,j], labels=group3, col=col.vec2[4], cex=0.4)
        #text(scores[group4,i], scores[group4,j], labels=group4, col=col.vec2[5], cex=0.4)
        text(scores[group5,i], scores[group5,j], labels=group5, col=col.vec2[6], cex=0.4)
        text(scores[group6,i], scores[group6,j], labels=group6, col=col.vec2[7], cex=0.4)
        
        #legend(group.legend.pos, legend=group.vec[1:7], col=col.vec[1:7], pch=16, title="Time points", cex=0.75, pt.cex=1.5, bg="transparent", bty="n")
        ##legend(idv.legend.pos, legend=idv.vec2[1:5], col="black", pch=pch.vec[1:5], title="Donors", cex=0.75, pt.cex=1.0, pt.lwd=1.25, pt.bg="grey90", bg="transparent", bty="n")
        dev.off()
      }
    }
  }  
  
}

#++++++++++++++++++++++ plotPCA2(...) +++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ plotPCA3(...) +++++++++++++++++++++++++++++++++++++++++++

plotPCA3 <- function(x=NULL, group.legend.pos=NULL, idv.legend.pos=NULL, filename.prefix=NULL, output.path=NULL){
  
  pca.dat <- x
  if(!is.null(output.path) && !dir.exists(paste0(output.path, "/PCA"))) {
    dir.create(paste0(output.path, "/PCA"))
  }
  output.path <- paste0(output.path, "/PCA")
  
  print(paste0("original feature number: ", nrow(pca.dat)))
  #pca.dat[is.na(pca.dat)] <- 1
  pca.dat <- na.omit(pca.dat)
  print(paste0("feature number after removing features containing NAs: ", nrow(pca.dat)))
  pca.dat <- log2(pca.dat)
  
  group1 <- grep("LEAN\\d+_0$", colnames(pca.dat), value=TRUE)
  group2 <- grep("ADIP\\d+_0$", colnames(pca.dat), value=TRUE)
  group3 <- grep("T2D\\d+_0$", colnames(pca.dat), value=TRUE)
  n1 <- length(group1)
  n2 <- length(group2)
  n3 <- length(group3)
  
  group.vec <- c(
    "lean-basal",
    "adipose-basal",
    "diabetic-basal"
  )
  
  col.vec <- c(
    adjustcolor("black", alpha=0.3),
    adjustcolor("red", alpha=0.3),
    adjustcolor("darkred", alpha=0.3),
    adjustcolor("deepskyblue4", alpha=0.3),
    adjustcolor("navy", alpha=0.3),
    adjustcolor("chartreuse3", alpha=0.3),
    adjustcolor("darkgreen", alpha=0.3),
    adjustcolor("magenta", alpha=0.3),
    adjustcolor("darkorange", alpha=0.3),
    adjustcolor("darkorchid", alpha=0.3),
    adjustcolor("darkolivegreen4", alpha=0.3),
    adjustcolor("gold", alpha=0.3),
    adjustcolor("plum", alpha=0.3),
    adjustcolor("yellow", alpha=0.3),
    adjustcolor("khaki", alpha=0.3),
    adjustcolor("seagreen", alpha=0.3),
    adjustcolor("saddlebrown", alpha=0.3),
    adjustcolor("lavender", alpha=0.3)
  )
  col.vec2 <- c(
    adjustcolor("black", alpha=0.8),
    adjustcolor("red", alpha=0.8),
    adjustcolor("darkred", alpha=0.8),
    adjustcolor("deepskyblue4", alpha=0.8),
    adjustcolor("navy", alpha=0.8),
    adjustcolor("chartreuse3", alpha=0.8),
    adjustcolor("darkgreen", alpha=0.8)
  )
  
  
  pcdat <- prcomp(t(pca.dat), center=TRUE, scale=TRUE)
  scores <- pcdat$x
  
  for (i in 1:2){
    for (j in i:2){
      if (i<j){
        XLIM <- c(-max(abs(scores[,i])), max(abs(scores[,i])))
        XLIM <- XLIM+(XLIM*0.1)
        YLIM <- c(-max(abs(scores[,j])), max(abs(scores[,j])))
        YLIM <- YLIM+(YLIM*0.1)
        #print(summary(pcdat)$importance)
        Xlab <- paste0("PC", i, " (", round(summary(pcdat)$importance[2,i]*100, 2), "%)")
        Ylab <- paste0("PC", j, " (", round(summary(pcdat)$importance[2,j]*100, 2), "%)")
        png(paste(output.path, "/", filename.prefix, "02_pca_", i, "_", j, "_filteredFeatures_BASAL.png", sep=""), width=3600, height=3600, pointsize=15, res=600)
        plot(scores[group1,i], scores[group1,j], type="n", xlab=Xlab, ylab=Ylab, xlim=XLIM, ylim=YLIM, pch=16, col=col.vec[1], cex=2)
        
        points(scores[group1,i], scores[group1,j], pch=21, col=col.vec[2], bg=col.vec[2], cex=2)
        points(scores[group2,i], scores[group2,j], pch=21, col=col.vec[3], bg=col.vec[3], cex=2)
        points(scores[group3,i], scores[group3,j], pch=21, col=col.vec[4], bg=col.vec[4], cex=2)
        
        text(scores[group1,i], scores[group1,j], labels=group1, col=col.vec2[2], cex=0.4)
        text(scores[group2,i], scores[group2,j], labels=group2, col=col.vec2[3], cex=0.4)
        text(scores[group3,i], scores[group3,j], labels=group3, col=col.vec2[4], cex=0.4)
        
        #legend(group.legend.pos, legend=group.vec[1:7], col=col.vec[1:7], pch=16, title="Time points", cex=0.75, pt.cex=1.5, bg="transparent", bty="n")
        ##legend(idv.legend.pos, legend=idv.vec2[1:5], col="black", pch=pch.vec[1:5], title="Donors", cex=0.75, pt.cex=1.0, pt.lwd=1.25, pt.bg="grey90", bg="transparent", bty="n")
        dev.off()
      }
    }
  }  
  
}

#++++++++++++++++++++++ plotPCA3(...) +++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ plotPCA4(...) +++++++++++++++++++++++++++++++++++++++++++

plotPCA4 <- function(x=NULL, group.legend.pos=NULL, idv.legend.pos=NULL, filename.prefix=NULL, output.path=NULL){
  
  pca.dat <- x
  if(!is.null(output.path) && !dir.exists(paste0(output.path, "/PCA"))) {
    dir.create(paste0(output.path, "/PCA"))
  }
  output.path <- paste0(output.path, "/PCA")
  
  print(paste0("original feature number: ", nrow(pca.dat)))
  #pca.dat[is.na(pca.dat)] <- 1
  pca.dat <- na.omit(pca.dat)
  print(paste0("feature number after removing features containing NAs: ", nrow(pca.dat)))
  pca.dat <- log2(pca.dat)
  
  group1 <- grep("LEAN\\d+_1$", colnames(pca.dat), value=TRUE)
  group2 <- grep("ADIP\\d+_1$", colnames(pca.dat), value=TRUE)
  group3 <- grep("T2D\\d+_1$", colnames(pca.dat), value=TRUE)
  n1 <- length(group1)
  n2 <- length(group2)
  n3 <- length(group3)
  
  group.vec <- c(
    "lean-clamp",
    "adipose-clamp",
    "diabetic-clamp"
  )
  
  col.vec <- c(
    adjustcolor("black", alpha=0.3),
    adjustcolor("red", alpha=0.3),
    adjustcolor("darkred", alpha=0.3),
    adjustcolor("deepskyblue4", alpha=0.3),
    adjustcolor("navy", alpha=0.3),
    adjustcolor("chartreuse3", alpha=0.3),
    adjustcolor("darkgreen", alpha=0.3),
    adjustcolor("magenta", alpha=0.3),
    adjustcolor("darkorange", alpha=0.3),
    adjustcolor("darkorchid", alpha=0.3),
    adjustcolor("darkolivegreen4", alpha=0.3),
    adjustcolor("gold", alpha=0.3),
    adjustcolor("plum", alpha=0.3),
    adjustcolor("yellow", alpha=0.3),
    adjustcolor("khaki", alpha=0.3),
    adjustcolor("seagreen", alpha=0.3),
    adjustcolor("saddlebrown", alpha=0.3),
    adjustcolor("lavender", alpha=0.3)
  )
  col.vec2 <- c(
    adjustcolor("black", alpha=0.8),
    adjustcolor("red", alpha=0.8),
    adjustcolor("darkred", alpha=0.8),
    adjustcolor("deepskyblue4", alpha=0.8),
    adjustcolor("navy", alpha=0.8),
    adjustcolor("chartreuse3", alpha=0.8),
    adjustcolor("darkgreen", alpha=0.8)
  )
  
  
  pcdat <- prcomp(t(pca.dat), center=TRUE, scale=TRUE)
  scores <- pcdat$x
  
  for (i in 1:2){
    for (j in i:2){
      if (i<j){
        XLIM <- c(-max(abs(scores[,i])), max(abs(scores[,i])))
        XLIM <- XLIM+(XLIM*0.1)
        YLIM <- c(-max(abs(scores[,j])), max(abs(scores[,j])))
        YLIM <- YLIM+(YLIM*0.1)
        #print(summary(pcdat)$importance)
        Xlab <- paste0("PC", i, " (", round(summary(pcdat)$importance[2,i]*100, 2), "%)")
        Ylab <- paste0("PC", j, " (", round(summary(pcdat)$importance[2,j]*100, 2), "%)")
        png(paste(output.path, "/", filename.prefix, "02_pca_", i, "_", j, "_filteredFeatures_CLAMP.png", sep=""), width=3600, height=3600, pointsize=15, res=600)
        plot(scores[group1,i], scores[group1,j], type="n", xlab=Xlab, ylab=Ylab, xlim=XLIM, ylim=YLIM, pch=16, col=col.vec[1], cex=2)
        
        points(scores[group1,i], scores[group1,j], pch=21, col=col.vec[2], bg=col.vec[2], cex=2)
        points(scores[group2,i], scores[group2,j], pch=21, col=col.vec[3], bg=col.vec[3], cex=2)
        points(scores[group3,i], scores[group3,j], pch=21, col=col.vec[4], bg=col.vec[4], cex=2)
        
        text(scores[group1,i], scores[group1,j], labels=group1, col=col.vec2[2], cex=0.4)
        text(scores[group2,i], scores[group2,j], labels=group2, col=col.vec2[3], cex=0.4)
        text(scores[group3,i], scores[group3,j], labels=group3, col=col.vec2[4], cex=0.4)
        
        #legend(group.legend.pos, legend=group.vec[1:7], col=col.vec[1:7], pch=16, title="Time points", cex=0.75, pt.cex=1.5, bg="transparent", bty="n")
        ##legend(idv.legend.pos, legend=idv.vec2[1:5], col="black", pch=pch.vec[1:5], title="Donors", cex=0.75, pt.cex=1.0, pt.lwd=1.25, pt.bg="grey90", bg="transparent", bty="n")
        dev.off()
      }
    }
  }  
  
}

#++++++++++++++++++++++ plotPCA4(...) +++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ plotProfiles(...) ++++++++++++++++++++++++++++++++++++++

plotProfiles <- function(x=NULL, info=NULL, prot.info=NULL, plot=TRUE, anonymize.ids=FALSE, output.path=NULL){
  
  time.pnt.vec <- c(0,1,2.5,5,15,30,60)
  # sample sequence --> "AR7", "AT13", "AT18", "AT22", "AT8"
  #pch.vec <- c(21,22,23,24,25,16)
  pch.vec <- c(21,22,23,24,25,4)
  col.vec <- c(
    adjustcolor("black", alpha=0.3),
    adjustcolor("darkgreen", alpha=0.3),
    adjustcolor("darkred", alpha=0.3),
    adjustcolor("navy", alpha=0.3),
    adjustcolor("magenta", alpha=0.3),
    "red",
    adjustcolor("deepskyblue4", alpha=0.3),
    adjustcolor("chartreuse3", alpha=0.3),
    adjustcolor("red", alpha=0.3),
    adjustcolor("darkorange", alpha=0.3),
    adjustcolor("darkorchid", alpha=0.3),
    adjustcolor("darkolivegreen4", alpha=0.3),
    adjustcolor("gold", alpha=0.3),
    adjustcolor("plum", alpha=0.3),
    adjustcolor("yellow", alpha=0.3),
    adjustcolor("khaki", alpha=0.3),
    adjustcolor("seagreen", alpha=0.3),
    adjustcolor("saddlebrown", alpha=0.3),
    adjustcolor("lavender", alpha=0.3)
  )
  
  print(paste0("original feature number: ", nrow(x)))
  x <- na.omit(x)
  print(paste0("feature number after removing features containing NAs: ", nrow(x)))
  
  pep.ids <- rownames(x)
  uniprot.ids <- gsub("\\_peptide\\d+", "", pep.ids)
  id.tab <- returnUniProtGeneNames(uniprot.ids=uniprot.ids,
                                   table=TRUE,
                                   iso.rm=TRUE,
                                   source="E:/Projekte/Phospho2.1/ID-Mapping/20220702_Phospho2.1_UniProt_to_GeneName.tsv")
  #print(id.tab[duplicated(id.tab$From),])
  #duplicate.idx <- which(duplicated(id.tab$From))
  #print(length(duplicate.idx))
  #print(nrow(id.tab))
  #id.tab <- id.tab[-duplicate.idx,]
  #print(nrow(id.tab))
  #print(id.tab[grep("O95819", id.tab$From),])
  rownames(id.tab) <- id.tab$From
  
  # ananoymization of human time course data:
  # "AR7" ---> "A"
  # "AT13" ---> "B"
  # "AT18" ---> "C"
  # "AT22" ---> "D"
  # "AT8" ---> "E
  idv.vec <- c("AR7", "AT13", "AT18", "AT22", "AT8")
  idv.vec2 <- c("A", "B", "C", "D", "E")
  group1 <- grep(idv.vec[1], colnames(x), value=T)
  group2 <- grep(idv.vec[2], colnames(x), value=T)
  group3 <- grep(idv.vec[3], colnames(x), value=T)
  group4 <- grep(idv.vec[4], colnames(x), value=T)
  group5 <- grep(idv.vec[5], colnames(x), value=T)
  groups <- rbind(group1, group2, group3, group4, group5)
  score.vec <- vector(mode="numeric", length=nrow(x))
  names(score.vec) <- rownames(x)
  cv.vec <- vector(mode="numeric", length=nrow(x))
  names(cv.vec) <- rownames(x)
  gene.name.vec <- vector(mode="numeric", length=nrow(x))
  names(gene.name.vec) <- rownames(x)
  mean.profile <- matrix(nrow=nrow(x), ncol=length(grep(idv.vec[1], colnames(x))))
  rownames(mean.profile) <- rownames(x)
  sd.profile <- matrix(nrow=nrow(x), ncol=length(grep(idv.vec[1], colnames(x))))
  rownames(sd.profile) <- rownames(x)
  cv.profile <- matrix(nrow=nrow(x), ncol=length(grep(idv.vec[1], colnames(x))))
  rownames(cv.profile) <- rownames(x)
  score.matrix <- array(NA, dim=c(nrow(x), choose(length(idv.vec),2)+1))
  rownames(score.matrix) <- rownames(x)
  for(i in 1:nrow(x)){
    pep.id <- names(score.vec)[i]
    #uniprot.id <- unlist(strsplit(pep.id, split="; "))
    #uniprot.id <- gsub("\\_peptide\\d+", "", uniprot.id)
    uniprot.id <- gsub("\\_peptide\\d+", "", pep.id)
    uniprot.id <- unique(gsub("\\-\\d+", "", x=uniprot.id))
    #gene.name.vec[i] <- paste0(unique(id.tab[uniprot.id,"To"]), collapse="; ")
    gene.name.vec[i] <- id.tab[uniprot.id,"To"]
    
    donor.profiles.tmp <- rbind(x[i, grep(idv.vec[1], colnames(x), value=TRUE)],
                            x[i, grep(idv.vec[2], colnames(x), value=TRUE)],
                            x[i, grep(idv.vec[3], colnames(x), value=TRUE)],
                            x[i, grep(idv.vec[4], colnames(x), value=TRUE)],
                            x[i, grep(idv.vec[5], colnames(x), value=TRUE)])
    #print(donor.profiles.tmp)
    mean.profile[i,] <- apply(donor.profiles.tmp, 2, mean)
    #print(mean.profile[i,])
    sd.profile[i,] <- apply(donor.profiles.tmp, 2, sd)
    #print(sd.profile[i,])
    cv.profile[i,] <- sd.profile[i,] / mean.profile[i,]
    #print(cv.profile[i,])
    
    tmp.score.vec <- vector(mode="numeric", length=choose(length(idv.vec),2))
    count <- 1
    for(j in 1:(length(idv.vec)-1)){
      idx1 <- grep(idv.vec[j], colnames(x), value=TRUE)
      for(k in (j+1):length(idv.vec)){
        idx2 <- grep(idv.vec[k], colnames(x), value=TRUE)
        
        tmp.score.vec[count] <- cor(log2(x[i,idx1]), log2(x[i,idx2]), method="pearson")
        #tmp.score.vec[count] <- cor(log2(x[i,idx1]), log2(x[i,idx2]), method="spearman") 
        count <- count+1
      }
    }
    score.vec[i] <- mean(tmp.score.vec)
    #score.vec[i] <- median(tmp.score.vec)
    cv.vec[i] <- mean(cv.profile[i,])
    
    score.matrix[i,] <- c(tmp.score.vec, score.vec[i])
  }
  output.score.matrix <- cbind(rownames(score.matrix), score.matrix)
  colnames(output.score.matrix) <- c("Peptide",
                                     "AR7_vs_AT13",
                                     "AR7_vs_AT18",
                                     "AR7_vs_AT22",
                                     "AR7_vs_AT8",
                                     "AT13_vs_AT18",
                                     "AT13_vs_AT22",
                                     "AT13_vs_AT8",
                                     "AT18_vs_AT22",
                                     "AT18_vs_AT8",
                                     "AT22_vs_AT8",
                                     "Mean")
  write.table(x=output.score.matrix,
              file=paste0(output.path, "/score.matrix.txt"),
              row.names=F,
              col.names=T,
              sep="\t")
  
  output <- cbind(rep(0,length(score.vec)),
                  names(score.vec),
                  info[names(score.vec),],
                  gene.name.vec,
                  prot.info[uniprot.ids,c("UniProt ID", "Description", "Gene Symbol")],
                  score.vec,
                  cv.vec)
  colnames(output) <- c("Pearsons r rank",
                        "Peptide",
                        colnames(info),
                        "Gene",
                        "UniProt ID",
                        "Description",
                        "Gene Symbol",
                        "Pearsons r",
                        "CV")
  rownames(output) <- names(score.vec)
  output <- output[order(as.numeric(output[,"Pearsons r"]), decreasing=TRUE),]
  score.vec <- sort(score.vec, decreasing=TRUE)
  cv.vec <- cv.vec[names(score.vec)]
  output[,"Pearsons r rank"] <- 1:length(score.vec) 
  #print(output)
  write.table(x=output,
              file=paste0(output.path, "/profile-scores.txt"),
              row.names=FALSE, col.names=TRUE, sep="\t")
  
  #-------------------- Begin: Plot the profiles?! -----------------------------
  if(plot==TRUE){
    if(anonymize.ids == TRUE){
      idv.vec3 <- idv.vec2  
    }else{
      idv.vec3 <- idv.vec  
    }
    for(i in 1:nrow(x)){
      pep.id <- names(score.vec)[i]
      #uniprot.id <- unlist(strsplit(pep.id, split="; "))
      #uniprot.id <- gsub("\\_peptide\\d+", "", uniprot.id)
      uniprot.id <- gsub("\\_peptide\\d+", "", pep.id)
      uniprot.id <- unique(gsub("\\-\\d+", "", x=uniprot.id))
      print(paste0("uniprot.id: ", uniprot.id))
      print(paste0("gene name: ", unique(id.tab[uniprot.id, "To"])))
      print(paste0("pep.id: ", pep.id))
      print(info[pep.id,"Modifications in Master Proteins"])
      if(length(unique(id.tab[uniprot.id, "To"])) > 1){
        gene.name.tmp <- paste(x=unique(id.tab[uniprot.id, "To"]), collapse=", ")
      }else{
        gene.name.tmp <- unique(id.tab[uniprot.id, "To"])
      }
      mods.tmp <- gsub(".*(\\[.*\\]).*", "\\1", info[pep.id,"Modifications in Master Proteins"])
      if(is.na(mods.tmp)) mods.tmp <- gsub(".*xPhospho (\\[.*\\]).*", "\\1", info[pep.id,"Modifications"])
      title <- paste0(x=uniprot.id,
                    " / ",
                    gene.name.tmp,
                    " ",
                    mods.tmp,
                    "\n avg. r = ",
                    round(score.vec[i], 4),
                    ", avg. CV = ",
                    round(cv.vec[i], 4)
                )
      print(title)
      print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
      ymin <- min(log10(x[pep.id,]))
      ymax <- max(log10(x[pep.id,]))
      ymax <- ymax + ymax*0.02
      png(filename=paste0(output.path, "/profile_", i, "_", gsub(", ","_",gene.name.tmp), "_", pep.id, ".png"), height=2000, width=2000, res=300)
        plot(x=time.pnt.vec, y=log10(x[pep.id, groups[1,]]), col=col.vec[1], bg=col.vec[1], type="b", lty=1, lwd=3, pch=pch.vec[1], cex=2, xaxt="n", main=title, xlab="Time (min)", ylab="log10(abundance)", ylim=c(ymin,ymax))
        points(x=time.pnt.vec, y=log10(x[pep.id, groups[2,]]), col=col.vec[2], bg=col.vec[2], type="b", lty=1, lwd=3, pch=pch.vec[2], cex=2)
        points(x=time.pnt.vec, y=log10(x[pep.id, groups[3,]]), col=col.vec[3], bg=col.vec[3], type="b", lty=1, lwd=3, pch=pch.vec[3], cex=2)
        points(x=time.pnt.vec, y=log10(x[pep.id, groups[4,]]), col=col.vec[4], bg=col.vec[4], type="b", lty=1, lwd=3, pch=pch.vec[4], cex=2)
        points(x=time.pnt.vec, y=log10(x[pep.id, groups[5,]]), col=col.vec[5], bg=col.vec[5], type="b", lty=1, lwd=3, pch=pch.vec[5], cex=2)
        points(x=time.pnt.vec, y=log10(mean.profile[pep.id,]), col="red", bg="red", type="b", lty=1, lwd=3, pch=pch.vec[6], cex=2)
        axis(1, at=time.pnt.vec, labels=c("0","1","2.5","5","15","30","60"), col.axis="black", las=2, cex.axis=1.0, tck=-0.01)
        legend("top",
             legend=c(idv.vec3, "mean"),
             col=col.vec[1:6],
             pch=pch.vec,
             lty=1,
             lwd=3,
             seg.len=3,
             pt.lwd=3,
             cex=0.75,
             pt.cex=1.5,
             pt.bg=col.vec[1:6],
             x.intersp=0.25,
             text.width = 5,
             text.font = 2,
             bg="transparent",
             bty="n",
             ncol=6)
      dev.off()
    }
  
  }
  #-------------------- End: Plot the profiles?! -----------------------------
  
  sel.n <- 33
  selected <- rev(rownames(output)[1:sel.n])
  
  # REMOVING KHDRBS1 PROTEOFORM WITH LOWER PEARSON'S r
  #selected <- selected[-grep("Q07666_peptide2", selected)]
  #sel.n <- sel.n-1
  
  ylab1 <- output[selected,"Gene"]
  ylab1 <- gsub("AKT2; AKT1; AKT3","AKT1/2/3",ylab1)
  ylab1 <- gsub("WASHC2C; WASHC2A","WASHC2C/A",ylab1)
  ylab1 <- gsub("ARHGEF6; ARHGEF7","ARHGEF6/7",ylab1)
  
  ylab2 <- gsub("\\[N-Term\\]", "", info[selected,"Modifications in Master Proteins"])
  ylab2 <- gsub("\\].*\\[", "/", ylab2)
  ylab2 <- gsub("\\(\\d+\\.*\\d*\\)", "", ylab2)
  ylab2 <- gsub("\\]", "", ylab2)
  ylab2 <- gsub("^.*\\[", "", ylab2)
  
  #FOR SHOWING ONLY FIRST SITE OF PEPTIDE GROUP:
  ylab2 <- gsub("([STY]\\d+).*", "\\1", ylab2)
  
  #png(filename=paste0(output.path, "/stripchart_donor_similarity_32.png"), height=2000, width=2000, res=300)
  png(filename=paste0(output.path, "/stripchart_donor_similarity.png"), height=2000, width=2000, res=300)
    mar.default <- c(5.1, 4.1, 4.1, 2.1)
    par(mar = mar.default + c(0, 6, 0, 6))
    stripchart(as.data.frame(t(score.matrix[selected,1:10])), axes=F, pch=21, cex=1.75, main="Pairwise similarity of donor profiles", xlab="Pearson's r")
    axis(side=1)
    axis(side=2, at=1:sel.n, labels=ylab1, las=1, cex.axis=0.7)
    axis(side=4, at=1:sel.n, labels=ylab2, las=1, cex.axis=0.7)
    abline(h=1:sel.n,lty=3,col="grey",lwd=2)
    stripchart(as.data.frame(t(score.matrix[selected,1:10])), pch=21, cex=1.75, add=T)
    stripchart(as.data.frame(t(score.matrix[selected,11])), pch=21, col="red", bg="red", cex=1.75, add=T) 
  dev.off()
}

#++++++++++++++++++++++ plotProfiles(...) ++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ plotVolcano(...) +++++++++++++++++++++++++++++++++++++++

plotVolcano <- function(fc=NULL,
                        p.val=NULL,
                        reg=NULL,
                        mark=NULL,
                        mark.txt=NULL,
                        move.txt.x=NULL,
                        move.txt.y=NULL,
                        col.idx=2,
                        mark2=NULL,
                        col.idx2=NULL,
                        main=NULL,
                        filename.tag=NULL,
                        output.path=NULL){
  
  if(!is.null(output.path) && !dir.exists(paste0(output.path, "/Volcanos"))) {
    dir.create(paste0(output.path, "/Volcanos"))
  }
  output.path <- paste0(output.path, "/Volcanos")
  
  
  col.vec <- c(
    adjustcolor("black", alpha=0.2),
    adjustcolor("red", alpha=0.8),
    adjustcolor("darkred", alpha=0.8),
    adjustcolor("deepskyblue4", alpha=0.8),
    adjustcolor("navy", alpha=0.8),
    adjustcolor("chartreuse3", alpha=0.8),
    adjustcolor("darkgreen", alpha=0.8),
    adjustcolor("magenta", alpha=0.8),
    adjustcolor("darkorange", alpha=0.8),
    adjustcolor("darkorchid", alpha=0.8),
    adjustcolor("darkolivegreen4", alpha=0.8),
    adjustcolor("gold", alpha=0.8),
    adjustcolor("plum", alpha=0.8),
    adjustcolor("yellow", alpha=0.8),
    adjustcolor("khaki", alpha=0.8),
    adjustcolor("seagreen", alpha=0.8),
    adjustcolor("saddlebrown", alpha=0.8),
    adjustcolor("lavender", alpha=0.8),
    adjustcolor("grey", alpha=0.75)
  )
  
  if(is.null(move.txt.x)){
    move.txt.x <- rep(0,length(mark.txt))
  }
  if(is.null(move.txt.y)){
    move.txt.y <- rep(0,length(mark.txt))
  }
  
  png(paste0(output.path, "/01_volcano_large", filename.tag, ".png"), 
      width=3600,
      height=3600,
      pointsize=15, res=600)
    #mar.default <- c(5.1, 4.1, 4.1, 2.1)
    #par(mar = mar.default + c(0, 0.5, 0, 0))
    mgp.default <- c(3, 1, 0)
    par(mgp = mgp.default + c(-0.25, 0, 0))
    plot(log2(fc),
       -log10(p.val),
       col=col.vec[19],
       #pch=16,
       cex=1.5,
       cex.lab=1.5,
       cex.axis=1.5, #---> tick-labels
       cex.main=1.75,
       #font.axis=2, ---> tick-labels
       main=main,
       ylab="-log10(p-values)",
       xlab="log2(fold change)")
    abline(h=-log10(0.05), lty=3, lwd=4, col="grey50")
    abline(v=log2(1.5), lty=3, lwd=4, col="grey50")
    abline(v=log2(2/3), lty=3, lwd=4, col="grey50")
    points(log2(fc[reg]),
         -log10(p.val[reg]),
         col=col.vec[1],
         pch=16,
         cex=1.5)
    points(log2(fc[mark]),
         -log10(p.val[mark]),
         col=col.vec[col.idx],
         pch=16,
         cex=1.5)
  dev.off() 
  
  
  if(!is.null(mark.txt)){
    png(paste0(output.path, "/02_volcano_large", filename.tag, ".png"), 
        width=3600,
        height=3600,
        pointsize=15, res=600)
      #mar.default <- c(5.1, 4.1, 4.1, 2.1)
      #par(mar = mar.default + c(0, 0.5, 0, 0))
      mgp.default <- c(3, 1, 0)
      par(mgp = mgp.default + c(-0.25, 0, 0))
      plot(log2(fc),
          -log10(p.val),
          col=col.vec[19],
          #pch=16,
          cex=1.5,
          cex.lab=1.5,
          cex.axis=1.5, #---> tick-labels
          cex.main=1.75,
          #font.axis=2, ---> tick-labels
          main=main,
          ylab="-log10(p-values)",
          xlab="log2(fold change)")
      abline(h=-log10(0.05), lty=3, lwd=4, col="grey50")
      abline(v=log2(1.5), lty=3, lwd=4, col="grey50")
      abline(v=log2(2/3), lty=3, lwd=4, col="grey50")
      points(log2(fc[reg]),
          -log10(p.val[reg]),
          col=col.vec[1],
          pch=16,
          cex=1.5)
      points(log2(fc[mark]),
          -log10(p.val[mark]),
          col=col.vec[col.idx],
          pch=16,
          cex=1.5)
      text(log2(fc[mark])+move.txt.x,
          -log10(p.val[mark])+move.txt.y,
          labels=mark.txt,
          pos=1,
          col=col.vec[col.idx],
          cex=0.6,
          font=2)
    dev.off()
  }
  
  
  png(paste0(output.path, "/03_volcano_small", filename.tag, ".png"), 
      width=3600,
      height=3600,
      pointsize=15, res=600)
    #mar.default <- c(5.1, 4.1, 4.1, 2.1)
    #par(mar = mar.default + c(0, 0.5, 0, 0))
    mgp.default <- c(3, 1, 0)
    par(mgp = mgp.default + c(-0.25, 0, 0))
    plot(log2(fc),
       -log10(p.val),
       col=col.vec[19],
       #pch=16,
       cex=2.5,
       cex.lab=1.8,
       cex.axis=1.8, #---> tick-labels
       cex.main=1.8,
       #font.axis=2, ---> tick-labels
       main=main,
       ylab="-log10(p-values)",
       xlab="log2(fold change)")
    abline(h=-log10(0.05), lty=3, lwd=4, col="grey50")
    abline(v=log2(1.5), lty=3, lwd=4, col="grey50")
    abline(v=log2(2/3), lty=3, lwd=4, col="grey50")
    points(log2(fc[reg]),
         -log10(p.val[reg]),
         col=col.vec[1],
         pch=16,
         cex=2.5)
    points(log2(fc[mark]),
         -log10(p.val[mark]),
         col=col.vec[col.idx],
         pch=16,
         cex=2.5)
    if(!is.null(mark2)){
      points(log2(fc[mark2]),
             -log10(p.val[mark2]),
             col=col.vec[col.idx2],
             pch=16,
             cex=2.5)  
    }
  dev.off() 
}

#++++++++++++++++++++++ plotVolcano(...) +++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ plotVolcanoCrowd(...) ++++++++++++++++++++++++++++++++++

plotVolcanoCrowd <- function(fc=NULL, p.val=NULL, filename.prefix=NULL, output.path=NULL){
  
  library(scatterplot3d)
  library(wesanderson)
  
  
  if(!is.null(output.path) && !dir.exists(paste0(output.path, "/VolcanoCrowd"))) {
    dir.create(paste0(output.path, "/VolcanoCrowd"))
  }
  output.path <- paste0(output.path, "/VolcanoCrowd")
  
  n <- nrow(fc)
  
  dat.3d.tmp1 <- cbind(rep(1,n), log2(fc[,1]), -log10(p.val[,1]))
  dat.3d.tmp2 <- cbind(rep(2,n), log2(fc[,2]), -log10(p.val[,2]))
  dat.3d.tmp3 <- cbind(rep(3,n), log2(fc[,3]), -log10(p.val[,3]))
  dat.3d.tmp4 <- cbind(rep(4,n), log2(fc[,4]), -log10(p.val[,4]))
  dat.3d.tmp5 <- cbind(rep(5,n), log2(fc[,5]), -log10(p.val[,5]))
  dat.3d.tmp6 <- cbind(rep(6,n), log2(fc[,6]), -log10(p.val[,6]))
  
  dat.3d <- rbind(dat.3d.tmp1, dat.3d.tmp2, dat.3d.tmp3, dat.3d.tmp4, dat.3d.tmp5, dat.3d.tmp6)
  
  col1 <- rep(adjustcolor("red", alpha=0.3), n)
  col2 <- rep(adjustcolor("darkred", alpha=0.3), n)
  col3 <- rep(adjustcolor("deepskyblue4", alpha=0.3), n)
  col4 <- rep(adjustcolor("navy", alpha=0.3), n)
  col5 <- rep(adjustcolor("chartreuse3", alpha=0.3), n)
  col6 <- rep(adjustcolor("darkgreen", alpha=0.3), n)
  #col1 <- rep(adjustcolor("grey80", alpha=1.0), n)
  #col2 <- rep(adjustcolor("grey70", alpha=1.0), n)
  #col3 <- rep(adjustcolor("grey60", alpha=1.0), n)
  #col4 <- rep(adjustcolor("grey45", alpha=1.0), n)
  #col5 <- rep(adjustcolor("grey30", alpha=1.0), n)
  #col6 <- rep(adjustcolor("grey0", alpha=1.0), n)
  
  #pal <- wes_palette(name="BottleRocket1", n=6, type="continuous")
  #pal <- wes_palette(name="BottleRocket2", n=6, type="continuous")
  #pal <- wes_palette(name="Rushmore", n=6, type="continuous")
  #pal <- wes_palette(name="Royal1", n=6, type="continuous")
  #pal <- wes_palette(name="Royal2", n=6, type="continuous")
  #pal <- wes_palette(name="Darjeeling1", n=6, type="continuous")
  #pal <- wes_palette(name="Darjeeling2", n=6, type="continuous")
  #pal <- wes_palette(name="Chevalier1", n=6, type="continuous")
  #pal <- wes_palette(name="FantasticFox1", n=6, type="continuous")
  #pal <- wes_palette(name="GrandBudapest1", n=6, type="continuous")
  #pal <- wes_palette(name="GrandBudapest2", n=6, type="continuous")
  #pal <- wes_palette(name="Zissou1", n=6, type="continuous")
  #pal <- wes_palette(name="Cavalcanti1", n=6, type="continuous")
  #pal <- wes_palette(name="Moonrise1", n=6, type="continuous")
  #pal <- wes_palette(name="Moonrise2", n=6, type="continuous")
  #pal <- wes_palette(name="Moonrise3", n=6, type="continuous")
  #col1 <- rep(pal[1], n)
  #col2 <- rep(pal[2], n)
  #col3 <- rep(pal[3], n)
  #col4 <- rep(pal[4], n)
  #col5 <- rep(pal[5], n)
  #col6 <- rep(pal[6], n)
  col <- c(col1, col2, col3, col4, col5, col6)
  
  time.vec <- c(
    "1 min",
    "2.5 min",
    "5 min",
    "15 min",
    "30 min",
    "60 min"
  )
  
  png(paste(output.path, "/", filename.prefix, "volcano_crowd.png", sep=""), width=7000, height=7000, pointsize=15, res=1200)
    scatterplot3d(
      x=dat.3d[,1:3],
      box=F,
      xlab="Time points",
      ylab="log2(fold change)",
      zlab="-log10(p-value)",
      scale.y = 0.6,
      ylim=c(-6.7,6.7),
      #x.ticklabs = time.vec,
      x.ticklabs = "",
      y.margin.add=0.25,
      color=col,
      bg=col,
      #pch=1,
      pch=21,
      angle=45
    )
    text(x=1:6,
         y=-0.1,
         #pos=1,
         labels=time.vec,
         srt=60,
         adj=1,
         xpd=TRUE,
         offset=5,
         cex=0.8)
  dev.off()
  
}

#++++++++++++++++++++++ plotVolcanoCrowd(...) ++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ prepareClueRDat(...) ++++++++++++++++++++++++++++++++++++

prepareClueRDat <- function(sig.reg.pept=NULL,
                               dat.ratio.used=NULL,
                               mods.master.used=NULL,
                               psm.number.used=NULL,
                               clue.pepforms=TRUE,
                               organism=NULL,
                               mode="kinase",
                               KSorPSanno=NULL,
                               log=TRUE,
                               ratios=NULL,
                               output.path=NULL,
                               file.name="clueR.dat"){
  
  if(mode=="kinase"){
    #devtools::install_github("PengyiYang/ClueR")
    #library(ClueR)
    #library(PhosR)
    #data(PhosphoSite) #PhosphoSitePlus-ClueR-Version: 1) overlap phos-sites with MS-data: 233 2) kinases: 206 3) substrates: 9830
    #data(PhosphoSitePlus) #PhosphoSitePlus-PhosR-Version: 1) overlap phos-sites with MS-data: 275 2) kinases: 379 3) substrates: 11998
    PhosphoSite.human <- preparePhosSitePlus(organism=organism) #own PhosphoSitePlus-Version: 1) overlap with MS-data: 319 2) kinases: 403 3) substrates: 13664
    #PhosphoSite.mouse <- preparePhosSitePlus(organism=organism)
    mode.folder <- "/ClueR/ClueR-K" 
  }else if(mode=="phosphatase"){
    PSanno <- KSorPSanno
    mode.folder <- "/ClueR/ClueR-P" 
  }
  
  
  dat.ratio.used.sig <- dat.ratio.used[sig.reg.pept,]
  mods.master.used.sig <- mods.master.used[sig.reg.pept]
  psm.number.used.sig <- psm.number.used[sig.reg.pept]
  
  clueR.dat.tmp <- matrix(ncol=ncol(dat.ratio.used.sig))
  psm.number.used.sig.clueR.tmp <- c()
  #test.idx <- c(85:100, 205:215)
  #test.idx <- 1100:1200
  #for(i in test.idx){
  for(i in 1:length(mods.master.used.sig)){
    print(mods.master.used.sig[i])
    if(clue.pepforms == TRUE && (length(grep("]; ", mods.master.used.sig[i])) > 0)){
      splitted.mods <- unlist(strsplit(mods.master.used.sig[i], "]; "))
      splitted.mods[1:(length(splitted.mods)-1)] <- paste0(splitted.mods[1:(length(splitted.mods)-1)], "]")
      
      last.uniprot.id <- ""
      for(j in 1:length(splitted.mods)){
        if(length(grep("^\\d+x", splitted.mods[j]) > 0)){
          splitted.mods[j] <- paste0(last.uniprot.id, " ", splitted.mods[j])
        }else{
          last.uniprot.id <- gsub("^([0-9A-Z]+) \\d+x.+", "\\1", splitted.mods[j])
        }  
      }
      
      for(j in 1:length(splitted.mods)){
        clueR.dat.tmp <- rbind(clueR.dat.tmp, dat.ratio.used.sig[i,]) 
        rownames(clueR.dat.tmp)[nrow(clueR.dat.tmp)] <- splitted.mods[j]
        print(rownames(clueR.dat.tmp)[nrow(clueR.dat.tmp)])
        psm.number.used.sig.clueR.tmp <- c(psm.number.used.sig.clueR.tmp, psm.number.used.sig[i])
        names(psm.number.used.sig.clueR.tmp)[length(psm.number.used.sig.clueR.tmp)] <- splitted.mods[j]
      }
    }else if(length(grep("]; ", mods.master.used.sig[i])) == 0 & !is.na(mods.master.used.sig[i])){
      clueR.dat.tmp <- rbind(clueR.dat.tmp, dat.ratio.used.sig[i,]) 
      rownames(clueR.dat.tmp)[nrow(clueR.dat.tmp)] <- mods.master.used.sig[i]
      print(rownames(clueR.dat.tmp)[nrow(clueR.dat.tmp)])
      psm.number.used.sig.clueR.tmp <- c(psm.number.used.sig.clueR.tmp, psm.number.used.sig[i])
      names(psm.number.used.sig.clueR.tmp)[length(psm.number.used.sig.clueR.tmp)] <- mods.master.used.sig[i]
    }
    print("****************************************************************")
  }
  if(sum(is.na(clueR.dat.tmp[1,])) == ncol(clueR.dat.tmp)) clueR.dat.tmp <- clueR.dat.tmp[-1,] # <---------- ????
  #print(clueR.dat.tmp)
  
  clueR.dat <- matrix(ncol=ncol(clueR.dat.tmp))
  psm.number.used.sig.clueR <- c()
  for(i in 1:nrow(clueR.dat.tmp)){
    print(paste0("original: ", rownames(clueR.dat.tmp)[i]))
    if(clue.pepforms == TRUE && (length(grep("); ", rownames(clueR.dat.tmp)[i])) > 0)){
      uniprot.id <- gsub("^([0-9A-Z]+)-*\\d* .*", "\\1", rownames(clueR.dat.tmp)[i], perl=TRUE)
      tmp <- gsub("^.* \\d+xPhospho ", "", rownames(clueR.dat.tmp)[i], perl=TRUE)
      tmp <- gsub("\\(\\d+\\.*\\d*\\)", "", tmp)
      tmp <- gsub("(\\[|\\])", "", tmp)
      splitted.sites <- unlist(strsplit(tmp, "; "))
      for(j in 1:length(splitted.sites)){
        print(paste0(uniprot.id, ":", splitted.sites[j]))
        clueR.dat <- rbind(clueR.dat, clueR.dat.tmp[i,])
        rownames(clueR.dat)[nrow(clueR.dat)] <- paste0(uniprot.id, ":", splitted.sites[j]) 
        print(paste0("PSM number: ", psm.number.used.sig.clueR.tmp[i]))
        psm.number.used.sig.clueR <- c(psm.number.used.sig.clueR, psm.number.used.sig.clueR.tmp[i])
        names(psm.number.used.sig.clueR)[length(psm.number.used.sig.clueR)] <- paste0(uniprot.id, ":", splitted.sites[j]) 
      }
    }else if(is.na(rownames(clueR.dat.tmp)[i])){
      print(paste0("NA rowname found for row: ", i)) 
      print(rownames(clueR.dat.tmp)[i])
      print(paste0("PSM number: ", psm.number.used.sig.clueR.tmp[i]))
    }else if(length(grep("); ", rownames(clueR.dat.tmp)[i])) == 0){
      uniprot.id <- gsub("^([0-9A-Z]+)-*\\d* .*", "\\1", rownames(clueR.dat.tmp)[i], perl=TRUE)
      site <- gsub(".*\\[([A-Z]{1}\\d+)\\(\\d+\\.*\\d*\\)\\]$", "\\1", rownames(clueR.dat.tmp)[i], perl=TRUE)
      print(paste0(uniprot.id, ":", site))
      clueR.dat <- rbind(clueR.dat, clueR.dat.tmp[i,])
      rownames(clueR.dat)[nrow(clueR.dat)] <- paste0(uniprot.id, ":", site)
      print(paste0("PSM number: ", psm.number.used.sig.clueR.tmp[i]))
      psm.number.used.sig.clueR <- c(psm.number.used.sig.clueR, psm.number.used.sig.clueR.tmp[i])
      names(psm.number.used.sig.clueR)[length(psm.number.used.sig.clueR)] <- paste0(uniprot.id, ":", site) 
    }
    print("------------------------------------------------------")
  }
  if(sum(is.na(clueR.dat[1,])) == ncol(clueR.dat)) clueR.dat <- clueR.dat[-1,]
  if(length(grep("\\[",rownames(clueR.dat))) > 0) clueR.dat <- clueR.dat[-grep("\\[",rownames(clueR.dat)), ]
  if(length(grep("/",rownames(clueR.dat))) > 0) clueR.dat <- clueR.dat[-grep("/",rownames(clueR.dat)), ]
  if(ratios==TRUE & anyNA(clueR.dat[,1])) clueR.dat <- clueR.dat[-which(is.na(clueR.dat[,1])), ]

  
  
  del.idx <- c()
  unique.rows <- unique(rownames(clueR.dat))
  for(i in 1:length(unique.rows)){
    if(length(grep(unique.rows[i], rownames(clueR.dat))) > 1){
      tmp.idx.vec <- grep(unique.rows[i], rownames(clueR.dat))
      print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
      print("tmp.idx.vec:")
      print(tmp.idx.vec)
      print("PSM-Numbers:")
      print(psm.number.used.sig.clueR[tmp.idx.vec])
      keep.idx <- which.max(psm.number.used.sig.clueR[tmp.idx.vec])
      print("keep.idx:")
      print(keep.idx)
      tmp.idx.vec <- tmp.idx.vec[-keep.idx]
      print("tmp.idx.vec:")
      print(tmp.idx.vec)
      del.idx <- c(del.idx, tmp.idx.vec)
    }
  }
  clueR.dat <- clueR.dat[-del.idx,]
  clueR.dat <- clueR.dat[!duplicated(clueR.dat),]
  
  
  
  id.vec <- gsub("\\:[A-Z]{1}\\d+$", "", rownames(clueR.dat))
  site.vec <- gsub("^[0-9A-Z]+\\:", "", rownames(clueR.dat))
  id.mapping <- returnUniProtGeneNames(uniprot.ids=id.vec,
                                       table=TRUE,
                                       source="E:/Projekte/Phospho2.1/ID-Mapping/20220702_Phospho2.1_UniProt_to_GeneName.tsv")
  #print("id.mapping: ")
  #print(id.mapping)
  for(i in 1:nrow(id.mapping)){
    id.vec[which(id.vec == id.mapping[i,"From"])] <- id.mapping[i,"To"]
    #id.vec[which(id.vec == id.mapping[i,1])] <- id.mapping[i,2]
  }
  rownames(clueR.dat) <- paste0(id.vec,";",site.vec,";")
  
  clueR.dat <- clueR.dat[!duplicated(rownames(clueR.dat)),]
  
  if(ratios==TRUE){
    baseline.col <- rep(1, nrow(clueR.dat))
    clueR.dat <- cbind("Abundance Ratio: (0) / (0)"=baseline.col, clueR.dat)
  }
  if(log==TRUE) clueR.dat <- log2(clueR.dat)
  colnames(clueR.dat) <- c("0min", "1min", "2.5min", "5min", "15min", "30min", "60min")
  
  if(organism == "human" && mode == "kinase") {
    phos.site.overlap <- length(intersect(unlist(PhosphoSite.human), rownames(clueR.dat)))
  }else if(organism == "mouse" && mode == "kinase"){
    phos.site.overlap <- length(intersect(unlist(PhosphoSite.mouse), rownames(clueR.dat)))  
  }else if(organism == "human" && mode == "phosphatase"){
    phos.site.overlap <- length(intersect(unlist(PSanno), rownames(clueR.dat)))
  }else stop("prepareClueRDat(...): Error - unknown organism and/or mode!!!")
  print(paste0("Overlap of measured phos-sites with database: ", phos.site.overlap))
  if(is.null(output.path) == FALSE){
      write.table(cbind(rownames(clueR.dat),clueR.dat), file=paste0(output.path, "/ClueR/", file.name, ".txt"), sep="\t", col.names = TRUE, row.names = FALSE)
      save(clueR.dat, file=paste0(output.path, mode.folder, "/", file.name, ".RData"), compress = "xz", compression_level = 9)
  }
  
  return(clueR.dat)
}

#++++++++++++++++++++++ prepareClueRDat(...) ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ preparePeptideInfo(...) ++++++++++++++++++++++++++++++++++++

preparePeptideInfo <- function(pep.dat=NULL, pep.cols1=NULL, prot.info=NULL, pep.cols2=NULL, output.path=NULL){
  
  #rownames(pep.cols1) <- rownames(pep.dat)
  prot.ids <- gsub("^([0-9A-Z]+-*\\d*)_peptide\\d+$", "\\1", rownames(pep.dat), perl=TRUE)
  
  output <- cbind(pep.cols1, prot.info[prot.ids,], pep.dat, pep.cols2)
  output <- cbind(rownames(output), output)
  colnames(output)[1] <- "Peptide ID"
  colnames(output)[2] <- "Sequence"
  colnames(output)[3] <- "Modifications" 
  
  if(!is.null(output.path)) write.table(x=output, file=paste0(output.path, "/phosphopeptides_info.txt"), sep="\t", row.names=FALSE, col.names=TRUE)
  
  return(output)
}

#++++++++++++++++++++++ preparePeptideInfo(...) ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ preparePhosSitePlus(...) +++++++++++++++++++++++++++++++

preparePhosSitePlus <- function(organism=NULL){
  
  if(organism=="human"){
    phosSitePlus.dat <- read.table(file="E:/Projekte/Phospho2.1-pd3.0/output_phospho2.1-workflow/ClueR/Kinase_Substrate_Dataset_2022-09-09_human.txt", header=T,sep="\t")
  }else if(organism=="mouse"){
    phosSitePlus.dat <- read.table(file="E:/Projekte/Phospho2.1-pd3.0/output_phospho2.1-workflow/ClueR/Kinase_Substrate_Dataset_2022-09-09_mouse.txt", header=T,sep="\t")
  }else{
    print("preparePhosSitePlus: Error - organism not defined!")
    stop()
  }
  
  n.kinases <- length(unique(phosSitePlus.dat$GENE))
  print(paste0("preparePhosSitePlus - organism: ", organism))
  print(paste0("preparePhosSitePlus - PhosphoSitePlus: ", n.kinases, " kinases"))
  print(paste0("preparePhosSitePlus - PhosphoSitePlus: ", nrow(phosSitePlus.dat), " phosphosites"))
  output <- vector(mode="list", length=n.kinases)
  names(output) <- sort(unique(phosSitePlus.dat$GENE)) 
  
  for(i in 1:nrow(phosSitePlus.dat)){
    kinase <- phosSitePlus.dat[i,"GENE"]  
    substrate <- phosSitePlus.dat[i,"SUB_GENE"] 
    phossite <- phosSitePlus.dat[i,"SUB_MOD_RSD"]
    
    output[[kinase]] <- c(output[[kinase]], paste0(substrate,";",phossite,";", collapse=""))
    #print(paste0(substrate,";",phossite,";", collapse=""))
  }
  
  return(output)
}

#++++++++++++++++++++++ preparePhosSitePlus(...) +++++++++++++++++++++++++++++++

#++++++++++++++++++++++ prepareProteinInfo(...) ++++++++++++++++++++++++++++++++++++

#Own PD p-values
prepareProteinInfo <- function(pep.names=NULL, prot.info=NULL){
  # pep.names: peptide names
  # prot.info: vector providing names of ratio columns
  
  col.idx <- c("Accession",
               "Gene Symbol",
               "Entrez Gene ID",
               "Description",
               "Reactome Pathways",
               "Master",
               "Contaminant",
               "# Protein Groups")
  prot.info <- prot.info[grep("Master Protein$", prot.info[,"Master"]), ]
  
  # remove peptide number substring to get protein accessions
  print(paste0("length(pep.names): ", length(pep.names)))
  prot.accessions <- gsub("\\_peptide\\d+", "", pep.names)
  print(paste0("length(unique(prot.accessions)): ", length(unique(prot.accessions))))
  
  # get protein group indices within vector of protein accessions, save them
  # in prot.groups vector and split protein groups to have complete list
  # containing separate accessions
  prot.groups.idx <- grep("\\; ", prot.accessions)
  print(paste0("number of protein groups: ", length(prot.groups.idx)))
  prot.groups <- grep("\\; ", prot.accessions, value=TRUE)
  splitted.accessions <- unlist(strsplit(grep("\\; ", prot.groups, value=TRUE), "\\; "))
  print(paste0("number of splitted accessions: ", length(splitted.accessions)))
  
  # replace protein groups by list of separate accessions 
  if(length(prot.groups.idx) > 0){
    prot.accessions <- prot.accessions[-prot.groups.idx]
    prot.accessions <- c(prot.accessions, splitted.accessions)
    print(paste0("prot.accessions: ", length(prot.accessions)))
    print(paste0("unique(prot.accessions): ", length(unique(prot.accessions))))
  }
  
  # finalize output
  output <- prot.info[prot.info[,"Accession"] %in% prot.accessions, col.idx]
  rownames(output) <- output[,"Accession"]
  #colnames(output) <- c("UniProt ID", "Description", "Gene Symbol", "Master", "Contaminant")
  colnames(output) <- c("UniProt ID",
                        "Gene Symbol",
                        "Entrez Gene ID",
                        "Description",
                        "Reactome Pathways",
                        "Master",
                        "Contaminant",
                        "# Protein Groups")
  print(paste0("nrow(output): ", nrow(output)))
  
  
  return(output)
}

#++++++++++++++++++++++ prepareProteinInfo(...) ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ regulatedPeptides(...) +++++++++++++++++++++++++++++++++

regulatedPeptides <- function(ratios=NULL,
                              p.values=NULL,
                              direction="positive",  #--> "positive", "negative", "both"
                              r.thresh.pos=1.5,
                              r.thresh.neg=2/3,
                              p.thresh=0.05,
                              output.type="boolean"  #--> "boolean", "row.names"
                            ){
  
  reg.idx <- rep(FALSE, length(ratios))
  na.idx <- which(is.na(p.values))
  print(paste0("length(na.idx): ", length(na.idx)))
  
  #print(table(reg.idx))
  if(direction == "positive"){
    for(i in 1:length(ratios)){
      if(is.na(ratios[i]) || is.na(p.values[i])) next
      if(ratios[i] >= r.thresh.pos && p.values[i] < p.thresh) reg.idx[i] <- TRUE
    }  
  }else if(direction == "negative"){
    for(i in 1:length(ratios)){
      if(is.na(ratios[i]) || is.na(p.values[i])) next
      if(ratios[i] <= r.thresh.neg && p.values[i] < p.thresh) reg.idx[i] <- TRUE
    }
  }else if(direction == "both"){
    for(i in 1:length(ratios)){
      if(is.na(ratios[i]) || is.na(p.values[i])) next
      if(ratios[i] >= r.thresh.pos && p.values[i] < p.thresh) reg.idx[i] <- TRUE
      if(ratios[i] <= r.thresh.neg && p.values[i] < p.thresh) reg.idx[i] <- TRUE
    }
  }else if(direction == "both2"){
    for(i in 1:length(ratios)){
      if(is.na(ratios[i]) || is.na(p.values[i])) next
      if((ratios[i] >= 1.5) | (ratios[i] <= 2/3) & p.values[i] <= 0.05) reg.idx[i] <- TRUE
    }  
  }else if(direction == "herwig"){
    for(i in 1:length(ratios)){
      if(is.na(p.values[i])) next
      if(p.values[i] < p.thresh) reg.idx[i] <- TRUE
    }  
  }else if(direction == "herwig-pd"){
    print(paste0("WARNING! The 'herwig-pd'-option simulates PD's p-value filte", 
          "ring, where rounding results in less number of matching peptides!"))
    for(i in 1:length(ratios)){
      if(is.na(p.values[i])) next
      if(round(p.values[i],2) < p.thresh) reg.idx[i] <- TRUE
    }  
  }

  print("table(reg.idx): ")
  print(table(reg.idx))
  print(paste0("length(unique(c(...)): ", length(unique(c(names(ratios)[reg.idx], names(ratios)[na.idx])))))
  print("*********************************************************************")
  
  if(length(na.idx) > 0) reg.idx[na.idx] <- TRUE
  
  if(output.type == "boolean"){
    return(reg.idx)      
  }else if(output.type == "row.names"){
    return(names(ratios)[reg.idx])
  }
    
}

#++++++++++++++++++++++ regulatedPeptides(...) +++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ returnUniProtGeneNames(...) ++++++++++++++++++++++++++++++++++++

returnUniProtGeneNames <- function(uniprot.ids=NULL, table=FALSE, iso.rm=TRUE, source=NULL, file.name=NULL, output.path=NULL){
  
  library(httr)
  
  if(iso.rm == TRUE){
      uniprot.ids <- unique(gsub("\\-\\d+", "", x=uniprot.ids))
      print(paste0("returnUniProtGeneNames() - uniprot-IDs: ", length(uniprot.ids)))
  }else if(iso.rm == FALSE){
      print(paste0("returnUniProtGeneNames() - uniprot-IDs: ", length(uniprot.ids)))  
  }
  
  if(is.null(source)){
    #chunk.size <- 820 #seems to be the query maximum of UniProt API... (in Nov. 2021)
    chunk.size <- 100 #seems to be the query maximum of UniProt API... (in Mai 2022)
    chunk.number <- ceiling(length(uniprot.ids) / chunk.size)
    res.string <- c()
    
    for(i in 1:chunk.number){
      
      idx1 <- (i-1)*chunk.size+1
      idx2 <- min(i*chunk.size, length(uniprot.ids))
      #print(paste0("idx1: ", idx1))
      #print(paste0("idx2: ", idx2))
      
      uniprot.ids.string <- paste0(uniprot.ids[idx1:idx2], collapse=" ")
      
      myurl <- "https://www.uniprot.org/uploadlists/"
      myquery <- list(from="ACC+ID", to="GENENAME", format="tab", query=uniprot.ids.string)
      res <- GET(url=myurl, query=myquery)
      #print(res)
      res.string.tmp <- httr::content(x=res, as="text", encoding="UTF-8")
      #res.string.tmp <- content(x=res)
      if(i > 1) res.string.tmp <- gsub("From\tTo\n", "", x=res.string.tmp) #remove column names for chunks > chunk 1
      #print(res.string.tmp)
      res.string <- c(res.string, res.string.tmp)
      print(res.string)
      
    }
    res.string <- paste0(res.string, collapse="")
    res.tab <- read.table(text=res.string, sep="\t", header=TRUE)
  }else{
    map.tab <- read.table(file=source, header=TRUE, sep="\t")
    map.tab[,"To"] <- gsub("-", "zzzzz", map.tab[,"To"])
    map.tab <- map.tab[order(map.tab[,"To"]),]
    map.tab[,"To"] <- gsub("zzzzz", "-", map.tab[,"To"])
    duplicate.idx <- which(duplicated(map.tab[,"From"]))
    map.tab <- map.tab[-duplicate.idx,]
    rownames(map.tab) <- map.tab[,"From"]
    
    res.tab <- matrix(ncol=2, nrow=length(uniprot.ids))
    colnames(res.tab) <- c("From","To")
    
    for(i in 1:length(uniprot.ids)){
      if(grepl(";", uniprot.ids[i], fixed=T)){
        splitted.group <- unlist(strsplit(uniprot.ids[i], "\\; "))
        splitted.group <- unique(splitted.group)
        for(j in 1:length(splitted.group)){
          splitted.group[j] <- map.tab[splitted.group[j],"To"]
        }
        res.tab[i,] <- c(uniprot.ids[i], paste0(splitted.group, collapse="; "))
      }else{
        res.tab[i,] <- c(uniprot.ids[i], map.tab[uniprot.ids[i],"To"])  
      }
    }
    
    res.tab <- as.data.frame(res.tab)
  }
  
  if(!is.null(output.path)) write.table(x=res.tab, file=paste0(output.path, "/", file.name, ".txt"), sep="\t", row.names=FALSE)
  
  if(table == TRUE){
    return(res.tab)  
  }else{
    return(res.tab$To)
  }
}

#++++++++++++++++++++++ returnUniProtGeneNames(...) ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ returnGeneNamesUniProtIDs(...) ++++++++++++++++++++++++++++++++++++

returnGeneNamesUniProtIDs <- function(gene.names=NULL, table=FALSE, source=NULL, file.name=NULL, output.path=NULL){
  
  if(!is.null(source)){
    map.tab <- read.table(file=source, header=TRUE, sep="\t")
    map.tab[,"To"] <- gsub("-", "zzzzz", map.tab[,"To"])
    map.tab <- map.tab[order(map.tab[,"To"]),]
    map.tab[,"To"] <- gsub("zzzzz", "-", map.tab[,"To"])
    duplicate.idx <- which(duplicated(map.tab[,"To"]))
    map.tab <- map.tab[-duplicate.idx,]
    rownames(map.tab) <- map.tab[,"To"]
    
    res.tab <- matrix(ncol=2, nrow=length(gene.names))
    colnames(res.tab) <- c("From","To")
    
    for(i in 1:length(gene.names)){
      res.tab[i,] <- c(gene.names[i], map.tab[gene.names[i],"From"]) 
    }
    
    res.tab <- as.data.frame(res.tab)  
  }else{
    stop("No source file available!")
  }
  
  if(!is.null(output.path)) write.table(x=res.tab, file=paste0(output.path, "/", file.name, ".txt"), sep="\t", row.names=FALSE)
  
  if(table == TRUE){
    return(res.tab)  
  }else{
    return(res.tab$To)
  }
  
}

#++++++++++++++++++++++ returnGeneNamesUniProtIDs(...) ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ returnUniProtEntrezGeneIDs(...) ++++++++++++++++++++++++++++++++++++

returnUniProtEntrezGeneIDs <- function(uniprot.ids=NULL, table=FALSE, iso.rm=TRUE, source=NULL, na.rm=FALSE, file.name=NULL, output.path=NULL){
  
  library(httr)
  
  if(iso.rm == TRUE){
    uniprot.ids <- unique(gsub("\\-\\d+", "", x=uniprot.ids))
    print(paste0("returnUniProtEntrezGeneIDs() - unique uniprot-IDs: ", length(uniprot.ids)))
  }else if (iso.rm == FALSE){
    print(paste0("returnUniProtEntrezGeneIDs() - uniprot-IDs: ", length(uniprot.ids)))  
  }
  
  if(is.null(source)){
    chunk.size <- 820 #seems to be the query maximum of UniProt API... (in Nov. 2021)
    chunk.number <- ceiling(length(uniprot.ids) / chunk.size)
    res.string <- c()
    
    for(i in 1:chunk.number){
      
      idx1 <- (i-1)*chunk.size+1
      idx2 <- min(i*chunk.size, length(uniprot.ids))
      #print(paste0("idx1: ", idx1))
      #print(paste0("idx2: ", idx2))
      
      uniprot.ids.string <- paste0(uniprot.ids[idx1:idx2], collapse=" ")
      
      myurl <- "https://www.uniprot.org/uploadlists/"
      myquery <- list(from="ACC+ID", to="P_ENTREZGENEID", format="tab", query=uniprot.ids.string)
      res <- GET(url=myurl, query=myquery)
      #print(res)
      #res.string.tmp <- httr::content(res, "text")
      res.string.tmp <- httr::content(x=res, as="text", encoding="UTF-8")
      #res.string.tmp <- content(res)
      if(i > 1) res.string.tmp <- gsub("From\tTo\n", "", x=res.string.tmp) #remove column names for chunks > chunk 1
      #print(res.string.tmp)
      res.string <- c(res.string, res.string.tmp)
        
    }
    res.string <- paste0(res.string, collapse="")
    res.tab <- read.table(text=res.string, sep="\t", header=TRUE)
  
  }else{
    map.tab <- read.table(file=source, header=TRUE, sep="\t")
    map.tab[,"To"] <- gsub("-", "zzzzz", map.tab[,"To"])
    map.tab <- map.tab[order(map.tab[,"To"]),]
    map.tab[,"To"] <- gsub("zzzzz", "-", map.tab[,"To"])
    duplicate.idx <- which(duplicated(map.tab[,"From"]))
    map.tab <- map.tab[-duplicate.idx,]
    rownames(map.tab) <- map.tab[,"From"]
    
    res.tab <- matrix(ncol=2, nrow=length(uniprot.ids))
    colnames(res.tab) <- c("From","To")
    
    for(i in 1:length(uniprot.ids)){
      if(grepl(";", uniprot.ids[i], fixed=T)){
        splitted.group <- unlist(strsplit(uniprot.ids[i], "\\; "))
        splitted.group <- unique(splitted.group)
        for(j in 1:length(splitted.group)){
          splitted.group[j] <- map.tab[splitted.group[j],"To"]
        }
        res.tab[i,] <- c(uniprot.ids[i], paste0(splitted.group, collapse="; "))
      }else{
        res.tab[i,] <- c(uniprot.ids[i], map.tab[uniprot.ids[i],"To"])  
      }
    }
    
    res.tab <- as.data.frame(res.tab)  
  }
  
  if(na.rm==TRUE){
    del.idx <- which(is.na(res.tab$To))
    if(length(res.tab > 0)) res.tab <- res.tab[-del.idx,]
  }
  
  if(!is.null(output.path)) write.table(x=res.tab, file=paste0(output.path, "/", file.name, ".txt"), sep="\t", row.names=FALSE)
  
  if(table == TRUE){
    return(res.tab)  
  }else{
    return(res.tab$To)
  }
}

#++++++++++++++++++++++ returnUniProtEntrezGeneIDs(...) ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ returnEntrezUniProtIDs(...) ++++++++++++++++++++++++++++++++++++

returnEntrezUniProtIDs <- function(entrez.ids=NULL, table=FALSE, source=NULL, swissprot=FALSE, na.rm=FALSE, file.name=NULL, output.path=NULL){
  
  library(httr)
  
  print(paste0("returnEntrezUniProtIDs() - unique Entrez Gene-IDs: ", length(entrez.ids)))
  
  if(is.null(source)){
    chunk.size <- 820 #seems to be the query maximum of UniProt API... (in Nov. 2021)
    chunk.number <- ceiling(length(entrez.ids) / chunk.size)
    res.string <- c()
    
    for(i in 1:chunk.number){
      
      idx1 <- (i-1)*chunk.size+1
      idx2 <- min(i*chunk.size, length(entrez.ids))
      #print(paste0("idx1: ", idx1))
      #print(paste0("idx2: ", idx2))
        
      entrez.ids.string <- paste0(entrez.ids[idx1:idx2], collapse=" ")
      
      myurl <- "https://www.uniprot.org/uploadlists/"
      if(swissprot == TRUE){
        myquery <- list(from="P_ENTREZGENEID", to="SWISSPROT", format="tab", query=entrez.ids.string)
      }else{
        myquery <- list(from="P_ENTREZGENEID", to="ACC", format="tab", query=entrez.ids.string)
      }
      res <- GET(url=myurl, query=myquery)
      #print(res)
      #res.string.tmp <- httr::content(res, "text")
      res.string.tmp <- httr::content(x=res, as="text", encoding="UTF-8")
      #res.string.tmp <- content(res)
      if(i > 1) res.string.tmp <- gsub("From\tTo\n", "", x=res.string.tmp) #remove column names for chunks > chunk 1
      #print(res.string.tmp)
      res.string <- c(res.string, res.string.tmp)
      
    }
    res.string <- paste0(res.string, collapse="")
    res.tab <- read.table(text=res.string, sep="\t", header=TRUE)
    
    if(!is.null(output.path)) write.table(x=res.tab, file=paste0(output.path, "/", file.name, ".txt"), sep="\t", row.names=FALSE)
    
    return(res.tab$To)
  }else{
    map.tab <- read.table(file=source, header=TRUE, sep="\t")
    del.idx <- which(duplicated(map.tab[,"To"]))
    entrez.duplicated <- map.tab[duplicated(map.tab[,"To"]),"To"]
    for(i in 1:length(entrez.duplicated)){
      map.tab[map.tab[,"To"] %in% entrez.duplicated[i],"From"] <- paste0(map.tab[map.tab[,"To"] %in% entrez.duplicated[i],"From"], collapse="/")
    }
    map.tab <- map.tab[-del.idx,]
    rownames(map.tab) <- map.tab[,"To"]
    
    res.tab <- matrix(ncol=2, nrow=length(entrez.ids))
    colnames(res.tab) <- c("From","To")
    
    for(i in 1:length(entrez.ids)){
      res.tab[i,] <- c(entrez.ids[i], map.tab[entrez.ids[i],"From"])
    }
    
    res.tab <- as.data.frame(res.tab)  
  
  
    if(na.rm==TRUE){
      del.idx <- which(is.na(res.tab$To))
      if(length(res.tab > 0)) res.tab <- res.tab[-del.idx,]
    }
  
    if(!is.null(output.path)) write.table(x=res.tab, file=paste0(output.path, "/", file.name, ".txt"), sep="\t", row.names=FALSE)
  
    if(table == TRUE){
      return(res.tab)  
    }else{
      return(res.tab$To)
    }
  }
}

#++++++++++++++++++++++ returnEntrezUniProtIDs(...) ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ timeCourseDotPlot(...) +++++++++++++++++++++++++++++++++

# ontologies: "BP", "MF" & "CC"

timeCourseDotPlot <- function(x=NULL,
                              time.course=1,
                              ontology="BP",
                              nodeSize=5,
                              fdr.thresh=0.01,
                              min.anno=5,
                              max.anno=700,
                              term.max=20,
                              cwd=NULL){
  
  #library(biomaRt) ---> NOTE: biomaRt gives only subset of gene names associated with entrez ids
  library(topGO)
  library(ggplot2)
  
  output.path <- paste0(cwd, "/Sankey")
  if(!is.null(cwd) && !dir.exists(output.path)) {
    dir.create(output.path)  
  }
  
  if(time.course==1){
    time.pnts <- c("Reg_1_1","Reg_2_2.5","Reg_3_5","Reg_4_15","Reg_5_30","Reg_6_60")
    time.pnts2 <- c("1 min", "2.5 min", "5 min", "15 min", "30 min", "60 min")
  }else if(time.course==2){
    time.pnts <- c("Reg_1_2.5", "Reg_2_5", "Reg_3_15", "Reg_4_30", "Reg_5_60")
    time.pnts2 <- c("2.5 min", "5 min", "15 min", "30 min", "60 min")  
  }else if(time.course==3){
    time.pnts <- c("Reg_1_5", "Reg_2_15", "Reg_3_30", "Reg_4_60")
    time.pnts2 <- c("5 min", "15 min", "30 min", "60 min")  
  }else if(time.course==4){
    time.pnts <- c("Reg_1_15", "Reg_2_30", "Reg_3_60")
    time.pnts2 <- c("15 min", "30 min", "60 min")  
  }else if(time.course==5){
    time.pnts <- c("Reg_1_30", "Reg_2_60")
    time.pnts2 <- c("30 min", "60 min")  
  }else if(time.course==6){
    time.pnts <- c("Reg_1_60")
    time.pnts2 <- c("60 min")  
  }
  
  dotplot.dat <- matrix(ncol=4)
  colnames(dotplot.dat) <- c("Term", "TimePoint", "Strength", "FDR")
  dotplot.dat <- as.data.frame(dotplot.dat)
  for(i in time.course:6){
    idx.tmp <- grep(time.pnts[i-time.course+1], x[,2+i])
    pep.id.tmp <- x[idx.tmp,1]
    
    uniprot.tmp <- unique(pep2prot(x[idx.tmp, 2]))
    print(paste0("uniprot.tmp: ", length(uniprot.tmp)))
    regulated.tab <- returnUniProtEntrezGeneIDs(uniprot.ids=uniprot.tmp, make.unique=FALSE, table=TRUE)
    ############################################################################
    # NOTE: Next rows are CRUCIAL, to fix UniProt-IDs with multiple Entrez-IDs #
    ############################################################################
    if(length(unique(regulated.tab$From)) < length(regulated.tab$To)){
      duplicate.idx <- which(duplicated(regulated.tab$From))
      del.idx.tmp <- c()
      #Find duplicated Entrez-IDs with minimal value
      for(j in 1:length(duplicate.idx)){
          idx.tmp <- grep(regulated.tab$From[duplicate.idx[j]], regulated.tab$From)
          keep.tmp <- min(as.numeric(regulated.tab$To[idx.tmp]))
          keep.tmp.idx <- grep(keep.tmp, regulated.tab$To)
          del.idx.tmp <- c(del.idx.tmp, idx.tmp[-grep(keep.tmp.idx, idx.tmp)])
      }
      regulated.tab <- regulated.tab[-del.idx.tmp,]
    }
    ############################################################################
    #NOTE: Next above are CRUCIAL, to fix UniProt-IDs with multiple Entrez-IDs #
    ############################################################################
    regulated <- regulated.tab$To
    print(paste0("regulated: ", length(regulated)))
    gene.tmp <- returnUniProtGeneNames(uniprot.ids=regulated.tab$From, make.unique=FALSE, table=FALSE)
    print(paste0("gene.tmp: ", length(gene.tmp)))
    gene.tmp.tab <- returnUniProtGeneNames(uniprot.ids=regulated.tab$From, make.unique=FALSE, table=TRUE)
    
    
    entrez.mapping <- cbind(regulated, regulated.tab$From, gene.tmp)
    colnames(entrez.mapping) <- c("Entrez", "UniProt", "Genes") 
    rownames(entrez.mapping) <- regulated
    
    if(i==1){
      write.table(x=gene.tmp.tab, file=paste0(output.path, "/last_gene_tmp_tab.txt"), col.names=TRUE, row.names=FALSE, sep="\t")
      write.table(x=regulated.tab, file=paste0(output.path, "/last_regulated_tmp_tab.txt"), col.names=TRUE, row.names=FALSE, sep="\t")
      write.table(x=entrez.mapping, file=paste0(output.path, "/last_entrez_mapping.txt"), col.names=TRUE, row.names=FALSE, sep="\t")
    }
    
    gene.annotation <- annFUN.org(whichOnto=ontology,
                                  feasibleGenes=NULL,
                                  mapping="org.Hs.eg.db",
                                  ID="entrez")
    
    allGenes <- unique(unlist(gene.annotation))
    geneList <- factor(as.integer(allGenes %in% regulated))
    names(geneList) <- allGenes
    
    sampleGOdata <- new("topGOdata",
                        ontology=ontology,
                        allGenes=geneList,
                        nodeSize=nodeSize,
                        annot=annFUN.org,
                        mapping="org.Hs.eg.db",
                        ID="entrez")
    
    resultFisher <- runTest(object=sampleGOdata,
                            algorithm="classic",
                            statistic = "fisher")
    
    score.adj <- p.adjust(resultFisher@score, method="fdr")
    topNodes <- sum(score.adj < fdr.thresh, na.rm=T)
    
    if(topNodes < 2) stop(paste0("timeCourseDotPlot - Error: There are/is only ", topNodes," term(s) with FDR smaller then ", fdr.thresh, "!!!"))
    if(topNodes < 10) warning(paste0("timeCourseDotPlot - Warning: There are only ", topNodes, " terms with FDR smaller then ", fdr.thresh, "!!!"))
    
    resTable <- GenTable(sampleGOdata,
                         Fisher=resultFisher,
                         topNodes=topNodes,
                         numChar=1000)
    print(resTable)
    
    if(nrow(resTable) < 10) warning(paste0("timeCourseDotPlot - Warning: 'resTab' has only ", nrow(resTable), " rows!!!"))
    
    FDR <- p.adjust(resTable$Fisher, method="fdr")
    resTable <- cbind(resTable, FDR)
    
    Strength <- vector(mode="numeric", length=topNodes)
    for(j in 1:topNodes){
      Strength[j] <- log10(resTable[j,"Significant"] / resTable[j, "Expected"]) 
    }
    resTable <- cbind(resTable, Strength)
    resTable <- resTable[order(resTable$Strength, decreasing=T),]
    
    # retrieve "all GO to genes" list from the "expanded" annotation in GOdata
    allGO2genes <- genesInTerm(sampleGOdata)
    GO2regGenes <- lapply(allGO2genes, function(x) x[x %in% regulated])
    
    #mart <- useMart(dataset="hsapiens_gene_ensembl", biomart='ensembl')
    Entrez <- vector(mode="character", length=topNodes)
    Genes <- vector(mode="character", length=topNodes)
    UniProt <- vector(mode="character", length=topNodes)
    for(j in 1:topNodes){
      entrez.id.tmp <- GO2regGenes[[resTable$GO.ID[j]]]
      Entrez[j] <- paste(entrez.id.tmp, collapse=", ")
      #genes.tmp <- biomaRt::select(mart, keys=GO2regGenes[[resTable$GO.ID[j]]], columns=c('hgnc_symbol'), keytype='entrezgene_id')
      uniprot.id.tmp <- entrez.mapping[entrez.id.tmp, "UniProt"]
      UniProt[j] <- paste(uniprot.id.tmp, collapse=", ")
      genes.tmp <- entrez.mapping[entrez.id.tmp, "Genes"]
      Genes[j] <- paste(sort(genes.tmp), collapse=", ")
    }
    resTable <- cbind(resTable, UniProt, Entrez, Genes)
    
    
    
    Filter <- rep(0, nrow(resTable))
    if(sum(resTable$Significant >= min.anno, na.rm=T) < 1){
      stop(paste0("timeCourseDotPlot - Error: In 'resTab' only ", sum(resTable$Significant >= min.anno, na.rm=T), " >= min.anno!!!"))
    }
    #resTable <- resTable[resTable$Significant >= min.anno,]
    Filter[resTable$Significant < min.anno] <- 1
    
    if(sum(resTable$Annotated <= max.anno, na.rm=T) < 1){
      stop(paste0("timeCourseDotPlot - Error: In 'resTab' only ", sum(resTable$Annotated <= max.anno, na.rm=T), " <= max.anno!!!"))
    }
    #resTable <- resTable[resTable$Annotated <= max.anno,]
    Filter[resTable$Annotated > max.anno] <- 1
    resTable <- cbind(resTable, Filter)

    
    
    Delete <- rep(0, nrow(resTable))
    for(j in 1:nrow(resTable)){
      if(resTable$Filter[j] == 1) next
      if(j == nrow(resTable)) break
      vec_j <- strsplit(resTable$Genes[j], split=", ")[[1]]
      n_j <- length(vec_j)
      for(k in (j+1):nrow(resTable)){
        vec_k <- strsplit(resTable$Genes[k], split=", ")[[1]]
        n_k <- length(vec_k)
        n_overlap <- length(intersect(vec_j, vec_k))
        if(n_overlap / n_j * 100 >= 90){
          Delete[k] <- 1
        }
      }  
    }
    resTable <- cbind(resTable, Delete)
    
    
    
    print(resTable[,1:8])
    write.table(x=resTable,
                file=paste0(output.path, "/resTable_", ontology, "_tc", time.course, "_tp", i, ".txt"),
                sep="\t",
                col.names=T,
                row.names=F)
    
    
    
    del.idx.resTable1 <- grep("1", resTable$Filter)
    del.idx.resTable2 <- grep("1", resTable$Delete)
    del.idx.resTable <- sort(unique(c(del.idx.resTable1, del.idx.resTable2)))
    resTable <- resTable[-del.idx.resTable,]
    
    
    
    if(nrow(resTable) >= term.max){
      term.n <- term.max
    }else{
      term.n <- nrow(resTable)
    }
    for(j in 1:term.n){
      row.tmp <- c(resTable$Term[j], time.pnts2[i-time.course+1], resTable$Strength[j], resTable$FDR[j])
      #print(row.tmp)
      dotplot.dat <- rbind(dotplot.dat, row.tmp)
    }
    
  }
  
  dotplot.dat$Strength <- as.numeric(dotplot.dat$Strength)
  dotplot.dat$FDR <- as.numeric(dotplot.dat$FDR)
  
  pl <- ggplot(data=dotplot.dat[-1,], aes(x=factor(TimePoint, level=time.pnts2), y=Term, 
                                          color=FDR, size=Strength)) + 
    geom_point() +
    scale_y_discrete(limits=rev) +
    scale_color_gradient(low = "red", high = "blue") +
    theme(
      plot.title = element_text(color="black", size=14, face="bold"),
      axis.text.x = element_text(size=12, face="bold", angle=90, vjust=0.5)
    ) + 
    #theme_bw() + 
    ylab("") + 
    xlab("") + 
    ggtitle(paste0("GO-based ORA (", ontology, ", tc", time.course, ")"))
  
  print(pl)
  
  png(filename=paste0(output.path, "/sankey_dotplot_", ontology, "_", time.course, ".png "), height=2600, width=2400, res=300) #3400 x 2400
      print(pl)
  dev.off()
  
}

#++++++++++++++++++++++ timeCourseDotPlot(...) +++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ timeCourseDotPlot2(...) +++++++++++++++++++++++++++++++++

# ontologies: "BP", "MF" & "CC"

timeCourseDotPlot2 <- function(x=NULL,
                              time.course=1,
                              ontology="BP",
                              nodeSize=5,
                              fdr.thresh=0.01,
                              min.anno=5,
                              max.anno=700,
                              term.max=20,
                              cwd=NULL){
  
  library(topGO)
  library(ggplot2)
  
  output.path <- paste0(cwd, "/Sankey/dereg-dotplots")
  if(!is.null(cwd) && !dir.exists(output.path)) {
    dir.create(output.path)  
  }
  
  if(time.course==1){
    time.pnts <- c("Reg_1_1","Reg_2_2.5","Reg_3_5","Reg_4_15","Reg_5_30","Reg_6_60")
    time.pnts2 <- c("2.5 min", "5 min", "15 min", "30 min", "60 min")
  }else if(time.course==2){
    time.pnts <- c("Reg_1_2.5", "Reg_2_5", "Reg_3_15", "Reg_4_30", "Reg_5_60")
    time.pnts2 <- c("5 min", "15 min", "30 min", "60 min")  
  }else if(time.course==3){
    time.pnts <- c("Reg_1_5", "Reg_2_15", "Reg_3_30", "Reg_4_60")
    time.pnts2 <- c("15 min", "30 min", "60 min")  
  }else if(time.course==4){
    time.pnts <- c("Reg_1_15", "Reg_2_30", "Reg_3_60")
    time.pnts2 <- c("30 min", "60 min")  
  }else if(time.course==5){
    time.pnts <- c("Reg_1_30", "Reg_2_60")
    time.pnts2 <- c("60 min")  
  }
  
  dotplot.dat <- matrix(ncol=4)
  colnames(dotplot.dat) <- c("Term", "TimePoint", "Strength", "FDR")
  dotplot.dat <- as.data.frame(dotplot.dat)
  for(i in time.course:5){
    idx.tmp1 <- grep(time.pnts[i-time.course+1], x[,2+i])
    idx.tmp2 <- grep(time.pnts[i-time.course+1+1], x[,2+i+1])
    idx.tmp <- setdiff(idx.tmp1, idx.tmp2)
    pep.id.tmp <- x[idx.tmp,1]
    
    x.tmp <- unique(pep2prot(x[idx.tmp, 2]))
    regulated <- returnUniProtEntrezGeneIDs(uniprot.ids=x.tmp)
    
    gene.annotation <- annFUN.org(whichOnto=ontology,
                                  feasibleGenes=NULL,
                                  mapping="org.Hs.eg.db",
                                  ID="entrez")
    
    allGenes <- unique(unlist(gene.annotation))
    geneList <- factor(as.integer(allGenes %in% regulated))
    names(geneList) <- allGenes
    
    sampleGOdata <- new("topGOdata",
                        ontology=ontology,
                        allGenes=geneList,
                        nodeSize=nodeSize,
                        annot=annFUN.org,
                        mapping="org.Hs.eg.db",
                        ID="entrez")
    
    resultFisher <- runTest(object=sampleGOdata,
                            algorithm="classic",
                            statistic = "fisher")
    
    score.adj <- p.adjust(resultFisher@score, method="fdr")
    topNodes <- sum(score.adj < fdr.thresh, na.rm=T)
    
    if(topNodes < 2) stop(paste0("timeCourseDotPlot2 - Error: There are/is only ", topNodes," term(s) with FDR smaller then ", fdr.thresh, "!!!"))
    if(topNodes < 10) warning(paste0("timeCourseDotPlot2 - Warning: There are only ", topNodes, " terms with FDR smaller then ", fdr.thresh, "!!!"))
    
    resTable <- GenTable(sampleGOdata,
                         Fisher=resultFisher,
                         topNodes=topNodes,
                         numChar=1000)
    
    print(resTable)
    
    if(nrow(resTable) < 10) warning(paste0("timeCourseDotPlot2 - Warning: 'resTab' has only ", nrow(resTable), " rows!!!"))
    
    FDR <- p.adjust(resTable$Fisher, method="fdr")
    resTable <- cbind(resTable, FDR)
    
    Strength <- vector(mode="numeric", length=topNodes)
    for(j in 1:topNodes){
      Strength[j] <- log10(resTable[j,"Significant"] / resTable[j, "Expected"]) 
    }
    resTable <- cbind(resTable, Strength)
    resTable <- resTable[order(resTable$Strength, decreasing=T),]
    
    print(resTable)
    write.table(x=resTable,
                file=paste0(output.path, "/dereg_resTable_", ontology, "_tc", time.course, "_tp", i, "_peps", length(idx.tmp), "_genes", length(regulated), ".txt"),
                sep="\t",
                col.names=T,
                row.names=F)
    
    
    
    #print(paste0("resTable$Significant >= min.anno: ", sum(resTable$Significant >= min.anno, na.rm=T)))
    if(sum(resTable$Significant >= min.anno, na.rm=T) < 1){
      stop(paste0("timeCourseDotPlot2 - Error: In 'resTab' only ", sum(resTable$Significant >= min.anno, na.rm=T), " >= min.anno!!!"))
    }
    resTable <- resTable[resTable$Significant >= min.anno,]
    
    #print(paste0("resTable$Annotated <= max.anno: ", sum(resTable$Annotated <= max.anno, na.rm=T)))
    if(sum(resTable$Annotated <= max.anno, na.rm=T) < 1){
      stop(paste0("timeCourseDotPlot2 - Error: In 'resTab' only ", sum(resTable$Annotated <= max.anno, na.rm=T), " <= max.anno!!!"))
    }
    resTable <- resTable[resTable$Annotated <= max.anno,]
    
    
    
    if(nrow(resTable) >= term.max){
      term.n <- term.max
    }else{
      term.n <- nrow(resTable)
    }
    for(j in 1:term.n){
      row.tmp <- c(resTable$Term[j], time.pnts2[i-time.course+1], resTable$Strength[j], resTable$FDR[j])
      #print(row.tmp)
      dotplot.dat <- rbind(dotplot.dat, row.tmp)
    }
    
  }
  
  dotplot.dat$Strength <- as.numeric(dotplot.dat$Strength)
  dotplot.dat$FDR <- as.numeric(dotplot.dat$FDR)
  
  pl <- ggplot(data=dotplot.dat[-1,], aes(x=factor(TimePoint, level=time.pnts2), y=Term, 
                                          color=FDR, size=Strength)) + 
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme(
      plot.title = element_text(color="black", size=14, face="bold"),
      axis.text.x = element_text(size=12, face="bold", angle=90, vjust=0.5)
    ) + 
    #theme_bw() + 
    ylab("") + 
    xlab("") + 
    ggtitle(paste0("GO-ORA dereg. peptides (", ontology, ", tc", time.course, ")"))
  
  print(pl)
  
  png(filename=paste0(output.path, "/dereg_sankey_dotplot_", ontology, "_", time.course, ".png "), height=3400, width=2400, res=300)
  print(pl)
  dev.off()
  
}

#++++++++++++++++++++++ timeCourseDotPlot2(...) +++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ updateUniProtIDs(...) ++++++++++++++++++++++++++++++++++++

updateUniProtIDs <- function(uniprot.ids=NULL, source=NULL, table=FALSE, make.unique=TRUE, file.name=NULL, output.path=NULL){
  
  library(httr)
  
  if(make.unique == TRUE){
    uniprot.ids <- unique(gsub("\\-\\d+", "", x=uniprot.ids))
    print(paste0("unique uniprot-IDs: ", length(uniprot.ids)))
  }else if (make.unique == FALSE){
    print(paste0("uniprot-IDs: ", length(uniprot.ids)))  
  }
  
  if(is.null(source)){
    chunk.size <- 820 #seems to be the query maximum of UniProt API... (in Nov. 2021)
    chunk.number <- ceiling(length(uniprot.ids) / chunk.size)
    res.string <- c()
    
    for(i in 1:chunk.number){
      
      idx1 <- (i-1)*chunk.size+1
      idx2 <- min(i*chunk.size, length(uniprot.ids))
      print(paste0("idx1: ", idx1))
      print(paste0("idx2: ", idx2))
      
      uniprot.ids.string <- paste0(uniprot.ids[idx1:idx2], collapse=" ")
      
      myurl <- "https://www.uniprot.org/uploadlists/"
      myquery <- list(from="ACC+ID", to="ACC", format="tab", query=uniprot.ids.string)
      res <- GET(url=myurl, query=myquery)
      #print(res)
      res.string.tmp <- httr::content(x=res, as="text", encoding="UTF-8")
      #res.string.tmp <- content(x=res)
      if(i > 1) res.string.tmp <- gsub("From\tTo\n", "", x=res.string.tmp) #remove column names for chunks > chunk 1
      print(res.string.tmp)
      res.string <- c(res.string, res.string.tmp)
      
    }
    res.string <- paste0(res.string, collapse="")
    res.tab <- read.table(text=res.string, sep="\t", header=TRUE)
    
    if(!is.null(output.path)) write.table(x=res.tab, file=paste0(output.path, "/", file.name, ".txt"), sep="\t", row.names=FALSE)
  
    if(table == TRUE){
      return(res.tab)  
    }else{
      return(res.tab$To)
    }
  }else{
    map.tab <- read.table(file=source, header=TRUE, sep="\t")
    rownames(map.tab) <- map.tab[,"From"]
    res.tab <- matrix(ncol=2, nrow=length(uniprot.ids))
    colnames(res.tab) <- c("From","To")
    
    for(i in 1:length(uniprot.ids)){
      res.tab[i,] <- c(uniprot.ids[i], map.tab[uniprot.ids[i],"To"])
    }
    res.tab <- as.data.frame(res.tab)  
    
    if(!is.null(output.path)) write.table(x=res.tab, file=paste0(output.path, "/", file.name, ".txt"), sep="\t", row.names=FALSE)
    
    if(table == TRUE){
      return(res.tab)  
    }else{
      return(res.tab$To)
    }
  }
}

#++++++++++++++++++++++ updateUniProtIDs(...) ++++++++++++++++++++++++++++++++++++
