
library(plyr)
library(tidyverse)
library(data.table)
library(magrittr)
library(readxl)

## To summarize all the PSMs from Proteome Discoverer

PSMSum <- function(PSM, ConfidenceLevel = "High"){
  SeqRep <- function(File){
    FL2 <- File[,c("Sequence", "Modifications")]%>%
      as.data.frame()
    
    C1 <- max(str_count(FL2$Sequence, "[a-z]"))
    
    FL2 <- cbind(FL2, str_split_fixed(FL2$Sequence, "[a-z]", C1))
    FL2 <- cbind(FL2, str_split_fixed(FL2$Modifications, "; ", C1))
    
    colnames(FL2) <- c("Sequence", "Modifications", 1: (2 * C1))
    
    NumRm <- function(Column){
      Column <- str_replace(Column, "[0-9]\\(", "\\(")
      Column <- str_replace(Column, "[0-9]\\(", "\\(")
    }
    
    FL2[,(length(FL2) - C1 +1): length(FL2)] <- apply(FL2[,(length(FL2) - C1 +1): length(FL2)], 2, NumRm)
    
    FL2 <- unite(FL2, col = Sequence, (rep(c(1, C1 + 1), C1) + rep(2:(C1+1), each = 2)), sep = "", remove = FALSE)
    File$ModifiedPeptide <- str_c("_", FL2$Sequence, "_", sep = "")
    
    File$ModifiedPeptide <- str_replace_all(File$ModifiedPeptide, "\\(", "\\[")
    File$ModifiedPeptide <- str_replace_all(File$ModifiedPeptide, "\\)", "\\]")
    File$ModifiedPeptide <- str_replace_all(File$ModifiedPeptide, "\\[1\\]", "\\(1\\)")
    
    return(File)
  }
  
  PSM <- data.table(PSM)%>%
    .[,Confidence := paste(A2, A4)]%>%
    .[Confidence %like% ConfidenceLevel]
  
  PSM <- SeqRep(PSM)%>%
    .[,Identifier := str_c(ModifiedPeptide, Charge)]%>%
    .[,SeqMz := str_c(toupper(Sequence), "_", round(`m/z [Da]`))]
  
  PSM.T <- PSM%>%
    .[Modifications %like% "Hex\\(1\\)HexNAc\\(1\\)",]%>%
    .[!Modifications %like% "HexNAc\\)"]
  
  p1 <- list(Peptides = PSM$Identifier, Glycopeptides = PSM.T$Identifier)
  
  p2 <- list(Peptides = PSM$SeqMz, Glycopeptides = PSM.T$SeqMz)
  
  return(list(PSM, PSM.T, p1, p2))
}

## To summarize all the entries in Glyco-DIA libraries

LbSum <- function(DIALibrary){
  DIALibrary[,Identifier := str_c(ModifiedPeptide, PrecursorCharge)]%>%
    .[,SeqMz := str_c(StrippedPeptide, "_", round(PrecursorMz))]
  
  DIALibrary.T <- DIALibrary%>%
    .[ModifiedPeptide %like% "Hex\\(1\\)HexNAc\\(1\\)",]%>%
    .[!ModifiedPeptide %like% "HexNAc\\)"]
  
  p1 <- list(Peptides = DIALibrary$Identifier, Glycopeptides = DIALibrary.T$Identifier)
  
  p2 <- list(Peptides = DIALibrary$SeqMz, Glycopeptides = DIALibrary.T$SeqMz)
  
  return(list(DIALibrary, DIALibrary.T, p1, p2))
  
}

## To summarize result files from Spectronaut

SpResSum <- function(SpRes){
  SpRes <- SpRes[EG.IsDecoy == FALSE & EG.Identified == TRUE]%>%
    .[, EG.ModifiedSequence := str_replace_all(`EG.ModifiedSequence`, "\\+365", "Hex\\(1\\)HexNAc\\(1\\)")]
  
  SpRes[,Identifier := str_c(EG.ModifiedSequence, FG.Charge)]%>%
    .[,SeqMz := str_c(EG.StrippedSequence, "_", round(FG.PrecMzCalibrated))]
  
  SpRes.T <- SpRes%>%
    .[EG.ModifiedSequence %like% "Hex\\(1\\)HexNAc\\(1\\)",]%>%
    .[!EG.ModifiedSequence %like% "203"]
  
  p1 <- list(Peptides = SpRes$Identifier, Glycopeptides = SpRes.T$Identifier)
  
  p2 <- list(Peptides = SpRes$SeqMz, Glycopeptides = SpRes.T$SeqMz)
  
  return(list(SpRes, SpRes.T, p1, p2))
  
}


## To calculate Isomers in the PSMs from Proteome Discoverer

IsomerCal <- function(ETD){
  SS <- data.table(table(ETD$Identifier))
  colnames(SS) <- c("Identifier", "SiteSpecificNum")
  
  ETD <- merge(ETD, SS, by = "Identifier", all.x = TRUE)
  
  SeqMz <- data.table(table(ETD$SeqMz))
  colnames(SeqMz) <- c("SeqMz", "SeqMzNum")
  
  ETD <- merge(ETD, SeqMz, by = "SeqMz", all.x = TRUE)
  
  ETD[,IsomerNum := SeqMzNum - SiteSpecificNum]
  
  Isomers <- ETD[IsomerNum > 0]
  
  MzRange <- Isomers[,.(SeqMz, Identifier, `RT [min]`)]%>%
    .[order(SeqMz, `RT [min]`)]
  
  tMzRange <- data.table(table(MzRange$SeqMz))%>%
    .[order(V1)]%>%
    .[,c(1:N),by = V1]
  
  MzRange[,Time := tMzRange[,2]]
  
  ggMzRange <- dcast(MzRange, value.var = "RT [min]", SeqMz ~ Time)
  
  ggMzRange <- melt(ggMzRange, id.vars = c("SeqMz"))%>%
    .[is.na(value) == FALSE]%>%
    .[order(SeqMz, variable)]
  
  ggMzRange$value2 <- ggMzRange$value[c(length(ggMzRange$value), 1:(length(ggMzRange$value)-1))]
  
  ggMzRange[,RTDiff := value - value2]
  
  ggMzRange[variable == 1]$RTDiff <- NA
  ggMzRange$value2 <- NULL
  
  return(list(MzRange, ggMzRange))
}

ETD.LB.Generator <- function(HCD.Index, ETD.HCD.ETD, ETD.HCD.ETD.LB){
  require(plyr)
  require(tidyverse)
  require(XML)
  require(data.table)
  
  wd <- getwd()
  
  RL <- data.frame(ETD.HCD.ETD.LB)
  
  ETD.HCD.ETD <- data.table(ETD.HCD.ETD)
  ETD.HCD.ETD <- ETD.HCD.ETD[,Confidence := paste(A2, A4)]%>%
    .[Confidence %like% "High"]
  
  ETD.HCD.ETD[, ID:= 1:dim(ETD.HCD.ETD)[1]]
  colnames(ETD.HCD.ETD)<- str_replace_all(colnames(ETD.HCD.ETD), " ", ".")
  
  
  # To subtract fragment ions from .txt files
  SubFrag <- function(Address){
    File1 <- read.delim(Address, stringsAsFactors=FALSE)
    File1$File <- Address
    FragList <- rbind(FragList, File1)
  }
  
  # To replace sequence
  SeqRep <- function(File){
    FL2 <- File[,c("Sequence", "Modifications")]
    
    C1 <- max(str_count(FL2$Sequence, "[a-z]"))
    
    FL2 <- cbind(FL2, str_split_fixed(FL2$Sequence, "[a-z]", C1))
    FL2 <- cbind(FL2, str_split_fixed(FL2$Modifications, "; ", C1))
    
    colnames(FL2) <- c("Sequence", "Modifications", 1: (2 * C1))
    
    NumRm <- function(Column){
      Column <- str_replace(Column, "[0-9]\\(", "\\(")
      Column <- str_replace(Column, "[0-9]\\(", "\\(")
    }
    
    FL2[,(length(FL2) - C1 +1): length(FL2)] <- apply(FL2[,(length(FL2) - C1 +1): length(FL2)], 2, NumRm)
    #FL2 <- FL2[1:50,]
    
    FL2 <- unite(FL2, col = Sequence, (rep(c(1, C1 + 1), C1) + rep(2:(C1+1), each = 2)), sep = "", remove = FALSE)
    File$ModifiedPeptide <- str_c("_", FL2$Sequence, "_", sep = "")
    return(File)
  }
  
  # To extract matched fragments, more advanced version
  FragMatch <- function(File){
    #File <- File[,c("Sequence", "matches")]
    File$FragmentType <- str_sub(File$matches, 1, 1)
    File$FragmentNumber <- str_split_fixed(File$matches, "[\\(\\)]", 3)[,2]
    File$FragmentCharge <- str_split_fixed(File$matches, "[ \\-]", 4)[,3]
    File$FragmentCharge <- str_sub(File$FragmentCharge, 2, 2)
    File$FragmentLossType <- str_split_fixed(File$matches, "-", 2)[,2]
    File$FragmentLossType <- str_replace_all(File$FragmentLossType, "H2O-1\\(\\+C8\\+H13\\+N\\+O5\\)", "1\\(\\+C8\\+H15\\+N\\+O6\\)")
    File$FragmentLossType <- str_replace_all(File$FragmentLossType, "NH3-1\\(\\+C8\\+H13\\+N\\+O5\\)", "1\\(\\+C8\\+H16\\+N2\\+O5\\)")
    File$FragmentLossType <- str_replace_all(File$FragmentLossType, "H2O-1\\(\\+C14\\+H23\\+N\\+O10\\)", "1\\(\\+C14\\+H25\\+N\\+O11\\)")
    File$FragmentLossType <- str_replace_all(File$FragmentLossType, "NH3-1\\(\\+C14\\+H23\\+N\\+O10\\)", "1\\(\\+C14\\+H26\\+N2\\+O10\\)")
    File$FragmentLossType <- str_replace_all(File$FragmentLossType, "NH3-H2O", "NH5O")
    return(File)
  }
  
  # To convert the index.html to data.table
  IndexCal <- function(Index){
    Index <- readLines(Index)
    Index <- str_replace_all(Index, "<a href=\"", "")
    Index <- str_replace_all(Index, "\">Peak List</a>", "")
    Index <- htmlParse(Index)
    Index <- readHTMLTable(Index)
    Index <- as.data.frame(Index[["NULL"]])
    Index <- data.frame(lapply(Index, as.character), stringsAsFactors = FALSE)
    colnames(Index) <- str_replace_all(Index[1,], "Â", "")
    Index <- subset.data.frame(Index, is.na(Index$Confidence) == FALSE)
    Index <- Index[-1,seq(2, 30, by = 2)]
    colnames(Index)<- str_replace_all(colnames(Index), " ", ".")
    colnames(Index) <- str_replace_all(colnames(Index), "Î”", "")
    Index$File <- Index$`Peak.List`
    Index[,c(5:7, 9:12)] <- data.frame(lapply(Index[,c(5:7, 9:12)], as.numeric), stringsAsFactors = FALSE)
    Index <- data.table(Index)
    Index[, ID:= 1:dim(Index)[1]]
    return(Index)
  }
  
  # To calculate the precursor mass
  MassCal <- function (sequence){
    for (sequence_number in 1:length(sequence)) {
      peptide_vector <- strsplit(sequence[sequence_number], 
                                 split = "")[[1]]
      peptide_length <- length(peptide_vector)
      C <- 12
      H <- 1.0078250321
      O <- 15.9949146221
      S <- 31.97207069
      N <- 14.0030740052
      proton <- 1.007276466
      electron <- 0.00054857990943
      residueMass <- function(residue) {
        if (residue == "A") 
          mass = C * 3 + H * 5 + N + O
        if (residue == "R") 
          mass = C * 6 + H * 12 + N * 4 + O
        if (residue == "N") 
          mass = C * 4 + H * 6 + N * 2 + O * 2
        if (residue == "D") 
          mass = C * 4 + H * 5 + N + O * 3
        if (residue == "E") 
          mass = C * 5 + H * 7 + N + O * 3
        if (residue == "Q") 
          mass = C * 5 + H * 8 + N * 2 + O * 2
        if (residue == "G") 
          mass = C * 2 + H * 3 + N + O
        if (residue == "H") 
          mass = C * 6 + H * 7 + N * 3 + O
        if (residue == "I") 
          mass = C * 6 + H * 11 + N + O
        if (residue == "L") 
          mass = C * 6 + H * 11 + N + O
        if (residue == "K") 
          mass = C * 6 + H * 12 + N * 2 + O
        if (residue == "M") 
          mass = C * 5 + H * 9 + N + O + S
        if (residue == "m") 
          mass = C * 5 + H * 9 + N + O * 2 + S
        if (residue == "F") 
          mass = C * 9 + H * 9 + N + O
        if (residue == "P") 
          mass = C * 5 + H * 7 + N + O
        if (residue == "S") 
          mass = C * 3 + H * 5 + N + O * 2
        if (residue == "T") 
          mass = C * 4 + H * 7 + N + O * 2
        if (residue == "W") 
          mass = C * 11 + H * 10 + N * 2 + O
        if (residue == "Y") 
          mass = C * 9 + H * 9 + N + O * 2
        if (residue == "V") 
          mass = C * 5 + H * 9 + N + O
        if (residue == "C") 
          mass = C * 5 + H * 8 + N * 2 + O * 2 + S
        return(mass)
      }
      masses <- sapply(peptide_vector, residueMass)
      pm <- sum(masses)
      return(pm)
    }
  }
  
  # To calculater the fragment ions of naked peptides
  FragCal <- function(File){
    #To make a sequences with all fragments
    PepSeq <- function(FileFile){
      Sequence <- FileFile$Sequence
      SeqLen <- str_length(Sequence[1])
      Seq1 <- rep(Sequence, times = SeqLen)
      File1 <- data.frame(Seq1, 1: SeqLen)
      File1$FragmentType <- ""
      File1$FragmentCharge <- ""
      File1$File <- FileFile$File
      File1 <- File1[,c(5,1:4)]
      colnames(File1) <- colnames(File)[1:5]
      bList <- File1
      yList <- File1
      bList$SubSeq <- str_sub(bList$Sequence, 1, bList$FragmentNumber)
      yList$SubSeq <- str_sub(yList$Sequence, - yList$FragmentNumber)
      bList2 <- as.data.frame(bList$SubSeq, stringsAsFactors = FALSE)
      bList$Mass <- apply(bList2, 1, MassCal)
      yList2 <- as.data.frame(yList$SubSeq, stringsAsFactors = FALSE)
      yList$Mass <- apply(yList2, 1, MassCal)
      File1 <- rbind.data.frame(bList, bList, bList, yList, yList, yList)
      File1$FragmentType <- rep(c("b", "y"), each = 3 * SeqLen)
      File1$FragmentCharge <- rep(c(1, 2, 3), each = SeqLen)
      File0 <- rbind.data.frame(File0, File1)
      return(File0)
    }
    #File <- Index
    File <- File[,c("File", "Sequence")]
    File$Sequence <- str_replace_all(File$Sequence, "s", "S")
    File$Sequence <- str_replace_all(File$Sequence, "t", "T")
    File$Sequence <- str_replace_all(File$Sequence, "y", "Y")
    File$Sequence <- str_replace_all(File$Sequence, "c", "C")
    
    #colnames(File) <- "Sequence"
    File$FragmentNumber <- str_length(File$Sequence)
    File$FragmentType <- ""
    File$FragmentCharge <- ""
    File$SubSeq <- ""
    File$Mass <- ""
    File$Sequence <- as.character(File$Sequence)
    File0 <- File[0,]
    File <- ddply(File, 1, PepSeq)
    
    File$FragmentLoss <- ""
    FileH2O <- File
    FileH2O$FragmentLoss <- "-H2O"
    FileH2O$Mass <- FileH2O$Mass - 18.01056
    
    FileNH3 <- File
    FileNH3$FragmentLoss <- "-NH3"
    FileNH3$Mass <- FileNH3$Mass - 17.02655
    
    File <- rbind.data.frame(File, FileH2O, FileNH3)
    
    colnames(File)[1] <- "FileID"
    
    #y fragments is the same as precursor fragments, so here keep y(length) as precursor ion
    File <- subset.data.frame(File, File$FragmentType == "y" | (File$Sequence != File$SubSeq))
    File456 <- subset.data.frame(File, File$FragmentType == "y" & (File$Sequence == File$SubSeq))
    File456$FragmentCharge <- File456$FragmentCharge + 3
    File <- rbind.data.frame(File, File456)
    
    FileT <- subset.data.frame(File, File$FragmentType == "y" & (File$Sequence == File$SubSeq))
    FileTn <- FileT
    
    FileT$Mass <- FileT$Mass + 365.132196
    FileTn$Mass <- FileTn$Mass + 203.079373
    
    FileT$FragmentLoss <- str_c(FileT$FragmentLoss, "-1(+C14+H23+N+O10)")
    FileTn$FragmentLoss <- str_c(FileTn$FragmentLoss, "-1(+C8+H13+N+O5)")
    
    File <- rbind.data.frame(File, FileT, FileTn)
    
    bIndex <- File$FragmentType == "b"
    yIndex <- File$FragmentType == "y"
    File[bIndex,]$Mass <- (File[bIndex,]$Mass + 1.007276466 * File[bIndex,]$FragmentCharge) / File[bIndex,]$FragmentCharge
    File[yIndex,]$Mass <- (File[yIndex,]$Mass + 18.010575 + 1.007276466 * File[yIndex,]$FragmentCharge) / File[yIndex,]$FragmentCharge
    
    File$NewMatches <- str_c(File$FragmentType, " ", "(", File$FragmentNumber, ")", " ", "(", File$FragmentCharge, "+)", File$FragmentLoss)
    
    File <- File[,c("FileID", "Sequence", "Mass", "NewMatches")]
    File$Identifier <- str_c(File$FileID, "_", round(File$Mass, digits = 0))
    return(File)
  }
  
  ## Import data files
  
  setwd(HCD.Index)
  
  HCD.Index <- IndexCal("index.html")
  
  HCD.Index[,RoundPre := round(`Precursor.m/z.[Da]`, digits = 1)]
  #HCD.Index$RoundPre <- round(HCD.Index$`Precursor.m/z.[Da]`, digits = 1)
  
  ETD.HCD.ETD[,RoundPre := round(`m/z.[Da]`, digits = 1)]
  
  Index <- merge.data.frame(ETD.HCD.ETD, HCD.Index, by = "RoundPre", all.x = TRUE)
  
  Index <- data.table(Index)
  
  colnames(Index) <- str_replace_all(colnames(Index), "\\.x", "")
  colnames(Index) <- str_replace_all(colnames(Index), "\\.y", "\\.HCD")
  
  Index[, ScanDiff := `First.Scan` - `First.Scan.HCD`]%>%
    .[,SameScan := abs(ScanDiff) < 5]
  
  Index <- Index[order(ID, abs(ScanDiff))]%>%
    .[,Duplicated := duplicated(ID)]
  
  Index <- Index[Duplicated == FALSE]
  
  Index <- Index[,c(4:34,50)]
  
  Index$Modifications <- str_replace_na(Index$Modifications)
  
  Index$Identifier <- str_c(Index$Sequence, Index$Modifications, Index$Charge)
  
  Index <- arrange(Index, desc(XCorr))
  
  Index <- Index[!duplicated(Index$Identifier),]
  
  Index <- data.frame(Index)
  
  FragmentIndex <- FragCal(Index)
  
  # Here to calculate glycan mass
  Index$TNum <- str_count(Index$Modifications, pattern = "\\(Hex\\(1\\)HexNAc\\(1\\)\\)")
  
  Index$TnNum <- str_count(Index$Modifications, pattern = "\\(HexNAc\\)")
  
  Index$SugarMass <- 365.132196 * Index$TNum + 203.079373 * Index$TnNum
  
  Index$SugarLossPre <- Index$m.z..Da. * Index$Charge - Index$SugarMass - Index$Charge * 1.007276466
  
  Index <- subset(Index, is.na(File) == FALSE)
  
  FragList <- read.delim(Index$File[1], stringsAsFactors=FALSE)
  FragList$File <- Index$File[1]
  FragList <- FragList[0,]
  
  FragList <- adply(Index$File,1,SubFrag)
  
  setwd(wd)
  
  FragList <- ddply(FragList, .(File), transform, rank = rank(desc(i), ties.method = "first"))
  
  FragList <- subset.data.frame(FragList, FragList$matches != "" | rank <= 15)
  
  FragList$Sequence <- str_split_fixed(FragList$File, "[_/]", 3)[,2]
  
  FragList$Identifier <- str_c(FragList$File, "_", round(FragList$m.z, digits = 0))
  
  FragList$X1 <- 1:length(FragList[,1])
  
  ComList <- merge.data.frame(FragList, FragmentIndex[,3:5], by = "Identifier", all.x = TRUE)
  
  ComList$deltaMass <- ComList$m.z - ComList$Mass
  
  ComList <- arrange(ComList, X1, abs(deltaMass))
  
  ComList <- ComList[!duplicated(ComList$X1),]
  
  massIndex <- is.na(ComList$deltaMass) == FALSE & abs(ComList$deltaMass) >= 0.01
  
  ComList[massIndex,c("Mass", "NewMatches", "deltaMass")] <- NA
  
  massIndex <- is.na(ComList$Mass) == FALSE
  
  ComList[massIndex, c("m.z", "matches")] <- ComList[massIndex, c("Mass", "NewMatches")]
  
  ComList$matches <- str_split_fixed(ComList$matches, ",", 2)[,1]
  
  ComList$matches <- str_replace_all(ComList$matches, "\\[M\\+[1-9]H\\]", str_c("y ", "(", str_length(ComList$Sequence), ")"))
  
  #To subset list with matches
  ComList <- subset.data.frame(ComList, matches != "")
  
  ComList <- ddply(ComList, .(File), transform, newrank = rank(desc(i), ties.method = "first"))
  
  #To normalize the prediction
  CLP <- ComList[, c("File", "newrank", "i")]
  
  CLP <- dcast(CLP, File ~ newrank, value.var = 'i')
  
  CLP[is.na(CLP)] <- 0
  
  CLP$max <- apply(CLP[,2:length(CLP)], 1, max)
  
  CLP <- CLP[,c("File", "max")]
  
  ComList <- merge(ComList, CLP, by = "File", all.x = TRUE)
  
  ComList$newPrediction <- 100 * ComList$i / ComList$max
  
  colnames(ComList)[8] <- "UpperSequence"
  
  ComList <- merge.data.frame(ComList, Index, by.x = "File", by.y = "File", all.x = TRUE)
  
  rm(CLP)
  
  RL <- RL[,c("PrecursorCharge", "ModifiedPeptide", "iRT", "UniProtIds", "Protein.Name", "ProteinDescription", "Organisms", "Genes", "FASTAName", "FASTAHeader", "Protein.Existence", "Sequence.Version")]
  
  RL <- unique(RL)
  
  RL$Identifier <- str_c(RL$ModifiedPeptide, RL$PrecursorCharge)
  
  ComList <- SeqRep(ComList)
  
  ComList$ModifiedPeptide <- str_replace_all(ComList$ModifiedPeptide, "\\(", "\\[")
  ComList$ModifiedPeptide <- str_replace_all(ComList$ModifiedPeptide, "\\)", "\\]")
  ComList$ModifiedPeptide <- str_replace_all(ComList$ModifiedPeptide, "\\[1\\]", "\\(1\\)")
  
  ComList$IntModifiedPeptide <- ComList$ModifiedPeptide
  
  ComList$IntModifiedPeptide <- str_replace_all(ComList$IntModifiedPeptide, "Hex\\(1\\)HexNAc\\(1\\)", "\\+365")
  ComList$IntModifiedPeptide <- str_replace_all(ComList$IntModifiedPeptide, "HexNAc2", "\\+203")
  ComList$IntModifiedPeptide <- str_replace_all(ComList$IntModifiedPeptide, "Carbamidomethyl", "\\+57")
  ComList$IntModifiedPeptide <- str_replace_all(ComList$IntModifiedPeptide, "Oxidation", "\\+16")
  
  ComList$LabeledPeptide <- ComList$ModifiedPeptide
  
  ComList <- FragMatch(ComList)
  
  ComList[str_detect(ComList$FragmentLossType, "") == FALSE,]$FragmentLossType <- "noloss"
  
  ComList$IsProteotypic <- "TRUE"
  
  ComList$ExcludeFromAssay <- FALSE
  
  ComList$Database <- "sp"
  
  ComList$Identifier <- str_c(ComList$ModifiedPeptide, ComList$Charge)
  
  ComList$Sequence <- toupper(ComList$Sequence)
  
  ComList <- merge.data.frame(ComList, RL, by = "Identifier", all.x = TRUE)
  
  
  ComList <- ComList[ ,c("Charge", "IntModifiedPeptide", "ModifiedPeptide.x", "Sequence", "iRT", "UniProtIds", "IsProteotypic", "LabeledPeptide", "m.z..Da.", "FragmentLossType", "FragmentNumber", "FragmentType", "FragmentCharge", "m.z", "newPrediction", "ExcludeFromAssay", "Database", "UniProtIds", "UniProtIds", "Protein.Name", "ProteinDescription", "Organisms", "Genes", "Protein.Existence", "Sequence.Version", "FASTAName", "FASTAHeader")]
  
  colnames(ComList) <- c("PrecursorCharge", "IntModifiedPeptide", "ModifiedPeptide", "StrippedPeptide", "iRT", "BGSInferenceId", "IsProteotypic", "LabeledPeptide", "PrecursorMz", "FragmentLossType", "FragmentNumber", "FragmentType", "FragmentCharge", "FragmentMz", "RelativeIntensity", "ExcludeFromAssay", "Database", "ProteinGroups", "UniProtIds", "Protein Name", "ProteinDescription", "Organisms", "Genes", "Protein Existence", "Sequence Version", "FASTAName", "FASTAHeader")
  
  ComList <- ddply(ComList, .(PrecursorCharge, StrippedPeptide), transform, rank = rank(desc(RelativeIntensity), ties.method = "first"))
  
  ComList <- subset.data.frame(ComList, is.na(ComList$iRT) == FALSE)
  
  return(ComList)
}

