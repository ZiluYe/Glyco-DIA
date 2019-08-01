
library(plyr)
library(tidyverse)
library(data.table)
library(magrittr)
library(readxl)

NTermLb <- function(OriginalLibrary){
  OL <- copy(OriginalLibrary)
  OL[,HexNAcNum := str_count(`ModifiedPeptide`, "HexNAc")]%>%
    .[,HexNum := str_count(`ModifiedPeptide`, "Hex") - `HexNAcNum`]
  
  OL[,IntModifiedPeptide := str_replace_all(`IntModifiedPeptide`, "\\[\\+203\\]", "")]%>%
    .[,IntModifiedPeptide := str_replace_all(`IntModifiedPeptide`, "\\[\\+365\\]", "")]
  
  OL[HexNAcNum != 0, IntModifiedPeptide := str_replace(`IntModifiedPeptide`, "_", 
                                                       str_c("_[+", (162 * `HexNum` + 203 * `HexNAcNum`), "]"))]
  
  OL[,ModifiedPeptide := str_replace_all(`ModifiedPeptide`, "\\[HexNAc]", "")]%>%
    .[,ModifiedPeptide := str_replace_all(`ModifiedPeptide`, "\\[Hex\\(1\\)HexNAc\\(1\\)]", "")]
  
  OL[HexNAcNum == 1 & HexNum == 0, ModifiedPeptide := str_replace(`ModifiedPeptide`, "_", "_[HexNAc]")]
  OL[HexNum != 0, ModifiedPeptide := str_replace(`ModifiedPeptide`, "_", 
                                                 str_c("_[Hex(", HexNum, ")HexNAc(", HexNAcNum, ")]"))]
  
  OL[, LabeledPeptide := ModifiedPeptide]
  
  return(OL[,1:27])
}

SiteSpecificLb <- function(ETD.HCD.HCD, ETD.HCD.ETD, ETD.HCD.HCD.LB){
  PSM.ETD.HCD.HCD <- PSMSum(ETD.HCD.HCD)
  PSM.ETD.HCD.ETD <- PSMSum(ETD.HCD.ETD)
  
  PSM.ETD.HCD.HCD.LB <- LbSum(ETD.HCD.HCD.LB)
  
  HCD <- PSM.ETD.HCD.HCD[[1]]
  ETD <- PSM.ETD.HCD.ETD[[1]]
  
  HCD[, ID:= 1:dim(HCD)[1]]
  ETD[, ID:= 1:dim(ETD)[1]]
  
  LB <- PSM.ETD.HCD.HCD.LB[[1]]
  
  LB <- LB[,.(iRT, Identifier, SeqMz)]%>%
    unique()
  
  LB[, LBID:= 1:dim(LB)[1]]
  
  MergeHCD <- merge.data.frame(HCD, ETD[,.(Sequence, `RT [min]`, `First Scan`, Identifier, SeqMz, ID)], by = "SeqMz", all.x = TRUE)
  colnames(MergeHCD) <- str_replace_all(colnames(MergeHCD), "\\.x", "")
  colnames(MergeHCD) <- str_replace_all(colnames(MergeHCD), "\\.y", "\\.ETD")
  MergeHCD <- data.table(MergeHCD)
  
  MergeHCD[, ScanDiff := `First Scan` - `First Scan.ETD`]%>%
    .[,SameIdentifier := Identifier == Identifier.ETD]%>%
    .[,RTDiff := `RT [min]` - `RT [min].ETD`]%>%
    .[,SameScan := abs(ScanDiff) < 5]
  
  
  MergeHCD <- MergeHCD[order(ID, abs(ScanDiff))]%>%
    .[,Duplicated := duplicated(ID)]
  
  
  MergeHCD <- MergeHCD[Duplicated == FALSE]
  
  MergeHCD <- merge(MergeHCD, LB, by = "Identifier", all.x = TRUE)
  
  MergeHCD[,iRTRank := rank(iRT)]%>%
    .[,RTRank := rank(`RT [min]`)]
  
  MatchList <- MergeHCD[,.(Identifier, `RT [min]`, ID.ETD, Identifier.ETD, SameIdentifier, iRT, LBID)]%>%
    unique()%>%
    .[is.na(LBID) == FALSE & is.na(Identifier.ETD) == FALSE]
  
  
  fit <- lm(iRT~`RT [min]`,data=MatchList)
  
  MatchList[, NewiRT := fit[["coefficients"]][["(Intercept)"]] + fit[["coefficients"]][["`RT [min]`"]] * `RT [min]`]%>%
    .[,Diff := NewiRT - iRT]
  
  newLB <- merge.data.frame(PSM.ETD.HCD.HCD.LB[[1]], MatchList, by = "Identifier", all.x = TRUE)
  colnames(newLB) <- str_replace_all(colnames(newLB), "\\.x", "")
  newLB <- data.table(newLB)
  
  newLB[LBID > 0, Identifier := Identifier.ETD]%>%
    .[LBID > 0, ModifiedPeptide := str_replace_all(Identifier.ETD, "_[2-7]", "_")]%>%
    .[LBID > 0, iRT := NewiRT]
  
  newLB$IntModifiedPeptide <- newLB$ModifiedPeptide
  
  newLB$IntModifiedPeptide <- str_replace_all(newLB$IntModifiedPeptide, "Hex\\(1\\)HexNAc\\(1\\)", "\\+365")
  newLB$IntModifiedPeptide <- str_replace_all(newLB$IntModifiedPeptide, "HexNAc2", "\\+203")
  newLB$IntModifiedPeptide <- str_replace_all(newLB$IntModifiedPeptide, "Carbamidomethyl", "\\+57")
  newLB$IntModifiedPeptide <- str_replace_all(newLB$IntModifiedPeptide, "Oxidation", "\\+16")
  
  newLB$LabeledPeptide <- newLB$ModifiedPeptide
  
  newLB1 <- newLB[LBID > 0, c(2:28)]
  
  newLB2 <- newLB[is.na(LBID), c(2:28)]%>%
    NTermLb()
  
  newLB <- rbind.data.frame(newLB1, newLB2)
  
  return(newLB)
}
