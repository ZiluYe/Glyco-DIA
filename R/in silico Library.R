
## To build in silico libraries from T-DIA and/or Tn-DIA library

library(stringr)
C <- 12
H <- 1.0078250321
O <- 15.9949146221
S <- 31.97207069
N <- 14.0030740052
proton <- 1.007276466
electron <- 0.00054857990943

## Use test data as Library

Library <- read.delim("test/test_Tn.xls", stringsAsFactors = FALSE)


Library$TnNum <- str_count(Library$IntModifiedPeptide, "203")
Library$TNum <- str_count(Library$IntModifiedPeptide, "365")
Library <- subset.data.frame(Library, TnNum + TNum != 0)

TnLb <- Library
TnLb$IntModifiedPeptide <- str_replace_all(TnLb$IntModifiedPeptide, "365", "203")
TnLb$ModifiedPeptide <- str_replace_all(TnLb$ModifiedPeptide, "Hex\\(1\\)HexNAc\\(1\\)", "HexNAc")
TnLb$LabeledPeptide <- TnLb$ModifiedPeptide
TnLb$FragmentLossType <- str_replace_all(TnLb$FragmentLossType, "\\+C14\\+H23\\+N\\+O10", "+C8+H13+N+O5")
TnLb$PrecursorMz <- TnLb$PrecursorMz - (C*6 + H*10 + O*5) * TnLb$TNum/TnLb$PrecursorCharge
TnLb$TnNum <- str_count(TnLb$IntModifiedPeptide, "203")


TLb <- TnLb
TLb$IntModifiedPeptide <- str_replace_all(TLb$IntModifiedPeptide, "203", "365")
TLb$ModifiedPeptide <- str_replace_all(TLb$ModifiedPeptide, "HexNAc", "Hex(1)HexNAc(1)")
TLb$LabeledPeptide <- TLb$ModifiedPeptide
TLb$FragmentLossType <- str_replace_all(TLb$FragmentLossType, "\\+C8\\+H13\\+N\\+O5", "+C14+H23+N+O10")
TLb$PrecursorMz <- TLb$PrecursorMz + (C*6 + H*10 + O*5) * TLb$TnNum/TLb$PrecursorCharge

STnLb <- TnLb
STnLb$IntModifiedPeptide <- str_replace_all(STnLb$IntModifiedPeptide, "203", "494")
STnLb$ModifiedPeptide <- str_replace_all(STnLb$ModifiedPeptide, "HexNAc", "HexNAc(1)NeuAc(1)")
STnLb$LabeledPeptide <- STnLb$ModifiedPeptide
STnLb$FragmentLossType <- str_replace_all(STnLb$FragmentLossType, "\\+C8\\+H13\\+N\\+O5", "+C19+H30+N2+O13")
STnLb$PrecursorMz <- STnLb$PrecursorMz + (C*11 + H*17 + N + O*8) * STnLb$TnNum/STnLb$PrecursorCharge

mSTLb <- TnLb
mSTLb$IntModifiedPeptide <- str_replace_all(mSTLb$IntModifiedPeptide, "203", "656")
mSTLb$ModifiedPeptide <- str_replace_all(mSTLb$ModifiedPeptide, "HexNAc", "Hex(1)HexNAc(1)NeuAc(1)")
mSTLb$LabeledPeptide <- mSTLb$ModifiedPeptide
mSTLb$FragmentLossType <- str_replace_all(mSTLb$FragmentLossType, "\\+C8\\+H13\\+N\\+O5", "+C25+H40+N2+O18")
mSTLb$PrecursorMz <- mSTLb$PrecursorMz + (C*17 + H*27 + N + O*13) * mSTLb$TnNum/mSTLb$PrecursorCharge

dSTLb <- TnLb
dSTLb$IntModifiedPeptide <- str_replace_all(dSTLb$IntModifiedPeptide, "203", "947")
dSTLb$ModifiedPeptide <- str_replace_all(dSTLb$ModifiedPeptide, "HexNAc", "Hex(1)HexNAc(1)NeuAc(2)")
dSTLb$LabeledPeptide <- dSTLb$ModifiedPeptide
dSTLb$FragmentLossType <- str_replace_all(dSTLb$FragmentLossType, "\\+C8\\+H13\\+N\\+O5", "+C36+H57+N3+O26")
dSTLb$PrecursorMz <- dSTLb$PrecursorMz + (C*28 + H*44 + N*2 + O*21) * dSTLb$TnNum/dSTLb$PrecursorCharge

NakedLb <- TnLb
NakedLb$IntModifiedPeptide <- str_replace_all(NakedLb$IntModifiedPeptide, "\\[203\\]", "")
NakedLb$ModifiedPeptide <- str_replace_all(NakedLb$ModifiedPeptide, "\\(HexNAc\\)", "")
NakedLb$LabeledPeptide <- NakedLb$ModifiedPeptide
NakedLb$FragmentLossType <- str_replace_all(NakedLb$FragmentLossType, "\\+C8\\+H13\\+N\\+O5", "noloss")
NakedLb$PrecursorMz <- NakedLb$PrecursorMz - (C*8 + H*13 + N + O*5) * NakedLb$TnNum/NakedLb$PrecursorCharge

TnLb <- TnLb[,1:27]
TLb <- TLb[,1:27]
STnLb <- STnLb[,1:27]
mSTLb <- mSTLb[,1:27]
dSTLb <- dSTLb[,1:27]
NakedLb <- NakedLb[,1:27]





