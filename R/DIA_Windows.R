library(tidyverse)
library(magrittr)
library(data.table)
library(readxl)
library(RColorBrewer)

DIA_windows <- read_excel("data/DIA_windows.xlsx")

DIA_windows$Height <- 100

DIA_windows$Rank <- rank(DIA_windows$mz)

ggplot(DIA_windows, aes(x = mz, y = Height, fill = factor(Width))) + 
  geom_bar(stat = "Identity", width = DIA_windows$Width - 1, alpha = 0.7) +  theme_bw() +
  scale_fill_manual(values = brewer.pal(3, "Accent")) + theme(legend.position = "none") + 
  geom_text(label = DIA_windows$mz, vjust = DIA_windows$Rank)


Hek_WT <- read.delim("data/HepG2_WT_Jacalin.xls", stringsAsFactors = FALSE)%>%
  data.table()%>%
  .[ModifiedPeptide %like% "Hex"]


Peptides <- Hek_WT[,.(ModifiedPeptide, PrecursorMz)]%>%
  unique()

ggplot(Peptides, aes(x = PrecursorMz)) + geom_histogram(binwidth = 10) + theme_bw() + 
  xlim(400, 1200)
