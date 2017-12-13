#! /usr/bin/env Rscript

rm(list=ls())

setwd("/ufrc/renne/sunantha.s/research/CLASH/Scripts")

#Read all CLASH data files with 5'miRNA-3'mRNA chimeras and assign column names

INDIR <- "../ExcelFiles/"
OUTDIR <- "../Binding/"

Human_mRNAs <- read.csv("/ufrc/renne/sunantha.s/research/Ensembl/HumanTrancriptome/Human_mRNAs.csv", stringsAsFactors = F)
Human_CDS <- Human_mRNAs[,c("Transcript_ID","cDNA_Length", "CDS_Start", "CDS_End")]

#Change the comparison column from factor to character (**check why**) and if loop to classify each mRNA fragment into 6 classes: 5'UTR,
#5'UTR-CDS, "CDS", "CDS-3'UTR" and "unknown". Write the output files to a different folder.

files <- list.files(path = INDIR ,pattern = "\\.csv$")

files_5miR <- c(grep("5mi3m", files, value=T), grep("5v3m", files, value=T))
files_3miR <- c(grep("5m3mi", files, value=T), grep("5m3v", files, value=T))

for (file in files_5miR) {
  
filename <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
filename$TranscriptID_3 <- as.character(filename$TranscriptID_3)

filename <- merge(filename, Human_CDS, by.x= "TranscriptID_3", by.y= "Transcript_ID", all.x=T)
  for (i in 1:nrow(filename)) {
    if(is.na(filename$Transcript_start_3[i]) || is.na(filename$Transcript_end_3[i]) || is.na(filename$CDS_Start[i]) || is.na(filename$CDS_End[i])) {
      filename$Binding_Region[i] <- "unknown"
    } else if (filename$Transcript_start_3[i] < filename$CDS_Start[i] && filename$Transcript_end_3[i] < filename$CDS_Start[i]) {
      filename$Binding_Region[i] <- "5'UTR"
    } else if (filename$Transcript_start_3[i] < filename$CDS_Start[i] && filename$Transcript_end_3[i] > filename$CDS_Start[i]) {
      filename$Binding_Region[i] <- "5'UTR-CDS"
    } else if (filename$Transcript_start_3[i] > filename$CDS_Start[i] && filename$Transcript_end_3[i] < filename$CDS_End[i]) {
      filename$Binding_Region[i] <- "CDS"
    } else if (filename$Transcript_start_3[i] < filename$CDS_End[i] && filename$Transcript_end_3[i] > filename$CDS_End[i]) {
      filename$Binding_Region[i] <- "CDS-3'UTR"
    } else if (filename$Transcript_start_3[i] > filename$CDS_End[i] && filename$Transcript_end_3[i] > filename$CDS_End[i]) {
      filename$Binding_Region[i] <- "3'UTR"
    }
  }
  write.csv(filename, paste0(OUTDIR, file), row.names=F)
}

for (file in files_3miR) {
  
  filename <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
  filename$TranscriptID_5 <- as.character(filename$TranscriptID_5)
  
  filename <- merge(filename, Human_CDS, by.x= "TranscriptID_5", by.y= "Transcript_ID", all.x=T)
  for (i in 1:nrow(filename)) {
    if(is.na(filename$Transcript_start_5[i]) || is.na(filename$Transcript_end_5[i]) || is.na(filename$CDS_Start[i]) || is.na(filename$CDS_End[i])) {
      filename$Binding_Region[i] <- "unknown"
    } else if (filename$Transcript_start_5[i] < filename$CDS_Start[i] && filename$Transcript_end_5[i] < filename$CDS_Start[i]) {
      filename$Binding_Region[i] <- "5'UTR"
    } else if (filename$Transcript_start_5[i] < filename$CDS_Start[i] && filename$Transcript_end_5[i] > filename$CDS_Start[i]) {
      filename$Binding_Region[i] <- "5'UTR-CDS"
    } else if (filename$Transcript_start_5[i] > filename$CDS_Start[i] && filename$Transcript_end_5[i] < filename$CDS_End[i]) {
      filename$Binding_Region[i] <- "CDS"
    } else if (filename$Transcript_start_5[i] < filename$CDS_End[i] && filename$Transcript_end_5[i] > filename$CDS_End[i]) {
      filename$Binding_Region[i] <- "CDS-3'UTR"
    } else if (filename$Transcript_start_5[i] > filename$CDS_End[i] && filename$Transcript_end_5[i] > filename$CDS_End[i]) {
      filename$Binding_Region[i] <- "3'UTR"
    }
  }
  write.csv(filename, paste0(OUTDIR, file), row.names=F)
}

#####################################################################################################################################

files <- list.files(OUTDIR, pattern= "\\.csv")

counting <- data.frame(matrix(NA, nrow=length(files), ncol=7))
colnames(counting) <- c("Sample","5'UTR", "5'UTR-CDS", "CDS", "CDS-3'UTR", "3'UTR", "Unknown")

count =1
for(file in files) {
  temp <- read.csv(paste0(OUTDIR,file))
  counting[count,1] <- file
  counting[count,2] <- sum(temp$Binding_Region == "5'UTR")
  counting[count,3] <- sum(temp$Binding_Region == "5'UTR-CDS")
  counting[count,4] <- sum(temp$Binding_Region == "CDS")
  counting[count,5] <- sum(temp$Binding_Region == "CDS-3'UTR")
  counting[count,6] <- sum(temp$Binding_Region == "3'UTR")
  counting[count,7] <- sum(temp$Binding_Region == "unknown")
  count= count+1
}
write.csv(counting, paste0(OUTDIR,"Counts/","BindingRegion.csv"), row.names=F)

#####################################################################################################################################