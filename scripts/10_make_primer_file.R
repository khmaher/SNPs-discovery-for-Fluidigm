rm(list = ls())

library(dplyr)
library(tidyr)
library(stringr)

# set working directory to where you have saved your files
setwd("")


# change to name of your input fasta file
fasta<-read.table("108.SNPs.fasta")

df<-data.frame(Name = fasta$V1[c(TRUE, FALSE)], 
           Sequence = fasta$V1[c(FALSE, TRUE)])

df2<-df

df2$Sequence<-toupper(df2$Sequence)

df2$Sequence<-gsub("^(.{250})(.*)$",         # Apply gsub
     "\\1]\\2",
     df2$Sequence)
df2$Sequence<-gsub("^(.{249})(.*)$",         # Apply gsub
                   "\\1[\\2",
                   df2$Sequence)

df2$Name<-gsub(">",         # Apply gsub
               "",
               df2$Name)

# change to name of your info file
info<-read.table("108.SNPs_info.txt")

info$V1021<-info$V2-250
info$V102<- info$V2+250
  
info$V103<-paste(info$V1,":", info$V1021, "-",info$V102, sep = "")
info$V104<-paste(info$V3, info$V4,sep = "")
info$V105<-paste(info$V1,"_", info$V2, "_", info$V104, sep = "")


info<-select(info, c('V103','V104','V105'))

df3<-merge(info, df2, by.x="V103", by.y = "Name")

test<-data.frame(do.call("rbind", strsplit(df3$Sequence, "[", 2)))
test2<-data.frame(do.call("rbind", strsplit(test$X2, "]", 2)))
test3<-as.data.frame(cbind(test$X1, df3$V104, test2$X2))
test3$V2<-gsub("^(.{1})(.*)$",         # Apply gsub
            "\\1/\\2",
            test3$V2)

df3$Sequence<-paste(test3$V1, "[", test3$V2, "]", test3$V3, sep = "")

df4<-df3[,3:4]

colnames(df4)<-c("Name", "Sequence")

# change to name of your output file
write.table(df4, "Fluidigm_SNP_Type_Assays_Target_Import_MonkParakeet_v1.txt", sep = "\t", quote=FALSE, row.names = F, col.names = T)


