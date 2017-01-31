# KT002_cleanup_all_timepoints

```r
rm(list=ls())
setwd("~/Desktop/Microarray Analyses/Kevin's array/151118 data/") #set working directory (folder where Apac array files are stored)
```

##########################################################################################
##load data

```r
print(load("./Data/Processed Data/0001.KTarray.Apac.cleaned_X1.RData"))
print(load("./Data/Processed Data/0001.KTarray.Apac.cleaned_X2.RData"))
print(load("./Data/Processed Data/0001.KTarray.Apac.cleaned_X3.RData"))

files <- sort(list.files("./Data/Processed Data")[grep("import_data.RData", list.files("./Data/Processed Data"))], decreasing=T)[1]
print(load(paste("./Data/Processed Data/", files, sep=""))) 
d1$spot[grep("std._", d1$spot)] <- gsub("std", "std0", d1$spot[grep("std._", d1$spot)])
```

##########################################################################################
##merge MEAN data 

```r
colnames(X1_raw_mean) <- paste("X1", colnames(X1_raw_mean), sep="_")
colnames(X2_raw_mean) <- paste("X2", colnames(X2_raw_mean), sep="_")
colnames(X3_raw_mean) <- paste("X3", colnames(X3_raw_mean), sep="_")

raw_mean <- cbind(X1_raw_mean,X2_raw_mean)
raw_mean <- cbind(raw_mean,X3_raw_mean)
```

##########################################################################################
##merge MEDIAN data 

```r
colnames(X1_raw_median) <- paste("X1", colnames(X1_raw_median), sep="_")
colnames(X2_raw_median) <- paste("X2", colnames(X2_raw_median), sep="_")
colnames(X3_raw_median) <- paste("X3", colnames(X3_raw_median), sep="_")

raw_median <- cbind(X1_raw_median,X2_raw_median)
raw_median <- cbind(raw_median,X3_raw_median)
```

##########################################################################################
##create pheSera dataframe with info on samples

```r
pheSera <- d1[duplicated(paste(d1$study, d1$sampleID, d1$slideN, sep="."))==F,c("study", "sampleID", "slideN", "subarrayN")]
names(pheSera) <- c("study", "Sample", "Slide", "Pad")
pheSera$Sample[pheSera$study=="Apac_X2" & pheSera$Slide==59] <- paste(pheSera$Sample[pheSera$study=="Apac_X2" & pheSera$Slide==59],"b", sep="")
```

##########################################################################################
##create annAnti dataframe with info on antigens

```r
antigens <- read.csv("./Data/Processed Data/150702.Apac.antigen_descriptions.csv")
antigens$antigen <- gsub("_", "-", antigens$antigen) 
antigens$antigen[grep("lank", antigens$antigen)] <- "blank"         # "Blank" and "blank" now are both "blank"
antigens <- antigens[duplicated(antigens$antigen)==F,]
antigens$antigen <- gsub(" ", "0", antigens$antigen)
antigens$antigen <- gsub("std010", "std10", antigens$antigen)

annAnti <- d1[duplicated(d1$spot)==F, c("INDEX", "spot", "MR", "MC", "SR", "SC")]
annAnti$MR <- 1
colnames(annAnti) <- c("INDEX", "Unique.ID", "Array.Row", "Array.Column", "Spot.Row", "Spot.Column")
annAnti$Unique.ID <- gsub("15-34", "1534", annAnti$Unique.ID)
annAnti$Unique.ID <- gsub("15-38", "1538", annAnti$Unique.ID)
annAnti$Unique.ID <- gsub("17-34", "1734", annAnti$Unique.ID)
annAnti$Unique.ID <- gsub("52-42", "5242", annAnti$Unique.ID)
annAnti$antigen <- gsub("_.+$", "", annAnti$Unique.ID)

annAnti <- merge(annAnti, antigens, by="antigen", all.x=T)
annAnti <- annAnti[order(annAnti$Unique.ID), c("Unique.ID", "antigen", "description", "INDEX", "Array.Row", "Array.Column", "Spot.Row", "Spot.Column")]
annAnti$ID <- annAnti$Unique.ID
```

##########################################################################################
##cleanup dataframes

```r
str(pheSera)
pheSera$Slide <- as.numeric(pheSera$Slide)
pheSera$Pad <- as.numeric(pheSera$Pad)

str(annAnti)
annAnti$description <- as.character(annAnti$description)
annAnti$INDEX <- as.numeric(annAnti$INDEX)
annAnti$Array.Column <- as.numeric(annAnti$Array.Column)
annAnti$Spot.Row <- as.numeric(annAnti$Spot.Row)
annAnti$Spot.Column <- as.numeric(annAnti$Spot.Column)
```

##########################################################################################
##cleanup workspace and save

```r
ls()
rm(antigens, d1, files,X1_raw_mean,X1_raw_median,X2_raw_mean,X2_raw_median,X3_raw_mean,X3_raw_median)
save.image("./Data/Processed Data/0002.KTarray.Apac.cleaned.RData")
```
