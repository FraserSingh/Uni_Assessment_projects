#!/usr/bin/R
#B217754
#Code taken from FGT tutorial files except where stated by ££ symbols,
#where code has been amended or adapted from other sources.

############################
#PREP
############################
#load the required libraries for the tutorials...
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager") # https://www.rdocumentation.org/packages/BiocManager/versions/1.30.20
BiocManager::install("affy") #https://www.bioconductor.org/packages/release/bioc/html/affy.html
BiocManager::install("limma") #https://bioconductor.org/packages/release/bioc/html/limma.html 
library(limma)
library(affy) #Bioconductor affy library including ReadAffy
library(mouse4302.db)# load chip-specific annotation, https://www.bioconductor.org/packages/release/data/annotation/html/mouse4302.db.html
install.packages("scatterplot3d",repo="http://cran.ma.imperial.ac.uk")#Install then load the library for PCA. here we request a specific package from a specific archive...
library(scatterplot3d)
library(annotate) #https://www.bioconductor.org/packages/release/bioc/html/annotate.html

#Unix commands to make new folder for files, download dataset from GEO, rename it, extract the tar file and unzip the zip folders. URL retrieved by manually searching GEO for accession number, and copying link for tar file 
system(
    "mkdir GSE18581_analysis; cd GSE18581_analysis; 
    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE18581&format=file';
    mv 'index.html?acc=GSE18581&format=file' GSE18581.tar;
    tar -xvf GSE18581.tar; 
    gzip -d *")

setwd("./GSE18581_analysis") #ensuring R is using the right directory

#NOTE view target.csv file for info on the samples
#target file was made as stated in essay

# Load the target file into an AnnotatedDataFrame object
#££ Targetfile is saved as targetfile.txt, attached in submission.
adf<-read.AnnotatedDataFrame("../targetfile_B217754.txt",header=TRUE,as.is=TRUE,row.names=NULL) #££ adjusted to have no rownames
rownames(adf)<-adf$Sample #££ set rownames


#sort the groups by cell group
order(adf$Group,decreasing=T)# This gives a set of indexes that defines the order we need to sort the table elements in decreasing order by Group. We can use these as indexes to re-order the adf file
adf <-adf[order(adf$Group,decreasing=T),]

#££ Loading expression data from local CEL files (ReadAffy generates an object of type AffyBatch in the R workspace (here called myICAdata)).
myICAdata <- ReadAffy(filenames=adf$Filename) 

############################
#Quality control of raw data
############################
#££ Manually gathering present-absent values. 
#Adapted from lecture 7 'Design'.
#Uses mas5, not log scaled expression data 
#so can be performed here before log converison later
    calls<-mas5calls(myICAdata)
    calls<-exprs(calls)
    colnames(calls)<-adf$Name #££ rename from filenames
    absent<-colSums(calls=='A')
    present<-colSums(calls=='P')

total_calls <-  # Define the total number of calls for each sample
percent_present <- round(((present / (absent + present)) * 100), digits=2) # Calculate the percentage of present calls for each sample
# Combine the original data and the new column and make a dataframe, adapted from https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/data.frame
table_data <- cbind(absent, present, percent_present)
presence_table <- data.frame(Absent=table_data[,1], Present=table_data[,2], Percent_Present=table_data[,3])
write.table(presence_table, file="Presence_metrics.csv", sep=",")

# Quality control plots
#££ plotted histogram to file, adapted from https://r-coder.com/add-legend-r/ and https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
plotting_colours=c("#FF0000", "#117733", "#FFA500", "#999933","#0000FF","#FF00FF")
png("GSE18581_QC_hist.png")
hist(myICAdata, main="Gene expression distribution in experiment GSE18581",col=plotting_colours)
legend("topright", legend = adf$Name, title = "Samples", fill=line_colours) #add specific location and parameters for the legend
dev.off()

# And a boxplot with different colour per sample group
#££ plotted with parameters to present nicely
#Code adapted from:
    #https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/mtext
    #https://www.statology.org/par-function-in-r/
    #https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/legend
png("GSE18581_QC_boxplot.png")
par(mar = c(5, 4, 4, 6) + 0.1) # adjust right margin to make room for legend
boxplot(myICAdata, col = adf$Colour, ylab = "Expression level", main = "Non-normalised expression values", names = adf$Name,las=2)
legend("topright", legend = unique(adf$Group), col = unique(adf$Colour), title = "Group", pch = 15, ncol = 2, xpd = TRUE, bty = "n",inset = c(-0.25,-0.12)) #add specific location and parameters for the legend
dev.off()


############################
# Data normalisation and plotting
############################

#In the overall Affymetrix analysis workflow AffyBatch files are commonly converted to ExpressionSet by normalisation with functions such as rma and mas5. 
#The normalisation of microarray data is a necessary step in order to account for variation that arises from the way the samples have been generated in the lab

# Normalise the data using RMA
eset <- rma(myICAdata)

# To obtain a matrix of the expression values, use exprs()
values <- exprs(eset)
colnames(values)<-rownames(adf)
#££ Boxplot to observe the results of normalisation
png("GSE18581_normalised_boxplot.png")
par(mar = c(5, 4, 4, 6) + 0.1) 
boxplot(values, col = adf$Colour, ylab = "Expression level (log2)", main = "Normalised expression values", las = 2, names = adf$Name)
legend("topright", legend = unique(adf$Group), col = unique(adf$Colour), title = "Group", pch = 15, ncol = 2, xpd = TRUE, bty = "n",inset = c(-0.25,-0.12)) #add specific location and parameters for the legend
dev.off()

#MA plot of all
# Note that in this code the pm() command is used to plot only the perfect match
# oligonucleotides. 
png("All_mvapairs_text.png", width=1000, height=1100)
mva.pairs(values,main='MVA for all samples (normalised)',label=adf$Name)
dev.off()

#MVA plot for further investigation of samples LCa.2 and HCa.2
png("MVALCa2HCa2.png")
mva.pairs(values[, c(2,5)],main='MVA for LCa.2 and HCa.2 (normalised)',label=adf$Name)
dev.off()

# # To facilitate interpretation, let’s replace the columns header,currently  displaying the filename, to show the name of each sample (if you have a targets file)
# colnames(values) <- first_order #££ using the object defined earlier

# Performs hierarchical clustering with average linkage based on
# Pearson’s Correlation Coefficient
hc<-hclust(as.dist(1-cor(values, method="pearson")), method="average")
png("GSE18581_heir_clust.png");plot(hc, label = adf$Name);dev.off() #££ neater plotting

#Plot PCA
# Perform PCA
pca <- prcomp(t(values), scale=T)

#Plotting PCA of normalised data, with PC3 as the focus
png("GSE18581_PCA_PC3.png")
s3d<-scatterplot3d(pca$x[,c(3,2,1)], pch=19,color=plotting_colours, main='PCA plot of all samples, PC3 as focus') #££ indexing changed
s3d.coords <- s3d$xyz.convert(pca$x[,c(3,2,1)])
text(s3d.coords$x, s3d.coords$y, labels = adf$Name,pos= 3,offset = 0.5)
dev.off()

############################
# Fold Filtering
############################

#££ Matrix of the expression values held in 'values' object
#RMA outputs log2 data while MAS5 outputs linear data
#To convert from log…
exprsvals10 <-2^values

# Check Sample Name Order
# For vectors of the names of the samples (essential for checking the order of the
# samples are read in) and to obtain the probe IDs of the genes use the functions
# sampleNames() and probeNames() respectively.
#check order of sample names
mysamples <- sampleNames(eset)
#it is useful to obtain a vector of ProbeIDs here
probesets <- probeNames(myICAdata)


# Build a Fold Change Vector
# Here we output the final summary data containing the expression means of each replicate
# group and the fold changes. We first use the apply method to apply the mean function to
# the selection of columns from each replicate group. We then build up each mean list.

# The cbind function is then used to add all the columns together into a table in order that
# the summary data can be conveniently output. Note that here the mean is calculated from
# the non-log transformed values. 

#Calculate the means
#Note mean of the log is not the same as the log of the mean!!
#££ changed for my sample names
SI1.mean <- apply(exprsvals10[,c("GSM213045","GSM462207")],1,mean) 
SI2.mean <- apply(exprsvals10[,c("GSM213051","GSM462208")],1,mean) 
SI3.mean <- apply(exprsvals10[,c("GSM213057","GSM462209")],1,mean)
HCa.mean <- apply(exprsvals10[,c("GSM213045","GSM213051","GSM213057")],1,mean)
LCa.mean <- apply(exprsvals10[,c("GSM462207","GSM462208","GSM462209")],1,mean)
#££ Calculate fold change between low calcium diet (control) and high calcium (treatment)
HiLoCa <-LCa.mean / HCa.mean

#Build a summary table to hold all the data
all.data= cbind(SI1.mean,SI2.mean,SI3.mean,HCa.mean,LCa.mean,HiLoCa)
#write the table of means as an output
write.table(all.data,file="Group_means.csv", quote=F,sep=",",col.names=NA)

############################
#Statistical analysis on RMA normalised expression data using Limma
############################

#We then generate list of differentially expressed genes which will be annotated
#Annotate the Results with Gene Names   
#Establish annotation for MOE430v2 
#modified from #http://gettinggeneticsdone.blogspot.co.uk/2012/01/annotating-limma-results-with-gene.html
#Mapping from http://genomicsclass.github.io/book/pages/mapping_features.ht#ml

#build an annotation table
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID, "mouse4302.db")
Name <- as.character(lookUp(ID, "mouse4302.db", "GENENAME"))
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name,stringsAsFactors=F)
tmp[tmp=="NA"] <- NA #fix padding with NA characters
#assign as feature data of the current Eset
fData(eset) <- tmp

#Uses linear models to analyse microarray data, requiring either one or two matrices to be specified. 
#Design Matrix indicates which samples have been applied to each array. Created using function model.matrix()
#££ This experiment contained 6 chips of 2 sample types, each element takes a value of 1 if element of the factor corresponding to the row of the matrix has the value correpsonding to the column of the matricx

#Build the design matrix
design <- model.matrix(~-1+factor(c(1,1,1,2,2,2))) #In this case the response are the gene profiles and the model the formulae we use to predict this response. We are just specifying the model part of this formulae
colnames(design) <- c("LCa","HCa") #setting column names to memorable names makes constructing contrast matrix easier.

#Contrasts Matrix
#Specifies which comparisons are to be made between samples. 
# ££ This instructs Limma which comparisons to make
contrastmatrix <- makeContrasts(HCa-LCa,levels=design)
#Fit lm for expression matrix exprs, design matrix design and contrast matrix contrast.matrix
#Fit the model...
fit <- lmFit(eset, design)
#...and make the contrasts, moderating the t-statistic using the borrowed variance approach
fit2 <- contrasts.fit(fit, contrastmatrix)
fit2 <- eBayes(fit2)
# ££ for toptreat to use. "topTreat assumes that the fit has been processed by treat."-https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/toptable
fit_treat<-treat(fit) 

# #Writing Summary Table from topTable
# #A list of the top differentially expressed genes for the first comparison specified by the contrast matrix can be obtained using the topTable() command:
#Here coef specifies which column of the contrast matrix to use for the comparison and adjust specifies what statistical adjustments to use (in this case ‘false discovery rate’).
myresults_topTable <-topTable(fit2,coef=1, adjust="fdr",number=nrow(eset)) #££ FIXME: adding lfc to filter by fold change values gave nothing
#Sort by logFC
myresults_topTable<-myresults_topTable[order(myresults_topTable$logFC),]
#Save topTable_results
write.table(myresults_topTable,"GSE18581results_topTable.csv", sep=",")

# #Writing Summary Table from topTreat
myresults_topTreat <-topTreat(fit_treat,coef=1, adjust="fdr",number=nrow(eset)) #££ Uses toptreat, 
#Rank by p -https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/toptable
myresults_topTreat<-myresults_topTreat[order(myresults_topTreat$adj.P.Val),]
#topTreat does not remove genes but ranks genes by evidence that their log-fold-change exceeds lfc.
write.table(myresults_topTreat,"GSE18581results_topTreat.csv", sep=",")

# #Classify each gene according to the pair-wise comparisons specified in the contrast matrix and display the results.
# decideTests() to # compare patterns of expression across the contrasts, because this gives control of FDR across the different contrasts. 
#££ adapted from https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/decideTests
png("GSE18581_vennDiagram_decideTests_labelled.png")
clas_decideTests <- decideTests(fit_treat, p.value=0.05)
vennDiagram(clas_decideTests)
title(main="Overlap of genes differentially expressed in both levels of comparison",cex.main=1, sub="FDR<0.05")
par(mar=c(2,2,2,2))
dev.off() #££

############################
#Functional Enrichment analysis
############################

system("cp /shared_files/MSigDB/Mm.h.all.v7.1.entrez.rds .") #Use Hallmark Gene signatures
Mm.H <- readRDS("Mm.h.all.v7.1.entrez.rds")

############################
#Annotate expression data
############################

#Here we select from the annotation a number of keys with the primary key being PROBEID
res <- select(mouse4302.db, keys = rownames(eset), columns =c("ENTREZID", "ENSEMBL","SYMBOL"), keytype="PROBEID")
#find the index of each row of the expression set in the annotation object res
idx <- match(rownames(eset), res$PROBEID)
#Use the index to set the phenotypic data in the ExpressionSet
fData(eset) <- res[idx, ]
#Find all rows that don’t have an EntrezID and remove then
eset_t<-eset[is.na(fData(eset)$ENTREZID)==0,]

#Enrichment analysis using MRoast or Camera or Romer

# We first need to convert the lists of IDs in the signatures database into an index into the
# expression table- this allows efficient compute of the results. Then we run the analysis and
# save the results. 
#To run these commands we run mroast or camera or romer.
#convert to indexes
H.indices <- ids2indices(Mm.H,fData(eset_t)$ENTREZID)

#Find enrichment signatures in the data, run mroast
results_mroast <-mroast(eset_t,
 index=H.indices,
 design=design,
 contrast=contrastmatrix[,1],
 adjust.method = "BH")

#View the results
head(results_mroast) #££ I think we want this one, self-contained experiment

#Writing a Summary Table
write.table(results_mroast,"enrichment_mroast.csv",sep=",")


#Part2

#Generate volcano plot using fit2, which moderates the t-statistic using the borrowed variance
#'highlight=10' prints the name of the 10 top genes on the plot
#Code adapted from 
#https://rdrr.io/bioc/limma/man/volcanoplot.html
#https://statisticsglobe.com/increase-font-size-in-plot-in-r
png("GSE18581_volcanoPlot2.png", width=800, height=800)
par(cex=1.25) #Set label size 
volcanoplot(fit2, coef = 1, style = "p-value", #choose data to plot
highlight = 10, names = fit2$genes$Symbol, hl.col="blue", #choose data to highlight
xlab = "Log2 Fold Change", ylab = NULL, pch=16, #Set axis labelling
main="Volcano plot of all data", sub="(Top 10 differentially expressed genes labelled in blue)", #set graph titles
cex=0.25,cex.lab = 1,cex.axis = 1,cex.main = 1,cex.sub = 0.75) #set all other text size
dev.off() #££
