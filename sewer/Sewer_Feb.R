##############################################
###  Corona Virus whole genome analysis    ###
###  Sewer in Israel in February           ###
###  February 2021                         ###
##############################################


#### Setup ####
setwd("/Users/netazuck/Documents/projects/Michal/CoronaVirus/data/Sewer/Sewer_routine/Feb2021")
library(msa); library(seqinr); library(ggplot2); library(stringdist); library(reshape2)
load("/Users/netazuck/Documents/projects/Michal/CoronaVirus/data/Sewer/refs/StartEndRegions.RData")
load("results/puList.RData")
load("results/puList_freqPos.RData")
load("results/mat_MutsByGene.RData")


#### align (augur's mafft)
cd /Users/netazuck/Documents/projects/Michal/CoronaVirus/data/Sewer/Sewer_routine/Feb2021

augur align \
--sequences alignment/Feb2021.fasta \
--reference-sequence /Users/netazuck/Documents/projects/Michal/CoronaVirus/data/Sewer/refs/REF_NC_045512.2.fasta \
--output alignment/Feb2021_aligned.fasta 



#### Import & arrange ALIGNED sequences ####
multAlign = read.fasta(file="alignment/Feb2021_aligned.fasta", seqtype="DNA")
multAlignNames = names(multAlign)

NucList = list()  #put sequences in list
for (i in 1:length(multAlign)){ NucList[[i]] = multAlign[[i]]}
names(NucList) = multAlignNames
save(NucList,file="results/NucList.RData")

NucMat = matrix(NA,nrow=length(NucList),ncol=length(NucList[[1]])); rownames(NucMat)=names(NucList); colnames(NucMat)=c(1:length(NucList[[1]]))
for (i in 1:length(NucList)){
  for (j in 1:length(NucList[[1]])){
    NucMat[i,j] = NucList[[i]][j]
  }
}
save(NucMat,file="results/NucMat.RData")


#### GET NUCLEOTIDE MIX IN CLADE POSITIONS (from BAM files - rsamtools) USING THIS! #####
library(Rsamtools)

## read bam files one by one into the pu_# object and the make a list of all (do this once and then load the puList.RData)
bamfile = "BAM/s-7356_ShafdanJan2021.mapped.sorted.bam"  
bf = BamFile(bamfile)
p_param = PileupParam(max_depth=10000,distinguish_strands=F,min_mapq=20)
#param = ScanBamParam(which=GRanges("chr17",IRanges(start=685640, end=685640)))
pu_7356 = pileup(bf,pileupParam=p_param)
#pu[pu$pos==11083,]

puList_new = list(puList[[1]],puList[[2]],puList[[3]],puList[[4]],puList[[5]],puList[[6]],      #beer-sheva
                  puList[[7]],puList[[8]],puList[[9]],puList[[10]],puList[[11]],puList[[12]],   #ashdod
                  puList[[13]],puList[[14]],puList[[15]],puList[[16]],puList[[17]],puList[[18]],#rahat
                  puList[[19]],puList[[20]],puList[[21]],puList[[22]],puList[[23]],             #jer-sor
                  puList[[24]],puList[[25]],puList[[26]],puList[[27]],puList[[28]],pu_7356,     #shafdan
                  puList[[30]],puList[[31]],puList[[32]],puList[[33]],puList[[34]],             #natanya
                  puList[[35]],puList[[36]],puList[[37]],puList[[38]],puList[[39]],puList[[40]],#haifa
                  puList[[41]],puList[[42]],puList[[43]],puList[[44]],puList[[45]],             #tzfat
                  puList[[46]],puList[[47]],puList[[48]],puList[[49]],puList[[50]])             #el-hamra
#puList[[51]])                                                                 #hebron

names(puList_new) = c("BeerSheva_aug","BeerSheva_sept","BeerSheva_oct","BeerSheva_nov","BeerSheva_dec","BeerSheva_jan",
                      "Ashdod_aug","Ashdod_sept","Ashdod_oct","Ashdod_nov","Ashdod_dec","Ashdod_jan",
                      "Rahat_aug","Rahat_sept","Rahat_oct","Rahat_nov","Rahat_dec","Rahat_jan",
                      "JerSor_aug","JerSor_sept","JerSor_oct","JerSor_dec","JerSor_jan",
                      "Shafdan_aug","Shafdan_sept","Shafdan_oct","Shafdan_nov","Shafdan_dec","Shafdan_jan",
                      "Natanya_aug","Natanya_sept","Natanya_oct","Natanya_nov","Natanya_jan",
                      "Haifa_aug","Haifa_sept","Haifa_oct","Haifa_nov","Haifa_dec","Haifa_jan",
                      "Tzfat_aug","Tzfat_sept","Tzfat_oct","Tzfat_dec","Tzfat_jan",
                      "ElHamra_aug","ElHamra_sept","ElHamra_oct","ElHamra_nov","ElHamra_dec")
puList = puList_new
save(puList,file="results/puList.RData")

# update puList w/ percentage for each position
load("results/NucMat.RData"); refSeq = NucMat[1,]; refSeq=toupper(refSeq)

for (pu in 1:length(puList)){  
  x=puList[[pu]]
  posTab = table(x$pos); max(posTab); table(table(x$pos))
  #pos1=as.integer(names(posTab[posTab==1])); pos2=as.integer(names(posTab[posTab==2]))
  #pos3=as.integer(names(posTab[posTab==3])); pos4=as.integer(names(posTab[posTab==4])); pos5=as.integer(names(posTab[posTab==5]))
  perc = rep(NA,nrow(x))
  
  for (i in 1:nrow(x)){
    xallpos = x[x$pos==x[i,]$pos,]  #get all positions
    sumpos = sum(xallpos$count) #get depth in this position
    perc[i] = (x[i,]$count / sumpos) * 100
  }
  x$perc = perc #add percentage
  
  #add ref seq to table
  ref = rep(NA,nrow(x)); diff = rep(0,nrow(x))
  for (i in 1:nrow(x)){ 
    ref[i] = refSeq[x[i,]$pos]  
    if (as.character(x[i,]$nucleotide) != ref[i]){ diff[i]=1 }
  }
  x$ref = ref
  x$diff = diff
  
  puList[[pu]] = x
  
}
save(puList,file="results/puList.RData")



#### GET FREQUENCY OF MUTAITONS ACROSS THE WHOLE GENOME  ####
load("results/puList.RData")
load("refs/StartEndRegions.RData")
StartEndRegions = StartEndRegions[-c(2:4),] #remove ORF1s
library(data.table)

puList_freqPos = puList
freq=1

for (i in 1:length(puList)){
  x = puList[[i]]
  
  #find frequently mutated positions
  freqPos = rep(NA,nrow(x)) 
  for (j in 1:nrow(x)){ if ((x[j,]$diff != 0) && (x[j,]$perc > freq) && (x[j,]$count >= 5)){ freqPos[j]=j }}
  freqPos=freqPos[!is.na(freqPos)]
  x = x[freqPos,]
  
  #add gene info
  gene = rep(NA,nrow(x))
  for (j in 1:nrow(x)){
    for (k in 1:nrow(StartEndRegions)){ 
      if (x[j,]$pos %inrange% c(StartEndRegions[k,]$start,StartEndRegions[k,]$end)) { gene[j]=as.character(StartEndRegions[k,]$segment); break }
    }
  }
  x$gene = gene
  
  puList_freqPos[[i]] = x
}
save(puList_freqPos,file="results/puList_freqPos.RData")



#### Search freq tables - manually 
posToSearch = 27972

for (i in 1:length(puList_freqPos)){
  pu = puList_freqPos[[i]]
  puName = names(puList_freqPos)[i]
  
  if (!is.na(match(posToSearch,pu$pos))){ 
    Sys.sleep(0.1)
    print(paste(puName,"  ",pu[match(posToSearch,pu$pos),]$perc))
    flush.console()
  }
}



#### HEATMAP FREQUENCY PLOT ####
library(gplots); library(RColorBrewer); Npalette = colorRampPalette(c("white","yellow","red"))(n = 299)

tab = read.delim2("tmp.txt"); rownames(tab)=tab[,1]; tab=tab[,-1] #import freq table from excel
hm = matrix(as.numeric(as.matrix(tab)),ncol=ncol(tab),nrow=nrow(tab)); rownames(hm)=rownames(tab); colnames(hm)=colnames(tab)
hm[hm<1]=0
#hm = hm[-5,]
hml = log2(hm);hml[hml==-Inf]=0  #transform to log2
hmmelt = reshape2::melt(hml)

pdf("figs/NovelVariants_sewer_heatmap.pdf")
ggplot(data = hmmelt, aes(x = Var1, y = Var2)) + 
  geom_tile(aes(fill = value), color = "white", size = 1) + 
  scale_fill_gradient(low = "cornsilk1", high = "tomato") + 
  theme_grey(base_size = 10) + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 8, colour = "gray50"),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size=8))
dev.off()



## PRINT FREQ TABLE (!!!)

#print all positions for all cities/months, sort by highest freq
for (i in 1:length(puList_freqPos)){
  write.table(puList_freqPos[[i]][,-1],file=paste("results/freqPos/",names(puList_freqPos)[i],".txt",sep=""),sep="\t",row.names=F,quote=F)
}




#### SEARCH FOR KNOWN VARIANTS ####

## B.1.1.7 - UK
UKmuts = c(21765:21770,21992:21994,23063,23271,23604,23709,24506,24914,27972,28048,28111,28280,28977,3267,5388,6954,11288,11289,11290,11291,11292,11293,11294,11295,11296)
res = matrix(NA,nrow=length(puList_freqPos),ncol=length(UKmuts))  
rownames(res) = names(puList_freqPos); colnames(res) = UKmuts
for (i  in 1:nrow(res)){ #for each entry
  p = puList_freqPos[[i]]
  for (j in 1:ncol(res)){  #for each mut pos
    if (!is.na(match(UKmuts[j],p$pos))){ res[i,j] = p[match(UKmuts[j],p$pos),]$perc }
  }
}
write.table(res,file="results/knownVariants/UKmuts_freq.txt",sep="\t")
#puList_freqPos[[1]][which(puList_freqPos[[1]]$pos=="23604"),]  #check P681H (if nuc is A)


## B.1.351 - SA
SAmuts = c(21801,22286:22294,22813,23012,23063,23664,25563,25904,26456,28887,1059,5230,10323,11288:11296)
res = matrix(NA,nrow=length(puList_freqPos),ncol=length(SAmuts))  
rownames(res) = names(puList_freqPos); colnames(res) = SAmuts
for (i  in 1:nrow(res)){ #for each entry
  p = puList_freqPos[[i]]
  for (j in 1:ncol(res)){  #for each mut pos
    if (!is.na(match(SAmuts[j],p$pos))){ res[i,j] = p[match(SAmuts[j],p$pos),]$perc }
  }
}
write.table(res,file="results/knownVariants/SAmuts_freq.txt",sep="\t")


## P.2 - MANAUS
MANAUSmuts = c(21614,21621,21638,21974,22132,22812,23012,23063,23525,24642,28167,28512,733,2749,3828,5648,11288:11296,12778,13860,17259)
res = matrix(NA,nrow=length(puList_freqPos),ncol=length(MANAUSmuts))  
rownames(res) = names(puList_freqPos); colnames(res) = MANAUSmuts
for (i  in 1:nrow(res)){ #for each entry
  p = puList_freqPos[[i]]
  for (j in 1:ncol(res)){  #for each mut pos
    if (!is.na(match(MANAUSmuts[j],p$pos))){ res[i,j] = p[match(MANAUSmuts[j],p$pos),]$perc }
  }
}
write.table(res,file="results/knownVariants/MANAUSmuts_freq.txt",sep="\t")

## B.1.525 - mystery variant
B1525muts = c(21717,21762,21765:21770,21992:21994,23012,23593,24224,24748,26767,27205:27207,28308,28699,28887,29543,1498,1807,2659,6285,8593,9565,11288:11296,14407,18171,20724)
res = matrix(NA,nrow=length(puList_freqPos),ncol=length(B1525muts))  
rownames(res) = names(puList_freqPos); colnames(res) = B1525muts
for (i  in 1:nrow(res)){ #for each entry
  p = puList_freqPos[[i]]
  for (j in 1:ncol(res)){  #for each mut pos
    if (!is.na(match(B1525muts[j],p$pos))){ res[i,j] = p[match(B1525muts[j],p$pos),]$perc }
  }
}
write.table(res,file="results/knownVariants/1525muts_freq.txt",sep="\t")

## A.23.1 - UGANDA
UGANDAmuts = c(21867,22033,22661,23401,23604,24097,28167,28378,4573,10747,11230,11266,11521,16575)
res = matrix(NA,nrow=length(puList_freqPos),ncol=length(UGANDAmuts))  
rownames(res) = names(puList_freqPos); colnames(res) = UGANDAmuts
for (i  in 1:nrow(res)){ #for each entry
  p = puList_freqPos[[i]]
  for (j in 1:ncol(res)){  #for each mut pos
    if (!is.na(match(UGANDAmuts[j],p$pos))){ res[i,j] = p[match(UGANDAmuts[j],p$pos),]$perc }
  }
}
write.table(res,file="results/knownVariants/UGANDAmuts_freq.txt",sep="\t")
#puList_freqPos[[1]][which(puList_freqPos[[1]]$pos=="23604"),]  #check P681R (if nuc is G)

## B.1.429 - CALIFORNIA
CALImuts = c(21600,22018,22917,12878,17014)
res = matrix(NA,nrow=length(puList_freqPos),ncol=length(CALImuts))  
rownames(res) = names(puList_freqPos); colnames(res) = CALImuts
for (i  in 1:nrow(res)){ #for each entry
  p = puList_freqPos[[i]]
  for (j in 1:ncol(res)){  #for each mut pos
    if (!is.na(match(CALImuts[j],p$pos))){ res[i,j] = p[match(CALImuts[j],p$pos),]$perc }
  }
}
write.table(res,file="results/knownVariants/CALImuts_freq.txt",sep="\t")

## B.1.526 - NEW YORK
NYmuts = c(21575,21846,22320,23012,23664,25517,27925,28869,28975,9867,11288:11296,16500,20262)
res = matrix(NA,nrow=length(puList_freqPos),ncol=length(NYmuts)) 
rownames(res) = names(puList_freqPos); colnames(res) = NYmuts
for (i  in 1:nrow(res)){ #for each entry
  p = puList_freqPos[[i]]
  for (j in 1:ncol(res)){  #for each mut pos
    if (!is.na(match(NYmuts[j],p$pos))){ res[i,j] = p[match(NYmuts[j],p$pos),]$perc }
  }
}
write.table(res,file="results/knownVariants/NYmuts_freq.txt",sep="\t")


## NEXTSTRAIN KNOWN CLADES
CLADEmuts = c(8782,14408,28144,14408,23403,28881,28882,28883,1059,25563)
res = matrix(NA,nrow=length(puList_freqPos),ncol=length(CLADEmuts))  
rownames(res) = names(puList_freqPos); colnames(res) = CLADEmuts
for (i  in 1:nrow(res)){ #for each city
  pu = puList_freqPos[[i]]
  for (j in 1:ncol(res)){  #for each mut pos
    if (!is.na(match(CLADEmuts[j],pu$pos))){ res[i,j] = pu[match(CLADEmuts[j],pu$pos),]$perc }
  }
}
write.table(res,file="results/knownVariants/CLADEmuts.txt",sep="\t")




#### COUNT # of mutations in each gene  ####

#load library size 
libSize = read.delim2("docs/LibrarySize.txt")

#calculate gene lengths
glengths = rep(NA,nrow(StartEndRegions)); names(glengths) = StartEndRegions$segment  
for (i in 1:length(glengths)) { glengths[i] = StartEndRegions[i,]$end - StartEndRegions[i,]$start }

mat_MutsByGene = matrix(NA,nrow=length(puList_freqPos),ncol=nrow(StartEndRegions))
rownames(mat_MutsByGene) = names(puList_freqPos); colnames(mat_MutsByGene) = StartEndRegions$segment
for (i  in 1:nrow(mat_MutsByGene)){ #for each city/month
  MutsByGene = rep(0,nrow(StartEndRegions)); names(MutsByGene) = StartEndRegions$segment
  x=table(puList_freqPos[[i]]$gene)
  MutsByGene[match(names(x),StartEndRegions$segment)] = x
  #MutsByGene = (MutsByGene/glengths)*10000  #normalize by gene length
  #MutsByGene = (MutsByGene / libSize[i,]$total.mapped.reads)*1000000
  #MutsByGene = (MutsByGene / libSize[i,]$total.reads)*1000000
  mat_MutsByGene[i,] = MutsByGene 
  
}

save(mat_MutsByGene,file="results/mat_MutsByGene.RData")
write.table(mat_MutsByGene,file="tmp.txt",sep="\t")

















### OLD
## Search freq tables for muts - automatically
UKmutsS = c(21765,21766,21767,21768,21769,21770,21991,21992,21993,23063,23271,23604,23709,24506,24914)
res = matrix(NA,nrow=length(puList),ncol=length(UKmutsS))  #UK, S gene only
rownames(res) = names(puList); colnames(res) = UKmutsS
for (i  in 1:nrow(res)){ #for each city
  pu = puList[[i]]
  if(!is.na(match("-",pu[pu$pos=="21765",]$nucleotide))) {res[i,1] = pu[pu$pos=="21765" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="21766",]$nucleotide))) {res[i,2] = pu[pu$pos=="21766" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="21767",]$nucleotide))) {res[i,3] = pu[pu$pos=="21767" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="21768",]$nucleotide))) {res[i,4] = pu[pu$pos=="21768" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="21769",]$nucleotide))) {res[i,5] = pu[pu$pos=="21769" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="21770",]$nucleotide))) {res[i,6] = pu[pu$pos=="21770" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="21991",]$nucleotide))) {res[i,7] = pu[pu$pos=="21991" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="21992",]$nucleotide))) {res[i,8] = pu[pu$pos=="21992" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="21993",]$nucleotide))) {res[i,9] = pu[pu$pos=="21993" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("T",pu[pu$pos=="23063",]$nucleotide))) {res[i,10] = pu[pu$pos=="23063" & pu$nucleotide=="T",]$perc}
  if(!is.na(match("A",pu[pu$pos=="23271",]$nucleotide))) {res[i,11] = pu[pu$pos=="23271" & pu$nucleotide=="A",]$perc}
  if(!is.na(match("A",pu[pu$pos=="23604",]$nucleotide))) {res[i,12] = pu[pu$pos=="23604" & pu$nucleotide=="A",]$perc}
  if(!is.na(match("T",pu[pu$pos=="23709",]$nucleotide))) {res[i,13] = pu[pu$pos=="23709" & pu$nucleotide=="T",]$perc}
  if(!is.na(match("G",pu[pu$pos=="24506",]$nucleotide))) {res[i,14] = pu[pu$pos=="24506" & pu$nucleotide=="G",]$perc}
  if(!is.na(match("C",pu[pu$pos=="24914",]$nucleotide))) {res[i,15] = pu[pu$pos=="24914" & pu$nucleotide=="C",]$perc}
}
write.table(res,file="results/freqPosMutTables/UKmutsS.txt",sep="\t")

SAmutsS = c(21801,22286:22294,22813,23012,23063,23664)
res = matrix(NA,nrow=length(puList),ncol=length(SAmutsS))  #SA, S gene only
rownames(res) = names(puList); colnames(res) = SAmutsS
for (i  in 1:nrow(res)){ #for each city
  pu = puList[[i]]
  if(!is.na(match("C",pu[pu$pos=="21801",]$nucleotide))) {res[i,1] = pu[pu$pos=="21801" & pu$nucleotide=="C",]$perc}
  if(!is.na(match("-",pu[pu$pos=="22286",]$nucleotide))) {res[i,2] = pu[pu$pos=="22286" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="22287",]$nucleotide))) {res[i,3] = pu[pu$pos=="22287" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="22288",]$nucleotide))) {res[i,4] = pu[pu$pos=="22288" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="22289",]$nucleotide))) {res[i,5] = pu[pu$pos=="22289" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="22290",]$nucleotide))) {res[i,6] = pu[pu$pos=="22290" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="22291",]$nucleotide))) {res[i,7] = pu[pu$pos=="22291" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="22292",]$nucleotide))) {res[i,8] = pu[pu$pos=="22292" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="22293",]$nucleotide))) {res[i,9] = pu[pu$pos=="22293" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("-",pu[pu$pos=="22294",]$nucleotide))) {res[i,10] = pu[pu$pos=="22294" & pu$nucleotide=="-",]$perc}
  if(!is.na(match("T",pu[pu$pos=="22813",]$nucleotide))) {res[i,11] = pu[pu$pos=="22813" & pu$nucleotide=="T",]$perc}
  if(!is.na(match("A",pu[pu$pos=="23012",]$nucleotide))) {res[i,12] = pu[pu$pos=="23012" & pu$nucleotide=="A",]$perc}
  if(!is.na(match("T",pu[pu$pos=="23063",]$nucleotide))) {res[i,13] = pu[pu$pos=="23063" & pu$nucleotide=="T",]$perc}
  if(!is.na(match("T",pu[pu$pos=="23664",]$nucleotide))) {res[i,14] = pu[pu$pos=="23664" & pu$nucleotide=="T",]$perc}
}
write.table(res,file="results/freqPosMutTables/SAmutsS.txt",sep="\t")





UKmuts = c(21765,21766,21767,21768,21769,21770,23063,23271,23604,23709,24506,24914,27972,28048,28111,28280,28977,3267,5388,6954,11288,11289,11290,11291,11292,11293,11294,11295,11296)
res = matrix(NA,nrow=length(puList_freqPos),ncol=length(UKmuts))  #UK
rownames(res) = names(puList_freqPos); colnames(res) = UKmuts
for (i  in 1:nrow(res)){ #for each city
  pu = puList_freqPos[[i]]
  for (j in 1:ncol(res)){  #for each mut pos
    if (!is.na(match(UKmuts[j],pu$pos))){ res[i,j] = pu[match(UKmuts[j],pu$pos),]$perc }
  }
}
write.table(res,file="results/freqPosMutTables/UKmuts.txt",sep="\t")

SAmuts = c(22289,22290,22291,22292,22293,22294,22295,22296,22297,23063,21801,22813,23012,23664,25563,25904,28887,26456,1059,5230,10323,11288,11289,11290,11291,11292,11293,11294,11295,11296)
res = matrix(NA,nrow=length(puList_freqPos),ncol=length(SAmuts))  #SA
rownames(res) = names(puList_freqPos); colnames(res) = SAmuts
for (i  in 1:nrow(res)){ #for each city
  pu = puList_freqPos[[i]]
  for (j in 1:ncol(res)){  #for each mut pos
    if (!is.na(match(SAmuts[j],pu$pos))){ res[i,j] = pu[match(SAmuts[j],pu$pos),]$perc }
  }
}
write.table(res,file="results/freqPosMutTables/SAmuts.txt",sep="\t")

BRZmuts = c(23012,28253,28628,28975,100,29754,11288,11289,11290,11291,11292,11293,11294,11295,11296)
res = matrix(NA,nrow=length(puList_freqPos),ncol=length(BRZmuts))  #BRZL
rownames(res) = names(puList_freqPos); colnames(res) = BRZmuts
for (i  in 1:nrow(res)){ #for each city
  pu = puList_freqPos[[i]]
  for (j in 1:ncol(res)){  #for each mut pos
    if (!is.na(match(BRZmuts[j],pu$pos))){ res[i,j] = pu[match(BRZmuts[j],pu$pos),]$perc }
  }
}
write.table(res,file="results/freqPosMutTables/BRZmuts.txt",sep="\t")

CLADEmuts = c(8782,14408,28144,14408,23403,28881,28882,28883,1059,25563)
res = matrix(NA,nrow=length(puList_freqPos),ncol=length(CLADEmuts))  #BRZL
rownames(res) = names(puList_freqPos); colnames(res) = CLADEmuts
for (i  in 1:nrow(res)){ #for each city
  pu = puList_freqPos[[i]]
  for (j in 1:ncol(res)){  #for each mut pos
    if (!is.na(match(CLADEmuts[j],pu$pos))){ res[i,j] = pu[match(CLADEmuts[j],pu$pos),]$perc }
  }
}
write.table(res,file="results/freqPosMutTables/CLADEmuts.txt",sep="\t")






## SEARCH FOR PATTERNS IN FREQ TABLES

sft = puList_freqPos[38:42] #puList_freqPos[34:37] #puList_freqPos[29:33] #puList_freqPos[25:28] #puList_freqPos[20:24] #puList_freqPos[11:15] #puList_freqPos[6:10] #puList_freqPos[1:5]
pos5 = Reduce(intersect, list(sft[[1]]$pos,sft[[2]]$pos,sft[[3]]$pos,sft[[4]]$pos,sft[[5]]$pos))  #all pos
pos4 = unique(c(Reduce(intersect, list(sft[[1]]$pos,sft[[2]]$pos,sft[[3]]$pos,sft[[4]]$pos)),
                Reduce(intersect, list(sft[[1]]$pos,sft[[2]]$pos,sft[[3]]$pos,sft[[5]]$pos)),
                Reduce(intersect, list(sft[[2]]$pos,sft[[3]]$pos,sft[[4]]$pos,sft[[5]]$pos)),
                Reduce(intersect, list(sft[[1]]$pos,sft[[2]]$pos,sft[[4]]$pos,sft[[5]]$pos))))
pos3 = unique(c(Reduce(intersect, list(sft[[1]]$pos,sft[[2]]$pos,sft[[3]]$pos)),Reduce(intersect, list(sft[[1]]$pos,sft[[2]]$pos,sft[[4]]$pos)),
                Reduce(intersect, list(sft[[1]]$pos,sft[[2]]$pos,sft[[5]]$pos)),Reduce(intersect, list(sft[[1]]$pos,sft[[3]]$pos,sft[[4]]$pos)),
                Reduce(intersect, list(sft[[1]]$pos,sft[[3]]$pos,sft[[5]]$pos)),Reduce(intersect, list(sft[[1]]$pos,sft[[4]]$pos,sft[[5]]$pos)),
                Reduce(intersect, list(sft[[2]]$pos,sft[[3]]$pos,sft[[4]]$pos)),Reduce(intersect, list(sft[[2]]$pos,sft[[3]]$pos,sft[[5]]$pos)),
                Reduce(intersect, list(sft[[2]]$pos,sft[[4]]$pos,sft[[5]]$pos)),Reduce(intersect, list(sft[[3]]$pos,sft[[4]]$pos,sft[[5]]$pos))))
posLAST = Reduce(intersect, list(sft[[4]]$pos,sft[[5]]$pos))
freqPos = unique(pos3,pos4,pos5,posLAST); write.table(freqPos,file="tmp.txt",sep="\t")

#freqPos = unique(c(Reduce(intersect, list(sft[[1]]$pos,sft[[2]]$pos,sft[[3]]$pos,sft[[4]]$pos)),
#                   Reduce(intersect, list(sft[[1]]$pos,sft[[2]]$pos,sft[[3]]$pos)),Reduce(intersect, list(sft[[1]]$pos,sft[[2]]$pos,sft[[4]]$pos)),
#                   Reduce(intersect, list(sft[[2]]$pos,sft[[3]]$pos,sft[[4]]$pos)),Reduce(intersect, list(sft[[1]]$pos,sft[[3]]$pos,sft[[4]]$pos)),
#                   Reduce(intersect, list(sft[[3]]$pos,sft[[4]]$pos))))


#### DETERMINE CLADE MIXES  ####  
load("results/puList.RData")


### ALL POSITIONS

## subset clade-positions
cladePos = c(8782,14408,28144,23403,28881,28882,28883,1059,25563,1397,11083,23929,13730,18877,15324,20268)
puList_cladesPos = puList

for (pu in 1:length(puList)){
  x =puList[[pu]]
  x = x[x$pos %in% cladePos,]
  puList_cladesPos[[pu]] = x
}
save(puList_cladesPos,file="results/puList_cladesPos.RData")



### ONE POSITION PER CLADE

## subset clade-positions
cladePos1 = c(8782,14408,23403,28881,1059); puList_cladesPos1 = puList
cladePos2 = c(8782,23403,28882,25563); puList_cladesPos2 = puList

for (pu in 1:length(puList)){
  x =puList[[pu]]
  x = x[x$pos %in% cladePos2,]
  puList_cladesPos2[[pu]] = x
}
save(puList_cladesPos2,file="results/puList_cladesPos2.RData")

# write to file
for (i in 1:length(puList_cladesPos2)){ 
  fname = paste("results/clades_pos2/",names(puList_cladesPos2)[[i]],".txt",sep="")
  write.table(puList_cladesPos2[[i]][,-1],file=fname,sep="\t",quote=F,row.names=F)
}


## plot summary  
#https://www.r-graph-gallery.com/48-grouped-barplot-with-ggplot2.html

cladetab = read.delim2("results/clades_pos2/clades2Summary.txt")
cladetab$per = as.numeric(as.matrix(cladetab$per)) #aug$per=as.numeric(levels(aug$per))[aug$per]

#by city
ggplot(cladetab[cladetab$city=="ElHamra",], aes(fill=clade, y=per, x=month)) + 
  geom_bar(position="fill", stat="identity") + ggtitle("ElHamra") +
  theme(panel.background = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=12),axis.text.x=element_text(size=12),
        axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.x=element_blank()) 

#by month
ggplot(cladetab[cladetab$month=="aug",], aes(fill=clade, y=per, x=city)) + 
  geom_bar(position="fill", stat="identity") +
  theme(panel.background = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=12),axis.text.x=element_text(size=12,angle=90,hjust=1),
        axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.x=element_blank()) 

ggplot(cladetab[cladetab$month=="sept",], aes(fill=clade, y=per, x=city)) + 
  geom_bar(position="fill", stat="identity") +
  theme(panel.background = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=12),axis.text.x=element_text(size=12,angle=90,hjust=1),
        axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.x=element_blank()) 

ggplot(cladetab[cladetab$month=="oct",], aes(fill=clade, y=per, x=city)) + 
  geom_bar(position="fill", stat="identity") +
  theme(panel.background = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=12),axis.text.x=element_text(size=12,angle=90,hjust=1),
        axis.title.x=element_blank(),axis.title.y=element_blank(),axis.ticks.x=element_blank()) 


























mat = matrix(NA,ncol=length(sft)*3,nrow=length(freqPos)); rownames(mat)=freqPos; colnames(mat)=rep(c("REF","nuc","perc"),length(sft))
counter = 1
for (i in 1:length(sft)){
  mat[,counter:(counter+2)] = as.matrix(sft[[i]][match(freqPos,sft[[i]]$pos),c(6,3,5)])
  counter = counter+3
}
write.table(mat,file="tmp.txt",sep="\t")


## PLOT HEATMAP of general patterns
library(gplots); library(RColorBrewer); Npalette = colorRampPalette(c("white","yellow","red"))(n = 299)

tab = read.delim2("tmp.txt"); rownames(tab)=tab[,1]; tab=tab[,-1] #import freq table from excel
hm = matrix(as.numeric(as.matrix(tab)),ncol=ncol(tab),nrow=nrow(tab)); rownames(hm)=rownames(tab); colnames(hm)=colnames(tab)
hml = log2(hm);hml[hml==-Inf]=0  #transform to log2
hmmelt = reshape2::melt(hml)

pdf("figs/GeneralVar_sewer_heatmap.pdf")
ggplot(data = hmmelt, aes(x = Var1, y = Var2)) + 
  geom_tile(aes(fill = value), color = "white", size = 1) + 
  scale_fill_gradient(low = "cornsilk1", high = "tomato") + 
  theme_grey(base_size = 10) + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 8, colour = "gray50"),
        axis.text.x = element_text(angle = 90), text = element_text(size = 8))
dev.off()













### PLOT

### per month, all cities
names(puList)
augList = puList[c(1,6,10,15,19,23,27,32,35)]
septList = puList[c(2,7,11,16,20,24,28,33,36)]
octList = puList[c(3,8,12,17,21,25,29,34,37)]
novList = puList[c(4,9,13,22,26,30,38)]
decList = puList[c(5,14,18,31)]
mList = decList

colors = rainbow(length(mList))
plot(1, type="n", xlab="", ylab="", xlim=c(0,30000), ylim=c(0,100))
for (i in 1:length(mList)){ #for each city
  x = mList[[i]]
  for (j in 1:nrow(x)){ #for each position in the city
    if ((x[j,]$diff != 0) && (x[j,]$perc > 10) && (x[j,]$count >= 5)) { 
      points(x[j,]$pos,x[j,]$perc, col=colors[i], pch=20) }
  }
}
#legend(30,80,legend=c("BeerSheva","Ashdod","Rahat","Jerusalem-Sorek","Shafdan","Natanya","Haifa","Tzfat","ElHamra"),col=colors,pch=20)





### per city, all months
names(puList)
cList = puList[c(25:27)]#puList[c(22:24)]#puList[c(19:21)]#puList[c(16:18)]#puList[c(13:15)]#puList[c(10:12)] #puList[c(7:9)]#puList[c(4:6)]#puList[c(1:3)]

colors = c("darkorange2","cyan3","purple") #rainbow(length(cList)) #barplot(1:length(cList), col=c("lightgreen","cyan3","purple"))
plot(1, type="n", xlab="", ylab="", xlim=c(0,30000), ylim=c(0,100))
for (i in 1:length(cList)){ #for each month
  x = cList[[i]]
  for (j in 1:nrow(x)){ #for each position in the month
    if ((x[j,]$diff != 0) && (x[j,]$perc > 10) && (x[j,]$count >= 5)) { 
      points(x[j,]$pos,x[j,]$perc, col=colors[i], pch=20) }
  }
}
#legend(30,80,legend=c("Aug","Sept","Oct"),col=colors,pch=20)

text(x=29260, y=101, offset=0.7, labels="x",cex=.7,col="blue")








### per city, all months, just spike
names(puList_freqPos)
cList = puList_freqPos[c(1:5)] #beerSheva
cList = puList_freqPos[c(27:31)] #Haifa
cList = puList_freqPos[c(23:26)] #Natanya

par(mfrow=c(5,1))
for (i in 1:length(cList)){
  x=cList[[i]]
  plot(x$pos,x$perc,pch=20,col="darkgrey") #draw all
  points(x[x$perc>60,]$pos,x[x$perc>60,]$perc,col="red",pch=20) #color high freq
  thigmophobe.labels(x[x$perc>60,]$pos,x[x$perc>60,]$perc,labels=x[x$perc>60,]$pos,cex=.8,col="red")
  
  points(x[x$perc>20 & x$perc<60,]$pos,x[x$perc>20 & x$perc<60,]$perc,col="darkorange2",pch=20) #color mid freq
  thigmophobe.labels(x[x$perc>20 & x$perc<60,]$pos,x[x$perc>20 & x$perc<60,]$perc,labels=x[x$perc>20 & x$perc<60,]$pos,cex=.7,col="darkorange2")
  
  thigmophobe.labels(x[x$perc<20,]$pos,x[x$perc<20,]$perc,labels=x[x$perc<20,]$pos,cex=.6,col="darkgrey")
  
}



### per mutation, all cities / months
mut = 23604#21768#21765#8782#1059#28881# 23604
mutarr = rep(NA,length(puList_freqPos)); names(mutarr)=names(puList_freqPos)

for (i in 1:length(puList_freqPos)){
  x=puList_freqPos[[i]]
  if (mut %in% x$pos) {
    mutarr[i]=x[match(mut,x$pos),]$perc}
}
barplot(mutarr,ylim=c(0,100),las=2,cex.names=.7,col="lightgrey")
box(lty=1)





####  COMPARE TO HUMAN MUTATIONS   - BRITISH MUTS ####
load("/Users/netazuck/Documents/projects/Michal/CoronaVirus/data/NGS_runs/refs/NovelMuts_new.RData") 
UKmuts = c(21765,21766,21767,21768,21769,21770,23063,23271,23604,23709,24506,24914,27972,28048,28111,28280,28977,3267,5388,6954,11288,11289,11290,11291,11292,11293,11294,11295,11296)


## AUG
multAlign = read.fasta(file="Humans/all_aligned.fasta", seqtype="DNA")  #aligned w/ augur as part of nextstrain pipeline
multAlignNames = names(multAlign)
NucList = list()  #put sequences in list
for (i in 1:length(multAlign)){ NucList[[i]] = multAlign[[i]]}
names(NucList) = multAlignNames

NucMat = matrix(NA,nrow=length(NucList),ncol=length(NucList[[1]])); rownames(NucMat)=names(NucList); colnames(NucMat)=c(1:length(NucList[[1]]))
for (i in 1:length(NucList)){
  for (j in 1:length(NucList[[1]])){
    NucMat[i,j] = NucList[[i]][j]
  }
}

augList = read.delim2("Humans/aug/aug.txt")
table(augList$district)
augMat = NucMat[c(1,match(augList$strain,rownames(NucMat))),]
augMat_UK = augMat[,UKmuts]

mutfreq = rep(0,length(UKmuts))
for (i in 1:length(mutfreq)){ if (length(unique(augMat_UK[,i])) > 1) { mutfreq[i]=UKmuts[i] }}
which(mutfreq!=0)


## SEPT
multAlign = read.fasta(file="Humans/all_aligned.fasta", seqtype="DNA")  #aligned w/ augur as part of nextstrain pipeline
multAlignNames = names(multAlign)
NucList = list()  #put sequences in list
for (i in 1:length(multAlign)){ NucList[[i]] = multAlign[[i]]}
names(NucList) = multAlignNames

NucMat = matrix(NA,nrow=length(NucList),ncol=length(NucList[[1]])); rownames(NucMat)=names(NucList); colnames(NucMat)=c(1:length(NucList[[1]]))
for (i in 1:length(NucList)){
  for (j in 1:length(NucList[[1]])){
    NucMat[i,j] = NucList[[i]][j]
  }
}

septList = read.delim2("Humans/sept/sept.txt")
table(septList$district)
septMat = NucMat[c(1,match(septList$strain,rownames(NucMat))),]
septMat_UK = septMat[,UKmuts]

mutfreq = rep(0,length(UKmuts))
for (i in 1:length(mutfreq)){ if (length(unique(septMat_UK[,i])) > 1) { mutfreq[i]=UKmuts[i] }}
which(mutfreq!=0)


## OCT
augur align \
--sequences Humans/oct/allOctober_TLV.fasta \
--reference-sequence refs/REF_NC_045512.2.fasta \
--output Humans/oct/allOctober_TLV_aligned.fasta 

multAlign = read.fasta(file="Humans/oct/allOctober_TLV_aligned.fasta", seqtype="DNA")  #aligned w/ augur as part of nextstrain pipeline
multAlignNames = names(multAlign)
NucList = list()  #put sequences in list
for (i in 1:length(multAlign)){ NucList[[i]] = multAlign[[i]]}
names(NucList) = multAlignNames

NucMat = matrix(NA,nrow=length(NucList),ncol=length(NucList[[1]])); rownames(NucMat)=names(NucList); colnames(NucMat)=c(1:length(NucList[[1]]))
for (i in 1:length(NucList)){
  for (j in 1:length(NucList[[1]])){
    NucMat[i,j] = NucList[[i]][j]
  }
}

octList = read.delim2("Humans/oct/oct.txt")
tmp = match(octList$Barcode,rownames(NucMat))
octMat = NucMat[c(1,tmp),]  
octMat_UK = octMat[,UKmuts]

mutfreq = rep(0,length(UKmuts))
for (i in 1:length(mutfreq)){ if (length(unique(octMat_UK[,i])) > 1) { mutfreq[i]=UKmuts[i] }}
which(mutfreq!=0)

table(octList$district)
tmp = octMat_UK[which(octList$district=="Central"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = octMat_UK[which(octList$district=="Haifa"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = octMat_UK[which(octList$district=="Jerusalem"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = octMat_UK[which(octList$district=="Northern"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = octMat_UK[which(octList$district=="Southern"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = octMat_UK[which(octList$district=="Tel aviv"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = octMat_UK[which(octList$district=="West bank"),]; apply(tmp[,which(mutfreq!=0)],2,table)


## NOV
augur align \
--sequences Humans/nov/allNovemberConsortium.fasta \
--reference-sequence refs/REF_NC_045512.2.fasta \
--output Humans/nov/allNovemberConsortium_aligned.fasta 

multAlign = read.fasta(file="Humans/nov/allNovemberConsortium_aligned.fasta", seqtype="DNA")  #aligned w/ augur as part of nextstrain pipeline
multAlignNames = names(multAlign)
NucList = list()  #put sequences in list
for (i in 1:length(multAlign)){ NucList[[i]] = multAlign[[i]]}
names(NucList) = multAlignNames

NucMat = matrix(NA,nrow=length(NucList),ncol=length(NucList[[1]])); rownames(NucMat)=names(NucList); colnames(NucMat)=c(1:length(NucList[[1]]))
for (i in 1:length(NucList)){
  for (j in 1:length(NucList[[1]])){
    NucMat[i,j] = NucList[[i]][j]
  }
}

novList = read.delim2("Humans/nov/nov.txt")
tmp = rownames(NucMat[-1,]); for(i in 1:length(tmp)){ tmp[i]=strsplit(tmp[i],"_")[[1]][1]}; tmp=c("NC_045512.2",tmp); rownames(NucMat)=tmp
tmp = match(novList$Barcode,rownames(NucMat))
novMat = NucMat[c(1,tmp),]  
novMat_UK = novMat[,UKmuts]

mutfreq = rep(0,length(UKmuts))
for (i in 1:length(mutfreq)){ if (length(unique(novMat_UK[,i])) > 1) { mutfreq[i]=UKmuts[i] }}
which(mutfreq!=0)

table(novList$district)
tmp = novMat_UK[which(octList$district=="Central"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = novMat_UK[which(octList$district=="Haifa"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = novMat_UK[which(octList$district=="Jerusalem"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = novMat_UK[which(octList$district=="Northern"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = novMat_UK[which(octList$district=="Southern"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = novMat_UK[which(octList$district=="Tel aviv"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = novMat_UK[which(octList$district=="West bank"),]; apply(tmp[,which(mutfreq!=0)],2,table)


## DEC
multAlign = read.fasta(file="Humans/all_aligned.fasta", seqtype="DNA")  #aligned w/ augur as part of nextstrain pipeline
multAlignNames = names(multAlign)
NucList = list()  #put sequences in list
for (i in 1:length(multAlign)){ NucList[[i]] = multAlign[[i]]}
names(NucList) = multAlignNames

NucMat = matrix(NA,nrow=length(NucList),ncol=length(NucList[[1]])); rownames(NucMat)=names(NucList); colnames(NucMat)=c(1:length(NucList[[1]]))
for (i in 1:length(NucList)){
  for (j in 1:length(NucList[[1]])){
    NucMat[i,j] = NucList[[i]][j]
  }
}

decList = read.delim2("Humans/dec/dec.txt")
tmp = match(decList$strain,rownames(NucMat))
decMat = NucMat[c(1,tmp),]; decMat_UK = decMat[,UKmuts]

mutfreq = rep(0,length(UKmuts))
for (i in 1:length(mutfreq)){ if (length(unique(decMat_UK[,i])) > 1) { mutfreq[i]=UKmuts[i] }}
which(mutfreq!=0)

table(decList$District)
tmp = decMat_UK[which(decList$District=="Central"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = decMat_UK[which(decList$District=="Haifa"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = decMat_UK[which(decList$District=="Jerusalem"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = decMat_UK[which(decList$District=="Northern"),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = decMat_UK[c(1,which(decList$District=="Southern ")),]; apply(tmp[,which(mutfreq!=0)],2,table)
tmp = decMat_UK[which(decList$District=="Tel Aviv"),]; apply(tmp[,which(mutfreq!=0)],2,table)















#### CHECK SPIKE MUTATIONS ####
#Smuts = c(21770,21771,21772,21773,21774,21775,23066,22922)
Omuts = c(3267,5388,6954,11288:11296,27972,28048,28111,28280,28977)
Smuts = c(21765:21770,21991:21993,23063,23271,23604,23709,24506,24914)
muts_brit = c(21765:21770,21991:21993,23063,23271,23604,23709,24506,24914,27972,28048,28111,28280,28977)
tmp = rep(NA,length(puList_freqPos))
for (i in 1:length(puList_freqPos)){
  if(length(which(puList_freqPos[[i]]$pos%in%muts_brit))>0){ tmp[i]=i }
  #if(length(which(puList_freqPos[[i]]$pos%in%Omuts))>0){ tmp[i]=i }
}
names(puList_freqPos)[tmp]

x=puList_freqPos[[3]]; x[x$pos%in%muts_brit,] #beer sheva, oct (deletion in S69/70 - 27%) **
x=puList_freqPos[[8]]; x[x$pos%in%muts_brit,] #Ashdod, oct (deletion in S69/70 - 5%) **
x=puList_freqPos[[10]]; x[x$pos%in%muts_brit,] #Rahat, aug (deletion in S69/70 - 26%) **
x=puList_freqPos[[13]]; x[x$pos%in%muts_brit,] #Rahat, nov (D1118H/S - 17%) **
x=puList_freqPos[[14]]; x[x$pos%in%muts_brit,] #Rahat, dec (P681H/S - 9%, D1118H - 7.5%) **
x=puList_freqPos[[17]]; x[x$pos%in%muts_brit,] #JerSor, oct (Y144del/S, 12%) **
x=puList_freqPos[[22]]; x[x$pos%in%muts_brit,] #Shafdan, nov (deletion in S69/70 - 46%) **
x=puList_freqPos[[25]]; x[x$pos%in%muts_brit,] #Natanya, oct (R52I/orf8 - 6%) **
x=puList_freqPos[[26]]; x[x$pos%in%muts_brit,] #Natanya, nov (P681H/S - 98%) **
x=puList_freqPos[[30]]; x[x$pos%in%muts_brit,] #Haifa, nov (P681H/S - 98%) **

















### GET NUCLEOTIDE MIX IN CLADE POSITIONS (from BAM files - bam-readcount) NOT USED  ####
export PATH=$PATH:/usr/local/bin/bam-readcount
bam-readcount -f refs/REF_NC_045512.2.fasta BAM/47-S-4354_S45_L001_001.mapped.sorted.bam --site-list site_list -b 20 > Depth/S-4354.txt
bam-readcount -f refs/REF_NC_045512.2.fasta BAM/42-S-3909_S40_L001_001.mapped.sorted.bam --site-list site_list -b 20 > Depth/S-3909.txt












### GET NUCLEOTIDE MIX IN CLADE POSITIONS (from BAM files - deepSNV) NOT USED ####
library(deepSNV)  #http://bioconductor.org/packages/release/bioc/html/deepSNV.html
bamfile = "bam/57-S-2366_S55_L001_001.mapped.sorted.bam"
b2r = bam2R(bamfile,q=20,chr="NC_045512.2",1,29300)
b2r[11083,]







#### DRAW COVERAGE FOR EACH SAMPLE ####

## create a list for all samples, each with w/ depth for the whole genome 
dfiles = list.files("QC/depth", pattern="*.txt", full.names=TRUE)
DepthList_WG = list() #empty list for each sample, each which will contain gdepth for the whole genome

for (i in 1:length(dfiles)) {  #for each sample
  sdepth = read.delim2(dfiles[i],header=F)
  DepthList_WG[[i]] = sdepth$V3
}
names(DepthList_WG) = rownames(NucMat[-1,])
save(DepthList_WG,file="results/DepthList_WG.RData")


## plot mean depths for all patients (line plot)
i=27
plot(DepthList_WG[[i]],type="l",
     main=paste(nameTable[match(names(DepthList_WG)[[i]],nameTable$sampleName),]$region,nameTable[match(names(DepthList_WG)[[i]],nameTable$sampleName),]$month,sep="_"))










































####  TRANSLATE TO AMINO ACIDS  ####

load("data/Strains/refs/StartEndRegions.RData") #CDS regions, NC_045512.2

NucMat_UTR5 = NucMat[,StartEndRegions[which(StartEndRegions$segment=="UTR5"),]$start:StartEndRegions[which(StartEndRegions$segment=="UTR5"),]$end]
NucMat_UTR3 = NucMat[,StartEndRegions[which(StartEndRegions$segment=="UTR3"),]$start:StartEndRegions[which(StartEndRegions$segment=="UTR3"),]$end]


# ORF1ab  (https://www.ncbi.nlm.nih.gov/protein/1802476803)
NucMat_ORF1ab = NucMat[,StartEndRegions[which(StartEndRegions$segment=="ORF1ab"),]$start:StartEndRegions[which(StartEndRegions$segment=="ORF1ab"),]$end]
AAmat_ORF1ab = matrix(NA,ncol=floor(ncol(NucMat_ORF1ab)/3),nrow=nrow(NucMat_ORF1ab))  
rownames(AAmat_ORF1ab)=rownames(NucMat_ORF1ab)
for (i in 1:nrow(AAmat_ORF1ab)) {  #translate
  tmp = translate(NucMat_ORF1ab[i,])
  for (j in 1:ncol(AAmat_ORF1ab)){  #put in matrix
    AAmat_ORF1ab[i,j] = tmp[j]
  }
}
AAmat_ORF1ab[,1:10]
AAmat_ORF1ab_iCVL = AAmat_ORF1ab


# NSP1 (https://www.ncbi.nlm.nih.gov/protein/1802476805, MESLV...)
NucMat_NSP1 = NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP1"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP1"),]$end]
AAmat_NSP1 = matrix(NA,ncol=floor(ncol(NucMat_NSP1)/3),nrow=nrow(NucMat_NSP1)); rownames(AAmat_NSP1)=rownames(NucMat_NSP1)
for (i in 1:nrow(AAmat_NSP1)) {  #translate
  tmp = translate(NucMat_NSP1[i,])
  for (j in 1:ncol(AAmat_NSP1)){  #put in matrix
    AAmat_NSP1[i,j] = tmp[j]
  }
}
AAmat_NSP1[,1:10]
AAmat_NSP1_iCVL = AAmat_NSP1


# NSP2 (https://www.ncbi.nlm.nih.gov/protein/1802476806, AYTRY...)
NucMat_NSP2 = NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP2"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP2"),]$end]
AAmat_NSP2 = matrix(NA,ncol=floor(ncol(NucMat_NSP2)/3),nrow=nrow(NucMat_NSP2)); rownames(AAmat_NSP2)=rownames(NucMat_NSP2)
for (i in 1:nrow(AAmat_NSP2)) {  #translate
  tmp = translate(NucMat_NSP2[i,])
  for (j in 1:ncol(AAmat_NSP2)){  #put in matrix
    AAmat_NSP2[i,j] = tmp[j]
  }
}
AAmat_NSP2[,1:10]
AAmat_NSP2_iCVL = AAmat_NSP2


# NSP3  (https://www.ncbi.nlm.nih.gov/protein/1802476807, APTKV...)
NucMat_NSP3 = NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP3"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP3"),]$end]
AAmat_NSP3 = matrix(NA,ncol=floor(ncol(NucMat_NSP3)/3),nrow=nrow(NucMat_NSP3)); rownames(AAmat_NSP3)=rownames(NucMat_NSP3)
for (i in 1:nrow(AAmat_NSP3)) {  #translate
  tmp = translate(NucMat_NSP3[i,])
  for (j in 1:ncol(AAmat_NSP3)){  #put in matrix
    AAmat_NSP3[i,j] = tmp[j]
  }
}
AAmat_NSP3[,1:10]
AAmat_NSP3_iCVL = AAmat_NSP3


# NSP4 (https://www.ncbi.nlm.nih.gov/protein/1802476808, KIVNN..)
NucMat_NSP4 = NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP4"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP4"),]$end]
AAmat_NSP4 = matrix(NA,ncol=floor(ncol(NucMat_NSP4)/3),nrow=nrow(NucMat_NSP4)); rownames(AAmat_NSP4)=rownames(NucMat_NSP4)
for (i in 1:nrow(AAmat_NSP4)) {  #translate
  tmp = translate(NucMat_NSP4[i,])
  for (j in 1:ncol(AAmat_NSP4)){  #put in matrix
    AAmat_NSP4[i,j] = tmp[j]
  }
}
AAmat_NSP4[,1:10]
AAmat_NSP4_iCVL = AAmat_NSP4


# NSP5 (https://www.ncbi.nlm.nih.gov/protein/1802476809, SGFRK...)
NucMat_NSP5 = NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP5"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP5"),]$end]
AAmat_NSP5 = matrix(NA,ncol=floor(ncol(NucMat_NSP5)/3),nrow=nrow(NucMat_NSP5)); rownames(AAmat_NSP5)=rownames(NucMat_NSP5)
for (i in 1:nrow(AAmat_NSP5)) {  #translate
  tmp = translate(NucMat_NSP5[i,])
  for (j in 1:ncol(AAmat_NSP5)){  #put in matrix
    AAmat_NSP5[i,j] = tmp[j]
  }
}
AAmat_NSP5[,1:10]
AAmat_NSP5_iCVL = AAmat_NSP5


# NSP6 (https://www.ncbi.nlm.nih.gov/protein/1802476810, SAVKRT...)
NucMat_NSP6 = NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP6"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP6"),]$end]
AAmat_NSP6 = matrix(NA,ncol=floor(ncol(NucMat_NSP6)/3),nrow=nrow(NucMat_NSP6)); rownames(AAmat_NSP6)=rownames(NucMat_NSP6)
for (i in 1:nrow(AAmat_NSP6)) {  #translate
  tmp = translate(NucMat_NSP6[i,])
  for (j in 1:ncol(AAmat_NSP6)){  #put in matrix
    AAmat_NSP6[i,j] = tmp[j]
  }
}
AAmat_NSP6[,1:10]
AAmat_NSP6_iCVL = AAmat_NSP6


# NSP7 (https://www.ncbi.nlm.nih.gov/protein/1802476811, SKMSD...)
NucMat_NSP7 = NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP7"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP7"),]$end]
AAmat_NSP7 = matrix(NA,ncol=floor(ncol(NucMat_NSP7)/3),nrow=nrow(NucMat_NSP7)); rownames(AAmat_NSP7)=rownames(NucMat_NSP7)
for (i in 1:nrow(AAmat_NSP7)) {  #translate
  tmp = translate(NucMat_NSP7[i,])
  for (j in 1:ncol(AAmat_NSP7)){  #put in matrix
    AAmat_NSP7[i,j] = tmp[j]
  }
}
AAmat_NSP7[,1:10]
AAmat_NSP7_iCVL = AAmat_NSP7


# NSP8 (https://www.ncbi.nlm.nih.gov/protein/1802476812, AIASEF...)
NucMat_NSP8 = NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP8"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP8"),]$end]
AAmat_NSP8 = matrix(NA,ncol=floor(ncol(NucMat_NSP8)/3),nrow=nrow(NucMat_NSP8)); rownames(AAmat_NSP8)=rownames(NucMat_NSP8)
for (i in 1:nrow(AAmat_NSP8)) {  #translate
  tmp = translate(NucMat_NSP8[i,])
  for (j in 1:ncol(AAmat_NSP8)){  #put in matrix
    AAmat_NSP8[i,j] = tmp[j]
  }
}
AAmat_NSP8[,1:10]
AAmat_NSP8_iCVL = AAmat_NSP8


# NSP9 (https://www.ncbi.nlm.nih.gov/protein/1802476813, NNELSP...)
NucMat_NSP9= NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP9"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP9"),]$end]
AAmat_NSP9 = matrix(NA,ncol=floor(ncol(NucMat_NSP9)/3),nrow=nrow(NucMat_NSP9)); rownames(AAmat_NSP9)=rownames(NucMat_NSP9)
for (i in 1:nrow(AAmat_NSP9)) {  #translate
  tmp = translate(NucMat_NSP9[i,])
  for (j in 1:ncol(AAmat_NSP9)){  #put in matrix
    AAmat_NSP9[i,j] = tmp[j]
  }
}
AAmat_NSP9[,1:10]
AAmat_NSP9_iCVL = AAmat_NSP9


# NSP10 (https://www.ncbi.nlm.nih.gov/protein/1802476814, AGNAT...)
NucMat_NSP10= NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP10"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP10"),]$end]
AAmat_NSP10 = matrix(NA,ncol=floor(ncol(NucMat_NSP10)/3),nrow=nrow(NucMat_NSP10)); rownames(AAmat_NSP10)=rownames(NucMat_NSP10)
for (i in 1:nrow(AAmat_NSP10)) {  #translate
  tmp = translate(NucMat_NSP10[i,])
  for (j in 1:ncol(AAmat_NSP10)){  #put in matrix
    AAmat_NSP10[i,j] = tmp[j]
  }
}
AAmat_NSP10[,1:10]
AAmat_NSP10_iCVL = AAmat_NSP10


# NSP11  (https://www.ncbi.nlm.nih.gov/protein/1802476820, SADAQS...)
NucMat_NSP11= NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP11"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP11"),]$end]
AAmat_NSP11 = matrix(NA,ncol=floor(ncol(NucMat_NSP11)/3),nrow=nrow(NucMat_NSP11)); rownames(AAmat_NSP11)=rownames(NucMat_NSP11)
for (i in 1:nrow(AAmat_NSP11)) {  #translate
  tmp = translate(NucMat_NSP11[i,])
  for (j in 1:ncol(AAmat_NSP11)){  #put in matrix
    AAmat_NSP11[i,j] = tmp[j]
  }
}
AAmat_NSP11[,1:10]
AAmat_NSP11_iCVL = AAmat_NSP11



# NSP12  (https://www.ncbi.nlm.nih.gov/protein/1802476815)   ###TRANSLATION WEIRD! note stem&loop elements for frameshifting
NucMat_NSP12= NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP12"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP12"),]$end]
AAmat_NSP12 = matrix(NA,ncol=floor(ncol(NucMat_NSP12)/3),nrow=nrow(NucMat_NSP12)); rownames(AAmat_NSP12)=rownames(NucMat_NSP12)
for (i in 1:nrow(AAmat_NSP12)) {  #translate
  tmp = translate(NucMat_NSP12[i,])
  for (j in 1:ncol(AAmat_NSP12)){  #put in matrix
    AAmat_NSP12[i,j] = tmp[j]
  }
}
AAmat_NSP12[,1:10]
AAmat_NSP12_iCVL = AAmat_NSP12


# NSP13  (https://www.ncbi.nlm.nih.gov/protein/1802476816, AVGAC...)
NucMat_NSP13= NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP13"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP13"),]$end]
AAmat_NSP13 = matrix(NA,ncol=floor(ncol(NucMat_NSP13)/3),nrow=nrow(NucMat_NSP13)); rownames(AAmat_NSP13)=rownames(NucMat_NSP13)
for (i in 1:nrow(AAmat_NSP13)) {  #translate
  tmp = translate(NucMat_NSP13[i,])
  for (j in 1:ncol(AAmat_NSP13)){  #put in matrix
    AAmat_NSP13[i,j] = tmp[j]
  }
}
AAmat_NSP13[,1:10]
AAmat_NSP13_iCVL = AAmat_NSP13


# NSP14  (https://www.ncbi.nlm.nih.gov/protein/1802476817, AENVT...)
NucMat_NSP14= NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP14"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP14"),]$end]
AAmat_NSP14 = matrix(NA,ncol=floor(ncol(NucMat_NSP14)/3),nrow=nrow(NucMat_NSP14)); rownames(AAmat_NSP14)=rownames(NucMat_NSP14)
for (i in 1:nrow(AAmat_NSP14)) {  #translate
  tmp = translate(NucMat_NSP14[i,])
  for (j in 1:ncol(AAmat_NSP14)){  #put in matrix
    AAmat_NSP14[i,j] = tmp[j]
  }
}
AAmat_NSP14[,1:10]
AAmat_NSP14_iCVL = AAmat_NSP14


# NSP15  (https://www.ncbi.nlm.nih.gov/protein/1802476818, SLENVAF...)
NucMat_NSP15= NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP15"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP15"),]$end]
AAmat_NSP15 = matrix(NA,ncol=floor(ncol(NucMat_NSP15)/3),nrow=nrow(NucMat_NSP15)); rownames(AAmat_NSP15)=rownames(NucMat_NSP15)
for (i in 1:nrow(AAmat_NSP15)) {  #translate
  tmp = translate(NucMat_NSP15[i,])
  for (j in 1:ncol(AAmat_NSP15)){  #put in matrix
    AAmat_NSP15[i,j] = tmp[j]
  }
}
AAmat_NSP15[,1:10]
AAmat_NSP15_iCVL = AAmat_NSP15


# NSP16  (https://www.ncbi.nlm.nih.gov/protein/1802476819, SSQAW...)
NucMat_NSP16= NucMat[,StartEndRegions[which(StartEndRegions$segment=="NSP16"),]$start:StartEndRegions[which(StartEndRegions$segment=="NSP16"),]$end]
AAmat_NSP16 = matrix(NA,ncol=floor(ncol(NucMat_NSP16)/3),nrow=nrow(NucMat_NSP16)); rownames(AAmat_NSP16)=rownames(NucMat_NSP16)
for (i in 1:nrow(AAmat_NSP16)) {  #translate
  tmp = translate(NucMat_NSP16[i,])
  for (j in 1:ncol(AAmat_NSP16)){  #put in matrix
    AAmat_NSP16[i,j] = tmp[j]
  }
}
AAmat_NSP16[,1:10]
AAmat_NSP16_iCVL = AAmat_NSP16


# SPIKE  (https://www.ncbi.nlm.nih.gov/protein/YP_009724390.1, MFVFL...)
NucMat_SPIKE= NucMat[,StartEndRegions[which(StartEndRegions$segment=="SPIKE"),]$start:StartEndRegions[which(StartEndRegions$segment=="SPIKE"),]$end]
AAmat_SPIKE = matrix(NA,ncol=floor(ncol(NucMat_SPIKE)/3),nrow=nrow(NucMat_SPIKE)); rownames(AAmat_SPIKE)=rownames(NucMat_SPIKE)
for (i in 1:nrow(AAmat_SPIKE)) {  #translate
  tmp = translate(NucMat_SPIKE[i,])
  for (j in 1:ncol(AAmat_SPIKE)){  #put in matrix
    AAmat_SPIKE[i,j] = tmp[j]
  }
}
AAmat_SPIKE[,1:10]
AAmat_SPIKE_iCVL = AAmat_SPIKE


# ORF3a  (https://www.ncbi.nlm.nih.gov/protein/1796318599, MDLF...)
NucMat_ORF3a= NucMat[,StartEndRegions[which(StartEndRegions$segment=="ORF3a"),]$start:StartEndRegions[which(StartEndRegions$segment=="ORF3a"),]$end]
AAmat_ORF3a = matrix(NA,ncol=floor(ncol(NucMat_ORF3a)/3),nrow=nrow(NucMat_ORF3a)); rownames(AAmat_ORF3a)=rownames(NucMat_ORF3a)
for (i in 1:nrow(AAmat_ORF3a)) {  #translate
  tmp = translate(NucMat_ORF3a[i,])
  for (j in 1:ncol(AAmat_ORF3a)){  #put in matrix
    AAmat_ORF3a[i,j] = tmp[j]
  }
}
AAmat_ORF3a[,1:10]
AAmat_ORF3a_iCVL = AAmat_ORF3a


# ENVELOPE  (https://www.ncbi.nlm.nih.gov/protein/1796318600, MYSFV...)
NucMat_ENVELOPE= NucMat[,StartEndRegions[which(StartEndRegions$segment=="ENVELOPE"),]$start:StartEndRegions[which(StartEndRegions$segment=="ENVELOPE"),]$end]
AAmat_ENVELOPE = matrix(NA,ncol=floor(ncol(NucMat_ENVELOPE)/3),nrow=nrow(NucMat_ENVELOPE)); rownames(AAmat_ENVELOPE)=rownames(NucMat_ENVELOPE)
for (i in 1:nrow(AAmat_ENVELOPE)) {  #translate
  tmp = translate(NucMat_ENVELOPE[i,])
  for (j in 1:ncol(AAmat_ENVELOPE)){  #put in matrix
    AAmat_ENVELOPE[i,j] = tmp[j]
  }
}
AAmat_ENVELOPE[,1:10]
AAmat_ENVELOPE_iCVL = AAmat_ENVELOPE


# MEMBRANE  (https://www.ncbi.nlm.nih.gov/protein/1796318601, MADSNG...)
NucMat_MEMBRANE= NucMat[,StartEndRegions[which(StartEndRegions$segment=="MEMBRANE"),]$start:StartEndRegions[which(StartEndRegions$segment=="MEMBRANE"),]$end]
AAmat_MEMBRANE = matrix(NA,ncol=floor(ncol(NucMat_MEMBRANE)/3),nrow=nrow(NucMat_MEMBRANE)); rownames(AAmat_MEMBRANE)=rownames(NucMat_MEMBRANE)
for (i in 1:nrow(AAmat_MEMBRANE)) {  #translate
  tmp = translate(NucMat_MEMBRANE[i,])
  for (j in 1:ncol(AAmat_MEMBRANE)){  #put in matrix
    AAmat_MEMBRANE[i,j] = tmp[j]
  }
}
AAmat_MEMBRANE[,1:10]
AAmat_MEMBRANE_iCVL = AAmat_MEMBRANE


# ORF6  (https://www.ncbi.nlm.nih.gov/protein/1796318602, MFHLV...)
NucMat_ORF6= NucMat[,StartEndRegions[which(StartEndRegions$segment=="ORF6"),]$start:StartEndRegions[which(StartEndRegions$segment=="ORF6"),]$end]
AAmat_ORF6 = matrix(NA,ncol=floor(ncol(NucMat_ORF6)/3),nrow=nrow(NucMat_ORF6)); rownames(AAmat_ORF6)=rownames(NucMat_ORF6)
for (i in 1:nrow(AAmat_ORF6)) {  #translate
  tmp = translate(NucMat_ORF6[i,])
  for (j in 1:ncol(AAmat_ORF6)){  #put in matrix
    AAmat_ORF6[i,j] = tmp[j]
  }
}
AAmat_ORF6[,1:10]
AAmat_ORF6_iCVL = AAmat_ORF6


# ORF7a  (https://www.ncbi.nlm.nih.gov/protein/1796318603, MKII...)
NucMat_ORF7a= NucMat[,StartEndRegions[which(StartEndRegions$segment=="ORF7a"),]$start:StartEndRegions[which(StartEndRegions$segment=="ORF7a"),]$end]
AAmat_ORF7a = matrix(NA,ncol=floor(ncol(NucMat_ORF7a)/3),nrow=nrow(NucMat_ORF7a)); rownames(AAmat_ORF7a)=rownames(NucMat_ORF7a)
for (i in 1:nrow(AAmat_ORF7a)) {  #translate
  tmp = translate(NucMat_ORF7a[i,])
  for (j in 1:ncol(AAmat_ORF7a)){  #put in matrix
    AAmat_ORF7a[i,j] = tmp[j]
  }
}
AAmat_ORF7a[,1:10]
AAmat_ORF7a_iCVL = AAmat_ORF7a


# ORF7b  (https://www.ncbi.nlm.nih.gov/protein/1820616061, MIELS...)
NucMat_ORF7b= NucMat[,StartEndRegions[which(StartEndRegions$segment=="ORF7b"),]$start:StartEndRegions[which(StartEndRegions$segment=="ORF7b"),]$end]
AAmat_ORF7b = matrix(NA,ncol=floor(ncol(NucMat_ORF7b)/3),nrow=nrow(NucMat_ORF7b)); rownames(AAmat_ORF7b)=rownames(NucMat_ORF7b)
for (i in 1:nrow(AAmat_ORF7b)) {  #translate
  tmp = translate(NucMat_ORF7b[i,])
  for (j in 1:ncol(AAmat_ORF7b)){  #put in matrix
    AAmat_ORF7b[i,j] = tmp[j]
  }
}
AAmat_ORF7b[,1:10]
AAmat_ORF7b_iCVL = AAmat_ORF7b


# ORF8  (https://www.ncbi.nlm.nih.gov/protein/1796318604, MKFLV...)
NucMat_ORF8= NucMat[,StartEndRegions[which(StartEndRegions$segment=="ORF8"),]$start:StartEndRegions[which(StartEndRegions$segment=="ORF8"),]$end]
AAmat_ORF8 = matrix(NA,ncol=floor(ncol(NucMat_ORF8)/3),nrow=nrow(NucMat_ORF8)); rownames(AAmat_ORF8)=rownames(NucMat_ORF8)
for (i in 1:nrow(AAmat_ORF8)) {  #translate
  tmp = translate(NucMat_ORF8[i,])
  for (j in 1:ncol(AAmat_ORF8)){  #put in matrix
    AAmat_ORF8[i,j] = tmp[j]
  }
}
AAmat_ORF8[,1:10]
AAmat_ORF8_iCVL = AAmat_ORF8


# NUCAP  (https://www.ncbi.nlm.nih.gov/protein/1798174255, MSDNG...)
NucMat_NUCAP= NucMat[,StartEndRegions[which(StartEndRegions$segment=="NUCAP"),]$start:StartEndRegions[which(StartEndRegions$segment=="NUCAP"),]$end]
AAmat_NUCAP = matrix(NA,ncol=floor(ncol(NucMat_NUCAP)/3),nrow=nrow(NucMat_NUCAP)); rownames(AAmat_NUCAP)=rownames(NucMat_NUCAP)
for (i in 1:nrow(AAmat_NUCAP)) {  #translate
  tmp = translate(NucMat_NUCAP[i,])
  for (j in 1:ncol(AAmat_NUCAP)){  #put in matrix
    AAmat_NUCAP[i,j] = tmp[j]
  }
}
AAmat_NUCAP[,1:10]
AAmat_NUCAP_iCVL = AAmat_NUCAP


# ORF10  (https://www.ncbi.nlm.nih.gov/protein/1798174256, MGYNIV...)
NucMat_ORF10= NucMat[,StartEndRegions[which(StartEndRegions$segment=="ORF10"),]$start:StartEndRegions[which(StartEndRegions$segment=="ORF10"),]$end]
AAmat_ORF10 = matrix(NA,ncol=floor(ncol(NucMat_ORF10)/3),nrow=nrow(NucMat_ORF10)); rownames(AAmat_ORF10)=rownames(NucMat_ORF10)
for (i in 1:nrow(AAmat_ORF10)) {  #translate
  tmp = translate(NucMat_ORF10[i,])
  for (j in 1:ncol(AAmat_ORF10)){  #put in matrix
    AAmat_ORF10[i,j] = tmp[j]
  }
}
AAmat_ORF10[,1:10]
AAmat_ORF10_iCVL = AAmat_ORF10




## FIND MUTATIONS - NUC / AA / RSmut
load("~/Documents/projects/Bioinformatics/AAproperties/AAproperties.RData")
load("data/Strains/refs/StartEndRegions.RData") #CDS regions, NC_045512.2

AAofChoice = list(AAmat_NSP1,AAmat_NSP2,AAmat_NSP3,AAmat_NSP4,AAmat_NSP5,AAmat_NSP6,AAmat_NSP7,AAmat_NSP8,
                  AAmat_NSP9,AAmat_NSP10,AAmat_NSP11,AAmat_NSP12,AAmat_NSP13,AAmat_NSP14,AAmat_NSP15,AAmat_NSP16,
                  AAmat_SPIKE,AAmat_ORF3a,AAmat_ENVELOPE,AAmat_MEMBRANE,AAmat_ORF6,AAmat_ORF7a,AAmat_ORF7b,AAmat_ORF8,
                  AAmat_NUCAP,AAmat_ORF10)
names(AAofChoice) = c("NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NSP11","NSP12","NSP13","NSP14",
                      "NSP15","NSP16","SPIKE","ORF3a","ENVELOPE","MEMBRANE","ORF6","ORF7a","ORF7b","ORF8","NUCAP","ORF10")
AAofChoice_name = names(AAofChoice)


mutTable_all = list()
for (AAofChoice_ind in 1:length(AAofChoice)){  #goes over all AA matrices
  #debug: AAofChoice_ind=12; AAofChoice_name[AAofChoice_ind]; dim(AAofChoice[[AAofChoice_ind]])
  
  AAstart = StartEndRegions[match(AAofChoice_name[AAofChoice_ind],StartEndRegions$segment),]$start
  AAend = StartEndRegions[match(AAofChoice_name[AAofChoice_ind],StartEndRegions$segment),]$end
  mat = NucMat[,AAstart:AAend]
  mutTrack_noNs = rep(NA,ncol(mat))
  AAcounter = 1; currAA = rep(NA,nrow(mat))
  
  mutTable = as.data.frame(matrix(NA,nrow=ncol(mat),ncol=3+nrow(mat)+nrow(mat)+1))
  colnames(mutTable)=c("nuc_count",rownames(mat),"aa_count",rownames(mat),"R/S","AAgroup")
  
  for (i in 1:ncol(mat)){  #find mutations
    nucpos = mat[,i]
    
    if (AAofChoice_name[AAofChoice_ind]=="NSP12" && AAcounter>931) { break }  #for some reason the AAcounter goes 1 extra...
    if (AAofChoice_name[AAofChoice_ind]=="ORF7b" && AAcounter>164) { break }  #for some reason the AAcounter goes 1 extra...
    
    else {
      if (i != ncol(mat))
        currAA = AAofChoice[[AAofChoice_ind]][,AAcounter]
      if (i%%3 == 0){ AAcounter = AAcounter+1 } #update AA only every 3 "i" iterations
      
      if (length(unique(nucpos)) > 1) { #if there are different nucleotides
        
        if (length(unique(nucpos)) == 2){ #if there's ONE difference, check if it's not "n" 
          if (!((any(unique(nucpos)=="n"))|(any(unique(nucpos)=="-")))) { #if the mutation is not "n" or "-"
            mutTrack_noNs[i] = i    #record mutational position
            mutTable[i,2:125] = nucpos
            mutTable[i,]$nuc_count=AAstart+i-1; mutTable[i,]$aa_count = AAcounter; mutTable[i,127:250]=currAA
            
            if (length(unique(currAA)) == 2 && !any(unique(currAA)=="X")) { 
              mutTable[i,]$`R/S`= "R"   #update to R mutation
              tmp=names(table(currAA))  #update if in the same AA group
              if (AAproperties[match(tmp[which.max(table(currAA))],AAproperties$Abbv1),]$properties_group !=
                  AAproperties[match(tmp[which.min(table(currAA))],AAproperties$Abbv1),]$properties_group)
              { mutTable[i,]$AAgroup = paste(AAproperties[match(tmp[which.min(table(currAA))],AAproperties$Abbv1),]$properties," (",tmp[which.min(table(currAA))],") , ",
                                             AAproperties[match(tmp[which.max(table(currAA))],AAproperties$Abbv1),]$properties," (",tmp[which.max(table(currAA))],")",sep="") }
            } 
            else if (length(unique(currAA)) > 2) { 
              mutTable[i,]$`R/S`= "R"   #update to R mutation
              tmp=names(table(currAA))  #update if in the same AA group
              if ("X"%in%tmp){ currAA = currAA[-which(currAA=="X")]; tmp=names(table(currAA)) }
              if (AAproperties[match(tmp[which.max(table(currAA))],AAproperties$Abbv1),]$properties_group !=
                  AAproperties[match(tmp[which.min(table(currAA))],AAproperties$Abbv1),]$properties_group) {
                mutTable[i,]$AAgroup = paste(AAproperties[match(tmp[which.min(table(currAA))],AAproperties$Abbv1),]$properties," (",tmp[which.min(table(currAA))],") , ",
                                             AAproperties[match(tmp[which.max(table(currAA))],AAproperties$Abbv1),]$properties," (",tmp[which.max(table(currAA))],")",sep="") }
            } 
            else { mutTable[i,]$`R/S`= "S" }
            
          }
        }
        
        if (length(unique(nucpos)) > 2){ #if there's more than one difference
          mutTrack_noNs[i] = i    #also record mutational position
          mutTable[i,2:125] = nucpos
          mutTable[i,]$nuc_count=AAstart+i-1; mutTable[i,]$aa_count=AAcounter; mutTable[i,127:250]=currAA
          
          if (length(unique(currAA)) == 2 && !any(unique(currAA)=="X")) { 
            mutTable[i,]$`R/S`= "R"   #update to R mutation
            tmp=names(table(currAA))  #update if in the same AA group
            if (AAproperties[match(tmp[which.max(table(currAA))],AAproperties$Abbv1),]$properties_group !=
                AAproperties[match(tmp[which.min(table(currAA))],AAproperties$Abbv1),]$properties_group)
            { mutTable[i,]$AAgroup = paste(AAproperties[match(tmp[which.min(table(currAA))],AAproperties$Abbv1),]$properties," (",tmp[which.min(table(currAA))],") , ",
                                           AAproperties[match(tmp[which.max(table(currAA))],AAproperties$Abbv1),]$properties," (",tmp[which.max(table(currAA))],")",sep="") }
          } 
          else if (length(unique(currAA)) > 2) { 
            mutTable[i,]$`R/S`= "R"   #update to R mutation
            tmp=names(table(currAA))  #update if in the same AA group
            if ("X"%in%tmp){ currAA = currAA[-which(currAA=="X")]; tmp=names(table(currAA)) }
            if (AAproperties[match(tmp[which.max(table(currAA))],AAproperties$Abbv1),]$properties_group !=
                AAproperties[match(tmp[which.min(table(currAA))],AAproperties$Abbv1),]$properties_group) {
              mutTable[i,]$AAgroup = paste(AAproperties[match(tmp[which.min(table(currAA))],AAproperties$Abbv1),]$properties," (",tmp[which.min(table(currAA))],") , ",
                                           AAproperties[match(tmp[which.max(table(currAA))],AAproperties$Abbv1),]$properties," (",tmp[which.max(table(currAA))],")",sep="") }
          } 
          else { mutTable[i,]$`R/S`= "S" }
          
        }
      }
    } #end else for "break"
  } #end "i" loop
  mutTrack_noNs=as.numeric(na.omit(mutTrack_noNs))
  mutTable = mutTable[mutTrack_noNs,]
  mutTable_all[[AAofChoice_ind]] = mutTable
  
}  #end loop for all AA matrices

names(mutTable_all) = names(AAofChoice)
mutTable_all_df = do.call(rbind.data.frame, mutTable_all)
write.table(mutTable_all_df,file="data/Eritreas/results/mutTable.txt",sep="\t")













