#Filtering oTSS with p-value > 0.05
library(ggplot2)
library(dplyr)
library(forcats)

setwd("D:/Dropbox/RandauLab/Publications/Small_proteins_sulfolobus/Data_analysis/oTSS/oTSS_filtered")

orphanTSS = read.table("Orphan.tsv", header=T)

filtered_orphanTSS = subset(orphanTSS, Difference >= 10 & Score >= 600)

ggplot(filtered_orphanTSS, aes(x=log10(Difference))) + geom_histogram()
ggplot(filtered_orphanTSS, aes(x=Score)) + geom_histogram()

write.table(filtered_orphanTSS, file="oTSS_diff10_score600_pval005.tab", sep="\t", quote = F, row.names = F)

#Matching oTSS with six-frame translation from Sulfolobus

split_orphan = split(filtered_orphanTSS, filtered_orphanTSS$Strand) #Split oTSS in two according to strand information


#Evaluating fwd oTSS
fwd_proteome = read.table("Six-frame-proteome/Protein_database/Start_codon_forward_proteome.txt")
fwd_otss = as.data.frame(split_orphan[2])
fwd_otss <- fwd_otss %>% mutate(otss_range = X..oTSS + 15)
t_fwd = fwd_otss
coding_otss = data.frame()

i=0
for(i in 1:nrow(fwd_proteome)){
  a = fwd_proteome[i,]
  t_fwd =  t_fwd %>% mutate(TranslationStart = ifelse(a >= X..oTSS & a <= otss_range,a,NA))
  t_fwd_a = na.omit(t_fwd)
  coding_otss = rbind(coding_otss,t_fwd_a)
  
}

#Evaluating rev oTSS
rev_otss  = as.data.frame(split_orphan[1])
rev_otss <- rev_otss %>% mutate(otss_range = X..oTSS - 15)
rev_proteome = read.table("Six-frame-proteome/Protein_database/Start_codon_reverse_proteome.txt")
t_rev = rev_otss

i=0
for(i in 1:nrow(rev_proteome)){
  a = rev_proteome[i,]
  t_rev =  t_rev %>% mutate(TranslationStart = ifelse(a <= X..oTSS & a >= otss_range,a,NA))
  t_rev_a = na.omit(t_rev)
  coding_otss = rbind(coding_otss,t_rev_a)
}

#Calculating distance from Translation start site
coding_otss = coding_otss %>% mutate(DistFromTranslation = abs(TranslationStart-X..oTSS))
which(duplicated(coding_otss$X..oTSS)==TRUE)
ggplot(coding_otss, aes(x=DistFromTranslation)) + geom_histogram(binwidth = 1)
write.table(coding_otss$TranslationStart, file="Six-frame_otss_match.txt",quote=F,row.names = F, col.names = F)


a = data.frame(Chr="NC_007181.1",Source="oTSS_study",Type="gene",Start=coding_otss$X..oTSS, Stop=coding_otss$otss_range,
               Value=".", Strand = coding_otss$X..Strand, Value2=".", Gene=coding_otss$X..Identifier)

write.table(a, file="coding_OTSS_r1.bed",sep="\t",quote = F, row.names = F, col.names = F)



#Filtering initial dataset
nc_otss_r1 = rbind(fwd_otss,rev_otss)
nc_otss_r1 = merge(nc_otss_r1, coding_otss, by="otss_range", all=TRUE)
which(duplicated(nc_otss_r1$X..oTSS.x)==TRUE)
nc_otss_r1 = subset(nc_otss_r1, is.na(nc_otss_r1$TranslationStart))
nc_otss_r1 = nc_otss_r1 %>% select_if(~ !any(is.na(.))) #List of RNAs without coding potential from the first test
nc_otss_r1 = nc_otss_r1 %>% mutate(Stop = case_when(X..Strand.x == "-" ~ X..oTSS.x+0,
                                                     X..Strand.x == "+" ~ X..oTSS.x+250))

nc_otss_r1 = nc_otss_r1 %>% mutate(Start = case_when(X..Strand.x == "-" ~ X..oTSS.x-250,
                                                    X..Strand.x == "+" ~ X..oTSS.x+0))

a = data.frame(Chr="NC_007181.1",Source="oTSS_study",Type="gene",Start=nc_otss_r1$Start, Stop=nc_otss_r1$Stop,
               Value=".", Strand = nc_otss_r1$X..Strand.x, Value2=".",Gene=nc_otss_r1$X..Identifier.x)

write.table(a, file="ncOTSS_r1.bed",sep="\t",quote = F, row.names = F, col.names = F)

#Second round of filtereing from codingOTSS
filter_tab = read.table("CodingOTSS_curatedR1.tab", header=T)
filter_df = merge(coding_otss,filter_tab, by="X..Identifier")
write.table(filter_df, file="coding_OTSS_r2.bed", col.names = T, row.names = F, sep="\t", quote = F)

#a = data.frame(Chr="NC_007181.1",Source="oTSS_study",Type="gene",Start=filter_df$X..oTSS, Stop=filter_df$otss_range,
#               Value=".", Strand = filter_df$X..Strand, Value2=".", Gene=filter_df$X..Identifier)

#write.table(a, file="coding_OTSS_r2.bed",sep="\t",quote = F, row.names = F, col.names = F)

#Loading spTSS after curation
curated_sptss = read.table(file="spTSS_info_r2.tab", header=T, sep="\t")

#UTR length distribution
ggplot(curated_sptss, aes(x=DistFromTranslation),colour="black") + 
  geom_histogram(aes(y=..density..),binwidth = 1, color = "black", fill="light grey") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12,
                                  face="bold"))+
  ggtitle('UTR Length of novel pTSS')+
  ylab("Counts")+
  xlab("Distance from Start Codon")+
  geom_density(lwd = 1, colour = "black")

#Protein size distribution
ggplot(curated_sptss, aes(x=Protein_size),colour="black") + 
  geom_histogram(binwidth = 1, color = "black", fill="light grey") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12,
                                  face="bold"))+
  ggtitle(paste('Number of new putative proteins:',nrow(curated_sptss),'            Size range :',min(curated_sptss$Protein_size),
'~',max(curated_sptss$Protein_size),'aa'))+
  ylab("Counts")+
  xlab("Protein length (aa)")+
  geom_vline(xintercept = 50,
             linetype="dashed", lwd = 1.2)+
  geom_vline(xintercept = 25,
             linetype="dashed", lwd = 1.2)

#Start codong usage:
codon_usage = read.table("start_codon_spTSSr2.tab",sep="\t")
start_codon_usage = codon_usage[order(codon_usage$V2, decreasing = T),]


ggplot(start_codon_usage, aes(x=fct_inorder(V1), y=V2)) + 
  geom_segment( aes(x=fct_inorder(V1), 
                    xend=fct_inorder(V1), 
                    y=0, yend=V2), 
                    color="grey",
                    lwd=1.05)+
  geom_point(size=4) + 
  theme_bw() +
  theme(
        plot.title = element_text(size=12,
                                  face="bold"))+
  ggtitle('Codon usage of putative proteins')+
  ylab("Counts")+
  xlab("Start Codon")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 55))
 

#Loading pHMMER results
#First check: Trembl analysis
phmmer_trembl_res = read.table("pHMMER-results/sulfolobus_spTSS/spTSS_trembl_tblout.tab", sep="\t", header=T)
phmmer_trembl_doublesig = subset(phmmer_trembl_res, E.value <= 0.05 & E.value.1 <= 0.05)
phmmer_trembl_singlesig = subset(phmmer_trembl_res, E.value <= 0.05)

#Second check: KB Analysis
phmmer_kb_res = read.table("pHMMER-results/sulfolobus_spTSS/spTSS_kb_tblout.tab", sep="\t", header=T)
phmmer_kb_doublesig = subset(phmmer_kb_res, E.value <= 0.05 & E.value.1 <= 0.05)
phmmer_kb_singlesig = subset(phmmer_kb_res, E.value <= 0.05)
