{
  library(ggplot2)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(ggpubr)
  library(purrr)
  library(plotly)
  library(ggrepel)
  library(ggalluvial)
  library(RIdeogram)
  library(ggtree)
  library(treeio)
  library(tibble)
  library(stringr) 
  library(ggnewscale)
  library(ggtreeExtra)
  library(ggstar)
  library(ape)
  library(patchwork)
  library(phytools)
  library(ggsignif)
  library(UpSetR)
  library(changepoint)
  library(ComplexHeatmap)
  library(circlize)
  library(rpart)
  library(changepoint)
  library(ggforce)
  library(ggmsa)
} # Library
{
scaf_ratio <- read.table("vac.scaf_mf.ratio") 
colnames(scaf_ratio)<-c("scaffold","length","sum_ratio","mean_ratio","count_ratio")
summary(scaf_ratio)
#filter the too short scaffolds
main_scaf <- scaf_ratio %>% filter(length>10e4,count_ratio<2)
main_scaf <- main_scaf %>% mutate(chromosomes=ifelse(mean_ratio>1.8,"Chr Z",
                                      ifelse((length>10e5),"Autosomes",
                                             ifelse((mean_ratio<0.2)&(count_ratio<0.2),"Chr W",""))))
main_scaf_sub <- main_scaf %>% filter(chromosomes !="")
#confirm HiC_scaffold_18 is the Z chromosome
ggplot(main_scaf_sub, aes(x = count_ratio, y = mean_ratio)) + 
  geom_point(aes(size = length/1e6,color = chromosomes))+
  scale_color_manual(name="",values = c("Autosomes"="#4daf4a","Chr Z"="#0000cd","Chr W"="#cd0000"))+
  scale_size_continuous(name ="Size/Mb",range = c(2, 10))+
  labs(x="\u2642/\u2640 mappable site",y="\u2642/\u2640 coverage")+
  coord_cartesian(xlim = c(0,1.2),ylim=c(0,2.2))+
  theme_bw()+theme(legend.position = "left",
                                 legend.key.size = unit(0.5, 'cm'),
                                 legend.title = element_text(size=10),
                                 legend.text = element_text(size=8))
{
matchinfo <- read.table("vac.zw.minimap.txt")
matchinfo <- matchinfo %>% mutate(strata=ifelse(V8 > 6.76e6,"S0",
                                          ifelse((V8>1e6),"S1",
                                                 ifelse((V8>0.57e6),"S2","PAR"))))
ggplot()+
  geom_segment(data = matchinfo, aes(x = (V3+V4) / 2000000, y = 1, 
                                       xend = (V8 + V9) / 2000000, yend = 2,color=strata)) +
  theme_minimal()
} # Z and ragtag_W minimap result



} # Sex chromosomes identification of Va
{
##data depth
scafdepth  <- read.table(file = "salvator_depth_scaf.bed")
scafdepth <-  scafdepth[order(-scafdepth$V3),]
scafdepth$V1 <- gsub("JAIXND010000","scaf",scafdepth$V1)
scafsub <- filter(scafdepth,V3>1e6) #filter 1mb scaffolds
ggplot(scafsub)+geom_density(aes(x=V5))+ theme_classic() + labs(x="depth",y="")
##Vsa pep homo by miniprot of Vac SC homo scafs
homoinfo  <- read.table(file = "vac2vsa.pep.homo.txt")
vsagenenum <- read.table(file = "vsa.geneNumber")
vsascaflen <- read.table(file = "vsa.scaflen.txt")
colnames(homoinfo) <- c("vac","s1","e1","gene","vsa","s2","e2")
homoinfosubSC <- data.frame(table(subset(homoinfo,vac=="chrW"|vac=="chrZ")$vsa))
homoinfosubSC$total <- vsagenenum$V1[match(homoinfosubSC$Var1,vsagenenum$V2)]
homoinfosubSC$perc <- homoinfosubSC$Freq/homoinfosubSC$total *100
homoinfosubSC <- homoinfosubSC[order(-homoinfosubSC$Freq),]
candiscaf <- head(homoinfosubSC,5)
##depth in 10k win of scafs
all10kdepth  <- read.table(file = "salvator_depth.10k.bed")
all10kdepth$V1 <- gsub("JAIXND010000","scaf",all10kdepth$V1)
all10kdepthsub <- all10kdepth[all10kdepth$V1 %in% candiscaf$Var1,]
all10kdepthsub2 <- all10kdepth %>% filter(V1=="scaf674.1"|V1=='scaf856.1')
ggplot(rbind(all10kdepthsub,all10kdepthsub2))+geom_boxplot(aes(x=factor(V1),y=V5),outliers = F)+
  theme(axis.text.x=element_text(angle=45, size=9,hjust = 1))+theme_classic()
##chr2 depth
all10kdepthsubChr2 <- all10kdepth %>% filter(V1 %in% chr2homo)
ggplot(all10kdepthsubChr2)+geom_boxplot(aes(x="",y=V5),outliers = F)+
  theme(axis.text.x=element_text(angle=45, size=9,hjust = 1))+theme_classic()
##Check best hits of these scaffolds on VAC
candiscaf <- c("scaf674.1","scaf802.1","scaf729.1","scaf715.1","scaf856.1","scaf746.1")
homoinfo <- homoinfo %>% 
  mutate(vac = ifelse(vac=="chr2" & e1 < 103e6, "chr2p",
                      ifelse(vac=="chr2" & s1 > 103e6,"chr2q",vac)))
create_pie <- function(df, candiscaf_item) {
  subset_df <- df[df$vsa == candiscaf_item, ]
  frequency <- table(subset_df$vac)
  freq_df <- as.data.frame(frequency)
  names(freq_df) <- c("vac", "count")
  freq_df$vac <- factor(freq_df$vac)
  
 pie_chart <- ggplot(freq_df, aes(x = "", y = count, fill = vac)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    ggtitle(paste("Pie Chart for", candiscaf_item)) +
    theme_void() +
    theme(legend.title = element_blank())
  
  return(pie_chart)
}
pie_charts <- lapply(candiscaf, function(scaf) create_pie(homoinfo, scaf))
ggarrange(plotlist = pie_charts, ncol = 2, nrow = 3)
##Check scaf729.1
homoinfo729 <- homoinfo %>% filter(vsa=="scaf729.1")
ggplot(homoinfo729)+
  geom_segment(aes(x = (s2 + e2) / 2000000, y = 1, 
                   xend = (s2 + e2) / 2000000 , yend = 2, col = vac))+theme_classic()+
  scale_color_manual(values = c(chr2p="#004a08",chr2q="#4daf4a",
                                chrZ="#0000cd",chrW="#cd0000"))
freq_df <- as.data.frame(table(homoinfo729$vac))
ggplot(freq_df, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.title = element_blank())
##Pie chart of all chr2 homologs scafs
chr2homo <- c("scaf827.1","scaf762.1","scaf001.1","scaf795.1","scaf789.1","scaf791.1","scaf703.1",
                           "scaf674.1","scaf792.1","scaf794.1","scaf685.1","scaf697.1")
chr2homo2tip <- c("scaf792.1","scaf794.1","scaf685.1","scaf697.1")
homosubchr2 <- homoinfo %>% filter(vsa %in% chr2homo2tip)
freq_df <- as.data.frame(table(homosubchr2$vac))
ggplot(freq_df, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.title = element_blank())
##Chiseq test
data <- matrix(c(520, 6, 62, 11), nrow = 2, byrow = TRUE,
               dimnames = list(Group = c("A", "B"),
                               Status = c("chr2", "chrW")))
print(data)
chi_test <- chisq.test(data)
print(chi_test)
##Back check Vac chrW to chr2 paralogs
wparainfo <- read.table("vac.chrWpep.besthit.info")  #paralog of w genes
zparainfo <- read.table("vac.chrZpep.besthit")  #paralog of z genes
targetdf <- rbind(wparainfo %>% dplyr::select(1:4),zparainfo) %>%  filter(V4=="chr2")

targetdf <- zparainfo %>%  filter(V4=="chr2")
targetdf <- wparainfo %>% dplyr::select(1:4) %>%  filter(V4=="chr2")


transloci <- read.table("vac.transcript.loci.txt",header = F,fill=T) #vac transcript location
colnames(targetdf) <- c("Gene","Chr","bestPara","bestParaLoci")
targetdf$zS <- transloci$V2[match(targetdf$Gene,transloci$V4)] #z start
targetdf$zE <- transloci$V3[match(targetdf$Gene,transloci$V4)] #z end
targetdf$paraS <- transloci$V2[match(targetdf$bestPara,transloci$V4)] #paralog start
targetdf$paraE <- transloci$V2[match(targetdf$bestPara,transloci$V4)] #paralog end
targetdf$bestParaLoci <- ifelse(targetdf$bestParaLoci=="chr2" & targetdf$paraE < 125e6, "chr2p",
                                ifelse(targetdf$bestParaLoci=="chr2" & targetdf$paraS > 125e6, "chr2q",targetdf$bestParaLoci)) #assign 2p 2q W and other chromosomes
##Density of VAC SC-genes to VSA scafs depth
homoinfoBackSC <- homoinfo %>% filter(gene %in% targetdf$Gene)
ggplot(data=(scafdepth %>% filter(V1 %in% homoinfoBack$vsa))) + geom_density(aes(x=V5),size=2)+xlim(0,180)
depth1 <- scafdepth %>% filter(V1 %in% homoinfoBackSC$vsa)

##Density of VAC Chr2-genes to VSA scafs depth
homoinfoBackChr2 <- homoinfo %>% filter(gene %in% targetdf$bestPara)
depth2 <- scafdepth %>% filter(V1 %in% homoinfoBackChr2$vsa)
ggplot(data=(scafdepth %>% filter(V1 %in% homoinfoBackChr2$vsa))) + geom_density(aes(x=V5),size=2)+xlim(0,180) 


test <- rbind(homoinfoBackSC,homoinfoBackChr2)


test2 <- subset(wparainfo,V4=="chr2")
test2$wscaf <- test$vsa[match(test2$V1,test$gene)]
test2$autoscaf <- test$vsa[match(test2$V3,test$gene)]
test2$wscafS <- test$s2[match(test2$V1,test$gene)]
test2$wscafE <- test$e2[match(test2$V1,test$gene)]
test2$autoscafS <- test$s2[match(test2$V3,test$gene)]
test2$autoscafE <- test$e2[match(test2$V3,test$gene)]
test2$wDEP <- scafdepth$V5[match(test2$wscaf,scafdepth$V1)]
test2$aDEP <- scafdepth$V5[match(test2$autoscaf,scafdepth$V1)]
test2$wtype <- ifelse(test2$wDEP<100,"w-w","a-w")
test2$atype <- ifelse(test2$aDEP<100,"w-a","a-a")

test2 <- test2 %>%
  mutate(type = ifelse(wscafE >= autoscafS & wscafS <= autoscafE, "T", "F"))



homoinfoscro <-read.table(file = "C:/Users/13719/Desktop/r_file/vac2scro.pep.homo.txt")


test2$Wscro <- homoinfoscro$V5[match(test2$V1,homoinfoscro$V4)]
test2$Ascro <- homoinfoscro$V5[match(test2$V3,homoinfoscro$V4)]

test2$WscroS <- homoinfoscro$V6[match(test2$V1,homoinfoscro$V4)]
test2$WscroE <- homoinfoscro$V7[match(test2$V1,homoinfoscro$V4)]

test2$AscroS <- homoinfoscro$V6[match(test2$V3,homoinfoscro$V4)]
test2$AscroE <- homoinfoscro$V7[match(test2$V3,homoinfoscro$V4)]

test2 <- test2 %>%
  mutate(type2 = ifelse(WscroE >= AscroS & WscroS <= AscroE, "T", "F"))

##data repeat
vsarepeat  <- read.table(file = "vsa.final.repeatlen.100k.txt")
target <- "scaf654.1"
vsarepeatsub <- subset(vsarepeat,V1==target)
ggplot(vsarepeatsub)+geom_step(aes(x=V2/1e6,y=V4))


} # Sex chromosomes identification of Vs
{
homoinfo  <- read.table(file = "vac2gga.eachBestHit.addloci.txt")
homoinfo <- homoinfo %>% dplyr::select(1:3,5:7)
homoinfo$fill <- "cccccc"
colnames(homoinfo) <- c("Species_1","Start_1","End_1", "Species_2", "Start_2", "End_2", "fill")
synteny_dual_comparison <- homoinfo 
##
synteny_dual_comparison <- subset(synteny_dual_comparison,(synteny_dual_comparison$Species_1=="chr2"|
                                                             synteny_dual_comparison$Species_1=="chrZ"|
                                                             synteny_dual_comparison$Species_1=="chrW")&(
                                                               synteny_dual_comparison$Species_2=="chrZ"|
                                                               synteny_dual_comparison$Species_2=="chr1"|
                                                              synteny_dual_comparison$Species_2=="chr12"|
                                                                synteny_dual_comparison$Species_2=="chr13"|
                                                                synteny_dual_comparison$Species_2=="chr16"|
                                                                synteny_dual_comparison$Species_2=="chr18"|
                                                                synteny_dual_comparison$Species_2=="chr28"
                                                             ))
##assign chr ID
synteny_dual_comparison$Species_1 <- gsub("\\bchr2\\b","1",synteny_dual_comparison$Species_1)
synteny_dual_comparison$Species_1 <- gsub("\\bchrZ\\b","2",synteny_dual_comparison$Species_1)
synteny_dual_comparison$Species_1 <- gsub("\\bchrW\\b","3",synteny_dual_comparison$Species_1)
synteny_dual_comparison$Species_2 <- gsub("\\bchrZ\\b","1",synteny_dual_comparison$Species_2)
synteny_dual_comparison$Species_2 <- gsub("\\bchr1\\b","2",synteny_dual_comparison$Species_2)
synteny_dual_comparison$Species_2 <- gsub("\\bchr12\\b","3",synteny_dual_comparison$Species_2)
synteny_dual_comparison$Species_2 <- gsub("\\bchr13\\b","4",synteny_dual_comparison$Species_2)
synteny_dual_comparison$Species_2 <- gsub("\\bchr16\\b","5",synteny_dual_comparison$Species_2)
synteny_dual_comparison$Species_2 <- gsub("\\bchr18\\b","6",synteny_dual_comparison$Species_2)
synteny_dual_comparison$Species_2 <- gsub("\\bchr28\\b","7",synteny_dual_comparison$Species_2)
synteny_dual_comparison$Species_1 <- as.numeric(synteny_dual_comparison$Species_1)
synteny_dual_comparison$Species_2 <- as.numeric(synteny_dual_comparison$Species_2)
##assign cols as vac chrs
cols <- data.frame(t(data.frame("1"="FDA8A8", "2"="F8F8D8","3"="F8F8D8")))
rownames(cols) <- 1:3
cols$V2 <- rownames(cols)
colnames(cols) <- c("col","chr")
synteny_dual_comparison$fill <- cols$col[match(synteny_dual_comparison$Species_1,cols$chr)]
##genome info
Tchrsize <- data.frame(V1=c("1","2","3"),
                       V2=c("1","1","1"),
                       V3=c("284531044","11831237","15022549"))
Tchrsize <- Tchrsize[order(as.numeric(Tchrsize$V1)),]
Tchrsize$V1 <- as.numeric(Tchrsize$V1)
Tchrsize$fill <- "969696"
Tchrsize$species <- "Va"
Tchrsize$size <- "12"
Tchrsize$color <- "252525"
colnames(Tchrsize) <- c("Chr", "Start", "End", "fill", "species", "size", "color")
#
Qchrsize <- data.frame(V1=c("1","2","3","4","5","6","7"),
                       V2=c("1","1","1","1","1","1","1"),
                       V3=c("86044486","196449156","20119077",
                            "17905061","2706039","11623896","5437364"))
Qchrsize$fill <- "969696"
Qchrsize$species <- "Gg"
Qchrsize$size <- "12"
Qchrsize$color <- "252525"
colnames(Qchrsize) <- c("Chr", "Start", "End", "fill", "species", "size", "color")
karyotype_dual_comparison <- rbind(Tchrsize,Qchrsize)
table(karyotype_dual_comparison$species)
karyotype_dual_comparison$Start <- as.numeric(karyotype_dual_comparison$Start)
karyotype_dual_comparison$End <- as.numeric(karyotype_dual_comparison$End)
##visualize
ideogram(karyotype = karyotype_dual_comparison, synteny = synteny_dual_comparison)
convertSVG("chromosome.svg", device = "png")

} # Synteny of Va to Gga7b
{
vacrpt  <- read.table(file = "vac.rpt.txt")
vacchrlen <- read.table(file = "vac.chrLen.txt")
vacrpt <-  vacrpt %>% 
  filter(!str_detect(V1, "scaf"))  
vacsub <- vacrpt  %>%
  setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
  mutate(
    class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
    family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
    length = end - start,
    window = floor(start / 100000) * 100000
  ) %>%
  select(-factor) %>% 
  filter(!str_detect(family, "\\?")) %>% 
  filter(!str_detect(class, "\\?")) %>% 
  filter(!(class == "SINE" & family == "tRNA"))
{
massbyfamily <- vacsub %>%
  group_by(chr, family) %>%
  summarise(total_length = sum(length)) %>%
  ungroup() %>% mutate(
    chrlen = vacchrlen$V2[match(chr,vacchrlen$V1)]
  )  %>% filter(!str_detect(family, "\\?"))
massbyfamily_trans<- massbyfamily %>%
  mutate(value = total_length / chrlen) %>%        
  select(chr, family, value) %>%                    
  pivot_wider(names_from = chr, values_from = value) %>%
  column_to_rownames(var = "family")  
massbyfamily_trans[is.na(massbyfamily_trans)] <- 0
chr_order <- data.frame(
  chr = c(paste0("chr", 1:20), "chrZ", "chrW") %>%
    keep(~ !str_detect(.,"chr18")),  # 去掉 chr18
  type = rep(c("marcro", "micro"), times = c(8, 13)) 
)
massbyfamily_trans <- massbyfamily_trans[, chr_order$chr]
#Class info
classInfo <- unique(vacsub %>% select(family, class)) %>% 
  filter(!str_detect(family, "\\?")) %>% 
  filter(!str_detect(class, "\\?")) %>% 
  filter(!(class == "SINE" & family == "tRNA")) %>% 
  column_to_rownames(var = "family") %>% 
  mutate(
    class = ifelse(class %in% c("LINE", "SINE", "LTR", "DNA"), class, "Other")
  ) 
#Split and get ordered heat map for each Class
massDNA <- massbyfamily_trans[rownames(massbyfamily_trans) %in% rownames(subset(classInfo, class == "DNA")),]
massLINE <- massbyfamily_trans[rownames(massbyfamily_trans) %in% rownames(subset(classInfo, class == "LINE")),]
massSINE <- massbyfamily_trans[rownames(massbyfamily_trans) %in% rownames(subset(classInfo, class == "SINE")),]
massLTR <- massbyfamily_trans[rownames(massbyfamily_trans) %in% rownames(subset(classInfo, class == "LTR")),]
massOther <- massbyfamily_trans[rownames(massbyfamily_trans) %in% rownames(subset(classInfo, class == "Other")),]
#Each heatmap
create_heatmap <- function(data) {
  zscaled_data <- t(scale(t(data)))  
  classAnno <- classInfo[match(rownames(zscaled_data), rownames(classInfo)), ]
  annotationsRow <- rowAnnotation(
    Class = classAnno, 
    col = list(class = c("DNA" = "#0000cd", "LINE" = "#cd0000", 
                         "SINE" = "#4daf4a", "LTR" = "#984ea3", 
                         "Other" = "#ff7f00"))
  )
  col_fun <- colorRamp2(c(-4, 0, 4), c("#0000cd", "white", "#cd0000"))
  heatmap <- Heatmap(zscaled_data,
                     name = "Repeat abundance",
                     col = col_fun,
                     cluster_rows = TRUE,  
                     cluster_columns = FALSE,  
                     show_row_names = TRUE, 
                     show_column_names = TRUE, 
                     right_annotation = annotationsRow,
                     column_title = " ",
                     column_names_rot = 45,
                     column_split = factor(chr_order$type,
                                           levels = c("marcro", "micro")),
                     border = TRUE)
  return(heatmap)
}
htDNA <- create_heatmap(massDNA)
htLINE <- create_heatmap(massLINE)
htSINE <- create_heatmap(massSINE)
htLTR <- create_heatmap(massLTR)
htOther <- create_heatmap(massOther)
ht_list <- htDNA %v% htLINE %v% htSINE %v% htLTR %v% htOther
draw(ht_list)
} # Heat map of all TE normalized by chromosomes by families split strata
{
vacnewsub <- vacsub %>%
  mutate(chrtype = case_when(
    chr %in% paste0("chr", 1:8) ~ "macro",
    chr %in% paste0("chr", 9:20) ~ "micro",
    chr == "chrZ" ~ "chrZ",
    chr == "chrW" ~ "chrW",
    chr == "chr2" & start >= 237e6 ~ "chr2tip",
    chr == "chr2" & start >= 237e6 ~ "chr2tip",
    TRUE ~ NA_character_  
  )) %>% na.omit() #split S0 S1
newChrlen <- data.frame(chrtype = c("macro","micro","chrZS1","chrZS0","chrW"),
                        length = c("1280481873","267739032","5760000","4351237","15022549"))
vacnewsub <- vacsub %>%
  mutate(chrtype = case_when(
    chr == "chr2" & start >= 237e6 ~ "chr2tip",
    chr %in% paste0("chr", 1:8) ~ "macro",
    chr %in% paste0("chr", 9:20) ~ "micro",
    chr == "chrZ" ~ "chrZ",
    chr == "chrW" ~ "chrW",
    TRUE ~ NA_character_  
  )) %>% na.omit() 
newChrlen <- data.frame(chrtype = c("macro","micro", "chr2tip", "chrZ","chrW"),
                        length = c("1232950829","267739032", "47531044", "11831237","15022549"))
massbyfamily <- vacnewsub %>%
  group_by(chrtype, family) %>%
  summarise(total_length = sum(length)) %>%
  ungroup() %>% mutate(
    chrlen = newChrlen$length[match(chrtype,newChrlen$chrtype)]
  )  %>% filter(!str_detect(family, "\\?")) 
#top10 mass repeat
massbyfamilyTOP20 <- massbyfamily %>%
  group_by(family) %>%
  summarise(total_length_all = sum(total_length)) %>%
  ungroup() %>% filter(family!="Unknown") %>%
  arrange(desc(total_length_all)) %>%
  slice_head(n = 20)
#transform
massbyfamily_trans <- massbyfamily %>% filter(family %in% massbyfamilyTOP20$family) %>%
mutate(value = total_length / as.numeric(chrlen)) %>%        
  select(chrtype, family, value) %>%                    
  pivot_wider(names_from = chrtype, values_from = value) %>%
  column_to_rownames(var = "family")  
massbyfamily_trans[is.na(massbyfamily_trans)] <- 0
#Class info
classInfo <- unique(vacsub %>% dplyr::select(family, class)) %>% 
  filter(!str_detect(family, "\\?")) %>% 
  filter(!str_detect(class, "\\?")) %>% 
  filter(!(class == "SINE" & family == "tRNA")) %>% 
  column_to_rownames(var = "family") %>% 
  mutate(
    class = ifelse(class %in% c("LINE", "SINE", "LTR", "DNA"), class, "Other")
  ) 
#write.table(classInfo, file = "vac.repeat.family2class.txt")
#Each heatmap
create_heatmap <- function(data) {
  zscaled_data <- t(scale(t(data)))  
  classAnno <- classInfo[match(rownames(zscaled_data), rownames(classInfo)), ]
  annotationsRow <- rowAnnotation(
    Class = classAnno, 
    col = list(class = c("DNA" = "#0000cd", "LINE" = "#cd0000", 
                         "SINE" = "#4daf4a", "LTR" = "#984ea3", 
                         "Other" = "#ff7f00"))
  )
  col_fun <- colorRamp2(c(-2, 0, 2), c("#0000cd", "white", "#cd0000"))
  heatmap <- Heatmap(zscaled_data,
                     name = "Repeat abundance",
                     col = col_fun,
                     cluster_rows = TRUE,  
                     cluster_columns = FALSE,  
                     show_row_names = TRUE, 
                     show_column_names = TRUE, 
                     right_annotation = annotationsRow,
                     column_title = " ",
                     column_names_rot = 45,
                     row_names_rot = -45,
                     border = TRUE)
  return(heatmap)
}
create_heatmap(massbyfamily_trans)
#Check pvalues
targetZte <- c("hAT-Charlie", "TcMar-Tc2", "L2", "CR1", "Simple_repeat", "tRNA-Core-RTE",
               "5S", "MIR", "5S-Sauria-RTE", "ID", "Maverick")
massbyWindow <- vacnewsub %>% 
  filter(family %in% targetZte) %>% 
  filter(chrtype %in% c("chrZS0", "chrZS1")) %>% 
  group_by(chrtype, chr, family, window) %>%
  summarise(total_length = sum(length)) %>%
  ungroup() %>% filter(!str_detect(family, "\\?")) 
comparisons <- c("chrZS0", "chrZS1")
ggplot(data = massbyWindow) +
  geom_boxplot(aes(x = chrtype, y = log2(total_length + 1))) +
  stat_compare_means(
    comparisons = comparisons, 
    method = "wilcox.test"    
  ) + facet_wrap( ~ family, scales = "free_y")
ggplot(massbyWindow, aes(x = chrtype, y = total_length)) +
  geom_boxplot() +
  theme_classic() +
  stat_compare_means(comparisons = list(
    c("chrZS0", "chrZS1")), 
  method = "wilcox.test", method.args = list(alternative = "greater"),
  label = "p.signif" ,) + facet_wrap( ~ family, scales = "free_y")
} # Simplified TE heatmap with strata boxplot

{
massbyWindow <- vacsub %>%
  group_by(chr, window) %>%
  summarise(total_length = sum(length), .groups = "drop")  %>%
  group_by(chr) %>%
  mutate(norm_length = scale(total_length)) %>%  # Zscale
  ungroup() %>%  
  mutate(type = "All")
L1massbyWindow <- vacsub %>%
  filter(family == "L1") %>%
  group_by(chr, window) %>%
  summarise(total_length = sum(length), .groups = "drop") %>%
  group_by(chr) %>%
  mutate(norm_length = scale(total_length)) %>%  
  ungroup() %>%
  mutate(type = "L1")
massMerge <- rbind(massbyWindow, L1massbyWindow) 
massMerge$chr <- factor(massMerge$chr, levels = chr_order$chr)
ggplot(massMerge) + 
  geom_smooth(aes(x = window/1E6, y = norm_length, color = type), 
              method = "loess", span = 0.2) + 
  facet_wrap(~ chr, scales = "free") +
  theme_classic()+
  labs(title = "Z-scaled repeat mass", x = "", y ="") +scale_color_manual(values = c("All"="black","L1"="#cd0000")) 
} # Line plot of all TE and L1 density by chromosomes

} # Total Repeat abundance by chromosomes in Va
{
#M & F coverage lineplot
read.table("female_Z_10depth.bed") %>% mutate_at(4:6,rev) -> f_Z_depth
read.table("male_Z_10depth.bed") %>% mutate_at(4:6,rev) -> m_Z_depth
rbind(data.frame(f_Z_depth,gender="female(ZW)"),data.frame(m_Z_depth,gender="male(ZZ)"))->fm_Z_depth
colnames(fm_Z_depth)<-c("scaffold","start","end","sum","mean","count","gender")
fm_Z_depth%>%mutate(gender=factor(gender))->fm_Z_depth
ggplot(fm_Z_depth%>%filter(mean>5,mean<25),aes(x=start/1e6,y=mean,col=gender))+geom_point(size=0.2)+
    geom_smooth(method = "loess",se = F,span=0.3)+
  labs(x="Position(Mbp)",y="Coverage")+
    #xlim(c(0,1))+
  scale_color_manual(values = c("male(ZZ)"="#2565B8","female(ZW)"="#C43F3D"),name="Gender")+theme_bw()
#ratio
fm_Z_depth%>%select(-4,-6)%>%tidyr::spread(key =gender,mean)%>%mutate(ratio=`male(ZZ)`/`female(ZW)`)%>%
  filter(ratio<4)-> fm_Z_depth1
fm_Z_depth1=fm_Z_depth1%>%mutate(ratio1=ifelse(start/1e6<0.57,"1.0","2.0"))

p1 <- ggplot(fm_Z_depth1,aes(x=start/1e6,y=ratio))+
    geom_point(size=0.2)+
    geom_smooth(aes(col=ratio1),method = "loess",se = F,span=1)+
  labs(x="Position(Mbp)",y="\u2642/\u2640 coverage")+
  scale_color_manual(values = c("1.0"="#4DD6CD","2.0"="#B16BE0"),name="Ratio")+theme_bw()
p1 + geom_vline(xintercept = 0.57,linetype=2)+theme(legend.position = "none")
#ZW similarity lastz
read.table("../data/Z/zw.simi.lastz.final10k",header = F)%>%mutate(V2=11830000-V2)->zw.simi
colnames(zw.simi)<-c("Chr","Pos","Mis","All","Similarity")

(zw.simi%>%
    filter(All>200)%>%
    ggplot(.,aes(x=Pos/1e6,y=Similarity))+geom_point(aes(size=All),color="#17B0AB")+
    geom_smooth(method = "loess",se = F,span=0.35,color="#17B0AB")+
    scale_radius(range=c(0.1,0.4))+
    labs(x="Position(Mbp)",y="Z-W Similarity")+theme_bw()+theme(legend.position = "none")->simi_p)
#change point analysis
cpa(zw.simi,2,5)->simi_cp
simi_p+geom_vline(xintercept = simi_cp$pos[1:20]/1e6,linetype="dashed")
#Heterozygosity SNP
read.table("../data/Z/single_filtermale18_10k.snpden",header = T)%>%filter(CHROM=="HiC_scaffold_18")%>%
  mutate(BIN_START=11830000-BIN_START)->m_snpden
read.table("../data/Z/single_filterfemale18_10k.snpden",header = T)%>%filter(CHROM=="HiC_scaffold_18")%>%
  mutate(BIN_START=11830000-BIN_START)->f_snpden
rbind(data.frame(f_snpden,gender="female(ZW)"),data.frame(m_snpden,gender="male(ZZ)"))->fm_Z_snpden
fm_Z_snpden%>%mutate(gender=factor(gender))->fm_Z_snpden
(snp_p=ggplot(fm_Z_snpden,aes(x=BIN_START/1e6,y=VARIANTS.KB,col=gender))+geom_point(size=0.2)+
    geom_smooth(method = "loess",se = F,span=0.08)+
    labs(x="Position(Mbp)",y="SNP density on ChrZ")+
    #xlim(c(0,1))+
    scale_color_manual(values = c("male(ZZ)"="#2565B8","female(ZW)"="#C43F3D"),name="Gender")+theme_bw())
snp_p=ggplot()+geom_point(data = fm_Z_snpden%>%filter(gender=="female(ZW)"),aes(x=BIN_START/1e6,y=VARIANTS.KB),size=0.2,col="#C43F3D")+
  geom_smooth(data = fm_Z_snpden%>%filter(gender=="female(ZW)"),aes(x=BIN_START/1e6,y=VARIANTS.KB),method = "loess",se = F,span=0.06,col="#C43F3D")+
  geom_point(data = fm_Z_snpden%>%filter(gender=="male(ZZ)"),aes(x=BIN_START/1e6,y=VARIANTS.KB),size=0.2,col="#2565B8")+
  geom_smooth(data = fm_Z_snpden%>%filter(gender=="male(ZZ)"),aes(x=BIN_START/1e6,y=VARIANTS.KB),method = "loess",se = F,span=0.4,col="#2565B8")+
  labs(x="Position(Mbp)",y="SNP density on ChrZ")+theme_bw()
#change point analysis
cpa(m_snpden,2,4)->snp_cp
snp_p+geom_vline(xintercept = snp_cp$pos[1:30]/1e6,linetype="dashed")
snp_p+geom_vline(xintercept = c(0.3,0.55,0.86,1,6.76),linetype="dashed")+
  geom_text(aes(x = c(0.3,0.55,0.86,1.2,6.76),y=c(23,25,23,25,24),
                label= c(0.3,0.55,0.86,"1.0",6.76)),size=3)
} # Strata identification
{
##data of input
wparainfo <- read.table("vac.chrWpep.besthit.info")  #paralog of w genes
wgenestrata <- read.table("vac.w.trans.strata.txt") #strata of w genes defined by minimap ragtagW
geneloci <- read.table("vac_zjuV2.final.geneloci",header = F,fill=T) #vac gene location
transloci <- read.table("vac.transcript.loci.txt",header = F,fill=T) #vac transcript location
worf <- read.table("vac.w.orfinfo",header = F,fill=T) #valid orf of wgenes
vactpm <- read.table("vac.tpmmeanFinalfilt.txt",sep = ",",header = T,row.names = 1)
vactpm <- vactpm %>% dplyr::filter_all(.,any_vars(.>=1)) #remove unexpressed genes
##assign strata to w-genes that have known zw gametologs
translocisub <-  subset(transloci, transloci$V1=="chrZ")
translocisub$strata <- ifelse(translocisub$V2 >= 0 & translocisub$V3 <= 5.7e5, "PAR",
                          ifelse(translocisub$V2 > 5.7e5 & translocisub$V3 <= 1e6, "S2",
                                 ifelse(translocisub$V2 > 1e6 & translocisub$V3 <= 6.76e6, "S1", "S0")))
wparainfoSub <- wparainfo %>%
  mutate(strata = translocisub$strata[match(wparainfo$V3,translocisub$V4)]) %>% filter(V4=="chrZ") %>% dplyr::select(3,6)
colnames(wparainfoSub) <- c("wgene","strata")
wgenestrata <- wgenestrata %>% dplyr::select(4,5)
colnames(wgenestrata) <- c("wgene","strata")
wgenestrata <- rbind(wgenestrata,wparainfoSub)
##Count of w genes with intact orf gametolog
disruptORF <- read.table(file = "vac.disruptORF.txt")
df <- wparainfo[grepl("Z", wparainfo$V5),] %>% mutate(
  orf = ifelse(V1 %in% disruptORF$V1, "Dis","Int"),
  tpm = if_else(gsub("-T[0-9]", "", V1) %in% rownames(vactpm), "T", "F")
)
##Assign status
wgenestrata <- wgenestrata %>%
  mutate(
    bestPara = wparainfo$V3[match(wgene, wparainfo$V1)],
    bestParaLoci = wparainfo$V4[match(wgene, wparainfo$V1)],
    orf = worf$V2[match(wgene, worf$V1)]
  ) %>%
  filter(orf != "") %>% filter(is.na(bestParaLoci) | bestParaLoci == "chrZ") %>%
  mutate(
    orf = if_else(orf == "True", "T", "F"),
    tpm = if_else(gsub("-T[0-9]", "", wgene) %in% rownames(vactpm), "T", "F"),
    game = if_else(is.na(bestParaLoci),"F", "T")
  )
##
wgenestrata <- wgenestrata %>% dplyr::select(1,2,5:7)
wgenestrata <- subset(wgenestrata,wgenestrata$strata=="S0"|wgenestrata$strata=="S1")
##
merged_df <- data.frame(table(factor(paste(wgenestrata$strata,wgenestrata$orf,wgenestrata$tpm,wgenestrata$game))))
merged_df <- merged_df %>%
  separate(Var1, into = c("strata", "orf", "tpm", "game"), sep = " ")
##reorder
strata_levels <- c("S1", "S0")
sorted_df <- merged_df %>%
  mutate(
    strata = factor(strata, levels = strata_levels)
  ) %>%
  arrange(
    strata,
    desc(game == "T"),
    desc(orf == "T"),
    desc(tpm == "T")
  )
df1mut <- sorted_df %>%
  mutate(Var1 = paste(strata, orf, tpm, game)) %>%
  arrange(factor(Var1, levels = unique(Var1))) %>% 
  mutate(ypos = cumsum(Freq) - 0.5 * Freq) 
#visualize
ggplot(df1mut, aes(x = "", y = Freq, fill = factor(Var1,levels =df1mut$Var1))) +
  geom_bar(stat = "identity", width = 1,color="white") +
  coord_polar("y", start = 0)+
  theme_void()  +
  theme(legend.position = "left") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())
} # W gene degeneration pie
{
dnds  <- read.table(file = "vac2php.dnds.txt")
colnames(dnds) <- c("Vgene","w","dN","dS")
loci <- read.table("vac2pp.rbh.addloci.txt",header = F,fill=T)
dnds <- dnds %>%
  mutate(
    Vchr = loci$V1[match(Vgene, loci$V4)],
    Vstart = loci$V2[match(Vgene, loci$V4)],
    Vend = loci$V3[match(Vgene, loci$V4)],
    Pchr = loci$V5[match(Vgene, loci$V4)],
    Pstart = loci$V6[match(Vgene, loci$V4)],
    Pend = loci$V7[match(Vgene, loci$V4)],
  ) %>%
  filter(!str_starts(Vchr, "scaf")) %>%
  filter(dS > 0.01 & dS < 3) %>% na.omit() %>%
  mutate(
    Pchr = case_when(
      Pchr %in% paste0("chr", 1:6) ~ "macro",
      Pchr %in% c("chr3_1","chr3_2") ~ "macro",
      Pchr == "chr15" ~ "X",
      Pchr %in% paste0("chr", 7:17) ~ "micro"
    )
  )  %>% 
  mutate(
    Vchr = case_when(
      Vchr %in% paste0("chr", 1:8) ~ "macro",
      Vchr == "chrZ" ~ "Z",
      Vchr == "chrW" ~ "W",
      Vchr %in% paste0("chr", 9:20) ~ "micro"
    )
  )  %>%
  mutate(
    type = case_when(
      Vchr == "Z" & Vstart > 6.76e6 & Pchr == "macro" ~ "macro-S0",
      Vchr == "Z" & Vstart > 1e6 & Vend < 6.76e6 & Pchr == "macro" ~ "macro-S1",
      Vchr == "Z" & Vstart > 6.76e6 & Pchr == "micro" ~ "micro-S0",
      Vchr == "Z" & Vstart > 1e6 & Vend < 6.76e6 & Pchr == "micro" ~ "micro-S1",
      Vchr == "macro" & Pchr == "macro" ~ "macro-macro",
      Vchr == "micro" & Pchr == "micro" ~ "micro-micro",
      Vchr == "macro" & Pchr == "micro" ~ "macro-micro",
      Vchr == "micro" & Pchr == "macro" ~ "micro-macro",
      Vchr == "macro" & Pchr == "X" ~ "X-macro",
      Vchr == "micro" & Pchr == "X" ~ "X-micro",
    )
  ) %>% na.omit() %>%
  filter(type %in% c("macro-macro","micro-micro","micro-S0","micro-S1","X-micro"))  
#psig
type_levels <- unique(dnds$type)
comparisons <- combn(type_levels, 2, simplify = FALSE)
#
ggplot(dnds, aes(x=type, y=w,color=type))+
  stat_boxplot(aes(type, w), 
                geom='errorbar', linetype=1, width=0.5)+ 
  geom_boxplot( aes(type, w),outlier.shape=NA) + 
  theme_classic() +
  labs(x="", y="dN/dS") +
  stat_compare_means(
    method = "wilcox.test", 
    comparisons = comparisons, 
    label = "p.signif",
  ) 
} # Faster-Z pairwise dNdS with Va-Ph
{
{
allcount <- read.table(file="gga.all.gene.count",header = T,row.names = 1)
txdb <- makeTxDbFromGFF("GRCg7b.gtf",format="gtf")
txname <- c(row.names(allcount))
exons.list.per.gene <- exonsBy(txdb,by="gene")[txname]
exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))
##TPM
len = exonic.gene.sizes
tpm3 <- function(dat,len) {
    x <- dat/len
    return(t(t(x)*1e6/colSums(x)))
}
tpm = matrix(0,nrow = 25466,ncol = 12)
colnames(tpm) = colnames(allcount)
rownames(tpm) = rownames(allcount)
for (i in 1:12) {
    tpm[,i] = tpm3(allcount[,i],len)
  }
tpm <- data.frame(tpm)
#get mean of 2 reps
for (i in seq(from=1, to=12, by=2)) {
    new <- rowMeans(tpm %>% dplyr::select(i,i+1))
    tpm[,ncol(tpm) + 1] <- new
    colnames(tpm)[ncol(tpm)] <- paste0("new", i)
}
tpmmean <- tpm[13:18]
colnames(tpmmean) <- gsub("_1","",colnames(tpm %>% dplyr::select(1,3,5,7,9,11)))
#write.csv(tpmmean,"gga.tpmmean.txt",quote = F)
} # Get Gga 7b tpm
vactpm <- read.table("vac.tpmmeanFinalfilt.txt",sep = ",",header = T,row.names = 1)
geneloci <- read.table("vac_zjuV2.final.geneloci",header = F,fill=T)
geneloci$name <- gsub("_[0-9]*","", geneloci$V5)
ggatpm <- read.table("gga.tpmmean.txt",header = T, sep = ",", row.names = 1)
#Zs0
vazs0sub <- subset(geneloci, V1=="chrZ" & V2 > 6.76e6) %>% filter(V5!="")
vazs0tpm <- vactpm[rownames(vactpm) %in% vazs0sub$V4,] %>% filter_all(.,any_vars(.>=1))
ggzs0tpm <- ggatpm[rownames(ggatpm) %in% vazs0sub$name,]
#Va
rownames(vazs0tpm) <- vazs0sub$V5[match(rownames(vazs0tpm),vazs0sub$V4)] 
col_fun <- colorRamp2(c(-2, 0, 2), c("#0000cd", "white", "#cd0000"))
htVa <- Heatmap(t(scale(t(vazs0tpm))),
                 name = "Repeat abundance",
                 col = col_fun,
                 cluster_rows = TRUE,
                 cluster_columns = FALSE,  
                 show_row_names = TRUE,
                 show_column_names = TRUE, 
                 column_title = " ",
                 column_names_rot = 45,
                 border = TRUE,
                 show_row_dend = TRUE,
                 row_km = 2)
#Gg
rownames1 <- rownames(vazs0tpm)[row_order(htVa)$`1`]
rownames2 <- rownames(vazs0tpm)[row_order(htVa)$`2`]
combined_rownames <- gsub("_[0-9]*","",unique(c(rownames1, rownames2)))

ggzs0tpm <- ggzs0tpm[match(combined_rownames,rownames(ggzs0tpm)),]
ggzs0tpm[is.na(ggzs0tpm)] <- 0
ggzs0tpm <- ggzs0tpm %>% filter_all(.,any_vars(.>=1))

Heatmap(t(scale(t(ggzs0tpm))),
                name = "Repeat abundance",
                col = col_fun,
                cluster_rows = FALSE,
                cluster_columns = FALSE,  
                show_row_names = TRUE,
                show_column_names = TRUE, 
                column_title = " ",
                column_names_rot = 45,
                border = TRUE,
                show_row_dend = TRUE)
htGg
} # Va gene expression with Gga7b as reference
{
##data
homoinfo  <- read.table(file = "vac2elgmul.chrLevel.homo.txt")
homoinfo <- homoinfo %>% dplyr::select(1:3,5:7)
homoinfo$fill <- "cccccc"
colnames(homoinfo) <- c("Species_1","Start_1","End_1", "Species_2", "Start_2", "End_2", "fill")
synteny_dual_comparison <- homoinfo
##
synteny_dual_comparison <- subset(synteny_dual_comparison,(synteny_dual_comparison$Species_1=="chr2"|
                                    synteny_dual_comparison$Species_1=="chrZ"|
                                  synteny_dual_comparison$Species_1=="chrW")&
                                    (synteny_dual_comparison$Species_2=="chr3"|
                                       synteny_dual_comparison$Species_2=="chr6"|
                                       synteny_dual_comparison$Species_2=="chr23"))
##assign chr ID
synteny_dual_comparison$Species_1 <- gsub("chr","",synteny_dual_comparison$Species_1)
synteny_dual_comparison$Species_1 <- gsub("2","1",synteny_dual_comparison$Species_1)
synteny_dual_comparison$Species_1 <- gsub("Z","2",synteny_dual_comparison$Species_1)
synteny_dual_comparison$Species_1 <- gsub("W","3",synteny_dual_comparison$Species_1)
synteny_dual_comparison$Species_2 <- gsub("chr","",synteny_dual_comparison$Species_2)
synteny_dual_comparison$Species_2 <- gsub("\\b3","1",synteny_dual_comparison$Species_2)
synteny_dual_comparison$Species_2 <- gsub("6","2",synteny_dual_comparison$Species_2)
synteny_dual_comparison$Species_2 <- gsub("23","3",synteny_dual_comparison$Species_2)
synteny_dual_comparison$Species_1 <- as.numeric(synteny_dual_comparison$Species_1)
synteny_dual_comparison$Species_2 <- as.numeric(synteny_dual_comparison$Species_2)
##assign cols as vac chrs
cols <- data.frame(t(data.frame("1"="FDA8A8", "2"="F8F8D8", "3"="B7D5A1")))
rownames(cols) <- 1:3
cols$V2 <- rownames(cols)
colnames(cols) <- c("col","chr")
synteny_dual_comparison$fill <- cols$col[match(synteny_dual_comparison$Species_1,cols$chr)]
##genome info
Tchrsize <- data.frame(V1=c("1","2","3"),
                       V2=c("1","1","1"),
                       V3=c("284531044","11831237","15022549"))
Tchrsize <- Tchrsize[order(as.numeric(Tchrsize$V1)),]
Tchrsize$V1 <- as.numeric(Tchrsize$V1)
Tchrsize$fill <- "969696"
Tchrsize$species <- "Vac"
Tchrsize$size <- "12"
Tchrsize$color <- "252525"
colnames(Tchrsize) <- c("Chr", "Start", "End", "fill", "species", "size", "color")
#
Qchrsize <- data.frame(V1=c("1","2","3"),
                       V2=c("1","1","1"),
                       V3=c("176344784","121985326","13854160"))
Qchrsize$fill <- "969696"
Qchrsize$species <- "Elgmul"
Qchrsize$size <- "12"
Qchrsize$color <- "252525"
colnames(Qchrsize) <- c("Chr", "Start", "End", "fill", "species", "size", "color")
karyotype_dual_comparison <- rbind(Tchrsize,Qchrsize)
table(karyotype_dual_comparison$species)
karyotype_dual_comparison$Start <- as.numeric(karyotype_dual_comparison$Start)
karyotype_dual_comparison$End <- as.numeric(karyotype_dual_comparison$End)

##visualize
ideogram(karyotype = karyotype_dual_comparison, synteny = synteny_dual_comparison)
convertSVG("chromosome.svg", device = "png")
} # Synteny of Va to Em
{
##data
homoinfo  <- read.table(file = "elgmul2phplat.txt")
homoinfo <- homoinfo %>% dplyr::select(1:3,5:7)
homoinfo$fill <- "cccccc"
colnames(homoinfo) <- c("Species_1","Start_1","End_1", "Species_2", "Start_2", "End_2", "fill")
synteny_dual_comparison <- homoinfo
##
synteny_dual_comparison <- subset(synteny_dual_comparison,(synteny_dual_comparison$Species_1=="chr3"|
                                    synteny_dual_comparison$Species_1=="chr6"|
                                    synteny_dual_comparison$Species_1=="chr23")&
                                    (synteny_dual_comparison$Species_2=="chr2"|
                                       synteny_dual_comparison$Species_2=="chr17"))
##assign chr ID
synteny_dual_comparison$Species_1 <- gsub("chr","",synteny_dual_comparison$Species_1)
synteny_dual_comparison$Species_1 <- gsub("\\b3","1",synteny_dual_comparison$Species_1)
synteny_dual_comparison$Species_1 <- gsub("\\b6","2",synteny_dual_comparison$Species_1)
synteny_dual_comparison$Species_1 <- gsub("\\b23","3",synteny_dual_comparison$Species_1)
synteny_dual_comparison$Species_2 <- gsub("chr","",synteny_dual_comparison$Species_2)
synteny_dual_comparison$Species_2 <- gsub("\\b2","1",synteny_dual_comparison$Species_2)
synteny_dual_comparison$Species_2 <- gsub("\\b17","2",synteny_dual_comparison$Species_2)
synteny_dual_comparison$Species_1 <- as.numeric(synteny_dual_comparison$Species_1)
synteny_dual_comparison$Species_2 <- as.numeric(synteny_dual_comparison$Species_2)
##assign cols as vac chrs
cols <- data.frame(t(data.frame("1"="FDA8A8", "2"="F8F8D8", "3"="B7D5A1")))
rownames(cols) <- 1:3
cols$V2 <- rownames(cols)
colnames(cols) <- c("col","chr")
synteny_dual_comparison$fill <- cols$col[match(synteny_dual_comparison$Species_1,cols$chr)]
##genome info
Tchrsize <- data.frame(V1=c("1","2","3"),
                       V2=c("1","1","1"),
                       V3=c("176344784","121985326","13854160"))
Tchrsize <- Tchrsize[order(as.numeric(Tchrsize$V1)),]
Tchrsize$V1 <- as.numeric(Tchrsize$V1)
Tchrsize$fill <- "969696"
Tchrsize$species <- "Em"
Tchrsize$size <- "12"
Tchrsize$color <- "252525"
colnames(Tchrsize) <- c("Chr", "Start", "End", "fill", "species", "size", "color")
#
Qchrsize <- data.frame(V1=c("1","2"),
                       V2=c("1","1"),
                       V3=c("336734412","8897685"))
Qchrsize$fill <- "969696"
Qchrsize$species <- "Ph"
Qchrsize$size <- "12"
Qchrsize$color <- "252525"
colnames(Qchrsize) <- c("Chr", "Start", "End", "fill", "species", "size", "color")
karyotype_dual_comparison <- rbind(Tchrsize,Qchrsize)
table(karyotype_dual_comparison$species)
karyotype_dual_comparison$Start <- as.numeric(karyotype_dual_comparison$Start)
karyotype_dual_comparison$End <- as.numeric(karyotype_dual_comparison$End)

##visualize
ideogram(karyotype = karyotype_dual_comparison, synteny = synteny_dual_comparison)
convertSVG("chromosome.svg", device = "png")
} # Synteny of Em to Ph
{
{
vacrpt  <- read.table(file = "vac.rpt.txt")
vacsub <- vacrpt %>% dplyr::filter(V1=="chr2") %>%
  setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
  mutate(
    class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
    family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
    length = end - start,
    window = floor(start / 100000) * 100000
  ) %>%
  select(-factor)
massbyclass<- vacsub %>%
  group_by(chr, window, class) %>%
  summarise(total_length = sum(length)) %>%
  ungroup()
##check class
ggplot(massbyclass, aes(x = as.numeric(window), y = total_length, fill = class)) +
  geom_bar(stat = "identity", position = "dodge",alpha=.6) +
  facet_wrap(~ class, scales = "free_y") +
  labs(x = "Window", y = "Total Length", title = "Total Length per 100Kb Window by Class") +
  theme_bw() +
  theme(axis.text.x = element_blank())
##check family
massbyfamily<- vacsub %>% dplyr::filter(class=="LINE") %>%
  group_by(chr, window, family) %>%
  summarise(total_length = sum(length)) %>%
  ungroup()
ggplot(massbyfamily, aes(x = as.numeric(window), y = total_length, fill = family)) +
  geom_bar(stat = "identity", position = "dodge",alpha=.6) +
  facet_wrap(~ family, scales = "free_y") +
  labs(x = "Window", y = "Total Length", title = "Total Length per 100Kb Window by Class") +
  theme_bw() +
  theme(axis.text.x = element_blank())
##check subfamily
massbysubfamily<- vacsub %>% dplyr::filter(family=="L1") %>%
  group_by(chr, window, subfamily) %>%
  summarise(total_length = sum(length)) %>%
  ungroup()
ggplot(massbysubfamily, aes(x = as.numeric(window), y = total_length, fill = subfamily)) +
  geom_bar(stat = "identity", position = "dodge",alpha=.6) +
  facet_wrap(~ subfamily, scales = "free_y") +
  labs(x = "Window", y = "Total Length", title = "Total Length per 100Kb Window by Class") +
  theme_bw() +
  theme(axis.text.x = element_blank())
##barcodes
vacsubForbar <- vacrpt %>% dplyr::filter(V1=="chr2"|V1=="chrZ"|V1=="chrW"|V1=="chr1") %>%
  setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
  mutate(
    class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
    family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
    length = end - start,
    window = floor(start / 100000) * 100000
  ) %>%
  select(-factor)
##major burst of line sine ltr
massofTar <- vacsubForbar %>%
  dplyr::filter(class=="LINE"|class=="LTR"|class=="SINE")  %>% 
  group_by(chr, window) %>%
  summarise(total_length = sum(length)) %>%
  ungroup() %>%
  mutate(
    mass = total_length / max(total_length),
    mid_point = window
  ) %>%
  ungroup()
massofTar$mass_level <- cut(massofTar$mass, breaks = quantile(massofTar$mass, probs = seq(0, 1, length.out = 101)),
                        include.lowest = TRUE, labels = FALSE)
##
massofTar <- massofTar %>%
  mutate(fill_color = case_when(
    mass_level >= 75 ~ "High",
    mass_level < 75 ~ "Other"
  ))
ggplot(massofTar, aes(fill = fill_color)) +
  geom_rect(aes(xmin = window-1e5, xmax = window, ymin = as.numeric(factor(chr))- 0.5, ymax = as.numeric(factor(chr)) + 0.5)) +
  theme_minimal() +
  scale_fill_manual(values = c("black", "white")) + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  guides(fill = "none")  
##repeat density as controls
ggplot(massofTar %>% filter(chr=="chr1" | chr=="chr2"), aes(x=window,y=mass))+
  geom_point(size=.1,alpha=.5)+ theme_classic()+
  facet_wrap(~ chr, scales = "free_x")
} # Va
{
emrpt  <- read.table(file = "em.rep.txt")
emsub <- emrpt %>% dplyr::filter(V1=="chr3") %>%
  setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
  mutate(
    class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
    family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
    length = end - start,
    window = floor(start / 100000) * 100000
  ) %>%
  select(-factor)
#check line step of repeat
changepoint <- emsub %>%
  dplyr::filter(class=="LINE"|class=="LTR"|class=="SINE")  %>% 
  group_by(chr, window) %>%
  summarise(total_length = sum(length)) %>%
  ungroup() %>%
  mutate(
    mass = total_length / max(total_length),
    mid_point = window
  ) %>%
  ungroup()
#changepoint
cpt_mass <- cpt.mean(changepoint$mass)
library(cpm)
df <- data.frame(x = changepoint$window/100000, y = changepoint$mass)
plot(df, col = "steelblue", lwd = 2)
shapiro.test(df$y)

cpm.res = processStream(df$y, cpmType = "Mann-Whitney")

plot(df, col = "steelblue", lwd = 2)
abline(v = cpm.res$changePoints, lwd = 3.5, col = "red")
print(cpm.res$changePoints)
##
changepoints <- 1494
ggplot(changepoint, aes(x = as.numeric(window), y = mass)) +
  geom_step() +
  geom_vline(xintercept = changepoint$window[changepoints], 
             color = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(x = "Window", y = "Mass", title = "Change Point Analysis with Step Plot")
##check class
massbyclass<- emsub %>%
  group_by(chr, window, class) %>%
  summarise(total_length = sum(length)) %>%
  ungroup()
ggplot(massbyclass, aes(x = as.numeric(window), y = total_length, fill = class)) +
  geom_bar(stat = "identity", position = "dodge",alpha=.6) +
  facet_wrap(~ class, scales = "free_y") +
  labs(x = "Window", y = "Total Length", title = "Total Length per 100Kb Window by Class") +
  theme_bw() +
  theme(axis.text.x = element_blank())
##check family
massbyfamily<- emsub %>% dplyr::filter(class=="LINE") %>%
  group_by(chr, window, family) %>%
  summarise(total_length = sum(length)) %>%
  ungroup()
ggplot(massbyfamily, aes(x = as.numeric(window), y = total_length, fill = family)) +
  geom_bar(stat = "identity", position = "dodge",alpha=.6) +
  facet_wrap(~ family, scales = "free_y") +
  labs(x = "Window", y = "Total Length", title = "Total Length per 100Kb Window by Class") +
  theme_bw() +
  theme(axis.text.x = element_blank())
##check subfamily
massbysubfamily<- emsub %>% dplyr::filter(family=="L1") %>%
  group_by(chr, window, subfamily) %>%
  summarise(total_length = sum(length)) %>%
  ungroup()
ggplot(massbysubfamily, aes(x = as.numeric(window), y = total_length, fill = subfamily)) +
  geom_bar(stat = "identity", position = "dodge",alpha=.6) +
  facet_wrap(~ subfamily, scales = "free_y") +
  labs(x = "Window", y = "Total Length", title = "Total Length per 100Kb Window by Class") +
  theme_bw() +
  theme(axis.text.x = element_blank())
##barcodes
emsubForbar <- emrpt %>% dplyr::filter(V1=="chr3"|V1=="chr6"|V1=="chr23"|V1=="chr2") %>%
  setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
  mutate(
    class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
    family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
    length = end - start,
    window = floor(start / 100000) * 100000
  ) %>%
  select(-factor)
##major burst of line sine ltr
massofTar <- emsubForbar %>%
  dplyr::filter(class=="LINE"|class=="LTR"|class=="SINE")  %>% 
  group_by(chr, window) %>%
  summarise(total_length = sum(length)) %>%
  ungroup() %>%
  mutate(
    mass = total_length / max(total_length),
    mid_point = window
  ) %>%
  ungroup()
massofTar$mass_level <- cut(massofTar$mass, breaks = quantile(massofTar$mass, probs = seq(0, 1, length.out = 101)),
                            include.lowest = TRUE, labels = FALSE)
##
massofTar <- massofTar %>%
  mutate(fill_color = case_when(
    mass_level > 75 ~ "High",
    mass_level <= 75 ~ "Other"
  ))
ggplot(massofTar, aes(fill = fill_color)) +
  geom_rect(aes(xmin = window-1e5, xmax = window, ymin = as.numeric(factor(chr))- 0.5, ymax = as.numeric(factor(chr)) + 0.5)) +
  theme_minimal() +
  scale_fill_manual(values = c("black", "white")) + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  guides(fill = "none")  
##repeat density as controls
ggplot(massofTar %>% filter(chr=="chr3" | chr=="chr2"), aes(x=window/1e6,y=mass))+
  geom_point(size=.1,alpha=.5)+ theme_classic()+
  facet_wrap(~ chr, scales = "free_x")
##Check telomeric region
emTELsub <- emrpt %>% dplyr::filter(V1=="chr3"|V1=="chr6"|V1=="chr23"|V1=="chr2") %>%
  setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
  mutate(
    class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
    family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
    length = end - start,
    window = floor(start / 100000) * 100000
  ) %>%
  select(-factor)
massofTEL <- emTELsub %>%
  dplyr::filter(subfamily=="(TTAGGG)n")  %>% 
  group_by(chr, window) %>%
  summarise(total_length = sum(length)) %>%
  ungroup() %>%
  mutate(
    mass = total_length / max(total_length),
    mid_point = window
  ) %>%
  ungroup() %>%
  dplyr::filter(chr=="chr3")
##
ggplot(massofTEL, aes(x = window,y=mass))+
  geom_bar(stat = "identity")
##
massofTEL$mass_level <- cut(massofTEL$mass, breaks = quantile(massofTEL$mass, probs = seq(0, 1, length.out = 101)),
                            include.lowest = TRUE, labels = FALSE)
massofTEL <- massofTEL %>%
  mutate(fill_color = case_when(
    mass_level > 75 ~ "High",
    mass_level <= 75 ~ "Other"
  ))
ggplot(massofTEL, aes(fill = fill_color)) +
  geom_rect(aes(xmin = window-1e5, xmax = window, ymin = as.numeric(factor(chr))- 0.5, ymax = as.numeric(factor(chr)) + 0.5)) +
  theme_minimal() +
  scale_fill_manual(values = c("black", "white")) + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  guides(fill = "none")  
} # Em
{
phrpt  <- read.table(file = "ph.rpt.txt")
phsub <- phrpt %>% dplyr::filter(V1=="chr2") %>%
  setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
  mutate(
    class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
    family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
    length = end - start,
    window = floor(start / 100000) * 100000
  ) %>%
  select(-factor)
massbyclass<- phsub %>%
  group_by(chr, window, class) %>%
  summarise(total_length = sum(length)) %>%
  ungroup()
##check class
ggplot(massbyclass, aes(x = as.numeric(window), y = total_length, fill = class)) +
  geom_bar(stat = "identity", position = "dodge",alpha=.6) +
  facet_wrap(~ class, scales = "free_y") +
  labs(x = "Window", y = "Total Length", title = "Total Length per 100Kb Window by Class") +
  theme_bw() +
  theme(axis.text.x = element_blank())
##check family
massbyfamily<- phsub %>% dplyr::filter(class=="Satellite") %>%
  group_by(chr, window, family) %>%
  summarise(total_length = sum(length)) %>%
  ungroup()
ggplot(massbyfamily, aes(x = as.numeric(window), y = total_length, fill = family)) +
  geom_bar(stat = "identity", position = "dodge",alpha=.6) +
  facet_wrap(~ family, scales = "free_y") +
  labs(x = "Window", y = "Total Length", title = "Total Length per 100Kb Window by Class") +
  theme_bw() +
  theme(axis.text.x = element_blank())
##barcodes
phsubForbar <- phrpt %>% dplyr::filter(V1=="chr2"|V1=="chr17"|V1=="chr1") %>%
  setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
  mutate(
    class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
    family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
    length = end - start,
    window = floor(start / 100000) * 100000
  ) %>%
  select(-factor)
##major burst of line sine ltr
massofTar <- phsubForbar %>%
  dplyr::filter(class=="LINE"|class=="LTR"|class=="SINE")  %>% 
  group_by(chr, window) %>%
  summarise(total_length = sum(length)) %>%
  ungroup() %>%
  mutate(
    mass = total_length / max(total_length),
    mid_point = window
  ) %>%
  ungroup()
massofTar$mass_level <- cut(massofTar$mass, breaks = quantile(massofTar$mass, probs = seq(0, 1, length.out = 101)),
                            include.lowest = TRUE, labels = FALSE)
##
massofTar <- massofTar %>%
  mutate(fill_color = case_when(
    mass_level > 75 ~ "High",
    mass_level <= 75 ~ "Other"
  ))
ggplot(massofTar, aes(fill = fill_color)) +
  geom_rect(aes(xmin = window-1e5, xmax = window, ymin = as.numeric(factor(chr))- 0.5, ymax = as.numeric(factor(chr)) + 0.5)) +
  theme_minimal() +
  scale_fill_manual(values = c("black", "white")) + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  guides(fill = "none")  
##repeat density as controls
ggplot(massofTar %>% filter(chr=="chr1" | chr=="chr2"), aes(x=window/1e6,y=mass))+
  geom_point(size=.1,alpha=.5)+ theme_classic()+
  facet_wrap(~ chr, scales = "free_x")
##Check telomeric region
phTELsub <- phrpt  %>%
  setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
  mutate(
    class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
    family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
    length = end - start,
    window = floor(start / 100000) * 100000
  ) %>%
  select(-factor)
massofTEL <- phTELsub %>%
  dplyr::filter(subfamily=="(TTAGGG)n")  %>% 
  group_by(chr, window) %>%
  summarise(total_length = sum(length)) %>%
  ungroup() %>%
  mutate(
    mass = total_length / max(total_length),
    mid_point = window
  ) %>%
  ungroup() %>% dplyr::filter(chr=="chr1") 
##
ggplot(massofTEL,aes(x=window,y=mass)) +
  geom_step()+ylim(0,max(massofTEL$mass))
  

massofTEL$mass_level <- cut(massofTEL$mass, breaks = quantile(massofTEL$mass, probs = seq(0, 1, length.out = 101)),
                            include.lowest = TRUE, labels = FALSE)
##
massofTEL <- massofTEL %>%
  mutate(fill_color = case_when(
    mass_level > 50 ~ "High",
    mass_level <= 50 ~ "Other"
  ))
ggplot(massofTEL, aes(fill = fill_color)) +
  geom_rect(aes(xmin = window-1e5, xmax = window, ymin = as.numeric(factor(chr))- 0.5, ymax = as.numeric(factor(chr)) + 0.5)) +
  theme_minimal() +
  scale_fill_manual(values = c("black", "white")) + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  guides(fill = "none")  
} # Ph
} # Barcode of TE bundency in Va, Em and Ph
{
  noNbed <- read.table("chr2.noN.100k.bed") # Sliding window cover 100Kb non-N nucleobase
  nuclAlign <- read.table("2toW.delta.filt.coords")
  ##P1. Higher identity of 2q-end region
  idCk <- nuclAlign %>% filter(V7 < 95) %>% 
    mutate(type = case_when(
      V2 < 237e6 ~ "Other region of Chr2",
      V1 >= 237e6 & V2 < 264e6 ~ "2q-end 2",
      V1 >= 264e6 ~ "2q-end 1"
    )) %>% select(V7,type)
  p1 <- ggplot(idCk, aes(x = type, y = V7)) +
    geom_boxplot(width=0.5) +
    stat_compare_means(method = "wilcox.test", 
                       ref.group =  c("2q-end 1"),
                       method.args = list(altanative = "greater")) +  # One-sided test, testing if values in one group are greater
    theme_bw() +
    labs(x = "", y = "Identity")+  
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),  
      axis.text = element_text(size = 12), 
      axis.title = element_text(size = 14, face = "bold"), 
      legend.title = element_text(size = 12), 
      legend.text = element_text(size = 10) 
    )
  ##P2. Line plots of 2q-end
  nuclAlignSub <- nuclAlign %>% filter(V7<95) %>% filter(V1>237e6) %>% select(12,1,2) #Remove too high identity
  grA <- GRanges(seqnames = nuclAlignSub$V12,
                 ranges = IRanges(start = nuclAlignSub$V1, end = nuclAlignSub$V2))
  grB <- GRanges(seqnames = noNbed$V1,
                 ranges = IRanges(start = noNbed$V2, end = noNbed$V3))
  overlaps <- findOverlaps(grB, grA) 
  noNbed$total_length <- 0
  for (i in 1:length(grB)) {
    overlap_indices <- subjectHits(overlaps[queryHits(overlaps) == i]) # Obtain index of overlapped regions
    overlap_length <- sum(width(grA[overlap_indices])) # Sum up length
    noNbed$total_length[i] <- overlap_length 
  }
  noNbed <- noNbed %>% filter(V2>237e6)
  pbar <- ggplot(noNbed, aes(x = V2, xend = V3, y = 1, yend = 1, color = total_length)) + 
    geom_segment(size = 10) +  # 使用 geom_segment 来绘制每个条形
    scale_color_gradient(low = "white", high = "#cd0000") +  # 设置颜色渐变
    labs(x = "Position", y = "", title = "Barcode plot based on total length") + 
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),  
      axis.ticks.y = element_blank(),  
      axis.title.y = element_blank()   
    )
  pline <- ggplot(noNbed, aes(x = V2,y= total_length)) + 
    geom_step() + 
    labs(x = "Position", y = "", title = "Barcode plot based on total length") + 
    theme_minimal() +  
    theme(
      axis.text = element_text(size = 12),  # Increase axis text size
      axis.title = element_text(size = 14, face = "bold"),  # Make axis titles bold
      legend.title = element_text(size = 12),  # Legend title size
      legend.text = element_text(size = 10))  # Legend text size
  p2 <- ggarrange(pbar,pline,ncol=1)
  ##P3. Line plots of repeats
  vacrpt  <- read.table(file = "vac.rpt.txt")
  vacrpt <-  vacrpt %>% 
    filter(V1 == "chr2" & V2 >=237e6) 
  vacsub <- vacrpt  %>%
    setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
    mutate(
      class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
      family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
      length = end - start,
      window = floor(start / 500000) * 500000
    ) %>%
    select(-factor) %>% 
    filter(!str_detect(family, "\\?")) %>% 
    filter(!str_detect(class, "\\?")) %>% 
    filter(class %in% c("SINE","LINE","LTR","DNA"))
  massbyfamily <- vacsub %>%
    group_by(family, window) %>%
    summarise(total_length = sum(length)) %>%
    ungroup() %>% filter(family %in% c("RTE-BovB","L2","RTE-X"))
  p3 <- ggplot(massbyfamily, aes(x = window,y= total_length)) + 
    geom_step(size = .5) + facet_wrap( ~ family, scale = "free_y", ncol =1) +# 设置颜色渐变
    labs(x = "Position", y = "", title = "Barcode plot based on total length") + 
    theme_minimal() +  
    theme(
      axis.text = element_text(size = 12),  # Increase axis text size
      axis.title = element_text(size = 14, face = "bold"),  # Make axis titles bold
      legend.title = element_text(size = 12),  # Legend title size
      legend.text = element_text(size = 10)  # Legend text size
    )
  bovbCk <- massbyfamily %>% filter(family == "RTE-BovB") %>% 
    mutate(type = case_when(
      window < 237e6 ~ "Other region of Chr2",
      window >= 237e6 & window < 264e6 ~ "2q-end 2",
      window >= 264e6 ~ "2q-end 1"
    )) %>% select(total_length, type)
  p4 <- ggplot(bovbCk, aes(x = type, y = total_length)) +
    geom_boxplot(width=0.5) +
    stat_compare_means(method = "wilcox.test", 
                       ref.group =  c("2q-end 1"),
                       method.args = list(altanative = "greater")) +  # One-sided test, testing if values in one group are greater
    theme_bw() +
    labs(x = "", y = "Length")+  
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),  
      axis.text = element_text(size = 12), 
      axis.title = element_text(size = 14, face = "bold"), 
      legend.title = element_text(size = 12), 
      legend.text = element_text(size = 10) 
    )
} # Alignments of zoom-in region of 2q end (The duplication may occured once)
{
  gcdf <- read.table("gc.100k.txt")
  gcdf <- subset(gcdf, V1 %in% c("HiC_scaffold_2","HiC_scaffold_4","HiC_scaffold_18","HiC_scaffold_W"))
  ggplot(gcdf,aes(x=V2, y=V4)) +
    geom_step() + facet_wrap(~ V1, scales = "free_x")
} # GC level see chr2 if non-recombining
{
  #all kimura
  alldiv <- read.table("vac.te.all.kimura")
  alldivSub <- alldiv %>% filter(V5=="LINE/RTE-BovB" & V1 %in% c("chr2","chrW")) %>%
    mutate(
      V1 = ifelse(V1 == "chr2" & V3 > 237e6, "chr2tip", V1),
      V6 = gsub("Kimura", "", V6),
      factor = factor(paste(V1, V4)),
      len = V3-V2
    ) %>% filter(V1 %in% c("chr2tip","chrW","chr2")) 
  
  densityDf <- alldivSub %>%
    mutate(V6_rounded = round(as.numeric(V6))) %>%
    group_by(V1, V6_rounded) %>%
    summarise(total_len = sum(len)) %>% 
    group_by(V1) %>% 
    mutate(total_len = total_len / max(total_len)) %>% ungroup() ##normalized 0-1 
  ##add colors
  allchr_colors <- c(chr2tip = "#4daf4a", chrW = "#cd0000",chr2 ="black")
  ##plot
  library(cpm)
  sub2tip <- densityDf %>% filter(V1 == "chr2tip") %>% select("V6_rounded","total_len") 
  subW <- densityDf %>% filter(V1 == "chrW") %>% select("V6_rounded","total_len") 
  sub2 <- densityDf %>% filter(V1 == "chr2") %>% select("V6_rounded","total_len") 
  
  detectChangePointBatch(sub2$total_len,"Student",alpha=0.05)
  
  ggplot(densityDf, aes(x = V6_rounded, y = total_len, color = V1)) + xlim(0,50) +
    geom_step(size =1) + facet_wrap(~ V1, ncol =1) + 
    theme_classic() +
    scale_color_manual(values = allchr_colors) +
    labs(x = "Kimura distance", y = "Density")
  
  group1 <- densityDf %>% filter(V1 == "chr2tip") %>% pull(total_len)
  group2 <- densityDf %>% filter(V1 == "chrW") %>% pull(total_len)
  
  ks.test(group1, group2)
  
  wilcox.test(group1, group2, alternative = "greater")
  
  
} # Kimura of RTE-BovB on chr2 and chrW
{
  vacrpt  <- read.table(file = "vac.rpt.txt")
  vacsub <- vacrpt %>% dplyr::filter(V1=="chr2") %>%
    setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
    mutate(
      class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
      family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
      length = end - start,
      window = floor(start / 100000) * 100000
    ) %>%
    select(-factor)
  massbyclass<- vacsub %>%
    group_by(chr, window, class) %>%
    summarise(total_length = sum(length)) %>%
    ungroup()
  ##check family
  massbyfamily<- vacsub %>% dplyr::filter(class=="LINE") %>%
    group_by(chr, window, family) %>%
    summarise(total_length = sum(length)) %>%
    ungroup()
  ggplot(massbyfamily, aes(x = as.numeric(window), y = total_length, fill = family)) +
    geom_smooth() +
    facet_wrap(~ family, scales = "free_y") +
    labs(x = "Window", y = "Total Length", title = "Total Length per 100Kb Window by Class") +
    theme_bw() 
  ##check subfamily
  massbysubfamily<- vacsub %>% dplyr::filter(family=="L1") %>%
    group_by(chr, window, subfamily) %>%
    summarise(total_length = sum(length)) %>%
    ungroup()
  ggplot(massbysubfamily, aes(x = as.numeric(window), y = total_length, fill = subfamily)) +
    geom_bar(stat = "identity", position = "dodge",alpha=.6) +
    facet_wrap(~ subfamily, scales = "free_y") +
    labs(x = "Window", y = "Total Length", title = "Total Length per 100Kb Window by Class") +
    theme_bw() +
    theme(axis.text.x = element_blank())
} # RTE-BovB density line on chr2 and chr4
{
{
vadf <- read.table(file = "vac.rpt.txt")
vsdf <- read.table(file = "vsa.rep.txt")
emdf <- read.table(file = "em.rep.txt")
phdf <- read.table(file = "ph.rpt.txt")
varpt <-  vadf %>% 
  filter(!str_detect(V1, "scaf")) %>%
  setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
  mutate(
    class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
    family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
    length = end - start,
    window = floor(start / 100000) * 100000
  ) %>%
  select(-factor) %>% 
  filter(!str_detect(family, "\\?")) %>% 
  filter(!str_detect(class, "\\?")) %>% 
  filter(class == "SINE" | class == "LINE" | class == "LTR") 
vsrpt <-  vsdf %>% 
  setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
  mutate(
    class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
    family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
    length = end - start,
    window = floor(start / 100000) * 100000
  ) %>%
  select(-factor) %>% 
  filter(!str_detect(family, "\\?")) %>% 
  filter(!str_detect(class, "\\?")) %>% 
  filter(class == "SINE" | class == "LINE" | class == "LTR") 
emrpt <-  emdf %>% 
  setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
  mutate(
    class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
    family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
    length = end - start,
    window = floor(start / 100000) * 100000
  ) %>%
  select(-factor) %>% 
  filter(!str_detect(family, "\\?")) %>% 
  filter(!str_detect(class, "\\?")) %>% 
  filter(class == "SINE" | class == "LINE" | class == "LTR")  
phrpt <-  phdf %>% 
  setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
  mutate(
    class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
    family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
    length = end - start,
    window = floor(start / 100000) * 100000
  ) %>%
  select(-factor) %>% 
  filter(!str_detect(family, "\\?")) %>% 
  filter(!str_detect(class, "\\?")) %>% 
  filter(class == "SINE" | class == "LINE" | class == "LTR") 
} # data import
{
#
vaMassZscale <- varpt %>% 
  group_by(chr, window) %>%
  summarise(total_length = sum(length), .groups = "drop")  %>%
  mutate(norm_length = scale(total_length)) %>%  # Zscale
  mutate(species = "Va") %>% filter(chr %in% c("chr1", "chr2")) %>% 
  group_by(chr)  %>%
  mutate(norm_window = window/max(window),
         new_chr = chr) %>%
  ungroup() %>% select(chr, norm_window, norm_length, new_chr, species) %>% 
    mutate(chrtype = ifelse(new_chr == "chr2", "Va Chr2", "Va Chr1"))
#
vaTreeNodel1 <- vaMassZscale %>% filter(new_chr == "chr1") %>% 
  select(new_chr, norm_window, norm_length) %>% 
  mutate(
    predict_tree = predict(rpart(norm_length ~ norm_window, data = .),
                           newdata = .))
vaTreeNodel2 <- vaMassZscale %>% filter(new_chr == "chr2") %>% 
  select(new_chr, norm_window, norm_length) %>% 
  mutate(
    predict_tree = predict(rpart(norm_length ~ norm_window, data = .),
                           newdata = .))
vaTreeNodel <- rbind(vaTreeNodel1, vaTreeNodel2)
ggplot(vaTreeNodel, aes(norm_window*100, norm_length, color = new_chr)) +
  geom_line(aes(y=predict_tree),size = 1)+
  theme_classic()+
  labs(title = "Z-scaled repeat mass", x = "", y ="") 
ggplot(vaMassZscale) + 
  geom_smooth(aes(x = norm_window*100, y = norm_length, color = new_chr), 
              method = "loess", span = 0.2) +
  theme_classic()+
  labs(title = "Z-scaled repeat mass", x = "", y ="") 
#
vsChr1scaf <- c("scaf765.1","scaf755.1","scaf624.1","scaf798.1","scaf779.1","scaf778.1")
vsChr2scaf <- c("scaf827.1","scaf762.1","scaf001.1","scaf795.1","scaf789.1","scaf791.1","scaf703.1","scaf674.1",
                "scaf792.1","scaf794.1","scaf685.1","scaf697.1")
#vsChr2tipscaf <- c("scaf792.1","scaf794.1","scaf685.1","scaf697.1") 
vsMassZscale1 <- vsrpt %>% 
  filter(chr %in% c(vsChr1scaf, vsChr2scaf)) %>%
  group_by(chr, window) %>%
  summarise(total_length = sum(length), .groups = "drop")  %>%
  mutate(norm_length = scale(total_length)) %>%  # Zscale
  mutate(species = "Vs") %>% 
  mutate(
    chr_group = case_when(
      chr %in% vsChr1scaf ~ "chr1",
      chr %in% vsChr2scaf ~ "chr2",
      TRUE ~ "other"
    ) # Separate deal
  ) %>% 
  group_by(chr_group, chr) %>%
  filter(chr_group == "chr1") %>%
  mutate(
    chr = factor(chr, levels = vsChr1scaf)
  ) %>%
  arrange(chr, window) %>%
  ungroup() %>%
  mutate(
    new_window = ((row_number() - 1) * 100000)
  ) 
vsMassZscale2 <- vsrpt %>% 
  filter(chr %in% c(vsChr1scaf, vsChr2scaf)) %>%
  group_by(chr, window) %>%
  summarise(total_length = sum(length), .groups = "drop")  %>%
  mutate(norm_length = scale(total_length)) %>%  # Zscale
  mutate(species = "Vs") %>% 
  mutate(
    chr_group = case_when(
      chr %in% vsChr1scaf ~ "chr1",
      chr %in% vsChr2scaf ~ "chr2",
      TRUE ~ "other"
    ) # Separate deal
  ) %>% 
  group_by(chr_group, chr) %>%
  filter(chr_group == "chr2") %>% 
  mutate(
    chr = factor(chr, levels = vsChr2scaf)
  ) %>%
  arrange(chr, window) %>%
  ungroup() %>%
  mutate(
    new_window = ((row_number() - 1) * 100000)
  )
vsMassZscale <- rbind(vsMassZscale1, vsMassZscale2) %>%
 group_by(chr_group) %>% 
  mutate(norm_window = new_window/max(new_window),
         new_chr = chr_group,
         species = "Vs") %>%
  ungroup() %>% select(chr, norm_window, norm_length, new_chr, species)  %>% 
  mutate(chrtype = ifelse(new_chr == "chr2", "Va Chr2", "Va Chr1"))
#
emMassZscale <- emrpt %>%
  group_by(chr, window) %>%
  summarise(total_length = sum(length), .groups = "drop") %>%
  ungroup() %>%  
  mutate(norm_length = scale(total_length)) %>%  # Zscale
  mutate(species = "Em") %>% 
  filter(chr %in% c("chr2", "chr3", "chr4", "chr6")) %>%
  # Create chr2_4 and chr6_3 combinations
  mutate(new_window = case_when(
    chr == "chr2" ~ window,  # Keep window of chr2 unchanged
    chr == "chr4" ~ window + max(window[chr == "chr2"]),  # Shift chr4 window by max of chr2
    chr == "chr6" ~ max(window[chr == "chr6"]) - window + 1,  # Reverse chr6 window
    chr == "chr3" ~ window + max(window[chr == "chr6"]),  # Shift chr3 window by max of reversed chr6
    TRUE ~ window  # For other cases (although not necessary here)
  )) %>%
  #Create new chr labels
  mutate(new_chr = case_when(
    chr == "chr2" ~ "chr2_4",  # chr2 stays chr2_4
    chr == "chr4" ~ "chr2_4",  # chr4 becomes chr2_4
    chr == "chr3" ~ "chr6_3",  # chr3 becomes chr6_3
    chr == "chr6" ~ "chr6_3",  # chr6 becomes chr6_3
    TRUE ~ chr  # For other cases (although not necessary here)
  )) %>%
  #Normalize the new_window values
  group_by(new_chr) %>%
  mutate(norm_window = new_window / max(new_window)) %>%
  ungroup() %>% select(chr, norm_window, norm_length, new_chr, species) %>% 
  mutate(chrtype = ifelse(new_chr == "chr6_3", "Va Chr2", "Va Chr1"))
#
emTreeNodel1 <- emMassZscale %>% filter(chrtype == "Va Chr1") %>% 
  select(chrtype, norm_window, norm_length) %>% 
  mutate(
    predict_tree = predict(rpart(norm_length ~ norm_window, data = .),
                           newdata = .))
emTreeNodel2 <- emMassZscale %>% filter(chrtype == "Va Chr2") %>% 
  select(chrtype, norm_window, norm_length) %>% 
  mutate(
    predict_tree = predict(rpart(norm_length ~ norm_window, data = .),
                           newdata = .))
emTreeNodel <- rbind(emTreeNodel1, emTreeNodel2)
ggplot(emTreeNodel, aes(norm_window*100, norm_length, color = chrtype)) +
  geom_line(aes(y=predict_tree),size = 1)+
  theme_classic()+
  labs(title = "Z-scaled repeat mass", x = "", y ="") 
ggplot(emMassZscale) + 
  geom_smooth(aes(x = norm_window*100, y = norm_length, color = new_chr), 
              method = "loess", span = 0.2) +
  theme_classic()+
  labs(title = "Z-scaled repeat mass", x = "", y ="") 

#
phMassZscale <- phrpt %>%
  group_by(chr, window) %>%
  summarise(total_length = sum(length), .groups = "drop")  %>%
  mutate(norm_length = scale(total_length)) %>%  # Zscale
  ungroup() %>%  
  mutate(species = "Ph") %>% filter(chr %in% c("chr1", "chr2")) %>% 
  group_by(chr)  %>%
  mutate(norm_window = window/max(window),
         new_chr = chr) %>%
  ungroup() %>% select(chr, norm_window, norm_length, new_chr, species) %>% 
  mutate(chrtype = ifelse(new_chr == "chr2", "Va Chr2", "Va Chr1"))
#
phTreeNodel1 <- phMassZscale %>% filter(chrtype == "Va Chr1") %>% 
  select(chrtype, norm_window, norm_length) %>% 
  mutate(
    predict_tree = predict(rpart(norm_length ~ norm_window, data = .),
                           newdata = .))
phTreeNodel2 <- phMassZscale %>% filter(chrtype == "Va Chr2") %>% 
  select(chrtype, norm_window, norm_length) %>% 
  mutate(
    predict_tree = predict(rpart(norm_length ~ norm_window, data = .),
                           newdata = .))
phTreeNodel <- rbind(phTreeNodel1, phTreeNodel2)
ggplot(phTreeNodel, aes(norm_window*100, norm_length, color = chrtype)) +
  geom_line(aes(y=predict_tree),size = 1)+
  theme_classic()+
  labs(title = "Z-scaled repeat mass", x = "", y ="") 
ggplot(phMassZscale) + 
  geom_smooth(aes(x = norm_window*100, y = norm_length, color = new_chr), 
              method = "loess", span = 0.2) +
  theme_classic()+
  labs(title = "Z-scaled repeat mass", x = "", y ="") 
#
tip2merge <- rbind(vaMassZscale, vsMassZscale, emMassZscale, phMassZscale)  %>% 
  select(norm_window, norm_length, chrtype, species)
tip2merge$species <- factor(tip2merge$species, levels = c("Va", "Vs", "Em", "Ph"))
ggplot(tip2merge) + 
  geom_smooth(aes(x = norm_window*100, y = norm_length, color = chrtype), 
              method = "loess", span = 0.2) + facet_wrap(~ species) +
  theme_classic()+
  labs(title = "Z-scaled repeat mass", x = "", y ="") + 
  scale_color_manual(name = "Homologous to",
                     values = c("Va Chr1"="#4daf4a","Va Chr2"="#ff7700"))
} # Scaled retroposon density in 2 Chr
{
tip2merge_split <- tip2merge %>%
  group_by(species) %>%
  group_split()
result <- data.frame() #init blank df
for (species_data in tip2merge_split) {
  Va_Chr1 <- species_data %>% filter(chrtype == "Va Chr1")
  Va_Chr2 <- species_data %>% filter(chrtype == "Va Chr2")
  fit_Va_Chr1 <- loess(norm_length ~ norm_window, data = Va_Chr1, span = 0.2)
  fit_Va_Chr2 <- loess(norm_length ~ norm_window, data = Va_Chr2, span = 0.2)
  pred_Va_Chr1 <- predict(fit_Va_Chr1) 
  pred_Va_Chr2 <- predict(fit_Va_Chr2) # get smooth data
  diff_pred <- pred_Va_Chr2 - pred_Va_Chr1 # get Chr2 repeat dens - Chr1 repeat dens
  #bind data
  temp_result <- data.frame(
    norm_window = Va_Chr1$norm_window, 
    diff_norm_length = diff_pred,
    species = Va_Chr1$species[1]  
  )
  result <- bind_rows(result, temp_result)
}
result$species <- factor(result$species, levels = c("Va", "Vs", "Em", "Ph"))
ggplot(result) +
  geom_smooth(aes(x = norm_window * 100, y = diff_norm_length, color = species), 
              method = "lm", 
              formula = y ~ poly(x, 20), 
              span = 0.1) +
  facet_wrap(~ species) + 
  theme_classic() +
  labs(title = "Repeat(Va Chr2) - Repeat(Va Chr1)", 
       x = "Relative Genomic position", 
       y = "Scaled repeat length Diff.") +
  scale_color_manual(values = c("Va" = "#cd0000",
                                "Vs" = "black",
                                "Em" = "#4daf4a", 
                                "Ph" = "#0000cd"),
                     name = "Species")


} # Difference of two Chr
{
subset1 <- subset(vaMassZscale, chrtype=="Va Chr2" & norm_window >0.83)
subset2 <- subset(vsMassZscale, chr %in% c("scaf792.1","scaf794.1","scaf685.1","scaf697.1"))
subset3 <- subset(emMassZscale, chrtype=="Va Chr2" & norm_window >0.83)
subset4 <- subset(phMassZscale, chrtype=="Va Chr2" & norm_window >0.83)

boxplotmergd <- rbind(subset1,subset2,subset3,subset4) %>% select(species,norm_length)
boxplotmergd$species <- factor(boxplotmergd$species, levels = c("Va", "Vs", "Em", "Ph"))
ggplot(boxplotmergd, aes(x = species, y = norm_length)) +
  geom_boxplot() +
  theme_classic() +
  stat_compare_means(comparisons = list(
    c("Va", "Vs"),
    c("Vs", "Em"),
    c("Em", "Ph")
  ), 
  method = "wilcox.test", 
  label = "p.signif") +  # Use 'label = "p.signif"' for significance markers (*, **, etc.)
  labs(x = "Species", y = "Scaled repeat mass", title = "")

} # Boxplot of scaled 2q region mass
} # Relative enriched Retroposon in Va 2q end
{
  v2rloci  <- read.table("v2rVerifiedNR.txt") # V2R loci
  v2eHomoinfo  <- read.table("vac2elgmul.eachBestVac.txt") # One-way Vac best protein hit on Emu
  emuV2R <- read.table("emu.v2r.pep") # Emu V2R protein IDs
  geneloci <- read.table("vac_zjuV2.final.geneloci",header = F,fill=T) # Vac gene loci
  disruptORF <- read.table("vac.disruptORF.txt") # Vac gene with disrupted ORF
  ohnologs <- read.table("ggallus.Pairs.Intermediate.2R.txt",fill = T) # Ohnologs of chicken
  vacchrlen <- read.table(file = "vac.chrLen.txt") # Chromosome length of Vac
  #Orthofinder output
  duploci  <- read.table("VaVsEmPpTi.orth.group.txt",fill = T) # with em and pp and ti
  dupCount <- read.table("Orthogroups.GeneCount.tsv",fill = T, header = T) # with em and pp and ti
  duploci <- duploci %>% filter(V4 !="")
  #Vsa scaffolds assignments
  vsaChr2scaf <- c("scaf827.1","scaf762.1","scaf001.1","scaf795.1","scaf789.1","scaf791.1","scaf703.1","scaf674.1",
                   "scaf792.1","scaf794.1","scaf685.1","scaf697.1")
  vsaChr4scaf <- c("scaf654.1","scaf782.1","scaf781.1","scaf753.1")
  vsaZscaf <- c("scaf802.1")
  vsaWscaf <- c("scaf729.1","scaf715.1","scaf856.1","scaf746.1")
} # Data input of paralogs analysis
{
scAndchr2Groups <- duploci %>%
  filter(V3 == "vac") %>%
  mutate(
    type = ifelse(V1 %in% v2rloci$V4, "V2R", "Other"),
    GeneProduct = gsub("_[0-9]*","",geneloci$V5[match(V1, paste(geneloci$V4,"-T1",sep = ""))]),
    GeneProduct = ifelse(GeneProduct == "" & type == "V2R", "V2R", GeneProduct),
    GeneProduct = ifelse(GeneProduct == "" ,"Unknown", GeneProduct)
  ) %>% filter(!(GeneProduct %in% toupper(ohnologs$V3) | GeneProduct %in% toupper(ohnologs$V4)))  %>% # Remove Gga Ohnologs
  group_by(V4) %>%
  filter(any(V2 == "chr2") & any(V2 == "chrW")) %>%
  ungroup() } # Ortho Groups with Vac chr2 and chrW dups 
{
dfForBarcodeChr2 <- scAndchr2Groups %>% filter(V2 == "chr2") %>%
  mutate(
    st = geneloci$V2[match(V1, paste(geneloci$V4,"-T1",sep = ""))],
    ed = geneloci$V3[match(V1, paste(geneloci$V4,"-T1",sep = ""))]
  ) 
noNbed <- read.table("chr2.noN.100k.bed") # Sliding window cover 100Kb non-N nucleobase
noNbed <- noNbed %>% mutate(window = paste(V2, V3, sep = "-"))
grA <- GRanges(seqnames = dfForBarcodeChr2$V2,
               ranges = IRanges(start = dfForBarcodeChr2$st, end = dfForBarcodeChr2$ed))
grB <- GRanges(seqnames = noNbed$V1,
               ranges = IRanges(start = noNbed$V2, end = noNbed$V3))
overlaps <- findOverlaps(grB, grA) 
dfForBarcodeChr2$window <- NA 
dfForBarcodeChr2$window[overlaps@to] <- paste0(noNbed$V2[overlaps@from], "-", noNbed$V3[overlaps@from])
dupfreq <- dfForBarcodeChr2  %>%
  group_by(window) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(noNbed, by = c("window" = "window")) %>% mutate(
    width = V3 - V2,
    mid = V2 + width / 2
  ) # Count gene frequency of each window 
breaks = unique(c(dupfreq$V2, dupfreq$V3))
ggplot(dupfreq, aes(x = V2/1e6, xend = V3/1e6, y = 1, yend = 1, color = pmin(count, 5))) + 
  geom_segment(size = 10) +  
  scale_color_gradient(low = "white", high = "#cd0000") +  
  labs(x = "Chr2", y = "", title = "Gene Frequency") + 
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(),  
    axis.title.y = element_blank()   
  ) } # Barcodes for duplication counts located at Chr2 in 100kb non-N window
{
dupOnChr2tip <- scAndchr2Groups %>% filter(V2 == "chr2") %>%
  mutate(
    st = geneloci$V2[match(V1, paste(geneloci$V4,"-T1",sep = ""))],
    ed = geneloci$V3[match(V1, paste(geneloci$V4,"-T1",sep = ""))]
  ) %>% filter(st > 237e6) %>% select(V2, st, ed, V1) %>%
  rename(
    "V2" = "chr",
    "V1" = "gene"
  )
dupOnChrW <- scAndchr2Groups %>% filter(V2 == "chrW") %>%
  mutate(
    st = geneloci$V2[match(V1, paste(geneloci$V4,"-T1",sep = ""))],
    ed = geneloci$V3[match(V1, paste(geneloci$V4,"-T1",sep = ""))]
  ) %>% select(V2, st, ed, V1) %>%
  rename(
    "V2" = "chr",
    "V1" = "gene"
  )
dupInfo <- rbind(dupOnChr2tip,dupOnChrW)
#write.table(dupInfo, file = "vac.dupChr2ChrW.txt", quote = F,row.names = F, col.names = F,sep = "\t")
tecount <- read.table(file = "vac.nearby.count.te.prop.txt") 
genecount <- read.table(file = "vac.nearby.count.gene.prop.txt") 
tecount <- tecount %>% filter(grepl("2tip|w",V4)) %>% filter(grepl("LINE|SINE|LTR",V2))
genecount <- genecount %>% filter(grepl("2tip|w",V4)) %>% filter(grepl("LINE|SINE|LTR",V2))
tecount <- tecount %>% mutate(
  factor = paste(V1,V2,V4,sep = ":"),
  ratioTE = (V3/V5)*100
)
genecount <- genecount %>% mutate(
  factor = paste(V1,V2,V4,sep = ":"),
  ratioGene = (V3/V5)*100
)
targetTE <- c("rnd-5_family-305","rnd-5_family-1095") #Define target TE sub-familiy
finalDf <- merge(tecount[, c("factor", "ratioTE")], 
                 genecount[, c("factor", "ratioGene")], 
                 by = "factor")
finalDf <- finalDf %>%
  mutate(factor = as.character(factor)) %>%
  separate(factor, into = c("subfamily", "info", "location"), sep = ":") %>%
  separate(info, into = c("class", "family"), sep = "/") %>%
  mutate(type = ifelse(subfamily %in% targetTE, "target","other"))
wdata <- finalDf %>% filter(location=="w.te"|location=="w.ctl.te") %>% mutate(
  data = "Wcopy"
)
twodata <- finalDf %>% filter(location=="2tip.te"|location=="2tip.ctl.te") %>% mutate(
  data = "2copy"
)
mergedf <- rbind(wdata,twodata) %>%
  mutate(origin = ifelse(location=="w.te"|location=="2tip.te", "target","control"))
ggplot(mergedf,aes(x=ratioTE,y=ratioGene,color=type,shape=origin,alpha = type))+
  geom_point() + theme_classic() +
  scale_color_manual(values = c("#bdbdbd", "black")) + 
  geom_text_repel(
    data = subset(mergedf, type == "target"), 
    aes(label = subfamily),  
    size = 3,  
    box.padding = 0.5, 
    point.padding = 0.5,
    max.overlaps = 100
  ) + labs(x = "% TE within 10k of genes", y = "% Gene with specific TE within 10k")+
  facet_wrap(~ data) } # Repeat nearby those duplicates (That may regulate the segmental duplication)
{
dupCountsub <- dupCount %>%
  filter(Orthogroup %in% scAndchr2Groups$V4) %>%
  filter(Emu > 0 | Ppl > 0 | Til > 0) # Select orthogroups which had at least one ortholog in out group species
duplociVacSub <- scAndchr2Groups %>% filter(V4  %in% dupCountsub$Orthogroup)
#Assign genes to each Group
propGene <- duplociVacSub %>% group_by(V4,GeneProduct) %>%
  summarise(freq = n(), .groups = "drop") 
PieDonut(propGene, aes(pies = V4, donuts = GeneProduct , count = freq)) 
groupAssign <- propGene %>%
  group_by(V4) %>%
  mutate(total_freq = sum(freq)) %>%  
  mutate(freq_ratio = freq / total_freq) %>%  
  arrange(V4, desc(freq_ratio)) %>%  
  slice_max(order_by = freq_ratio, n = 1) %>% filter(V4!="OG0000067") %>% # Remove OG0000067, it was NFI genes which is onholog
  filter(GeneProduct != "Unknown")  # Remove OG0000010 and OG0000218, it was unknown function gene in out group species also
OGassign <- groupAssign %>% select(1,2)
#Gene Locations and production in Vac
dfForBubble <- duplociVacSub %>% 
  filter(V4 %in% OGassign$V4) %>% 
  mutate(
    Gene = OGassign$GeneProduct[match(V4, OGassign$V4)]
  ) %>% 
  mutate(chr = ifelse(V2 %in% c("chr2","chrZ","chrW","chr4"), V2, "Unassigned scaffolds")) %>%
  group_by(chr,Gene) %>%
  summarise(freq = n(), .groups = "drop") 
ggplot(dfForBubble, aes(x = chr, y = Gene)) +
  geom_point(aes(size = freq, color = Gene))+ 
  theme_bw() +
  labs(x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic"),  
    legend.title = element_text("Paralogs Count") 
  ) + 
  scale_size_continuous(
    name = "Paralogs Count", 
    breaks = c(50, 100, 200),  
    limits = c(0, max(dfForBubble$freq)),  
    trans = "identity"  
  )
#Expression of paralogs
dfForExp <- duplociVacSub %>% 
  filter(V4 %in% OGassign$V4) %>% 
  mutate(
    Gene = OGassign$GeneProduct[match(V4, OGassign$V4)]
  ) %>% 
  mutate(chr = ifelse(V2 %in% c("chr2","chrZ","chrW","chr4"), V2, "Unassigned scaffolds"))
vactpm <- read.table("vac.tpmmeanFinalfilt.txt",sep = ",",header = T,row.names = 1)
geneloci <- read.table("vac_zjuV2.final.geneloci",header = F,fill=T)
v2rexp <- vactpm[rownames(vactpm) %in% 
                   gsub("-T[0-9]*","",subset(dfForExp, GeneProduct == "V2R")$V1),] %>% 
  filter_all(.,any_vars(.>=1))
hlaexp <- vactpm[rownames(vactpm) %in% 
                   gsub("-T[0-9]*","",subset(dfForExp, GeneProduct == "HLA-DPB1")$V1),] %>% 
  filter_all(.,any_vars(.>=1))
rnfexp <- vactpm[rownames(vactpm) %in% 
                   gsub("-T[0-9]*","",subset(dfForExp, GeneProduct == "EEF2")$V1),] %>% 
  filter_all(.,any_vars(.>=1))
library(pheatmap)
tarexp <- rnfexp
anno_row <- data.frame(chr = geneloci$V1[match(rownames(tarexp), geneloci$V4)])
rownames(anno_row) <- rownames(tarexp)
pheatmap(log2(tarexp+1),
         cluster_rows = TRUE,       
         cluster_cols = F,       
         color = colorRampPalette(c("white", "#cd0000"))(50),  
         fontsize = 10,
         annotation_row = anno_row) } # What are those duplicates and expressions
{
dfForChr2GeneBarcode <- duplociVacSub %>% 
  filter(V4 %in% OGassign$V4) %>% 
  mutate(
    Gene = OGassign$GeneProduct[match(V4, OGassign$V4)]
  ) %>%  select(1,7) %>% 
  mutate(
    chr = geneloci$V1[match(V1, paste(geneloci$V4,"-T1",sep = ""))],
    st = geneloci$V2[match(V1, paste(geneloci$V4,"-T1",sep = ""))],
    ed = geneloci$V3[match(V1, paste(geneloci$V4,"-T1",sep = ""))]
  ) %>% filter(!grepl("^scaf", chr)) %>% # Remove genes on unasssigned scaffolds
  mutate(
    chrLen = vacchrlen$V2[match(chr, vacchrlen$V1)],
    stNorm = st/chrLen,
    edNorm = ed/chrLen,
    Gene = ifelse(V1 %in% disruptORF$V1, "Disrupted", Gene)
  )
ggplot(dfForChr2GeneBarcode) +
  geom_segment(aes(x = (stNorm + edNorm)/2, xend = (stNorm + edNorm)/2, y = 1, yend = 2, color = Gene)) +
  geom_vline(xintercept = 237e6/289531044) + facet_wrap( ~ chr, ncol = 1) +
  theme_bw() + 
  labs(x = "Relative position", y = "Chromosome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_manual(values = c("V2R" = "#ff7f00", "RNF39" = "#CD0000", "HLA-DPB1" = "#4DAF4A", "EEF2" = "#0000CD",
                                "Disrupted" = "black"))
#How about the scaffolds with duplicates, do they highly-enriched with RTE-BovB?
scaffoldsCandi <-  duplociVacSub %>% 
  filter(V4 %in% OGassign$V4) %>% 
  mutate(
    Gene = OGassign$GeneProduct[match(V4, OGassign$V4)]
  ) %>% filter(grepl("^scaf", V2)) 
vacrpt  <- read.table(file = "vac.rpt.txt")
vacsub <- vacrpt %>% filter(grepl("^scaf", V1)) %>%
  setNames(c("chr", "start", "end", "subfamily", "factor")) %>%  
  mutate(
    class = ifelse(grepl("/", factor), sub("/.*", "", factor), factor),
    family = ifelse(grepl("/", factor), sub(".*/", "", factor), factor),
    length = end - start,
    window = floor(start / 100000) * 100000
  ) %>%
  select(-factor)
##check family
massbyfamily<- vacsub %>% dplyr::filter(class=="LINE") %>%
  group_by(chr, window, family) %>%
  summarise(total_length = sum(length)) %>%
  ungroup() %>% mutate(
    scafType = ifelse(chr %in% scaffoldsCandi$V2, "Scaffolds with duplicates","Other")
  ) %>% filter(family %in% c("RTE-BovB", "L2", "RTE-X"))

ggplot(massbyfamily, aes(x = scafType, y = total_length)) + 
  geom_boxplot() + facet_wrap( ~ family) + 
  stat_compare_means(comparisons = list(
    c("Scaffolds with duplicates", "Other")), 
    method = "t.test", 
    label = "p.signif") + 
  labs(x = "Unassigned Scaffolds", y = "Repeat length in 100k window (bp)", title = "") + theme_bw()
} # Specific location of involved orthogroups in Vac
{
#How about Vs
duplociVsaSub <- duploci %>%
  filter(V3 == "vsa") %>% 
  filter(V4 %in% OGassign$V4) %>% mutate(
    Gene = OGassign$GeneProduct[match(V4, OGassign$V4)],
    Chr = case_when(
      V2 %in% vsaChr2scaf ~ "Chr2",
      V2 %in% vsaChr4scaf ~ "Chr4",
      V2 %in% vsaZscaf ~ "ChrZ",
      V2 %in% vsaWscaf ~ "ChrW",
      !(V2 %in% c(vsaChr2scaf, vsaChr4scaf, vsaZscaf, vsaWscaf)) ~ "Other")) %>%
  group_by(Chr, Gene) %>%
  summarise(freq = n(), .groups = "drop") 
ggplot(duplociVsaSub, aes(x = Chr, y = Gene)) +
  geom_point(aes(size = freq, color = Gene))+ 
  theme_bw() +
  labs(x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic"),  
    legend.title = element_text("Paralogs Count") 
  ) + 
  scale_size_continuous(
    name = "Paralogs Count", 
    breaks = c(50, 100, 200),  
    limits = c(0, max(dfForBubble$freq)),  
    trans = "identity"  
  )
#See chromosome location of out group species
outSpeDup <- duploci %>%
  filter(V3 != "vac" & V3 != "vsa") %>% 
  filter(V4 %in% OGassign$V4) %>% mutate(
    Gene = OGassign$GeneProduct[match(V4, OGassign$V4)]
  )
dfForBubbleOut <- outSpeDup %>%
  group_by(V2, V3, Gene) %>%
  summarise(freq = n(), .groups = "drop") %>% 
  mutate(
    VacHomo = case_when(
      V3 == "em" & V2 == "chr23" ~ "Sex",
      V3 == "em" & V2 == "chr11" ~ "Chr4",
      V3 == "em" & (V2 == "chr3" | V2 == "chr6") ~ "Chr2",
      V3 == "em" & !(V2 %in% c("chr11", "chr23", "chr3", "chr6")) ~ "Other Chr",
      V3 == "pp" & V2 == "chr6" ~ "Chr4",
      V3 == "pp" & V2 == "chr2" ~ "Chr2",
      V3 == "pp" & V2 == "chr17" ~ "Sex",
      V3 == "til" & V2 == "chr8" ~ "Sex",
      V3 == "til" & V2 == "chr5" ~ "Chr4",
      V3 == "til" & V2 == "chr2" ~ "Chr2",
      TRUE ~ NA_character_  
    )
  )
#For RNF39 and HLA-DPB1 all from homologous chromosome to Vac Chr2. So the Vac-W copy comes from Chr2
#For EEF2, all species both have Sex chromosome copy and Chr2 copy. So we can't define the orientarion.
#For V2R, Two cases, if W and Chr2 copies duplicate independently, they will form 2 clades in gene tree, if recent translocated, that would be a chimeric tree.
ggplot(dfForBubbleOut, aes(x = VacHomo, y = Gene)) +
  geom_point(aes(size = freq, color = Gene))+ facet_wrap( ~ V3) +
  theme_bw() +
  labs(x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic"),  
    legend.title = element_text("Paralogs Count") 
  ) + 
  scale_size_continuous(
    name = "Paralogs Count", 
    breaks = c(50, 100, 200),  
    limits = c(0, max(dfForBubbleOut$freq)),  
    trans = "identity"  
  )
#Merged all species
dfvac <- dfForBubble %>% mutate(chr = ifelse(chr == "Unassigned scaffolds","chr2", chr)) %>% 
  filter(grepl("^chr", chr)) %>% mutate(spe = "Vac") %>% 
  group_by(chr, Gene, spe) %>%
  summarize(freq = sum(freq), .groups = "drop")
dfvsa <- duplociVsaSub %>% mutate(Chr = ifelse(Chr == "Other","Chr2", Chr)) %>% 
  filter(grepl("^Chr", Chr)) %>% rename("Chr" = "chr") %>% mutate(chr = gsub("Chr", "chr", chr), spe = "Vsa") %>% 
  group_by(chr, Gene, spe) %>%
  summarize(freq = sum(freq), .groups = "drop")
dfoutSpe <- dfForBubbleOut %>% select(5, 3, 4, 2) %>%
  rename("VacHomo" = "chr", "V3" = "spe") %>% 
  mutate(chr = gsub("Chr", "chr", chr)) %>% 
  mutate(chr = gsub("Sex", "chrZ", chr)) %>% 
  filter(grepl("^chr", chr)) 
dfForAllbubble <- rbind(dfvac, dfvsa, dfoutSpe) 
ggplot(dfForAllbubble, aes(x = factor(spe, levels = c("til", "pp", "em", "Vsa", "Vac")), y = freq)) +
  geom_bar(aes(fill = chr), stat = "identity")+ facet_wrap( ~ Gene, nrow =1 , scales = "free_y") + theme_bw()+ 
  scale_fill_manual(values = c("chr4" = "#4daf4a", "chr2" = "#ff7f00", "chrW" = "#cd0000", "chrZ" = "#0000cd"))
ggplot(dfForAllbubble, aes(x = chr, y = Gene)) +
  geom_point(aes(size = freq, color = Gene))+ facet_wrap( ~ spe, nrow =1) +
  theme_bw() +
  labs(x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic"),  
    legend.title = element_text("Paralogs Count") 
  ) + 
  scale_size_continuous(
    name = "Paralogs Count", 
    breaks = c(50, 100, 200),  
    limits = c(0, max(dfForBubbleOut$freq)),  
    trans = "identity"  
  )
#Get duplication OG list
VaOGlist <- dupCount %>% filter(Vac > 0) %>% 
  filter(Vsa > 0 | Emu > 0 | Ppl > 0 | Til > 0) # Select orthogroups which had at least one ortholog in out group species
VacGeneMultiOG <- dupCount %>% filter(Vac >= 3) %>% 
  filter(Vsa > 0 | Emu > 0 | Ppl > 0 | Til > 0)
wGeneInOG <- duploci %>% filter(V4 %in% VaOGlist$Orthogroup) %>%
  filter(V3 == "vac" & V2 =="chrW") %>% mutate(
    copy = ifelse(V4 %in% VacGeneMultiOG$Orthogroup, "Multi", "Single")
  ) 
wGeneAll <- geneloci %>% filter(V1 == "chrW") %>% 
  mutate(V4 = paste0(V4, "-T1", sep ="")) %>% select(V4) %>% mutate(
    OG = duploci$V4[match(V4, duploci$V1)],
    type = wGeneInOG$copy[match(V4, wGeneInOG$V1)],
    type = ifelse(is.na(type), "Ab initio", type),
    type = ifelse(OG %in% c("OG0000000","OG0000042"), "MultiV2R", type),
    orf = ifelse(V4 %in% disruptORF$V1, "Disrupt", "Intact")
  ) %>% group_by(type, orf) %>%
  summarise(freq = n(), .groups = "drop") 
PieDonut(wGeneAll, aes(pies = type,
                       donuts = orf , count = freq)) } # Paralogs duplications in Vsa, Emu, Ppl and Til
{
library(ggtree)
tr1 = read.tree("v2r.fasttree.OG.nwk") # all V2Rs from OG0000000
#fasttree
treeanno <- read.table("v2rOG.anno.txt",fill = T)
treeanno <- treeanno %>% mutate(
  label =  V1,
  species = gsub("[0-9].*","", V1),
  chr = case_when(
    species == "vac" & V3 == "chr2" ~ "Chr2",
    species == "vac" & V3 == "chr4" ~ "Chr4",
    species == "vac" & (V3 == "chrW" | V3 == "chrZ")  ~ "Sex",
    species == "vac" & !(V3 %in% c("chr2", "chr4", "chrZ", "chrW")) ~ "Other",
    species == "vsa" & V3 %in% vsaChr2scaf ~ "Chr2",
    species == "vsa" & V3 %in% vsaChr4scaf ~ "Chr4",
    species == "vsa" & (V3 %in% vsaZscaf | V3 %in% vsaWscaf ) ~ "Sex",
    species == "vsa" & !(V3 %in% c(vsaChr2scaf, vsaChr4scaf, vsaZscaf, vsaWscaf)) ~ "Other",
    species == "em" & V3 == "chr23" ~ "Sex",
    species == "em" & V3 == "chr11" ~ "Chr4",
    species == "em" & (V3 == "chr3" | V3 == "chr6") ~ "Chr2",
    species == "em" & !(V3 %in% c("chr11", "chr23", "chr3", "chr6")) ~ "Other",
    species == "pp" & V3 == "chr6" ~ "Chr4",
    species == "pp" & V3 == "chr2" ~ "Chr2",
    species == "pp" & V3 == "chr17" ~ "Sex",
    species == "til" & V3 == "chr8" ~ "Sex",
    species == "til" & V3 == "chr5" ~ "Chr4",
    species == "til" & V3 == "chr2" ~ "Chr2",
    TRUE ~ NA_character_ ),
  orf = ifelse(V2 %in% disruptORF$V1, "Disrupted", NA)
) %>% select(label, species, chr, orf)  
rownames(treeanno) <- treeanno$label
treeanno <- treeanno %>% select(-label) %>% select(-orf)
p1 <- ggtree(tr1, layout = "fan", size=.2, open.angle= 10, branch.length = "none") 
gheatmap(p1, treeanno, offset = 1, width=.3) + scale_fill_manual(
  values = c("Sex" = "#cd0000", "Chr4" = "#4daf4a", "Chr2" = "#ff7f00", "Other" = "grey",
             "vac" = "#cd0000", "vsa"="#cd0000", "em" = "#0000cd",
             "pp"="#4daf4a", "til" = "black"))
#Other three genes
trHla <- read.tree("hla.txt")
trRnf39 <- read.tree("rnf39.txt")
trEef2 <- read.tree("eef2.txt")
allGeneAnno <- read.table("allGeneOG.anno.txt",fill = T)
allGeneAnno <- allGeneAnno %>% mutate(
  species = gsub("[0-9].*","", V2),
  species = case_when(
    species == "em" ~ "Emu",
    species == "pp" ~ "Ppl",
    species == "til" ~ "Til",
    species == "vac" ~ "Vac",
    species == "vsa" ~ "Vsa"),
  label = paste(species, V2, sep = "_"),
  chr = case_when(
    species == "Vac" & V3 == "chr2" ~ "Chr2",
    species == "Vac" & V3 == "chr4" ~ "Chr4",
    species == "Vac" & (V3 == "chrW" | V3 == "chrZ")  ~ "Sex",
    species == "Vac" & !(V3 %in% c("chr2", "chr4", "chrZ", "chrW")) ~ "Other",
    species == "Vsa" & V3 %in% vsaChr2scaf ~ "Chr2",
    species == "Vsa" & V3 %in% vsaChr4scaf ~ "Chr4",
    species == "Vsa" & (V3 %in% vsaZscaf | V3 %in% vsaWscaf ) ~ "Sex",
    species == "Vsa" & !(V3 %in% c(vsaChr2scaf, vsaChr4scaf, vsaZscaf, vsaWscaf)) ~ "Other",
    species == "Emu" & V3 == "chr23" ~ "Sex",
    species == "Emu" & V3 == "chr11" ~ "Chr4",
    species == "Emu" & (V3 == "chr3" | V3 == "chr6") ~ "Chr2",
    species == "Emu" & !(V3 %in% c("chr11", "chr23", "chr3", "chr6")) ~ "Other",
    species == "Ppl" & V3 == "chr6" ~ "Chr4",
    species == "Ppl" & V3 == "chr2" ~ "Chr2",
    species == "Ppl" & V3 == "chr17" ~ "Sex",
    species == "Til" & V3 == "chr8" ~ "Sex",
    species == "Til" & V3 == "chr5" ~ "Chr4",
    species == "Til" & V3 == "chr2" ~ "Chr2",
    TRUE ~ NA_character_ )) %>% select(label, species, chr)  
rownames(allGeneAnno) <- allGeneAnno$label
allGeneAnno <- allGeneAnno %>% select(-label) 
p1 <- ggtree(trEef2, layout = "fan", size=.2, open.angle= 10, branch.length = "none") 
gheatmap(p1, allGeneAnno, offset = 1, width=.3) + scale_fill_manual(
  values = c("Sex" = "#cd0000",  
             "Chr4" = "#4daf4a", "Chr2" = "#ff7f00", "Other" = "grey",
             "Vac" = "#cd0000", "Vsa"="#cd0000", "Emu" = "#0000cd",
             "Ppl"="#4daf4a", "Til" = "black",
             "Disrupted" = "black"))
} # Gene tree for V2R, HLA-DPB1, RNF39 and EEF2
{
#all kimura
alldiv <- read.table("vac.te.all.kimura")
alldiv <- alldiv %>% 
  mutate(
    V1 = ifelse(V1 == "chr2" & V3 > 237e6, "chr2tip", V1),
    V6 = gsub("Kimura", "", V6),
    factor = factor(paste(V1, V4)),
    len = V3-V2
  )
#plots by class
plot_function <- function(type) {
  targetdf <- subset(alldiv, V5 == type) %>% 
    filter(V1 == "chr1" | V1 == "chr2tip" | V1 == "chrZ" |
           V1 == "chr2" | V1 == "chrW")
  ##
  targetdf <- targetdf %>%
    mutate(V6_rounded = round(as.numeric(V6)))
  ##
  targetdf_filtered <- targetdf %>%
    group_by(V1, V6_rounded) %>%
    summarise(total_len = sum(len)) %>%
    ungroup()
  ##normalized by chr len
  targetdf_filtered <- targetdf_filtered %>%
    mutate(total_len = case_when(
      V1 == "chr1" ~ total_len / 305726166,
      V1 == "chr2" ~ total_len / 237000000,
      V1 == "chr2tip" ~ total_len / (284531044 - 237000000),
      V1 == "chrW" ~ total_len / 15022549,
      V1 == "chrZ" ~ total_len / 11831237,
      TRUE ~ total_len
    ))
  ##add colors
  allchr_colors <- c("chr2tip" = "#4daf4a", "chrW" = "#cd0000", 
                     "chr2" = "#004a08", "chrZ" = "#0000cd",
                     "chr1" = "black")
  ##plot
  ggplot(targetdf_filtered, aes(x = V6_rounded, y = total_len, fill = V1)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = allchr_colors) +
    theme_minimal() +
    ggtitle(type) +
    theme(legend.position = "none")
}
plot_list <- lapply(unique(alldiv$V5), plot_function)
ggarrange(plotlist = plot_list, ncol = 1, nrow = length(plot_list))
##BovB
test <- subset(alldiv, V4=="rnd-5_family-305")  %>% 
  filter(V1 == "chr1" | V1 == "chr2tip" | V1 == "chrZ" |
           V1 == "chr2" | V1 == "chrW")
test <- test %>%
  mutate(V6_rounded = round(as.numeric(V6)))
##
test_filtered <- test %>%
  group_by(V1, V6_rounded) %>%
  summarise(total_len = sum(len)) %>%
  ungroup()
##normalized by chr len
test_filtered <- test_filtered %>%
  mutate(total_len = case_when(
    V1 == "chr1" ~ total_len / 305726166,
    V1 == "chr2" ~ total_len / 237000000,
    V1 == "chr2tip" ~ total_len / (284531044 - 237000000),
    V1 == "chrW" ~ total_len / 15022549,
    V1 == "chrZ" ~ total_len / 11831237,
    TRUE ~ total_len
  ))
##add colors
allchr_colors <- c("chr2tip" = "#4daf4a", "chrW" = "#cd0000", 
                   "chr2" = "#004a08", "chrZ" = "#0000cd",
                   "chr1" = "black")
##plot
ggplot(test_filtered, aes(x = V6_rounded, y = total_len, fill = V1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = allchr_colors) +
  theme_minimal() +
  theme(legend.position = "none")
} # TE kimura and copy length
{
allcountFam <- read.table(file="mergeByfamily.txt",header = T,row.names = 1) #2tip
allcountSubFam <- read.table(file="mergeBysubfamily.txt",header = T,row.names = 1) #2tip
allcountFam <- read.table(file="vac.te.wholeGeno.mergeByfamily.txt",header = T,row.names = 1) #2tip
allcountSubFam <- read.table(file="vac.te.wholeGeno.mergeBysubfamily.txt",header = T,row.names = 1) #2tip

family2class <- read.table(file = "vac.repeat.family2class.txt")
#For all family
targetct <- allcountFam %>% filter(rownames(.) != "Unknown") %>%
  filter(!grepl("\\?", rownames(.)))
#CPM
sample_sums <- colSums(targetct)
cpm <- as.data.frame(t(t(targetct) / sample_sums * 1e6))
#get mean of 2 reps
for (i in seq(from=1, to=8, by=2)) {
  new <- rowMeans(cpm %>% dplyr::select(i,i+1))
  cpm[,ncol(cpm) + 1] <- new
  colnames(cpm)[ncol(cpm)] <- paste0("new", i)
}
cpmmean <- cpm[9:12]
colnames(cpmmean) <- gsub("rep1","",colnames(cpm %>% dplyr::select(1,3,5,7)))
colnames(cpmmean) <- c("Female-Brain","Ovary","Male-Brain","Testis")
cpmmean <- cpmmean %>% dplyr::filter_all(.,any_vars(.>=1))
create_heatmap <- function(data) {
  log2data <- log2(data+1)
  classAnno <- family2class[match(rownames(log2data), rownames(family2class)), ]
  annotationsRow <- rowAnnotation(
    Class = classAnno, 
    col = list(class = c("DNA" = "#0000cd", "LINE" = "#cd0000", 
                         "SINE" = "#4daf4a", "LTR" = "#984ea3", 
                         "Other" = "#ff7f00"))
  )
  col_fun <-  c(colorRampPalette(c("#0000cd", "white"))(10),
                colorRampPalette(c("white", "#cd0000"))(10))
  heatmap <- Heatmap(log2data,
                     name = "Repeat abundance",
                     col = col_fun,
                     cluster_rows = TRUE,  
                     cluster_columns = FALSE,  
                     show_row_names = TRUE, 
                     show_column_names = TRUE, 
                     right_annotation = annotationsRow,
                     column_title = " ",
                     column_names_rot = 45,
                     row_names_rot = -45,
                     border = TRUE)
  return(heatmap)
}
create_heatmap(cpmmean)
##Only RTE BovB
agetint <- data.frame(
  family = c("rnd-4_family-45", "rnd-5_family-305", "rnd-2_family-52", "rnd-6_family-57", 
             "rnd-5_family-1095", "rnd-5_family-2055", "rnd-1_family-43", "rnd-1_family-102", 
             "rnd-3_family-21", "rnd-3_family-796", "rnd-6_family-2990", "rnd-4_family-292", 
             "rnd-5_family-127", "rnd-6_family-54", "rnd-1_family-163", "rnd-2_family-44", 
             "rnd-1_family-71", "rnd-3_family-457", "rnd-1_family-45"),
  age = c(rep("young", 10), rep("old", 9)) 
)
targetctSub <- allcountSubFam
# CPM
sample_sums <- colSums(targetctSub)
cpm <- as.data.frame(t(t(targetctSub) / sample_sums * 1e6))
#get mean of 2 reps
for (i in seq(from=1, to=8, by=2)) {
  new <- rowMeans(cpm %>% dplyr::select(i,i+1))
  cpm[,ncol(cpm) + 1] <- new
  colnames(cpm)[ncol(cpm)] <- paste0("new", i)
}
cpmmean <- cpm[9:12]
colnames(cpmmean) <- gsub("rep1","",colnames(cpm %>% dplyr::select(1,3,5,7)))
colnames(cpmmean) <- c("Female-Brain","Ovary","Male-Brain","Testis")
cpmmean <- cpmmean %>% filter(rownames(.) %in% agetint$family) 
create_heatmap <- function(data) {
  log2data <- log2(data+1)
  classAnno <- family2class[match(rownames(log2data), rownames(family2class)), ]
  annotationsRow <- rowAnnotation(
    Class = classAnno, 
    col = list(class = c("DNA" = "#0000cd", "LINE" = "#cd0000", 
                         "SINE" = "#4daf4a", "LTR" = "#984ea3", 
                         "Other" = "#ff7f00"))
  )
  col_fun <-  c(colorRampPalette(c("#0000cd", "white"))(10),
                colorRampPalette(c("white", "#cd0000"))(10))
  heatmap <- Heatmap(log2data,
                     name = "Repeat abundance",
                     col = col_fun,
                     cluster_rows = TRUE,  
                     cluster_columns = FALSE,  
                     show_row_names = TRUE, 
                     show_column_names = TRUE, 
                     right_annotation = annotationsRow,
                     column_title = " ",
                     column_names_rot = 45,
                     row_names_rot = -45,
                     border = TRUE)
  return(heatmap)
}
create_heatmap(cpmmean)
} # TE expression
{
{
geneloci <- read.table("vac.transcript.loci.txt",header = F,fill=T)
disruptORF <- read.table(file = "vac.disruptORF.txt")
v2rloci  <- read.table(file = "v2rVerifiedNR.txt")
v2rFalseORF <- data.frame(
  chr = c("chr2", "chr2", "chr2", "chr4", "chr4", "chr4", "chr4", "chrW", "chrW", "chrW", "chrW", "chrW", "chrW", "chrW", "chrW", "chrW"),
  start = c(255092986, 274581408, 280468408, 334636, 436331, 2758805, 10478135, 4130074, 5758187, 8857316, 9878671, 9919866, 9995389, 10119379, 10920924, 14855038),
  end = c(255093453, 274582195, 280469287, 394531, 438129, 2759119, 10482293, 4133725, 5758531, 8871161, 9879839, 9923711, 9999097, 10120915, 10922001, 14861092)
)  %>%  mutate(
  factor = paste(chr,start,end,sep = "_")
)
v2rloci <- v2rloci %>% mutate(
  chr = case_when(
    V1 %in% c("chr2", "chr4", "chrZ", "chrW") ~ V1,
    grepl("scaf", V1) ~ "Unassigned",
    TRUE ~ "Other Chr"
  ),
  factor = paste(V1,V2,V3,sep = "_")
) %>% filter(V6>100) %>% mutate(
  orf = ifelse(factor %in% v2rFalseORF$factor, "Disrupted", "Intact")
)
genelocisub <- geneloci %>% mutate(
  v2r = ifelse(V4 %in% v2rloci$V4, "v2r", "other"),
  orf = ifelse(V4 %in% disruptORF$V1, "Disrupt", "Intact"),
  chr = case_when(
    V1 %in% c("chr2", "chr4", "chrZ", "chrW") ~ V1,
    TRUE ~ "Other Chr"),
  factor = paste(chr, v2r)) %>% 
  group_by(chr, orf, factor) %>%
  summarise(freq = n(), .groups = "drop")
ggplot(genelocisub, aes(x = factor, y = freq, fill = orf)) +
  geom_bar(position="fill", stat="identity") +   theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "Chr_V2R", y = "Frequency", fill = "ORF Status", 
       title = "Stacked Barplot of Chr_V2R by ORF") + 
  scale_fill_manual(values = c(Intact = "#4daf4a", Disrupt = "#bdbdbd")) 
#
df <- genelocisub %>% filter(!grepl("Other Chr", chr))
results <- df %>%
  group_by(chr) %>%
  summarize(
    chi_test = list({
      counts <- freq
      row_labels <- unique(factor)
      contingency_table <- matrix(
        counts,
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
          c(row_labels[1], row_labels[2]),
          c("Disrupted", "Intact")
        )
      )
      chisq.test(contingency_table)
    })
  )
results <- results %>%
  mutate(
    p_value = sapply(chi_test, function(x) x$p.value),
    X_squared = sapply(chi_test, function(x) x$statistic)
  )
print(results)

} # Va
{
vsageneloci <- read.table("vsa.pep.loci")
vsadisrupt <-  read.table("vsa.disrupt.list")
vsachr2tipscaf <- c("scaf792.1","scaf794.1","scaf685.1","scaf697.1")
vsachr4scaf <- c("scaf654.1","scaf782.1","scaf781.1","scaf753.1")
vsaZscaf <- c("scaf802.1")
vsaWscaf <- c("scaf729.1","scaf715.1","scaf856.1","scaf746.1")  
vsv2rloci <- read.table(file = "vsa.v2r.txt")
vsv2rsub <- vsv2rloci %>% filter(V1 %in% c(vsachr2tipscaf,vsachr4scaf,vsaZscaf,vsaWscaf)) %>% mutate(
  type = case_when(V1 %in% vsachr2tipscaf ~ "VaChr2",
                   V1 %in% vsachr4scaf ~ "VaChr4",
                   V1 %in% vsaZscaf ~ "VaChrZ",
                   V1 %in% vsaWscaf ~ "VaChrW"))
#
vsainfodf <- vsageneloci %>% mutate(
  v2r = ifelse(V4 %in% vsv2rsub$V4, "V2R", "Other"),
  orf = ifelse(V4 %in% vsadisrupt$V1, "Disrupted", "Intact"),
  chr = case_when(V1 %in% vsachr2tipscaf ~ "Chr2",
                  V1 %in% vsachr4scaf ~ "Chr4",
                  V1 %in% vsaZscaf ~ "ChrZ",
                  V1 %in% vsaWscaf ~ "ChrW")
) %>% replace_na(list(v2r = "Other", orf = "Other", chr = "Other")) 
summary_table <- vsainfodf %>%
  mutate(chr_v2r = paste(chr, v2r, sep = "_")) %>%  
  count(chr_v2r, orf, name = "count")
  
ggplot(summary_table, aes(x = chr_v2r, y = count, fill = orf)) +
  geom_bar(position="fill", stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "Chr_V2R", y = "Frequency", fill = "ORF Status", 
       title = "Stacked Barplot of Chr_V2R by ORF") + 
  scale_fill_manual(values = c("Disrupted"="#bdbdbd","Intact"="#4daf4a"))
#
df <- summary_table %>%
  mutate(chr = sub("_.*", "", chr_v2r)) %>% filter(!grepl("Other", chr))
results <- df %>%
  group_by(chr) %>%
  summarize(
    chi_test = list({
      counts <- count
      row_labels <- unique(chr_v2r)
      contingency_table <- matrix(
        counts,
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
          c(row_labels[1], row_labels[2]),
          c("Disrupted", "Intact")
        )
      )
      chisq.test(contingency_table)
    })
  )
results <- results %>%
  mutate(
    p_value = sapply(chi_test, function(x) x$p.value),
    X_squared = sapply(chi_test, function(x) x$statistic)
  )
print(results)
#
} # Vs
} # Proportion of V2R with disrupted ORFs
{
geneloci <- read.table("vac_zjuV2.final.geneloci",header = F,fill=T)
disruptORF <- read.table(file = "vac.disruptORF.txt")
#define single copy of multicopy genes
copynumber <- geneloci %>% select(V1, V5) %>% mutate(
  product = toupper(gsub("_[0-9]*","", V5))
) %>% filter(!grepl("scaf", V1)) %>% filter(product != "") %>% 
  group_by(product) %>%
  summarise(allV1 = paste(V1, collapse = ",")) %>%
  ungroup() %>%
  mutate(
    type = case_when(
      str_count(allV1, ",") == 0 ~ "singlecopy", 
      allV1 %in% c("chrW,chrZ", "chrZ,chrW") ~ "singlecopy", 
      TRUE ~ "multicopy"
    )
  )
#V2R genes
v2rloci  <- read.table(file = "v2rVerifiedNR.txt")
v2rFalseORF <- data.frame(
  chr = c("chr2", "chr2", "chr2", "chr4", "chr4", "chr4", "chr4", "chrW", "chrW", "chrW", "chrW", "chrW", "chrW", "chrW", "chrW", "chrW"),
  start = c(255092986, 274581408, 280468408, 334636, 436331, 2758805, 10478135, 4130074, 5758187, 8857316, 9878671, 9919866, 9995389, 10119379, 10920924, 14855038),
  end = c(255093453, 274582195, 280469287, 394531, 438129, 2759119, 10482293, 4133725, 5758531, 8871161, 9879839, 9923711, 9999097, 10120915, 10922001, 14861092)
)  %>%  mutate(
  factor = paste(chr,start,end,sep = "_")
)
v2rloci <- v2rloci %>% mutate(
  chr = case_when(
    V1 %in% c("chr2", "chr4", "chrZ", "chrW") ~ V1,
    grepl("scaf", V1) ~ "Unassigned",
    TRUE ~ "Other Chr"
  ),
  factor = paste(V1,V2,V3,sep = "_")
) %>% filter(V6>100) %>% mutate(
  orf = ifelse(factor %in% v2rFalseORF$factor, "Disrupted", "Intact")
) 
v2rgeneW <- subset(v2rloci, V1=="chrW") %>% mutate(
  gene = gsub("-T[0-9]*","",V4)) %>% mutate(
    type = "multicopy",
    product = "V2R"
  )

#W gene numbers
wgenes <- geneloci %>% filter(V1 == "chrW") %>% mutate(
  orf = ifelse(V4 %in% gsub("-T[0-9]*", "", disruptORF$V1), "Disrupted", "Intact"),
  product = ifelse(V4 %in% v2rgeneW$gene, "V2R", toupper(gsub("_[0-9]*", "", V5))),
  type = ifelse(V4 %in% v2rgeneW$gene, "V2R", 
                coalesce(copynumber$type[match(product, copynumber$product)], "Ab initio"))) 
t1 <- subset(wgenes, orf == "Disrupted") 
t2 <- subset(wgenes, orf == "Intact")
data1 <- data.frame(table(t1$type)
)
data2 <- data.frame(table(t2$type)
)
pie1 <- ggplot(data1, aes(x = "", y = Freq, fill = 
                            factor(Var1, levels = c("Ab initio","V2R", "multicopy","singlecopy")))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y")  +
  geom_text(aes(label = Freq), 
            position = position_stack(vjust = 0.5))+
  theme_void() + scale_fill_manual(values = c("singlecopy" = "#CD0000", "multicopy"="#0000cd",
                                              "V2R"="#4faf4a","Ab initio" = "#bdbdbd")) +
  labs(title = "Disrupted ORF") 
pie2 <- ggplot(data2, aes(x = "", y = Freq, fill = 
                            factor(Var1, levels = c("Ab initio","V2R", "multicopy","singlecopy")))) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = Freq), 
            position = position_stack(vjust = 0.5))+
  coord_polar(theta = "y") +
  theme_void() + scale_fill_manual(values = c("singlecopy" = "#CD0000", "multicopy"="#0000cd",
                                              "V2R"="#4faf4a","Ab initio" = "#bdbdbd")) +
  labs(title = "Intact ORF")
ggarrange(pie1,pie2, legend = "bottom")
    
} # W gene numbers donuts
{
  dnds <- read.table(file = "va2em.all.dnds.txt")
  dnds <- dnds %>% mutate(
    gene = gsub("-T[0-9]*","", V1)
  )
  geneproduct <- read.table("vac_zjuV2.final.geneloci",header = F,fill=T)
  HLADPB1 <- duplociVacSub %>% filter(GeneProduct == "HLA-DPB1")
  RNF39 <- duplociVacSub %>% filter(GeneProduct == "RNF39")
  MR1 <- geneproduct %>% filter(grepl("\\bMR1", geneproduct$V5)) %>% filter(!grepl("scaf",V1))
  singlecopy <- geneproduct %>% filter(!grepl("_[0-9]", geneproduct$V5)) %>% 
    filter(!grepl("scaf",V1)) %>% filter(V5 !="")
  V2R <-  duplociVacSub %>% filter(GeneProduct == "V2R")
  ##
  dndssub <- dnds %>% mutate(
    type = case_when(
      V1 %in% V2R$V1 ~ "V2R",
      V1 %in% HLADPB1$V1 ~ "HLADPB1",
      V1 %in% RNF39$V1 ~ "RNF39",
      gene %in% MR1$V4 ~ "MR1",
      gene %in% singlecopy$V4 ~ "SingleCopy",
      TRUE ~ "Other"
    )) %>% filter(type != "Other") %>% mutate(
      chr = geneloci$V1[match(gsub("-T[0-9]*","",V1), geneloci$V4)],
      orf = ifelse(V1 %in% disruptORF$V1, "Disrupted", "Intact"),
    ) %>% filter(!grepl("scaf", chr)) %>% filter(!(type == "V2R" & !(chr %in% c("chr2", "chr4", "chrZ", "chrW")))) %>% mutate(
      type = case_when(
        type == "V2R" & chr == "chrW" ~ "V2R_W",
        type == "V2R" & chr == "chr2" ~ "V2R_2",
        type == "V2R" & chr == "chr4" ~ "V2R_4",
        type == "V2R" & chr == "chrZ" ~ "V2R_Z", 
        TRUE ~ type
      )
    ) %>% filter(orf == "Intact")
  ggplot(dndssub, aes(x = type, y = V2, color = type)) +
    geom_boxplot() + theme_classic() + coord_cartesian(ylim = c(0,1)) +
    labs(x="", y= "dN/dS", title = "") 
  wilcox.test(subset(dndssub, type == "V2R_4")$V2, subset(dndssub, type == "MR1")$V2)
  wilcox.test(subset(dndssub, type == "V2R_4")$V2, subset(dndssub, type == "HLADPB1")$V2)
  wilcox.test(subset(dndssub, type == "V2R_4")$V2, subset(dndssub, type == "RNF39")$V2)
  wilcox.test(subset(dndssub, type == "V2R_4")$V2, subset(dndssub, type == "SingleCopy")$V2)
} # V2R genes dnds





