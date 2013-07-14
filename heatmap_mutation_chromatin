library(gdata)
library(GenomicRanges)
library(rtracklayer)
library(RColorBrewer)
library(gplots)


NCBI37=read.table("NCBI37.txt")
names(NCBI37)=c("chr","coordinate")
total=sum(as.numeric(NCBI37$coordinate))

# import cancer data
Gto_lift_over_COLO_829_result =  import.bed("Gto_lift_over_COLO_829_result.bed")
Gto_lift_over_NCI_H209_result = import.bed("Gto_lift_over_NCI_H209_result.bed")       
Gto_lift_over_LEUKEMIA_1= import.bed( "Gto_lift_over_LEUKEMIA_1"   )                   
Gto_lift_over_LEUKEMIA_2=  import.bed( "Gto_lift_over_LEUKEMIA_2"  )                    
Gto_lift_over_LEUKEMIA_3=  import.bed( "Gto_lift_over_LEUKEMIA_3"  )                    
Gto_lift_over_LEUKEMIA_4=  import.bed( "Gto_lift_over_LEUKEMIA_4" )                     
GPM_P_1_result = import.bed("GPM_P_1_result.bed")                      
GPM_P_2_result =import.bed( "GPM_P_2_result.bed")                           
GPM_P_3_result= import.bed( "GPM_P_3_result.bed")                           
GPM_P_4_result= import.bed( "GPM_P_4_result.bed")                            
GPM_P_5_result =import.bed( "GPM_P_5_result.bed")                           
GPM_P_6_result =import.bed( "GPM_P_6_result.bed")                           
GPM_P_7_result= import.bed( "GPM_P_7_result.bed")                           

XPoint_Mutations_list=list()
XPoint_Mutations_list=list(Gto_lift_over_COLO_829_result,Gto_lift_over_NCI_H209_result,
Gto_lift_over_LEUKEMIA_1,Gto_lift_over_LEUKEMIA_2,Gto_lift_over_LEUKEMIA_3,Gto_lift_over_LEUKEMIA_4,
GPM_P_1_result,GPM_P_2_result ,GPM_P_3_result,GPM_P_4_result,GPM_P_5_result,GPM_P_6_result,GPM_P_7_result)
GPoint_Mutations_list=list()
for(i in 1:13){ GPoint_Mutations_list[[i]]=lapply(XPoint_Mutations_list,function(x) as(x,"GRanges"))[[i]]}
LEUKEMIA=c(GPoint_Mutations_list[[3]],GPoint_Mutations_list[[4]],GPoint_Mutations_list[[5]],GPoint_Mutations_list[[6]])
PROSTATE=c(GPoint_Mutations_list[[7]],GPoint_Mutations_list[[8]],GPoint_Mutations_list[[9]],GPoint_Mutations_list[[10]],
GPoint_Mutations_list[[11]],GPoint_Mutations_list[[12]],GPoint_Mutations_list[[13]])

GPoint_Mutations_list_pulled=list()
GPoint_Mutations_list_pulled= list(GPoint_Mutations_list[[1]],GPoint_Mutations_list[[2]],LEUKEMIA,PROSTATE)

#checks
unlist(lapply(GPoint_Mutations_list_pulled, function(x) length(x)))
unlist(lapply(GPoint_Mutations_list, function(x) length(x)))
sum(unlist(lapply(GPoint_Mutations_list_pulled, function(x) length(x))))
GBP_P= GPoint_Mutations_list_pulled 
 
#change to ensembl annotation
for ( i in 1: length(GBP_P)){
seqlevels(GBP_P[[i]]) <- sub("chr23", "chrX", seqlevels(GBP_P[[i]]))
seqlevels(GBP_P[[i]]) <- sub("chr24", "chrY", seqlevels(GBP_P[[i]]))}


#load chromatin regions (`colors of chromatin`). colors.rda contains GRanges objects that represent each chromatin region ( 8 in total). In addition this was calculated for 6 different cell lines (CL1-CL6). For clarity purposes, CL5 is represented here.
load("colors.rda")


colors=list(CL5_T,CL5_E, CL5_R, CL5_CTCF, CL5_PF, CL5_TSS,CL5_WE)
df_chromothripsis_color_5=list()

for ( paciente in 1:length(GBP_P)){

snps_1000_wo_colors=c();snps_1000_wo_colors=unlist(lapply(colors, function(x) length(queryHits(findOverlaps( GBP_P[[paciente]],x,ignore.strand=TRUE)))))

a_w=c(snps_1000_wo_colors,length(GBP_P[[paciente]])- sum(snps_1000_wo_colors))

b_w=length(GBP_P[[paciente]])-a_w

c_w=c(unlist(lapply(colors,function(x) sum(as.numeric(width(x))))),total-sum(unlist(lapply(colors,function(x) sum(as.numeric(width(x)))))))

d_w=total-c_w

matrix_snps_colors=list()

df_snps_colors1=c()

df_snps_colors2=c()

for( i in 1:8){

df_snps_colors1[i]=as.data.frame(matrix(c(a_w[i],b_w[i],sum(a_w[i],b_w[i]))))

df_snps_colors2[i]=as.data.frame(matrix(c(c_w[i],d_w[i],sum(c_w[i], d_w[i]))))

matrix_snps_colors[[i]]=matrix(c(a_w[i],b_w[i],c_w[i],d_w[i]),ncol=2)}

v_matrix_snps_colors=list();
chisq.value=list()

for( i in 1:8){

v_matrix_snps_colors[[i]]=c( chisq.test(matrix_snps_colors[[i]])$p.value,0,0)         

chisq.value[[i]]=c(c( chisq.test(matrix_snps_colors[[i]]))$observed[1,1]-c( chisq.test(matrix_snps_colors[[i]]))$expected[1,1],0,0)}

cat("PACIENTE!!!",paciente, "\n")

df_chromothripsis_color_5[[paciente]]=(as.data.frame(cbind(unlist(df_snps_colors1),unlist(df_snps_colors2),unlist(v_matrix_snps_colors), unlist(chisq.value)))) }

for( i in 1:paciente){ df_chromothripsis_color_5[[i]][,5]=NA } ; empty_vector=list()
for( i in 1:paciente){ empty_vector[[i]]=rep(NA,paciente)}
for( i in 1:paciente){
empty_vector[[i]][which(p.adjust(df_chromothripsis_color_5[[i]][,3][seq(1,24,3)],"BH")<0.05 & df_chromothripsis_color_5[[i]][,4][seq(1,24,3)]>0)]=4 # red 0 This gonna be green
empty_vector[[i]][which(p.adjust(df_chromothripsis_color_5[[i]][,3][seq(1,24,3)],"BH")<0.05 & df_chromothripsis_color_5[[i]][,4][seq(1,24,3)]<0)]=0 # black 2 This gonna be red
empty_vector[[i]][which(p.adjust(df_chromothripsis_color_5[[i]][,3][seq(1,24,3)],"BH")>0.05)]=2} # green 4. This gonna be black
Mdf_chromothripsis_color_5 =matrix(unlist(empty_vector)[is.na(unlist(empty_vector))!="TRUE"],ncol=paciente)



bindedM=rbind(Mdf_chromothripsis_color_1[1:8,]  , Mdf_chromothripsis_color_2[1:8,],Mdf_chromothripsis_color_3[1:8,],
Mdf_chromothripsis_color_4[1:8,]  , Mdf_chromothripsis_color_5[1:8,],Mdf_chromothripsis_color_6[1:8,])

M1=matrix(unlist(lapply(df_chromothripsis_color_1, function(x) x[,1][seq(1,24,3)])),ncol=length(GBP_P))
M2=matrix(unlist(lapply(df_chromothripsis_color_2, function(x) x[,1][seq(1,24,3)])),ncol=length(GBP_P))
M3=matrix(unlist(lapply(df_chromothripsis_color_3, function(x) x[,1][seq(1,24,3)])),ncol=length(GBP_P))
M4=matrix(unlist(lapply(df_chromothripsis_color_4, function(x) x[,1][seq(1,24,3)])),ncol=length(GBP_P))
M5=matrix(unlist(lapply(df_chromothripsis_color_5, function(x) x[,1][seq(1,24,3)])),ncol=length(GBP_P))
M6=matrix(unlist(lapply(df_chromothripsis_color_6, function(x) x[,1][seq(1,24,3)])),ncol=length(GBP_P))


#MT=rbind(M1)
MT=rbind(M1,M2,M3,M4,M5,M6)
MT=matrix(MT,ncol=length(GBP_P))

bindedM[MT<5]=2

chi_M1= matrix(unlist(lapply(df_chromothripsis_color_1, function(x) round(x[,1][seq(1,24,3)]/(x[,1][seq(1,24,3)]-x[,4][seq(1,24,3)]),3))),ncol=length(GBP_P))
chi_M2= matrix(unlist(lapply(df_chromothripsis_color_2, function(x) round(x[,1][seq(1,24,3)]/(x[,1][seq(1,24,3)]-x[,4][seq(1,24,3)]),3))),ncol=length(GBP_P))
chi_M3= matrix(unlist(lapply(df_chromothripsis_color_3, function(x) round(x[,1][seq(1,24,3)]/(x[,1][seq(1,24,3)]-x[,4][seq(1,24,3)]),3))),ncol=length(GBP_P))
chi_M4= matrix(unlist(lapply(df_chromothripsis_color_4, function(x) round(x[,1][seq(1,24,3)]/(x[,1][seq(1,24,3)]-x[,4][seq(1,24,3)]),3))),ncol=length(GBP_P))
chi_M5= matrix(unlist(lapply(df_chromothripsis_color_5, function(x) round(x[,1][seq(1,24,3)]/(x[,1][seq(1,24,3)]-x[,4][seq(1,24,3)]),3))),ncol=length(GBP_P))
chi_M6= matrix(unlist(lapply(df_chromothripsis_color_6, function(x) round(x[,1][seq(1,24,3)]/(x[,1][seq(1,24,3)]-x[,4][seq(1,24,3)]),3))),ncol=length(GBP_P))

M1=paste(paste("R=",chi_M1,sep=""), paste(M1,")",sep=""),sep=" (")
M2=paste(paste("R=",chi_M2,sep=""), paste(M2,")",sep=""),sep=" (")
M3=paste(paste("R=",chi_M3,sep=""), paste(M3,")",sep=""),sep=" (")
M4=paste(paste("R=",chi_M4,sep=""), paste(M4,")",sep=""),sep=" (")
M5=paste(paste("R=",chi_M5,sep=""), paste(M5,")",sep=""),sep=" (")
M6=paste(paste("R=",chi_M6,sep=""), paste(M6,")",sep=""),sep=" (")

MT=rbind(M1,M2,M3,M4,M5,M6)

#MT=rbind(M1)
MT=matrix(MT,ncol=length(GBP_P))

chromo_names=c("Melanoma","Lung cancer ","Leukemia","Prostate cancer")

save(bindedM,MT,file="last.check.colors.rda")

chromatinsegmentation=c(

paste(c("Transcribed","Enhancer", "Repressed", "CTCF", "Promoter","TSS","WeakEnhancer","Black"),"Gm12878",sep="_"),

paste(c("Transcribed","Enhancer", "Repressed", "CTCF", "Promoter","TSS","WeakEnhancer","Black"),"H1hESC",sep="_"),

paste(c("Transcribed","Enhancer", "Repressed", "CTCF", "Promoter","TSS","WeakEnhancer","Black"),"Helas3",sep="_"),

paste(c("Transcribed","Enhancer", "Repressed", "CTCF", "Promoter","TSS","WeakEnhancer","Black"),"Hepg2",sep="_"),

paste(c("Transcribed","Enhancer", "Repressed", "CTCF", "Promoter","TSS","WeakEnhancer","Black"),"Huvec",sep="_"),

paste(c("Transcribed","Enhancer", "Repressed", "CTCF", "Promoter","TSS","WeakEnhancer","Black"),"K562",sep="_") )

chromatinsegmentation_copy=c(

paste(c("1Transcribed","2Enhancer", "3Repressed", "4CTCF", "5Promoter","6TSS","7WeakEnhancer", "8Black"),"Gm12878",sep="_"),

paste(c("1Transcribed","2Enhancer", "3Repressed", "4CTCF", "5Promoter","6TSS","7WeakEnhancer","8Black"),"H1hESC",sep="_"),

paste(c("1Transcribed","2Enhancer", "3Repressed", "4CTCF", "5Promoter","6TSS","7WeakEnhancer","8Black"),"Helas3",sep="_"),

paste(c("1Transcribed","2Enhancer", "3Repressed", "4CTCF", "5Promoter","6TSS","7WeakEnhancer", "8Black"),"Hepg2",sep="_"),

paste(c("1Transcribed","2Enhancer", "3Repressed", "4CTCF", "5Promoter","6TSS","7WeakEnhancer","8Black"),"Huvec",sep="_"),

paste(c("1Transcribed","2Enhancer", "3Repressed", "4CTCF", "5Promoter","6TSS","7WeakEnhancer","8Black"),"K562",sep="_"))


pdf("colors.pdf",width = 8, height = 8)
par(oma=c(0, 0, 0, 5)) #c(bottom, left, top, right)’
heatmap.2(bindedM[order(chromatinsegmentation_copy),] , cexCol = 0.75, labRow= chromatinsegmentation[order(chromatinsegmentation_copy)],
labCol= chromo_names,cexRow = 0.5, cellnote= MT, notecex=0.3, vline=NULL, hline=NULL,tracecol=NULL ,key=FALSE,
notecol="yellow",col=redgreen, breaks=c(-1,1,3,5),Rowv=FALSE,Colv= FALSE,lhei = c(0.05,0.95),lwid = c(0.05,0.95))
dev.off()
