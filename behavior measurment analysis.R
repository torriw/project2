#####################conners_lp_se data##########################

#reading in neuro file containing all behaviors
neuro<-read.csv("04_neuro_all.csv", header= TRUE)

#adding a sample filter of samples used in LD analysis
samples<-gsub("[[:blank:]]+", ",",samples)
samples<-c('CCA209','CCA236','CCA237','CCA238','CCA239','CCA240','CCA242','CCA246','CCA250','CCA255','CCA256','CCA258','CCA259','CCA262','CCA265','CCA267','CCA268','CCA270','CCA271','CCA275','CCA277','CCA302','CCA303','CCA306','CCA308','CCA309','CCA310','CCA311','CCA312','CCA313','CCA317','CCA319','CCA320','CCA321','CCA323','CCA324','CCA325','CCA328','CCA332','CCA333','JWM1044','JWM1054','JWM3324','JWM3442','JWM3583','JWM6201','JWM6204','JWM6206','JWM6208','JWM6209','JWM6234','JWM6236','JWM6238','JWM6240','JWM6241','JWM6242','JWM6243','JWM6248','JWM6249','JWM6251','JWM6255','JWM6257','JWM6266','JWM6267','JWM6270','JWM6285','JWM6286','JWM6301','JWM6314','JWM6332','JWM6337','JWM6343','JWM6344','JWM6346','JWM6348','JWM6350','JWM6354','JWM6375','JWM6380','JWM6400','JWM6405','JWM6410','JWM6414','JWM6415','JWM6431','JWM6494','JWM6495','JWM6500','JWM6503','JWM6513','JWM6525','JWM6527','JWM6534','JWM6535','JWM6547','JWM6553','JWM6556','JWM6566','JWM6578','JWM6579','JWM6621','JWM6649','JWM6652','JWM6665','JWM6670','JWM6676','JWM6679','JWM6681','JWM6682','JWM6688','JWM6691','JWM6718','JWM6719','JWM6790','SMS115','SMS119','SMS123','SMS128','SMS137','SMS139','SMS142','SMS146','SMS154','SMS155','SMS158','SMS168','SMS170','SMS176','SMS178','SMS179','SMS186','SMS203','SMS212','SMS221','SMS227','SMS229','SMS238','SMS239','SMS240','SMS241','SMS244','SMS247','SMS265','SMS277','SMS279','SMS282','SMS283','SMS290','SMS291','SMS292','SMS293','SMS294','SMS295','SMS298','SMS308','SMS310','SMS311','SMS317','SMS322','SMS323','SMS325','SMS326','SMS328','SMS329','SMS332','SMS334','SMS335','SMS343','SMS346','SMS348')

sample_subset<-subset(neuro,IID %in% samples ,select = c(IID,conners_lp_se))

#reading in genotypes that have been converted to # values & creating a subset for matching samples
genotypes<-read.csv("converted genotypes.csv")
genotype_subset<-subset(genotypes,ID %in% sample_subset$IID)


#adding the snp genotypes to the conners behavior measurement 
sample_subset$rs2771040<-(genotype_subset$rs2771040)
sample_subset$rs3199966<-(genotype_subset$rs3199966)
sample_subset$rs10991629<-(genotype_subset$rs10991629)

#omitting all values with na
genotype_subset_na_omit<-na.omit(sample_subset)

#scatter plots for 3 key snps (actual points/vs predicted points)
library(ggplot2)
rs277_plot_point<-ggplot(genotype_subset_na_omit,aes(x=rs2771040,y=conners_lp_se))+geom_point()+geom_smooth()
rs277_plot_jitter<-ggplot(genotype_subset_na_omit,aes(x=rs2771040,y=conners_lp_se))+geom_jitter()+geom_smooth()

rs319_plot_point<-ggplot(genotype_subset_na_omit,aes(x=rs3199966,y=conners_lp_se))+geom_point()+geom_smooth()
rs319_plot_jitter<-ggplot(genotype_subset_na_omit,aes(x=rs3199966,y=conners_lp_se))+geom_jitter()+geom_smooth()

rs109_plot_point<-ggplot(genotype_subset_na_omit,aes(x=rs10991629,y=conners_lp_se))+geom_point()+geom_smooth()
rs109_plot_jitter<-ggplot(genotype_subset_na_omit,aes(x=rs10991629,y=conners_lp_se))+geom_jitter()+geom_smooth()

##################dasii#######################

genotypes_switched<-read.csv("converted genotypes.csv")
genotype_subset<-subset(genotypes_switched,ID %in% sample_subset$IID)

dasii_sample_subset<-subset(neuro,IID %in% c('CCA209','CCA236','CCA237','CCA238','CCA239','CCA240','CCA242','CCA246','CCA250','CCA255','CCA256','CCA258','CCA259','CCA262','CCA265','CCA267','CCA268','CCA270','CCA271','CCA275','CCA277','CCA302','CCA303','CCA306','CCA308','CCA309','CCA310','CCA311','CCA312','CCA313','CCA317','CCA319','CCA320','CCA321','CCA323','CCA324','CCA325','CCA328','CCA332','CCA333','JWM1044','JWM1054','JWM3324','JWM3442','JWM3583','JWM6201','JWM6204','JWM6206','JWM6208','JWM6209','JWM6234','JWM6236','JWM6238','JWM6240','JWM6241','JWM6242','JWM6243','JWM6248','JWM6249','JWM6251','JWM6255','JWM6257','JWM6266','JWM6267','JWM6270','JWM6285','JWM6286','JWM6301','JWM6314','JWM6332','JWM6337','JWM6343','JWM6344','JWM6346','JWM6348','JWM6350','JWM6354','JWM6375','JWM6380','JWM6400','JWM6405','JWM6410','JWM6414','JWM6415','JWM6431','JWM6494','JWM6495','JWM6500','JWM6503','JWM6513','JWM6525','JWM6527','JWM6534','JWM6535','JWM6547','JWM6553','JWM6556','JWM6566','JWM6578','JWM6579','JWM6621','JWM6649','JWM6652','JWM6665','JWM6670','JWM6676','JWM6679','JWM6681','JWM6682','JWM6688','JWM6691','JWM6718','JWM6719','JWM6790','SMS115','SMS119','SMS123','SMS128','SMS137','SMS139','SMS142','SMS146','SMS154','SMS155','SMS158','SMS168','SMS170','SMS176','SMS178','SMS179','SMS186','SMS203','SMS212','SMS221','SMS227','SMS229','SMS238','SMS239','SMS240','SMS241','SMS244','SMS247','SMS265','SMS277','SMS279','SMS282','SMS283','SMS290','SMS291','SMS292','SMS293','SMS294','SMS295','SMS298','SMS308','SMS310','SMS311','SMS317','SMS322','SMS323','SMS325','SMS326','SMS328','SMS329','SMS332','SMS334','SMS335','SMS343','SMS346','SMS348')
,select = c(IID,dasiiseqquanreasperc))

dasii_genotype_subset<-subset(genotypes_switched,ID %in% dasii_sample_subset$IID)

dasii_sample_subset$rs10991629<-(dasii_genotype_subset$rs10991629)
dasii_sample_subset$rs2771040<-(dasii_genotype_subset$rs2771040)
dasii_sample_subset$rs3199966<-(dasii_genotype_subset$rs3199966)
dasii_sample_subset$rs75106836<-(dasii_genotype_subset$rs75106836)
dasii_sample_subset$rs143438338<-(dasii_genotype_subset$rs143438338)
dasii_sample_subset$rs59370172<-(dasii_genotype_subset$rs59370172)



dasii_subset_na_omit<-na.omit(dasii_sample_subset)

#exporting and converting file to csv
write.csv(dasii_subset_na_omit,"C:/Users/torriw/Desktop/dasii_subset_na_omit.csv",row.names = TRUE)

#rs277
png("dasii_rs277_plot_point.png", units = "in", width = 7, height = 7, res = 600)
ggplot(dasii_subset_na_omit,aes(x=rs2771040,y=dasiiseqquanreasperc))+geom_point()+geom_smooth()+ggtitle("Behavioral Measurement \n by Genotype")+theme(plot.title = element_text(hjust = 0.5))+xlab("Genotype (rs2771040)")
dev.off()


png("dasii_rs277_plot_jitter.png", units = "in", width = 7, height = 7, res = 600)
ggplot(dasii_subset_na_omit,aes(x=rs2771040,y=dasiiseqquanreasperc))+geom_jitter(width = 0.1, height = 0.1)+geom_smooth()+ggtitle("Behavioral Measurement \n by Genotype")+theme(plot.title = element_text(hjust = 0.5))+xlab("Genotype (rs2771040)")
dev.off()

#rs319
png("dasii_rs319_plot_point.png", units = "in", width = 7, height = 7, res = 600)
ggplot(dasii_subset_na_omit,aes(x=rs3199966,y=dasiiseqquanreasperc))+geom_point()+geom_smooth()+ggtitle("Behavioral Measurement \n by Genotype")+theme(plot.title = element_text(hjust = 0.5))+xlab("Genotype (rs3199966)")
dev.off()

png("dasii_rs319_plot_jitter.png", units = "in", width = 7, height = 7, res = 600)
ggplot(dasii_subset_na_omit,aes(x=rs3199966,y=dasiiseqquanreasperc))+geom_jitter(width = 0.1, height = 0.1)+geom_smooth()+ggtitle("Behavioral Measurement \n by Genotype")+theme(plot.title = element_text(hjust = 0.5))+xlab("Genotype (rs3199966)")
dev.off()

#rs109
png("dasii_rs109_plot_point.png", units = "in", width = 7, height = 7, res = 600)
ggplot(dasii_subset_na_omit,aes(x=rs10991629,y=dasiiseqquanreasperc))+geom_point()+geom_smooth()+ggtitle("Behavioral Measurement \n by Genotype")+theme(plot.title = element_text(hjust = 0.5))+xlab("Genotype (rs10991629)")
dev.off()

png("dasii_rs109_plot_jitter.png", units = "in", width = 7, height = 7, res = 600)
ggplot(dasii_subset_na_omit,aes(x=rs10991629,y=dasiiseqquanreasperc))+geom_jitter(width = 0.1, height = 0.1)+geom_smooth()+ggtitle("Behavioral Measurement \n by Genotype")+theme(plot.title = element_text(hjust = 0.5))+xlab("Genotype (rs10991629)")
dev.off()

#rs751
png("dasii_rs751_plot_point.png", units = "in", width = 7, height = 7, res = 600)
ggplot(dasii_subset_na_omit,aes(x=rs75106836,y=dasiiseqquanreasperc))+geom_point()+geom_smooth()+ggtitle("Behavioral Measurement \n by Genotype")+theme(plot.title = element_text(hjust = 0.5))+xlab("Genotype (rs75106836)")
dev.off()

png("dasii_rs751_plot_jitter.png", units = "in", width = 7, height = 7, res = 600)
ggplot(dasii_subset_na_omit,aes(x=rs75106836,y=dasiiseqquanreasperc))+geom_jitter(width = 0.1, height = 0.1)+geom_smooth()+ggtitle("Behavioral Measurement \n by Genotype")+theme(plot.title = element_text(hjust = 0.5))+xlab("Genotype (rs75106836)")
dev.off()

#rs143
png("dasii_rs143_plot_point.png", units = "in", width = 7, height = 7, res = 600)
ggplot(dasii_subset_na_omit,aes(x=rs143438338,y=dasiiseqquanreasperc))+geom_point()+geom_smooth()+ggtitle("Behavioral Measurement \n by Genotype")+theme(plot.title = element_text(hjust = 0.5))+xlab("Genotype (rs143438338)")
dev.off()

png("dasii_rs143_plot_jitter.png", units = "in", width = 7, height = 7, res = 600)
ggplot(dasii_subset_na_omit,aes(x=rs143438338,y=dasiiseqquanreasperc))+geom_jitter(width = 0.1, height = 0.1)+geom_smooth()+ggtitle("Behavioral Measurement \n by Genotype")+theme(plot.title = element_text(hjust = 0.5))+xlab("Genotype (rs143438338)")
dev.off()

#rs593
png("dasii_rs593_plot_point.png", units = "in", width = 7, height = 7, res = 600)
ggplot(dasii_subset_na_omit,aes(x=rs59370172,y=dasiiseqquanreasperc))+geom_point()+geom_smooth()+ggtitle("Behavioral Measurement \n by Genotype")+theme(plot.title = element_text(hjust = 0.5))+xlab("Genotype (rs59370172)")
dev.off()

png("dasii_rs593_plot_jitter.png", units = "in", width = 7, height = 7, res = 600)
ggplot(dasii_subset_na_omit,aes(x=rs59370172,y=dasiiseqquanreasperc))+geom_jitter(width = 0.1, height = 0.1)+geom_smooth()+ggtitle("Behavioral Measurement \n by Genotype")+theme(plot.title = element_text(hjust = 0.5))+xlab("Genotype (rs59370172)")
dev.off()

#regression line values
a<-dasii_subset_na_omit$rs10991629
b<-dasii_subset_na_omit$dasiiseqquanreasperc
rs109_reg = lm(a~b)
rs109_summary<-summary(rs109_reg)

c<-dasii_subset_na_omit$rs2771040
d<-dasii_subset_na_omit$dasiiseqquanreasperc
rs277_reg = lm(c~d)
rs277_summary<-summary(rs277_reg)

e<-dasii_subset_na_omit$rs3199966
f<-dasii_subset_na_omit$dasiiseqquanreasperc
rs319_reg = lm(e~f)
rs319_summary<-summary(rs319_reg)

g<-dasii_subset_na_omit$rs143438338
h<-dasii_subset_na_omit$dasiiseqquanreasperc
rs143_reg = lm(g~h)
rs143_summary<-summary(rs143_reg)

i<-dasii_subset_na_omit$rs75106836
j<-dasii_subset_na_omit$dasiiseqquanreasperc
rs751_reg = lm(i~j)
rs751_summary<-summary(rs751_reg)

k<-dasii_subset_na_omit$rs59370172
l<-dasii_subset_na_omit$dasiiseqquanreasperc
rs593_reg = lm(k~l)
rs593_summary<-summary(rs593_reg)

#finding ethnic proportions
samples_ethnicity<-read.csv("Ethnic Proportion by Site sample_260.csv", header = TRUE)
samples_ethnicity_subset<-subset(samples_ethnicity,IID.1 %in% dasii_subset_na_omit$IID)

#writing out file
write.csv(samples_ethnicity_subset,"C:/Users/torriw/Desktop/samples_ethnicity_subset.csv",row.names = FALSE)

#finding FAS vs Control 
expvscon<-read.csv("CIFASD_3_demgroupclass.csv", header = TRUE)
expvscon_subset<-subset(expvscon,subjectid %in% dasii_subset_na_omit$IID)

write.csv(expvscon_subset,"C:/Users/torriw/Desktop/expvscon_subset.csv",row.names = FALSE )
