setwd("C:/Users/zywlm/Documents/02_my project/02_gbm_CMT3/version3.0/06_gbM_featu/2_work/11.ML/")

library(lubridate)
library(caret)
library(caretEnsemble)

#########################################
load(file="./modelselect.RData")

########################################prediction perfermance by test data

pred2=predict(model_rf,g.test,type = "raw")
confusionMatrix(pred2, g.test$status)

pdf(
  "2.feature.usefullness.pdf",
  width = 3,
  height = 3,
)
varimp_rf=varImp(model_rf)
#varimp_rf$importance$Overall=c(100.00,0.45,35.10,47.65,35.89,54.28,51.14)
plot(varimp_rf,xlab="importance of features",cex=1,cex.axis=5)

dev.off()

pred=predict(model_rf,g.test,type = "prob")
predr=as.data.frame(predict(model_rf,g.test,type = "raw"))
row.names(predr)=row.names(pred)
colnames(predr)="predict"

######################################correlation between predict gbM probablity and gbM conservation?
tr=as.data.frame(trainr$gbM.freq)
colnames(tr)="gbM.freq";
row.names(tr)=row.names(trainr)
res <- Reduce(function(a,b){
  ans <- merge(a,b,by="row.names",all=T)
  row.names(ans) <- ans[,"Row.names"]
  ans[,!names(ans) %in% "Row.names"]
}, list(g.test,pred,predr,tr))


m=na.omit(res)###combined gbmfreq and predict probability together for test dataset with original gbM status

write.table(m,file="testingresult.table",sep="\t",quote=F)
m2=m[m[,1]=='X1',]
m9=m2[m2$gbM.freq>0.9,]
tp=length(m9[m9$X1>=0.5,1])
fn=length(m9[m9$X1<0.5,1])


#################from 0.01-0.1
lim=1
m0=m2[m2$gbM.freq<=lim,]
hist(m0$gbM.freq)
l=length(m0[,1])
d=20
step=(lim)/d
up=c()##record the ranges 
down=c()
start=0
ratio=c()
m3=matrix()
m4=matrix()
freq=c()
for (i in 1:d) {
  up[i]=start
  down[i]=up[i]+step
  start=down[i]
  m3=m0[m0$gbM.freq<=down[i],]
  m4=m3[m3$gbM.freq>up[i],]
  tp=length(m4[m4$X1>=0.5,1])
  r=tp/length(m4[,1])
  ratio[i]=r
  freq[i]=length(m4[,1])
}
c=data.frame(up=up,down=down,ratio=ratio,freq=freq)
plot(c[,1],c[,3])
plot(c[,1],c[,4])
write.table(c,file="lowgbMfreq.table",sep="\t",quote=F)
####gglot2 dr
# png(
#   "2.low.gbMgene.prediction.png",
#   width = 9,
#   height = 6,
#   units = "in",
#   bg = "white",
#   res = 600
# )

pdf(
  "2.low.gbMgene.prediction.v2.0.pdf",
  width = 6,
  height = 4
)
library(ggplot2)
ggplot(c,aes(x=up*100,y=ratio,size=freq,col=freq))+
  geom_point()+
  scale_color_gradient(name="Number of genes",low = "black",high = "#1D70B7")+
  scale_size(guide = F)+
  theme_classic(base_size=15) +
  ylab("Sensitivity")+
  xlab("gbM epiallele frequency (%)")+
  xlim(c(0,100))+
  ylim(c(0.55,1))
  
dev.off()

################################correlation between gbM freq and classification accurary

#xrg <- ggplot_build(p)$layout$panel_ranges[[1]]$x.range
# pdf(
#   "1.judgement.pdf",
#   width = 9,
#   height = 6,
# )

pdf(
  "1.judgement.pdf",
  width = 9,
  height = 6,
)
cbbPalette <- c("#b30000", "#56b4e9" )
p=ggplot(m2,aes(x=X1,y=gbM.freq))+
  geom_point(aes(col=predict),pch = 21, alpha=0.5,size=1) +
  stat_density_2d(n=200,alpha=0.5,contour = TRUE)+
  scale_color_manual(name="Sensitivity",labels=c("TP","FN"),values = cbbPalette) +
  theme_classic(base_size=25) +
  xlab("predict probability of gbM")+
  ylab("gbM epiallele frequency") 

p + annotate("text", x = 0.5, y = 1.07, label = "Prediction by random forest", size = 5) +
  annotate("text", x = 0.12, y = 1, label = "Acc: 0.79", size = 5,col='black',alpha=0.6) +
  annotate("text", x = 0.12, y = 0.95, label = "Sen: 0.77", size = 5,col='red',alpha=0.6) +
  annotate("text", x = 0.12, y = 0.90, label = "Spe: 0.80", size = 5,col='black',alpha=0.6)

dev.off()

##############################plot model accuracy rankings

t=read.table("./data/2.model.accu.eva.txt",header = T,stringsAsFactors = FALSE)

# t$models=factor(t$models,levels=c("rf","gbm","svmRadial","avNNet","treebag","earth",
#                                     "LMT","J48","kknn","rpart","regLogistic","glm"))

t$models=factor(t$models,levels=c("glm","regLogistic","rpart","kknn","J48","LMT",
                                  "earth","treebag","avNNet","svmRadial","gbm","rf"))

# cbbpale<- c("#CC79A7", "#0072B2", "#D55E00", "#009E73","#E69F00", "#F0E442", "#56B4E9",
#                 "#d556e9","#5f1cd7","#94d71c","#1cd7bd","#d71c36")
cbbpale2<- c("#CC79A7", "#CC79A7", "#0072B2", "#D55E00","#0072B2", "#009E73", "#CC79A7",
            "#009E73","#D55E00","#D55E00","#009E73","#009E73")
pdf(
   "2.modelrank.pdf",
   width = 8,
   height = 6,
)
p=ggplot(t,aes(x=models,y=Mean,ymin=X1st,ymax=X3rd,col=models)) +
  geom_pointrange(size=2,alpha=0.8,shape=20) +
  scale_color_manual(values = cbbpale2,guide=F ) +
  theme_classic(base_size=25) +
  ylab("Accuracy")+
  xlab(" ") +
  coord_flip()
p

dev.off()

############################plot a example of 
