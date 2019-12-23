#source("dge_plot_funcs.r")
#install.package(xx,lib="~/R/lib")
#save.image(file="*.RData")
#load(file="")

library("caret")

trainr=read.table("./1.combined.feature.txt",header = T,row.names = 1,stringsAsFactors = FALSE)

#######################select a subset of data
train=subset(trainr,status!=4)
train=subset(train,status!=2)
train$frqCG=train$frqCGT+train$frqCGA+train$frqCGC+train$frqCGG
features <- c("status", "genelength", "transnum","expression",
                "DnDs", "Dis_Cen", "frqCWG", "frqCG")

train=train[,features]
str(train)
#############################impute data
pre.pro=preProcess(train[,-1],method="knnImpute")
train.imp=predict(pre.pro,train[,-1])
train$DnDs=train.imp$DnDs
train$Dis_Cen=train.imp$Dis_Cen

############SCALE the data
pr.train=preProcess(train[,-1],method=c('center','scale'))
f.train=predict(pr.train,train[,-1])
train=cbind(train$status,f.train)
names(train)=c("status","genelen","transnum","expr","dnds",
                  "discen","fCWG","fCG")

##############refine the name of status
levels=unique(c(train[,1]))
train[,1]=factor(train[,1],labels=make.names(levels))

set.seed(10)
indexes <- createDataPartition(train$status,
                               times = 1,
                               p = 0.5,
                               list = FALSE)
g.train <- train[indexes,]
g.test <- train[-indexes,]

prop.table(table(train$status))
prop.table(table(g.train$status))
prop.table(table(g.test$status))

###########################################model training

fitControl <- trainControl(
  method = 'repeatedcv',                   # k-fold cross validation
  number = 10,  
  repeats = 10,
  search="grid" 
) 

######################################test different models
#modelnames <- paste(names(getModelInfo()), collapse=',  ')
#modelnames
#write.table(modelnames,file = "modelnames.txt")

##1random forest good
model_rf=train(status ~. ,data=g.train, method='rf',trControl=fitControl,tuneLength=5)
varimp_rf=varImp(model_rf)
#plot(varimp_rf,xlab="importance of features",cex=2)
##2svm good slow
#model_svm=train(status ~. ,data=g.train, method='svmRadial',trControl=fitControl,tuneLength=5)
#varimp_svm=varImp(model_svm)


##3glm so so
#model_glm=train(status ~. ,data=g.train, method='glm',trControl=fitControl,tuneLength=5)
#varimp_glm=varImp(model_glm)


##4 Bagged cart ok
#model_bc=train(status ~. ,data=g.train, method='treebag',trControl=fitControl,tuneLength=5)
#varimp_bc=varImp(model_bc)


##5 generalized booseted model good
#model_gbm=train(status ~. ,data=g.train, method='gbm',trControl=fitControl,tuneLength=5)

## 6 earth good
#model_earth=train(status ~. ,data=g.train, method='earth',trControl=fitControl,tuneLength=5)
#varimp_earth=varImp(model_earth)


## 7 neural network good
#model_nnet=train(status ~. ,data=g.train, method='nnet',trControl=fitControl,tuneLength=5)
#varimp_nnet=varImp(model_nnet)


## 8 averageed nnet good but slow
#model_annet=train(status ~. ,data=g.train, method='avNNet',trControl=fitControl,tuneLength=5)
#varimp_annet=varImp(model_annet)


## 9 logic regression so so

#model_lgrg=train(status ~. ,data=g.train, method='regLogistic',trControl=fitControl,tuneLength=5)
#varimp_lgrg=varImp(model_lgrg)


## 10 J48 so so
#library("RWeka")
#model_J48=train(status ~. ,data=g.train, method='J48',trControl=fitControl,tuneLength=5)
#varimp_J48=varImp(model_J48)


## 11 cart so so
#library("rpart")
#model_cart=train(status ~. ,data=g.train, method='rpart',trControl=fitControl,tuneLength=5)
#varimp_cart=varImp(model_cart)


## 12 knn ok
#model_knn=train(status ~. ,data=g.train, method='kknn',trControl=fitControl,tuneLength=5)
#varimp_knn=varImp(model_knn)


## 13 Logistic model trees good
#model_LMT=train(status ~. ,data=g.train, method='LMT',trControl=fitControl,tuneLength=5)
#varimp_LMT=varImp(model_LMT)

save.image(file="./modelselect.RData")

## 14 tree model from genetic algorithms 

#model_evtree=train(status ~. ,data=g.train, method='evtree',trControl=fitControl,tuneLength=5)
#varimp_evtree=varImp(model_evtree)


## 15 extreme gradient boosting good very slow
#model_xgb=train(status ~. ,data=g.train, method='xgbTree',trControl=fitControl,tuneLength=5)
#varimp_xgb=varImp(model_xgb)

#save.image(file="./modelselect.RData")
