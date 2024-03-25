

library(PerFit)
library(pROC)
library(ggplot2)
library(randomForest)
library(pracma)

#----------------------------------------
#save each condition with 100 iterations

dcv = NULL
for (cs in c(1, 2, 3)){
  for (cl in c(0.1, 0.3, 0.5)){
    for (i in c(10, 20, 30)) {
      for (n in c(500, 750, 1000)) {
        for (r in 1:100) {
          
          dataorg = read.table(paste("sdb-",cs,"-",cl,"-",i,"-",n,"-",r,".dat",sep=""), header=FALSE)
          temp=data.frame(dataorg[,1:i],dataorg[,(2*i+1)])
          dcv=rbind(dcv,temp)

        }
        write.table(dcv, file = paste("cv-",cs,"-",cl,"-",i,"-",n,".csv",sep=""), row.names=FALSE, na="", col.names=T, sep=",")
        dcv = NULL
        }}}}

#------------------
# cross validation

start = proc.time()
d.cv = NULL
for (cs in c(1, 2, 3)){
  for (cl in c(0.1, 0.3, 0.5)){
    for (i in c(10, 20, 30)) {
      for (n in c(500, 750, 1000)) {
        
        dataorg = read.table(paste("cv-",cs,"-",cl,"-",i,"-",n,".csv",sep=""), header=TRUE, sep = ",")
        colnames(dataorg)[i+1] = "class"
        #shuffle the data
        d.shuff = dataorg[sample(1:nrow(dataorg)), ]
        
        #create 5 folds
        nrFolds <- 5
        
        # generate array containing fold-number for each sample (row)
        folds <- rep_len(1:nrFolds, nrow(d.shuff))
        
        # actual cross validation
        for(k in 1:nrFolds) {
          # actual split of the data
          fold <- which(folds == k)
          train <- d.shuff[-fold,]
          test <- d.shuff[fold,]
          
          ############################################################# train
          #14 person fit indices
          rpb.tr <- r.pbis(train[,1:i], IRT.PModel = "2PL")
          csato.tr <- C.Sato(train[,1:i], IRT.PModel = "2PL")
          cstar.tr <- Cstar(train[,1:i], IRT.PModel = "2PL")
          g.tr <- G(train[,1:i], IRT.PModel = "2PL")
          gn.tr <- Gnormed(train[,1:i], IRT.PModel = "2PL")
          akb.tr <- A.KB(train[,1:i], IRT.PModel = "2PL")
          dkb.tr <- D.KB(train[,1:i], IRT.PModel = "2PL")
          ekb.tr <- E.KB(train[,1:i], IRT.PModel = "2PL")
          u3.tr <- U3(train[,1:i], IRT.PModel = "2PL")
          zu3.tr <- ZU3(train[,1:i], IRT.PModel = "2PL")
          nci.tr <- NCI(train[,1:i], IRT.PModel = "2PL")
          ht.tr <- Ht(train[,1:i], IRT.PModel = "2PL")
          lz.tr <- lz(train[,1:i], IRT.PModel = "2PL")
          lzstar.tr <- lzstar(train[,1:i], IRT.PModel = "2PL")
          
          temp.tr = data.frame(rpb.tr[["PFscores"]], csato.tr[["PFscores"]], cstar.tr[["PFscores"]], g.tr[["PFscores"]], gn.tr[["PFscores"]], 
                               akb.tr[["PFscores"]], dkb.tr[["PFscores"]], ekb.tr[["PFscores"]], u3.tr[["PFscores"]], zu3.tr[["PFscores"]], 
                               nci.tr[["PFscores"]], ht.tr[["PFscores"]], lz.tr[["PFscores"]], lzstar.tr[["PFscores"]],train[,i+1])
          
          ################################################################### test
          
          #14 person fit indices
          rpb.ts <- r.pbis(test[,1:i], IRT.PModel = "2PL")
          csato.ts <- C.Sato(test[,1:i], IRT.PModel = "2PL")
          cstar.ts <- Cstar(test[,1:i], IRT.PModel = "2PL")
          g.ts <- G(test[,1:i], IRT.PModel = "2PL")
          gn.ts <- Gnormed(test[,1:i], IRT.PModel = "2PL")
          akb.ts <- A.KB(test[,1:i], IRT.PModel = "2PL")
          dkb.ts <- D.KB(test[,1:i], IRT.PModel = "2PL")
          ekb.ts <- E.KB(test[,1:i], IRT.PModel = "2PL")
          u3.ts <- U3(test[,1:i], IRT.PModel = "2PL")
          zu3.ts <- ZU3(test[,1:i], IRT.PModel = "2PL")
          nci.ts <- NCI(test[,1:i], IRT.PModel = "2PL")
          ht.ts <- Ht(test[,1:i], IRT.PModel = "2PL")
          lz.ts <- lz(test[,1:i], IRT.PModel = "2PL")
          lzstar.ts <- lzstar(test[,1:i], IRT.PModel = "2PL")
          
          temp.ts = data.frame(rpb.ts[["PFscores"]], csato.ts[["PFscores"]], cstar.ts[["PFscores"]], g.ts[["PFscores"]], gn.ts[["PFscores"]], 
                               akb.ts[["PFscores"]], dkb.ts[["PFscores"]], ekb.ts[["PFscores"]], u3.ts[["PFscores"]], zu3.ts[["PFscores"]], 
                               nci.ts[["PFscores"]], ht.ts[["PFscores"]], lz.ts[["PFscores"]], lzstar.ts[["PFscores"]],test[,i+1])
          
          
          #calculate AUC
          temp3 = NULL
          for (p in 1:14) {
            
            d.train = data.frame(temp.tr[,p], temp.tr[,15])
            d.train <- na.omit(d.train)
            colnames(d.train)=c("pfscr","class")
            
            # c1=biased, c2=unbiased
            d.train$class[d.train$class==1] <- 1
            d.train$class[d.train$class==2] <- 0
            glm.fit.tr = glm(class ~ pfscr, family=binomial, data=d.train)
            
            ## ROC by random forest
            rf.model.tr <- randomForest(factor(class) ~ pfscr, data = d.train, keep.forest=TRUE)
            #roc.rf.tr = roc(d.train$class, rf.model.tr$votes[,2])
            
            #----------------------------------------------------------------------------
            
            d.test = data.frame(temp.ts[,p], temp.ts[,15])
            d.test <- na.omit(d.test)
            colnames(d.test)=c("pfscr","class")
            
            #ROC by logistic regression
            d.test$class[d.test$class==1] <- 1
            d.test$class[d.test$class==2] <- 0
            glm.pred = predict(glm.fit.tr, newdata = d.test, family=binomial)
            roc.glm.ts = roc(d.test$class, glm.pred)
            
            ## ROC by random forest
            rf.pred = predict(rf.model.tr, newdata = d.test, type="prob")
            roc.rf.ts = roc(d.test$class, rf.pred[,2])
            
            temp2 = cbind(roc.glm.ts[["auc"]], roc.rf.ts[["auc"]])
            temp3 = cbind(temp3, temp2)
          } #p
          
          temp3 = data.frame(cs,cl,i,n,k,temp3)
          d.cv = rbind(d.cv,temp3)
          
          print(paste("k=",k))
        } #k
        
        print(paste("cs=",cs,", cl=",cl,", i=",i,", n=",n))
          
      } #n
    } #i
  } #cl
} #cs

write.table(d.cv, file = "d.cv.csv", row.names=FALSE, na="", col.names=T, sep=",")

end = proc.time()
end - start

#-----------------------------------
#means of 14 pf AUCs across 5 folds

data = read.csv("d.cv.csv", header=TRUE)

count = seq(1,405, by=5)

temp3 = NULL

for (i in count) {
  temp2 = NULL
  for (j in 6:33) {
    temp <- mean(data[i:(i+4), j])
    temp2 = cbind(temp2,temp)
  }
  temp3=rbind(temp3, temp2) 
}

temp3 = as.data.frame(temp3)

colnames(temp3)[1:28]=c("lr.rpb","rf.rpb", "lr.csato","rf.csato", "lr.cstar","rf.cstar", "lr.g","rf.g", "lr.gn","rf.gn",
                        "lr.akb","rf.akb", "lr.dkb","rf.dkb", "lr.ekb","rf.ekb", "lr.u3","rf.u3", "lr.zu3","rf.zu3",
                        "lr.nci","rf.nci", "lr.ht","rf.ht", "lr.lz","rf.lz", "lr.lzstar","rf.lzstar")   

write.csv(temp3, "cv means.csv")

#--------------------------------
#sd of 14 pf AUCs across 5 folds

data = read.csv("d.cv.csv", header=TRUE)

count = seq(1,405, by=5)

temp3 = NULL

for (i in count) {
  temp2 = NULL
  for (j in 6:33) {
    temp <- sd(data[i:(i+4), j])
    temp2 = cbind(temp2,temp)
  }
  temp3=rbind(temp3, temp2) 
}

temp3 = as.data.frame(temp3)

colnames(temp3)[1:28]=c("lr.rpb","rf.rpb", "lr.csato","rf.csato", "lr.cstar","rf.cstar", "lr.g","rf.g", "lr.gn","rf.gn",
                        "lr.akb","rf.akb", "lr.dkb","rf.dkb", "lr.ekb","rf.ekb", "lr.u3","rf.u3", "lr.zu3","rf.zu3",
                        "lr.nci","rf.nci", "lr.ht","rf.ht", "lr.lz","rf.lz", "lr.lzstar","rf.lzstar")   

write.csv(temp3, "cv sd.csv")

#---------
#sd range

data = read.csv("cv sd lr.csv", header=TRUE)
range(data[,2:15])

data = read.csv("cv sd rf.csv", header=TRUE)
range(data[,2:15])

