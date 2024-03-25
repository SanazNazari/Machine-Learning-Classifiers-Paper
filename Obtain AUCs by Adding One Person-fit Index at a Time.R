

library(PerFit)
library(pROC)
library(ggplot2)
library(randomForest)
library(pracma)

#-------------------------------------------------------------
#compute 81 pf-score datasets for 81 conditions as features

start = proc.time()

for (cs in c(1, 2, 3)){
  for (cl in c(0.1, 0.3, 0.5)){
    for (i in c(10, 20, 30)) {
      for (n in c(500, 750, 1000)) {
          
          dataorg = read.table(paste("cv-",cs,"-",cl,"-",i,"-",n,".csv",sep=""), header=TRUE, sep = ",")
          colnames(dataorg)[i+1] = "class"
          data=dataorg[,1:i]
          
          #compute pf score
          rpb.out <- r.pbis(data, IRT.PModel = "2PL")
          csato.out <- C.Sato(data, IRT.PModel = "2PL")
          cstar.out <- Cstar(data, IRT.PModel = "2PL")
          g.out <- G(data, IRT.PModel = "2PL")
          gn.out <- Gnormed(data, IRT.PModel = "2PL")
          akb.out <- A.KB(data, IRT.PModel = "2PL")
          dkb.out <- D.KB(data, IRT.PModel = "2PL")
          ekb.out <- E.KB(data, IRT.PModel = "2PL")
          u3.out <- U3(data, IRT.PModel = "2PL")
          zu3.out <- ZU3(data, IRT.PModel = "2PL")
          nci.out <- NCI(data, IRT.PModel = "2PL")
          ht.out <- Ht(data, IRT.PModel = "2PL")
          lz.out <- lz(data, IRT.PModel = "2PL")
          lzstar.out <- lzstar(data, IRT.PModel = "2PL")
          
          temp = data.frame(rpb.out[["PFscores"]], csato.out[["PFscores"]], cstar.out[["PFscores"]], g.out[["PFscores"]], gn.out[["PFscores"]], 
                            akb.out[["PFscores"]], dkb.out[["PFscores"]], ekb.out[["PFscores"]], u3.out[["PFscores"]], zu3.out[["PFscores"]], 
                            nci.out[["PFscores"]], ht.out[["PFscores"]], lz.out[["PFscores"]], lzstar.out[["PFscores"]],dataorg[,i+1])
          
          temp = na.omit(temp)
          write.table(temp, file = paste("pf-",cs,"-",cl,"-",i,"-",n,".csv",sep=""), row.names=FALSE, na="", col.names=T, sep=",")
          temp = NULL 

      } #n
    } #i
  } #cl
} #cs

end = proc.time()
end - start

#---------------------------------------------------------
#increase the number of pfs and get aucs by random forest

d.means = read.csv("auc means rf.csv", header=TRUE)

start = proc.time()

d.inc = NULL
counter = 0  
for (cs in c(1, 2, 3)){
  for (cl in c(0.1, 0.3, 0.5)){
    for (i in c(10, 20, 30)) {
      for (n in c(500, 750, 1000)) {
        
        counter = counter + 1
        d.temp = d.means[counter,]
        index = order(d.temp)
        
        data = read.table(paste("pf-",cs,"-",cl,"-",i,"-",n,".csv",sep=""), header=TRUE, sep = ",")
        colnames(data)[1:15] = c("rpb", "c", "mc", "g", "ng", "a", "d", "e", "u3", "zu3", "nci", "ht", "lz", "clz", "class")
        
        # c1=biased, c2=unbiased
        data$class[data$class==1] <- 1
        data$class[data$class==2] <- 0

        data.ordered = cbind(data[,index],data$class)
        
        temp3 = NULL #create storage
        for (u in 1:14) {
          
          pfdata = data.frame(data.ordered[,15], data.ordered[,1:u])
          colnames(pfdata)[1]= "class"
          
          ## ROC by random forest
          rf.model <- randomForest(factor(pfdata$class) ~ ., pfdata)
          roc.rf = roc(pfdata$class, rf.model$votes[,2])
          
          temp2 = roc.rf[["auc"]]
          temp3 = cbind(temp3, temp2)
          
        } #u
        
        print(paste("counter=",counter,"cs=",cs,", cl=",cl,", i=",i,", n=",n))
        temp3 = data.frame(cs,cl,i,n,temp3)
        d.inc=rbind(d.inc,temp3)
        
      } #n
    } #i
  } #cl
} #cs

write.table(d.inc, file = "d.rf.increase.csv", row.names=FALSE, na="", col.names=T, sep=",")

end = proc.time()
end - start

#---------------------------------------------------------------
#increase the number of pfs and get aucs by logistic regression

d.means = read.csv("auc means lr.csv", header=TRUE)

start = proc.time()

d.inc = NULL
counter = 0  
for (cs in c(1, 2, 3)){
  for (cl in c(0.1, 0.3, 0.5)){
    for (i in c(10, 20, 30)) {
      for (n in c(500, 750, 1000)) {
        
        counter = counter + 1
        d.temp = d.means[counter,]
        index = order(d.temp)
        
        data = read.table(paste("pf-",cs,"-",cl,"-",i,"-",n,".csv",sep=""), header=TRUE, sep = ",")
        colnames(data)[1:15] = c("rpb", "c", "mc", "g", "ng", "a", "d", "e", "u3", "zu3", "nci", "ht", "lz", "clz", "class")
        
        # c1=biased, c2=unbiased
        data$class[data$class==1] <- 1
        data$class[data$class==2] <- 0
        
        data.ordered = cbind(data[,index],data$class)
        
        temp3 = NULL #create storage
        for (u in 1:14) {
          
          pfdata = data.frame(data.ordered[,15], data.ordered[,1:u])
          colnames(pfdata)[1]= "class"
          
          #ROC by logistic regression
          glm.fit = glm(pfdata$class ~ ., family=binomial,data = pfdata)
          roc.glm = roc(pfdata$class, glm.fit$fitted.values)
          
          temp2 = roc.glm[["auc"]]
          temp3 = cbind(temp3, temp2)
          
        } #u
        
        print(paste("counter=",counter,"cs=",cs,", cl=",cl,", i=",i,", n=",n))
        temp3 = data.frame(cs,cl,i,n,temp3)
        d.inc=rbind(d.inc,temp3)
        
      } #n
    } #i
  } #cl
} #cs

write.table(d.inc, file = "d.lr.increase.csv", row.names=FALSE, na="", col.names=T, sep=",")

end = proc.time()
end - start





