
library(PerFit)
library(pROC)
library(ggplot2)
library(randomForest)
library(pracma)

#--------------------------------------------------------------------
#Create a large dataset including 100 iterations for each condition.

#r = replication
#n = sample size
#i = number of items
#cl = cross loading
#cs = class separation

start = proc.time()

danova = NULL
for (cs in c(1, 2, 3)){
  for (cl in c(0.1, 0.3, 0.5)){
    for (i in c(10, 20, 30)) {
      for (n in c(500, 750, 1000)) {
        for (r in 1:100) {
          
          dataorg = read.table(paste("sdb-",cs,"-",cl,"-",i,"-",n,"-",r,".dat",sep=""), header=FALSE)
          data=dataorg[,1:i]
          
          #Calculate 14 person fit indices
          
          ####### 1 #######
          #personal point-biserial correlation
          rpb.out <- r.pbis(data, IRT.PModel = "2PL")
          
          ####### 2 #######
          #caution statistic C.Sato
          csato.out <- C.Sato(data, IRT.PModel = "2PL")
          
          ####### 3 #######
          #modified caution statistic Cstar c*
          cstar.out <- Cstar(data, IRT.PModel = "2PL")
          
          ####### 4 #######
          #Number of Guttman errors
          g.out <- G(data, IRT.PModel = "2PL")
          
          ####### 5 #######
          #Gnormed
          gn.out <- Gnormed(data, IRT.PModel = "2PL")
          
          ####### 6 #######
          #Agreement statistic
          akb.out <- A.KB(data, IRT.PModel = "2PL")
          
          ####### 7 #######
          #Disagreement
          dkb.out <- D.KB(data, IRT.PModel = "2PL")
          
          ####### 8 #######
          #Dependability
          ekb.out <- E.KB(data, IRT.PModel = "2PL")
          
          ####### 9 #######
          #U3 van der Flier's statistics
          u3.out <- U3(data, IRT.PModel = "2PL")
          
          ####### 10 #######
          #ZU3
          zu3.out <- ZU3(data, IRT.PModel = "2PL")
          
          ####### 11 #######
          #NCI
          nci.out <- NCI(data, IRT.PModel = "2PL")
          
          ####### 12 #######
          #Ht
          ht.out <- Ht(data, IRT.PModel = "2PL")
          
          ####### 13 #######
          #Iz
          lz.out <- lz(data, IRT.PModel = "2PL")
          
          ####### 14 #######
          #Iz*
          lzstar.out <- lzstar(data, IRT.PModel = "2PL")
          
          temp = data.frame(rpb.out[["PFscores"]], csato.out[["PFscores"]], cstar.out[["PFscores"]], g.out[["PFscores"]], gn.out[["PFscores"]], 
                            akb.out[["PFscores"]], dkb.out[["PFscores"]], ekb.out[["PFscores"]], u3.out[["PFscores"]], zu3.out[["PFscores"]], 
                            nci.out[["PFscores"]], ht.out[["PFscores"]], lz.out[["PFscores"]], lzstar.out[["PFscores"]],dataorg[,2*i+1])
          
          #calculate AUC
          temp3 = NULL #create storage
          for (p in 1:14) {
            
            pfdata = data.frame(temp[,p], temp[,15])
            pfdata <-na.omit(pfdata)
            colnames(pfdata)=c("pfscr","class")
            
            # c1=biased, c2=unbiased
            pfdata$class[pfdata$class==1] <- 1
            pfdata$class[pfdata$class==2] <- 0
            
            #ROC by logistic regression
            glm.fit = glm(pfdata$class ~ pfdata$pfscr, family=binomial)
            roc.glm = roc(pfdata$class, glm.fit$fitted.values)
            
            ## ROC by random forest
            rf.model <- randomForest(factor(pfdata$class) ~ pfscr, pfdata)
            roc.rf = roc(pfdata$class, rf.model$votes[,2])
            
            #ROC by person-fit scores
            roc.pf = roc(pfdata$class, pfdata$pfscr)
            
            temp2 = cbind(roc.glm[["auc"]], roc.rf[["auc"]], roc.pf[["auc"]])
            temp3 = cbind(temp3, temp2)
          }
          
          temp3 = data.frame(cs,cl,i,n,r,temp3)
          danova=rbind(danova,temp3)  
        } #r
      } #n
    } #i
  } #cl
} #cs

end = proc.time()
end - start

write.table(danova, file = "danova.csv", row.names=FALSE, na="", col.names=T, sep=",")


