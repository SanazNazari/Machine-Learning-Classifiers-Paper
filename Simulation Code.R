
library(glue)
library(MplusAutomation)

#r = replication
#n = sample size
#i = number of items
#cl = cross loading
#cs = class separation

start = proc.time()

for (r in 1:100) {
  for (n in c(500, 750, 1000)) {
    for (i in c(10, 20, 30)) {
      for (cl in c(0.1, 0.3, 0.5)){
        for (cs in c(1, 2, 3)){
          
    fmm  <- mplusObject(
            
    TITLE = glue("SDB 2 classes"), 
            
    MONTECARLO = glue(
    "names are f1-f{i} s1-s{i};
    !threshold number = 1
    generate = f1-f{i}(1) s1-s{i}(1);
    categorical = f1-f{i} s1-s{i};
    genclasses = c(2);
    classes = c(2);
    nobs = {n};
    nrep = 1;
    REPSAVE = ALL;
    save = sdb-{cs}-{cl}-{i}-{n}-{r}.dat;"),
            
    ANALYSIS = 
    "type = mixture;
    ALGORITHM=INTEGRATION;",
            
    MODELPOPULATION = glue(
    "%overall%
    
    sdb BY s1@1 s2-s{i}*0.8;
    foc BY f1@1 f2-f{i}*0.8;
    sdb by f1@1 f2-f{i}*{cl};
    
    [sdb@0];
    [foc@0];
    sdb*1;
    foc*1;
    
    sdb with foc@0;
    [c#1*0];
     
    %c#1%
       
    %c#2%
    sdb by f1-f{i}@0;
    [sdb@{cs}];
    [foc@{cs}];"),
            
    MODEL = glue(
    "%overall%
    sdb BY s1@1 s2-s{i};
    foc BY f1@1 f2-f{i};
    sdb by f1@1 f2-f{i};
    
    sdb;
    foc;
    sdb with foc@0;
    
    %c#1%
      
    %c#2%
    sdb by f1-f{i}@0;
    sdb;
    foc;")
    )
          
    mplusModeler(fmm,dataout=glue("dataout.dat"),
                modelout=glue("modelout-{cs}-{cl}-{i}-{n}-{r}.inp") ,run = TRUE)  
          
        }
      }
    }
  }
}
    
end = proc.time()
end - start

        
        
        
        
        