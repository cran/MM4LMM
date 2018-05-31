.MM_ML <-
  function(Y , ListX , VarList , Init=NULL , MaxIter = 10 , Crit = 10e-3 , NbCores=1){
    
    if (is.null(Init)) Init <- rep(1,length(VarList))
    if (is.null(names(VarList))){
      NamesSigma <- paste0("Var",1:length(VarList))
    }else{
      NamesSigma <- names(VarList)
    }
    
    
    Res <-  mclapply(ListX , function(X) {
      if (is.null(colnames(X))){
        colnames(X) <- paste0("X",1:ncol(X))
      }
      
      NamesFixed <- colnames(X)
      
      res <- .MLMM(Y, X , VarList , Init , MaxIter, Crit)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      
      return(res)
    } ,mc.cores=NbCores)
    
    return(Res)
  }

.MM_ML2Mat <-
  function(Y , ListX , VarList , Init=NULL , MaxIter = 100 , Crit = 10e-3 , NbCores=1){
    
    if (is.null(Init)) Init <- rep(1,2)
    
    K1 <- VarList[[1]]
    K2 <- VarList[[2]]
    List <- .PrepMat(Y , K1 , K2)
    Ytilde <- List$Ytilde
    Passage <- t(List$U)
    
    if (is.null(names(VarList))){
      NamesSigma <- paste0("Var",1:length(VarList))
    }else{
      NamesSigma <- names(VarList)
    }
    
    Res <- mclapply(ListX , function(X) {
      if (is.null(colnames(X))){
        colnames(X) <- paste0("X",1:ncol(X))
      }
      NamesFixed <- colnames(X)
      res <- .MM_ML2MatRcpp(Ytilde , X , Passage , List$Diag , Init ,MaxIter, Crit )
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- colnames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      
      return(res)}
      ,mc.cores=NbCores)
    
    return(Res)
  }

.MM_Reml <-
  function(Y , ListX , VarList , Init=NULL , MaxIter = 10 , Crit = 10e-3 , NbCores=1){
    
    if (is.null(Init)) Init <- rep(1,length(VarList))
    if (is.null(names(VarList))){
      NamesSigma <- paste0("Var",1:length(VarList))
    }else{
      NamesSigma <- names(VarList)
    }
    
    Res <-  mclapply(ListX , function(X) {
      if (is.null(colnames(X))){
        colnames(X) <- paste0("X",1:ncol(X))
      }
      NamesFixed <- colnames(X)
      res <- .RemlMM(Y, X , VarList , Init , MaxIter, Crit)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- colnames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      
      return(res)
    } ,mc.cores=NbCores)
    
    return(Res)
  }

.MM_Reml2Mat <-
  function(Y , ListX , VarList , Init=NULL , MaxIter = 100 , Crit = 10e-3 , NbCores=1){
    
    if (is.null(Init)) Init <- rep(1,1)
    
    K1 <- VarList[[1]]
    K2 <- VarList[[2]]
    List <- .PrepMat(Y , K1 , K2)
    Ytilde <- List$Ytilde
    Passage <- t(List$U)
    
    
    
    
    if (is.null(names(VarList))){
      NamesSigma <- paste0("Var",1:length(VarList))
    }else{
      NamesSigma <- names(VarList)
    }
    
    Res <- mclapply(ListX , function(X) {
      if (is.null(colnames(X))){
        colnames(X) <- paste0("X",1:ncol(X))
      }
      
      NamesFixed <- colnames(X)
      
      res <- .MM_Reml2MatRcpp(Ytilde , X , Passage , c(List$Diag) , Init ,MaxIter, Crit )
      
      
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- colnames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      
      return(res)}
      ,mc.cores=NbCores) 
    return(Res)
  }

.MM_RemlHen <-
  function(Y , ListX , Z , GList , GinvList , Rinv , logdetV , NameVar=NULL , Init=NULL , MaxIter = 10 , Crit = 10e-3 , NbCores=1){
    
    if (is.null(NameVar)){
      NamesSigma <- paste0("Var",1:length(NameVar))
    }else{
      NamesSigma <- NameVar
    }
    
    
    Res <-  mclapply(ListX , function(X) {
      if (is.null(colnames(X))){
        colnames(X) <- paste0("X",1:ncol(X))
      }
      
      NamesFixed <- colnames(X)
      
      res <- .RemlMMHen(Y, X , Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, Crit)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- colnames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      
      return(res)
    } ,mc.cores=NbCores)
    
    return(Res)
  }

.MM_RemlHenDiag <-
  function(Y , ListX , Z , GList , GinvList , Rinv , logdetV , NameVar=NULL , Init=NULL , MaxIter = 10 , Crit = 10e-3 , NbCores=1){
    
    if (is.null(NameVar)){
      NamesSigma <- paste0("Var",1:length(NameVar))
    }else{
      NamesSigma <- NameVar
    }
    
    
    Res <-  mclapply(ListX , function(X) {
      if (is.null(colnames(X))){
        colnames(X) <- paste0("X",1:ncol(X))
      }
      
      NamesFixed <- colnames(X)
      res <- .RemlMMHenDiag(Y, X , Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, Crit)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- colnames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      
      return(res)
    } ,mc.cores=NbCores)
    
    return(Res)
  }

