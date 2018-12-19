.MM_ML <-
  function(Y , Cofactor , X , formula , Factors , VarList , Init=NULL , MaxIter = 10 , CritVar = 1e-3 , CritLogLik = 1e-3 , NbCores=1){
    
    if (is.null(Init)) Init <- rep(1,length(VarList))
    if (is.null(names(VarList))){
      NamesSigma <- paste0("Var",1:length(VarList))
    }else{
      NamesSigma <- names(VarList)
    }
    
    if (is.null(X)){
      
      
      Mat <- Cofactor
      Tmp <- model.matrix(formula,data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      res <- .MLMM(Y, Fixed , VarList , Init , MaxIter, CritVar, CritLogLik)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
      res$Factors <- Factors
      Res <- list(res)
      names(Res) <- "NullModel"
    }else{
      if (is.list(X)){
        Res <-  mclapply(X , function(x) {
          if (is.null(colnames(x))){
            colnames(x) <- paste0("X",1:ncol(X))
          }
          
          Mat <- cbind(Cofactor,x)
          Tmp <- model.matrix(formula,data=as.data.frame(Mat))
          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .MLMM(Y, Fixed , VarList , Init , MaxIter, CritVar, CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- names(X)
      }else{
        Res <-  mclapply(1:ncol(X) , function(i) {
          Mat <- Cofactor
          Mat$Xeffect <- X[,i]
          Tmp <- model.matrix(formula,data=as.data.frame(Mat))
          #colnames(Tmp) <- gsub("Xeffect",i,colnames(Tmp))
          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .MLMM(Y, Fixed , VarList , Init , MaxIter, CritVar, CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- colnames(X)
      }
    }
    return(Res)
  }

.MM_ML2Mat <-
  function(Y , Cofactor , X , formula , Factors , VarList , Init=NULL , MaxIter = 100 , CritVar = 1e-3 , CritLogLik = 1e-3 , NbCores=1){
    
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
    
    if (is.null(X)){
      Mat <- Cofactor
      Tmp <- model.matrix(formula,data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      
      res <- .MM_ML2MatRcpp(Ytilde , Fixed , Passage , List$Diag , Init ,MaxIter, CritVar , CritLogLik)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
      res$Factors <- Factors
      Res <- list(res)
      names(Res) <- "NullModel"
    }else{
      if (is.list(X)){
        Res <-  mclapply(X , function(x) {
          if (is.null(colnames(x))){
            colnames(x) <- paste0("X",1:ncol(X))
          }
          
          Mat <- cbind(Cofactor,x)
          Tmp <- model.matrix(formula,data=as.data.frame(Mat))
          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .MM_ML2MatRcpp(Ytilde , Fixed , Passage , List$Diag , Init ,MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- names(X)
      }else{
        Res <-  mclapply(1:ncol(X) , function(i) {
          Mat <- Cofactor
          Mat$Xeffect <- X[,i]
          Tmp <- model.matrix(formula,data=as.data.frame(Mat))
          #colnames(Tmp) <- gsub("Xeffect",i,colnames(Tmp))
          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .MM_ML2MatRcpp(Ytilde , Fixed , Passage , List$Diag , Init ,MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- colnames(X)
      }
    }
    return(Res)
  }

.MM_Reml <-
  function(Y , Cofactor , X , formula , Factors , VarList , Init=NULL , MaxIter = 100 , CritVar = 1e-3 , CritLogLik = 1e-3 , NbCores=1){
    
    if (is.null(Init)) Init <- rep(1,length(VarList))
    if (is.null(names(VarList))){
      NamesSigma <- paste0("Var",1:length(VarList))
    }else{
      NamesSigma <- names(VarList)
    }
    if (is.null(X)){
      Mat <- Cofactor
      Tmp <- model.matrix(formula,data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      
      res <- .RemlMM(Y, Fixed , VarList , Init , MaxIter, CritVar , CritLogLik)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
      res$Factors <- Factors
      Res <- list(res)
      names(Res) <- "NullModel"
    }else{
      if (is.list(X)){
        Res <-  mclapply(X , function(x) {
          if (is.null(colnames(x))){
            colnames(x) <- paste0("X",1:ncol(X))
          }
          
          Mat <- cbind(Cofactor,x)
          Tmp <- model.matrix(formula,data=as.data.frame(Mat))
          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMM(Y, Fixed , VarList , Init , MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- names(X)
      }else{
        Res <-  mclapply(1:ncol(X) , function(i) {
          Mat <- Cofactor
          Mat$Xeffect <- X[,i]
          Tmp <- model.matrix(formula,data=as.data.frame(Mat))
          #colnames(Tmp) <- gsub("Xeffect",i,colnames(Tmp))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMM(Y, Fixed , VarList , Init , MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- colnames(X)
      }
    }
    return(Res)
  }

.MM_Reml1Mat <-
  function(Y , Cofactor , X , formula , Factors , VarInv , logdetVar , NamesSigma , NbCores=1){
    
    if (is.null(NamesSigma)){
      NamesSigma <- "Var"
    }
    if (is.null(X)){
      Mat <- Cofactor
      Tmp <- model.matrix(formula,data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      
      res <- .RemlMM1Mat(Y, Fixed , VarInv , logdetVar)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
      res$Factors <- Factors
      Res <- list(res)
      names(Res) <- "NullModel"
    }else{
      if (is.list(X)){
        Res <-  mclapply(X , function(x) {
          if (is.null(colnames(x))){
            colnames(x) <- paste0("X",1:ncol(X))
          }
          
          Mat <- cbind(Cofactor,x)
          Tmp <- model.matrix(formula,data=as.data.frame(Mat))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMM1Mat(Y, Fixed , VarInv , logdetVar)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- names(X)
      }else{
        Res <-  mclapply(1:ncol(X) , function(i) {
          Mat <- Cofactor
          Mat$Xeffect <- X[,i]
          Tmp <- model.matrix(formula,data=as.data.frame(Mat))
          #colnames(Tmp) <- gsub("Xeffect",i,colnames(Tmp))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMM1Mat(Y, Fixed , VarInv , logdetVar)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- colnames(X)
      }
    }
    return(Res)
  }


.MM_Reml2Mat <-
  function(Y , Cofactor , X , formula , Factors , VarList , Init=NULL , MaxIter = 100 , CritVar = 1e-3 , CritLogLik = 1e-3 , NbCores=1){
    
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
    if (is.null(X)){
      Mat <- Cofactor
      Tmp <- model.matrix(formula,data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      
      res <- .MM_Reml2MatRcpp(Ytilde , Fixed , Passage , c(List$Diag) , Init ,MaxIter, CritVar , CritLogLik)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]
      res$Factors <- Factors
      Res <- list(res)
      names(Res) <- "NullModel"
    }else{
      if (is.list(X)){
        Res <-  mclapply(X , function(x) {
          if (is.null(colnames(x))){
            colnames(x) <- paste0("X",1:ncol(X))
          }
          
          Mat <- cbind(Cofactor,x)
          Tmp <- model.matrix(formula,data=as.data.frame(Mat))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .MM_Reml2MatRcpp(Ytilde , Fixed , Passage , c(List$Diag) , Init ,MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]          
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- names(X)
      }else{
        Res <-  mclapply(1:ncol(X) , function(i) {
          Mat <- Cofactor
          Mat$Xeffect <- X[,i]
          Tmp <- model.matrix(formula,data=as.data.frame(Mat))
          #colnames(Tmp) <- gsub("Xeffect",i,colnames(Tmp))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .MM_Reml2MatRcpp(Ytilde , Fixed , Passage , c(List$Diag) , Init ,MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]      
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- colnames(X)
      }
    }
    return(Res)
  }

.MM_RemlHen <-
  function(Y , Cofactor , X , formula , Factors , Z , GList , GinvList , Rinv , logdetV , NameVar=NULL , Init=NULL , MaxIter = 100 , CritVar = 1e-3 , CritLogLik = 1e-3 , NbCores=1){
    
    if (is.null(NameVar)){
      NamesSigma <- paste0("Var",1:length(NameVar))
    }else{
      NamesSigma <- NameVar
    }
    if (is.null(X)){
      Mat <- Cofactor
      Tmp <- model.matrix(formula,data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      
      res <- .RemlMMHen(Y, Fixed , Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, CritVar , CritLogLik)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]      
      res$Factors <- Factors
      Res <- list(res)
      names(Res) <- "NullModel"
    }else{
      if (is.list(X)){
        Res <-  mclapply(X , function(x) {
          if (is.null(colnames(x))){
            colnames(x) <- paste0("X",1:ncol(X))
          }
          
          Mat <- cbind(Cofactor,x)
          Tmp <- model.matrix(formula,data=as.data.frame(Mat))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMMHen(Y, Fixed , Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]    
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- names(X)
      }else{
        Res <-  mclapply(1:ncol(X) , function(i) {
          Mat <- Cofactor
          Mat$Xeffect <- X[,i]
          Tmp <- model.matrix(formula,data=as.data.frame(Mat))
          #colnames(Tmp) <- gsub("Xeffect",i,colnames(Tmp))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMMHen(Y, Fixed , Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]    
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- colnames(X)
      }
    }
    return(Res)
  }

.MM_RemlHenDiag <-
  function(Y , Cofactor , X , formula , Factors , Z , GList , GinvList , Rinv , logdetV , NameVar=NULL , Init=NULL , MaxIter = 10 , CritVar = 1e-3 , CritLogLik = 1e-3 , NbCores=1){
    
    if (is.null(NameVar)){
      NamesSigma <- paste0("Var",1:length(NameVar))
    }else{
      NamesSigma <- NameVar
    }
    if (is.null(X)){
      Mat <- Cofactor
      Tmp <- model.matrix(formula,data=as.data.frame(Mat))
      QR <- qr(Tmp)
      Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
      NamesFixed <- colnames(Fixed)
      
      res <- .RemlMMHenDiag(Y, Fixed , Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, CritVar , CritLogLik)
      
      names(res$Beta) <- NamesFixed
      rownames(res$VarBeta) <- NamesFixed
      names(res$Sigma2) <- NamesSigma
      res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]      
      res$Factors <- Factors
      Res <- list(res)
      names(Res) <- "NullModel"
    }else{
      if (is.list(X)){
        Res <-  mclapply(X , function(x) {
          if (is.null(colnames(x))){
            colnames(x) <- paste0("X",1:ncol(X))
          }
          
          Mat <- cbind(Cofactor,x)
          Tmp <- model.matrix(formula,data=as.data.frame(Mat))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMMHenDiag(Y, Fixed , Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]      
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- names(X)
      }else{
        Res <-  mclapply(1:ncol(X) , function(i) {
          Mat <- Cofactor
          Mat$Xeffect <- X[,i]
          Tmp <- model.matrix(formula,data=as.data.frame(Mat))
          #colnames(Tmp) <- gsub("Xeffect",i,colnames(Tmp))

          QR <- qr(Tmp)
          Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
          NamesFixed <- colnames(Fixed)
          
          res <- .RemlMMHenDiag(Y, Fixed , Z , GList , GinvList , Rinv , logdetV , Init , MaxIter, CritVar , CritLogLik)
          
          names(res$Beta) <- NamesFixed
          rownames(res$VarBeta) <- NamesFixed
          names(res$Sigma2) <- NamesSigma
          res$attr <- attr(Tmp,"assign")[QR$pivot[1:QR$rank]]      
          res$Factors <- Factors
          return(res)
        } ,mc.cores=NbCores)
        names(Res) <- colnames(X)
      }
    }
    return(Res)
  }

