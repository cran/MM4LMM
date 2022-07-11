MMVcov <- function(ResMM,Y,Cofactor=NULL,formula=NULL,ZList=NULL,VarList,information="Expected"){
  NbVar <- length(VarList)
  NamesVar <- names(VarList)
  
  if (!is.null(ZList)){
    ## Check ZList and VarList dimension
    if (length(ZList) != NbVar) stop("Length of ZList and VarList should be the same.")
    VarList <- lapply(1:NbVar , function(k) {
      if (ncol(ZList[[k]]) != nrow(VarList[[k]])) stop(paste0("Dimensions of ZList and VarList element ",k," do not match."))
      if (nrow(ZList[[k]]) != length(Y)) stop(paste0("Dimensions of Y and ZList element ",k," do not match."))
      return(tcrossprod(ZList[[k]] , tcrossprod(ZList[[k]], VarList[[k]])))
    })
  }else{
    invisible(sapply(1:NbVar , function(k) {
      if (nrow(VarList[[k]]) != length(Y)) stop(paste0("Dimensions of Y and VarList element ",k," do not match."))
    }))
  }
  
  if (is.null(Cofactor)){
    Cofactor <- cbind(rep(1,length(Y)))
  }
  
  ## Check on Cofactor
  if (!is.data.frame(Cofactor)){
    if (!is.matrix(Cofactor)){
      stop('Cofactor should be either a matrix or a data.frame')        
    } else {
      Cofactor <- as.data.frame(Cofactor)
      names(Cofactor) <- colnames(Cofactor)
    }
    Names <- names(Cofactor)
    if (is.null(Names)){
      names(Cofactor) <- paste0("Cof",1:ncol(Cofactor))
    } else {
      names(Cofactor) <- gsub('([[:punct:]])|\\s+','',Names)
    }
  }
  CofName <- names(Cofactor)
  
  if (is.null(formula)){
    formulaComp <- as.formula(paste0("~",paste0(CofName,collapse="+")))
    Factors <- c(CofName)
  }else{
    SplitForm <- strsplit(as.character(formula)," \\+ ")
    Factors <- SplitForm[[length(SplitForm)]]
    formulaComp <- formula
  }
  
  Tmp <- model.matrix(terms(formulaComp,keep.order=TRUE),data=as.data.frame(Cofactor))
  QR <- qr(Tmp)
  Fixed <- cbind(Tmp[,QR$pivot[1:QR$rank]])
  NamesFixed <- colnames(Fixed)
  
  ## Check the matrix used for the variance covariance approximation matrix
  if (!information %in% c("AI","Expected")) stop("information can be equal to AI or Expected only.")
  
  
  X <- Fixed
  ResSig <- ResMM[[1]]$Sigma2
  Sigma <- Reduce('+',lapply(1:length(VarList) , function(i) VarList[[i]]*ResSig[i]))
  Sigma_delta <- Sigma/tail(ResSig,n=1)
  Sigma_inv <- solve(Sigma_delta)  
  
  if (information=="AI"){
    P <- Sigma_inv - crossprod(Sigma_inv,crossprod(t(X),solve(crossprod(X,crossprod(Sigma_inv,X)),crossprod(X,Sigma_inv))))
    P <- P/tail(ResSig,n=1)
    PV <- lapply(VarList , function(k) crossprod(k,P))
    PY <- crossprod(P,Y)
    AI <- matrix(0,NbVar,NbVar)
    AI_NoResidual <- AI[-NbVar,-NbVar]
    #if (NbVar>2){
    AI[lower.tri(AI,diag=F)] <- sapply(combn(1:(NbVar),2,simplify=F) , function(x) {
      return(crossprod(Y,crossprod(PV[[x[1]]],crossprod(PV[[x[2]]],PY)))/2)
    })
    AI <- AI + t(AI)
    diag(AI) <- sapply(1:(NbVar) , function(x) {
      return(crossprod(Y,crossprod(PV[[x]],crossprod(PV[[x]],PY)))/2)
    })
    vcovMat <- solve(AI)
    
    return(vcovMat)
  }else{
    if (ResMM[[1]]$Method=="Reml"){
      P <- Sigma_inv - crossprod(Sigma_inv,crossprod(t(X),solve(crossprod(X,crossprod(Sigma_inv,X)),crossprod(X,Sigma_inv))))
      P_all <- P/tail(ResSig,n=1)
      V_P <- lapply(VarList , function(k) crossprod(k,P_all))
      Fisher <- matrix(0,NbVar,NbVar)
      Fisher[lower.tri(Fisher,diag=F)] <- sapply(combn(1:NbVar , 2 , simplify=F) , function(x) {
        return(1/2 * sum(diag(crossprod(t(V_P[[x[1]]]),V_P[[x[2]]]))))
      })
      Fisher <- Fisher + t(Fisher)
      diag(Fisher) <- sapply(1:NbVar , function(x) {
        return(1/2 * sum(diag(crossprod(V_P[[x]]))))
      })
    }else{
      Fisher <- matrix(0,NbVar,NbVar)
      V_Sinv <- lapply(VarList , function(k) crossprod(k,Sigma_inv)/tail(ResSig,n=1))
      Fisher[lower.tri(Fisher,diag=F)] <- sapply(combn(1:NbVar , 2 , simplify=F) , function(x) {
        return(1/2 * sum(diag(crossprod(t(V_Sinv[[x[1]]]),V_Sinv[[x[2]]]))))
      })
      Fisher <- Fisher + t(Fisher)
      diag(Fisher) <- sapply(1:NbVar , function(x) {
        return(1/2 * sum(diag(crossprod(V_Sinv[[x]]))))
      })
    }
    vcovMat <- solve(Fisher)
  }

  colnames(vcovMat) <- rownames(vcovMat) <- NamesVar

return(list(vcovMat = vcovMat , SE = sqrt(diag(vcovMat))))
}
