MMBlup=function(Y,Cofactor = NULL, X = NULL, fmla = NULL,ZList=NULL,VarList,ResMM){
  
  ### Checkings 
  
  ## On Y
  if (!is.numeric(Y)){
    stop('Y should be numeric')
  }
  
  ## On X
  if(!is.null(X)){
    stop('This function does not handle cases where X is not NULL...')
  }
  ## On covariance matrices
  CheckSym <- VarList %>% 
    map_lgl(isSymmetric) %>% 
    all
  if(!CheckSym){
    stop("All covariance matrices should be symetric")
  }
  
  ## On cov and incidence matrices
  CheckZandC <- map2_lgl(ZList, VarList, ~ ncol(.x)==ncol(.y)) %>% 
    all
  if(!CheckZandC){
    stop("Incompatible dimensions between Z and Var")
  }
  
  ## Checks on names
  
  if(any(names(VarList) %in% names(ResMM$NullModel$Sigma2))==FALSE)
  {"Not all random effects are shared by VarList and ResMM"}
  
  
  ### Build the fixed effect matrix
  
  #assign("fmla",fmla,.GlobalEnv)
  Fixed <- .GetFixedEffects(Cofactor=Cofactor,fmla=fmla,Nind=length(Y))
  
  ### Perform Blup estimation
  
  VarComp <- ResMM$NullModel$Sigma
  VarY <- pmap(list(ZList,VarList,VarComp), ~ ..3*tcrossprod(..1,tcrossprod(..1,..2))) %>% 
    Reduce('+',.)
  NormY <- solve(VarY,Y-Fixed%*%ResMM$NullModel$Beta)
  Blup <- pmap(list(VarComp,ZList,VarList), ~ ..1*..3 %*%crossprod(..2,NormY))
  
  ### Output the results
  
  return(Blup)
}