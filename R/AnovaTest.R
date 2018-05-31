AnovaTest <- 
function(ResMMEst , TestedCombination=NULL , Type = "TypeIII" , NbCores=1){

	if (length(ResMMEst[[1]]) != 6) stop("ResMMEst should be the output of the function MMEst")

	if (!is.null(TestedCombination)){
		message("Wald tests about the combination are computed")
		if (is.matrix(TestedCombination)) TestedCombination <- list(TestedCombination)
		if (!is.list(TestedCombination)) stop("TestedCombination should be either a list of matrices or a matrix")

		Res <- mclapply(ResMMEst , function(x) {
			Beta <- x$Beta
			VarBeta <- x$VarBeta

			Test <- t(sapply(TestedCombination , function(C) {
				if (ncol(C)!=length(Beta)){
					Stat <- NA
					Pval <- NA
				}else{
					rgC <- qr(C)$rank
					CBeta <- tcrossprod(C,t(Beta))
					Stat <- crossprod(CBeta, solve(tcrossprod(C,tcrossprod(C,VarBeta)), CBeta))

					Pval <- pchisq(Stat , df = rgC , lower.tail=FALSE)
				}
				return(c(Stat,Pval))			
			}))
			colnames(Test) <- c("Wald","pval")
			rownames(Test) <- names(TestedCombination)

			return(Test)
		},mc.cores=NbCores)
		names(Res) <- names(ResMMEst)	
	}else{
		if (Type=="TypeI"){
			Res <- mclapply(ResMMEst , function(x) {
				Beta <- x$Beta
				Names <- sapply(strsplit(names(Beta),"_:_"), function(y) y[[1]])
				CommonName <- lapply(unique(Names) , function(y) which(Names==y))
				names(CommonName) <- unique(Names)
				MatTI <- lapply(unique(Names) , function(y){
					mat <- as.numeric(1:length(Beta) %in% CommonName[[y]])
					return(mat)				
				})

				VarBeta <- x$VarBeta
				Chol <- chol(solve(VarBeta))
				TypeIcomp <- tcrossprod(Chol,t(Beta))^2
				Test <- t(sapply(MatTI , function(y) {
					if (length(y)!=length(Beta)){
						Stat <- NA
						Pval <- NA
					}else{
						Stat <- sum(y*TypeIcomp)
						Pval <- pchisq(Stat , df=sum(y) , lower.tail=FALSE)
					}
					return(c(Stat,Pval))
				}))

				colnames(Test) <- c("Wald (Type I)","pval")
				rownames(Test) <- unique(Names)
				return(Test)
			},mc.cores=NbCores)
			names(Res) <- names(ResMMEst)
		}
		if (Type=="TypeIII"){
			Res <- mclapply(ResMMEst , function(x) {
				Beta <- x$Beta
				Names <- sapply(strsplit(names(Beta),"_:_"), function(y) y[[1]])
				CommonName <- lapply(unique(Names) , function(y) which(Names==y))
				names(CommonName) <- unique(Names)
				MatTIII <- lapply(unique(Names) , function(y){
					EffectSelected <- CommonName[[y]]
					mat <- t(sapply(EffectSelected , function(z) as.numeric(1:length(Beta)==z)))
					return(mat)				
				})

				VarBeta <- x$VarBeta

				Test <- t(sapply(MatTIII , function(C) {
					if (ncol(C)!=length(Beta)){
						Stat <- NA
						Pval <- NA
					}else{
						rgC <- qr(C)$rank
						CBeta <- tcrossprod(C,t(Beta))
						Stat <- crossprod(CBeta, solve(tcrossprod(C,tcrossprod(C,VarBeta)), CBeta))
						Pval <- pchisq(Stat , df = rgC , lower.tail=FALSE)
					}
					return(c(Stat,Pval))
				}))

				colnames(Test) <- c("Wald (Type III)","pval")
				rownames(Test) <- unique(Names)
				return(Test)
			},mc.cores=NbCores)
			names(Res) <- names(ResMMEst)
		}
		if ((Type!="TypeI")&&(Type!="TypeIII")) warning("AnovaTest computes TypeI or TypeIII test when TestedCombination is NULL")
	}
	CompNA <- sum(unlist(mclapply(Res,function(x) sum(is.na(x[,2])) , mc.cores=NbCores)))
	if (CompNA > 0) message(paste0(CompNA," tests were not compute because of incompatible dimension"))
	return(Res)
}



