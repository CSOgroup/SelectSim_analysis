# Author: Debora Sesia
library(glmnet)

ENanalysis=function(y,X,n.rep=10,type_regression){
  
  X=as.matrix(X)
  X=X[,which(colSums(X)>=min(nrow(X)*0.002,5))]    
  
  # Run EN
  ENCoeff=list()
  bestalphavector=c()
  
  for (i in 1:n.rep){
    print(i)
    alpha.grid=seq(0.1,0.9,0.1)
    
    min.cvm=c()
    lambda.min=c()
    n.a=1
    for (a in alpha.grid){
      cv=cv.glmnet(X, y, alpha=a, family=type_regression)
      min.cvm=c(min.cvm,min(cv$cvm)); names(min.cvm)[n.a]=a
      lambda.min=c(lambda.min,cv$lambda.min); names(lambda.min)[n.a]=a
      n.a=n.a+1
    }
    
    bestalpha=names(which.min(min.cvm))
    bestalphavector[i]=bestalpha
    bestlam=lambda.min[bestalpha]
    
    fit=glmnet(X,y,alpha=bestalpha,family=type_regression)
    coefficients=coef(fit,s=bestlam)[-1,]
    ENCoeff[[i]]=coefficients
  }
  
  coeff.matrix=matrix(NA,ncol=n.rep,nrow = length(coefficients)); colnames(coeff.matrix)=paste0("rep",1:n.rep); rownames(coeff.matrix)=names(coefficients)
  for(i in 1:n.rep){
    for(j in 1:length(ENCoeff[[i]])){
      gene.temp=names(ENCoeff[[i]][j])
      coeff.matrix[gene.temp,i]=as.numeric(ENCoeff[[i]][j])
    }
  }
  
  coeff.summary=data.frame(matrix(NA,ncol=3,nrow = nrow(coeff.matrix)))
  colnames(coeff.summary)=c("gene","mean coeff","times selected")
  coeff.summary$gene=rownames(coeff.matrix)
  for(i in 1:nrow(coeff.matrix)){
    coeff.summary[i,"mean coeff"]=mean(coeff.matrix[i,])
    coeff.summary[i,"times selected"]=length(which((coeff.matrix[i,]!=0)))
  }
  
  results=list(data.frame(coeff.matrix),coeff.summary,bestalphavector)
  names(results)=c("coeff.matrix","summary","alpha")
  
  return(results)
}