
##############################################################
### Simula matriz binaria con diferentes grados de desbalance

# n: Número de filas
# p: Número de columnas
# k: Número de dimensiones subyacentes
# D: Proporción esperada de 1's en la matriz, desbalanceo.
# C: Ampliar escala

library(data.table)
library(pracma)
library(MASS)
library(BiplotML)

## Usando agrupaciones
simBin.Clus <- function(n, p, k, D, C = 1){
  Bp <- matrix(rnorm(p*k), p)
  B <- gramSchmidt(Bp)$Q  #De esta manera t(B)B = I_k
  
  S <-  diag(k)
  mu <- c(rep(qlogis(D), p))
  
  res <- sapply(centros, 
                function(x)
                  mvtnorm::rmvnorm(n/k, mean = x, sigma = S), simplify = F)
  A <- as.matrix(rbindlist(lapply(res, as.data.frame)))
  
  logOdds <- rep(1,n)%*%t(mu) + C *(A %*% t(B))
  P <- plogis(logOdds)
  
  M <- matrix(NA, n, p)
  X <- t(sapply(1:n, function(i) sapply(1:p, function(j) M[i, j] = rbinom(1, 1, P[i,j]))))
  
  desb <- sum(X)/length(X)
  
  out <- list(X = X, P = P, Theta = logOdds, A = A, B = B, mu = mu, D_t = D, D_r = desb, n = n, p = p)  
}


## Sin agrupaciones
simBin <- function(n, p, k, D, C = 1){
  Bp <- matrix(rnorm(p*k), p)
  B <- pracma::gramSchmidt(Bp)$Q  #De esta manera t(B)B = I_k
  
  S <-  diag(k)
  mu <- c(rep(qlogis(D), p))
  
  centro <- c(0, 0, 0) 
  A <- mvtnorm::rmvnorm(n, mean = centro, sigma = S)

  logOdds <- rep(1,n)%*%t(mu) + C *(A %*% t(B))
  P <- plogis(logOdds)
  
  M <- matrix(NA, n, p)
  X <- t(sapply(1:n, function(i) sapply(1:p, function(j) M[i, j] = rbinom(1, 1, P[i,j]))))
  
  desb <- sum(X)/length(X)
  
  out <- list(X = X, P = P, Theta = logOdds, A = A, B = B, mu = mu, D_t = D, D_r = desb, n = n, p = p)  
  
}

#====================================================
# Selección optima del umbral

thresholds <- function(x, P, ncuts = 100){
  P <- as.matrix(P)
  if(is.null(colnames(x))) colnames(x) <- paste0('V', 1:ncol(x))
  
  N1 <- ceiling(sum(x, na.rm = T))
  N0 <- length(as.matrix(x)) - N1
  TE <- lapply(seq(0, 1, length.out = ncuts), function(z){
    Pr <- ifelse(P>=z, 1, 0)
    c1 <- 1 - apply((Pr == 1) & (x == 1), 2, sum, na.rm=TRUE)/apply(x == 1, 2, sum)
    c2 <- 1 - apply((Pr == 0) & (x == 0), 2, sum, na.rm=TRUE)/apply(x == 0, 2, sum)
    TE <- 100/2 * (c1 + c2)
    lista <- list(TE)
    
  })
  TEp <-  data.frame(dplyr::bind_rows(TE), threshold = seq(0, 1, length.out = ncuts))
  
  thresholds <- TEp %>%
    tidyr::pivot_longer(-threshold, names_to = "variable", values_to = "BACC") %>%
    group_by(variable) %>%
    mutate(merror = min(BACC)) %>%
    dplyr::filter(BACC == merror) %>%
    mutate(row = dplyr::row_number()) %>%
    filter(row == 1) %>% ungroup %>%
    dplyr::select(variable, threshold, BACC)
  
  thresholds <- thresholds[match(colnames(x), thresholds$variable),]
  Pr <- matrix(NA, nrow(P), ncol(P))
  for(p in 1:ncol(P)){
    Pr[,p] <- ifelse( P[,p] >= thresholds$threshold[p], 1, 0)
  }
  c1 <- 1 - sum(apply((Pr == 1) & (x == 1), 2, sum, na.rm=TRUE))/N1
  c2 <- 1 - sum(apply((Pr == 0) & (x == 0), 2, sum, na.rm=TRUE))/N0
  BACC <- round(100/2 * (c1 + c2), 2)
  
  out <- list(pred = Pr, thres = thresholds, BACC = BACC)
}


#======================================================
# selección de pliegues - patron de valores perdidos Wold
# Bro et. al (2008)...k, k+K, k+2K, etc

train_miss_pattern <- function(data, K = 7){
  X <- as.matrix(data)
  out <- sapply(1:K, function(x) seq(x, length(X), by = K))
  train <- list()
  for(i in 1:K){
    temp <- X
    temp[out[[i]]] <- NA
    train[[i]] <- temp
  }
  
  lout <- list(missWold = out, Xtrain = train)
  return(lout)
}

#=======================================================.
# Cross-Validate

cv_LogBip <- function(data, k = 2, K = 7, method = NULL, type = NULL){
  
  x <- as.matrix(data)
  min <- min(k)
  j <- k
  
  cvT <- matrix(NA, length(k), 2)
  cvD <- matrix(NA, length(k), 2)
  for(k in j){

  if(method == "CG" & k > 0){
    bip <- BiplotML::LogBip(x, k = k, method = method, type = type, plot = FALSE)
    thres <- BiplotML::pred_LB(bip, x, ncuts = 50)$thresholds
  }else if(method == "BFGS" & k > 0){
    bip <- BiplotML::LogBip(x, k = k, method = method, plot = FALSE)
    thres <- BiplotML::pred_LB(bip, x, ncuts = 50)$thresholds
  }else if(method == "MM" & k > 0){
     bip <- BiplotML::LogBip(x, k = k, method = MM, plot = FALSE)
     thres <- BiplotML::pred_LB(bip, x, ncuts = 50)$thresholds
  }else{
    P <- rep(1, nrow(x)) %*% t(as.matrix(colMeans(x, na.rm=TRUE)))
    thres <- thresholds(x = x, P = P, ncuts = 50)$thres
  }
  
  folds <- train_miss_pattern(x, K = K)
  missVal <- folds$missWold
  train <- folds$Xtrain
 
  cv_errT = list()
  cv_errD = list()
  for(i in 1:K){  
  if(method == "CG" & k > 0){
    missBip <- BiplotML::LogBip(train[[i]], k = k, 
                                method = method, type = type, plot = FALSE)
    P <- fitted_LB(missBip, type = "response")
  }else if(method == "BFGS" & k > 0){
    missBip <- BiplotML::LogBip(train[[i]], k = k,
                                method = method, plot = FALSE)
    P <- fitted_LB(missBip, type = "response")
  }else if(method == "MM" & k > 0){
    missBip <- BiplotML::LogBip(train[[i]], k = k, method="MM", plot = FALSE)
    P <- fitted_LB(missBip, type = "response")    
  }else{
    P <- rep(1, nrow(x)) %*% t(as.matrix(colMeans(train[[i]], na.rm=TRUE)))
  }
  
  Xhat <- matrix(NA, nrow(P), ncol(P))
  for(p in 1:ncol(P)){
      Xhat[,p] <- ifelse(P[,p] >= thres$threshold[p], 1, 0)
  }
  Xhat[is.na(Xhat)] <- 0
  
  Xhat_pred <- Xhat[missVal[[i]]]
  xReal <- x[missVal[[i]]]
  
  n1 <- sum(xReal)  
  n0 <- length(xReal) - n1
  cv_errT[[i]] = 100 - 100*sum(ifelse((xReal==1 & Xhat_pred==1) | (xReal==0 & Xhat_pred==0), 1, 0))/length(xReal)
  
  err0 <- sum(ifelse(xReal == 0 & Xhat_pred == 1, 1, 0))/n0 
  err1 <- sum(ifelse(xReal == 1 & Xhat_pred == 0, 1, 0))/n1
  
  cv_errD[[i]] <- 100/2 * (err0 + err1)
  }
  row <- k - min + 1 
  cvT[row, 1] <- k; cvD[row, 1] <- k
  cvT[row, 2] <- round(mean(sapply(cv_errT, mean, na.rm = TRUE), na.rm = TRUE), 2)   
  cvD[row, 2] <- round(mean(sapply(cv_errD, mean, na.rm = TRUE), na.rm = TRUE), 2)  
  
  }
  cvT <- as.data.frame(cvT)
  colnames(cvT) <- c("k", "BACC")
  
  cvD <- as.data.frame(cvD)
  colnames(cvD) <- c("k", "BACC")
  
  out <- list(cvT = cvT, cvD = cvD)
  return(out)
}
  

###----  RMS(\theta) = ||\theta - \hat{\theta}||^2/||\theta||^2

rmse <- function(X, Xhat){
  X <- as.matrix(X)
  Xhat <- as.matrix(Xhat)
  
  #--- Norma de Frobenius.
  rmse <- norm(X - Xhat, type = "F")^2 / norm(X, type = "F")^2
  return(rmse)
}



### Función del Cv especial para esta simulación

crossval <- function(x, k = 2, K = 7, thres = NULL, method = NULL, type = NULL){
  
    folds <- train_miss_pattern(x, K = K)
    missVal <- folds$missWold
    train <- folds$Xtrain
    
    cv_errT = list()
    cv_errD = list()
    for(i in 1:K){  
      if(method == "CG" & k > 0){
        missBip <- BiplotML::LogBip(train[[i]], k = k, 
                                        method = method, type = type, plot = FALSE)
        P <- fitted_LB(missBip, type = "response")
      }else if(method == "BFGS" & k > 0){
        missBip <- BiplotML::LogBip(train[[i]], k = k,
                                        method = method, plot = FALSE)
        P <- fitted_LB(missBip, type = "response")
      }else if(method == "MM" & k > 0){
        missBip <- BiplotML::LogBip(train[[i]], k = k, method = "MM", plot = FALSE)
        P <- fitted_LB(missBip, type = "response")    
      }else{
        P <- rep(1, nrow(x)) %*% t(as.matrix(colMeans(train[[i]], na.rm=TRUE)))
      }
      
      Xhat <- matrix(NA, nrow(P), ncol(P))
      for(p in 1:ncol(P)){
        Xhat[,p] <- ifelse(P[,p] >= thres$threshold[p], 1, 0)
      }
      Xhat[is.na(Xhat)] <- 0
      
      Xhat_pred <- Xhat[missVal[[i]]]
      xReal <- x[missVal[[i]]]
      
      n1 <- sum(xReal)  
      n0 <- length(xReal) - n1
      cv_errT[[i]] = 100 - 100*sum(ifelse((xReal==1 & Xhat_pred==1) | (xReal==0 & Xhat_pred==0), 1, 0))/length(xReal)
      
      err0 <- sum(ifelse(xReal == 0 & Xhat_pred == 1, 1, 0))/n0 
      err1 <- sum(ifelse(xReal == 1 & Xhat_pred == 0, 1, 0))/n1
      
      cv_errD[[i]] <- 100/2 * (err0 + err1)
    }
    cvT <- round(mean(sapply(cv_errT, mean, na.rm = TRUE), na.rm = TRUE), 2)   
    cvD <- round(mean(sapply(cv_errD, mean, na.rm = TRUE), na.rm = TRUE), 2)  
    
  out <- list(cvT = cvT, cvD = cvD)
  return(out)
}
