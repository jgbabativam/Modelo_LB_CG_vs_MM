rm(list = ls())
library(here)
library(tidyverse)
library(bcv)
library(BiplotML)
library(data.table)
library(pracma)
library(MASS)

source(here::here("./R/0. functions.R"))

params <- tribble(
  ~n,  ~p, ~k,   ~D,  ~C,
 100,  50,  3,  0.5,  20
 300,  50,  3,  0.5,  20,
 500,  50,  3,  0.5,  20,
1000,  50,  3,  0.5,  20,
 100, 100,  3,  0.5,  20,
 300, 100,  3,  0.5,  20,
 500, 100,  3,  0.5,  20,
1000, 100,  3,  0.5,  20,
 100, 200,  3,  0.5,  20,
 300, 200,  3,  0.5,  20,
 500, 200,  3,  0.5,  20,
 100,  50,  3,  0.1,  20,  #Desbalance 0.3
 300,  50,  3,  0.1,  20,  #Desbalance 0.3
 500,  50,  3,  0.1,  20,  #Desbalance 0.3
1000,  50,  3,  0.1,  20,  #Desbalance 0.3
 100, 100,  3, 0.15,  20,  #Desbalance 0.3
 300, 100,  3, 0.15,  20,  #Desbalance 0.3
 500, 100,  3, 0.15,  20,  #Desbalance 0.3
1000, 100,  3, 0.15,  20,  #Desbalance 0.3
 100, 200,  3,  0.2,  20,
 300, 200,  3,  0.2,  20,
 500, 200,  3,  0.2,  20,
 100,  50,  3, 0.02,  20,  #Desbalance 0.2
 300,  50,  3, 0.02,  20,  #Desbalance 0.2
 500,  50,  3, 0.02,  20,  #Desbalance 0.2
1000,  50,  3, 0.02,  20,  #Desbalance 0.2
 100, 100,  3, 0.055, 20,  #Desbalance 0.2
 300, 100,  3, 0.055, 20,  #Desbalance 0.2
 500, 100,  3, 0.055, 20,  #Desbalance 0.2
1000, 100,  3, 0.055, 20,  #Desbalance 0.2
 100, 200,  3, 0.085, 20,
 300, 200,  3, 0.085, 20,
 500, 200,  3, 0.085, 20,
 100,  50,  3, 0.002, 20,  #Desbalance 0.1
 300,  50,  3, 0.002, 20,  #Desbalance 0.1
 500,  50,  3, 0.002, 20,  #Desbalance 0.1
1000,  50,  3, 0.002, 20,  #Desbalance 0.1
 100, 100,  3, 0.009, 20,  #Desbalance 0.1
 300, 100,  3, 0.009, 20,  #Desbalance 0.1
 500, 100,  3, 0.009, 20,  #Desbalance 0.1
1000, 100,  3, 0.009, 20,  #Desbalance 0.1
 100, 200,  3, 0.025, 20,
 300, 200,  3, 0.025, 20,
 500, 200,  3, 0.025, 20,
)

###-----------------Simular X, \Theta, \Pi.

set.seed(12345)
samples <- function(n, p, k, D, C){
  xs <- simBin(n = n, p = p, k = k, D = D, C = C)
  return(xs)
}

###-----------------------------------.
### Monte Carlo
library(future)
library(furrr)
library(future.apply)
library(parallel)
library(beepr)

plan(multiprocess)

inicio <- Sys.time()
### Estimo los errores para R remuestras usando programación funcional
sale <- furrr::future_map(1:2, function(x){
          lsamples = params %>% purrr::pmap(samples)
          out <- future_mapply(function(x, dimen){
            #--- BFGS Algorithm
            BFGS <- BiplotML::LogBip(x$X, k = dimen, method = "BFGS", plot = FALSE, random_start = FALSE)
            ThetaBFGS <- BiplotML::fitted_LB(BFGS, type = "link")
            predBFGS <- BiplotML::pred_LB(BFGS, x = x$X, ncuts = 50)
            rBFGS <- predBFGS$BACC; rmse_BFGS <- rmse(x$Theta, ThetaBFGS)
            cvBFGS <- crossval(x$X, k = dimen, thres = predBFGS$thresholds, method = "BFGS")
            #--- Fletcher-Reeves Algorithm     
            FR <- BiplotML::LogBip(x$X, k = dimen, method = "CG", type = 1, plot = FALSE, random_start = FALSE)
            ThetaFR <- BiplotML::fitted_LB(FR, type = "link")
            predFR <- BiplotML::pred_LB(FR, x = x$X, ncuts = 50)
            fr_CG <- predFR$BACC; rmse_FR <- rmse(x$Theta, ThetaFR)
            cvFR <- crossval(x$X, k = dimen, thres = predFR$thresholds, method = "CG", type = 1)
            #--- Polak--Ribiere Algorithm    
            PR <- BiplotML::LogBip(x$X, k = dimen, method = "CG", type = 2, plot = FALSE, random_start = FALSE)
            ThetaPR <- BiplotML::fitted_LB(PR, type = "link")
            predPR <- BiplotML::pred_LB(PR, x = x$X, ncuts = 50)
            pr_CG <- predPR$BACC; rmse_PR <- rmse(x$Theta, ThetaPR)
            cvPR <- crossval(x$X, k = dimen, thres = predPR$thresholds, method = "CG", type = 2)
            #--- Beale--Sorenson Algorithm    
            BS <- BiplotML::LogBip(x$X, k = dimen, method = "CG", type = 3, plot = FALSE, random_start = FALSE)
            ThetaBS <- BiplotML::fitted_LB(BS, type = "link")
            predBS <- BiplotML::pred_LB(BS, x = x$X, ncuts = 50)
            bs_CG <- predBS$BACC; rmse_BS <- rmse(x$Theta, ThetaBS)
            cvBS <- crossval(x$X, k = dimen, thres = predBS$thresholds, method = "CG", type = 3)
            #--- Dai--Yuan Algorithm    
            DY <- BiplotML::LogBip(x$X, k = dimen, method = "CG", type = 4, plot = FALSE, random_start = FALSE)
            ThetaDY <- BiplotML::fitted_LB(DY, type = "link")
            predDY <- BiplotML::pred_LB(DY, x = x$X, ncuts = 50)
            DY_CG <- predDY$BACC; rmse_DY <- rmse(x$Theta, ThetaDY)
            cvDY <- crossval(x$X, k = dimen, thres = predDY$thresholds, method = "CG", type = 4)
            #--- Majorization Minimization algorithm    
            tMM <- BiplotML::LogBip(x$X, k = dimen, method = "MM", plot = FALSE)
            ThetaMM <- BiplotML::fitted_LB(tMM, type = "link")
            predMM <- BiplotML::pred_LB(tMM, x = x$X, ncuts = 50)
            CD_MM <- predMM$BACC; rmse_MM <- rmse(x$Theta, ThetaMM)
            cvMM <- crossval(x$X, k = dimen, thres = predMM$thresholds, method = "MM")
            #----- Salida.
            outr <- list(n = x$n, p = x$p, k = dimen,
                         D_t = x$D_t, D_r = round(x$D_r, 2),  # n.filas, n.columnas, Desb.teorico, Desb.Real   
                         BFGS = rBFGS, fr_CG = fr_CG, pr_CG = pr_CG, bs_CG=bs_CG, DY_CG = DY_CG, MM = CD_MM, 
                         rmse_BFGS = rmse_BFGS, rmse_FR = rmse_FR, rmse_PR = rmse_PR, rmse_BS = rmse_BS, rmse_DY = rmse_DY, rmse_MM = rmse_MM,
                         cvTBFGS = cvBFGS$cvT, cvTFR = cvFR$cvT, cvTPR = cvPR$cvT, cvTBS = cvBS$cvT, cvTDY = cvDY$cvT, cvTMM = cvMM$cvT,
                         cvDBFGS = cvBFGS$cvD, cvDFR = cvFR$cvD, cvDPR = cvPR$cvD, cvDBS = cvBS$cvD, cvDDY = cvDY$cvD, cvDMM = cvMM$cvD)
            return(outr)
                 }, lsamples, 1:5) 
        }, .progress = TRUE)



resulta <- as.data.frame(matrix(unlist(sale), ncol=29, byrow=TRUE))
colnames(resulta) <- c("n", "p", "k", "desb.teorico", "desb.real", 
                       "BFGS", "CG_Fletcher", "CG_Polak", "CG_Beale", "CG_DY", "MM",
                       "RMSE_BFGS", "RMSE_FR", "RMSE_PR", "RMSE_BS", "RMSE_DY", "RMSE_MM",
                       "cvT.BFGS", "cvT.FR", "cvT.PR", "cvT.BS", "cvT.DY", "cvT.MM",
                       "cvD.BFGS", "cvD.FR", "cvD.PR", "cvD.BS", "cvD.DY", "cvD.MM")

#save(resulta, file=here::here("data/results.rda"))
fin <- Sys.time()
fin - inicio
beepr::beep(8)

###---- Medidas finales de precisión y error
resumen <- resulta %>% group_by(n, p, k, desb.teorico) %>% 
           summarise_all(mean) %>% ungroup()

resumen %>% 
  dplyr::select(n, p, k, desb.teorico, desb.real, starts_with("cvD")) %>% 
  pivot_longer(-c("n", "p", "k", "desb.teorico", "desb.real"), names_to = "Algoritmo", values_to = "Error") %>%
  ggplot(aes(x = k, y = Error, group = Algoritmo)) +
  geom_line(aes(linetype=Algoritmo, color=Algoritmo)) +
  geom_point(aes(color=Algoritmo))+
  scale_y_continuous(breaks = seq(0,100,5)) +
  scale_x_continuous(breaks = seq(0, 5, 1)) +
  labs(y = "BACC (%)", x = "k") + theme_bw() +
  facet_grid( n ~ p)
  
####---- Fin

