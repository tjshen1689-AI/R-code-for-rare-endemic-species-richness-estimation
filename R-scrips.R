########################################################################################################
########################################################################################################
#R code for conducting area-based estimation of total and rare endemic species richness
#Using different estimators
#2025-11-28
#Written by Tsung-Jen Shen, Chao-Sui Lu, Nai-Tzu Cheng, Youhua Chen
########################################################################################################
########################################################################################################



########################################################################################################
#basic functions
f.to.X = function(f){
  X = c()
  k=1
  while(k<=length(f)){
    if(k>1) X = c(X,rep(k,f[k])) else {
      X = rep(k,f[k])
    }
    k = k+1;
  }
  return(X);
}

########################################################################################################
#basic functions
X.to.f = function(X) {
  f = factor(X, levels = 0:max(X))
  f = table(f, exclude = 0) ## frequency counts
  dimnames(f) = NULL
  return(f)
}

########################################################################################################
#improved estimation of species spatial density parameter
PoiLam <- function(x, zero=FALSE){
  x <- unlist(x)
  # print(x)
  n <- sum(x)  
  f1 <- sum(x==1)
  f2 <- sum(x==2)
  if(zero==FALSE) x <- x[x>0]
  b.solve <- function(b){
    a <- f1 / sum(x*exp(-b*x))
    obj <- sum((x * (1 - a * exp(-b*x)))^2) - sum(x*(x-1)) + 2*f2 
    abs(obj)
  }
  b <- tryCatch(optimize(b.solve, c(0,1))$min, error = function(a) {1})
  a <- f1 / sum(x*exp(-b*x))
  
  return(c(a,b))
}



########################################################################################################
#Chao1 estimator of species richness
Chao1.Pois = function(f) {
  est = sum(f)
  if (length(f) > 1) {
    if (f[2] > 0) {
      A = f[1] / f[2]
    } else {
      if (f[1] > 0) {
        A = f[1] - 1
      } else {
        A = 0
      }
    }
  } else
    A = f[1] - 1
  est = est + A * f[1] / 2
  
  return(est)
}




########################################################################################################
#theoretical true value of estimating rare species richness for a given abundance threshold
#真数值
E_tau_h <- function(h, S, tau_max, lambda) {
  results <- numeric(tau_max)
  for (tau in 1:tau_max) {
    result <- 0
    
    for (i in 1:S) {
      term <- (1 / factorial(tau)) * exp(-h * lambda[i]) * (h * lambda[i])^tau
      result <- result + term
    }
    results[tau] <- result  # 儲存每個 tau 的結果
  }
  return(results)
}

########################################################################################################
#Bayesian weighted estimator for abundance threshold-specific rare endemic species richness estimation
BW.Rh <- function(A, b, h, f,tau_max, alpha_hat, beta_hat) {
  results <- numeric(tau_max)  # 初始化儲存結果的向量
  for (tau in 1:tau_max) {
    result <- 0
    for (n in 1:tau) {
      # 計算公式中的每一項
      term <- f[n] * exp(-(A - b) * n / b * (1 - alpha_hat * exp(-beta_hat * n))) *
        ((A - b) * (n / b) * (1 - alpha_hat * exp(-beta_hat * n)))^(tau - n) /
        factorial(tau - n)
      # 累積結果
      result <- result + term
    }
    # 計算 \hat{R}_tau(h) 的值
    results[tau] <- (h^tau / (A^tau - (A - b)^tau)) * result
  }
  return(results)
}


########################################################################################################
#Unweighted estimator for abundance threshold-specific rare endemic species richness estimation
UW.Rh <- function(A, b, h, f, tau_max, alpha_hat, beta_hat) {
  results <- numeric(tau_max)  # 初始化儲存結果的向量
  for (tau in 1:tau_max) {
    # 計算 UW.Rh 的值
    result <- (h / b)^tau * f[tau] * exp(-(A - b) * (tau / b) * (1 - alpha_hat * exp(-beta_hat * tau)))
    results[tau] <- result  # 儲存每個 tau 的結果
  }
  return(results)
}



########################################################################################################
#Chao-derived estimator for abundance threshold-specific rare species richness estimation
Chao.Eh = function(f=NULL,xi=NULL,h, b,alpha_hat,beta_hat, f0, k.show=3){
  
  if(is.null(f)*is.null(xi)) {
    print("Please input either species frequency counts data or species abundance data!!")
    return(NULL);
  }
  if(is.null(f)) {
    Xi = xi
    f = X.to.f(Xi)
  }
  if(is.null(xi)) {
    # f=Q
    Xi = f.to.X(f)
  }
  
  kmax = length(f)
  
  # cut.pts = max(10,kmax)
  # cut.pts = min(k.show, kmax)
  cut.pts = k.show
  
  # n = sum(f*(1:kmax))
  ## modified f
  # a.f = rep(0,cut.pts)
  ## delta i
  # d = rep(0, length(Xi)) 
  ## estimated fk
  est.Rh = rep(0, cut.pts) 
  # print(d)  
  if(length(f)==1 && f[1]>0) {f=c(f[1],0)}
  if(f[2]==0){
    f1 <- max(f[1] - 1, 0)
    f2 <- 1
  } else {
    f1 = f[1]
    f2 = f[2]
  }
  
  d0 = 0
  if(f0>0) d0 = f1/f0
  
  
  ### detection intensities estimation
  d = Xi/b *(1-alpha_hat *exp(-beta_hat*Xi))
  
  
  for(k in 1:cut.pts){
    est.Rh[k] = sum(1/factorial(k)*exp(-h*d)*(h*d)^k);
    est.Rh[k] = est.Rh[k] + f0*1/factorial(k)*exp(-h*d0)*(h*d0)^k
  }# loop: k
  
  
  return(est.Rh)
}


########################################################################################################
#Horvitz-Thompson estimator for abundance threshold-specific rare species richness estimation
HT.Eh <- function(b, h, d, tau_max) {
  results <- numeric(tau_max)
  for (tau in 1:tau_max) {
    result <- 0  # 初始化結果為 0
    for (i in 1:length(d)) {
      # 計算公式中的各個部分
      term <- (1 / factorial(tau)) * exp(-h * d[i]) * (h * d[i])^tau
      denominator <- 1 - exp(-b * d[i])  # 分母計算
      result <- result + (term / denominator)
    }  
    results[tau] <- result
  }
  return(results)
}


########################################################################################################
#Two-step estimator for abundance threshold-specific rare endemic species richness estimation
#new2(A,b,d,h,3,lambda_hat_star_1)
TS.Eh=function(b,d,h,tau_max,lambda_hat_star_1,Xi){
  results <- numeric(tau_max)
  for (tau in 1:tau_max) {
    seen <- sum((1 / factorial(tau)) * exp(-h * d) * (h * d)^tau)
    unseen <- sum((Xi == 1) * (1 / factorial(tau)) * exp(-h * lambda_hat_star_1) * (h * lambda_hat_star_1)^tau /(b*lambda_hat_star_1))
    result <- unseen+seen
    results[tau] <- result
  }
  return(results)
}



########################################################################################################
calculate_rbias <- function(estimate, true_value) {
  mean(estimate - true_value) / true_value
}
########################################################################################################
calculate_rrmse <- function(estimate, true_value) {
  sqrt(mean((estimate - true_value)^2)) / true_value
}



########################################################################################################
### Generating species frequency counts (fk's) using the bootstrapping method
# f is frequency counts
# b is the estimated parameters (a, b) in dual equations for estimating true detection intensities of species
boot.abundance.Pois.fun = function(S.hat, f, alpha_hat ,beta_hat) {
  D = sum(f)  #已觀測的物種數
  ind = 1:length(f)  #計算已觀測lambda時的Xi
  kmax = length(f)  #最大的已觀測lambda的Xi
  d = rep(0, kmax)  #複製最大frequency count次
  n = sum(ind * f)  #總共物種量
  S.hat = ceiling(S.hat)  #向上取整
  
  ## Estimated lambda of a species with Xi > 0
  boot.lambda = ind /b * (1 - alpha_hat * exp(-beta_hat * ind))
  ## Generate a bootstrap sample
  boot.Xi = unlist(sapply(1:kmax,
                          function(x) {
                            rpois(f[x], b*boot.lambda[x])
                          }))
  
  
  if (S.hat > D) {
    ## Number of unseen species in the local sample
    f0 = S.hat - D
    ## Estimated lambda of a species with Xi = 0
    g0 = f[1] / f0
    
    boot.Xi = c(rpois(f0, g0), boot.Xi)
  }
  
  ## frequency counts
  f.count = X.to.f(boot.Xi)
  
  if(length(f.count) == 1) f.count =  c(f.count[1]-1, 1)
  
  return(f.count)
}#回傳值是一frequency count
#############

########################################################################################################
##### main function for estimating number of rare species with a specific abundance in an additional sample of area size 'h'
## f: species frequency counts data
## xi: species abundance data
Pred.abundance.Pois.rare = function(boot.rep = 50,
                                    A = 1,
                                    f = NULL,
                                    xi = NULL,
                                    h,
                                    k.show = 3
                                    #b.seed = 1234
) {#此處應帶入生成完的Xi
  if (is.null(f) * is.null(xi)) {
    print("Please input either species frequency counts data or species abundance data!!")
    return(NULL)
    
  }
  if (is.null(f)) {
    Xi = xi
    f = X.to.f(Xi)
  }
  if (is.null(xi)) {
    Xi = f.to.X(f)
  }
  # #### Point estimates
  a.p.1 = PoiLam(Xi, zero = FALSE)
  est.R = hat_R_tau_h(A, b, h, f,tau_max, alpha_hat, beta_hat)
  est.f0 = Chao1.Pois(f) - sum(f)
  est.R2 = UW.Rh(A, b, h, f, tau_max, alpha_hat, beta_hat)
  est.R3 = Chao.Eh(f=NULL,xi=NULL,h, b,alpha_hat,beta_hat, f0, k.show=3)  #Chao.Eh(f=NULL,xi=NULL,h, b,alpha_hat,beta_hat, f0, k.show=3)
  est.R4 = HT.Eh(b,h,d,tau_max = 3)  #HT.Eh(b, h, d, tau_max)
  est.R5 = TS.Eh(b,d,h,tau_max=3,lambda_hat_star_1,Xi)
  #### Calculating bootstrap SEs and CIs
  #### Create a space for storing bootstrap samples
  boot.output = NULL
  #set.seed(b.seed)
  
  for (i in 1:boot.rep) {
    b.f = boot.abundance.Pois.fun(S.hat = est.f0 + sum(f),
                                  f = f,
                                  alpha_hat = alpha_hat,
                                  beta_hat = beta_hat)
    
    b.Xi = f.to.X(b.f)
    #####################################
    #below are common for all estimators
    b.a.p.1 = PoiLam(b.Xi, zero = FALSE)
    b.est.f0 = Chao1.Pois(b.f) - sum(b.f)
    b.d = b.Xi/b*(1-b.a.p.1[1]*exp(-b.Xi*b.a.p.1[2]) )
    b.lambda_hat_star_1 <- 1/b*(1- b.a.p.1[1]* exp(-b.a.p.1[2]))
    #####################################
    b.pred.fk.BW = hat_R_tau_h(A, b, h, b.f, tau_max, b.a.p.1[1], b.a.p.1[2])
    b.pred.fk.unweighted = UW.Rh(A, b, h, b.f, tau_max, b.a.p.1[1], b.a.p.1[2])
    b.pred.fk.Chao_derived = new1(b.f,NULL,A,h,b,b.a.p.1[1],b.a.p.1[2],f0 = b.est.f0,k.show = 3)
    b.pred.fk.HT = HT_func(A,b,h,b.d,3)
    b.pred.fk.proposed = new2(A,b,b.d,h,3,b.lambda_hat_star_1)
    boot.output = rbind(boot.output,
                        c(bayesian = b.pred.fk.BW[1:k.show], unweighted = b.pred.fk.unweighted[1:k.show],
                          Chao_derived = b.pred.fk.Chao_derived[1:k.show],HT = b.pred.fk.HT[1:k.show],
                          proposed = b.pred.fk.proposed[1:k.show])
    )
  }### loop: boot.rep
  
  point.est = cbind(bayesian = est.R, unweighted = est.R2, Chao_derived = est.R3, HT = est.R4,proposed = est.R5)
  boot.sd = apply(boot.output, 2, sd, na.rm = T)
  boot.sd = matrix(boot.sd, ncol = 5, byrow = F)
  #boot.ci = apply(boot.output, 2, quantile, probs = c(0.025, 0.975))
  
  Normal.ci = cbind(point.est - 1.96 * boot.sd, point.est + 1.96 * boot.sd)
  for (i in 1:nrow(Normal.ci)) {
    for (j in 1:ncol(Normal.ci)) {
      if (Normal.ci[i, j] < 0)
        Normal.ci[i, j] = 0
    }
  }
  
  output = list()
  output[["Data information"]] = as.matrix(c(A, h),ncol = 1)
  
  rownames(output[["Data information"]]) = c(
    "  Area size of the original sample (A):    ",
    "  Area size of an additional sample (h):   "  )
  colnames(output[["Data information"]]) = c(" ")
  
  output[["Bayesian-weight method"]] = round(cbind(1:k.show, point.est[, 1], boot.sd[, 1], Normal.ci[, c(1, 6)]), 3)
  output[["Unweighted method"]] = round(cbind(1:k.show, point.est[, 2], boot.sd[, 2], Normal.ci[, c(2, 7)]), 3)
  output[["Chao-derived method"]] = round(cbind(1:k.show, point.est[, 3], boot.sd[, 3], Normal.ci[,c(3,8)]),3)
  output[["Horvitz-Thompson method"]] = round(cbind(1:k.show, point.est[, 4], boot.sd[, 4], Normal.ci[,c(4,9)]),3)
  output[["proposed method"]] = round(cbind(1:k.show, point.est[, 5], boot.sd[, 5], Normal.ci[,c(5,10)]),3)
  
  colnames(output[["Bayesian-weight method"]]) =
    colnames(output[["Unweighted method"]]) =
    colnames(output[["Chao-derived method"]]) =
    colnames(output[["Horvitz-Thompson method"]]) =
    colnames(output[["proposed method"]]) =
    c("  k",
      "  Estimate",
      "  Estimated SE",
      "  95% lower limit",
      "  95% upper limit")
  
  output
} ### end of R function "Pred.abundance.Pois.rare"
# 定義參數


########################################################################################################
#for conducting spatially implicit numerical simulations for estimating abundance threshold-specific
#rare species richness with BOOTSTRAPPING
Numerical.Simulation<-function(lambdas,A, S, b_values,h_values,num_runs,Bruns,k.show=3)
{
  # 創建一個空的 data.frame 用來存儲結果
  results_df <- data.frame()
  
  for (b in b_values) {
    #  h_values <- b * multipliers
    for (h in h_values) {
      
      # 初始化儲存估計值的矩陣
      #abundance class=1,2,3, ncol=3!
      BW_estimates_matrix <- matrix(NA, nrow = num_runs, ncol = 3)
      UW_estimates_matrix <- matrix(NA, nrow = num_runs, ncol = 3)
      new1_estimates_matrix <- matrix(NA, nrow = num_runs, ncol = 3)
      HT_estimates_matrix <- matrix(NA, nrow = num_runs, ncol = 3)
      new2_estimates_matrix <- matrix(NA, nrow = num_runs, ncol = 3)
      
      ####################################
      OUTPUT=list()
      for (i in 1:num_runs) {
        Xi = apply(as.matrix(lambdas, ncol=1), 1, function(x) { rpois(1, b * x) })
        f = X.to.f(Xi)
        a.p.1 <- PoiLam(Xi, zero = FALSE)
        alpha_hat <- a.p.1[1]   # 預估參數 alpha
        beta_hat <- a.p.1[2]    # 預估參數 beta
        d = Xi[Xi>0]/b*(1-alpha_hat*exp(-beta_hat*Xi[Xi>0]))
        est.f0 = Chao1.Pois(f)-sum(f)
        lambda_hat_star_1 <- 1/b*(1- alpha_hat* exp(-beta_hat))
        
        # 計算 tau 從 1 到 3 的估計值
        tau_max <- 3
        hat_R_tau_h_values <- BW.Rh(A, b, h,f,tau_max,alpha_hat, beta_hat)
        UW.Rh_values <- UW.Rh(A,b,h,f,tau_max,alpha_hat,beta_hat)
        new1_values <- Chao.Eh(f,NULL,h,b,alpha_hat ,beta_hat,est.f0,3)
        HT.est <- HT.Eh(b,h,d,tau_max)
        new2_values <- TS.Eh(b,d,h,tau_max,lambda_hat_star_1,Xi)
        # 將結果存儲到矩陣的第 i 行
        BW_estimates_matrix[i, ] <- hat_R_tau_h_values
        UW_estimates_matrix[i, ] <- UW.Rh_values
        new1_estimates_matrix[i, ] <- new1_values
        HT_estimates_matrix[i, ] <- HT.est
        new2_estimates_matrix[i, ] <- new2_values
        
        ####################################
        boot.output=NULL
        for (ii in 1:Bruns) {
          b.f = boot.abundance.Pois.fun(S.hat = est.f0 + sum(f),
                                        f = f,
                                        alpha_hat = alpha_hat,
                                        beta_hat = beta_hat)
          
          b.Xi = f.to.X(b.f)
          #####################################
          #below are common for all estimators
          b.a.p.1 = PoiLam(b.Xi, zero = FALSE)
          b.est.f0 = Chao1.Pois(b.f) - sum(b.f)
          b.d = b.Xi/b*(1-b.a.p.1[1]*exp(-b.Xi*b.a.p.1[2]) )
          b.lambda_hat_star_1 <- 1/b*(1- b.a.p.1[1]* exp(-b.a.p.1[2]))
          #####################################
          b.pred.fk.BW = BW.Rh(A, b, h, b.f, tau_max, b.a.p.1[1], b.a.p.1[2])
          b.pred.fk.unweighted = UW.Rh(A, b, h, b.f, tau_max, b.a.p.1[1], b.a.p.1[2])
          b.pred.fk.Chao_derived = Chao.Eh(b.f,NULL,h,b,b.a.p.1[1],b.a.p.1[2],f0 = b.est.f0,k.show = 3)
          #Chao.Eh(f,NULL,h,b,alpha_hat ,beta_hat,est.f0,3)
          b.pred.fk.HT = HT.Eh(b,h,b.d,tau_max)
          b.pred.fk.proposed = TS.Eh(b,b.d,h,tau_max,b.lambda_hat_star_1,b.Xi)
          boot.output = rbind(boot.output,
                              c(bayesian = b.pred.fk.BW[1:k.show], unweighted = b.pred.fk.unweighted[1:k.show],
                                Chao_derived = b.pred.fk.Chao_derived[1:k.show],HT = b.pred.fk.HT[1:k.show],
                                proposed = b.pred.fk.proposed[1:k.show])
          )
        }### loop: boot.rep
        ####################################   
        point.est = cbind(bayesian = hat_R_tau_h_values, unweighted = UW.Rh_values, Chao_derived = new1_values, HT = HT.est,proposed = new2_values)
        boot.sd = apply(boot.output, 2, sd, na.rm = T)
        boot.sd = matrix(boot.sd, ncol = 5, byrow = F)
        #boot.ci = apply(boot.output, 2, quantile, probs = c(0.025, 0.975))
        
        Normal.ci = cbind(point.est - 1.96 * boot.sd, point.est + 1.96 * boot.sd)
        for (k in 1:nrow(Normal.ci)) {
          for (l in 1:ncol(Normal.ci)) {
            if (Normal.ci[k, l] < 0)
              Normal.ci[k, l] = 0
          }
        }
        ################
        output = list()
        output[["Data information"]] = as.matrix(c(A, h),ncol = 1)
        
        rownames(output[["Data information"]]) = c(
          "  Area size of the original sample (A):    ",
          "  Area size of an additional sample (h):   "  )
        colnames(output[["Data information"]]) = c(" ")
        
        output[["Bayesian-weight method"]] = round(cbind(1:k.show, point.est[, 1], boot.sd[, 1], Normal.ci[, c(1, 6)]), 3)
        output[["Unweighted method"]] = round(cbind(1:k.show, point.est[, 2], boot.sd[, 2], Normal.ci[, c(2, 7)]), 3)
        output[["Chao-derived method"]] = round(cbind(1:k.show, point.est[, 3], boot.sd[, 3], Normal.ci[,c(3,8)]),3)
        output[["Horvitz-Thompson method"]] = round(cbind(1:k.show, point.est[, 4], boot.sd[, 4], Normal.ci[,c(4,9)]),3)
        output[["proposed method"]] = round(cbind(1:k.show, point.est[, 5], boot.sd[, 5], Normal.ci[,c(5,10)]),3)
        
        colnames(output[["Bayesian-weight method"]]) =
          colnames(output[["Unweighted method"]]) =
          colnames(output[["Chao-derived method"]]) =
          colnames(output[["Horvitz-Thompson method"]]) =
          colnames(output[["proposed method"]]) =
          c("  k",
            "  Estimate",
            "  Estimated SE",
            "  95% lower limit",
            "  95% upper limit")
        
        #####################################
        OUTPUT[[i]]=output
      }#i numerical simulation replicate
      ####################################   
      # 計算每一行的平均值
      #BSE?! 2026-7-16
      BW_row_means <- apply(BW_estimates_matrix, 2, mean)
      UW_row_means <- apply(UW_estimates_matrix, 2, mean)
      new1_row_means <- apply(new1_estimates_matrix, 2, mean)
      HT_row_means <- apply(HT_estimates_matrix, 2, mean)
      new2_row_means <- apply(new2_estimates_matrix, 2, mean)
      
      # 計算每一行的樣本標準差
      BW_row_sd <- apply(BW_estimates_matrix, 2, sd)
      UW_row_sd <- apply(UW_estimates_matrix, 2, sd)
      new1_row_sd <- apply(new1_estimates_matrix, 2, sd)
      HT_row_sd <- apply(HT_estimates_matrix, 2, sd)
      new2_row_sd <- apply(new2_estimates_matrix, 2, sd)
      
      # 計算 Rbias 和 Rrmse
      R_tau_true <- E_tau_h(h, S, tau_max, lambdas)
      rbias_BW <- sapply(1:length(R_tau_true), function(i) calculate_rbias(BW_estimates_matrix[, i], R_tau_true[i]))
      rrmse_BW <- sapply(1:length(R_tau_true), function(i) calculate_rrmse(BW_estimates_matrix[, i], R_tau_true[i]))
      rbias_unweighted <- sapply(1:length(R_tau_true), function(i) calculate_rbias(UW_estimates_matrix[, i], R_tau_true[i]))
      rrmse_unweighted <- sapply(1:length(R_tau_true), function(i) calculate_rrmse(UW_estimates_matrix[, i], R_tau_true[i]))
      rbias_new1 <- sapply(1:length(R_tau_true), function(i) calculate_rbias(new1_estimates_matrix[, i], R_tau_true[i]))
      rrmse_new1 <- sapply(1:length(R_tau_true), function(i) calculate_rrmse(new1_estimates_matrix[, i], R_tau_true[i]))
      rbias_HT <- sapply(1:length(R_tau_true), function(i) calculate_rbias(HT_estimates_matrix[, i], R_tau_true[i]))
      rrmse_HT <- sapply(1:length(R_tau_true), function(i) calculate_rrmse(HT_estimates_matrix[, i], R_tau_true[i]))
      rbias_new2 <- sapply(1:length(R_tau_true), function(i) calculate_rbias(new2_estimates_matrix[, i], R_tau_true[i]))
      rrmse_new2 <- sapply(1:length(R_tau_true), function(i) calculate_rrmse(new2_estimates_matrix[, i], R_tau_true[i]))
      
      tau_values=c(1:k.show)
      # 將結果添加到 data.frame，避免重複的行
      temp_df <- data.frame(
        b = rep(b, each = length(tau_values)),
        h = rep(h, times = length(tau_values)),
        Tau = tau_values,
        Theoretical_value = round(R_tau_true, 3),
        Avg_Unweighted = round(UW_row_means,3),
        Rbias_Unweighted = round(rbias_unweighted,3),
        Rrmse_Unweighted = round(rrmse_unweighted,3),
        SSE_Unweighted = round(UW_row_sd,3),
        Avg_Bayesian = round(BW_row_means,3),
        Rbias_Bayesian = round(rbias_BW,3),
        Rrmse_Bayesian = round(rrmse_BW,3),
        SSE_Bayesian = round(BW_row_sd,3),
        Avg_new1 = round(new1_row_means,3),
        Rbias_new1 = round(rbias_new1,3),
        Rrmse_new1 = round(rrmse_new1,3),
        SSE_new1 = round(new1_row_sd,3),
        Avg_HT = round(HT_row_means,3),
        Rbias_HT = round(rbias_HT,3),
        Rrmse_HT = round(rrmse_HT,3),
        SSE_HT = round(HT_row_sd,3),
        Avg_new2 = round(new2_row_means,3),
        Rbias_new2 = round(rbias_new2,3),
        Rrmse_new2 = round(rrmse_new2,3),
        SSE_new2 = round(new2_row_sd,3)
      )
      
      # 將當前的 temp_df 合併到 results_df
      #results_df <- rbind(results_df, temp_df)
      ###############################################
      #coverage percent
      Nestimator=5
      CP=matrix(0,nrow=Nestimator,ncol=k.show)
      rownames(CP)=c("BW","UW","Chao","HT","TS")
      for(i in 1:length(tau_values))
      {
        for(j in 1:Nestimator)
        {
          cv=0
          for(k in 1:num_runs)
          {
            if(OUTPUT[[k]][[j+1]][i,4]<=R_tau_true[i] & OUTPUT[[k]][[j+1]][i,5]>=R_tau_true[i])
            {
              cv=cv+1
            }
          }#k
          CP[j,i]=cv/num_runs
        }#j
      }#i
      #############
      temp_df=cbind(temp_df,t(CP))
      results_df=rbind(results_df,temp_df)
      ###############################################
    }#h
  }#b
  ########
  return(results_df)
}#
############








########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
############################A simulated example#########################################################
#a hypothetical sampled biodiversity information vector x
x=c(82,63,48,40,42,37,29,22,15,16,18,18,14,15,17,11,16,3,9,11,11,7,5,8,
    5,5,3,5,5,1,5,3,3,5,4,1,2,2,1,2,2,3,1,3,4,3,1,1,2,1,2,1,1,2,1,2,1,2,1,
    4,2,1,2,1,1,1,1,3,1,1,1,1,1,1,1,2,1,2,1,1,1,2,3,1,1,1,1,2,1,1,1,1,1,1,
    1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
##############################################

##############################################
#numerical simulations for testing the performance of the proposed rare endemic richness estimators
A <- 1       # total region area size, usually set to 1 without explicitly stated and without losing generality
S <- length(x)        # total species richness number based on the sampled biodiversity information vector x

#possible combinations of targeted area size h and sampled area size b for evaluating the performance of estimators
#both values are 0~1
b_values <- c(0.2,0.4,0.6,0.8)
h_values <- c(0.2,0.4,0.6,0.8)
#number of runing times for numerical simulations and estimator peformance check
num_runs <- 5 
Bruns=5
###############
res=Numerical.Simulation(x,A,S, b_values,h_values,num_runs,Bruns)
##############################################


