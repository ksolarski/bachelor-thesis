#Kacper Solarski i6153241
set.seed(6153241)
options(scipen=999)

#PACKAGES####
require("plm")
require("pglm")
require("bife")
require("parallel")
require("base")
require("extraDistr")
require("boot")
require("data.table")
require("xtable")

#GLOBAL VARIABLES####
t <- 8 #time (number of years)
n <- 250 #number of individuals
m <- 10 #Initial number of replications
k <- m #This is the number of succesfull replications (will be updated in monte carlo)
R <- 2 #number of bootstrap replications (for APE RE)
beta_x <- 1 #parameter for x
beta_y <- 1 #parameter for lagged y
beta_vector <- c(beta_y,beta_x)

#Variables for RE
alpha_0 <- -1
alpha_1 <- 1
alpha_3 <- 1
alpha <- c(alpha_0,alpha_1,alpha_3)
sigma <- 1

#FUNCTIONS#### 

#Function which calculates APE for probit model
APE_probit_RE <- function(df,beta_vector,alpha,sigma,dynamic=TRUE){
  sigma_sq <- sigma^2
  x <- as.matrix(subset(df, select = c(x))) #selecting our regressors
  temp <- (1+sigma_sq)^(-1/2) #we multiply all the parameters with this value
  if (dynamic==TRUE){
    alpha <- alpha[1]*temp + alpha[2]*df$y_i0*temp + df$x_bar*alpha[3]*temp
  }
  else alpha <- alpha[1]*temp + df$x_bar*alpha[3]*temp
  beta_vector <- beta_vector*temp
  beta_y <- beta_vector[1]
  beta_x <- beta_vector[2]
  APE <- vector("numeric",length=2)
  APE[1] <- mean(pnorm(beta_y + x %*% beta_x + alpha)-pnorm(x %*% beta_x + alpha))
  APE[2] <- mean(dnorm(beta_y + x %*% beta_x + alpha)*beta_x)
  return(APE)
}

#True APE for Probit model
true_APE_probit_RE <- function(df){
  x <- as.matrix(subset(df, select = c(x))) #selecting our regressor
  c <- as.matrix(subset(df, select = c(c))) #c is the individual effect
  APE <- vector("numeric",length=2)
  APE[1] <- mean(pnorm(beta_y + x %*% beta_x + c)-pnorm(x %*% beta_x + c))
  APE[2] <- mean(dnorm(beta_y + x %*% beta_x + c)*beta_x)
  return(APE)
}

#Function which checks rejections (t-test)
rejection_t_test <- function(coeff,vector=beta_vector,df=(n=1)*(t-1),dynamic=TRUE){
  rejection_matrix <- matrix(NA,ncol=2,nrow=2)
  colnames(rejection_matrix)=c("0.05","0.1")
  rownames(rejection_matrix)=c("y_lag","x")
  if (dynamic==TRUE){
    rejection_matrix["y_lag","0.05"] <- (abs((coeff["y_lag",1]-vector[1])/coeff["y_lag",2])>qt(p=0.975, df=df))
    rejection_matrix["y_lag","0.1"] <- (abs((coeff["y_lag",1]-vector[1])/coeff["y_lag",2])>qt(p=0.95, df=df))
  }
  rejection_matrix["x","0.05"] <- (abs((coeff["x",1]-vector[2])/coeff["x",2])>qt(p=0.975, df=df))
  rejection_matrix["x","0.1"] <- (abs((coeff["x",1]-vector[2])/coeff["x",2])>qt(p=0.95, df=df))
  return(rejection_matrix)
}

#Function which checks rejections (Z-test)
rejection_z_test <- function(coeff,vector=beta_vector,dynamic=TRUE){
  rejection_matrix <- matrix(NA,ncol=2,nrow=2)
  colnames(rejection_matrix)=c("0.05","0.1")
  rownames(rejection_matrix)=c("y_lag","x")
  if (dynamic==TRUE){
    rejection_matrix["y_lag","0.05"] <- (abs((coeff["y_lag",1]-vector[1])/coeff["y_lag",2])>qnorm(p=0.975, mean=0,sd=1))
    rejection_matrix["y_lag","0.1"] <- (abs((coeff["y_lag",1]-vector[1])/coeff["y_lag",2])>qnorm(p=0.95, mean=0,sd=1))
  }
  rejection_matrix["x","0.05"] <- (abs((coeff["x",1]-vector[2])/coeff["x",2])>qnorm(p=0.975, mean=0,sd=1))
  rejection_matrix["x","0.1"] <- (abs((coeff["x",1]-vector[2])/coeff["x",2])>qnorm(p=0.95, mean=0,sd=1))
  return(rejection_matrix)
}

#Function which checks rejection for RE APE using bootstrap
rejection_APE_RE_Bootstrap <- function(data, true_APE_vector,df=(n=1)*(t-1),dynamic=TRUE){
  rejection_matrix <- matrix(NA,ncol=2,nrow=2)
  colnames(rejection_matrix)=c("0.05","0.1")
  rownames(rejection_matrix)=c("y_lag","x")
  boot_est <- boot(data=data,statistic=APE_bootstrap,R=R,parallel = "multicore",ncpus = 2,dynamic=dynamic)
  if (dynamic==TRUE){
    rejection_matrix["y_lag","0.05"] <- (abs((boot_est$t0[1]-true_APE_vector[1])/sd(boot_est$t[,1])) > qt(p=0.975, df=df))
    rejection_matrix["y_lag","0.1"] <- (abs((boot_est$t0[1]-true_APE_vector[1])/sd(boot_est$t[,1])) > qt(p=0.95, df=df))
  }
  rejection_matrix["x","0.05"] <- (abs((boot_est$t0[2]-true_APE_vector[2])/sd(boot_est$t[,2])) > qt(p=0.975, df=df))
  rejection_matrix["x","0.1"] <- (abs((boot_est$t0[2]-true_APE_vector[2])/sd(boot_est$t[,2])) > qt(p=0.95, df=df))
  return(rejection_matrix)
}

#Bootstrap function for APE RE
APE_bootstrap <- function(df,i,dynamic){
  data <- df[i,]
  if (dynamic ==TRUE){
    op_RE <- pglm(y ~ y_lag + x + x_bar + y_i0,data = data,family = binomial(link='probit'),
                  method = 'nr', effect = "individual", model = "random")
  }
  else (op_RE <- pglm(y ~ x + x_bar,data = data,family = binomial(link='probit'),
                      method = 'nr', effect = "individual", model = "random"))
  #Getting latent variables RE
  beta_y <- as.numeric(coef(op_RE)["y_lag"])
  beta_x <- as.numeric(coef(op_RE)["x"])
  
  #Creating a vector of coefficients
  beta_vector <- c(beta_y,beta_x)
  
  #Replacing NAs with 0: needed for static model
  beta_vector[is.na(beta_vector)] <- 0
  
  #Getting estimated alphas and sigma
  alpha_0 <- as.numeric(coef(op_RE)["(Intercept)"])
  alpha_1 <- as.numeric(coef(op_RE)["y_i0"])
  alpha_3 <- as.numeric(coef(op_RE)["x_bar"])
  sigma <- as.numeric(coef(op_RE)["sigma"])
  
  #Putting alphas in a vector
  alpha <- c(alpha_0,alpha_1,alpha_3)
  
  #Replacing NAs with 0: needed for static model
  alpha[is.na(alpha)] <- 0
  
  APE <- APE_probit_RE(data, beta_vector,alpha,sigma,dynamic)
  return(APE)
}



#Function which deletes all the observations where y doesn't change (consequently I can apply LPM and RE on reduced dataset)
reduce_data <- function(df){
  new_df <- data.frame(matrix(ncol = ncol(df), nrow = 0)) #Creating new matrix with the same number of columns
  names <- names(df)
  colnames(new_df) <- names
  y <- df$y  #Vetor y
  for (i in 1:n){
    one_individual <- y[(i*t-(t-1)):(i*t)] #vector of y observations for one individual
    if (nrow(as.matrix(unique(one_individual))) > 1){ #or ==2
      temp <- df[(i*t-(t-1)):(i*t),]
      new_df <- rbind(new_df, temp)
    }
  }
  new_df <- as.data.frame(new_df)
  return(new_df)
}

#Function which puts estimated alphas into matrix (not used)
alpha_into_matrix <- function(est){
  alphas <- as.vector(unlist(est["alpha"]))
  alphas_matrix <- matrix(nrow=0,ncol=1)
  for (i in 1:(length(alphas))){
    temp <- as.matrix(rep(alphas[i],t))
    alphas_matrix <- rbind(alphas_matrix,temp)
  }
  return(alphas_matrix)
}
#STATIC PROBIT MODEL #####

#Simulating X vector
simulate_x_static_probit <- function(){
  x <- matrix(0, nrow = t, ncol=1)
  x_0 <- runif(1, -0.5, 0.5) #Initialize the first value (X at time 0 which we don't include in matrix)
  x[1] <- 0.1 + x_0/2 + runif(1, -0.5, 0.5) #That's the first value which is relevant for us
  for (i in 2:t){ #We need a for loop since X at time t depends on X at time t-1
    x[i] <- i/10 + (x[i-1])/2 + runif(1, -0.5, 0.5)
  }
  return(list(x,x_0))
}

#Simulating the individual effects from N(0,sigma)
simulate_effects_static_probit <- function(x){
  alpha_i <- rnorm(1,0,sigma)
  c_i <- alpha_0 + mean(x)*alpha_3+alpha_i
  indiv_effects <- matrix(rep(c_i,t), nrow = t, ncol = 1) #putting in it matrix
  return(indiv_effects)
}

#Simulating the disturbancess from N(0,1)
simulate_errors_normal <- function(){
  errors <- matrix(0,nrow = t,ncol = 1)
  errors[] <- rnorm(t,0,1)
  return(errors)
}

#Generating observed variable (Y)
#I use a for loop here - some package might make this part faster
generate_y_static_probit <- function(x_list,c,u){
  y <- matrix(0, nrow = t, ncol = 1)
  y_latent <- y <- matrix(0, nrow = t, ncol = 1)
  x <- x_list[[1]]
  x_0 <- x_list[[2]]
  y_0 <- x_0*beta_x + c[1] + runif(1,-0.5,0.5)
  if (y_0>0){
    y_0 <-1
  }
  else y_0 <-0
  for (i in 1:t){
    y_latent[i] <- x[i]*beta_x + c[i] + u[i]
    if (x[i]*beta_x + c[i] + u[i] > 0){
      y[i] <- 1
    }
  }
  return(list(y,y_latent))
}

#Function which generates data for static probit model and puts in into dataframe
generate_static_probit <- function(){
  df <- data.frame(ID=numeric(),Time=numeric(), y=numeric(),y_latent=numeric(), x=numeric(), x_bar=numeric(), c=numeric(), u=numeric())
  for (i in 1:n){
    x_list <- simulate_x_static_probit()
    x <- x_list[[1]]
    x_0 <-  x_list[[2]]
    c <- simulate_effects_static_probit(x)
    u <- simulate_errors_normal()
    y_list <- generate_y_static_probit(x_list,c,u)
    y <- y_list[[1]]
    y_latent <- y_list[[2]]
    data <- data.frame(ID = rep(i,t),Time = c(1:t), y=y, y_latent=y_latent, x=x, x_bar=mean(x), c=c, u=u)
    df <- rbind(df, data)
  }
  df <- pdata.frame(df, index=c("ID","Time"))
  return(df)
}

#DYNAMIC LOGIT MODEL (NOT RELEVANT) #####

#Simulating X vector
simulate_x_dynamic_logit <- function(){
  x <- matrix(0, nrow = t, ncol = 1)
  x[] <- rnorm(t,0,pi^2/3)
  return(x)
}

#Simulating the individual effects which are correlated with regressor
simulate_effects_dynamic_logit <- function(x){
  temp <- sum(x[(1:4)])/4 # individual effect is an average of the first 4 values of X
  indiv_effects <- matrix(rep(temp,t), nrow = t, ncol = 1)
  return(indiv_effects)
}

#Simulating the disturbances from Logistic(0,pi^2/3)
simulate_errors_logistic <- function(){
  errors <- matrix(0,nrow = t,ncol = 1)
  errors[] <- rlogis(t,0, pi^2/3)
  return(errors)
}

#Generating observed variable (Y)
generate_y_dynamic_logit <- function(x,c,u){
  y <- matrix(0, nrow = t, ncol = 1) #matrix of variable y_t
  y_lag <- matrix(0, nrow = t, ncol = 1) #matrix of variable y_t-1 (which we will use as a regressor)
  y_0 <- rnorm(1,0,pi^2/3)*beta_x + c[1] - rlogis(1,0, pi^2/3) #The latent starting value (y_0)
  if (y_0 > 0){ #Converting latent starting value to observed starting value
    y_0 <- 1
  }
  else y_0 <- 0
  temp <- y_0 #this will be our lagged variable which will keep being updated
  for (i in 1:t){
    y_lag[i] <- temp
    latent_y <- y_lag[i]*beta_y + x[i]*beta_x + c[i] - u[i]
    if (latent_y > 0){
      y[i] <- 1
    }
    temp <- y[i] #we update temp because this will be used as a regressor in the next interation
  }
  combined_matrix <-  list(y,y_lag,y_0) #we need to combine the two matrices
  return(combined_matrix)
}

#Function which generates data for dynamic logit model and puts in into dataframe
generate_dynamic_logit <- function(){
  df <- data.frame(ID=numeric(), Time=numeric(), y=numeric(), y_lag=numeric(),
                   x=numeric(), y_i0=numeric(),x_bar=numeric(), c=numeric())
  for (i in 1:n){
    x <- simulate_x_dynamic_logit()
    u <- simulate_errors_logistic()
    c <- simulate_effects_dynamic_logit(x)
    y_list <- generate_y_dynamic_logit(x,c,u)
    y <- y_list[[1]]
    y_lag <- y_list[[2]]
    y_i0 <- y_list[[3]]
    temp <- data.frame(Time = c(1:t), ID = rep(i,t), y=y, y_lag=y_lag, x=x, x_bar=mean(x), y_i0=y_i0, c=c)
    df <- rbind(df, temp)
  }
  df <- pdata.frame(df, index=c("ID","Time"))
  return(df)
}

#DYNAMIC PROBIT MODEL#####

#Function for errors is already in static probit model
simulate_x_dynamic_probit_random <- function(){
  #set.seed(6153241)
  x <- matrix(0, nrow = t, ncol = 1)
  x[] <- rnorm(t,0,1)
  #rm(.Random.seed, envir=.GlobalEnv)
  return(x)
}

simulate_effects_dynamic_probit_random <- function(y_0,x,i){
  #set.seed(i)
  alpha_i <- rnorm(1,0,sigma) #Change here for misspecification
  c_i <- alpha_0+y_0*alpha_1+mean(x)*alpha_3+alpha_i #Or here #+ alpha_1*y_0
  indiv_effects <- matrix(rep(c_i,t), nrow = t, ncol = 1) #putting in it matrix
  #rm(.Random.seed, envir=.GlobalEnv)
  return(indiv_effects)
}

generate_y_dynamic_probit_random <- function(x,c,u,y_0){
  y <- matrix(0, nrow = t, ncol = 1) #matrix of variable y_t
  y_lag <- matrix(0, nrow = t, ncol = 1) #matrix of variable y_t-1 (which we will use as a regressor)
  temp <- y_0 #this will be our lagged variable which will keep being updated
  for (i in 1:t){
    y_lag[i] <- temp
    latent_y <- y_lag[i]*beta_y + x[i]*beta_x + c[i] + u[i]
    if (latent_y > 0){
      y[i] <- 1
    }
    temp <- y[i] #we update temp because this will be used as a regressor in the next interation
  }
  combined_matrix <-  matrix(c(y,y_lag), nrow = t, ncol = 2) #we need to combine the two matrices
  return(combined_matrix)
}

intomatrix <- function(y_0){ #Not used
  y_i0 <- as.matrix(rep(y_0,t))
  return(y_i0)
}

generate_dynamic_probit_random <- function(){
  df <- data.frame(ID=numeric(), Time=numeric(), y=numeric(), y_lag=numeric(),
                   x=numeric(), y_i0=numeric(),x_bar=numeric(), c=numeric())
  for (i in 1:n){
    x <- simulate_x_dynamic_probit_random()
    u <- simulate_errors_normal()
    y_0 <- rnorm(1,0,sigma) #Simulating latent variable for y_i0
    if (y_0 > 0){ #Converting latent starting value to observed starting value
      y_0 <- 1
    }
    else y_0 <- 0
    c <- simulate_effects_dynamic_probit_random(y_0,x,i)
    matrix_y <- generate_y_dynamic_probit_random(x,c,u,y_0)
    y <- matrix_y[,1]
    y_lag <- matrix_y[,2]
    temp <- data.frame(ID = rep(i,t),Time = c(1:t), y=y, y_lag=y_lag, x=x, y_i0=y_0,x_bar=mean(x), c=c)
    df <- rbind(df, temp)
  }
  df <- pdata.frame(df, index=c("ID","Time"))
  return(df)
}

#MONTE CARLO DYNAMIC PROBIT####
#I cannot use tryCatch outside the whole function because then bootstrap doesn't work
run_dynamic_probit <- function(){
  names <- c("RE - Beta_y","RE(R) - Beta_y","FE - Beta_y","FE(BC) - Beta_y",
             "RE - Beta_x","RE(R) - Beta_x","FE - Beta_x","FE(BC) - Beta_x",
             "RE - APE(y)","RE(R) - APE(y)","FE - APE(y)","FE(BC) - APE(y)","LPM - APE(y)","LPM(R) - APE(y)",
             "RE - APE(x)","RE(R) - APE(x)","FE - APE(x)","FE(BC) - APE(x)","LPM - APE(x)","LPM(R) - APE(x)")
  
  # Running Monte Carlo
  montecarlo_dynamic_RE <- function(data){
    estimations <- rep(NA,length(names)) #rows - replications; columns - different methods
    names(estimations) <- names
    rejections_0.05 <- rep(NA,length(names)) #rows - replications; columns - different methods
    names(rejections_0.05) <- names
    rejections_0.1 <- rep(NA,length(names)) #rows - replications; columns - different methods
    names(rejections_0.1) <- names
    
    #True APE
    true_APE_vector <- true_APE_probit_RE(data)
    true_APE_y <- true_APE_vector[1]
    true_APE_x <- true_APE_vector[2]
    
    #RE####
    #RE coefficients
    tryCatch(
      expr = {
        op_RE <- pglm(y ~ y_lag + x + x_bar + y_i0,data = data,family = binomial(link='probit'),
                      method = 'nr', effect = "individual", model = "random")
      },
      error = function(e){
        return(NA) #We return NA and then we'll delete them
      },
      warning = function(w){
      },
      finally = {}
    )    
    
    #Getting latent variables RE
    beta_y_hat <- estimations["RE - Beta_y"] <- as.numeric(coef(op_RE)["y_lag"])
    beta_x_hat <- estimations["RE - Beta_x"] <- as.numeric(coef(op_RE)["x"])
    
    #Creating a vector of coefficients
    beta_vector_hat <- c(beta_y_hat,beta_x_hat)
    
    #Getting estimated alphas and sigma
    alpha_0_hat <- as.numeric(coef(op_RE)["(Intercept)"])
    alpha_1_hat <- as.numeric(coef(op_RE)["y_i0"])
    alpha_3_hat <- as.numeric(coef(op_RE)["x_bar"])
    sigma_hat <- as.numeric(coef(op_RE)["sigma"])  
    
    #Putting alphas in a vector
    alpha_hat <- c(alpha_0_hat,alpha_1_hat,alpha_3_hat)
    
    #Getting APE for RE
    APE_hat <- APE_probit_RE(data,beta_vector_hat,alpha_hat,sigma_hat)
    estimations["RE - APE(y)"] <- APE_hat[1]/true_APE_y
    estimations["RE - APE(x)"] <- APE_hat[2]/true_APE_x
    
    #Checking rejection for latent RE
    coeff <- summary(op_RE)$estimate
    rejection_matrix <- rejection_t_test(coeff)
    rejections_0.05["RE - Beta_y"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["RE - Beta_y"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["RE - Beta_x"] <- rejection_matrix["x","0.05"]
    rejections_0.1["RE - Beta_x"] <- rejection_matrix["x","0.1"]
    
    #Checking rejection for APE RE (Bootstrap)
    #rejection_matrix <- rejection_APE_RE_Bootstrap(data,true_APE_vector)
    #rejections_0.05["RE - APE(y)"] <- rejection_matrix["y_lag","0.05"]
    #rejections_0.1["RE - APE(y)"] <- rejection_matrix["y_lag","0.1"]
    #rejections_0.05["RE - APE(x)"] <- rejection_matrix["x","0.05"]
    #rejections_0.1["RE - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #FE####
    op_FE <- bife(y ~ y_lag + x | ID,
                  data = data, model = "probit")
    
    #FE coefficients
    estimations["FE - Beta_y"] <- as.numeric(coef(op_FE)["y_lag"])
    estimations["FE - Beta_x"] <- as.numeric(coef(op_FE)["x"])
    
    #FE - APE
    estimations["FE - APE(y)"] <- get_APEs(op_FE)$delta["y_lag"]/true_APE_y
    estimations["FE - APE(x)"] <- get_APEs(op_FE)$delta["x"]/true_APE_x
    
    #Checking rejection for latent FE
    coeff <- summary(op_FE)$cm
    rejection_matrix <- rejection_z_test(coeff)
    rejections_0.05["FE - Beta_y"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["FE - Beta_y"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["FE - Beta_x"] <- rejection_matrix["x","0.05"]
    rejections_0.1["FE - Beta_x"] <- rejection_matrix["x","0.1"]
    
    #Checking rejection for APE FE
    coeff <- summary(get_APEs(op_FE))
    rejection_matrix <- rejection_z_test(coeff,vector=true_APE_vector)
    rejections_0.05["FE - APE(y)"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["FE - APE(y)"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["FE - APE(x)"] <- rejection_matrix["x","0.05"]
    rejections_0.1["FE - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #FE(BC)####
    tryCatch(
      expr = {
        op_FE <- bias_corr(op_FE, L = 1L)
      },
      error = function(e){
        return(NA) #We return NA and then we'll delete them
      },
      warning = function(w){
      },
      finally = {}
    )    
    
    #FE(BC) coefficients
    estimations["FE(BC) - Beta_y"] <- as.numeric(coef(op_FE)["y_lag"])
    estimations["FE(BC) - Beta_x"] <- as.numeric(coef(op_FE)["x"])
    
    #FE(BC) - APe
    estimations["FE(BC) - APE(y)"] <- get_APEs(op_FE)$delta["y_lag"]/true_APE_y
    estimations["FE(BC) - APE(x)"] <- get_APEs(op_FE)$delta["x"]/true_APE_x
    
    #Checking rejection for latent FE(BC)
    coeff <- summary(op_FE)$cm
    rejection_matrix <- rejection_z_test(coeff)
    rejections_0.05["FE(BC) - Beta_y"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["FE(BC) - Beta_y"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["FE(BC) - Beta_x"] <- rejection_matrix["x","0.05"]
    rejections_0.1["FE(BC) - Beta_x"] <- rejection_matrix["x","0.1"]
    
    #Checking rejection for APE FE(BC)
    coeff <- summary(get_APEs(op_FE))
    rejection_matrix <- rejection_z_test(coeff,vector=true_APE_vector)
    rejections_0.05["FE(BC) - APE(y)"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["FE(BC) - APE(y)"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["FE(BC) - APE(x)"] <- rejection_matrix["x","0.05"]
    rejections_0.1["FE(BC) - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #LPM####
    #LPM with all the observations
    op_LPM <- plm(y ~ y_lag + x,data=data,effect="individual",model="within")
    
    #Coefficients
    estimations["LPM - APE(y)"] <- coef(op_LPM)["y_lag"]/true_APE_y
    estimations["LPM - APE(x)"] <- coef(op_LPM)["x"]/true_APE_x
    
    #Checking rejections for LPM
    coeff <- summary(op_LPM)$coefficients
    df <-summary(op_LPM)$df[2]
    rejection_matrix <- rejection_t_test(coeff=coeff,vector=true_APE_vector,df=df)
    rejections_0.05["LPM - APE(y)"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["LPM - APE(y)"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["LPM - APE(x)"] <- rejection_matrix["x","0.05"]
    rejections_0.1["LPM - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #Reducing the dataset to obervations varying over time
    data <- reduce_data(data)
    
    #RE with reduces dataset
    tryCatch(
      expr = {
        op_RE <- pglm(y ~ y_lag + x + x_bar + y_i0,data = data,family = binomial(link='probit'),
                      method = 'nr', effect = "individual", model = "random")
      },
      error = function(e){
        return(NA) #We return NA and then we'll delete them
      },
      warning = function(w){
      },
      finally = {}
    )    
    
    #Getting latent variables RE
    beta_y_hat <- estimations["RE(R) - Beta_y"] <- as.numeric(coef(op_RE)["y_lag"])
    beta_x_hat <- estimations["RE(R) - Beta_x"] <- as.numeric(coef(op_RE)["x"])
    
    #Creating a vector of coefficients
    beta_vector_hat <- c(beta_y_hat,beta_x_hat)
    
    #Getting estimated alphas and sigma
    alpha_0_hat <- as.numeric(coef(op_RE)["(Intercept)"])
    alpha_1_hat <- as.numeric(coef(op_RE)["y_i0"])
    alpha_3_hat <- as.numeric(coef(op_RE)["x_bar"])
    sigma_hat <- as.numeric(coef(op_RE)["sigma"])  
    
    #Putting alphas in a vector
    alpha_hat <- c(alpha_0_hat,alpha_1_hat,alpha_3_hat)
    
    #Getting APE for RE
    APE_hat <- APE_probit_RE(data,beta_vector_hat,alpha_hat,sigma_hat)
    estimations["RE(R) - APE(y)"] <- APE_hat[1]/true_APE_y
    estimations["RE(R) - APE(x)"] <- APE_hat[2]/true_APE_x
    
    #Checking rejection for latent RE
    coeff <- summary(op_RE)$estimate
    rejection_matrix <- rejection_t_test(coeff,df=(nrow(data)/t-1)*(t-1))
    rejections_0.05["RE(R) - Beta_y"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["RE(R) - Beta_y"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["RE(R) - Beta_x"] <- rejection_matrix["x","0.05"]
    rejections_0.1["RE(R) - Beta_x"] <- rejection_matrix["x","0.1"]
    
    #Checking rejection for APE RE (Bootstrap)
    #rejection_matrix <- rejection_APE_RE_Bootstrap(data,true_APE_vector)
    #rejections_0.05["RE(R) - APE(y)"] <- rejection_matrix["y_lag","0.05"]
    #rejections_0.1["RE(R) - APE(y)"] <- rejection_matrix["y_lag","0.1"]
    #rejections_0.05["RE(R) - APE(x)"] <- rejection_matrix["x","0.05"]
    #rejections_0.1["RE(R) - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #LPM with reduced dataset
    op_LPM <- plm(y ~ y_lag + x,data=data,effect="individual",model="within")
    
    #Coefficients
    estimations["LPM(R) - APE(y)"] <- coef(op_LPM)["y_lag"]/true_APE_y
    estimations["LPM(R) - APE(x)"] <- coef(op_LPM)["x"]/true_APE_x
    
    #Checking rejections for reduced LPM
    coeff <- summary(op_LPM)$coefficients
    df <-summary(op_LPM)$df[2]
    rejection_matrix <- rejection_t_test(coeff=coeff,vector=true_APE_vector,df=df)
    rejections_0.05["LPM(R) - APE(y)"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["LPM(R) - APE(y)"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["LPM(R) - APE(x)"] <- rejection_matrix["x","0.05"]
    rejections_0.1["LPM(R) - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #Returning esitmations and rejection data
    return(c(estimations,rejections_0.05,rejections_0.1))
  }
  #Generating list of dataframes
  df_list <- list()
  for (i in 1:m){
    df_list <- c(list(generate_dynamic_probit_random()),df_list)
  }
  
  #Applying estimations to each dataframe
  estimations_rejections <-mclapply(df_list, FUN = montecarlo_dynamic_RE, mc.cores=2)
  
  #Transforming the results
  t_matrix <- matrix(NA, ncol = length(names)*3, nrow = 0)
  for (i in 1:length(estimations_rejections)){
    t_matrix <- rbind(t(unlist(estimations_rejections[i])),t_matrix)
  }
  
  #Splitting the results
  estimations <- as.matrix(t_matrix[,1:20])
  estimations_p_0.05 <- as.matrix(t_matrix[,21:40])
  estimations_p_0.1 <- as.matrix(t_matrix[,41:60])
  
  #Deleting NA rows from estimations matrices
  estimations <- na.omit(estimations)
  
  #Updating the number of succesful replications
  k <- nrow(estimations)
  
  #Here we obtain results of interests from estimation matrix
  obtain_results <- function (estimations){
    results_dynamic_probit <- matrix(0, ncol=ncol(estimations), nrow=7, 
                                     dimnames = list(c("Mean Bias","Median Bias", "MSE", "MAE", "SD", "p_0.05","p_0.1"), names))
    
    #Filling in mean,SD,p0.05,p0.1 (this is the same for all the esitmators)
    for (i in 1:ncol(estimations)){
      results_dynamic_probit["SD",i] <- sd(estimations[,i])
      results_dynamic_probit["p_0.05",i] <- mean(estimations_p_0.05[,i],na.rm=TRUE)
      results_dynamic_probit["p_0.1",i] <- mean(estimations_p_0.1[,i],na.rm=TRUE)
    }
    
    #MSE and MAE for estimators of beta_y
    for (i in 1:3){
      results_dynamic_probit["Mean Bias",i] <- mean(estimations[,i]-beta_y)
      results_dynamic_probit["Median Bias",i] <- median(estimations[,i]-beta_y)
      results_dynamic_probit["MSE",i] <- mean((estimations[,i]-beta_y)^2) 
      results_dynamic_probit["MAE",i] <- mean(abs(estimations[,i]-beta_y))
    }
    
    #MSE and MAE for estimators of beta_x
    for (i in 4:6){
      results_dynamic_probit["Mean Bias",i] <- mean(estimations[,i]-beta_x)
      results_dynamic_probit["Median Bias",i] <- median(estimations[,i]-beta_x)
      results_dynamic_probit["MSE",i] <- mean((estimations[,i]-beta_x)^2)
      results_dynamic_probit["MAE",i] <- mean(abs(estimations[,i]-beta_x))
    }
    
    #MSE and MAE for estimators of APE (1 is the true value since it's the ratio)
    for (i in 7:ncol(estimations)){
      results_dynamic_probit["Mean Bias",i] <- mean(estimations[,i]-1)
      results_dynamic_probit["Median Bias",i] <- median(estimations[,i]-1)
      results_dynamic_probit["MSE",i] <- mean((estimations[,i]-1)^2)
      results_dynamic_probit["MAE",i] <- mean(abs(estimations[,i]-1))
    }
    
    return(t(results_dynamic_probit)) #Transpose is easier to read for me
  }
  
  results_dynamic_probit <- obtain_results(estimations)
  
  #Printing results
  printing <- function(){
    closeAllConnections()
    sink("Results_Dynamic_Probit.txt")
    cat("Number of indivuduals (n) =",n)
    cat(sep="\n")
    cat("Time (t) =",t)
    cat(sep="\n")
    cat("Number of initial iterations =",m)
    cat(sep="\n")
    cat("Number of actual iterations =",k)
    cat(sep="\n")
    cat("Success rare of simulation =",k/m)
    cat(sep="\n")
    cat("Bootstrap Replications =",R)
    cat(sep="\n")
    cat("Beta = ",beta_vector)
    cat(sep="\n")
    cat("Alpha = ", alpha)
    cat(sep="\n")
    cat("Sigma = ", sigma)
    cat(sep="\n")
    print(results_dynamic_probit,sep="\n")
    sink()
  } 
  printing()
  return(results_dynamic_probit)
}
#Execution Dynamic Probit
results_dynamic_probit <- run_dynamic_probit()
#MONTE CARLO STATIC PROBIT####

beta_y <- beta_vector[1] <- 0 #This way we can still use the functions for dynamic probit model

run_static_probit <- function(){
  names <- c("RE - Beta_x","RE(R) - Beta_x","FE - Beta_x","FE(BC) - Beta_x",
             "RE - APE(x)","RE(R) - APE(x)","FE - APE(x)","FE(BC) - APE(x)","LPM - APE(x)","LPM(R) - APE(x)")
  
  # Running Monte Carlo
  montecarlo_static_RE <- function(data){
    estimations <- rep(NA,length(names)) #rows - replications; columns - different methods
    names(estimations) <- names
    rejections_0.05 <- rep(NA,length(names)) #rows - replications; columns - different methods
    names(rejections_0.05) <- names
    rejections_0.1 <- rep(NA,length(names)) #rows - replications; columns - different methods
    names(rejections_0.1) <- names
    
    #True APE
    true_APE_vector <- true_APE_probit_RE(data)
    true_APE_x <- true_APE_vector[2]
    
    #RE####
    #RE coefficients
    tryCatch(
      expr = {
        op_RE <- pglm(y ~ x + x_bar,data = data,family = binomial(link='probit'),
                      method = 'nr', effect = "individual", model = "random")
      },
      error = function(e){
        return(NA) #We return NA and then we'll delete them
      },
      warning = function(w){
      },
      finally = {}
    )    
    
    #Getting latent variables RE
    beta_x_hat <- estimations["RE - Beta_x"] <- as.numeric(coef(op_RE)["x"])
    beta_vector_hat <- c(0,beta_x_hat)
    
    #Getting estimated alphas and sigma
    alpha_0_hat <- as.numeric(coef(op_RE)["(Intercept)"])
    alpha_3_hat <- as.numeric(coef(op_RE)["x_bar"])
    sigma_hat <- as.numeric(coef(op_RE)["sigma"])
    
    #Putting alphas in a vector
    alpha_hat <- c(alpha_0_hat,0,alpha_3_hat)
    
    #Getting APE for RE
    APE_hat <- APE_probit_RE(data,beta_vector_hat,alpha_hat,sigma_hat,dynamic=FALSE)
    estimations["RE - APE(x)"] <- APE_hat[2]/true_APE_x
    
    #Checking rejection for latent RE
    coeff <- summary(op_RE)$estimate
    rejection_matrix <- rejection_t_test(coeff,dynamic=FALSE)
    rejections_0.05["RE - Beta_x"] <- rejection_matrix["x","0.05"]
    rejections_0.1["RE - Beta_x"] <- rejection_matrix["x","0.1"]
    
    #Checking rejection for APE RE (Bootstrap)
    #rejection_matrix <- rejection_APE_RE_Bootstrap(data,true_APE_vector,dynamic=FALSE)
    #rejections_0.05["RE - APE(x)"] <- rejection_matrix["x","0.05"]
    #rejections_0.1["RE - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #FE####
    op_FE <- bife(y ~ x | ID,
                  data = data, model = "probit")
    
    #FE coefficients
    estimations["FE - Beta_x"] <- as.numeric(coef(op_FE)["x"])
    
    #FE - APE
    estimations["FE - APE(x)"] <- get_APEs(op_FE)$delta["x"]/true_APE_x
    
    #Checking rejection for latent FE
    coeff <- summary(op_FE)$cm
    rejection_matrix <- rejection_z_test(coeff,dynamic=FALSE)
    rejections_0.05["FE - Beta_x"] <- rejection_matrix["x","0.05"]
    rejections_0.1["FE - Beta_x"] <- rejection_matrix["x","0.1"]
    
    #Checking rejection for APE FE
    coeff <- summary(get_APEs(op_FE))
    rejection_matrix <- rejection_z_test(coeff,vector=true_APE_vector,dynamic=FALSE)
    rejections_0.05["FE - APE(x)"] <- rejection_matrix["x","0.05"]
    rejections_0.1["FE - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #FE(BC)####
    tryCatch(
      expr = {
        op_FE <- bias_corr(op_FE)
      },
      error = function(e){
        return(NA) #We return NA and then we'll delete them
      },
      warning = function(w){
      },
      finally = {}
    )
    
    #FE(BC) coefficients
    estimations["FE(BC) - Beta_x"] <- as.numeric(coef(op_FE)["x"])
    
    #FE(BC) - APe
    estimations["FE(BC) - APE(x)"] <- get_APEs(op_FE)$delta["x"]/true_APE_x
    
    #Checking rejection for latent FE(BC)
    coeff <- summary(op_FE)$cm
    rejection_matrix <- rejection_z_test(coeff,dynamic=FALSE)
    rejections_0.05["FE(BC) - Beta_x"] <- rejection_matrix["x","0.05"]
    rejections_0.1["FE(BC) - Beta_x"] <- rejection_matrix["x","0.1"]
    
    #Checking rejection for APE FE(BC)
    coeff <- summary(get_APEs(op_FE))
    rejection_matrix <- rejection_z_test(coeff,vector=true_APE_vector,dynamic=FALSE)
    rejections_0.05["FE(BC) - APE(x)"] <- rejection_matrix["x","0.05"]
    rejections_0.1["FE(BC) - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #LPM####
    #LPM with all the observations
    op_LPM <- plm(y ~ x,data=data,effect="individual",model="within")
    
    #Coefficients
    estimations["LPM - APE(x)"] <- coef(op_LPM)["x"]/true_APE_x
    
    #Checking rejections for LPM
    coeff <- summary(op_LPM)$coefficients
    df <-summary(op_LPM)$df[2]
    rejection_matrix <- rejection_t_test(coeff=coeff,vector=true_APE_vector,df=df,dynamic=FALSE)
    rejections_0.05["LPM - APE(x)"] <- rejection_matrix["x","0.05"]
    rejections_0.1["LPM - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #Reducing data to observations that vary over time
    data <- reduce_data(data)
    
    #RE####
    #RE coefficients
    tryCatch(
      expr = {
        op_RE <- pglm(y ~ x + x_bar,data = data,family = binomial(link='probit'),
                      method = 'nr', effect = "individual", model = "random")
      },
      error = function(e){
        return(NA) #We return NA and then we'll delete them
      },
      warning = function(w){
      },
      finally = {}
    )    
    
    #Getting latent variables RE
    beta_x_hat <- estimations["RE(R) - Beta_x"] <- as.numeric(coef(op_RE)["x"])
    beta_vector_hat <- c(0,beta_x_hat)
    
    #Getting estimated alphas and sigma
    alpha_0_hat <- as.numeric(coef(op_RE)["(Intercept)"])
    alpha_3_hat <- as.numeric(coef(op_RE)["x_bar"])
    sigma_hat <- as.numeric(coef(op_RE)["sigma"])  
    
    #Putting alphas in a vector
    alpha_hat <- c(alpha_0_hat,0,alpha_3_hat)
    
    #Getting APE for RE
    APE_hat <- APE_probit_RE(data,beta_vector_hat,alpha_hat,sigma_hat,dynamic=FALSE)
    estimations["RE(R) - APE(x)"] <- APE_hat[2]/true_APE_x
    
    #Checking rejection for latent RE
    coeff <- summary(op_RE)$estimate
    rejection_matrix <- rejection_t_test(coeff,dynamic=FALSE)
    rejections_0.05["RE(R) - Beta_x"] <- rejection_matrix["x","0.05"]
    rejections_0.1["RE(R) - Beta_x"] <- rejection_matrix["x","0.1"]
    
    #Checking rejection for APE RE (Bootstrap)
    #rejection_matrix <- rejection_APE_RE_Bootstrap(data,true_APE_vector,df=(nrow(data)/t-1)*(t-1),dynamic=FALSE)
    #rejections_0.05["RE(R) - APE(x)"] <- rejection_matrix["x","0.05"]
    #rejections_0.1["RE(R) - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #LPM with reduced dataset
    op_LPM <- plm(y ~ x,data=data,effect="individual",model="within")
    
    #Coefficients
    estimations["LPM(R) - APE(x)"] <- coef(op_LPM)["x"]/true_APE_x
    
    #Checking rejections for reduced LPM
    coeff <- summary(op_LPM)$coefficients
    df <-summary(op_LPM)$df[2]
    rejection_matrix <- rejection_t_test(coeff=coeff,vector=true_APE_vector,df=df,dynamic=FALSE)
    rejections_0.05["LPM(R) - APE(x)"] <- rejection_matrix["x","0.05"]
    rejections_0.1["LPM(R) - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #Returning esitmations and rejection data
    return(c(estimations,rejections_0.05,rejections_0.1))
  }
  
  #Generating list of dataframes
  df_list <- list()
  for (i in 1:m){
    df_list <- c(list(generate_static_probit()),df_list)
  }
  
  #Applying estimations to each dataframe
  estimations_rejections <-mclapply(df_list, FUN = montecarlo_static_RE, mc.cores=2)
  
  #Transforming the results
  t_matrix <- matrix(NA, ncol = length(names)*3, nrow = 0)
  for (i in 1:length(estimations_rejections)){
    t_matrix <- rbind(t(unlist(estimations_rejections[i])),t_matrix)
  }
  
  #Splitting the results
  estimations <- as.matrix(t_matrix[,1:10])
  estimations_p_0.05 <- as.matrix(t_matrix[,11:20])
  estimations_p_0.1 <- as.matrix(t_matrix[,21:30])
  
  #Deleting NA rows from estimations matrices
  estimations <- na.omit(estimations)
  
  #Updating the number of succesful replications
  k <- nrow(estimations)
  
  #Here we obtain results of interests from estimation matrix
  obtain_results <- function (estimations){
    results <- matrix(0, ncol=ncol(estimations), nrow=7, 
                      dimnames = list(c("Mean Bias","Median Bias", "MSE", "MAE", "SD", "p_0.05","p_0.1"), names))
    
    #Filling in SD,p0.05,p0.1 (this is the same for all the esitmators)
    for (i in 1:ncol(estimations)){
      results["SD",i] <- sd(estimations[,i])
      results["p_0.05",i] <- mean(estimations_p_0.05[,i])
      results["p_0.1",i] <- mean(estimations_p_0.1[,i])
    }
    
    #MSE and MAE for estimators of beta_x
    for (i in 1:3){
      results["Mean Bias",i] <- mean(estimations[,i]-beta_x)
      results["Median Bias",i] <- median(estimations[,i]-beta_x)
      results["MSE",i] <- mean((estimations[,i]-beta_x)^2)
      results["MAE",i] <- mean(abs(estimations[,i]-beta_x))
    }
    
    #MSE and MAE for estimators of APE (1 is the true value since it's the ratio)
    for (i in 4:ncol(estimations)){
      results["Mean Bias",i] <- mean(estimations[,i]-1)
      results["Median Bias",i] <- median(estimations[,i]-1)
      results["MSE",i] <- mean((estimations[,i]-1)^2)
      results["MAE",i] <- mean(abs(estimations[,i]-1))
    }
    
    return(t(results)) #Transpose is easier to read for me
  }
  results_static_probit <- obtain_results(estimations)
  
  #Printing results
  printing <- function(){
    closeAllConnections()
    sink("Results_Static_Probit.txt")
    cat("Number of indivuduals (n) =",n)
    cat(sep="\n")
    cat("Time (t) =",t)
    cat(sep="\n")
    cat("Number of initial iterations =",m)
    cat(sep="\n")
    cat("Number of actual iterations =",k)
    cat(sep="\n")
    cat("Success rare of simulation =",k/m)
    cat(sep="\n")
    cat("Bootstrap Replications =",R)
    cat(sep="\n")
    cat("Beta = ",beta_vector)
    cat(sep="\n")
    cat("Alpha = ", alpha)
    cat(sep="\n")
    cat("Sigma = ", sigma)
    cat(sep="\n")
    print(results_static_probit,sep="\n")
    sink()
  } 
  printing()
  return(results_static_probit)
}
#Execution Static Probit
results_static_probit <- run_static_probit()

#MONTE CARLO DYNAMIC LOGIT (NOT RELEVANT)####

beta_y <- beta_vector[1] <- 1

run_dynamic_logit <- function(){
  names <- c("RE - Beta_y","FE - Beta_y","FE(BC) - Beta_y",
             "RE - Beta_x","FE - Beta_x","FE(BC) - Beta_x",
             "RE - APE(y)","FE - APE(y)","FE(BC) - APE(y)","LPM-FS - APE(y)","LPM - APE(y)",
             "RE - APE(x)","FE - APE(x)","FE(BC) - APE(x)","LPM-FS - APE(x)","LPM - APE(x)")
  
  # Running Monte Carlo
  montecarlo_dynamic_logit <- function(data){
    estimations <- rep(NA,length(names)) #rows - replications; columns - different methods
    names(estimations) <- names
    rejections_0.05 <- rep(NA,length(names)) #rows - replications; columns - different methods
    names(rejections_0.05) <- names
    rejections_0.1 <- rep(NA,length(names)) #rows - replications; columns - different methods
    names(rejections_0.1) <- names
    
    #True APE
    APE <- true_APE_logit_RE(data)
    true_APE_y <- APE[1]
    true_APE_x <- APE[2]
    true_APE_vector <- c(true_APE_y,true_APE_x)
    
    #RE####
    #RE coefficients
    op_RE <- pglm(y ~ y_lag + x + x_bar + y_i0,data = data,family = binomial(link='logit'),
                  method = 'nr', effect = "individual", model = "random")
    
    #Getting latent variables RE
    beta_y_hat <- estimations["RE - Beta_y"] <- as.numeric(coef(op_RE)["y_lag"])
    beta_x_hat <- estimations["RE - Beta_x"] <- as.numeric(coef(op_RE)["x"])
    
    #Creating a vector of coefficients
    beta_vector_hat <- c(beta_y_hat,beta_x_hat)
    
    #Getting estimated alphas and sigma
    alpha_0_hat <- as.numeric(coef(op_RE)["(Intercept)"])
    alpha_1_hat <- as.numeric(coef(op_RE)["y_i0"])
    alpha_3_hat <- as.numeric(coef(op_RE)["x_bar"])
    sigma_hat <- as.numeric(coef(op_RE)["sigma"])  
    
    #Putting alphas in a vector
    alpha_hat <- c(alpha_0_hat,alpha_1_hat,alpha_3_hat)
    
    #Getting APE for RE
    APE_hat <- APE_logit_RE(data,beta_vector_hat,alpha_hat,sigma_hat)
    estimations["RE - APE(y)"] <- APE_hat[1]/true_APE_y
    estimations["RE - APE(x)"] <- APE_hat[2]/true_APE_x
    
    #Checking rejection for latent RE
    coeff <- summary(op_RE)$estimate
    rejection_matrix <- rejection_t_test(coeff)
    rejections_0.05["RE - Beta_y"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["RE - Beta_y"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["RE - Beta_x"] <- rejection_matrix["x","0.05"]
    rejections_0.1["RE - Beta_x"] <- rejection_matrix["x","0.1"]
    
    #Checking rejection for APE RE (Bootstrap)
    #rejection_matrix <- rejection_APE_RE_Bootstrap(data,beta_vector_hat,alpha_hat,sigma_hat,true_APE_vector)
    #rejections_0.05["RE - APE(y)"] <- rejection_matrix["y_lag","0.05"]
    #rejections_0.1["RE - APE(y)"] <- rejection_matrix["y_lag","0.1"]
    #rejections_0.05["RE - APE(x)"] <- rejection_matrix["x","0.05"]
    #rejections_0.1["RE - APE(x)"] <- rejection_matrix["x","0.1"]
    #FE####
    op_FE <- bife(y ~ y_lag + x | ID,
                  data = data, model = "logit")
    
    #FE coefficients
    estimations["FE - Beta_y"] <- as.numeric(coef(op_FE)["y_lag"])
    estimations["FE - Beta_x"] <- as.numeric(coef(op_FE)["x"])
    
    #FE - APE
    estimations["FE - APE(y)"] <- get_APEs(op_FE)$delta["y_lag"]/true_APE_y
    estimations["FE - APE(x)"] <- get_APEs(op_FE)$delta["x"]/true_APE_x
    
    #Checking rejection for latent FE
    coeff <- summary(op_FE)$cm
    rejection_matrix <- rejection_z_test(coeff)
    rejections_0.05["FE - Beta_y"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["FE - Beta_y"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["FE - Beta_x"] <- rejection_matrix["x","0.05"]
    rejections_0.1["FE - Beta_x"] <- rejection_matrix["x","0.1"]
    
    #Checking rejection for APE FE
    coeff <- summary(get_APEs(op_FE))
    rejection_matrix <- rejection_z_test(coeff,vector=true_APE_vector)
    rejections_0.05["FE - APE(y)"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["FE - APE(y)"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["FE - APE(x)"] <- rejection_matrix["x","0.05"]
    rejections_0.1["FE - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #FE(BC)####
    op_FE <- bias_corr(op_FE, L = 1L)
    
    #FE(BC) coefficients
    estimations["FE(BC) - Beta_y"] <- as.numeric(coef(op_FE)["y_lag"])
    estimations["FE(BC) - Beta_x"] <- as.numeric(coef(op_FE)["x"])
    
    #FE(BC) - APe
    estimations["FE(BC) - APE(y)"] <- get_APEs(op_FE)$delta["y_lag"]/true_APE_y
    estimations["FE(BC) - APE(x)"] <- get_APEs(op_FE)$delta["x"]/true_APE_x
    
    #Checking rejection for latent FE(BC)
    coeff <- summary(op_FE)$cm
    rejection_matrix <- rejection_z_test(coeff)
    rejections_0.05["FE(BC) - Beta_y"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["FE(BC) - Beta_y"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["FE(BC) - Beta_x"] <- rejection_matrix["x","0.05"]
    rejections_0.1["FE(BC) - Beta_x"] <- rejection_matrix["x","0.1"]
    
    #Checking rejection for APE FE(BC)
    coeff <- summary(get_APEs(op_FE))
    rejection_matrix <- rejection_z_test(coeff,vector=true_APE_vector)
    rejections_0.05["FE(BC) - APE(y)"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["FE(BC) - APE(y)"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["FE(BC) - APE(x)"] <- rejection_matrix["x","0.05"]
    rejections_0.1["FE(BC) - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #LPM####
    #LPM with all the observations
    op_LPM <- plm(y ~ y_lag + x,data=data,effect="individual",model="within")
    
    #Coefficients
    estimations["LPM-FS - APE(y)"] <- coef(op_LPM)["y_lag"]/true_APE_y
    estimations["LPM-FS - APE(x)"] <- coef(op_LPM)["x"]/true_APE_x
    
    #Checking rejections for LPM
    coeff <- summary(op_LPM)$coefficients
    df <-summary(op_LPM)$df[2]
    rejection_matrix <- rejection_t_test(coeff=coeff,vector=true_APE_vector,df=df)
    rejections_0.05["LPM-FS - APE(y)"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["LPM-FS - APE(y)"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["LPM-FS - APE(x)"] <- rejection_matrix["x","0.05"]
    rejections_0.1["LPM-FS - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #LPM only with observations which vary over time
    data_reduced <- reduce_data(data)
    op_LPM <- plm(y ~ y_lag + x,data=data_reduced,effect="individual",model="within")
    
    #Coefficients
    estimations["LPM - APE(y)"] <- coef(op_LPM)["y_lag"]/true_APE_y
    estimations["LPM - APE(x)"] <- coef(op_LPM)["x"]/true_APE_x
    
    #Checking rejections for reduced LPM
    coeff <- summary(op_LPM)$coefficients
    df <-summary(op_LPM)$df[2]
    rejection_matrix <- rejection_t_test(coeff=coeff,vector=true_APE_vector,df=df)
    rejections_0.05["LPM - APE(y)"] <- rejection_matrix["y_lag","0.05"]
    rejections_0.1["LPM - APE(y)"] <- rejection_matrix["y_lag","0.1"]
    rejections_0.05["LPM - APE(x)"] <- rejection_matrix["x","0.05"]
    rejections_0.1["LPM - APE(x)"] <- rejection_matrix["x","0.1"]
    
    #Returning esitmations and rejection data
    return(c(estimations,rejections_0.05,rejections_0.1))
  }
  
  #Generating list of dataframes
  df_list <- list()
  for (i in 1:m){
    df_list <- c(list(generate_dynamic_probit_random()),df_list)
  }
  
  #Applying estimations to each dataframe
  estimations_rejections <-mclapply(df_list, FUN = montecarlo_dynamic_logit, mc.cores=2)
  
  #Transforming the results
  t_matrix <- matrix(NA, ncol = length(names)*3, nrow = 0)
  for (i in 1:length(estimations_rejections)){
    if (!is.na(estimations_rejections[i])){
      t_matrix <- rbind(t(unlist(estimations_rejections[i])),t_matrix)
    }
  }
  
  #Splitting the results
  estimations <- t_matrix[,1:16]
  estimations_p_0.05 <- t_matrix[,17:32]
  estimations_p_0.1 <- t_matrix[,33:48]
  
  #Updating the number of succesful replications
  k <- nrow(estimations)
  
  #Here we obtain results of interests from estimation matrix
  
  obtain_results <- function (estimations){
    results <- matrix(0, ncol=ncol(estimations), nrow=6, 
                      dimnames = list(c("Bias", "MSE", "MAE", "SD", "p_0.05","p_0.1"), names))
    
    #Filling in mean,SD,p0.05,p0.1 (this is the same for all the esitmators)
    for (i in 1:ncol(estimations)){
      results["SD",i] <- sd(estimations[,i])
      results["p_0.05",i] <- mean(estimations_p_0.05[,i])
      results["p_0.1",i] <- mean(estimations_p_0.1[,i])
    }
    
    #MSE and MAE for estimators of beta_y
    for (i in 1:3){
      results["Bias",i] <- mean(estimations[,i]-beta_y)
      results[2,i] <- mean((estimations[,i]-beta_y)^2) 
      results[3,i] <- mean(abs(estimations[,i]-beta_y))
    }
    
    #MSE and MAE for estimators of beta_x
    for (i in 4:6){
      results["Bias",i] <- mean(estimations[,i]-beta_x)
      results[2,i] <- mean((estimations[,i]-beta_x)^2)
      results[3,i] <- mean(abs(estimations[,i]-beta_x))
    }
    
    #MSE and MAE for estimators of APE (1 is the true value since it's the ratio)
    for (i in 7:ncol(estimations)){
      results["Bias",i] <- mean(estimations[,i]-1)
      results[2,i] <- mean((estimations[,i]-1)^2)
      results[3,i] <- mean(abs(estimations[,i]-1))
    }
    return(t(results)) #Transpose is easier to read for me
  }
  
  results <- obtain_results(estimations)
  
  #Printing results
  printing <- function(){
    closeAllConnections()
    sink("Results_Dynamic_Logit.txt")
    cat("Number of indivuduals (n) =",n)
    cat(sep="\n")
    cat("Time (t) =",t)
    cat(sep="\n")
    cat("Number of initial iterations =",m)
    cat(sep="\n")
    cat("Number of actual iterations =",k)
    cat(sep="\n")
    cat("Success rare of simulation =",k/m)
    cat(sep="\n")
    cat("Beta = ",beta_vector)
    cat(sep="\n")
    cat("Alpha = ", alpha)
    cat(sep="\n")
    cat("Sigma = ", sigma)
    cat(sep="\n")
    print(results,sep="\n")
    sink()
  }
  printing()
}

#PROVING THAT SIGMA IS WRONG####
t <- 8 #time (number of years)
n <- 250 #number of individuals
beta_x <- 1 #parameter for x
beta_y <- 0 #parameter for lagged y
alpha_0 <- 1
alpha_1 <- 0
alpha_3 <- 0
alpha <- c(alpha_0,alpha_1,alpha_3)
sigma <- 1

pglm_est <- function(m){
  est_coefficients <- matrix(NA,nrow=m,ncol=3)
  for (i in 1:m){
    data <- generate_static_probit()
    op_RE <- pglm(y ~ x,data = data,family = binomial(link='probit'),
                  method = 'nr', effect = "individual", model = "random")
    est_coefficients[i,1] <- as.numeric(coef(op_RE)["(Intercept)"])
    est_coefficients[i,2] <- as.numeric(coef(op_RE)["x"])
    est_coefficients[i,3] <- as.numeric(coef(op_RE)["sigma"])
  }
  return(est_coefficients)
}
#est_coefficients <- pglm_est(100) #Uncomment to check why PGLM doesn't esitmate sigma correctly
#coeff_means <- c(mean(est_coefficients[,1]),mean(est_coefficients[,2]),mean(est_coefficients[,3]))