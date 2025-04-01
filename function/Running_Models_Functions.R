run_lm <- function(formula, df, num){
  
  # Purpose: Fitting a lm model of the given formula, 
  #          and returning required statistics
  
  # Arguments: formula - (string) the formula of the model
  #            df - (dataframe) the dataframe containing the data
  #            num - (int) the number of iteration (only used for return value)
  
  # Returns: a dataframe contains columns necessary for simsum
  
  model <- eval(bquote( lm(.(as.formula(formula)), data=df) ))
  summ <- summary(model)
  point_estimate = summ$coefficients[2,1]
  se = summ$coefficients[2,2]
  p_value <- summ$coefficients[2,4]
  
  data.frame(
    num = i,
    method = "LR",
    estimand = "mean difference",
    b = point_estimate,
    se = se,
    p_value = p_value
  )
  
}

run_log_lm <- function(formula, df, num){
  
  # Purpose: Fitting a lm model on the log-transformed outcome data 
  #          using the given formula, and returning required statistics
  
  # Arguments: formula - (string) the formula of the model
  #            df - (dataframe) the dataframe containing the data
  #            num - (int) the number of iteration (only used for return value)
  
  # Returns: a dataframe contains columns necessary for simsum
  
  df$outcome <- log(df$outcome)
  model <- eval(bquote( lm(.(as.formula(formula)), data=df)))
  summ <- summary(model)
  point_estimate = summ$coefficients[2,1]
  se = summ$coefficients[2,2]
  p_value <- summ$coefficients[2,4]
  
  data.frame(
    num = i,
    method = "Log_LR",
    estimand = "log median ratio",
    b = point_estimate,
    se = se,
    p_value = p_value
  )
}

run_gaussian_log_link <- function(formula, df, num){
  
  # Purpose: Fitting a Gaussian glm model with log link using the given 
  #          formula, and returning the required statistics
  
  # Arguments: formula - (string) the formula of the model
  #            df - (dataframe) the dataframe containing the data
  #            num - (int) the number of iteration (only used for return value)
  
  # Returns: a dataframe contains columns necessary for simsum
  
  model <- eval(bquote( glm(.(as.formula(formula)), data=df, family=gaussian(link="log"))))
  summ <- summary(model)
  point_estimate = summ$coefficients[2,1]
  se = summ$coefficients[2,2]
  p_value <- summ$coefficients[2,4]
  
  data.frame(
    num = i,
    method = "GLM_Gaussian",
    estimand = "log mean ratio",
    b = point_estimate,
    se = se,
    p_value = p_value
  )
  
}

run_gamma_log_link <- function(formula, df, num){
  
  # Purpose: Fitting a Gamma glm model with log link using the given 
  #          formula, and returning the required statistics
  
  # Arguments: formula - (string) the formula of the model
  #            df - (dataframe) the dataframe containing the data
  #            num - (int) the number of iteration (only used for return value)
  
  # Returns: a dataframe contains columns necessary for simsum
  
  model <- eval(bquote( glm(.(as.formula(formula)), data=df, family=Gamma(link="log"))))
  summ <- summary(model)
  point_estimate = summ$coefficients[2,1]
  se = summ$coefficients[2,2]
  p_value <- summ$coefficients[2,4]
  
  data.frame(
    num = i,
    method = "GLM_Gamma",
    estimand = "log mean ratio",
    b = point_estimate,
    se = se,
    p_value = p_value
  )
  
}

run_olr <- function(formula, df, num){
  
  # Purpose: Fitting a ordinal logistic regression model using the given 
  #          formula, and returning the required statistics
  
  # Arguments: formula - (string) the formula of the model
  #            df - (dataframe) the dataframe containing the data
  #            num - (int) the number of iteration (only used for return value)
  
  # Returns: a dataframe contains columns necessary for simsum
  
  dd <- datadist(df)
  options(datadist = dd)
  model <- eval(bquote( orm(.(as.formula(formula)), data=df, se.fit=TRUE)))
  summ <- summary(model)
  
  point_estimate = summ['treatment','Effect']
  se = summ['treatment', "S.E."]
  z_score <- point_estimate / se
  p_value <- 2 * (1 - pnorm(abs(z_score)))
  
  data.frame(
    num = i,
    method = "OLR",
    estimand = "log odds?",
    b = point_estimate,
    se = se,
    p_value = p_value
  )
}


