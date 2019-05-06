
# Modified Gale-Shapley algorithm #
# This script is used to obtain Counterfactuals based on dating data in combination with census data from IPUMS #
# Sections:
#   0. Parameters: Sets the parameters for the counterfactual
#   1. Input data: This section is used to input the dating data and the census data used throughout the code
#   2. Function design: This section builds the functions needed to obtain the counterfactuals:
#                     - Match_fun: Matching function - matches individuals based on their preferences and education attainment
#                     - ToOpt: Used when we want to estimated the number of voluntary singles (singles by choice)
#                     - ObjectiveF: Objective function to be minimized in order to estimate the shares of singles.
#   3. Matching Algorithm: This section wraps all the code and runs the counterfactuals.
#   4. Results


#########################################################################################################################################
#                                                   0.PARAMETERS                                                                        #
#########################################################################################################################################

# rm(list = setdiff(ls(), lsf.str()))
#set.seed(0)

modified_gale_shapley = function(year_edu, dim, counterfactual_year, estimate_vol_singles = "yes",
                                 penalized, supress_results = F, supress_graphics = F){
  
  
  # If counterfactual_year = 9999 then no counterfactual is estimated.
  # This counterfactual uses the estimated shares to build a new counterfactual with fixed shares and
  #   preferences from year_edu and census data from counterfactual_year
  
  
  dim = dim
  is_with_singles = 1 #(alternatives are 1 and 0)
  yearOfEduDistr = year_edu #(acepted values: 1980, 2000, 2011)
  yearOfEduDistr_Block3 = counterfactual_year #(alternatives are 1980 or 2000)
  estimate_volunt_singles = estimate_vol_singles #(alternatives are "yes" and "no")
  optim_penalty = penalized
  set.seed(12345)
  
  print(paste0("Year of education distribution is ", yearOfEduDistr))
  #########################################################################################################################################
  #                                                   1.INPUT DATA                                                                        #
  #########################################################################################################################################
  
  Distr_F_Mat = numeric() # Female vector
  NDistr_M_Mat = numeric() # Male vector
  
  #################################################
  #               CENSUS DATA                     #
  #################################################
  if (yearOfEduDistr == 2011) {
    NDistr_Couples = matrix(
      c(46080, 34640, 2460,
        31140, 237380, 78460,
        560, 28540, 95640),
      nrow = dim,
      ncol = dim,
      byrow = T)
  } else {
    if (yearOfEduDistr == 2000){
      NDistr_Couples = matrix(
        c(78480, 54100, 2620,
          78580, 312580, 44480,
          1960, 31140, 52420),
        nrow = dim,
        ncol = dim,
        byrow = T)
    } else {
      NDistr_Couples = matrix(
        c(231120, 67340, 360,
          194800, 332940, 8580,
          1840, 30100, 17080),
        nrow = dim,
        ncol = dim,
        byrow = T)
    }
  }
  
  
  # Census data metrics:
  # 1. Share of homogamous couples (L, M and H educated couples)
  census_shc_l = NDistr_Couples[1,1]/sum(NDistr_Couples)
  census_shc_m = NDistr_Couples[2,2]/sum(NDistr_Couples)
  census_shc_h = NDistr_Couples[3,3]/sum(NDistr_Couples)
  census_shc_total = sum(diag(NDistr_Couples))/sum(NDistr_Couples)
  
  print(paste0("Share of homogamous couple in census data ", round(census_shc_total*100,2),"%"))
  
  # Married Females
  NDistr_F_inrelation = matrix(colSums(NDistr_Couples), nrow = 1, ncol = 3, byrow = T)
  # Married Males
  NDistr_M_inrelation = matrix(colSums(t(NDistr_Couples)), nrow = 3, ncol = 1, byrow = T)
  
  
  # SINGLES
  # Single Females
  if (yearOfEduDistr == 2011) {
    NDistr_F_singles = matrix(
      c(136000, 315480, 61420),
      nrow = 1,
      ncol = 3,
      byrow = T)
  } else {
    if (yearOfEduDistr == 2000){
      NDistr_F_singles = matrix(
        c(97800, 221480, 74780),
        nrow = 1,
        ncol = 3,
        byrow = T)
    } else {
      NDistr_F_singles = matrix(
        c(78100, 95680, 12440),
        nrow = 1,
        ncol = 3,
        byrow = T)
    }
  }
  
  
  # Single Males
  if (yearOfEduDistr == 2011) {
    NDistr_M_singles = matrix(
      c(128060,
        471540,
        145860),
      nrow = 3,
      ncol = 1,
      byrow = T
    )
    
  } else {
    if (yearOfEduDistr == 2000){
      NDistr_M_singles = matrix(
        c(136000,
          315480,
          61420),
        nrow = 3,
        ncol = 1,
        byrow = T)
    } else {
      NDistr_M_singles = matrix(
        c(107420,
          128140,
          11100),
        nrow = 3,
        ncol = 1,
        byrow = T)
    }
  }
  
  
  # Marginal distributions of Females and Males
  # Females
  if (is_with_singles == 1) {
    NDistr_F = NDistr_F_inrelation + NDistr_F_singles
  } else {
    NDistr_F = NDistr_F_inrelation
  }
  
  # Males
  if (is_with_singles == 1) {
    NDistr_M = NDistr_M_inrelation + NDistr_M_singles
  } else {
    NDistr_M = NDistr_M_inrelation
  }
  
  #################################################
  #               DATING DATA                     #
  #################################################
  nPref_F_Mat = matrix(c(19, 22, 2,
                         17, 156, 83,
                         2, 7, 47), nrow = dim, ncol = dim, byrow = T)
  
  
  nPref_M_Mat = matrix(c(46, 23, 0,
                         57, 149, 3,
                         9, 61, 32), nrow = dim, ncol = dim, byrow = T)
  
  # Adding dummy row (column) for voluntary females (males)
  if (is_with_singles == 1) {
    nPref_F_Mat = rbind(nPref_F_Mat, 0)
  }
  
  if (is_with_singles == 1) {
    nPref_M_Mat = cbind(nPref_M_Mat, 0)
  }
  
  
  # Calculating frequencies for each cell in the matrices of preferences
  rPref_F_Mat = t(t(nPref_F_Mat) / rowSums(t(nPref_F_Mat)))
  rPref_M_Mat = nPref_M_Mat / rowSums(nPref_M_Mat)
  
  
  # Updating census marginals with dating data information
  NPref_F_Mat = t(t(rPref_F_Mat) * as.vector(t(NDistr_F)))
  NPref_M_Mat = rPref_M_Mat * as.vector(NDistr_M)
  
  
  
  
  #########################################################################################################################################
  #                                                   2.FUNCTION DESIGN                                                                   #
  #########################################################################################################################################
  
  
  #################################################
  #            2.1. Match function                #
  #################################################
  
  Match_fun = function(Female_preferences, Male_preferences, include_singles = 1){
    
    # Defining parameters in case the analysis includes singles
    if (include_singles == 1) {
      dim_withSingle = ncol(Male_preferences)
      dim = dim_withSingle - 1
      
      NPref_F_Mat_NoSingle = Female_preferences[1:dim, ]   # dim = (3x3)
      NPref_M_Mat_NoSingle = Male_preferences[, 1:dim]   # dim = (3x3)
      
      NPref_F_Single = Female_preferences[dim_withSingle, ]  # Is vector. Length = 3
      NPref_M_Single = Male_preferences[, dim_withSingle]  # Is vector. Length = 3
      
    } else {
      dim = nrow(Male_preferences)
      
      NPref_F_Mat_NoSingle = Female_preferences   # dim = (3x3)
      NPref_M_Mat_NoSingle = Male_preferences   # dim = (3x3)
      
      NPref_F_Single = numeric()
      NPref_M_Single = numeric()
    }
    
    
    Vec_M = numeric()
    Vec_F = numeric()
    
    #Index i identifies males and j identifies females
    for (i in seq(dim, 1, by = -1)) {
      #Loops from most educated to less educated males
      for (j in seq(dim, 1, by = -1)) {
        #Loops from most educated to less educated women
        
        Vec_M = rbind(Vec_M, NPref_M_Mat_NoSingle[i, j])
        Vec_F = rbind(Vec_F, NPref_F_Mat_NoSingle[j, i])
      }
    }
    
    
    #First proposals are made to the best available female
    index_F = 0
    Vec_F_old = Vec_F
    Vec_M_old = Vec_M
    
    NPref_M_Mat_old = NPref_M_Mat_NoSingle
    NPref_F_Mat_old = NPref_F_Mat_NoSingle
    Cumul_M = 0
    
    
    #Output matrices
    
    Matched_Mat = matrix(0, nrow = dim, ncol = dim)
    Resid_F_Mat = matrix(0, nrow = dim, ncol = dim)
    Resid_M_Mat = matrix(0, nrow = dim, ncol = dim)
    
    
    # Matching Loop. Index i identifies males and index j identifies females
    # Matching will be made on the male side.
    # Ordered males choose from the females in educational decreasing order
    
    for (i in seq(dim, 1, by = -1)) {
      for (j in seq(dim, 1, by = -1)) {
        index_F = index_F + 1   #Index of Vec_F
        
        Cumul_M = sum(NPref_M_Mat_old[j:dim, 1:i])  #cumulative sum of males by education level
        # i=3, j=3 sums H educ males across L,M & H females
        # i=3, j=2 sums H & M educ males across L,M & H females
        # i=3, j=1 sums H,M & L educ males across L,M & H females
        
        # i=2, j=3 sums H educ males across L & M females
        # ...
        
        matched = min(Vec_F_old[index_F], Cumul_M)
        Vec_F_new = Vec_F_old
        Vec_F_new[index_F] = Vec_F_old[index_F] - matched
        Resid_F_Mat[j, i] = NPref_F_Mat_NoSingle[j, i] - matched
        
        NPref_M_Mat_new = NPref_M_Mat_old
        
        for (v in seq(dim, j, by = -1)) {
          #Rows
          for (m in seq(i, 1, by = -1)) {
            #Columns
            
            matchedwith = min(NPref_M_Mat_old[v, m], matched)
            NPref_M_Mat_new[v, m] = NPref_M_Mat_old[v, m] - matchedwith
            matched = matched - matchedwith
            Matched_Mat[v, i] = Matched_Mat[v, i] + matchedwith
            
          }
        }
        
        Vec_F_old = Vec_F_new
        NPref_M_Mat_old = NPref_M_Mat_new
        
      }
    }
    
    
    # Output matrices
    Resid_M_Mat = NPref_M_Mat_new
    
    Resid_M_Mat = matrix(cbind(Resid_M_Mat, NPref_M_Single), nrow = dim)
    Resid_F_Mat = matrix(rbind(Resid_F_Mat, NPref_F_Single), ncol = dim)
    
    
    
    return(list(Matched_Mat,
                Resid_M_Mat,
                Resid_F_Mat))
    
    
  }
  
  
  #################################################
  #            2.2. ToOpt function                #
  #################################################
  
  ToOpt = function(Vec_param_init, penalty = 1){
    
    NDistr_F_singles = NDistr_F_singles
    NDistr_M_singles = NDistr_M_singles
    
    # Running objective function minimization to find vector of voluntary singles
    
    value_ = matrix(NA, nrow = 100, ncol = 6)
    estimates = rep(0, 6)
    bounds = c(NDistr_F, NDistr_M)
    
    for (j in 1:6){
      param = estimates
      for (i in 0:99){
        param[j] = i/100 * bounds[j]
        value_[(i+1),j] = ObjectiveF(param, penalty = optim_penalty)
      }
      idx = which.min(value_[,j])
      estimates[j] = (idx-1)/100*bounds[j]
    }
    
    print("Initital values")
    print(estimates)
    
    vec_par = estimates

    for (i in 1:20){
    
      vec_par = Sim.Annealing(vec_par, penalty, cycles = 500, trials = 20)
      #print(paste("Annealing cycle:", i," Estimated vector"))
      #print(vec_par)
       
    }
    
    value_ = matrix(NA, nrow = 100, ncol = 6)
    estimates = vec_par
    bounds = c(NDistr_F, NDistr_M)
    
    for (j in 1:6){
      param = estimates
      for (i in 0:99){
        param[j] = i/100 * bounds[j]
        value_[(i+1),j] = ObjectiveF(param, penalty = optim_penalty)
      }
      idx = which.min(value_[,j])
      estimates[j] = (idx-1)/100*bounds[j]
    }
    
    
    if (supress_graphics == F){
      par(mfrow=c(2,3))
      for(j in 1:6){
        ts.plot(value_[,j])
      }
      par(mfrow=c(1,1))
    }
    

    # Use the information of the boundaries on the number of singles
    # to create an algorithm which separately finds the minimum for the objective function.
    
    #Optimize = Sim.Annealing(initial_values, penalty = penalty)
    
    # Estimates
    Vec_param_est = estimates
    
    # Objective function value
    Value = ObjectiveF(Vec_param_est, penalty = penalty)
    
    # Normalization
    Vec_param_est = pmax(Vec_param_est,0)
    
    # Female and male voluntary singles
    VoluntSingle_F = t(Vec_param_est[1:3])
    VoluntSingle_M = t(t(Vec_param_est[4:6]))
    
    
    # Frequencies
    rVoluntSingle_F = VoluntSingle_F / NDistr_F
    rVoluntSingle_M = VoluntSingle_M / NDistr_M
    
    # Updating frequencies with voluntary singles estimates
    rPref_F_Mat=t(t(nPref_F_Mat) / rowSums(t(nPref_F_Mat))) * matrix((1 - rVoluntSingle_F),nrow = 4, ncol = 3, byrow = T)
    rPref_F_Mat = rbind(rPref_F_Mat[1:3, 1:3], rVoluntSingle_F)
    
    
    rPref_M_Mat = nPref_M_Mat / rowSums(nPref_M_Mat) * as.vector(1 - rVoluntSingle_M)
    rPref_M_Mat = cbind(rPref_M_Mat[1:3, 1:3], rVoluntSingle_M)
    
    
    # Updating contingency table
    NPref_F_Mat_new = t(t(rPref_F_Mat) * as.vector(t(NDistr_F)))
    
    NPref_M_Mat_new = rPref_M_Mat * as.vector(NDistr_M)
    
    
    # Running matching algorithm
    include_singles = 1
    
    MatchFun=Match_fun(NPref_F_Mat_new, NPref_M_Mat_new, include_singles = 1)
    Matched_Mat = MatchFun[[1]]
    Resid_M_Mat = MatchFun[[2]]
    Resid_F_Mat = MatchFun[[3]]
    
    # Output
    return(list(Value,
                Matched_Mat,
                Resid_M_Mat,
                Resid_F_Mat,
                Vec_param_est))
    
  }
  
  #################################################
  #            2.3. Objective function            #
  #################################################
  
  ObjectiveF = function(x, penalty = 1) {
    
    NDistr_F_singles = NDistr_F_singles
    NDistr_M_singles = NDistr_M_singles
    
    # Setting bounds on parameters
    bound=c(NDistr_F,NDistr_M)
    x = pmax(pmin(x, bound), 0)
    
    # Female and Male vector of parameters
    VoluntSingle_F = t(x[1:3])
    VoluntSingle_M = t(t(x[4:6]))
    
    # Setting bounds of voluntary singles
    VoluntSingle_F=pmax(pmin(VoluntSingle_F,NDistr_F),0)
    VoluntSingle_M=pmax(pmin(VoluntSingle_M,NDistr_M),0)
    
    # Frequencies
    rVoluntSingle_F = VoluntSingle_F / NDistr_F
    rVoluntSingle_M = VoluntSingle_M / NDistr_M
    
    
    # Updating preferences
    rPref_F_Mat=t(t(nPref_F_Mat) / rowSums(t(nPref_F_Mat))) * matrix((1 - rVoluntSingle_F),nrow = 4, ncol = 3, byrow = T)
    rPref_F_Mat = rbind(rPref_F_Mat[1:3, 1:3], rVoluntSingle_F)
    
    rPref_M_Mat = nPref_M_Mat / rowSums(nPref_M_Mat) * as.vector(1 - rVoluntSingle_M)
    rPref_M_Mat = cbind(rPref_M_Mat[1:3, 1:3], rVoluntSingle_M)
    
    NPref_F_Mat_new = t(t(rPref_F_Mat) * as.vector(t(NDistr_F)))
    NPref_M_Mat_new = rPref_M_Mat * as.vector(NDistr_M)
    
    # Run Matching algorithm
    include_singles = 1
    MatchFun=Match_fun(NPref_F_Mat_new, NPref_M_Mat_new, include_singles = 1)
    
    Matched_Mat = MatchFun[[1]]
    Resid_M_Mat = MatchFun[[2]]
    Resid_F_Mat = MatchFun[[3]]
    
    sumResid_F = colSums(pmax(Resid_F_Mat,0))
    sumResid_M = rowSums(pmax(Resid_M_Mat,0))
    
    sumResid_F=pmax(pmin(sumResid_F,NDistr_F),0)
    sumResid_M=pmax(pmin(sumResid_M,NDistr_M),0)
    
    est_shc_total <<- sum(diag(Matched_Mat))/sum(Matched_Mat)
    
    h = matrix(c(t((NDistr_F_singles - sumResid_F) / NDistr_F),
                 (NDistr_M_singles - sumResid_M) / NDistr_M),
               nrow = 6,
               ncol = 1,
               byrow = T)
    
    # Objective function value - Sums of square deviations
    Q = t(h) %*% h + penalty*(census_shc_total - est_shc_total)^2
    
    return(Q)
    
  }
  
  #################################################
  #            2.4. Solver Function            #
  #################################################
  
  Sim.Annealing=function(par, penalty, cycles = 100, trials = 100){
    # Number of cycles
    n = cycles
    
    # Number of trials per cycle
    m = trials
    
    # Number of accepted solutions
    na = 0.0
    
    # Probability of accepting worse solution at the start
    p1 = 0.8
    
    # Probability of accepting worse solution at the end
    p50 = 0.01
    
    # Initial temperature
    t1 = -1.0 / log(p1) 
    
    # Final temperature
    t50 = -1.0 / log(p50)
    
    # Fractional reduction every cycle
    frac = (t50 / t1) ^ (1.0 / (n - 1.0))
    
    # Initialize x
    
    x_start = par
    x = matrix(0, nrow = n + 1, ncol = 6)
    x[1, ] = x_start
    xi = c(0.0, 0.0,0.0,0.0,0.0,0.0)
    xi = x_start
    na = na + 1.0
    
    # Current best result so far
    xc = c(0.0, 0.0,0.0,0.0,0.0,0.0)
    xc = x[1,]
    fc = ObjectiveF(xi, penalty = penalty)
    fs = rep(0, n + 1)
    fs[1] = fc
    
    # Current temperature
    t = t1
    
    # DeltaE Average
    DeltaE_avg = 0.0
    
    for (i in 1:n) {
      # if (i%%50 == 0){
      #   msg_1 = paste("Cycle", i, "with temperature", t)
      #   print(msg_1)
      # }

      for (j in 1:m) {
        # Generate candidate points
        xi = xc + (runif(6, 0, 1) - 0.5) * 1000/j
        
        
        # Constrain search bounds
        xi = pmax(pmin(xi, c(NDistr_F, NDistr_M)), 0)
        
        func_val = ObjectiveF(xi, penalty = penalty)
        DeltaE = abs(func_val - fc)
        
        if (func_val > fc) {
          if (i == 0 & j == 0) {
            DeltaE_avg = DeltaE
          }
          
          p = exp(-DeltaE / (DeltaE_avg * t))
          
          if (runif(1, 0, 1) < p) {
            accept = TRUE
          } else {
            accept = FALSE
          }
          
        } else {
          accept = TRUE
        }
        
        if (accept == TRUE) {
          xc = xi
          
          fc = ObjectiveF(xc, penalty = penalty)
          
          na = na + 1.0
          
          DeltaE_avg = (DeltaE_avg * (na - 1.0) + DeltaE) / na
        }
        
        
      }
      
      x[i + 1, ] = xc
      
      fs[i+1]=fc
      
      if (runif(1,0,1)>0.99){
        t = t1
      } else {
        t = frac * t
      }

      
    }
    
    xc=x[which.min(fs),]
    
    return(xc)
  }
  
  
  #########################################################################################################################################
  #                                                   3. MATCHING ALGORITHM                                                               #
  #########################################################################################################################################
  
  if (estimate_volunt_singles == "no"){
    
    # Run male-optimum matching function: No voluntary singles estimation
    MatchFun = Match_fun(NPref_F_Mat, NPref_M_Mat, is_with_singles)
    Matched_Mat = MatchFun[[1]]
    Resid_M_Mat = MatchFun[[2]]
    Resid_F_Mat = MatchFun[[3]]
    
    contingency_1 = matrix(0, nrow = 5, ncol = 5)
    
    contingency_1[1:3, 1:3] = Matched_Mat
    contingency_1[4, 1:3] = Resid_F_Mat[4,]
    contingency_1[1:3, 4] = Resid_M_Mat[, 4]
    contingency_1[5, 1:3] = colSums(Resid_F_Mat[1:3,])
    contingency_1[1:3, 5] = rowSums(Resid_M_Mat[, 1:3])
    
  } else {
    
    # Run male-optimum matching function: Voluntary singles estimation
    print("Initializing optimization procedure")
    print(paste0("Running constrained penalized optimization. Penalty term equal to ", optim_penalty))
    
    
    Vec_param_init = 0.3*c(NDistr_F, NDistr_M)
    # Running ToOpt fun to obtain vector of voluntary singles
    ToOptFun = ToOpt(Vec_param_init, penalty = optim_penalty)
    
    # Results from optimization
    Value = ToOptFun[[1]]
    Matched_Mat = ToOptFun[[2]]
    Resid_M_Mat = ToOptFun[[3]]
    Resid_F_Mat = ToOptFun[[4]]
    Vec_param_est = ToOptFun[[5]]
    
    
    # Contingency table
    contingency = matrix(0, nrow = 5, ncol = 5)
    
    contingency[1:3, 1:3] = Matched_Mat
    contingency[4, 1:3] = Resid_F_Mat[4,]
    contingency[1:3, 4] = Resid_M_Mat[, 4]
    contingency[5, 1:3] = colSums(Resid_F_Mat[1:3,])
    contingency[1:3, 5] = rowSums(Resid_M_Mat[, 1:3])
    
    
    # Allocate share of voluntrary and involuntary singles to their own vectors
    Shares_Mat_Singles = matrix(c(t(t(Resid_F_Mat[4,])/rowSums(t(NDistr_F))),
                                  Resid_M_Mat[, 4] / rowSums(NDistr_M)),nrow = 6, ncol = 1, byrow = T)
    
    
    
    Shares_Mat_Coupled = matrix(c(t(t(t(t(nPref_F_Mat)/rowSums(t(nPref_F_Mat))))*(1 - Shares_Mat_Singles[1:3,]))[1:3,],
                                  t((nPref_M_Mat / rowSums(nPref_M_Mat) * (1 - Shares_Mat_Singles[4:6,]))[, 1:3])), nrow = 18, ncol = 1, byrow = T)
    
    
    
    rownames(Shares_Mat_Coupled) = c("f_LL",
                                     "f_LM",
                                     "f_LH",
                                     "f_ML",
                                     "f_MM",
                                     "f_MH",
                                     "f_HL",
                                     "f_HM",
                                     "f_HH",
                                     "m_LL",
                                     "m_LM",
                                     "m_LH",
                                     "m_ML",
                                     "m_MM",
                                     "m_MH",
                                     "m_HL",
                                     "m_HM",
                                     "m_HH")
    
    rownames(Shares_Mat_Singles) = c("f_LS", "f_MS", "f_HS",
                                     "m_LS", "m_MS", "m_HS")
    
    
    
    SHC = list()
    SHC$L = Matched_Mat[1, 1] / sum(Matched_Mat)
    SHC$M = Matched_Mat[2, 2] / sum(Matched_Mat)
    SHC$H = Matched_Mat[3, 3] / sum(Matched_Mat)
    SHC$Agg = sum(diag(Matched_Mat)) / sum(Matched_Mat)
    
    
    # Allocate singles to each own vector
    # Singles by change
    InvoluntarySingles = matrix(c(colSums(Resid_F_Mat[1:3, ]), rowSums(Resid_M_Mat[, 1:3])), nrow = 6, ncol = 1, byrow = T)
    
    # Singles by choice
    Voluntary_Singles = matrix(Vec_param_est, nrow = 6, ncol = 1, byrow = T)
    
    
    rownames(Matched_Mat) = c("L", "M", "H")
    colnames(Matched_Mat) = c("L", "M", "H")
    rownames(InvoluntarySingles) = c("F_L", "F_M", "F_H", "M_L", "M_M", "M_H")
    rownames(Voluntary_Singles) = c("F_L", "F_M", "F_H", "M_L", "M_M", "M_H")
    constraint_bind = census_shc_total - SHC$Agg
    
    
  }
  #########################################################################################################################################
  #                                                   4. RESULTS                                                                          #
  #########################################################################################################################################
  if (supress_results == F){
    
    
    print(paste("Objective Function Value:", Value, sep = " "))
  
    print("Estimated Number of Involuntary Singles:")
    print(InvoluntarySingles)
    
    print("Estimated Number of Voluntary Singles:")
    print(Voluntary_Singles)
    
    print("Match matrix")
    print(Matched_Mat)
    
    print("Estimated shares - coupled")
    print(Shares_Mat_Coupled)
    
    print("Estimated shares - singles")
    print(Shares_Mat_Singles)
    
    print(paste("Share of homogamous couples- Low Educated:", SHC$L, sep =" "))
    print(paste("Census Share of homogamous couples - Low Educated:", census_shc_l, sep=" "))
    print(paste("Share of homogamous couples - Middle Educated:", SHC$M, sep =" "))
    print(paste("Census Share of homogamous couples - Middle Educated:", census_shc_m, sep=" "))
    print(paste("Share of homogamous couples - High Educated:", SHC$H, sep =" "))
    print(paste("Census Share of homogamous couples - High Educated:", census_shc_h, sep=" "))
    print(paste("Share of homogamous couples - Aggregated:", SHC$Agg, sep =" "))
    print(paste("Census Share of homogamous couples - Aggregated:", census_shc_total, sep =" "))
    print(paste("Constraint check:", census_shc_total - SHC$Agg, sep = " "))
    print("Counterfactual contingency table")
    print(contingency)
  }
  
  
  #Objective function graphs
  # N = 100
  # par(mfrow = c(2, 3))
  # Vec_param = Vec_param_est
  # 
  # for (i in 1:6) {
  #   NSingles = seq(1, c(NDistr_F, NDistr_M)[i], length = N)
  #   B = matrix(Vec_param,
  #              nrow = N,
  #              ncol = 6,
  #              byrow = T)
  #   
  #   B[, i] = NSingles
  #   FunEvals = rep(NA, N)
  #   FunEvals = apply(B, 1, ObjectiveF)
  #   FunEvals = ifelse(FunEvals == 9e10, NA, FunEvals)
  #   plot(NSingles, FunEvals, type = "l")
  # }
  
  
  
  if (yearOfEduDistr == 2011 & counterfactual_year!=9999){
    print("Running Counterfactual 3 with estimated shares")
    
    # Estimated shares of coupled individuals
    f_LL= Shares_Mat_Coupled[1]
    f_LM= Shares_Mat_Coupled[2]
    f_LH= Shares_Mat_Coupled[3]
    f_ML= Shares_Mat_Coupled[4]
    f_MM= Shares_Mat_Coupled[5]
    f_MH= Shares_Mat_Coupled[6]
    f_HL= Shares_Mat_Coupled[7]
    f_HM= Shares_Mat_Coupled[8]
    f_HH= Shares_Mat_Coupled[9]
    m_LL= Shares_Mat_Coupled[10]
    m_LM= Shares_Mat_Coupled[11]
    m_LH= Shares_Mat_Coupled[12]
    m_ML= Shares_Mat_Coupled[13]
    m_MM= Shares_Mat_Coupled[14]
    m_MH= Shares_Mat_Coupled[15]
    m_HL= Shares_Mat_Coupled[16]
    m_HM= Shares_Mat_Coupled[17]
    m_HH= Shares_Mat_Coupled[18]
    
    
    # Estimated shares of single individuals
    f_LS = Shares_Mat_Singles[1]
    f_MS = Shares_Mat_Singles[2]
    f_HS = Shares_Mat_Singles[3]
    m_LS = Shares_Mat_Singles[4]
    m_MS = Shares_Mat_Singles[5]
    m_HS = Shares_Mat_Singles[6]
    
    
    # Building matrices of shares for females
    Shares_F=matrix(c(f_LL,f_ML,f_HL,
                      f_LM,f_MM,f_HM,
                      f_LH,f_MH,f_HH,
                      f_LS,f_MS,f_HS),nrow = 4,ncol = 3,byrow = T)
    
    
    # Building matrices of shares for males
    Shares_M=matrix(c(m_LL,m_LM,m_LH,m_LS,
                      m_ML,m_MM,m_MH,m_MS,
                      m_HL,m_HM,m_HH,m_HS),nrow = 3,ncol = 4,byrow = T)
    
    
    
    
    # Census data
    if (yearOfEduDistr_Block3 == 2000){
      NDistr_Couples = matrix(
        c(78480, 54100, 2620,
          78580, 312580, 44480,
          1960, 31140, 52420),
        nrow = dim,
        ncol = dim,
        byrow = T)
    } else {
      NDistr_Couples = matrix(
        c(231120, 67340, 360,
          194800, 332940, 8580,
          1840, 30100, 17080),
        nrow = dim,
        ncol = dim,
        byrow = T)
    }
    
    
    # Individuals in a relationship
    NDistr_F_inrelation = matrix(colSums(NDistr_Couples), nrow = 1, ncol = 3, byrow = T)
    NDistr_M_inrelation = matrix(colSums(t(NDistr_Couples)), nrow = 3, ncol = 1, byrow = T)
    
    
    # Singles
    # Females
    # Single Females
    if (yearOfEduDistr_Block3 == 2000){
      NDistr_F_singles = matrix(
        c(97800, 221480, 74780),
        nrow = 1,
        ncol = 3,
        byrow = T)
    } else {
      NDistr_F_singles = matrix(
        c(78100, 95680, 12440),
        nrow = 1,
        ncol = 3,
        byrow = T)
    }
    
    
    # Single Males
    if (yearOfEduDistr_Block3 == 2000){
      NDistr_M_singles = matrix(
        c(136000,
          315480,
          61420),
        nrow = 3,
        ncol = 1,
        byrow = T)
    } else {
      NDistr_M_singles = matrix(
        c(107420,
          128140,
          11100),
        nrow = 3,
        ncol = 1,
        byrow = T)
    }
    
    
    # Marginal distributions of Females and Males
    # Females
    if (is_with_singles == 1) {
      NDistr_F = NDistr_F_inrelation + NDistr_F_singles
    } else {
      NDistr_F = NDistr_F_inrelation
    }
    
    # Males
    if (is_with_singles == 1) {
      NDistr_M = NDistr_M_inrelation + NDistr_M_singles
    } else {
      NDistr_M = NDistr_M_inrelation
    }
    
    
    
    # Updating distributions with estimated shares
    # Females
    NPref_F_Mat=t(t(Shares_F)*as.vector(NDistr_F))
    
    # Males
    NPref_M_Mat=Shares_M*as.vector(NDistr_M)
    
    
    
    # Run Matching algorithm
    Match = Match_fun(NPref_F_Mat, NPref_M_Mat, is_with_singles)
    Matched_Mat = Match[[1]]
    Resid_M_Mat = Match[[2]]
    Resid_F_Mat = Match[[3]]
    
    contingency_C3 = matrix(0, nrow = 5, ncol = 5)
    
    contingency_C3[1:3, 1:3] = Matched_Mat
    contingency_C3[4, 1:3] = Resid_F_Mat[4, ]
    contingency_C3[1:3, 4] = Resid_M_Mat[, 4]
    contingency_C3[5, 1:3] = colSums(Resid_F_Mat[1:3, ])
    contingency_C3[1:3, 5] = rowSums(Resid_M_Mat[, 1:3])
    
    print("Counterfactual 3 contingency table")
    print(contingency_C3)
    

  }
  
  return(list(Value,
              constraint_bind,
              InvoluntarySingles,
              Voluntary_Singles,
              Matched_Mat,
              Shares_Mat_Coupled,
              Shares_Mat_Singles,
              contingency))
  
}



results = modified_gale_shapley(year_edu = 2011, dim = 3, counterfactual_year = 9999,
                      estimate_vol_singles = "yes", penalized = 40, supress_results = F,
                      supress_graphics = F)



# 2011 - penalty = 40
# 2000 - penalty = 120
# 1980 - penalty = 130


