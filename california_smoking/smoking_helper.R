#######################
##### Helpers #########
#######################
library(optimx)
library(rgenoud)
library(LowRankQP)

# helper function to get the data matrix for the desired SC.
smoking.sc.dataprep <- function(smoking.data, treatment.unit, control.units){
  
  preds = c("lnincome","beer","age15to24","retprice")
  spec.preds =  list(
    list("beer", 1984:1988, "mean"),
    list("cigsale", 1975, "mean"),
    list("cigsale", 1980, "mean"),
    list("cigsale", 1988, "mean"))
  
  dataprep.out <- dataprep(
    foo = smoking.data,
    predictors = preds,
    predictors.op = "mean",
    time.predictors.prior = 1970:1988,
    special.predictors = spec.preds,
    dependent = "cigsale",
    unit.variable = "stateno",
    #unit.names.variable = "state",
    time.variable = "year",
    treatment.identifier = treatment.unit,
    controls.identifier = control.units,
    time.optimize.ssr = 1970:1988,
    time.plot = 1970:2000)
  
  return( dataprep.out)
}

#######################
##### Functions from Synth Package #########
#######################


synth <-
  function(           data.prep.obj = NULL,
                      X1 = NULL,
                      X0 = NULL,
                      Z0 = NULL,
                      Z1 = NULL,
                      custom.v = NULL,
                      optimxmethod = c("Nelder-Mead","BFGS"),
                      genoud = FALSE,
                      quadopt = "ipop",
                      Margin.ipop = 0.0005,
                      Sigf.ipop = 5,
                      Bound.ipop = 10,
                      verbose = FALSE,
                      ...
  )
  { 
    # Retrieve dataprep objects
    if(is.null(data.prep.obj) == FALSE)
    {
      cat("\nX1, X0, Z1, Z0 all come directly from dataprep object.\n\n")
      X1 <- data.prep.obj$X1
      Z1 <- data.prep.obj$Z1
      X0 <- data.prep.obj$X0
      Z0 <- data.prep.obj$Z0
    } else {
      cat("X1,X0,Z1,Z0 were individually input (not dataprep object.)\n\n")
    }
    
    # routine checks
    store <- list(X1=X1,X0=X0,Z1=Z1,Z0=Z0)
    for(i in 1:4){
      if(is.null(store[[i]]))
      {stop(paste("\n",names(store)[i],"is missing \n"))}
      if(sum(is.na(store[[i]]))>0)
      {stop(paste("\n NAs in",names(store)[i],"\n"))}
      if(is.matrix(store[[i]]) == FALSE)
      {stop(paste("\n",names(store)[i],"is not a matrix object\n"))}
    }
    
    # geometry checks
    if(ncol(X1)!=1){stop("\n Please specify only one treated unit: X1 has to have ncol= 1")}
    if(ncol(Z1)!=1){stop("\n Please specify only one treated unit: Z1 has to have ncol= 1")}
    
    if(ncol(X0)<2){stop("\n Please specify at least two control units: X0 has to have ncol >= 2 ")}
    if(ncol(Z0)<2){stop("\n Please specify only one treated unit: Z0 has to have ncol >= 2")}
    
    if(nrow(Z0)!=nrow(Z1)){stop("\n Different number of periods for treated and controls: nrow(Z0) unequal nrow(Z1)")}
    if(nrow(X0)!=nrow(X1)){stop("\n Different number of predictors for treated and controls: nrow(X0) unequal nrow(X1)")}
    
    if(nrow(X0)==0){stop("No predictors specified. Please specify at least on predictor")}
    if(nrow(Z0)==0){stop("No periods specified for Z1 and Z0. Please specify at least on period")}
    
    if(0 %in% apply(X0,1,sd))
    {stop("\n At least one predictor in X0 has no variation across control units. Please remove this predictor.")}
    
    # collinearity check
    # check <- try(solve(t(X0)%*%X0),silent=TRUE)
    # if(class(check)=="try-error")
    #  {stop("\n Some of the predictors in X0 are collinear (t(X0)%*%X0) not invertible")}
    
    # Normalize X
    nvarsV <- dim(X0)[1]
    big.dataframe <- cbind(X0, X1)
    divisor <- sqrt(apply(big.dataframe, 1, var))
    scaled.matrix <-
      t(t(big.dataframe) %*% ( 1/(divisor) *
                                 diag(rep(dim(big.dataframe)[1], 1)) ))
    
    X0.scaled <- scaled.matrix[,c(1:(dim(X0)[2]))]
    if(is.vector(X0.scaled)==TRUE)
    {X0.scaled <- t(as.matrix(X0.scaled))}
    X1.scaled <- scaled.matrix[,dim(scaled.matrix)[2]]
    
    
    # check if custom v weights are supplied or
    # if only on predictor is specified,
    # we jump to quadratic optimization over W weights
    # if not start optimization over V and W
    if(is.null(custom.v) & nrow(X0) != 1)
    {
      
      # two attemps for best V are made: 
      # equal weights and regression based starting values 
      cat("\n****************",
          "\n searching for synthetic control unit  \n","\n"
      )
      
      if(genoud == TRUE) # if user wants genoud as well
      {
        # we run genoud first
        cat("\n****************",
            "\n genoud() requested for optimization\n","\n"
        )
        
        rgV.genoud <- rgenoud::genoud(
          fn.V, 
          nvarsV, 
          X0.scaled = X0.scaled,
          X1.scaled = X1.scaled,
          Z0 = Z0,
          Z1 = Z1,
          quadopt = quadopt,
          margin.ipop = Margin.ipop,
          sigf.ipop = Sigf.ipop,
          bound.ipop = Bound.ipop
        )
        SV1 <- rgV.genoud$par  # and use these as starting values
        
        cat("\n****************",
            "\n genoud() finished, now running local optimization using optim()\n","\n"
        )
        
      } else {
        # if we don't use genoud first: set of starting values: equal weights
        SV1 <- rep(1/nvarsV,nvarsV)
      }
      
      # now we run optimization
      all.methods <- FALSE
      if(sum(optimxmethod %in% c("All"))==1){ all.methods <- TRUE }
      rgV.optim.1 <- optimx(par=SV1, fn=fn.V,
                            gr=NULL, hess=NULL, 
                            method=optimxmethod, itnmax=NULL, hessian=FALSE,
                            control=list(kkt=FALSE,
                                         starttests=FALSE,
                                         dowarn=FALSE,
                                         all.methods=all.methods),
                            X0.scaled = X0.scaled,
                            X1.scaled = X1.scaled,
                            Z0 = Z0,
                            Z1 = Z1,
                            quadopt = quadopt,
                            margin.ipop = Margin.ipop,
                            sigf.ipop = Sigf.ipop,
                            bound.ipop = Bound.ipop
      )
      # get minimum
      if(verbose==TRUE){print(rgV.optim.1)}
      rgV.optim.1 <- collect.optimx(rgV.optim.1,"min")
      
      # second set of starting values: regression method 
      # will sometimes not work because of collinear Xs
      # so it's wrapped in a try command
      Xall <- cbind(X1.scaled,X0.scaled)
      Xall <- cbind(rep(1,ncol(Xall)),t(Xall))
      Zall <- cbind(Z1,Z0)
      Beta <- try(solve(t(Xall)%*%Xall)%*%t(Xall)%*%t(Zall),silent=TRUE)
      
      # if inverses did not work, we
      # stick with first results    
      if(inherits(Beta,"try-error")) 
      {
        rgV.optim <- rgV.optim.1
      } else {
        # otherwise we run a second optimization with regression starting values
        Beta <- Beta[-1,]
        V    <- Beta%*%t(Beta)
        SV2  <- diag(V)
        SV2 <- SV2 / sum(SV2)
        
        rgV.optim.2 <- optimx(par=SV2, fn=fn.V,
                              gr=NULL, hess=NULL, 
                              method=optimxmethod, itnmax=NULL, hessian=FALSE,
                              control=list(kkt=FALSE,
                                           starttests=FALSE,
                                           dowarn=FALSE,
                                           all.methods=all.methods),
                              X0.scaled = X0.scaled,
                              X1.scaled = X1.scaled,
                              Z0 = Z0,
                              Z1 = Z1,
                              quadopt = quadopt,
                              margin.ipop = Margin.ipop,
                              sigf.ipop = Sigf.ipop,
                              bound.ipop = Bound.ipop
        )
        # get minimum
        if(verbose==TRUE){print(rgV.optim.2)}
        rgV.optim.2 <- collect.optimx(rgV.optim.2,"min")
        
        # ouput
        if(verbose == TRUE){
          cat("\n Equal weight loss is:",rgV.optim.1$value,"\n")
          cat("\n Regression Loss is:",rgV.optim.2$value,"\n")
        }       
        # and keep the better optim results    
        if(rgV.optim.1$value < rgV.optim.2$value) 
        {
          rgV.optim <- rgV.optim.1
        } else {
          rgV.optim <- rgV.optim.2
        }
      } # close if statement for second regression based optimization attempt
      
      # final V weights from optimization
      solution.v   <- abs(rgV.optim$par)/sum(abs(rgV.optim$par))
    } else { # jump here if only optimize over W
      
      
      cat("\n****************",
          "\n optimization over w weights: computing synthtic control unit \n","\n\n"
      )
      
      if(nrow(X0)==1)
      {
        custom.v <- 1 # this is the case where only one predictor is specified V is the identity matrix
      } else {
        # if user choose to supply v manually:
        cat("\n****************",
            "\n v weights supplied manually: computing synthtic control unit \n","\n\n"
        )
        if(length(custom.v) != nvarsV){stop("custom.V misspecified: length(custom.V) != nrow(X1)")}
        if(mode(custom.v)!="numeric"){stop("custom.V must be numeric")}
      }
      
      # enter solution.V
      rgV.optim  <- NULL
      solution.v <- abs(custom.v)/sum(custom.v)
      
    } # close else statment for by-passing V optimization
    
    
    # last step: now recover solution.w
    V <- diag(x=as.numeric(solution.v),nrow=nvarsV,ncol=nvarsV)
    H <- t(X0.scaled) %*% V %*% (X0.scaled)
    a <- X1.scaled
    c <- -1*c(t(a) %*% V %*% (X0.scaled) )
    A <- t(rep(1, length(c)))
    b <- 1
    l <- rep(0, length(c))
    u <- rep(1, length(c))
    r <- 0
    
    if(quadopt=="ipop"){
      res <- ipop(c = c, H = H, A = A, b = b, l = l, u = u, r = r,
                  margin = Margin.ipop, maxiter = 1000, sigf = Sigf.ipop, bound = Bound.ipop)
      solution.w <- as.matrix(primal(res))
    } else {
      # LowRankQP
      if(quadopt=="LowRankQP"){
        res <- LowRankQP::LowRankQP(Vmat=H,dvec=c,Amat=A,bvec=1,uvec=rep(1,length(c)),method="LU")
        solution.w <- as.matrix(res$alpha)
      } 
    }
    
    rownames(solution.w) <- colnames(X0)
    colnames(solution.w) <- "w.weight"
    names(solution.v) <- rownames(X0)
    
    loss.w <- t(X1.scaled - X0.scaled %*% solution.w) %*%
      V %*% (X1.scaled - X0.scaled %*% solution.w)
    
    loss.v <-
      t(Z1 - Z0 %*% as.matrix(solution.w)) %*%
      (Z1 - Z0 %*% as.matrix(solution.w)) 
    loss.v <- loss.v/nrow(Z0)      
    
    # produce viewable output
    cat("\n****************",
        "\n****************",
        "\n****************",
        "\n\nMSPE (LOSS V):", loss.v,
        #      "\n\nLOSS (W):", loss.w,
        "\n\nsolution.v:\n", round(as.numeric(solution.v), 10),
        "\n\nsolution.w:\n", round(as.numeric(solution.w), 10),
        "\n\n"
    )
    
    optimize.out <- list(
      solution.v = solution.v,
      solution.w = solution.w,
      loss.v = loss.v,
      loss.w = loss.w,
      custom.v = custom.v,
      rgV.optim = rgV.optim
    )
    
    return(invisible(optimize.out))
    
  }



dataprep <-
  function(
    foo = NULL,
    predictors = NULL,
    predictors.op = "mean",
    special.predictors = NULL,
    dependent = NULL,
    unit.variable = NULL,
    time.variable = NULL,
    treatment.identifier = NULL,
    controls.identifier = NULL,
    time.predictors.prior = NULL,
    time.optimize.ssr = NULL,
    time.plot = time.optimize.ssr,
    unit.names.variable = NA
  )

  {
    
    if(is.data.frame(foo) == FALSE ){stop("\n No data.frame supplied in foo.\n")}
    
    # identify unit.variable
    if(mode(unit.variable) == "character"){unit.variable<-which(names(foo)==unit.variable)}
    if(is.null(unit.variable) == TRUE ||  mode(foo[,unit.variable])!="numeric")
    { stop("\n unit.variable not found as numeric variable in foo.\n")}
    if(length(unit.variable)!=1){stop("\ Please specify only one unit.variable\n")}
    
    # identify time.variable
    if(mode(time.variable) == "character"){time.variable<-which(names(foo)==time.variable)}
    if(is.null(time.variable) == TRUE || mode(foo[,time.variable])!="numeric" )
    {
      stop("\n time.variable not found as numeric variable in foo.\n")
    }
    if(length(time.variable)!=1){stop("\ Please specify only one time.variable\n")}
    # identify units
    # check if unit.name var is there (required if any identifier is given as character)
    if(
      mode(treatment.identifier) == "character" ||
      mode(controls.identifier)  == "character" ||
      is.na(unit.names.variable) == FALSE
    )
    {
      if(mode(unit.names.variable) == "character")
      {unit.names.variable<-which(names(foo)==unit.names.variable)}
      if(is.na(unit.names.variable) == TRUE || mode(foo[,unit.names.variable])!="character")
      {stop("\n unit.names.variable not found as character variable in foo.\n")}
    }
    # check treated
    if(length(treatment.identifier) != 1){stop("\n please specify a single treated unit\n")}
    if(mode(treatment.identifier) == "character")
    {
      if(treatment.identifier %in% foo[,unit.names.variable] == FALSE  )
      {stop("\n treated unit not found in unit.names.variable.\n")}
      tr.id = unique(foo[foo[,unit.names.variable] == treatment.identifier,unit.variable])
      if(length(tr.id)!=1)
      {stop(paste("\n treatment unit name",treatment.identifier,"refers to multiple numbers in unit.variable"))}
      treatment.identifier.name = treatment.identifier
      treatment.identifier = tr.id
    } else {
      if(treatment.identifier %in% foo[,unit.variable] == FALSE )
      {stop("\n treated unit not found in unit.variable.\n")}
      if(is.na(unit.names.variable) == TRUE )
      {
        treatment.identifier.name = treatment.identifier
      } else {
        treatment.identifier.name = unique(foo[treatment.identifier==foo[,unit.variable],unit.names.variable])
        if(length(treatment.identifier.name)>1)
        {stop("\n treatment.identifier refers to multiple names in unit.names.variable")}
      }
    }
    
    # check controls
    if(length(controls.identifier) < 2 ){stop("\n please specify at least two control units\n")}
    if(sum(duplicated(controls.identifier))>0)
    {stop("\n duplicate control units in controls.identifier\n")}
    if(mode(controls.identifier) == "character")
    { co.store <- c()
    for(i in controls.identifier)
    {
      if(i %in% foo[,unit.names.variable] == FALSE )
      {stop(paste("\n control unit",i,"not found in unit.names.variable"))}
      co.id = unique(foo[foo[,unit.names.variable] == i,unit.variable])
      if(length(co.id) != 1)
      {stop(paste("\n control unit name",i," refers to multiple numbers in unit.variable"))}
      co.store <- c(co.store,co.id)
    }
    controls.identifier.name = controls.identifier
    controls.identifier = co.store
    } else {
      co.store <- c()
      for(i in controls.identifier)
      {
        if(i %in% foo[,unit.variable] == FALSE )
        {stop(paste("\n control unit",i,"not found in unit.variable"))}
        if(is.na(unit.names.variable) == FALSE )
        {
          co.id = unique(foo[foo[,unit.variable] == i,unit.names.variable])
          if(length(co.id) != 1)
          {stop(paste("\n control unit number",i," refers to multiple names in unit.names.variable"))}
          co.store <- c(co.store,co.id)
        }
      }
      if(is.na(unit.names.variable) == FALSE )
      {
        controls.identifier.name = co.store
      } else {
        controls.identifier.name = controls.identifier
      }
    }
    
    # more checks
    if(treatment.identifier.name %in% controls.identifier.name)
    {stop("\n treated unit among controls\n")}
    
    if(sum(duplicated(c(controls.identifier.name,treatment.identifier.name))) > 0)
    {stop("n duplicate unit.variable.names across units\n")}
    
    # sort first by unit, then by time variable
    foo[,time.variable] <- as.numeric(as.character(foo[,time.variable]))
    foo[,unit.variable] <- as.numeric(as.character(foo[,unit.variable]))
    foo <- foo[order(foo[,unit.variable], foo[,time.variable]),]
    
    # Get rows
    treatment.rows <- which(foo[,unit.variable] %in% treatment.identifier)
    control.rows   <- which(foo[,unit.variable] %in% controls.identifier)
    
    # check if panel is unbalanced  (possible cases where this too restrictive, but imposes discipline)
    balcheck <-       table(    foo[c(control.rows,treatment.rows),unit.variable],
                                foo[c(control.rows,treatment.rows),time.variable])         
    if( length(unique(balcheck)) != 1 || unique(balcheck)!= 1)
    {
      stop("\n Your panel, as described by unit.variable and time.variable, is unbalanced. Balance it and run again.\n")
    }
    
    # Now check and get time identifiers
    if(sum(is.null(time.predictors.prior))> 0)
    {stop("time.predictors.prior missing")}
    if(sum(is.null(time.optimize.ssr))>0)
    {stop("time.optimize.ssr missing")}
    if(sum(is.null(time.plot))>0)
    {stop("time.plot missing")}
    
    t.list <- list(time.predictors.prior,time.optimize.ssr,time.plot)
    names(t.list) <- c("predictors.prior","optimize.ssr","plot")
    for(i in 1:3)
    {
      if(mode(t.list[[i]]) != "numeric")
      {stop(paste("\n time",names(t.list[i]),"must be numeric\n"))}
      if(sum(duplicated(t.list[[i]]))>0)
      {stop(paste("\n duplicates in time",names(t.list[i]),"\n"))}
      if(length(t.list[[i]]) < 1)
      {stop(paste("\n specificy at least one period in time",names(t.list[i]),"\n"))}
      for(p in t.list[[i]])
      {
        if(p %in% unique(foo[,time.variable]) == FALSE)
          stop(paste("\n time period ",p," from time.",names(t.list[i])," not found in time.variable\n",sep=""))
      }
    }
    
    # get time rows
    time.predictors.prior.rows <-
      which(foo[,time.variable] %in% time.predictors.prior)
    
    time.optimize.ssr.rows <-
      which(foo[,time.variable] %in% time.optimize.ssr)
    
    # Build X1 & X0 Predictor matrices
    
    # first, run check on regular predictors, if specified
    if(is.null(predictors)==FALSE)
    {
      # predictor checks
      if(sum(duplicated(predictors))>0)
      {stop("\n duplicates found in predictors \n")}
      
      for(i in predictors)
      {
        if(mode(foo[,i]) != "numeric")
        {stop(paste("\n predictor",i,"not found as numeric variable in foo \n"))}
      }
      
      if(mode(predictors)=="character")
      {
        pred.no <- c()
        for(i in 1:length(predictors)){pred.no <- c(pred.no,which(names(foo) == predictors[i]))}
        predictors <- pred.no
      }
      # X1 matrix for treated
      X1 <-
        data.frame(foo[
          intersect(
            treatment.rows,
            time.predictors.prior.rows
          ),
          predictors
        ]
        )
      colnames(X1) <- names(foo)[predictors]
      rownames(X1) <- time.predictors.prior
      
      # deep missing checker:
      for(i in 1:ncol(X1))
      {
        if(sum(is.na(X1[,i])) == nrow(X1))
        {stop(paste("\n predictor",names(X1)[i],"has missing data for all periods in time.predictors.prior\n"))}
        
        for(j in 1:nrow(X1))
        {
          if(is.na(X1[j,i])){
            cat(paste("\n Missing data- treated unit; predictor:",names(X1)[i],"; for period:",rownames(X1)[j],
                      "\n","We ignore (na.rm = TRUE) all missing values for predictors.op.\n"))}
        }
      }
      
      # aggregate
      X1 <- apply(X1, 2, paste(predictors.op), na.rm = TRUE)
      X1 <- as.matrix(X1)
      rownames(X1) <- names(foo)[predictors]
      colnames(X1) <- treatment.identifier
      
    } else { # if no regular predictors are specified, pass a void matrix to special predictors
      
      X1 <- matrix(NA,0,1)
      colnames(X1) <- treatment.identifier 
    }
    
    # X0 matrix for controls
    if(is.null(predictors)==FALSE)
    {
      X0 <- data.frame(foo[intersect(control.rows,
                                     time.predictors.prior.rows
      ),
      c(predictors, unit.variable)
      ]
      )
      checknames <- rep(time.predictors.prior,length(controls.identifier))
      
      # deep missing checker:
      for(i in 1:(ncol(X0)-1))
      {
        
        for(p in controls.identifier)
        {
          if(sum(is.na(X0[X0[,ncol(X0)] == p,i])) == length(is.na(X0[X0[,ncol(X0)] == p,i])))
          {stop(paste("\n control unit:",p,"; predictor:",names(X0)[i],"has missing data for all periods in time.predictors.prior\n"))}
        }
        for(j in 1:nrow(X1))
        {
          if(is.na(X0[j,i])){
            cat(paste("\n Missing data - control unit:",X0[j,ncol(X0)],"; predictor:",names(X0)[i],"; for period:",checknames[j],
                      "\n","We ignore (na.rm = TRUE) all missing values for predictors.op.\n"))}
        }
      }
      
      X0 <- split(X0, X0[,dim(X0)[2]])
      X0 <- sapply(X0, apply, 2, mean, na.rm = TRUE, simplify = TRUE)
      X0 <- as.matrix(X0[-dim(X0)[1],])
      
      
      # Take transpose in presence of a single predictor only
      if(dim(X0)[2]==1)
      {
        X0 <- t(X0)
        rownames(X0) <- names(foo)[predictors]
        colnames(X0) <- controls.identifier
      }
      
    } else { # if no regular predictors are specified, pass a void matrix to special predictors
      X0 <- matrix(NA,0,length(controls.identifier))
      #rownames(X0) <- names(foo)[predictors]
      colnames(X0) <- controls.identifier
    }
    
    # Add Special Predictors to X1 and X0
    
    if(is.null(special.predictors) == FALSE)
    {
      if(is.list(special.predictors) == FALSE)
      {stop("\nspecial.predictors is not a list object\n")}
      
      for(i in 1:length(special.predictors))
      {
        
        # checks for special predictors
        if(is.list(special.predictors[[i]]) == FALSE)
        {stop(paste("\n special.predictor number",i,"is not a list object\n"))}
        
        if(length(special.predictors[[i]]) != 3)
        {stop(paste("\n special.predictor number",i,"is misspecified (requires predictor name, time periods, and name of operation\n"))}
        
        # name check
        sp.name <- special.predictors[[i]][[1]]
        if(is.na(sp.name) || length(sp.name) != 1)
        {stop(paste("\n predictor name",sp.name,"of special.predictor number",i,"misspecified\n"))}
        
        # availability check
        if(mode(foo[,sp.name]) != "numeric")
        {stop(paste("\n special predictor named ",sp.name,"not found as numeric variable in foo \n"))}
        
        # time check
        sp.time <- special.predictors[[i]][[2]]
        if(mode(sp.time) != "numeric")
        {stop(paste("\n for special predictor:",sp.name," (which is special.predictor number",i,") the time period is misspecified\n"))}
        if(sum(duplicated(sp.time))>0)
        {stop(paste("\n for special predictor:",sp.name," (which is special.predictor number",i,") time period contains duplicates\n"))}
        if(length(sp.time)<1)
        {stop(paste("\n for special predictor:",sp.name," (which is special.predictor number",i,") specify at least on time period\n"))}
        
        # time availability check
        for(p in sp.time)
        {
          if(p %in% unique(foo[,time.variable]) == FALSE)
          {stop(paste("\n for special predictor:",sp.name," (which is special.predictor number",i,") time period",p,"not found in time.variable\n"))}
        }
        
        # operator check
        sp.op <- special.predictors[[i]][[3]]
        if(mode(sp.op) != "character" || length(sp.op) !=1 )
        {stop(paste("\n for special predictor:",sp.name," (which is special.predictor number",i,") the operator is mispecified\n"))}
        
        # now go and built matrices
        spf <- spec.pred.func(list.object = special.predictors[[i]],
                              tr.numb = treatment.identifier,
                              co.numb = controls.identifier,
                              unit.var = unit.variable,
                              time.var = time.variable,
                              foo.object = foo,
                              X0.inner = X0,
                              X1.inner = X1
        )
        
        X0 <- spf[[1]]
        X1 <- spf[[2]]
        
      }
      
    }
    
    # no predictors check
    if(nrow(X0)==0)
    {stop("No predictors specified. Please specify at least on predictor")}
    
    # Dependent Variable Pretreatment
    
    if(is.null(dependent))
    {stop("\n dependent variabale is missing")}
    
    if(mode(foo[,dependent]) != "numeric")
    {stop(paste("\n dependent variable",dependent,"not found as numeric variable in foo \n"))}
    
    if(mode(dependent) == "character")
    {dependent <- which(is.element(names(foo), dependent))}
    #foo[,dependent] <- as.numeric(as.character(foo[,dependent]))
    
    Z1 <- foo[
      intersect(
        treatment.rows,
        time.optimize.ssr.rows
      ),
      dependent
    ]
    
    Z1 <- as.matrix(Z1)
    
    # missing check
    if(sum(is.na(Z1)) == length(Z1))
    {stop("\n treated unit: dependent variable is missing for all periods in time.optimize.ssr.rows \n")}
    
    rownames(Z1) <- time.optimize.ssr
    colnames(Z1) <- treatment.identifier
    
    for(i in 1:length(Z1))
    {
      if(is.na(Z1[i]))
      {stop(paste("\n treated unit: dependent variable is missing for period",rownames(Z1)[i],"in time.optimize.ssr.rows \n"))}
    }
    
    Z0 <-
      as.matrix(foo[
        intersect(
          control.rows,
          time.optimize.ssr.rows
        ),
        c(dependent)
      ]
      )
    
    Z0 <- matrix(
      Z0,
      byrow = FALSE,
      nrow = length(time.optimize.ssr),
      ncol = length(controls.identifier),
    )
    
    colnames(Z0)<- controls.identifier
    rownames(Z0)<- time.optimize.ssr
    
    for(i in 1:ncol(Z0))
    {
      for(j in 1:nrow(Z0))
      {
        if(is.na(Z0[j,i]))
        {stop(paste("\n control unit:",colnames(Z0)[i],"; dependent variable",dependent,"is missing for period",rownames(Z0)[j],"in time.optimize.ssr.rows \n"))}
      }
    }
    
    
    # dependent Variable Plot
    time.plot.rows <-
      which(is.element(foo[,time.variable], time.plot ) == TRUE)
    
    Y1plot <- foo[intersect(treatment.rows, time.plot.rows), dependent]
    Y1plot <- as.matrix(Y1plot)
    rownames(Y1plot) <- time.plot
    colnames(Y1plot) <- treatment.identifier
    
    if(sum(is.na(Y1plot)) == length(Y1plot))
    {stop("\n treated unit: dependent variable is missing for all periods in time.plot \n")}
    
    for(i in 1:length(Y1plot))
    {
      if(is.na(Y1plot[i]))
      {stop(paste("\n treated unit: dependent variable is missing for period",rownames(Y1plot)[i],"in time.plot \n"))}
    }
    
    Y0plot <- as.matrix(foo[
      intersect(control.rows, time.plot.rows),
      c(dependent)
    ]
    )
    
    Y0plot  <- matrix(
      Y0plot,
      byrow = FALSE,
      nrow = length(time.plot),
      ncol = length(controls.identifier),
    )
    
    
    rownames(Y0plot) <- time.plot
    colnames(Y0plot) <- controls.identifier
    
    
    # table with unit names
    names.and.numbers <-
      data.frame(c(treatment.identifier.name,controls.identifier.name),
                 c(treatment.identifier,controls.identifier))
    
    
    
    names(names.and.numbers) <- c(
      "unit.names",
      "unit.numbers"
    )
    
    
    
    #######################
    
    tag <- list(
      foo = as.character(foo),
      predictors = predictors,
      predictors.op = predictors.op,
      special.predictors = special.predictors,
      dependent = dependent,
      unit.variable = unit.variable,
      time.variable = time.variable,
      treatment.identifier = treatment.identifier,
      controls.identifier = controls.identifier,
      time.predictors.prior = time.predictors.prior,
      time.optimize.ssr = time.optimize.ssr,
      time.plot = time.plot,
      unit.names.variable = unit.names.variable
    )
    
    ######################################
    
    output <- list(
      X0 = X0,
      X1 = X1,
      Z0 = Z0,
      Z1 = Z1,
      Y0plot = Y0plot,
      Y1plot = Y1plot,
      names.and.numbers = names.and.numbers,
      tag = tag
    )
    
    return(invisible(output))
    
  }


spec.pred.func <-
  function(
    list.object = NULL,
    tr.numb = NULL,
    co.numb = NULL,
    unit.var = NULL,
    time.var = NULL,
    foo.object = NULL,
    X0.inner = NULL,
    X1.inner = NULL
  )
  {
    
    special.units.tr <-
      which(is.element(foo.object[,unit.var], tr.numb) == TRUE)
    
    special.units.co <-
      which(is.element(foo.object[,unit.var], co.numb) == TRUE)
    
    special.times <-
      which(is.element(foo.object[,time.var], list.object[[2]]) == TRUE)
    
    if(mode(list.object[[1]]) == "character")
    {
      list.object[[1]] <- which(names(foo.object) == list.object[[1]])
    }
    
    # create predictor name:
    name.predictor <- names(foo.object)[list.object[[1]]]
    if(length(list.object[[2]])>1)
    {
      name.predictor <-
        paste("special.", name.predictor, ".",
              paste(list.object[[2]][1],".",
                    list.object[[2]][length(list.object[[2]])]
                    , sep = ""), sep = "")
      
    } else {
      
      name.predictor <-
        paste("special.", name.predictor, ".",
              paste(list.object[[2]]
                    , sep = ""), sep = "")
    }
    
    X1.special <- as.matrix(foo.object[
      intersect(
        special.units.tr,
        special.times
      ),
      list.object[[1]]
    ]
    )
    
    rownames(X1.special) <-
      foo.object[
        intersect(special.units.tr, special.times),
        time.var
      ]
    colnames(X1.special) <- tr.numb
    
    # deep missing checker:
    for(i in 1:ncol(X1.special))
    {
      if(sum(is.na(X1.special[,i])) == nrow(X1.special))
      {stop(paste("\n special predictor",name.predictor,"has missing data for all periods in time.predictors.prior\n"))}
      
      for(j in 1:nrow(X1.special))
      {
        if(is.na(X1.special[j,i])){
          cat(paste("\n Missing data: treated unit; special predictor:",name.predictor,"; for period:",rownames(X1.special)[j],
                    "\n","We ignore (na.rm = TRUE) all missing values for predictors.op.\n"))}
      }
    }
    
    # aggregate if only one time period is there
    if(length(list.object[[2]]) > 1)
    {
      X1.special <- apply(X1.special,
                          2,
                          paste(list.object[[3]]),
                          na.rm = TRUE
      )
    }
    X1.special <- t(as.matrix(X1.special))
    X1.inner <- rbind(X1.inner, X1.special)
    rownames(X1.inner)[nrow(X1.inner)] <- name.predictor
    
    # controls
    X0.special <-
      as.matrix(foo.object[intersect(
        special.units.co,
        special.times
      ),
      c(list.object[[1]])
      ]
      )
    
    X0.special <-
      matrix(X0.special[,1], byrow = FALSE, ncol = length(co.numb))
    
    # Define row and column names
    
    row.names(X0.special) <-
      foo.object[intersect(special.units.tr, special.times), time.var]
    
    colnames(X0.special) <- co.numb
    
    # deep missing checker:
    for(i in 1:ncol(X0.special))
    {
      if(sum(is.na(X0.special[,i])) == nrow(X0.special))
      {stop(paste("\n control unit:",co.numb[i],"; special predictor:",name.predictor,"has missing data for all periods specified:",list.object[[2]],"\n"))}
      
      for(j in 1:nrow(X0.special))
      {
        if(is.na(X0.special[j,i])){
          cat(paste("\n Missing data - control unit:",co.numb[i],"; special predictor:",name.predictor,"; for period:",rownames(X0.special)[j],
                    "\n","We ignore (na.rm = TRUE) all missing values for predictors.op.\n"))}
      }
    }
    
    # Continue with object creation
    if(length(list.object[[2]]) > 1)
    {
      X0.special <- apply(
        X0.special,
        2,
        paste(list.object[[3]]),
        na.rm = TRUE
      )
      
      X0.special <- t(as.matrix(X0.special))
    }
    
    X0.inner <- rbind(X0.inner,X0.special)
    rownames(X0.inner)[nrow(X0.inner)] <- name.predictor
    
    
    
    
    # Prepare output
    special.output <- list(X0.inner = X0.inner,
                           X1.inner = X1.inner
    )
    
    return(invisible(special.output))
    
  }

fn.V <-
  function(
    variables.v = stop("variables.v missing"),
    X0.scaled = stop("X0.scaled missing"),
    X1.scaled = stop("X1.scaled missing"),
    Z0 = stop("Z0 missing"),
    Z1 = stop("Z1 missing"),
    margin.ipop = 0.0005,
    sigf.ipop = 5,
    bound.ipop = 10,
    quadopt = "ipop"
  )
    
  {
    
    # check quadopt
    Check <- sum(quadopt %in% c("ipop","LowRankQP"))
    if(Check!=1){
      stop("option quadopt must be one of ipop or LowRankQP") 
    }
    
    # rescale par
    #V <- diag( abs(variables.v)/sum(abs(variables.v)) )
    V <- diag(x=as.numeric(abs(variables.v)/sum(abs(variables.v))),
              nrow=length(variables.v),ncol=length(variables.v))
    
    # set up QP problem
    H <- t(X0.scaled) %*% V %*% (X0.scaled)
    a <- X1.scaled
    c <- -1*c(t(a) %*% V %*% (X0.scaled) )
    A <- t(rep(1, length(c)))
    b <- 1
    l <- rep(0, length(c))
    u <- rep(1, length(c))
    r <- 0
    
    # run QP and obtain w weights
    # ipop
    if(quadopt=="ipop"){
      res <- ipop(c = c, H = H, A = A, b = b, l = l, u = u, r = r, bound = bound.ipop,
                  margin = margin.ipop, maxiter = 1000, sigf = sigf.ipop)
      solution.w <- as.matrix(primal(res))
    } else {
      # LowRankQP
      if(quadopt=="LowRankQP"){
        res <- LowRankQP::LowRankQP(Vmat=H,dvec=c,Amat=A,bvec=1,uvec=rep(1,length(c)),method="LU")
        solution.w <- as.matrix(res$alpha)
      } 
    }
    
    # compute losses
    loss.w <- as.numeric(t(X1.scaled - X0.scaled %*% solution.w) %*%
                           (V) %*% (X1.scaled - X0.scaled %*% solution.w))
    
    loss.v <- as.numeric(t(Z1 - Z0 %*% solution.w) %*%
                           ( Z1 - Z0 %*% solution.w ))
    loss.v <- loss.v/nrow(Z0)
    
    return(invisible(loss.v))
  }

collect.optimx <-
  function(res,opt="min"){ 
    
    # reorder dataframe for best solution
    if(opt=="min"){
      res <- res[order(res$value,decreasing=FALSE),]
    }
    if(opt=="max"){
      res <- res[order(res$value,decreasing=TRUE),]
    }
    
    # index for value
    val <- which(colnames(res)=="value") 
    # index for pars
    pars <- 1:(val-1)
    
    # extract best solution
    out <- list(out.list=res,
                par=res[1,pars],
                value=res[1,val]
    )
    return(invisible(out))
    
  }

##ipop solves the quadratic programming problem
##minimize   c' * primal + 1/2 primal' * H * primal
##subject to b <= A*primal <= b + r
##           l <= x <= u
##           d is the optimizer itself
##returns primal and dual variables (i.e. x and the Lagrange
##multipliers for b <= A * primal <= b + r)
##for additional documentation see
##     R. Vanderbei
##     LOQO: an Interior Point Code for Quadratic Programming, 1992
## Author:      R version Alexandros Karatzoglou, orig. matlab Alex J. Smola
## Created:     12/12/97
## R Version:   12/08/03
## Updated:     13/10/05
## This code is released under the GNU Public License



setGeneric("ipop",function(c, H, A, b, l, u, r, sigf=7, maxiter=40, margin=0.05, bound=10, verb=0) standardGeneric("ipop"))
setMethod("ipop",signature(H="matrix"),
          function(c, H, A, b, l, u, r, sigf=7, maxiter=40, margin=0.05, bound=10, verb=0)
          {
            
            if(!is.matrix(H)) stop("H must be a matrix")
            if(!is.matrix(A)&&!is.vector(A)) stop("A must be a matrix or a vector")
            if(!is.matrix(c)&&!is.vector(c)) stop("c must be a matrix or a vector")
            if(!is.matrix(l)&&!is.vector(l)) stop("l must be a matrix or a vector")
            if(!is.matrix(u)&&!is.vector(u)) stop("u must be a matrix or a vector")
            
            n <- dim(H)[1]
            
            ## check for a decomposed H matrix
            if(n == dim(H)[2])
              smw <- 0
            if(n > dim(H)[2])
              smw <- 1
            if(n < dim(H)[2])
            {
              smw <- 1
              n <- dim(H)[2]
              H <- t(H)
            }
            
            if (is.vector(A)) A <- matrix(A,1)
            m <- dim(A)[1]
            primal <- rep(0,n)
            if (missing(b))
              bvec <- rep(0, m)
            ## if(n !=nrow(H))
            ##   stop("H matrix is not symmetric")
            if (n != length(c))
              stop("H and c are incompatible!")
            if (n != ncol(A))
              stop("A and c are incompatible!")
            if (m != length(b))
              stop("A and b are incompatible!")
            if(n !=length(u))
              stop("u is incopatible with H")
            if(n !=length(l))
              stop("l is incopatible with H")
            
            c <- matrix(c)
            l <- matrix(l)
            u <- matrix(u)
            
            m <- nrow(A)
            n <- ncol(A)
            H.diag <- diag(H)
            if(smw == 0)
              H.x <- H
            else if (smw == 1)
              H.x <- t(H)
            b.plus.1 <- max(svd(b)$d) + 1
            c.plus.1 <- max(svd(c)$d) + 1
            one.x <- -matrix(1,n,1)
            one.y <- -matrix(1,m,1)
            ## starting point
            if(smw == 0)
              diag(H.x) <- H.diag + 1
            else
              smwn <- dim(H)[2]
            H.y <- diag(1,m)
            c.x <- c
            c.y <- b
            ## solve the system [-H.x A' A H.y] [x, y] = [c.x c.y]
            if(smw == 0)
            {
              AP <- matrix(0,m+n,m+n)
              xp <- 1:(m+n) <= n
              AP[xp,xp] <- -H.x
              AP[xp == FALSE,xp] <- A
              AP[xp,xp == FALSE] <- t(A)
              AP[xp == FALSE, xp== FALSE] <- H.y
              s.tmp <- solve(AP,c(c.x,c.y))
              x <- s.tmp[1:n]
              y <- s.tmp[-(1:n)]
            }
            else
            {
              V <- diag(smwn)
              smwinner <- chol(V + crossprod(H))
              smwa1 <- t(A)
              smwc1 <- c.x
              smwa2 <- smwa1 - (H %*% solve(smwinner,solve(t(smwinner),crossprod(H,smwa1))))
              smwc2 <- smwc1 - (H %*% solve(smwinner,solve(t(smwinner),crossprod(H,smwc1)))) 
              y <- solve(A %*% smwa2 + H.y , c.y + A %*% smwc2)
              x <- smwa2 %*% y - smwc2
            }
            
            g <- pmax(abs(x - l), bound)
            z <- pmax(abs(x), bound)
            t <- pmax(abs(u - x), bound)
            s <- pmax(abs(x), bound)
            v <- pmax(abs(y), bound)
            w <- pmax(abs(y), bound)
            p <- pmax(abs(r - w), bound)
            q <- pmax(abs(y), bound)
            mu <- as.vector(crossprod(z,g) + crossprod(v,w) + crossprod(s,t) + crossprod(p,q))/(2 * (m + n))
            sigfig <- 0
            counter <- 0
            alfa <- 1
            if (verb > 0)	                       # print at least one status report
              cat("Iter    PrimalInf  DualInf  SigFigs  Rescale  PrimalObj  DualObj","\n")
            
            while (counter < maxiter)
            {
              ## update the iteration counter
              counter <- counter + 1
              ## central path (predictor)
              if(smw == 0)
                H.dot.x <- H %*% x
              else if (smw == 1)
                H.dot.x <- H %*% crossprod(H,x)
              rho <- b - A %*% x + w
              nu <- l - x + g
              tau <- u - x - t
              alpha <- r - w - p
              sigma <- c  - crossprod(A, y) - z + s + H.dot.x
              beta <- y + q - v
              gamma.z <- - z
              gamma.w <- - w
              gamma.s <- - s
              gamma.q <- - q
              ## instrumentation
              x.dot.H.dot.x <-  crossprod(x, H.dot.x)
              primal.infeasibility <- max(svd(rbind(rho, tau, matrix(alpha), nu))$d)/ b.plus.1
              dual.infeasibility <- max(svd(rbind(sigma,t(t(beta))))$d) / c.plus.1
              primal.obj <- crossprod(c,x) + 0.5 * x.dot.H.dot.x
              dual.obj <- crossprod(b,y) - 0.5 * x.dot.H.dot.x + crossprod(l, z) - crossprod(u,s) - crossprod(r,q)
              old.sigfig <- sigfig
              sigfig <- max(-log10(abs(primal.obj - dual.obj)/(abs(primal.obj) + 1)), 0)
              if (sigfig >= sigf) break
              if (verb > 0)		      	# final report
                cat( counter, "\t", signif(primal.infeasibility,6), signif(dual.infeasibility,6), sigfig, alfa, primal.obj, dual.obj,"\n")
              ## some more intermediate variables (the hat section)
              hat.beta <- beta - v * gamma.w / w
              hat.alpha <- alpha - p * gamma.q / q
              hat.nu <- nu + g * gamma.z / z
              hat.tau <- tau - t * gamma.s / s
              ## the diagonal terms
              d <- z / g + s / t
              e <- 1 / (v / w + q / p)
              ## initialization before the big cholesky
              if (smw == 0)
                diag(H.x) <- H.diag + d
              diag(H.y) <- e
              c.x <- sigma - z * hat.nu / g - s * hat.tau / t
              c.y <- rho - e * (hat.beta - q * hat.alpha / p)
              ## and solve the system [-H.x A' A H.y] [delta.x, delta.y] <- [c.x c.y]
              if(smw == 0){
                AP[xp,xp] <- -H.x
                AP[xp == FALSE, xp== FALSE] <- H.y
                s1.tmp <- solve(AP,c(c.x,c.y))
                delta.x<-s1.tmp[1:n] ; delta.y <- s1.tmp[-(1:n)]
              }
              else
              {
                V <- diag(smwn)
                smwinner <- chol(V + chunkmult(t(H),2000,d))
                smwa1 <- t(A)
                smwa1 <- smwa1 / d
                smwc1 <- c.x / d
                smwa2 <- t(A) - (H %*% solve(smwinner,solve(t(smwinner),crossprod(H,smwa1))))
                smwa2 <- smwa2 / d 
                smwc2 <- (c.x - (H %*% solve(smwinner,solve(t(smwinner),crossprod(H,smwc1)))))/d 
                delta.y <- solve(A %*% smwa2 + H.y , c.y + A %*% smwc2)
                delta.x <- smwa2 %*% delta.y - smwc2
              }
              
              ## backsubstitution
              delta.w <- - e * (hat.beta - q * hat.alpha / p + delta.y)
              delta.s <- s * (delta.x - hat.tau) / t
              delta.z <- z * (hat.nu - delta.x) / g
              delta.q <- q * (delta.w - hat.alpha) / p
              delta.v <- v * (gamma.w - delta.w) / w
              delta.p <- p * (gamma.q - delta.q) / q
              delta.g <- g * (gamma.z - delta.z) / z
              delta.t <- t * (gamma.s - delta.s) / s
              ## compute update step now (sebastian's trick)
              alfa <- - (1 - margin) / min(c(delta.g / g, delta.w / w, delta.t / t, delta.p / p, delta.z / z, delta.v / v, delta.s / s, delta.q / q, -1))
              newmu <- (crossprod(z,g) + crossprod(v,w) + crossprod(s,t) + crossprod(p,q))/(2 * (m + n))
              newmu <- mu * ((alfa - 1) / (alfa + 10))^2
              gamma.z <- mu / g - z - delta.z * delta.g / g
              gamma.w <- mu / v - w - delta.w * delta.v / v
              gamma.s <- mu / t - s - delta.s * delta.t / t
              gamma.q <- mu / p - q - delta.q * delta.p / p
              ## some more intermediate variables (the hat section)
              hat.beta <- beta - v * gamma.w / w
              hat.alpha <- alpha - p * gamma.q / q
              hat.nu <- nu + g * gamma.z / z
              hat.tau <- tau - t * gamma.s / s
              ## initialization before the big cholesky
              ##for (  i  in  1 : n H.x(i,i) <- H.diag(i) + d(i) ) {
              ##H.y <- diag(e)
              c.x <- sigma - z * hat.nu / g - s * hat.tau / t
              c.y <- rho - e * (hat.beta - q * hat.alpha / p)
              
              ## and solve the system [-H.x A' A H.y] [delta.x, delta.y] <- [c.x c.y]
              if (smw == 0)
              {
                AP[xp,xp] <- -H.x
                AP[xp == FALSE, xp== FALSE] <- H.y
                s1.tmp <- solve(AP,c(c.x,c.y))
                delta.x<-s1.tmp[1:n] ; delta.y<-s1.tmp[-(1:n)]
              }
              else if (smw == 1)
              {
                smwc1 <- c.x / d
                smwc2 <- (c.x - (H %*% solve(smwinner,solve(t(smwinner),crossprod(H,smwc1))))) / d
                delta.y <- solve(A %*% smwa2 + H.y , c.y + A %*% smwc2)
                delta.x <- smwa2 %*% delta.y - smwc2
              }
              ## backsubstitution
              delta.w <- - e * (hat.beta - q * hat.alpha / p + delta.y)
              delta.s <- s * (delta.x - hat.tau) / t
              delta.z <- z * (hat.nu - delta.x) / g
              delta.q <- q * (delta.w - hat.alpha) / p
              delta.v <- v * (gamma.w - delta.w) / w
              delta.p <- p * (gamma.q - delta.q) / q
              delta.g <- g * (gamma.z - delta.z) / z
              delta.t <- t * (gamma.s - delta.s) / s
              ## compute the updates
              alfa <- - (1 - margin) / min(c(delta.g / g, delta.w / w, delta.t / t, delta.p / p, delta.z / z, delta.v / v, delta.s / s, delta.q / q, -1))
              x <- x + delta.x * alfa
              g <- g + delta.g * alfa
              w <- w + delta.w * alfa
              t <- t + delta.t * alfa
              p <- p + delta.p * alfa
              y <- y + delta.y * alfa
              z <- z + delta.z * alfa
              v <- v + delta.v * alfa
              s <- s + delta.s * alfa
              q <- q + delta.q * alfa
              ## these two lines put back in ?
              ## mu <- (crossprod(z,g) + crossprod(v,w) + crossprod(s,t) + crossprod(p,q))/(2 * (m + n))
              ## mu <- mu * ((alfa - 1) / (alfa + 10))^2
              mu <- newmu
            }
            if (verb > 0)		      	## final report
              cat( counter, primal.infeasibility, dual.infeasibility, sigfig, alfa, primal.obj, dual.obj)
            
            ret <- new("ipop")               ## repackage the results
            primal(ret) <- x
            dual(ret)   <- drop(y)
            if ((sigfig > sigf) & (counter < maxiter))
              how(ret)    <- 'converged'
            else
            {					## must have run out of counts
              if ((primal.infeasibility > 10e5) & (dual.infeasibility > 10e5))
                how(ret)    <- 'primal and dual infeasible'
              if (primal.infeasibility > 10e5)
                how(ret)    <- 'primal infeasible'
              if (dual.infeasibility > 10e5)
                how(ret)    <- 'dual infeasible'
              else					## don't really know
                how(ret)    <- 'slow convergence, change bound?'
            }
            ret
          })


setGeneric("chunkmult",function(Z, csize, colscale) standardGeneric("chunkmult"))
setMethod("chunkmult",signature(Z="matrix"),
          function(Z, csize, colscale)
          {
            n <- dim(Z)[1]
            m <- dim(Z)[2]
            d <- sqrt(colscale)
            nchunks <- ceiling(m/csize)
            res <- matrix(0,n,n)
            
            for( i in 1:nchunks)
            {
              lowerb <- (i - 1) * csize + 1 
              upperb <- min(i * csize, m)
              buffer <- t(Z[,lowerb:upperb,drop = FALSE])
              bufferd <- d[lowerb:upperb]
              buffer <- buffer / bufferd
              res <- res + crossprod(buffer)
            }
            return(res)
          })

