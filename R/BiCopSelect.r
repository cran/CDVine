BiCopSelect <- function(u1,u2,familyset=NA,selectioncrit="AIC",indeptest=FALSE,level=0.05) 
{
  if(is.null(u1)==TRUE || is.null(u2)==TRUE) stop("u1 and/or u2 are not set or have length zero.")
  if(length(u1)!=length(u2)) stop("Lengths of 'u1' and 'u2' do not match.")
  if(length(u1)<2) stop("Number of observations has to be at least 2.")	
  if(!is.na(familyset[1])) for(i in 1:length(familyset)) if(!(familyset[i] %in% c(0,1,2,3,4,5,6,7,9,13,14,16,23,24,26,33,34,36))) stop("Copula family not implemented.")  
  if(selectioncrit != "AIC" && selectioncrit != "BIC") stop("Selection criterion not implemented.")
  if(level < 0 & level > 1) stop("Significance level has to be between 0 and 1.")
  
  out=list()

  data1 = u1
  data2 = u2

  if(!is.na(familyset[1]) & familyset[1] == 0){

    out$p.value.indeptest = NA
    out$family = 0
    out$par = c(0,0)

  }else{
  
    if(!is.na(familyset[1]) && (!any(c(1,2,5,23,24,26,33,34,36) %in% familyset) || !any(c(1:7,9,13,14,16) %in% familyset))) stop("'familyset' has to include at least one bivariate copula family for positive and one for negative dependence.")

    #emp_tau = cor(data1,data2,method="kendall")
    emp_tau = fasttau(data1,data2)

    if(indeptest == TRUE){
      out$p.value.indeptest = BiCopIndTest(data1,data2)$p.value
    }else{
      out$p.value.indeptest = NA
    }

    if(!is.na(out$p.value.indeptest) & out$p.value.indeptest >= level){

      out$family = 0
      out$par = c(0,0)      

    }else{

      start = list()
	    start[[1]] = sin(pi * emp_tau/2)
	    start[[2]] = start[[1]]
	    start[[3]] = 2*abs(emp_tau)/(1-abs(emp_tau))
	    start[[4]] = 1/(1-abs(emp_tau))
	    start[[5]] = Frank.itau.JJ(emp_tau)
	    start[[6]] = Joe.itau.JJ(abs(emp_tau))
	    start[[7]] = c(0.5, 1.5)
	    start[[8]] = c(0.5, 0.5)
	    start[[9]] = c(1.5, 0.5)
	    start[[13]] = 2*abs(emp_tau)/(1-abs(emp_tau))
	    start[[14]] = 1/(1-abs(emp_tau))
	    start[[16]] = Joe.itau.JJ(abs(emp_tau))
	    start[[23]] = -2*abs(emp_tau)/(1-abs(emp_tau))
	    start[[24]] = -1/(1-abs(emp_tau))
	    start[[26]] = -Joe.itau.JJ(abs(emp_tau))
	    start[[33]] = -2*abs(emp_tau)/(1-abs(emp_tau))
	    start[[34]] = -1/(1-abs(emp_tau))
	    start[[36]] = -Joe.itau.JJ(abs(emp_tau))

	    lb = list()
	    lb[[1]] = -0.9999
	    lb[[2]] = -0.9999
	    lb[[3]] = 0.0001
	    lb[[4]] = 1.0001
	    lb[[5]] = -Inf
	    lb[[6]] = 1
	    lb[[7]] = c(0.001,1.001)
	    lb[[8]] = c(0.001,0.001)
	    lb[[9]] = c(1.001,0.001)
	    lb[[13]] = 0.0001
	    lb[[14]] = 1.0001
	    lb[[16]] = 1.0001
	    lb[[23]] = -Inf
	    lb[[24]] = -Inf
	    lb[[26]] = -Inf
	    lb[[33]] = -Inf
	    lb[[34]] = -Inf
	    lb[[36]] = -Inf

	    ub = list()
	    ub[[1]] = 0.9999
	    ub[[2]] = 0.9999
	    ub[[3]] = Inf
	    ub[[4]] = Inf
	    ub[[5]] = Inf
	    ub[[6]] = Inf
	    ub[[7]] = c(5,6)
	    ub[[8]] = c(5,5)
	    ub[[9]] = c(5,6)
	    ub[[13]] = Inf
	    ub[[14]] = Inf
	    ub[[16]] = Inf
	    ub[[23]] = -0.0001
	    ub[[24]] = -1.0001
	    ub[[26]] = -1.0001
	    ub[[33]] = -0.0001
	    ub[[34]] = -1.0001
	    ub[[36]] = -1.0001

    	if(emp_tau < 0)
    	{
    		todo = c(1,2,5,23,24,26,33,34,36)
    	}
    	else
    	{
    		todo = c(1:7,9,13,14,16)
    	}
    	if(!is.na(familyset[1])) todo = todo[which(todo %in% familyset)]

      optiout = list()

      if(any(todo == 2)){
    		t_LL = function(param){
    				ll = .C("LL",as.integer(2),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
    		}
    		optiout[[2]] = optim(par=c(start[[2]], 5),fn=t_LL,
    				method="L-BFGS-B",
    				lower=c(-0.9999, 1.0001),
    				upper=c(0.9999, 30),
    				control=list(fnscale=-1,maxit = 500))

    		if(optiout[[2]]$par[2] >= 30){
    		  todo[todo==2] = 1
          todo = unique(todo)
          optiout[[2]] = list()
    		}
    	}

    	if(any(todo == 7)){
    	  t_LL = function(param){
    				ll = .C("LL",as.integer(7),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[7]] = optim(par=start[[7]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[7]],
    				upper=ub[[7]],
    				control=list(fnscale=-1,maxit = 500))
        if(optiout[[7]]$par[1] <= 0.1 | optiout[[7]]$par[2] <= 1.1){
    		  if(optiout[[7]]$par[1] <= 0.1){
            todo[todo==7] = 4
            todo = unique(todo)
    		  }else if(optiout[[7]]$par[2] <= 1.1){
      		  todo[todo==7] = 3
            todo = unique(todo)
      		}
      		optiout[[7]] = list()
        }
      }

    	if(any(todo == 8)){
    	  t_LL = function(param){
    				ll = .C("LL",as.integer(8),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[8]] = optim(par=start[[8]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[8]],
    				upper=ub[[8]],
    				control=list(fnscale=-1,maxit = 500))
        if(optiout[[8]]$par[2] <= 0.1){
    		  todo[todo==8] = 3
          todo = unique(todo)
          optiout[[8]] = list()
        }
      }

    	if(any(todo == 9)){
    	  t_LL = function(param){
    				ll = .C("LL",as.integer(9),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[9]] = optim(par=start[[9]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[9]],
    				upper=ub[[9]],
    				control=list(fnscale=-1,maxit = 500))
        if(optiout[[9]]$par[1] <= 1.1 | optiout[[9]]$par[2] <= 0.1){
          if(optiout[[9]]$par[1] <= 1.1){
      		  todo[todo==9] = 3
            todo = unique(todo)
      		}else if(optiout[[9]]$par[2] <= 0.1){
      		  todo[todo==9] = 6
            todo = unique(todo)
      		}
          optiout[[9]] = list()
        }
      }

      for(i in todo[todo!=2 & todo!=7 & todo!=8 & todo!=9]){

    			t_LL = function(param){
    				ll = .C("LL_mod2",as.integer(i),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param),as.double(0),as.double(0),PACKAGE='CDVine')[[7]]
   			    if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)

        }
        optiout[[i]] = optim(par=start[[i]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[i]],
    				upper=ub[[i]],
    				control=list(fnscale=-1,maxit = 500))
    	}

      if(selectioncrit == "AIC"){

        AICs = rep(Inf,36)

  		  for(i in todo){
  		    if(i %in% c(2,7,8,9)){
  		      AICs[i] = -2*optiout[[i]]$value + 4
          }else{
            AICs[i] = -2*optiout[[i]]$value + 2
          }
        }

        #out$AICs = AICs
  		  out$family = todo[which.min(AICs[todo])]

 		  }else if(selectioncrit == "BIC"){

        BICs = rep(Inf,36)

  		  for(i in todo){
  		    if(i %in% c(2,7,8,9)){
  		      BICs[i] = -2*optiout[[i]]$value + 2*length(data1)
          }else{
            BICs[i] = -2*optiout[[i]]$value + length(data1)
          }
        }

        #out$BICs = BICs
  		  out$family = todo[which.min(BICs[todo])]

 		  }

      #out$optimout = optiout
  		out$par = optiout[[out$family]]$par
   	  if(out$family != 2 & out$family != 7 & out$family != 8 & out$family != 9) out$par[2] = 0

    	#out$emp_tau = emp_tau

      #out$CondOn.1 = .C("Hfunc2",as.integer(out$family),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(out$par[1]),as.double(out$par[2]),as.double(rep(0,length(data1))),PACKAGE='CDVine')[[7]]
      #out$CondOn.2 = .C("Hfunc1",as.integer(out$family),as.integer(length(data1)),as.double(data2),as.double(data1),as.double(out$par[1]),as.double(out$par[2]),as.double(rep(0,length(data1))),PACKAGE='CDVine')[[7]]

    }
  }
  #neu
  out$par2=out$par[2]
  out$par=out$par[1]

  return(out)

}