BiCopSelect <- function(u1,u2,familyset=NA,selectioncrit="AIC",indeptest=FALSE,level=0.05) 
{
  if(is.null(u1)==TRUE || is.null(u2)==TRUE) stop("u1 and/or u2 are not set or have length zero.")
  if(length(u1)!=length(u2)) stop("Lengths of 'u1' and 'u2' do not match.")
  if(length(u1)<2) stop("Number of observations has to be at least 2.")	
  if(!is.na(familyset[1])) for(i in 1:length(familyset)) if(!(familyset[i] %in% c(0,1,2,3,4,5,6,7,8,9,10,13,14,16,17,18,19,20,23,24,26,27,28,29,30,33,34,36,37,38,39,40))) stop("Copula family not implemented.")  
  if(selectioncrit != "AIC" && selectioncrit != "BIC") stop("Selection criterion not implemented.")
  if(level < 0 & level > 1) stop("Significance level has to be between 0 and 1.")
  
  out=list()

  data1 = u1
  data2 = u2

  if(!is.na(familyset[1]) & any(familyset == 0)){

    out$p.value.indeptest = NA
    out$family = 0
    out$par = c(0,0)

  }else{
  
    if(!is.na(familyset[1]) && (!any(c(1,2,5,23,24,26:30,33,34,36:40) %in% familyset) || !any(c(1:10,13,14,16:20) %in% familyset))) stop("'familyset' has to include at least one bivariate copula family for positive and one for negative dependence.")

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
	    start[[8]] = c(1.5, 1.5)
	    start[[9]] = c(1.5, 0.5)
	    start[[10]] = c(1.5,0.5)
	    start[[13]] = 2*abs(emp_tau)/(1-abs(emp_tau))
	    start[[14]] = 1/(1-abs(emp_tau))
	    start[[16]] = Joe.itau.JJ(abs(emp_tau))
	    start[[17]] = c(0.5, 1.5)
	    start[[18]] = c(1.5, 1.5)
	    start[[19]] = c(1.5, 0.5)
	    start[[20]] = c(1.5,0.5)
	    start[[23]] = -2*abs(emp_tau)/(1-abs(emp_tau))
	    start[[24]] = -1/(1-abs(emp_tau))
	    start[[26]] = -Joe.itau.JJ(abs(emp_tau))
	    start[[27]] = c(-0.5, -1.5)
	    start[[28]] = c(-1.5, -1.5)
	    start[[29]] = c(-1.5, -0.5)
	    start[[30]] = c(-1.5,-0.5)
	    start[[33]] = -2*abs(emp_tau)/(1-abs(emp_tau))
	    start[[34]] = -1/(1-abs(emp_tau))
	    start[[36]] = -Joe.itau.JJ(abs(emp_tau))
	    start[[37]] = c(-0.5, -1.5)
	    start[[38]] = c(-1.5, -1.5)
	    start[[39]] = c(-1.5, -0.5)
	    start[[40]] = c(-1.5,-0.5)

	    lb = list()
	    lb[[1]] = -0.9999
	    lb[[2]] = -0.9999
	    lb[[3]] = 0.0001
	    lb[[4]] = 1.0001
	    lb[[5]] = -Inf
	    lb[[6]] = 1
	    lb[[7]] = c(0.001,1.001)
	    lb[[8]] = c(1.001,1.001)
	    lb[[9]] = c(1.001,0.001)
	    lb[[10]] = c(1.001,0.001)
	    lb[[13]] = 0.0001
	    lb[[14]] = 1.0001
	    lb[[16]] = 1.0001
	    lb[[17]] = c(0.001,1.001)
	    lb[[18]] = c(1.001,1.001)
	    lb[[19]] = c(1.001,0.001)
	    lb[[20]] = c(1.001,0.001)
	    lb[[23]] = -Inf
	    lb[[24]] = -Inf
	    lb[[26]] = -Inf
	    lb[[27]] = c(-5,-6)
	    lb[[28]] = c(-6,-6)
	    lb[[29]] = c(-5,-6)
	    lb[[30]] = c(-6,-0.999)
	    lb[[33]] = -Inf
	    lb[[34]] = -Inf
	    lb[[36]] = -Inf
	    lb[[37]] = c(-5,-6)
	    lb[[38]] = c(-6,-6)
	    lb[[39]] = c(-5,-6)
	    lb[[40]] = c(-6,-0.999)
	    
	    ub = list()
	    ub[[1]] = 0.9999
	    ub[[2]] = 0.9999
	    ub[[3]] = Inf
	    ub[[4]] = Inf
	    ub[[5]] = Inf
	    ub[[6]] = Inf
	    ub[[7]] = c(5,6)
	    ub[[8]] = c(6,6)
	    ub[[9]] = c(5,6)
	    ub[[10]] = c(6,0.999)
	    ub[[13]] = Inf
	    ub[[14]] = Inf
	    ub[[16]] = Inf
	    ub[[17]] = c(5,6)
	    ub[[18]] = c(6,6)
	    ub[[19]] = c(5,6)
	    ub[[20]] = c(6,0.999)
	    ub[[23]] = -0.0001
	    ub[[24]] = -1.0001
	    ub[[26]] = -1.0001
	    ub[[27]] = c(-0.001,-1.001)
	    ub[[28]] = c(-1.001,-1.001)
	    ub[[29]] = c(-1.001,-0.001)
	    ub[[30]] = c(-1.001,-0.001)
	    ub[[33]] = -0.0001
	    ub[[34]] = -1.0001
	    ub[[36]] = -1.0001
	    ub[[37]] = c(-0.001,-1.001)
	    ub[[38]] = c(-1.001,-1.001)
	    ub[[39]] = c(-1.001,-0.001)
	    ub[[40]] = c(-1.001,-0.001)
	    

    	if(emp_tau < 0)
    	{
    		todo = c(1,2,5,23,24,26:30,33,34,36:40)
    	}
    	else
    	{
    		todo = c(1:10,13,14,16:20)
    	}
    	if(!is.na(familyset[1])) todo = todo[which(todo %in% familyset)]

      optiout = list()

     if(any(todo == 2)){
			t_LL_start = function(param){
                     ll = .C("LL",as.integer(2),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(start[[2]]),as.double(param),as.double(0),PACKAGE='CDVine')[[7]]
                     if(is.infinite(ll) || is.na(ll))
                         ll = -10^300
                     return(ll)
             }
         optim_start = optim(par=5,fn=t_LL_start,method="L-BFGS-B",lower=1.0001,upper=30,control=list(fnscale=-1,maxit= 500))

             t_LL = function(param){
                     ll = .C("LL",as.integer(2),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
                     if(is.infinite(ll) || is.na(ll))
                         ll = -10^300
                     return(ll)
             }
             optiout[[2]] = optim(par=c(start[[2]], 
		     optim_start$par),fn=t_LL,
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
        if(optiout[[8]]$par[1] <= 1.1 | optiout[[8]]$par[2] <= 1.1){
    		  if(optiout[[8]]$par[1] <= 1.1){
            todo[todo==8] = 4
            todo = unique(todo)
    		  }else if(optiout[[8]]$par[2] <= 1.1){
      		  todo[todo==8] = 6
            todo = unique(todo)
      		}
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

      if(any(todo == 10)){
    	  t_LL = function(param){
    				ll = .C("LL",as.integer(10),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[10]] = optim(par=start[[10]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[10]],
    				upper=ub[[10]],
    				control=list(fnscale=-1,maxit = 500))
		    if(optiout[[10]]$par[2] >= 0.99){
      		  todo[todo==10] = 6
            todo = unique(todo)
          optiout[[10]] = list()
        }
      }

      if(any(todo == 17)){
    	  t_LL = function(param){
    				ll = .C("LL",as.integer(17),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[17]] = optim(par=start[[17]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[17]],
    				upper=ub[[17]],
    				control=list(fnscale=-1,maxit = 500))
        if(optiout[[17]]$par[1] <= 0.1 | optiout[[17]]$par[2] <= 1.1){
    		  if(optiout[[17]]$par[1] <= 0.1){
            todo[todo==17] = 14
            todo = unique(todo)
    		  }else if(optiout[[17]]$par[2] <= 1.1){
      		  todo[todo==17] = 13
            todo = unique(todo)
      		}
      		optiout[[17]] = list()
        }
      }

      if(any(todo == 18)){
    	  t_LL = function(param){
    				ll = .C("LL",as.integer(18),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[18]] = optim(par=start[[18]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[18]],
    				upper=ub[[18]],
    				control=list(fnscale=-1,maxit = 500))
        if(optiout[[18]]$par[1] <= 1.1 | optiout[[18]]$par[2] <= 1.1){
    		  if(optiout[[18]]$par[1] <= 1.1){
            todo[todo==18] = 14
            todo = unique(todo)
    		  }else if(optiout[[18]]$par[2] <= 1.1){
      		  todo[todo==18] = 16
            todo = unique(todo)
      		}
      		optiout[[18]] = list()
        }        
      }

      if(any(todo == 19)){
    	  t_LL = function(param){
    				ll = .C("LL",as.integer(19),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[19]] = optim(par=start[[19]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[19]],
    				upper=ub[[19]],
    				control=list(fnscale=-1,maxit = 500))
        if(optiout[[19]]$par[1] <= 1.1 | optiout[[19]]$par[2] <= 0.1){
          if(optiout[[19]]$par[1] <= 1.1){
      		  todo[todo==19] = 13
            todo = unique(todo)
      		}else if(optiout[[19]]$par[2] <= 0.1){
      		  todo[todo==19] = 16
            todo = unique(todo)
      		}
          optiout[[19]] = list()
        }
      }

      if(any(todo == 20)){
    	  t_LL = function(param){
    				ll = .C("LL",as.integer(20),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[20]] = optim(par=start[[20]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[20]],
    				upper=ub[[20]],
    				control=list(fnscale=-1,maxit = 500))
		    if(optiout[[20]]$par[2] >= 0.99){
      		  todo[todo==20] = 16
            todo = unique(todo)
          optiout[[20]] = list()
        }
      }

      if(any(todo == 27)){
    	  t_LL = function(param){
    				ll = .C("LL_mod2",as.integer(27),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[27]] = optim(par=start[[27]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[27]],
    				upper=ub[[27]],
    				control=list(fnscale=-1,maxit = 500))
        if(optiout[[27]]$par[1] >= -0.1 | optiout[[27]]$par[2] >= -1.1){
    		  if(optiout[[27]]$par[1] >= -0.1){
            todo[todo==27] = 24
            todo = unique(todo)
    		  }else if(optiout[[27]]$par[2] >= -1.1){
      		  todo[todo==27] = 23
            todo = unique(todo)
      		}
      		optiout[[27]] = list()
        }
      }

      if(any(todo == 28)){
    	  t_LL = function(param){
    				ll = .C("LL_mod2",as.integer(28),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[28]] = optim(par=start[[28]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[28]],
    				upper=ub[[28]],
    				control=list(fnscale=-1,maxit = 500))
        if(optiout[[28]]$par[1] >= -1.1 | optiout[[28]]$par[2] >= -1.1){
    		  if(optiout[[28]]$par[1] >= -1.1){
            todo[todo==28] = 24
            todo = unique(todo)
    		  }else if(optiout[[28]]$par[2] >= -1.1){
      		  todo[todo==28] = 26
            todo = unique(todo)
      		}
      		optiout[[28]] = list()
        }         
      }

      if(any(todo == 29)){
    	  t_LL = function(param){
    				ll = .C("LL_mod2",as.integer(29),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[29]] = optim(par=start[[29]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[29]],
    				upper=ub[[29]],
    				control=list(fnscale=-1,maxit = 500))
        if(optiout[[29]]$par[1] >= -1.1 | optiout[[29]]$par[2] >= -0.1){
          if(optiout[[29]]$par[1] >= -1.1){
      		  todo[todo==29] = 23
            todo = unique(todo)
      		}else if(optiout[[29]]$par[2] >= -0.1){
      		  todo[todo==29] = 26
            todo = unique(todo)
      		}
          optiout[[29]] = list()
        }
      }

      if(any(todo == 30)){
    	  t_LL = function(param){
    				ll = .C("LL_mod2",as.integer(30),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[30]] = optim(par=start[[30]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[30]],
    				upper=ub[[30]],
    				control=list(fnscale=-1,maxit = 500))
		    if(optiout[[30]]$par[2] <= -0.99){
      		  todo[todo==30] = 26
            todo = unique(todo)
          optiout[[30]] = list()
        }
      }

      if(any(todo == 37)){
    	  t_LL = function(param){
    				ll = .C("LL_mod2",as.integer(37),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[37]] = optim(par=start[[37]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[37]],
    				upper=ub[[37]],
    				control=list(fnscale=-1,maxit = 500))
        if(optiout[[37]]$par[1] >= -0.1 | optiout[[37]]$par[2] >= -1.1){
    		  if(optiout[[37]]$par[1] >= -0.1){
            todo[todo==37] = 24
            todo = unique(todo)
    		  }else if(optiout[[37]]$par[2] >= -1.1){
      		  todo[todo==37] = 33
            todo = unique(todo)
      		}
      		optiout[[37]] = list()
        }
      }

      if(any(todo == 38)){
    	  t_LL = function(param){
    				ll = .C("LL_mod2",as.integer(38),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[38]] = optim(par=start[[38]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[38]],
    				upper=ub[[38]],
    				control=list(fnscale=-1,maxit = 500))
        if(optiout[[38]]$par[1] >= -1.1 | optiout[[38]]$par[2] >= -1.1){
    		  if(optiout[[38]]$par[1] >= -1.1){
            todo[todo==38] = 24
            todo = unique(todo)
    		  }else if(optiout[[38]]$par[2] >= -1.1){
      		  todo[todo==38] = 26
            todo = unique(todo)
      		}
      		optiout[[38]] = list()
        } 
      }

      if(any(todo == 39)){
    	  t_LL = function(param){
    				ll = .C("LL_mod2",as.integer(39),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[39]] = optim(par=start[[39]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[39]],
    				upper=ub[[39]],
    				control=list(fnscale=-1,maxit = 500))
        if(optiout[[39]]$par[1] >= -1.1 | optiout[[39]]$par[2] >= -0.1){
          if(optiout[[39]]$par[1] >= -1.1){
      		  todo[todo==39] = 33
            todo = unique(todo)
      		}else if(optiout[[39]]$par[2] >= -0.1){
      		  todo[todo==39] = 36
            todo = unique(todo)
      		}
          optiout[[39]] = list()
        }
      }

      if(any(todo == 40)){
    	  t_LL = function(param){
    				ll = .C("LL_mod2",as.integer(40),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(param[1]),as.double(param[2]),as.double(0),PACKAGE='CDVine')[[7]]
    				if(is.infinite(ll) || is.na(ll))
    					ll = -10^300
    				return(ll)
	      }
	      optiout[[40]] = optim(par=start[[40]],fn=t_LL,
    				method="L-BFGS-B",
    				lower=lb[[40]],
    				upper=ub[[40]],
    				control=list(fnscale=-1,maxit = 500))
		    if(optiout[[40]]$par[2] <= -0.99){
      		  todo[todo==40] = 36
            todo = unique(todo)
          optiout[[40]] = list()
        }
      }

      for(i in todo[!(todo%in%c(2,7:10,17:20,27:30,37:40))]){

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

        AICs = rep(Inf,40)

  		  for(i in todo){
  		    if(i %in% c(2,7:10,17:20,27:30,37:40)){
  		      AICs[i] = -2*optiout[[i]]$value + 4
          }else{
            AICs[i] = -2*optiout[[i]]$value + 2
          }
        }

        #out$AICs = AICs
  		  out$family = todo[which.min(AICs[todo])]

 		  }else if(selectioncrit == "BIC"){

        BICs = rep(Inf,40)

  		  for(i in todo){
  		    if(i %in% c(2,7:10,17:20,27:30,37:40)){
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
   	  if(!(out$family%in%c(2,7:10,17:20,27:30,37:40)) ) out$par[2] = 0

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