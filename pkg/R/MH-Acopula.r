#Acopula

#====general setup====

remove(list=objects())#remove previously declared objects
library(fBasics)#load fBasics for generating Stable random variates
eps=1e-8#for Frank's and Joe's family, F is evaluated until F>=1-eps
maxnumofFvalues=500000#maximal number of computed F values for Frank's and Joe's family
numofsummandsjoe=#number of summands for calculating Kendall's tau for Joe's family

#====general functions====

#wrapper for root finding
rootfindingwrapper=function(theta,fun,value) return(fun(theta)-value)

#compute random variates from an exponentially-tilted Stable distribution \tilde{S}(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),V0*Indicator(alpha==1),h*Indicator(alpha!=1);1)
#with corresponding Laplace-Stieltjes transform exp(-V0*((h+t)^alpha-h^alpha))
retstable=function(alpha,V0,h){
	n=length(V0)
	variates=numeric(n)
	#case alpha==1
	if(alpha==1){
		return(V0)#sample from S(1,1,0,V0;1) with Laplace-Stieltjes transform exp(-V0*t)
	}
	#case alpha!=1 and h==0
	if(h==0){
		for(i in 1:n){
			variates[i]=rstable(1,alpha,1,(cos(alpha*pi/2)*V0[i])^(1/alpha),0,pm=1)#sample from S(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),0;1) with Laplace-Stieltjes transform exp(-V0*t^alpha)
		}	
		return(variates)
	}
	#case alpha!=1 and h!=0
	#call fast rejection algorithm
	for(i in 1:n){
		#determine optimal constant m for the fast rejection algorithm
		logc=V0[i]*h^alpha
		floorlogc=floor(logc)
		ceilinglogc=ceiling(logc)
		c=exp(logc)
		floorvalue=floorlogc*c^(1/floorlogc)
		ceilingvalue=ceilinglogc*c^(1/ceilinglogc)
		if(logc<=1){
			m=1
		}else if(logc>1&&floorvalue<=ceilingvalue){
			m=floorlogc
		}else{
			m=ceilinglogc
		}
		#apply standard rejection algorithm to sample the m-fold sum of variates from the distribution with Laplace-Stieltjes transform exp(-(V0/m)*((h+t)^alpha-h^alpha)) 
		variates[i]=0
		for(j in 1:m){
			repeat{
				Vtilde=rstable(1,alpha,1,(cos(alpha*pi/2)*V0[i]/m)^(1/alpha),0,pm=1)#sample from S(alpha,1,(cos(alpha*pi/2)V0/m)^(1/alpha),0;1) with Laplace-Stieltjes transform exp(-(V0/m)*t^alpha)
				u=runif(1)
				if(u<=exp(-h*Vtilde)){
					variates[i]=variates[i]+Vtilde
					break			
				}
			}
		}
	}
	return(variates)
}

#====Ali-Mikhail-Haq, see Nelsen (2007) p. 116, #3====

#generator
psiA=function(t,theta){
	return((1-theta)/(exp(t)-theta))
}

#generator inverse
psiinvA=function(t,theta){
	return(log((1-theta*(1-t))/t))
}

#parameter constraint
constrA=function(theta){
	return(ifelse(theta>=0&&theta<1,TRUE,FALSE))
}

#nesting constraint
nconstrA=function(theta0,theta1){
	return(ifelse(constrA(theta0)&&constrA(theta1)&&theta1>=theta0,TRUE,FALSE))
}

#V_0
V0A=function(n,theta){
	return(rgeom(n,1-theta)+1)
}

#V_{01}
V01A=function(V0,theta0,theta1){
	variates=sapply(V0,rgeom,prob=(1-theta1)/(1-theta0))
	variates=sapply(variates,sum)+V0
	return(variates)
}

#tau
tauA=function(theta){
	return(1-2*((1-theta)*(1-theta)*log(1-theta)+theta)/(3*theta*theta))
}

#tau inverse
tauinvA=function(tau){
	if(tau>0.33333333333){
		return("Error: It is not possible for an Ali-Mikhail-Haq copula to attain such Kendall's tau")
	}else{
		return(uniroot(f=rootfindingwrapper,interval=c(1e-12,1-1e-12),fun=tauA,value=tau)$root)
	}
}

#lambda_l
lambdalA=function(theta){
	return(0*theta)
}

#lambda_l inverse
lambdalinvA=function(lambda){
	return("Error: Any parameter choice for an Ali-Mikhail-Haq copula leads to a lower tail dependence coefficient equal to zero")
}

#lambda_u
lambdauA=function(theta){
	return(0*theta)
}

#lambda_u inverse
lambdauinvA=function(lambda){
	return("Error: Any parameter choice for an Ali-Mikhail-Haq copula leads to an upper tail dependence coefficient equal to zero")
}

#====Clayton, see Nelsen (2007) p. 116, #1 but we use a slightly simpler form of the generator====

#generator
psiC=function(t,theta){
	return((1+t)^(-1/theta))
}

#generator inverse
psiinvC=function(t,theta){
	return(t^(-theta)-1)
}

#parameter constraint
constrC=function(theta){
	return(ifelse(theta>0,TRUE,FALSE))
}

#nesting constraint
nconstrC=function(theta0,theta1){
	return(ifelse(constrC(theta0)&&constrC(theta1)&&theta1>=theta0,TRUE,FALSE))
}

#V_0
V0C=function(n,theta){
	return(rgamma(n,shape=1/theta))
}

#V_{01}
V01C=function(V0,theta0,theta1){
	alpha=theta0/theta1
	return(retstable(alpha,V0,1))
}

#tau
tauC=function(theta){
	return(theta/(theta+2))
}

#tau inverse
tauinvC=function(tau){
	return(2*tau/(1-tau))
}

#lambda_l
lambdalC=function(theta){
	return(2^(-1/theta))
}

#lambda_l inverse
lambdalinvC=function(lambda){
	return(-log(2)/log(lambda))
}

#lambda_u
lambdauC=function(theta){
	return(0*theta)
}

#lambda_u inverse
lambdauinvC=function(lambda){
	return("Error: Any parameter choice for a Clayton copula leads to an upper tail dependence coefficient equal to zero")
}

#====Frank, see Nelsen (2007) p. 116, #5====

#generator
psiF=function(t,theta){
	return(-log(1-(1-exp(-theta))*exp(-t))/theta)
}

#generator inverse
psiinvF=function(t,theta){
	return(-log((exp(-theta*t)-1)/(exp(-theta)-1)))
}

#parameter constraint
constrF=function(theta){
	return(ifelse(theta>0,TRUE,FALSE))
}

#nesting constraint
nconstrF=function(theta0,theta1){
	return(ifelse(constrF(theta0)&&constrF(theta1)&&theta1>=theta0,TRUE,FALSE))
}

#V_0
V0F=function(n,theta){
	W=runif(1)
	if(W>=1-exp(-theta)){
		return(1)
	}else{
		q=1-exp(theta*runif(1))
		qsquared=q*q
		if(W<=q^2){
			return(floor(1+log(W)/log(q)))
		}else if(qsquared<W&&W<=q){
			return(1)
		}else{
			return(2)
		}
	}
}

#V_{01}
V01F=function(V0,theta0,theta1){
	#compute the function values of the discrete distribution
	Fvalues=numeric(maxnumofFvalues)
	alpha=theta0/theta1
	c0=1-exp(-theta0)
	c1=1-exp(-theta1)
	ysum=(alpha/c0)*c1#y_1
	expsum=log(alpha)
	logc1=log(c1)
	k=1
	while(ysum<1-eps&&k<=maxnumofFvalues){
		Fvalues[k]=ysum
		k=k+1
		#compute y_k
		expsum=expsum+log(k-1-alpha)
		ysum=ysum+exp(expsum+k*logc1-log(factorial(k)))/c0
	}
	m=length(Fvalues)
	#sample by finding quantiles 
	n=length(V0)
	variates=numeric(n)
	#consider warnings
	numoftruncations=0
	for(i in 1:n){
		uniformvariates=runif(V0[i])
		variatesfromF=apply(outer(uniformvariates,Fvalues,">"),1,sum)+1#sample the summands of the sums involved
		numoftruncations=numoftruncations+length(variatesfromF[variatesfromF==m+1])
		variates[i]=sum(variatesfromF)
	}
	if(numoftruncations>0){
		print(paste("Warning: The distribution function F involved in sampling the inner distribution function for nested Frank copulas is truncated",numoftruncations,"times at",m+1,"as the largest computed value of F is",Fvalues[m],sep=" "))
	}
	return(variates)
}

#tau
tauF=function(theta){
	return(1+(4/theta)*((1/theta)*integrate(f=function(t) return(t/(exp(t)-1)),lower=0,upper=theta)[[1]]-1))
}

#tau inverse
tauinvF=function(tau){
	return(uniroot(f=rootfindingwrapper,interval=c(0.001,100),fun=tauA,value=tau)$root)
}

#lambda_l
lambdalF=function(theta){
	return(0*theta)
}

#lambda_l inverse
lambdalinvF=function(lambda){
	return("Error: Any parameter choice for a Frank copula leads to a lower tail dependence coefficient equal to zero")
}

#lambda_u
lambdauF=function(theta){
	return(0*theta)
}

#lambda_u inverse
lambdauinvF=function(lambda){
	return("Error: Any parameter choice for a Frank copula leads to an upper tail dependence coefficient equal to zero")
}

#====Gumbel, see Nelsen (2007) p. 116, #4====

#generator
psiG=function(t,theta){
	return(exp(-t^(1/theta)))
}

#generator inverse
psiinvG=function(t,theta){
	return((-log(t))^theta)
}

#parameter constraint
constrG=function(theta){
	return(ifelse(theta>=1,TRUE,FALSE))
}

#nesting constraint
nconstrG=function(theta0,theta1){
	return(ifelse(constrG(theta0)&&constrG(theta1)&&theta1>=theta0,TRUE,FALSE))
}

#V_0
V0G=function(n,theta){
	if(theta==1){
		return(rep(1,n))#sample from S(1,1,0,1;1) with Laplace-Stieltjes transform exp(-t)
	}else{
		alpha=1/theta
		return(rstable(n,alpha,1,(cos(alpha*pi/2))^(1/alpha),0,pm=1))#sample from S(alpha,1,(cos(alpha*pi/2))^(1/alpha),0;1) with Laplace-Stieltjes transform exp(-t^alpha)
	}
}

#V_{01}
V01G=function(V0,theta0,theta1){
	n=length(V0)
	variates=numeric(n)
	alpha=theta0/theta1
	if(alpha==1){
		return(V0)#sample from S(1,1,0,V0;1) with Laplace-Stieltjes transform exp(-V0*t)
	}else{
		for(i in 1:n){
			variates[i]=rstable(1,alpha,1,(cos(alpha*pi/2)*V0[i])^(1/alpha),0,pm=1)#sample from S(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),0;1) with Laplace-Stieltjes transform exp(-V0*t^alpha)
		}
	}
	return(variates)
}

#tau
tauG=function(theta){
	return((theta-1)/theta)
}

#tau inverse
tauinvG=function(tau){
	return(1/(1-tau))
}

#lambda_l
lambdalG=function(theta){
	return(0*theta)
}

#lambda_l inverse
lambdalinvG=function(lambda){
	return("Error: Any parameter choice for a Gumbel copula leads to a lower tail dependence coefficient equal to zero")
}

#lambda_u
lambdauG=function(theta){
	return(2-2^(1/theta))
}

#lambda_u inverse
lambdauinvG=function(lambda){
	return(log(2)/log(2-lambda))
}

#====Joe, see Nelsen (2007) p. 116, #6====

#generator
psiJ=function(t,theta){
	return(1-(1-exp(-t))^(1/theta))
}

#generator inverse
psiinvJ=function(t,theta){
	return(-log(1-(1-t)^theta))
}

#parameter constraint
constrJ=function(theta){
	return(ifelse(theta>=1,TRUE,FALSE))
}

#nesting constraint
nconstrJ=function(theta0,theta1){
	return(ifelse(constrJ(theta0)&&constrJ(theta1)&&theta1>=theta0,TRUE,FALSE))
}

#auxiliary function for computing the mass functions involved in V0J and V01J
auxJ=function(alpha){
	Fvalues=numeric(maxnumofFvalues)
	ysum=alpha#y_1
	expsum=log(alpha)
	k=1
	while(ysum<1-eps&&k<=maxnumofFvalues){
		Fvalues[k]=ysum
		k=k+1
		#compute y_k
		expsum=expsum+log(k-1-alpha)
		ysum=ysum+exp(expsum-log(factorial(k)))
	}	
	return(Fvalues)
}

#V_0
V0J=function(n,theta){
	#compute the function values of the discrete distribution
	Fvalues=auxJ(1/theta)
	m=length(Fvalues)
	#sample by finding quantiles
	uniformvariates=runif(n)
	variates=apply(outer(uniformvariates,Fvalues,">"),1,sum)+1
	#warning
	numoftruncations=length(variates[variates==m+1])
	if(numoftruncations>0){
		print(paste("Warning: The distribution function corresponding to Joe's family is truncated",numoftruncations,"times at",m+1,"as the largest computed value of F is",Fvalues[m],sep=" "))
	}
	return(variates)
}

#V_{01}
V01J=function(V0,theta0,theta1){
	#compute the function values of the discrete distribution
	alpha=theta0/theta1
	Fvalues=auxJ(alpha)
	m=length(Fvalues)
	#sample by finding quantiles 
	n=length(V0)
	variates=numeric(n)
	#consider warnings
	numoftruncations=0
	for(i in 1:n){
		uniformvariates=runif(V0[i])
		variatesfromF=apply(outer(uniformvariates,Fvalues,">"),1,sum)+1#sample the summands of the sums involved
		numoftruncations=numoftruncations+length(variatesfromF[variatesfromF==m+1])
		variates[i]=sum(variatesfromF)
	}
	if(numoftruncations>0){
		print(paste("Warning: The distribution function F involved in sampling the inner distribution function for nested Joe copulas is truncated",numoftruncations,"times at",m+1,"as the largest computed value of F is",Fvalues[m],sep=" "))
	}
	return(variates)
}

#tau
tauJ=function(theta){
	mysum=0
	for(k in 1:numofsummandsjoe){
		mysum=mysum+1/(k*(theta*k+2)*(theta*(k-1)+2))
	}
	return(1-4*mysum)
}

#tau inverse
tauinvJ=function(tau){
	return(uniroot(f=rootfindingwrapper,interval=c(1.001,100),fun=tauJ,value=tau)$root)
}

#lambda_l
lambdalJ=function(theta){
	return(0*theta)
}

#lambda_l inverse
lambdalinvJ=function(lambda){
	return("Error: Any parameter choice for a Joe copula leads to a lower tail dependence coefficient equal to zero")
}

#lambda_u
lambdauJ=function(theta){
	return(2-2^(1/theta))
}

#lambda_u inverse
lambdauinvJ=function(lambda){
	return(log(2)/log(2-lambda))
}

#====Tests====

psiA(1,0.5)
psiinvA(0.5,0.5)
nconstrA(0.5,0.2)
V0A(10,0.5)
V01A(10,0.5,0.8)
tauA(0.5)
tauinvA(0.2)
lambdalA(0.5)

psiC(1,0.5)
psiinvC(0.5,0.5)
nconstrC(0.2,0.5)
V0=V0C(10,0.5)
V01C(V0,0.5,0.8)
tauC(0.5)
tauinvC(0.2)
lambdalC(0.5)
lambdalinvC(0.5)

n=1000
theta0=tauinvC(0.05)
theta1=tauinvC(0.1)
startclock=proc.time()
V0=V0C(n,theta0)
V01=V01C(V0,theta0,theta1)
runtime=proc.time()-startclock
runtime