data one;
input id drink time trt;
drinkyes=drink>0;
cards;
1 0 1 0
1 2.5 2 0
1 3 3 0
2 0 1 1
2 5 2 1
2 4.5 3 1
...
;
run;

* Part I model, with random intercept only;
proc nlmixed data=one;
parms alpha0=3 alpha1=0 var1=2;
bounds var1>=0;
teta=alpha0 + a + alpha1 * time +alpha2 * trt;
expteta=exp(teta);
p=expteta / (1+expteta);
model drinkyes ~ binomial(1, p);
random a ~ normal(0, var1) subject=id;
run;
* Part I model, with random intercept and slope;
proc nlmixed data=one;
parms alpha0=3 alpha1=0 alpha2=1 var1=2 var2=1 cov12=0;
bounds var1 var2>=0;
teta=alpha0 + a + b* time + alpha1 * time +alpha2 * trt;
expteta=exp(teta);
p=expteta / (1+expteta);
model drinkyes ~ binomial(1, p);
random a b~ normal([0, 0], [var1, cov12, var2]) subject=id;
run;

* log skew normal model with heteroscedasticity, random intercept and slope;
proc nlmixed data=liulei.all qpoints=5;
parms  alpha0=3 alpha1=0 alpha2=1 
		beta0=1.5 beta1=.2 beta2=0
		gamma0=-1 gamma1=0 gamma2=0
		var1=1 var2=1 var3=1 var4=1 cov12=0 cov13=0 cov14=0 cov23=0 cov24=0 cov34=0 
		lambda=0;

bounds var1 var2 var3 var4>=0;


teta=alpha0 + a + b* time + alpha1 * time + alpha2 * trt  ;
expteta=exp(teta);
p=expteta / (1+expteta);


if drink=0 then loglik=log(1-p);
if drink>0 then do;

	mu=beta0 + c + d* time + beta1 * time + beta2 * trt  ;

	theta=gamma0 + gamma1 * time + gamma2 * trt;
	sigma=exp(theta/2);


	value1=sigma**2 + lambda **2;
	value2=(log(drink)-mu) / sqrt(value1);
	value3=value2 * lambda/sigma;

	logden=log(2/drink) -.5 * log(value1) - .5 *log(2*3.1415926) -.5 *value2**2 + log(CDF('NORMAL',value3)) ;

	loglik=log(p) + logden;
end;
model drink ~ general(loglik);
random a b c d~ normal([0, 0, 0, 0], [var1, cov12, var2, cov13, cov23, var3, cov14, cov24, cov34, var4]) subject=id;
estimate 'lognormal' lambda;
estimate 'rho12' cov12/sqrt(var1*var2);
estimate 'rho13' cov13/sqrt(var1*var3);
estimate 'rho14' cov14/sqrt(var1*var4);
estimate 'rho23' cov23/sqrt(var2*var3);
estimate 'rho24' cov24/sqrt(var2*var4);
estimate 'rho34' cov34/sqrt(var3*var4);
run;

