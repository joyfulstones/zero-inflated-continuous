libname liulei "s:\preventive_medicine\faculty\liulei\paper\paper36";

data one;
set Mylib.myibd_baseline_copy;
Al=0; Ba=0; Bi=0; Clo=0; col=0; cop=0; di=0; do=0; es=0; eu=0; fa=0; h=0; l=0; p=0; ro=0; ru=0; s=0; v=0;
if Alistipes>0 then Al=1;
if Bacteroides>0 then Ba=1;
if Bifidobacterium>0 then Bi=1;
if Clostridium>0 then clo=1;
if Collinsella>0 then col=1;
if Coprobacillus>0 then cop=1;
if Dialister>0 then di=1;
if Dorea>0 then do=1;
if Escherichia>0 then es=1;
if Eubacterium>0 then eu=1;
if Faecalibacterium>0 then fa=1;
if Haemophilus>0 then h=1;
if Lactobacillus>0 then l=1;
if Parabacteroides>0 then p=1;
if Roseburia>0 then ro=1;
if Ruminococcus>0 then ru=1;
if Streptococcus>0 then s=1;
if Veillonella>0 then v=1;
run;

%macro TPM(var, var2);

PROC nlmixed DATA=one METHOD=GAUSS;
bounds vara>0;
/* log-likelihood */
teta = alpha0 + a  + alpha1*Time + alpha2*Treatment + alpha3*Baseline&var2;
expteta = exp(teta);
p = expteta / (1 + expteta);

model &var ~ binary(p);
random a ~ normal(0, vara) subject=Subject;
ods output ParameterEstimates=est1 FitStatistics=fit1; 

run;

data two;
set one;
if &var>0;
run; 

* Part II for Beta distribution;

PROC nlmixed DATA=two METHOD=GAUSS;
bounds varb phi>0;

/* define initial values and bounds */

/* log-likelihood */
	eta=beta0 + b + beta1*Time + beta2*Treatment + beta3*Baseline&var2;
	mu = exp(eta)/(1+exp(eta));
	loglik = (mu*phi-1)*log(&var2) + ((1-mu)*phi-1)*log(1-&var2) + lgamma(phi) - lgamma(mu*phi) - lgamma((1-mu)*phi);


model &var2 ~ general(loglik);
random b ~ normal(0, varb) subject=Subject;
ods output ParameterEstimates=est2 FitStatistics=fit2; 

run;

data initialvalue1;
set est1 est2;
run;

PROC nlmixed DATA=one METHOD=GAUSS;

/* define initial values and bounds */
parms/data=initialvalue1;

bounds vara varb phi>0;

/* log-likelihood */
teta = alpha0 + a  + alpha1*Time + alpha2*Treatment + alpha3*Baseline&var2;
expteta = exp(teta);
p = expteta / (1 + expteta);


/* For zero Faecalibacterium */
if &var=0 then loglik=log(1-p);

/* For positive Faecalibacterium */
if &var>0 then do;
	mu = exp(beta0 + b + beta1*Time + beta2*Treatment + beta3*Baseline&var2)/(1+exp(beta0 + b + beta1*Time + beta2*Treatment + beta3*Baseline&var2));
	loglik = log(p) + (mu*phi-1)*log(&var2) + ((1-mu)*phi-1)*log(1-&var2) + lgamma(phi) - lgamma(mu*phi) - lgamma((1-mu)*phi);
end;

/* fit the above modle */
model &var ~ general(loglik);
random a b~ normal([0,0],[vara,0,varb]) subject=Subject;
ods output ParameterEstimates=est3 FitStatistics=fit3; 

run;

data parcov12;
parameter="cov12";
estimate=0;
run;

data initialvalue2;
set est3 parcov12;
run;
/* covariance */


PROC nlmixed DATA=one METHOD=GAUSS;
parms/data=initialvalue2;

bounds vara varb phi>0;

/* log-likelihood */
teta = alpha0 + a  + alpha1*Time + alpha2*Treatment + alpha3*Baseline&var2;
expteta = exp(teta);
p = expteta / (1 + expteta);


/* For zero Faecalibacterium */
if &var=0 then loglik=log(1-p);

/* For positive Faecalibacterium */
if &var>0 then do;
	mu = exp(beta0 + b + beta1*Time + beta2*Treatment + beta3*Baseline&var2)/(1+exp(beta0 + b + beta1*Time + beta2*Treatment + beta3*Baseline&var2));
	loglik = log(p) + (mu*phi-1)*log(&var2) + ((1-mu)*phi-1)*log(1-&var2) + lgamma(phi) - lgamma(mu*phi) - lgamma((1-mu)*phi);
end;

/* fit the above modle */
model &var ~ general(loglik);
random a b~ normal([0,0],[vara,cov12,varb]) subject=Subject;
ods output ParameterEstimates=est3 FitStatistics=fit3; 

run;

%mend;

%TPM(al, Alistipes); %TPM(ba, Bacteroides); %TPM(bi, Bifidobacterium); %TPM(clo, Clostridium); %TPM(col, Collinsella); %TPM(cop, Coprobacillus); 
%TPM(di, Dialister); %TPM(do, Dorea); %TPM(es, Escherichia); %TPM(eu, Eubacterium); %TPM(fa, Faecalibacterium); %TPM(h, Haemophilus); 
%TPM(l, Lactobacillus); %TPM(p, Parabacteroides); %TPM(ro, Roseburia); %TPM(ru, Ruminococcus); %TPM(s, Streptococcus); %TPM(v, Veillonella); 

