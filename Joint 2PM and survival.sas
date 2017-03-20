data cost;
input id anycost z1 month cost;
cards;
1 0 1 0 0
1 1 1 1 15
1 1 0 2 35
2 1 0 0 24
2 0 1 1 0
2 1 0 2 56
2 1 1 3 32
...
;
data cost2;
set cost;
bb=0;
run;
data fup;
input id z follow_time death;
cards;
1 1 12.3 0
2 0 24 1
3 1 6.7 1
...
;
run;
data fup;
set fup;
bb=1;
run;
* Obtain the 10 quantiles;
proc univariate data=fup noprint;
var follow_time; 
output out=quant pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=q; 
where death=1;
run;
data quant;
set quant;
bb=1;
run;
* Merge data with the quantiles;
data fup2;
merge fup quant;
by bb;
run;
data fup3;
set fup2;
array quant {11} q0 q10 q20 q30 q40 q50 q60 q70 q80 q90 q100;

array dur {10} dur1-dur10;

array event {10} event1-event10;

do i=1 to 10;
	dur{i}=0;
end;

do i=1 to 10;
	event{i}=0;
end;

do i=2 to 11;		
	if follow_time<=quant{i} then do;
		dur{i-1}=follow_time-quant{i-1};	/* duration in each death event quantiles */
		event{i-1}=death;					/* indicator of death event in each interval */
		i=11;
	end;
	else dur{i-1}=quant{i}-quant{i-1};
end;
run;
* Prepare the final data for analysis;
data cost_all;
set cost2 fup3;
if cost=. then cost=follow_time;
run;
proc sort data=cost_all;
by id bb month;
run;
proc nlmixed data=cost_all  qpoints=5;
parms alpha0=-1 alpha1=1 alpha2=.1 beta0=-1 beta1=0.5 beta2=.1  
 var1=1 var2=1 cov12=0.5 sigma=2 gamma1=.5 gamma2=.5 delta1=1
r01=0.01 r02=0.01 r03=0.015 r04=0.02 r05=0.02 r06=0.03 r07=0.03 r08=0.03 r09=0.04 r10=0.04;

bounds r01 r02 r03 r04 r05 r06 r07 r08 r09 r10 var1 var2 sigma >=0;

base_haz=r01 * event1 + r02 * event2 + r03 * event3 + r04 * event4 + r05 * event5 
+ r06 * event6 + r07 * event7 + r08* event8 +r09 * event9 + r10 * event10;

cum_base_haz=r01 * dur1 + r02 * dur2 + r03 * dur3 + r04 * dur4 + r05 * dur5 
+ r06 * dur6 + r07 * dur7 + r08* dur8 +r09 * dur9 + r10 * dur10;

teta=alpha0 + a + alpha1*z1 + alpha2 * month;
expteta=exp(teta);
p=expteta / (1+expteta);

mu1=beta0 + beta1 * z1 + beta2 * month + b ;

mu2= delta1 * z + gamma1 * a + gamma2 *b ;

* Cost data;
if bb=0 then do;
	if anycost=0 then loglik=log(1-p);
	if anycost=1 then loglik=log(p)-.5*((log(cost)-mu1)/sigma)**2-.5*log(sigma**2);
end;

* Survival data;
if bb=1 then do;
	if death=1 then loglik=log(base_haz) + mu2 -exp(mu2) * cum_base_haz;
	if death=0 then loglik=-exp(mu2) * cum_base_haz;
end;

model cost ~ general(loglik);
random a b ~ normal([0, 0], [var1, cov12, var2]) subject=id;
run;
