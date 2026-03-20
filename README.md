# Jolly-Seber-N-Prior-DA
Jolly-Seber MCMC samplers using N-prior data augmentation

This repository contains efficient MCMC samplers for Jolly-Seber models using N-prior data augmentation where we put Poisson priors directly on the
entries--the starting population size and yearly recruits, rather than Bernoulli data augmentation where Poisson assumptions only hold as
the data augmentation size M goes to infinity. For recruits specifically, this means the variance in the number of yearly recruits declines across primary occasions.
Poisson priors on entries avoids this effect and allows for substantially more efficient MCMC, especially compared to existing approaches
that allow for per capita recruitment rate to be estimated (Chandler and Clark 2014). 

N-prior data augmentation for closed models is explained/demonstrated here: 
https://github.com/benaug/SCR-N-Prior-Data-Augmentation

and for multisession, see here:
https://github.com/benaug/SCR_Dcov_Multisession

In both of these, we assume N ~ Poisson(lambda) and then we allocate the N real individuals and N0=M-N augmented individuals to 
the vector z[1:M] at random where order does not matter. There are (M choose N) ways to do this allocation, so the prior for 
[z[1:M] | N,M] is 1/(M choose N).

Moving to open populations, we assume 

N[1] ~ Poisson(lambda) and 

N.recruit[g] ~ Poisson(N[g]*gamma[g]) for g = 1,..., n.primary - 1.

Now, the prior for z.super[1:M], which indicates if an individual is ever in the population, must consider the number of ways to allocate
each Poisson RV to the indices of z.super[1:M], which is the multinomial coefficient of size M with partition sizes

N[1], N.recruit[1:(n.primary-1)], 

and N0=M-N[1]-N.recruit[1:(n.primary-1)]. 

So the prior [z.super[1:M] | N[1],N.recruit[1:(n.primary-1)],M] is


(N[1]!N.recruit[1]!...N.recruit[n.primary-1]!N0!)/ M!, the inverse multinomial coefficient.

This repository also has sex-specific versions with sex-specific population dynamics. For these models, the prior is analogous to that
described above except we now have sex-specific partitions for N[1] and N.recruit. There is a single N0 off class, not sex-specific.

My first attempt at this approach (Jolly-Seber Github repo) mistakenly used the same
prior as for closed populations, which appeared to work as expected with the larger population sizes used to test the code, but not with smaller population sizes,
especially for sex-specific models where there are more possible partitions in the z.super prior. 

More notes: individuals with z.super[i]=0 do not have a z vector describing entry year and survival (set to all 0). These z states are simulated when proposing to update
N.super/z.super. In SCR versions, the z.super[i] individuals maintain activity centers for convenience like all other DA approaches in the past,
but these could be removed and simulated when proposing to add as well. Finally, for sex-specific versions, I keep sex in the model for z.super=1 individuals
but they do not mean anything and this is done for convenience. The BUGS code needs sexes of z.super=0 individuals to compile and to compute
the correct survival probabilities when turning on individuals. I simulate a new sex when turning individuals on because there is no model for
individual sex when z.super=0.

Model versions:

1. JS: nonspatial, 1 continuous survival covariate
2. JS-SCR: same as 1 but spatial, fixed activity centers
3. JS-SexPopDy: sex-specific population dynamics, male and female survival and recruitment parameters
4. JS-SCR-SexPopDy: same as 3 but spatial, fixed activity centers.
5. JS-SCR-Dcov: same as 2 with habitat mask and/or density covariates
6: JS-SCR-Dcov-SexPopDy: same as 4 with habitat mask and/or density covariates
7. JS-Typical: This is a nonspatial version of the Jolly-Seber approach of Chandler and Clark (2014) that 
considers per capita recruitment for comparison.

These SCR version below consider activity center movement. A notable difference between the ones below and those above
is that I gate the activity center likelihood by z.super so that they are not in the model when z.super=0. They are set
to "0" to indicate they are turned off. Then, a new s trajectory is proposed when turning on a z.super index. This
improved mixing of the between year movement scale parameter for a few simulated data sets where I made the comparison.
Another note is that for mobile activity centers, you generally need pretty good SCR data to estimate the movement
parameter well--many individuals, many survival events that are documented (same inds captured in consecutive years),
larger trapping arrays, etc. You will also have to run the model for more MCMC iterations to get a good effective sample
size for the movement parameter. 

8. JS-SCR-mobileAC: same as 2 with BVN Markov activity center movement (truncated by state space boundary)
9. JS-SCR-SexPopDy-mobileAC: same as 8 but sex-specific population dynamics, detection, and movement parameters

I will upload mobile AC versions that work with inhomogenous density later.
