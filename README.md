# Jolly–Seber N-Prior Data Augmentation

Efficient MCMC samplers for Jolly–Seber models using N-prior data augmentation.

## Overview

This repository contains efficient MCMC samplers for Jolly–Seber models using **$N$-prior data augmentation**. Poisson priors are placed directly on population entry counts: the initial population size and the number of recruits entering before each subsequent primary occasion. Recruitment is per capita as a function of realized abundance in the previous primary occasion. This provides an exact Poisson recruitment specification for the per-capita recruitment model of Chandler and Clark (2014).

The repository also contains implementations of several existing Jolly–Seber data augmentation formulations. For each, one can add custom whole-history block Gibbs samplers based on the Type I and Type II updates of Wu et al. (2021), which build on the blocked latent-history sampling approach of Dupuis and Schwarz (2007). The N-prior models use more generic Gibbs updates of the entry and exit components of each individual's latent history from their full conditional distributions. This more general approach is less computationally efficient than a highly tailored Wu et al.-type update, but was only moderately slower in the limited comparisons I did and is considerably easier to adapt when the model structure or observation process is changes. These other Jolly-Seber formulations are described at the end of the readme below.

## Related repositories

N-prior data augmentation for closed-population models is explained and demonstrated here:

- [SCR N-Prior Data Augmentation](https://github.com/benaug/SCR-N-Prior-Data-Augmentation)

A multisession implementation is available here:

- [SCR_Dcov_Multisession](https://github.com/benaug/SCR_Dcov_Multisession)

If you do not have individual covariates, you can fit Jolly-Seber models without using data augmentation:

 - [Jolly-Seber-SS](https://github.com/benaug/Jolly-Seber-SS)

Here is a model that uses N-prior data augmentation as the backbone of a spatial integrated population model that allows telemetry-informed survival. You can mix and match two different SCR observation processes in each primary occasion, and/or have primary occasions with no observations.

- [Spatial-IPM-Telemetry](https://github.com/benaug/Spatial-IPM-Telemetry)

Here is open population spatial mark resight using the marginalized observation model. It is a spatial IPM where you can mix and match marking and sighting observation processes in each primary occasion and/or have primary occasions with no observations.

- [Spatial-Mark-Resight-Open-Marginal](https://github.com/benaug/Spatial-Mark-Resight-Open-Marginal)

Here is the conditional observation model version of the above. Mixes worse, but allows non-Poisson observation models for sighting data.

- [Spatial-Mark-Resight-Open-Conditional](https://github.com/benaug/Spatial-Mark-Resight-Open-Conditional)

Here is an open population version of the 2-flank spatial partial identity model

- [2-flank-SPIM-Jolly-Seber](https://github.com/benaug/2-flank-SPIM_Jolly_Seber)

## Closed-population N-prior data augmentation

For a closed population, assume

$$
N \sim \mathrm{Poisson}(\lambda),
$$

where $N$ is the number of real individuals and $M$ is the data augmentation size. The number of augmented individuals is

$$
N_0 = M-N.
$$

The $N$ real individuals and $N_0$ augmented individuals are randomly allocated among the elements of

$$
\mathbf{z} = (z_1,\ldots,z_M).
$$

Because of this random allocation, the ordering of individuals within the real and augmented classes does not matter. Therefore, 
there are

$$
\binom{M}{N}
$$

possible allocations, and the conditional prior for $\mathbf{z}$ is

$$
\Pr(\mathbf{z}\mid N,M) = \binom{M}{N}^{-1} = \frac{N!N_0!}{M!}.
$$

A custom Metropolis-Hastings update is required to jointly update $N$ and $\mathbf{z}$.

## Open-population formulation

Let $G$ denote the number of primary occasions. The initial population size is assigned the prior

$$
N_1 \sim \mathrm{Poisson}(\lambda).
$$

The number of recruits entering between primary occasions $g$ and $g+1$ is assigned the prior

$$
B_g \sim \mathrm{Poisson}(N_g\gamma_g),
\qquad
g=1,\ldots,G-1,
$$

where $B_g$ corresponds to `N.recruit[g]` in the code and $\gamma_g$ is the per capita recruitment rate. The total number of individuals that ever enter the population is

$$
N^{\mathrm{super}} = N_1+\sum_{g=1}^{G-1}B_g.
$$

The number of augmented individuals that never enter is

$$
N_0=M-N^{\mathrm{super}}.
$$

To describe the allocation of the $M$ indices, define the entry-occasion variable

$$
e_i \in \{0,1,\ldots,G\},
$$

where

$$
e_i=
\begin{cases}
0, & \text{individual } i \text{ never enters},\\
1, & \text{individual } i \text{ is present during the first primary occasion},\\
g, & \text{individual } i \text{ enters at primary occasion } g,
\quad g=2,\ldots,G.
\end{cases}
$$

The corresponding cohort counts are

$$
N_0=\sum_{i=1}^{M} I(e_i=0),
$$

$$
N_1=\sum_{i=1}^{M} I(e_i=1),
$$

and

$$
B_g=\sum_{i=1}^{M} I(e_i=g+1),
\qquad
g=1,\ldots,G-1.
$$

Thus, $e_i=1$ indicates membership in the initial population, and $e_i=g$ indicates recruitment before primary occasion $g$, for $g=2,\ldots,G$. The value $e_i=0$ indicates that the augmented individual never enters the population. In the code, `z.start[i]` is the entry-occasion variable $e_i$.

The binary superpopulation indicator is then a derived variable:

$$
z_i^{\mathrm{super}}=I(e_i\gt 0)
$$

Conditional on the cohort counts, the number of ways to allocate the $M$ augmented indices among the initial population, the $G-1$ recruitment cohorts, and the never-entered class is

$$
\frac{
M!
}{
N_0!N_1!\displaystyle\prod_{g=1}^{G-1}B_g!
}.
$$

Therefore, the prior probability of a particular entry-cohort allocation is

$$
\Pr\left(
\mathbf{e}
\mid
N_1,B_1,\ldots,B_{G-1},M
\right)=\frac{
N_0!N_1!\displaystyle\prod_{g=1}^{G-1}B_g!
}{
M!
}.
$$

This is the inverse multinomial coefficient. The entry occasion also determines the initial portion of the latent population-state history. For an individual with $e_i\gt 0$,

$$
z_{i,g}=0\qquad\text{for } g\lt e_i,
$$

and

$$
z_{i,e_i}=1.
$$

Subsequent states are governed by the survival model. For example,

$$
\mathbf{z}_i=(0,0,1,1,0)
$$

implies that $e_i=3$, or equivalently `z.start[i] = 3`: individual $i$ entered at primary occasion 3, remained alive through occasion 4, and was no longer alive at occasion 5.

The entry-occasion variable describes only when an individual enters. Survival after entry is governed separately by the survival model.

To update the latent population states, I use three custom samplers:

1. **Detected individuals:** Entry and exit occasions are updated separately using Gibbs samplers.

2. **Undetected individuals currently in the superpopulation:** For individuals with $z^{\mathrm{super}}_i = 1$ but no detections, an entry occasion and complete survival history are proposed from the process model and accepted or rejected jointly using a Metropolis-Hastings update.

3. **Superpopulation size:** To update the number of individuals that ever enter the population, the sampler proposes adding or removing one undetected individual at a time. An addition changes $e_i=0$ to $e_i\gt 0$ and jointly proposes the individual's entry occasion and survival history; consequently, $z_i^{\mathrm{super}}$ changes from 0 to 1. A removal changes $e_i \gt 0$ to $e_i=0$ and removes the associated entry and survival history; consequently, $z^{\mathrm{super}}_i$ changes from 1 to 0.


## Sex-specific population dynamics

Some model versions include sex-specific population dynamics, with separate initial population sizes and recruitment processes for females and males. Let $N_1^F$ and $N_1^M$ denote the numbers of females and males in the initial population, respectively. Let $B_g^F$ and $B_g^M$ denote the numbers of female and male recruits entering before primary occasion $g+1$, for $g=1,\ldots,G-1$. Define the recruitment vectors as

$$
\mathbf{B}^F=(B_1^F,\ldots,B_{G-1}^F)
\qquad\text{and}\qquad
\mathbf{B}^M=(B_1^M,\ldots,B_{G-1}^M).
$$

The total number of individuals that ever enter the population is

$$
N^{\mathrm{super}}=N_1^F+N_1^M+\sum_{g=1}^{G-1}(B_g^F+B_g^M),
$$

and the size of the single shared never-entered class is

$$
N_0=M-N_{\mathrm{super}}.
$$

Thus, the $M$ augmented indices are allocated among the female and male initial populations, the female and male recruitment cohorts, and the single never-entered class. The number of possible allocations is

$$
\frac{M!}{N_0!\left(N_1^F\right)!\left(N_1^M\right)!\displaystyle\prod_{g=1}^{G-1}\left(B_g^F\right)!\left(B_g^M\right)!}.
$$

Therefore, the prior probability of a particular joint allocation of entry cohort and sex is

$$
\Pr\left(\mathbf{e},\mathbf{sex}\mid N_1^F,N_1^M,\mathbf{B}^F,\mathbf{B}^M,M\right)=\frac{N_0!\left(N_1^F\right)!\left(N_1^M\right)!\displaystyle\prod_{g=1}^{G-1}\left(B_g^F\right)!\left(B_g^M\right)!}{M!},
$$

where $\mathbf{sex}$ is the vector of individual sex states, with $\text{sex}_i=1$ for male, $\text{sex}_i=2$ for female, and $\text{sex}_i=0$ for individuals that never enter the population ($e_i=0$ and $z_i^{\mathrm{super}}=0$). This is the inverse multinomial coefficient for the joint allocation of the $M$ augmented indices among entry cohorts, sex classes, and the single never-entered class.

The three custom samplers operate similarly in the sex-specific models, except that samplers 2 and 3 also jointly propose the individual’s sex.

## Model versions

### Fixed activity centers

1. **JS**

   Nonspatial model with one continuous survival covariate.

2. **JS-SCR**

   Spatial version of model 1 with fixed activity centers.

3. **JS-SexPopDy**

   Nonspatial model with sex-specific population dynamics, including male- and female-specific survival and recruitment parameters.

4. **JS-SCR-SexPopDy**

   Spatial version of model 3 with fixed activity centers.

5. **JS-SCR-Dcov**

   Spatial version of model 2 with a habitat mask, density covariates, or both.

6. **JS-SCR-Dcov-SexPopDy**

   Spatial version of model 4 with a habitat mask, density covariates, or both.


  
### Mobile activity centers

The SCR models below allow activity centers to move among primary occasions.

For these models, the activity center model is gated by $z_i^{\mathrm{super}}$. For individuals with $z^{\mathrm{super}}_i=0$, the activity centers (and associated availability and use distributions in RSF movement models) are set to zero and do not contribute to the model density. When an augmented individual is activated, a complete activity center trajectory is proposed jointly with its entry and survival history. This avoids evaluating the movement model for inactive augmented individuals and, in several simulated datasets, improved mixing of the between primary occasion movement parameter relative to retaining active activity center trajectories for all augmented individuals; however, I have not conducted a formal comparison.

Mobile activity center models generally require informative SCR data to estimate movement parameters reliably, e.g.,

- many observed individuals
- many documented survival transitions, with the same individuals detected in consecutive primary occasions
- sufficiently large trapping arrays
- adequate spatial recaptures within and among primary occasions
- telemetry data

These models may also require substantially more MCMC iterations to obtain adequate effective sample sizes for movement parameters and possibly survival and recruitment parameters if they are sensitive to the magnitude of movement.

7. **JS-SCR-mobileAC**

   Spatial version of model 2 with bivariate normal Markov activity center movement. The movement distribution is truncated by the state-space boundary.

8. **JS-SCR-SexPopDy-mobileAC**

   Version of model 8 with sex-specific population dynamics, detection parameters, and movement parameters.

9. **JS-SCR-Dcov-mobileAC**

   Version of model 8 with an inhomogeneous density model for activity centers during the first primary occasion and resource selection during subsequent activity center movement.

Let $I_c$ denote the habitat-mask indicator, where $I_c=1$ if cell $c$ is included in the state space and $I_c=0$ otherwise. Let $x_c$ denote the spatial covariate for cell $c$, and let $C$ denote the total number of grid cells.

With a single habitat covariate, during the first primary occasion, the relative intensity of activity centers in cell $c$ is

$$
\lambda_c=I_c\exp(\beta^Dx_c),
$$

where $\beta^D$ is the density covariate coefficient. The probability that an individual's initial activity center occurs in cell $c$ is

$$
\pi_c=\frac{\lambda_c}{\displaystyle\sum_{k=1}^{C}\lambda_k}.
$$

An initial activity center cell is drawn from this categorical distribution, after which the continuous activity center coordinates are drawn uniformly within the selected cell. During each subsequent primary occasion, activity centers follow a Markov movement model with resource selection. Conditional on the previous activity center $\mathbf{s}_{i,g-1}$, the availability probability for cell $c$ is

$$
a_{i,g,c}=\int_{\mathcal{A}_c}
\mathrm{BVN}\left(
\mathbf{v}
\mid
\mathbf{s}_{i,g-1},
\sigma_{\mathrm{move}}^2\mathbf{I}_2
\right)
\,d\mathbf{v},
$$

where $\mathcal{A}_c$ is the area of cell $c$, $\mathbf{v}=(v_x,v_y)$ is a possible continuous activity center location within that cell, and $\mathrm{BVN}(\mathbf{v}\mid\boldsymbol{\mu},\sigma^2\mathbf{I}_2)$ denotes an isotropic bivariate normal density with mean $\boldsymbol{\mu}$, equal coordinate variances $\sigma^2$, and zero covariance. The matrix $`\mathbf{I}_{2}`$ is the $`2\times2`$ identity matrix, and $`\sigma_{\mathrm{move}}`$ is the scale of activity center movement between primary occasions.

The resource selection weight for cell $c$ is

$$
r_c=I_c\exp(\beta^{\mathrm{RSF}}x_c),
$$

where $\beta^{\mathrm{RSF}}$ is the resource selection coefficient. The resulting use probability, or probability that individual $i$ selects cell $c$ during primary occasion $g$, is

$$
u_{i,g,c}=\frac{a_{i,g,c}r_c}
{\displaystyle\sum_{k=1}^{C}a_{i,g,k}r_k}.
$$

Thus, cell selection depends jointly on proximity to the previous activity center and selection for the spatial covariate. Because $r_c=0$ when $I_c=0$, only cells within the habitat state space can be selected. After a cell is selected, each activity center coordinate is drawn from its normal movement distribution centered on the previous activity center and truncated to the boundaries of the selected cell. This produces a continuous activity center location rather than assigning the individual to the cell centroid.

This movement formulation is a cell-integrated version of the movement model described by [Bischof et al. (2020)](https://doi.org/10.1073/pnas.2011383117). The initial point process formulations differ: Bischof et al. use an inhomogeneous binomial point process conditional on a fixed augmented population size, whereas this model assigns a Poisson distribution to initial abundance, so the activity centers of individuals present during the first primary occasion form an inhomogeneous Poisson point process.

The movement models are otherwise equivalent when the spatial covariate is constant within each grid cell. Bischof et al. define the transition density as the normalized product of an isotropic Gaussian movement kernel centered on the previous activity center and a resource selection weight. Integrating that density over cell $c$ gives a cell-selection probability proportional to $a_{i,g,c}r_c$, matching the cell-level probability used in this implementation. Conditional on selecting cell $c$, drawing the new activity center from the Gaussian movement distribution truncated to that cell recovers the same continuous transition density as Bischof et al. Thus, the two approaches differ in computational representation rather than in the underlying movement model.

This normalized product of a movement kernel and resource selection weight was used in the spatially explicit habitat-selection model of [Rhodes et al. (2005)](https://doi.org/10.1890/04-0912) and is also the standard structure of step selection functions, in which a resource-independent movement kernel is multiplied by a resource selection function and normalized over possible endpoints ([Forester et al. 2009](https://doi.org/10.1890/08-0874.1)).

In this implementation, the bivariate normal probability mass within each rectangular cell is calculated analytically from differences of univariate normal cumulative distribution functions, and calculations are restricted to cells with non-negligible probability mass. Storing the availability distributions allows them to be reused when updating $\beta^{\mathrm{RSF}}$. These computational savings come at the cost of increased RAM use.


10. **JS-SCR-Dcov-mobileAC-patchy**

    Modification of `JS-SCR-Dcov-mobileAC` that applies the same initial density and activity center movement model to a patchy habitat state space.

This version differs from model 10 in two respects:

1. The simulator uses a patchy habitat state space containing gaps that make local activity center proposals more difficult.
2. A third activity center sampler is included to improve proposals across gaps in the habitat state space.

A discrete proposal based directly on each individual's cell-level use distribution would provide a more direct way to move across habitat gaps but would be computationally expensive. The continuous proposal used here should work well in many applications, but its performance should be evaluated through simulation for each general model configuration.


For models 10 and 11, I need to add code that checks whether all initial activity center log probabilities are finite. If any are nonfinite, the sampler may not be able to recover and may cause R to crash. I did not encounter this problem in simulated datasets when initializing  $\sigma_{\mathrm{move}}$ at its true value. It is more likely to occur with real datasets in which the model assumptions are not perfectly met, such as when a small number of observed movements are much larger than the others. If this occurs, increase the initial value of $\sigma_{\mathrm{move}}$.


### Alternative Data Augmentation Approaches

Here are several other data augmentation approaches for Jolly-Seber models. This first set is set up with individual survival covariates. Except for the marginalized formulation, these implementations use custom block Gibbs samplers based on the Type I and Type II recursions of Wu et al. (2021). Rather than updating latent states one occasion at a time or explicitly evaluating every possible latent history, the Type I and Type II recursions propagate the relevant conditional likelihoods through time and update the entry and exit portions of an individual's latent history as separate blocks. The Schwarz-Arnason implementations directly adapt the recursions of Wu et al., whereas other formulations require model-specific modifications to accommodate their different process model specifications. Where applicable, superpopulation membership is also updated using the corresponding marginalized block probabilities. The resulting updates are hard coded for the observation models used in these implementations: changing an observation model, or adding an additional observation process, generally requires modifying the block sampler so that the appropriate observation-likelihood contributions are included in the recursions. The specific latent-state representation also differs among models, and the Chandler and Clark (2014) formulation requires an adapted Type I/Type II update because recruitment probabilities depend on realized abundance and therefore couple individual histories. These recursive block Gibbs samplers improve computational efficiency, with especially large gains for the Chandler and Clark model. Finally, the marginalized formulation integrates over the latent entry and survival states directly and therefore does not require a latent-state sampler.

11. **JS-RO**

    Nonspatial implementation of the Restricted Occupancy Jolly–Seber approach of Royle and Dorazio (2008). Can choose original priors or priors from Dorazio 2020. There are two sets of files, one with "Dorazio" in file names. I encountered mixing problems with the Dorazio priors, with realized recruits becoming stuck at 0. The corresponding Beta priors on recruitment probabilities have infinite density at 0 and induce substantial prior probability for zero realized recruitment.

12. **JS-CC**

    Nonspatial implementation of the Jolly–Seber approach of Chandler and Clark (2014), an extension of Restricted Occupancy that estimates per capita recruitment as a function of realized abundance. The N-prior data augmentation models in this repository use the corresponding exact Poisson recruitment specification.
   
13. **JS-SA-sequential**

    Nonspatial implementation of the Schwarz-Arnason Jolly–Seber approach of Royle and Dorazio (2008). This version uses the conditional entry probabilities and sequential latent-state parameterization of Royle and Dorazio. This testscript and the next two are set up with the same parameter values and seed for comparison.
   
14. **JS-SA-entryOccasion**

    Nonspatial implementation of the Schwarz–Arnason Jolly–Seber approach of Royle and Dorazio (2008). This is the same model as JS-SA-sequential, but with a different latent-state parameterization. Instead of representing entry through conditional entry probabilities in a sequential latent-state model, entry occasion is represented directly as `e[i] ~ dcat(pi[])`, where `pi` has a Dirichlet prior, and survival is modeled forward from the entry occasion. Augmented individuals outside the superpopulation have potential entry and survival trajectories that do not contribute to the realized population. This can be viewed as a data-augmentation implementation of the original Schwarz and Arnason (1996) entry-time formulation. Related data-augmentation approaches that explicitly represent latent entry/birth and exit/death processes have been used by Schofield and Barker (2011), Matechou et al. (2016), and Wu et al. (2021).

15. **JS-SA-marginal**

    Nonspatial implementation of the Schwarz-Arnason Jolly–Seber approach of Royle and Dorazio (2008). Unlike the other implementations here, this version marginalizes the latent states out of the likelihood with a 3-state forward algorithm (not yet entered, alive, dead) rather than using a latent-trajectory Gibbs sampler. It retains the M augmented individuals to accommodate an individual survival covariate without numerical integration. N, B, and N.super are recovered each iteration by forward-filtering backward-sampling the trajectories from their full conditionals. Marginalizing the latent states prevents use of the conjugate updates available for some parameters in the explicitly augmented versions, and in the comparisons I have made this version appears less computationally efficient. Marginalization can be considerably faster when capture histories can be aggregated, for example when there are no individual random effects. Marginalization is not computationally feasible when recruitment is per capita as a function of realized abundance because individuals are then coupled through the population process. This model and marginalization approach have also been implemented in Stan, e.g., https://github.com/mhollanders/cr-in-stan/blob/main/stan/js-rng.stanfunctions, where marginalization is required because Stan cannot sample discrete parameters. HMC may be more efficient, but I have not made this comparison.

Below are models that consider per capita recruitment as a function of expected abundance, rather than realized abundance as assumed in Chandler and Clark (2014) and the N-prior data augmentation models here (I prefer the use of realized abundance conceptually). For efficient MCMC, I removed the individual survival covariate from these implementations. Individual survival covariates can be included, but computing expected abundance then requires numerical integration, which can substantially slow down the MCMC.

Expected recruitment during an interval is proportional to expected abundance at the beginning of the interval. Two parameterizations are provided for the Restricted Occupancy and sequential Schwarz-Arnason implementations. Files without `V2` use a forward parameterization in which the initial entry probability (`psi`) supplies the population scale and expected abundance and recruitment are propagated forward through time. This substantially reduces dependencies in the NIMBLE model graph and therefore compiles much faster. A disadvantage is that the posteriors of `psi` and the per capita recruitment rate (`gamma`) can be strongly correlated, which can reduce MCMC efficiency. However, an AF slice block sampler for `psi` and `gamma` is very efficient, at least when a single `gamma` is shared across intervals.

Files with `V2` use the original parameterization, in which a relative expected demographic trajectory is constructed and normalized across entry occasions, with a separate superpopulation-scale parameter supplying the absolute population scale. This creates more dependencies among demographic parameters, entry probabilities, and latent states, so these models can take substantially longer to compile. Once compiled, however, they run quickly and generally avoid the strong posterior correlation between `psi` and `gamma` that can occur in the forward parameterization.

Even though V1 compiles faster than V2, these models compile much slower than the realized abundance per capita recruitment models (Chandler and Clark and NPDA models here), particularly as the number of primary occasions grow. This is a second reason to just use the realized abundance per capita models!

In limited comparisons (V1, not V2), the Restricted Occupancy implementation was more efficient than the Schwarz-Arnason

16. **JS-RO-perCapitaExpectedN**

    Restricted-occupancy implementation of the per capita expected-abundance model. Expected recruitment determines the distribution of entry across primary occasions and the initial and conditional entry probabilities of the restricted-occupancy formulation. There is no separate superpopulation membership indicator: an augmented individual outside the superpopulation is represented directly by a trajectory that never enters.

17. **JS-SA-sequential-perCapitaExpectedN**

    Schwarz-Arnason implementation of the per capita expected-abundance model using a sequential latent-state formulation. Relative expected recruitment determines the distribution of entry across primary occasions, which is represented through the conditional entry probabilities of the sequential state process.
    
