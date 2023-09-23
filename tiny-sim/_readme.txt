Near-simplest possible multispecies coalescent simulation: 
- 2 samples from each of 2 species in just one gene
- theta and lambda chosen so that deep coalescence is, 
    for all intents and purposes, impossible

Used to demonstrate that marginal likelihoods computed using SMC
are the same (excepting small sampling error) as those estimated
using the LoRaD method. LoRaD requires all gene trees and the 
species tree to be identical for all particles used in the marginal
likelihood estimation; hence, the choice of theta and lambda.

Instructions:

1. Compile dub with SAVE_PARAMS_FOR_LORAD defined (in conditionals.hpp)
2. Navigate into tiny-sim/simulate and run dub to generate simulated.nex
    and files for the true gene tree and true species tree
3. Navigate into tiny-sim/gene-marglike and run dub and then loradML
    This tests estimation of the marginal Felsenstein likelihood 
    conditioning on the species tree and theta: p(D|S,theta)
4. Navigate into tiny-sim/species-marglike and run dub and then loradML
    This tests estimation of the marginal coalescent likelihood
    conditioning on theta and lambda: p(G|theta, lambda)


