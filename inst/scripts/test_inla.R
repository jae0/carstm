
# some tests of INLA and CARSTM in case internals of INLA change 

# how INLA treats offsets are not always clear or consistent
# CARSTM tries to handle them: the following tests should succeed so long as INLA does not change base behaviour

library(carstm)
# loadfunctions("carstm")

# these are simple regression tests
carstm_test_inla("gaussian")

carstm_test_inla("poisson")

carstm_test_inla("binomial")
 
# ---- 
# testing of slightly more complex, "bym" models

# NOTE key difference:  .. trying to account for this divergent behaviour causes a lot of the complexity within carstm
# classical mode: predictions on user scale, incorporating offsets 
# experimental mode:  predictions are on user scale, NOT incorporating offsets  


carstm_test_inla("besag_poisson_nbinomial")



carstm_test_inla("bym2_poisson_nbinomial_random_effect")


