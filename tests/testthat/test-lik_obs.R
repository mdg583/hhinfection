# test_that("bcoef generates correct coefficients", {
#   expect_equal(.test_bcoef(),1)
# })

# test_that("count_bits calculates bits correctly", {
#   expect_equal(.test_count_bits(),1)
# })

# test_that("rbool times", {
#   expect_equal(.test_rbool_times(),1)
# })

# lprob_sec_obs_pri(2,0,1,2,0.02,Sn=.99,Sp=.99)

# test_that("prob of household observation based on actual is computed correctly", {
#   expect_equal(.test_lprob_obs_act(),1)
# })
#
# test_that("Test1: brute force likelihood on a single household without Sn Sp", {
#   # Worked out by hand for testing
#   obs = c(1,0)
#   prim = c(0.4,0.7)
#   hh = c(1,1)
#   q = 0.1
#   #message("Testing prob_obs_hhs")
#   lik = prob_obs_hhs(obs,prim,hh,q)
#   #message(lik)
#   expect_equal(lik, log(0.108))
# })
#
# test_that("Test2: brute force likelihood on a single household with Sn Sp", {
#   # Worked out using a computer for testing
#   obs = c(1,0)
#   prim = c(0.4,0.7)
#   hh = c(1,1)
#   q = 0.1
#   Sn = 0.9
#   Sp = 0.8
#   #message("Testing prob_obs_hhs with Sn Sp")
#   lik = prob_obs_hhs(obs,prim,hh,q,Sn,Sp)
#   message(lik)
#   expect_equal(lik, -1.93669276)
# })
#
# test_that("Test3: sim likelihood on a single household with Sn Sp", {
#   # Worked out using a computer for testing
#   obs = c(1,1,0,0,0,0)
#   prim = c(0.4,0.7,0.1,0.08,0.12,0.2)
#   hh = c(1,1,1,1,1,1)
#   q = 0.1
#   Sn = 0.9
#   Sp = 0.8
#   #message("Testing prob_obs_hhs with Sn Sp")
#   x = vapply(1:1000,function(i){
#     prob_obs_hhs(obs,prim,hh,q,Sn,Sp,M1=100,T1=50)
#   },0)
#   x_bf = prob_obs_hhs(obs,prim,hh,q,Sn,Sp)
#   tpval = t.test(x,mu=x_bf)$p.value
#   expect_gt(tpval, 0.001)
# })
#
# test_that("Test4: sim likelihood on a single household with Sn Sp", {
#   n = 15;
#   prim = exp(runif(n,log(0.01),log(0.3)))
#   obs = rbinom(n,1,prim)
#   hh = rep(1,n)
#   q = 0.1
#   Sn = 0.9
#   Sp = 0.8
#   x = vapply(1:200,function(i){
#     prob_obs_hhs(obs,prim,hh,q,Sn,Sp,M1=30,T1=1)
#   },0)
#   x_bf = prob_obs_hhs(obs,prim,hh,q,Sn,Sp)
#   exp(x_bf)
#   summary(exp(x))
#   hist(exp(x),breaks=50)
#   abline(v=exp(x_bf))
#   tpval = t.test(exp(x),mu=exp(x_bf))$p.value
#   tpval
#   expect_gt(tpval, 0.001)
# })
#
# test_that("Test5: Secondary Simulation", {
#   # Worked out using a computer for testing
#   obs = c(1,1,0,0,0,0)
#   prim = c(0.4,0.7,0.1,0.08,0.12,0.2)
#   hh = c(1,1,1,1,1,1)
#   q = 0.1
#   Sn = 0.9
#   Sp = 0.8
#   #message("Testing prob_obs_hhs with Sn Sp")
#   x = vapply(1:500,function(i){
#     prob_obs_hhs(obs,prim,hh,q,Sn,Sp,M1=0,M2=4)
#   },0)
#   x_bf = prob_obs_hhs(obs,prim,hh,q,Sn,Sp)
#   hist(exp(x),breaks=50)
#   abline(v=exp(x_bf))
#   tpval = t.test(exp(x),mu=exp(x_bf))$p.value
#   expect_gt(tpval, 0.001)
# })
#
# test_that("Test6: sim likelihood on a single household with Sn Sp", {
#   # Worked out using a computer for testing
#   n = 15;
#   prim = exp(runif(n,log(0.01),log(0.3)))
#   obs = rbinom(n,1,prim)
#   hh = rep(1,n)
#   q = 0.1
#   Sn = 0.9
#   Sp = 0.8
#   system.time({
#     x = vapply(1:50,function(i){
#       prob_obs_hhs(obs,prim,hh,q,Sn,Sp,M1=5,M2=5,T1=1,T2=1)
#     },0)
#   })
#   x_bf = prob_obs_hhs(obs,prim,hh,q,Sn,Sp)
#   exp(x_bf)
#   summary(exp(x))
#   hist(exp(x),breaks=50)
#   abline(v=exp(x_bf))
#   tpval = t.test(exp(x),mu=exp(x_bf))$p.value
#   tpval
#   expect_gt(tpval, 0.001)
# })
#
# test_that("pri prob is unbiased", {
#   expect_equal(.test_unbiased_pri(),1)
# })
#
# test_that("prob obs conditioned on pri is unbiased", {
#   expect_equal(.test_unbiased_obs_pri(),1)
# })
#
# test_that("secondary infection probabilities are unbiased", {
#   expect_equal(.test_unbiased_obs_pri_sim(),1)
# })
