logit = function(x) log(x/(1-x))
ilogit = function(x) exp(x)/(1+exp(x))

# Questions about fitting model to data
# - Does infection rate diminish with household size? By how much?
# - Does the presence of young children increase infection rates in the household?
# - Does the infection rate change with average age in the household?
# - Does the infection rate change with blood suguar (???)
# - Can I affirm evidence of camel contact risk for primary infection?
# - How much risk is not accounted for by modelling? (ie like community exposure)
# - Is error clusterred by household my models? (Implies unaccounted for clustering. How do I measure this?)

# hh: household identifier
# primary: boolean, whether a person is a primary infection
# sar: secondary attack rate
# Returns: simulated final infections
simulate.sar = function(primary, sar){
  n = length(primary)
  if(sum(primary)==0){
    rep(0,n)
  }else{
    M = randomgraph(n,sar)
    as.vector(clusters(list(M),primary))*1
  }
}

multi.fact = function(n,levels,probs=NULL){
  ln = length(levels)
  if(is.null(probs)) probs = rep(1/ln,ln)
  if(length(probs) == 1) probs = rep(probs[[1]],ln)
  x = sample(1:length(probs),n,replace=TRUE,probs)
  factor(levels[x],levels=levels)
}

gen.member.mat = function(hhs,
                          sex.p1=0,
                          age.min=5,age.max=100,age.p1=1,age.p2=3,
                          camel.activity.beta=c(logit(0.2),0,0)
                          ){
  do.call(rbind, lapply(unique(hhs),function(hh_id){
    n = sum(hhs==hh_id)
    sex = multi.fact(n,c("male","female"),c(1-ilogit(sex.p1),ilogit(sex.p1)))
    age = floor(rbeta(n,age.p1,age.p2)*(age.max-age.min)+age.min)
    p.camel.activity = ilogit(model.matrix(~1+sex+age,data.frame("sex"=sex,"age"=age)) %*% camel.activity.beta)
    camel.activity = factor(ifelse(rbinom(n,1,p.camel.activity)==1,"yes","no"),levels=c("no","yes"))
    data.frame(
      hh=hh_id,
      sex=sex,
      age=age,
      camel.activity=camel.activity
    )
  }))
}

# Three simulation methods to generate a dataset
# s1: primary infections follow binom(par)
gen.sim.data = function(member.df,beta.par=c(0,0,0,0),sar=0.02,icc=0,Sn=1,Sp=1){
  # get standard deviation accross households using icc
  # icc = sig0^2 / (sig0^2 + pi^2/3)
  sig0 = pi * sqrt(icc) / sqrt(3-3*icc)
  sig.dev.hh = rnorm(length(unique(member.df$hh)),0,sig0)
  # match(member.df$hh,unique(member.df$hh))

  y.df = do.call(rbind, lapply(unique(member.df$hh),function(hh_id){
    hh.df = member.df[member.df$hh == hh_id,]
    hh.matrix = model.matrix(~1+sex+age+camel.activity,data=hh.df[,c("sex","age","camel.activity")])
    hh.dev = sig.dev.hh[which(unique(member.df$hh)==hh_id)]
    lpprim = hh.matrix %*% beta.par + hh.dev
    # Simulate primary infections
    y.prim = rbinom(nrow(hh.df),1,ilogit(lpprim))
    # Simulate secondary infections
    y = simulate.sar(y.prim,sar)
    y.obs = rbinom(length(y),1,ifelse(
      y==0,1-Sp,Sn
    ))
    data.frame(y.prim=y.prim,y.act=y,y=y.obs)
  }))
  cbind(y.df,member.df)
}

sim1 = function(){
  set.seed(999)
  num.hh = 200    # number of households
  max.hh = 25     # max household size
  hh.rgamma = rgamma(num.hh,2,10)
  hh.sizes = ceiling(hh.rgamma*max.hh/max(hh.rgamma))

  hhs = do.call(c,lapply(1:length(hh.sizes),function(i){
    rep(i,hh.sizes[i])
  }))
  member.df = gen.member.mat(
    hhs,
    sex.p1=logit(0), # all female
    age.min = 25, age.max = 25, # all 25
    camel.activity.beta = c(logit(0),0,0) # no camel activity
  )
  member.df = gen.sim.data(member.df,beta.par = c(logit(0.1),0,0,0),sar = 0.02,icc = 0,Sn = 1, Sp=1)
  member.df
}

