#--- basic test that old version produces same output as new version ------------

require(msde) # stable version
require(msde2) # development version
pkg1 <- function(fun) {
  eval(parse(text = paste0("msde::", fun)))
}
pkg2 <- function(fun) {
  eval(parse(text = paste0("msde2::", fun)))
}
SEED <- 906

# 1. build with clang: failed (not supported by rtools?)
# 2. against -O3 -ffast-math: pass
# 3. Switch to Box-Muller: failed (slower than inversion)
# 4. Collapse lmvn double for-loop: pass (5-10

# make models
hmod <- pkg1("sde.make.model")(list = msde::hestList,
                               showOutput = TRUE, rebuild = TRUE)
hmod2 <- pkg2("sde.make.model")(list = msde2::hestList,
                                showOutput = TRUE, rebuild = TRUE)

identical(hmod, hmod2) # c++ codes now different

# simulate data

N <- 100
dT <- 1/252
theta <- c(alpha = .1, gamma = 5, beta = .8, sigma = .6, rho = -.7)
Y0 <- c(X = log(1e3), Z = 2*sqrt(.1))

set.seed(SEED)
hsim <- pkg1("sde.sim")(model = hmod, init.data = Y0,
                      params = theta, dt = dT, dt.sim = dT/100,
                      N = N, burn = 0, nreps = 5)

set.seed(SEED)
hsim2 <- pkg2("sde.sim")(model = hmod2, init.data = Y0,
                         params = theta, dt = dT, dt.sim = dT/100,
                         N = N, burn = 0, nreps = 5)

identical(hsim, hsim2)


# inference

k <- 1
par.index <- 1
init <- pkg1("sde.init")(data = hsim$data[1,,], dt = dT, k = k,
                         par.index = par.index, params = theta)
init2 <- pkg2("sde.init")(data = hsim2$data[1,,], dt = dT, k = k,
                          par.index = par.index, params = theta)


identical(init, init2)

nsamples <- 1e5
burn <- 10
rw.jump.sd <- rep(.1, 6)

set.seed(SEED)
tm.hp1 <- system.time({
  hpost <- pkg1("sde.post")(model = hmod, init = init,
                            nsamples = nsamples, burn = burn,
                            data.out.ind = 100,
                            rw.jump.sd = rw.jump.sd,
                            prior = NULL)
})

set.seed(SEED)
tm.hp2 <- system.time({
  hpost2 <- pkg2("sde.post")(model = hmod2, init = init2,
                             nsamples = nsamples, burn = burn,
                             data.out.ind = 100,
                             rw.jump.sd = rw.jump.sd,
                             prior = NULL)
})

identical(hpost, hpost2)

tm.hp2[3]/tm.hp1[3] # 15% faster

#--- just look at RNG speed -----------------------------------------------------

RNGkind(normal.kind="Inversion")
system.time({
  rnorm(1e7)
})

# slower for me...
RNGkind(normal.kind="Box-Muller")
system.time({
  rnorm(1e7)
})
