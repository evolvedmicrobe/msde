#--- basic test that old version produces same output as new version ------------

require(msde)
require(msde.devel)

# make models

hest.model <- msde::sde.make.model(list = msde::hestList)
hest.model2 <- msde.devel::sde.make.model(list = msde.devel::hestList)

identical(hest.model, hest.model2)

# simulate data

N <- 100
dT <- 1/252
theta <- c(alpha = .1, gamma = 5, beta = .8, sigma = .6, rho = -.7)
Y0 <- c(X = log(1e3), Z = 2*sqrt(.1))

set.seed(1)
hsim <- msde::sde.sim(model = hest.model, init.data = Y0,
                      params = theta, dt = dT, dt.sim = dT/100,
                      N = N, burn = 0, nreps = 5)

set.seed(1)
hsim2 <- msde.devel::sde.sim(model = hest.model2, init.data = Y0,
                             params = theta, dt = dT, dt.sim = dT/100,
                             N = N, burn = 0, nreps = 5)

identical(hsim, hsim2)


# inference

k <- 1
par.index <- 1
init <- msde::sde.init(data = hsim$data[1,,], dt = dT, k = k,
                       par.index = par.index, params = theta)
init2 <- msde.devel::sde.init(data = hsim2$data[1,,], dt = dT, k = k,
                              par.index = par.index, params = theta)


identical(init, init2)

nsamples <- 100
burn <- 10
rw.jump.sd <- rep(.1, 6)

set.seed(1)
hpost <- msde::sde.post(model = hest.model, init = init,
                        nsamples = nsamples, burn = burn,
                        data.out.ind = nsamples, rw.jump.sd = rw.jump.sd,
                        prior = NULL)

set.seed(1)
hpost2 <- msde.devel::sde.post(model = hest.model2, init = init2,
                               nsamples = nsamples, burn = burn,
                               data.out.ind = nsamples,
                               rw.jump.sd = rw.jump.sd,
                               prior = NULL)

identical(hpost, hpost2)
