#! /usr/local/R/2.15/bin/Rscript --vanilla

# Some tests for the mpmi package

# Got the idea for assert() from Julia language comparison with R
# https://github.com/JuliaLang/julia/blob/master/test/perf/micro/perf.R
#
assert = function(bool) {
    if (!isTRUE(bool)) stop('Assertion failed')
} 

# N.B., all.equal() uses tolerance .Machine$double.eps ^ 0.5 by default, 
# this may not be precise enough so it's overwritten below.
# (On my machine this tolerance is 1.490116e-08.)
#
all.equal <- function(...) 
{
    base:::all.equal(..., tolerance = .Machine$double.eps)
}
# All the tests for equality of the procedures should pass with 
# tolerance = .Machine$double.eps. The specific values are only
# stored by dput() to an accuracy of 1e-14.
#
all.equal.weak <- function(...) 
{
    base:::all.equal(..., tolerance = 1e-14)
}

library(mpmi)
data(mpmidata)

c1 <- cmi(cts)
d1 <- dmi(disc)
m1 <- mmi(cts, disc)

# Should be no missing values
assert(all(unlist(lapply(c1, function(x) all(!is.na(x))))))
assert(all(unlist(lapply(d1, function(x) all(!is.na(x))))))
assert(all(unlist(lapply(m1, function(x) all(!is.na(x))))))

# length of lists (bad code)
for (i in c("c1", "d1", "m1"))
{
    assert(length(eval(parse(text = i))) == 3)
}

# Dimensions
for (i in 1:3)
{
    assert(all.equal(dim(c1[[1]]), c(100L, 100L)))
    assert(all.equal(dim(d1[[1]]), c(75L, 75L)))
    assert(all.equal(dim(m1[[1]]), c(100L, 75L)))
}

# Test that the pairwise functions equal the implicitly parallel ones
# Diagonal:
c1pw1 <- cmi.pw(cts[, 1], cts[, 1])
d1pw1 <- dmi.pw(disc[, 1], disc[, 1])
m1pw1 <- mmi.pw(cts[, 1], disc[, 1])

# Off diagonal:
c1pw2 <- cmi.pw(cts[, 1], cts[, 2])
d1pw2 <- dmi.pw(disc[, 1], disc[, 2])
m1pw2 <- mmi.pw(cts[, 1], disc[, 2])

for (i in 1:3)
{
    assert(c1[[i]][1, 1] == c1pw1[[i]])
    assert(c1[[i]][1, 2] == c1pw2[[i]])

    assert(d1[[i]][1, 1] == d1pw1[[i]])
    assert(d1[[i]][1, 2] == d1pw2[[i]])

    assert(m1[[i]][1, 1] == m1pw1[[i]])
    assert(m1[[i]][1, 2] == m1pw2[[i]])
}

# Test that the no-jackknife versions get the same answers
c1n <- cminjk(cts)
d1n <- dminjk(disc)
m1n <- mminjk(cts, disc)

# Should be no missing values
assert(all(unlist(lapply(c1n, function(x) all(!is.na(x))))))
assert(all(unlist(lapply(d1n, function(x) all(!is.na(x))))))
assert(all(unlist(lapply(m1n, function(x) all(!is.na(x))))))

assert(all.equal(c1$mi, c1n))
assert(all.equal(d1$mi, d1n))
assert(all.equal(m1$mi, m1n))

# Pairwise no-jackknife equal to implicitly parallel version
# Diagonal
c1pw1n <- cminjk.pw(cts[, 100], cts[, 100])
d1pw1n <- dminjk.pw(disc[, 75], disc[, 75])
m1pw1n <- mminjk.pw(cts[, 100], disc[, 75])

# Off diagonal:
c1pw2n <- cminjk.pw(cts[, 99], cts[, 100])
d1pw2n <- dminjk.pw(disc[, 74], disc[, 75])
m1pw2n <- mminjk.pw(cts[, 99], disc[, 75])

assert(c1n[100, 100] == c1pw1n)
assert(c1n[99, 100] == c1pw2n)

assert(d1n[75, 75] == d1pw1n)
assert(d1n[74, 75] == d1pw2n)

assert(m1n[100, 75] == m1pw1n)
assert(m1n[99, 75] == m1pw2n)

# Check explicit paralellisation from the vignette 
# - Not portable
# - Makes some above tests redundant
# - For BCMI only
library(parallel)

# Mixed comparisons
hs <- apply(cts, 2, dpik, level = 3L, kernel = "epanech")
fi <- function(i)
{
    bcmis <- rep(NaN, 100)
    for (j in 1:100)
    {
        bcmis[j] <- mmi.pw(cts[,j], disc[,i], h = hs[j])$bcmi
    }
    return(bcmis)
}
parmmi <- mcmapply(fi, 1:75)

assert(all.equal(dim(parmmi), c(100L, 75L)))
assert(all.equal(m1$bcmi, parmmi))

# Continuous comparisons
hs2 <- apply(cts, 2, dpik, level = 3L, kernel = "epanech")
fi <- function(i)
{
    bcmis <- rep(NaN, 100)
    for (j in i:100)
    {
        bcmis[j] <- cmi.pw(cts[,i], cts[,j], h = hs2[c(i,j)])$bcmi
    }
    return(bcmis)
}

parcmi <- mcmapply(fi, 1:100)

assert(all.equal(dim(parcmi), c(100L, 100L)))
lt <- function(x) x[lower.tri(x, diag = TRUE)]
assert(all.equal(lt(c1$bcmi), lt(parcmi)))

# Discrete comparisons
fi <- function(i)
{
    bcmis <- rep(NaN, 75)
    for (j in i:75)
    {
        bcmis[j] <- dmi.pw(disc[,i], disc[,j])$bcmi
    }
    return(bcmis)
}

pardmi <- mcmapply(fi, 1:75)

assert(all.equal(dim(pardmi), c(75, 75L)))
lt <- function(x) x[lower.tri(x, diag = TRUE)]
assert(all.equal(lt(d1$bcmi), lt(pardmi)))

##################################################
# Check for some specific values

# Check using same dataset
cts2 <- 
structure(c(-0.0440282882314467, -0.669646072430327, -0.579628381117376, 
-0.90775085125724, 0.623125796191207, 0.382034585982775, -0.701974366236039, 
0.591820220745806, -0.279211766892763, 0.0871603688458017, 0.452041612791284, 
0.0105965395548507, -0.534345392411734, -0.172695440132611, -0.664825928667404, 
1.13256669577108, 0.990417318347305, -0.0279586343809654, -0.616254145262044, 
1.30194933160412, -1.70320152805629, -0.45162203955167, 2.31293070562161, 
0.503349078981106, -1.27900469979384, -0.7639311964385, -1.60937978957289, 
0.845827187579307, -0.522828471354191, -0.221120095161766, 1.00480218064388, 
1.09865177634465, 2.07279035633363, -1.72615781608062, 0.908906555169858, 
1.00165447371485, -0.438090965900502, -0.594426654913884, 1.88020912657578, 
0.0823343513612791, -0.308320823799523, -0.425178268228927, -0.743302611444874, 
0.520509405348235, -2.55501228067179, -0.206419487142881, 0.169094895070225, 
1.59388751937683, -0.617288172416216, -0.20305591440745, -0.0618448322549006, 
-0.242127299362674, -0.430127698279181, -0.974662750111274, 0.717782640854738, 
0.523615680259834, -0.824486350242906, 0.62522469629629, -0.22338124334855, 
0.246995232021837, 0.524009578138266, -0.142503405810888, -0.548809630577351, 
-0.189168863677349, -0.675749557414567, 1.11384917006589, 1.10178329153206, 
-0.104640348275532, -0.535181992286635, 1.15020310048664, -1.56312078784279, 
-0.297724013291145, 2.34107618703512, 0.525246896584568, -1.79942580058113, 
-0.69529425406869, -1.63166832869106, 0.592980509866606, -0.229288894288224, 
-0.201739037588088, 0.881045481001252, 0.738463622947662, 1.86590787918122, 
-1.80851978854533, 0.68963695729991, 1.32456687955228, -0.615687635935685, 
-0.552722441939986, 1.73307205954959, 0.350568591957014, -0.266447270865806, 
-0.712953770515152, -0.839854570959711, 0.943014706966472, -2.37997646058715, 
-0.360823429002887, 0.0878810182782329, 1.73226209948494, -0.613480644265789, 
-0.287775178749668), .Dim = c(50L, 2L), .Dimnames = list(NULL, 
    NULL))

disc2 <-
structure(c("H", "B", "H", "A", "B", "H", "A", "B", "A", "H", 
"H", "A", "H", "A", "B", "B", "A", "B", "A", "B", "A", "B", "B", 
"H", "A", "B", "H", "H", "H", "H", "B", "H", "A", "H", "B", "A", 
"A", "H", "B", "H", "H", "A", "A", "B", "B", "B", "B", "A", "A", 
"A", "B", "B", "A", "H", "B", "H", "A", "A", "H", "A", "B", "H", 
"H", "A", "B", "A", "A", "H", "H", "A", "H", "A", "A", "A", "H", 
"B", "B", "H", "A", "B", "B", "B", "H", "H", "B", "H", "A", "A", 
"A", "B", "H", "H", "H", "H", "B", "A", "B", "B", "B", "B"), .Dim = c(50L, 
2L))

assert(all.equal.weak(cts[, 99:100], cts2))
assert(all.equal.weak(disc[, 74:75], disc2))

# Calculated values:
cvold <- structure(list(mi = structure(c(1.27911130190712, 1.17834767387782,
                                1.17834767387782, 1.2794570792097), .Dim = c(2L,
                                2L)), bcmi = structure(c(1.33493636934022,
1.21512085444343, 1.21512085444343, 1.32510402507595), .Dim = c(2L, 2L)),
          zvalues = structure(c(10.8763971921941, 10.4618869120025,
                                10.4618869120025, 12.0155520566252), .Dim =
          c(2L, 2L))), .Names = c("mi", "bcmi", "zvalues")) 

dvold <- structure(list(mi = structure(c(1.09820954035319, 0.0914960536398212,
                                         0.0914960536398212, 1.09820954035319),
                                       .Dim = c(2L, 2L)), bcmi =
structure(c(1.11876141208849, 0.0446909702419578, 0.0446909702419578,
            1.11876141208849), .Dim = c(2L, 2L)), zvalues =
structure(c(268.437361450642, 0.69409221694838, 0.69409221694838,
            268.437361450642), .Dim = c(2L, 2L))), .Names = c("mi", "bcmi",
                        "zvalues"))

mvold <- structure(list(mi = structure(c(0.160514384030384, 0.246279973012949,
                                         0.0964862368145258,
                                         0.0885486467758484), .Dim = c(2L, 2L)),
                        bcmi = structure(c(0.0711043323121212,
                                           0.168974993559808,
                                           -0.0135059546908078,
                                           -0.036685028283886), .Dim = c(2L,
                                           2L)), zvalues =
                        structure(c(1.19099250708778, 2.83088622685523,
                                    -0.253782518304232, -0.646000734892397),
                                  .Dim = c(2L, 2L))), .Names = c("mi", "bcmi",
                        "zvalues"))

cvnew <- cmi(cts[,99:100])
dvnew <- dmi(disc[,74:75])
mvnew <- mmi(cts[,99:100], disc[,74:75])

for (i in 1:3)
{
    assert(all.equal.weak(cvold[[i]], cvnew[[i]]))
    assert(all.equal.weak(dvold[[i]], dvnew[[i]]))
    assert(all.equal.weak(mvold[[i]], mvnew[[i]]))
}
