dmi.pw <- function(disc1, disc2)
{
    disc1 <- as.integer(factor(disc1))
    disc2 <- as.integer(factor(disc2))
    ans <- as.double(0)
    return(.Fortran("dmi",
                    disc1,
                    as.integer(length(disc1)),
                    disc2,
                    as.integer(length(disc2)),
                    result = ans, DUP = FALSE)$result)
}

dmi <- function(dmat)
{
    # Convert to integers
    dint <- matrix(0L, nrow = nrow(dmat), ncol = ncol(dmat))
    for (i in 1: ncol(dmat))
    {
        dint[,i] <- as.integer(factor(dmat[,i]))
    }

    ans <- matrix(0.0, nrow = ncol(dint), ncol = ncol(dint))
    return(.Fortran("dmim",
                    dint,
                    as.integer(nrow(dint)),
                    as.integer(ncol(dint)),
                    result = ans, NAOK = TRUE, DUP = FALSE)$result)
}

