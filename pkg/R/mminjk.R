mminjk <-
function(cts, disc, level = 3L, na.rm = FALSE, h, ...)
{
    cts <- as.matrix(cts)
    mode(cts) <- "double"

    mi <- matrix(0.0, nrow = ncol(cts), ncol = ncol(disc))
    # bcmi <- matrix(0.0, nrow = ncol(cts), ncol = ncol(disc))
    # zans <- matrix(0.0, nrow = ncol(cts), ncol = ncol(disc))

    # Convert to integers (probably slow)
    dint <- matrix(0L, nrow = nrow(disc), ncol = ncol(disc))
    for (i in 1: ncol(disc))
    {
        dint[,i] <- as.integer(factor(disc[,i]))
    }

    # Calculate bandwidths
    if (missing(h))
    {
        if(na.rm)
        {
            h2 <- function(x)
            {
                return(dpik(x[is.finite(x)], level = level, kernel = "epanech", ...))
            }
            h <- apply(cts, 2, h2)
        } else 
        {
            h <- apply(cts, 2, dpik, level = level, kernel = "epanech", ...)
        }
    }
 
    result <- .Fortran("mmimnjk", 
                    cts, 
                    as.integer(nrow(cts)), 
                    as.integer(ncol(cts)),
                    dint,
                    as.integer(nrow(dint)), 
                    as.integer(ncol(dint)),
                    mi = mi, 
                    # bcmi = bcmi, 
                    # zvalues = zans,
                    as.double(h), 
                    NAOK = TRUE, 
                    DUP = TRUE)

    return(result$mi)
}

mminjk.pw <- function(cts, disc, h = dpik(cts[!is.na(cts)], level = 3L, kernel = "epanech"))
{
    if (length(cts) != length(disc)) stop("Input vectors must be the same length")

    # Remove missing values pairwise
    ok <- !is.na(disc) & !is.na(cts)
    disc <- disc[ok]
    cts <- cts[ok]

    lf <- length(cts)

    mi <- 0.0

    result <- .Fortran("mmipwnjk", cts = as.double(cts), 
                        lc = as.integer(lf), 
                        disc = as.integer(factor(disc)),
                        h = as.double(h),
                        mi = as.double(mi), 
                        # bcmi = as.double(0),
                        # zvalue = as.double(0),
                        DUP = TRUE)

    return(result$mi)
}
