chromRange <- function(i, chr="") {
  eval(parse(text=paste('RangesList("', chr, i,'"=IRanges(start=0, end=268435456L))',sep="")))
}

HUMAN.BUILDS=c("hg19")
MOUSE.BUILDS=c("mm10")

autosomeMaxValue <- function(gbuild) {
    if(gbuild %in% HUMAN.BUILDS) {
        return(22)
    } else if(gbuild %in% MOUSE.BUILDS) {
        return(19)
    } else {
        cat("\n\n    seqDNAcopy:bam2counts:autosomesMaxValue:: FATAL ERROR\n")
        cat("    Unknown genome build: ",paste0("[",gbuild,"]"),"\n\n\n")
        stop("FATAL-ERROR[seqDNAcopy:bam2counts:autosomesMaxValue(16)]")
    }
}

bam2fragments <- function(bamFile, X=FALSE, mapq=20, gbuild="hg19") {
    # get position, mate position and insert size (fragment length)
    what <- c("pos","mpos","isize")
    # get first mate from properly paired reads which pass QC
    flag=scanBamFlag(isNotPassingQualityControls=FALSE, isPaired=TRUE, isFirstMateRead=TRUE, hasUnmappedMate=FALSE, isDuplicate=FALSE, isSecondaryAlignment=FALSE)
    bam <- list()
    # autosomes

    maxAutosomeNumber=autosomeMaxValue(gbuild)

    for(i in 1:maxAutosomeNumber) {
        which <- chromRange(i)
        param <- ScanBamParam(flag=flag, which = which, what = what, mapqFilter=mapq)
        # scan the data
        bam[[i]] <- scanBam(bamFile, param=param)[[1]]
    }
    if (X) {
        # chromosome X
        which <- chromRange("X")
        param <- ScanBamParam(flag=flag, which = which, what = what)
        # scan the data
        bam[[maxAutosomeNumber+1]] <- scanBam(bamFile, param=param)[[1]]
    }
    bam
}

fragments2dataframe <- function(nbam, tbam, iSizeLim=c(75,750), gbuild="hg19") {
    # normal fragment midpoints
    nfmid <- lapply(nbam, function(x, ll, ul) {
        x$isize <- abs(x$isize)
        sort((pmin(x$pos, x$mpos) + x$isize/2)[x$isize >= ll & x$isize <= ul])
    }, iSizeLim[1], iSizeLim[2])
    # tumor fragment midpoints
    tfmid <- lapply(tbam, function(x, ll, ul) {
        x$isize <- abs(x$isize)
        sort((pmin(x$pos, x$mpos) + x$isize/2)[x$isize >= ll & x$isize <= ul])
    }, iSizeLim[1], iSizeLim[2])
    # bin the mid points into interval of size 100
    binnedcounts <- list()
    nchr <- length(nbam)
    # load blacklisted regions data
    data(hgBlacklist, package="seqDNAcopy",envir=environment())
    blgbuild <- paste(gbuild, "bl", sep="")
    gbl <- get(blgbuild, envir=environment())
    for(i in 1:nchr) binnedcounts[[i]] <- frags2counts(nfmid[[i]], tfmid[[i]], gbl[[i]][,1], gbl[[i]][,2])
    # total number of bins with nonzero count
    nbins <- unlist(lapply(binnedcounts, function(x) {length(x$pos)}))
    # convert to matrix
    fcounts <- matrix(0, sum(nbins), 4)
    colnames(fcounts) <- c("chrom","pos","normal","tumor")
    fcounts[,1] <- rep(1:nchr, nbins)
    for(i in 1:nchr) {
        ii <- fcounts[,1]==i
        fcounts[ii,2] <- binnedcounts[[i]]$pos
        fcounts[ii,3] <- binnedcounts[[i]]$normal
        fcounts[ii,4] <- binnedcounts[[i]]$tumor
    }
    # return data
    as.data.frame(fcounts)
}

bams2counts <- function(nBamFile, tBamFile, GCcorrect=TRUE, gbuild="hg19", mapq=20, iSizeLim=c(75,750), X=FALSE) {
    # normal bam read data ("pos","mpos","isize")
    nbam <- bam2fragments(nBamFile, X, mapq, gbuild)
    # tumor bam read data ("pos","mpos","isize")
    tbam <- bam2fragments(tBamFile, X, mapq, gbuild)
    # get the fragment count in bins centered every 100 bases
    out <- fragments2dataframe(nbam, tbam, iSizeLim, gbuild)
    if (GCcorrect) {
        # gc percentage
        gcpct <- rep(NA_real_, nrow(out))
        # get GC percentages from pctGCdata package
        # loop thru chromosomes
        nchr <- max(out$chrom) # IMPACT doesn't have X so only 22
        for (i in 1:nchr) {
            ii <- which(out$chrom==i)
            # allow for chromosomes with no SNPs i.e. not targeted
            if (length(ii) > 0) {
                gcpct[ii] <- getGCpct(i, out$pos[ii], gbuild)
            }
        }
        # GC correction
        tscl <- sum(out$normal)/sum(out$tumor)
        # normal and tumor fragment counts by GC percent level
        ncount <- tapply(out$normal, gcpct, sum)
        tcount <- tapply(out$tumor, gcpct, sum)
        # log-ratio of tumor to normal counts as function of GC percentages
        # standardized by the ratio of library sizes
        logtn <- log2(tcount*tscl) - log2(ncount)
        # percent GC
        pctgc <- as.numeric(names(logtn))
        # weights for each GC percent level
        wts <- log2(tcount+ncount)
        # bins that have non-zero count for both tumor and normal
        ii <- which(is.finite(logtn))
        # gc bias estimated using loess
        gcb <- loess(logtn[ii] ~ pctgc[ii], weights=wts[ii], span=0.25)
        # GC corrected library scaled tumor counts
        jj <- match(gcpct, gcb$x)
        gcscl <- 2^{-gcb$fitted[jj]}
        # gcscl is NA if a bin has NA for gcpct; make it 1
        gcscl[is.na(gcscl)] <- 1
        # scale tumor counts by gcscl
        out$tumor <- gcscl*out$tumor
    }
    out
}
