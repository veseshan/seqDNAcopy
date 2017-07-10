seqsegment <-
function(fcounts, sampleid="SeqSample", minBinCount=15, binSize=1000, ...) {
    # collapse the 100 base bins according to binsize
    fcounts$bins <- fcounts$chrom + floor(fcounts$pos/binSize)*binSize/2^28
    # data corresponding to the new bins
    zz0 <- list()
    zz0$chrom <- tapply(fcounts$chrom, fcounts$bins, function(x){x[1]})
    zz0$pos <- tapply(fcounts$pos, fcounts$bins, function(x){x[1]})
    zz0$pos <- binSize*(floor(zz0$pos/binSize) + 0.5)/1e6
    zz0$normal <- tapply(fcounts$normal, fcounts$bins, sum)
    zz0$tumor <- tapply(fcounts$tumor, fcounts$bins, sum)
    zz0 <- as.data.frame(zz0)
    # counts to segments
    zz <- zz0[zz0[,3]>minBinCount,]
    # scale the read counts
    tscl <- sum(zz[,"normal"])/sum(zz[,"tumor"])
    zchr <- zz[,"chrom"]
    zpos <- zz[,"pos"]
    # log-ratio
    zlr <- log2(zz[,"tumor"]*tscl+1) - log2(zz[,"normal"]+1)
    # weights
    zwts <- log2(zz[,"normal"]+1-minBinCount)

    zcna <- CNA(as.matrix(zlr), zchr, zpos, sampleid=sampleid)
    zout <- segment(zcna, weights=zwts, ...)
    zout
}
