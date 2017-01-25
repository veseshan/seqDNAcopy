\name{seqDNAcopy-internal}
\alias{bam2fragments}
\alias{bams2counts}
\alias{chromRange}
\alias{fragments2counts}
\alias{fragments2dataframe}
\alias{seqsegment}
\title{Internal seqDNAcopy functions}
\description{
  Internal functions of package seqDNAcopy.
}
\usage{
chromRange(i, chr="")
bam2fragments(bamFile, X=FALSE, mapq=20)
frags2counts(nfmid, tfmid, rstart, rend)
fragments2dataframe(nbam, tbam, iSizeLim=c(75,750), gbuild="hg19")
bams2counts(nBamFile, tBamFile, GCcorrect=TRUE, gbuild="hg19", mapq=20,
              iSizeLim=c(75,750), X=FALSE)
seqsegment(fcounts, sampleid="SeqSample", minBinCount=15, binSize=1000)
}
\details{These are not to be called directly by the user}
\keyword{internal}
