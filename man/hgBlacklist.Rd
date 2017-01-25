\name{hgBlacklist}
\alias{hgBlacklist}
\alias{hg18bl}
\alias{hg19bl}
\alias{hg38bl}
\title{Blacklisted regions for human genome}
\description{
  Regions in the genome that give sequencing artifacts as listed in
  UCSC Genome browser track wgEncodeDacMapabilityConsensusExcludable
}
\usage{
  data(hgBlacklist)
}
\format{
  A list containing 25 matrices corresponding to chromosomes 1-22, X, Y &
  MT in that order. Each row is an interval that is to be discarded.
}
\details{
  The intervals for hg19 were obtained from UCSC Genome browser track
  \url{https://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability}
  wgEncodeDacMapabilityConsensusExcludable. The intervals were padded
  by 100 bases on both ends and an extra row (2^28, 2^28) added to each
  chromosome for programming convenience. The hg18 and hg38 data were
  obtained using liftover \url{https://genome.ucsc.edu/cgi-bin/hgLiftOver}.

  A fragment whose mid-point lies in an interval is discarded from analysis.
}
\references{
  \url{https://sites.google.com/site/anshulkundaje/projects/blacklists}
}
\keyword{datasets}
