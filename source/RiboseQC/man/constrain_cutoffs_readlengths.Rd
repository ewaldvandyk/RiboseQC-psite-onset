% Please edit documentation in R/riboseqc.R
\name{constrain_cutoffs_readlengths}
\alias{constrain_cutoffs_readlengths}
\title{Choose cutoffs jointly across read lengths}
\usage{
constrain_cutoffs_readlengths(cutoffs_list, max_diff = 1)
}
\arguments{
\item{cutoffs_list}{list with one element per read length (names), each the output of 
\code{calc_cutoffs_from_profiles} with \code{method="edge"}. An element named "all" is left as it is.}

\item{max_diff}{Maximum difference in cutoff (nt) between two read lengths, per nt of difference in length. Defaults to 1}
}
\value{
\code{cutoffs_list}, where for each read length \code{final_cutoff} is the cutoff from the chain. 
The cutoff chosen for the read length on its own is kept in \code{final_cutoff_independent}; 
\code{cutoff_constraint} is "chain" or "imputed" (left out of the chain, cutoff taken from the neighbours); 
\code{cutoff_chain_margin} is the difference in score between the best chain and the best chain with another 
cutoff for the read length with the highest weight, as a fraction of the best score (0: undecided, 1: no alternative).
}
\description{
The distance between the 5'end of a footprint and the P-site changes slowly with read length: 
read lengths that differ by one nucleotide cannot have cutoffs that differ by more than about one nucleotide.
This function uses that to choose the cutoffs of all read lengths together, 
instead of each read length on its own.
}
\details{
Each read length has a set of possible cutoffs, 3 nt apart (the in-frame positions), each with a
\code{support} value calculated by \code{calc_cutoffs_from_profiles} (\code{method="edge"}). 
A chain is one cutoff per read length such that two read lengths never differ by more than 
\code{max_diff} nt per nt of difference in length. The chain with the highest sum of weighted support is chosen. 
The weight of a read length is \code{log10(1+n)}, with n its number of 5'ends around start codons, multiplied by 
its frame preference rescaled so that no preference (1/3) gives 0 and a perfect one gives 1: 
read lengths with many reads and a clear frame count the most, read lengths without frame preference do not count. \cr \cr
Because the possible cutoffs of one read length are 3 nt apart, a neighbouring read length leaves only one of them possible 
(with \code{max_diff=1}): the support of all read lengths is thus pooled into one decision. \cr \cr
A read length that cannot be part of the best chain (e.g. a sparse read length for which the wrong frame was found) 
is left out of it, and then gets the possible cutoff closest to that of its neighbours in the chain. 
On ties, the larger cutoff is chosen, as in \code{calc_cutoffs_from_profiles}.
}
\seealso{
\code{\link{calc_cutoffs_from_profiles}}, \code{\link{RiboseQC_analysis}}
}
\keyword{Ribo-seQC}
