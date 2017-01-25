#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List frags2counts(NumericVector nfmid, NumericVector tfmid, NumericVector rstart, NumericVector rend) {

  // normal fragments that are discarded
  int n = nfmid.size();
  NumericVector ndiscard(n);

  int i = 0;                           // first blacklisted region
  for(int j=0; j < n; ++j) {           // loop through fragments
    if (nfmid[j] < rstart[i]) {        // if fragmid before region start
      ndiscard[j] = 0;                 // keep fragment
    } else {                           // else check fragmid
      if (nfmid[j] <= rend[i]) {       // if fragmid within region
	ndiscard[j] = 1;               // discard fragment
      } else {                         // find the next region at or past it
	while (nfmid[j] > rend[i]) {   // while fragmid exceeds region end
	  i = i+1;                     // go to next region
	}                              // region change complete
	if (nfmid[j] < rstart[i]) {    // if fragmid before region start
	  ndiscard[j] = 0;             // keep fragment
	} else {                       // else discard it
	  ndiscard[j] = 1;             // since we know nfmid[j] <= rend[i]
	}
      }
    }
  }

  // tumor fragments that are discarded
  n = tfmid.size();
  NumericVector tdiscard(n);

  i = 0;                               // first blacklisted region
  for(int j=0; j < n; ++j) {           // loop through fragments
    if (tfmid[j] < rstart[i]) {        // if fragmid before region start
      tdiscard[j] = 0;                 // keep fragment
    } else {                           // else check fragmid
      if (tfmid[j] <= rend[i]) {       // if fragmid within region
	tdiscard[j] = 1;               // discard fragment
      } else {                         // find the next region at or past it
	while (tfmid[j] > rend[i]) {   // while fragmid exceeds region end
	  i = i+1;                     // go to next region
	}                              // region change complete
	if (tfmid[j] < rstart[i]) {    // if fragmid before region start
	  tdiscard[j] = 0;             // keep fragment
	} else {                       // else discard it
	  tdiscard[j] = 1;             // since we know tfmid[j] <= rend[i]
	}
      }
    }
  }

  // chromosome < 2^28 bases; so at most 2^28/100 = 2684355 bins
  // use the valid fragments to create a count matrix

  int ll;
  NumericVector ncount(2684355);
  NumericVector tcount(2684355);
  NumericVector keep(2684355);      // bins with non-zero counts

  for(int j=0; j < 2684355; ++j) {
    ncount[j] = 0;                  // initialize normal count to 0
    tcount[j] = 0;                  // initialize tumor count to 0
    keep[j] = 0;                    // initialize keep to 0
  }

  n = nfmid.size();
  for(int j=0; j < n; ++j) {        // loop through normal fragments
    if (ndiscard[j] == 0) {
      ll = floor(nfmid[j]/100);
      ncount[ll] += 1;
      keep[ll] = 1;
    }
  }

  n = tfmid.size();
  for(int j=0; j < n; ++j) {        // loop through tumor fragments
    if (tdiscard[j] == 0) {
      ll = floor(tfmid[j]/100);
      tcount[ll] += 1;
      keep[ll] = 1;
    }
  }

  int nnzbins = 0;                      // non-zero bin counter
  for(int j=0; j < 2684355; ++j) {
    nnzbins += keep[j];
  }

  NumericVector pos(nnzbins);
  NumericVector normal(nnzbins);
  NumericVector tumor(nnzbins);
  i = 0;
  for(int j=0; j < 2684355; ++j) {
    if (keep[j] == 1) {
      pos[i] = 50 + 100*j;
      normal[i] = ncount[j];
      tumor[i] = tcount[j];
      i += 1;
    }
  }

  return List::create(_["pos"] = pos, _["normal"] = normal, _["tumor"] = tumor);
}
