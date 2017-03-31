#include <Rcpp.h>

#include "rcpp_binode.h"
#include "splitter.h"
#include "ijw_rand.h" 


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::export]]
SEXP simple_split(Rcpp::IntegerMatrix d, Rcpp::IntegerVector positions) {
    splitter *s = new splitter(d);                                     // define the splitter object s
    for (int i=0; i<positions.size(); i++) s->split(positions[i]-1);   // split at positions
    Rcpp::XPtr< splitter > pt(s, true);                                // get pointer as SEXP
    return pt;
}

// [[Rcpp::export]]
int nleaves(SEXP ptr) {
  Rcpp::XPtr< splitter > s(ptr);
  return s->nleaves();
}

// [[Rcpp::export]]
Rcpp::List node_labels(SEXP ptr) {
  Rcpp::XPtr< splitter > s(ptr);
  Rcpp::List L;
  NLRIterator<binode> ii(s->root());
  ii.nextLeaf();   // the root can never be a leaf
  while (!ii.isend()) {
    L.push_back(Rcpp::IntegerVector((*ii)->labels.begin(), (*ii)->labels.end())+1);  // sugar takes care
    ii.nextLeaf();                                                                   // of adding one
  } 
  return L;
}

// [[Rcpp::export]]
Rcpp::List get_phylo(SEXP ptr) {
  Rcpp::XPtr< splitter > s(ptr);

  std::vector<std::pair<int,int> > edges;            // set up date structure for edges
  std::vector<std::vector<int> > labels;             // and for labels  
 
  s->apesplit(edges, labels);
  Rcpp::NumericMatrix e(edges.size(), 2);
  
  for (size_t i=0; i< edges.size(); i++) {
    e(i, 0) = edges[i].first;
    e(i, 1) = edges[i].second;
  }
  
  Rcpp::List L =  Rcpp::List::create(
    Rcpp::Named("edge") = e,
    Rcpp::Named("Nnode") = s->nleaves()-1,
    Rcpp::Named("tip.label") = Rcpp::seq_len(s->nleaves())
  );
  L.attr("class") ="phylo";
  return L;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix case_control_leaves(SEXP ptr, Rcpp::IntegerVector cases) {
  Rcpp::XPtr< splitter > s(ptr);
  Rcpp::IntegerMatrix xxx(s->nleaves(), 2);
  s->getCaseControlLeaves(xxx, cases);
  return xxx;
}
/********************************************************************************************/
/** Get the count of haplotypeat at each leaf                                               */
/********************************************************************************************/
// [[Rcpp::export]]
Rcpp::IntegerVector  leaf_count(SEXP ptr) {
  Rcpp::XPtr< splitter > s(ptr);
  Rcpp::IntegerVector counts(s->nleaves());
  NLRIterator<binode> ii(s->root());
  ii.nextLeaf();   // the root can never be a leaf
  int index=0;
  while (!ii.isend()) {
    counts[index++] = (*ii)->labels.size();
    ii.nextLeaf();
  }
  return counts;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix  leaf_positions(SEXP ptr, double gap=1) {
  Rcpp::XPtr< splitter > s(ptr);
  int leaves = s->nleaves();
  Rcpp::NumericMatrix ypositions(leaves);
  NLRIterator<binode> ii(s->root());
  ii.nextLeaf();   // the root can never be a leaf
  ypositions(0, 0) = 0;
  ypositions(0, 1) = (*ii)->labels.size();
  ii.nextLeaf();
  int index=1;
  while (!ii.isend()) {
    ypositions(index, 0) = gap + ypositions(index-1, 1);
    ypositions(index, 1) = ypositions(index, 0) + (*ii)->labels.size();
    index++;
    ii.nextLeaf();
  }
  return ypositions;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix node_positions(SEXP ptr, double gap=1) {
  Rcpp::XPtr< splitter > s(ptr);
  int leaves = s->nleaves();
  Rcpp::NumericMatrix ypositions(leaves);
  NLRIterator<binode> ii(s->root());
  ii.nextLeaf();   // the root can never be a leaf
  ypositions(0, 0) = 0;
  ypositions(0, 1) = (*ii)->labels.size();
  ii.nextLeaf();
  int index=1;
  while (!ii.isend()) {
    ypositions(index, 0) = gap + ypositions(index-1, 1);
    ypositions(index, 1) = ypositions(index, 0) + (*ii)->labels.size();
    index++;
    ii.nextLeaf();
  }
  return ypositions;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix  stumps(SEXP ptr, int pos1, int pos2, double gap=1) {
  Rcpp::XPtr< splitter > s(ptr);
  int leaves = s->nleaves();
  Rcpp::NumericMatrix ypositions(leaves);
  NLRIterator<binode> ii(s->root());
  ii.nextLeaf();   // the root can never be a leaf
  ypositions(0, 0) = 0;
  ypositions(0, 1) = (*ii)->labels.size();
  ii.nextLeaf();
  int index=1;
  while (!ii.isend()) {
    ypositions(index, 0) = gap + ypositions(index-1, 1);
    ypositions(index, 1) = ypositions(index, 0) + (*ii)->labels.size();
    index++;
    ii.nextLeaf();
  }
  return ypositions;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix  get_coordinates(SEXP ptr, double gap=1) {
  Rcpp::XPtr< splitter > s(ptr);
  s->calculate_top_bottom(gap);   // gets the tops and bottoms 
                                  // for the leaves and the rest of the tree
  int leaves = s->nleaves();
  Rcpp::NumericMatrix coords(2*leaves+3*(leaves-1), 3);
  
  s->get_coordinates(coords);
  for (int i=0; i<coords.nrow(); i++) {
    coords(i, 0) +=1; 
  }
  
  return coords;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_splitTestCC(Rcpp::IntegerMatrix data, 
                 Rcpp::IntegerVector cases,
                 Rcpp::IntegerVector positions,
                 int maxk,
                 int reps,
                 const std::string pickStat)
{
     splitter s(data);                           // define the splitter object s
     for (int i=0; i<positions.size(); i++) s.split(positions[i]-1);   // split at positions
     int ncases=cases.size();
     int samplesize = data.nrow();
  
     std::vector<int> caseVec(cases.begin(), cases.end());
     Rcpp::NumericMatrix randteststats(reps+1, maxk);
    
     std::vector<double> stat=s.getStat(caseVec, maxk, pickStat.c_str());   // use default test statistic
     for (int jj=0;jj < maxk;jj++) randteststats(0, jj) = stat[jj];
  
     rng r;
     std::vector<int> cc(samplesize);
     for (int i=0;i< samplesize;i++) cc[i]=i;
    
     for (int i=0;i< reps;i++) {
        permute(cc, r);
        caseVec.assign(cc.begin(),cc.begin()+ ncases);
        std::sort(caseVec.begin(),caseVec.begin()+ ncases);
        std::vector<double> rstat =  s.getStat(caseVec, maxk, pickStat.c_str());
        for (int jj=0; jj< maxk; jj++) {
            randteststats(i+1, jj) = rstat[jj];
      }
    }
    return randteststats;
}



// [[Rcpp::export]]
Rcpp::NumericMatrix rcppsplittestqtrait( Rcpp::IntegerMatrix data, 
                                      Rcpp::NumericVector qtl,
                                      Rcpp::IntegerVector positions, 
                                      int reps, 
                                      int maxk,
                                      const std::string statPick)
{
  splitter s(data);
  for (int i=0; i< positions.size(); i++) s.split(positions[i]-1);
  
  std::vector<double> myqtl(qtl.begin(), qtl.end());
  Rcpp::NumericMatrix randteststats(reps+1, maxk);

  std::vector<double> stat=s.qtlStat(myqtl, maxk, statPick.c_str());
  for (int jj=0; jj< maxk; jj++) randteststats(0, jj) = stat[jj];
  
  rng r;
  for (int i=0; i< reps; i++) {
    permute(myqtl, r);
    std::vector<double> rstat =  s.qtlStat(myqtl, maxk, statPick.c_str());
    for (int jj=0;jj< maxk;jj++) {
      randteststats(i+1, jj) = rstat[jj];
    }
  }
  
  return randteststats;
}

/*******************************************************************************************************/
/** Get a simple tree from a set of data doing all splitting within c++ and return an object of class 
* phylo 
*/
/******************************************************************************************************/
// [[Rcpp::export]]
Rcpp::List rcpp_split_simple( Rcpp::IntegerMatrix data, 
                   Rcpp::IntegerVector positions)  { 

  splitter s(data);                                                // define the splitter object s
  for (int i=0;i< positions.size();i++) s.split(positions[i]-1);   // split at positions

  std::vector<std::pair<int,int> > edges;                          // set up date structure for edges
  std::vector<std::vector<int> > labels;                           // and for labels  
  
  s.apesplit(edges, labels);
  
  // get the lengths of the labels to allow us to pass this information back to R
  // the information is returned in the correct order for ape which is root->left->right
  
  Rcpp::IntegerVector leafcount(s.nleaves());
  Rcpp::List comblabels;
  
  for (size_t ii=0; ii < labels.size(); ii++) {
    comblabels.push_back(Rcpp::IntegerVector(labels[ii].begin(), labels[ii].end())+1);
    leafcount[ii] = static_cast<int>(labels[ii].size());
  }
  
  int nedges=static_cast<int>(edges.size());
  if (nedges != 2*(s.nleaves()-1))
    throw std::range_error("problem in C++ code split_simple\n");
  
  Rcpp::IntegerMatrix edge(nedges, 2);
  for (int i=0; i<nedges; i++) {
    edge(i, 0) = edges[i].first;
    edge(i, 1) = edges[i].second;
  }
  // now try to get the lengths.  Note that the centre of this split is positions[0].
  std::vector<int> nodePos;
  s.getNodesPositions(nodePos);
  Rcpp::IntegerVector nodepos(nodePos.begin(), nodePos.end());
  for (size_t ii=0;ii<nodePos.size();ii++) nodepos[ii]=nodePos[ii];
  
  Rcpp::List L = Rcpp::List::create(
      Rcpp::Named("edge") = edge,
      Rcpp::Named("tip.label") = Rcpp::seq_len(s.nleaves()),
      Rcpp::Named("Nnode") = s.nleaves()-1,
      Rcpp::Named("labels") = comblabels,
      Rcpp::Named("nodepos") = nodepos,
      Rcpp::Named("leafcount") = leafcount
  );
  L.attr("class") ="phylo";
  return L;
}


/** Get the split in a form that is suitable for using within ape                    
* The nodepos and edgepos give the left hand and right hand 
* positions of the minimum and maximum range of the haplotypes at each 
* node and tip (note both are -1 for tips with only a single haplotype 
* i.e. those that have been successfully separated)                        
* The positions are calculated for all data, not just those positions that have
* been separated out, but they always include the starting position (which should guarantee
* that all positions that the haplotypes have been split on are included).   
*/
// [[Rcpp::export]]
Rcpp::NumericVector qtrait_statistics(SEXP ptr, Rcpp::NumericVector qtrait, const std::string statPick) {
  Rcpp::XPtr< splitter > s(ptr);
  int len=s->nleaves();
  
  std::vector<double> qv(qtrait.begin(), qtrait.end());
  double xbar = mean(qtrait);
  double sample_s = sd(qtrait);
  double s_squared=sample_s*sample_s;
  Rcpp::NumericVector sumx(len, 0.0);
  Rcpp::NumericVector n(len, 0.0);
  
  NLRIterator<binode> ii(s->root());
  ii.nextLeaf();   // the root can never be a leaf
  int index=0;
  while (!ii.isend()) {
    for (size_t jj=0; jj<(*ii)->labels.size(); jj++) 
      sumx[index] += qtrait[(*ii)->labels[jj]];
    n[index++] = (*ii)->labels.size();
    ii.nextLeaf();
  }
  switch(statPick.c_str()[0]) {
      case 'Z': return  n*(sumx/n - xbar)*(sumx/n - xbar)/s_squared;
      case 'A': return abs(sqrt(n)*(sumx/n - xbar)/sample_s);
      case 'P': return sqrt(n)*(sumx/n - xbar)/sample_s;
      case 'N': return -sqrt(n)*(sumx/n - xbar)/sample_s;
  default: 
    Rcpp::stop("pickstat must be one of 'Z', 'A', 'N' or 'P'");
  }
}
  
