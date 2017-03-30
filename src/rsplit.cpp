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
SEXP  simple_split(Rcpp::IntegerMatrix d, Rcpp::IntegerVector positions) {
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
Rcpp::List get_phylo(SEXP ptr) {
  std::vector<std::pair<int,int> > edges;            // set up date structure for edges
  std::vector<std::vector<int> > labels;             // and for labels  
  s.apesplit(edges, labels);
  NumericMatrix e(s->nleaves(), 2);
  
  
  
  
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix  case_control_leaves(SEXP ptr, Rcpp::IntegerVector cases) {
  Rcpp::XPtr< splitter > s(ptr);
  Rcpp::IntegerMatrix xxx(s->nleaves(), 2);
  s->getCaseControlLeaves(xxx, cases);
  return xxx;
}

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
    for (int jj=0;jj < maxk;jj++) randteststats(0,jj) = stat[jj];
  
    rng r;
    std::vector<int> cc(samplesize);
    for (int i=0;i< samplesize;i++) cc[i]=i;
    
    
    
    for (int i=0;i< reps;i++) {
        permute(cc,r);
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
Rcpp::NumericMatrix rcppsplittestQTL( Rcpp::IntegerMatrix data, 
                                      Rcpp::IntegerVector positions, 
                                      Rcpp::NumericVector qtl,
                                      int reps, 
                                      int maxk,
                                      const std::string statPick)
{
  Rcpp::NumericMatrix randteststats(reps+1, maxk);
  
  splitter s(data);
  for (int i=0; i< positions.size(); i++) s.split(positions[i]);
  
  std::vector<double> myqtl(qtl.begin(), qtl.end());
  std::vector<double> stat=s.qtlStat(myqtl, maxk, statPick.c_str());
  for (int jj=0;jj< maxk;jj++) randteststats(0, jj) = stat[jj];
  
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


// [[Rcpp::export]]
Rcpp::List rcpp_split_simple( Rcpp::IntegerMatrix data, 
                   Rcpp::IntegerVector positions)  { 

  splitter s(data);                                  // define the splitter object s
  for (int i=0;i< positions.size();i++) s.split(positions[i]);   // split at positions
  int len=s.nleaves();                               // how many leaves do we have on the tree
  
  std::vector<std::pair<int,int> > edges;            // set up date structure for edges
  std::vector<std::vector<int> > labels;             // and for labels  
  
  s.apesplit(edges, labels);
  
  // get the lengths of the labels to allow us to pass this information back to R
  // the information is returned in the correct order for ape which is root->left->right
  
  Rcpp::IntegerVector leafcount(len);
  Rcpp::IntegerVector comblabels(data.nrow());
  
  int count=0;
  for (size_t ii=0; ii < labels.size(); ii++) {
    for (size_t jj=0; jj < labels[ii].size(); jj++) 
      comblabels[count++] = labels[ii][jj]+1;  // get 1 offset not 0
    leafcount[ii] = static_cast<int>(labels[ii].size());
  }
  
  int nedges=static_cast<int>(edges.size());
  if (nedges != 2*(len-1))
    throw std::range_error("problem in C++ code split_simple\n");
  
  Rcpp::IntegerMatrix edge(nedges, 2);
  for (int i=0; i<nedges; i++) {
    edge(i,0) = edges[i].first;
    edge(i,1) = edges[i].second;
  }
  // now try to get the lengths.  Note that the centre of this split is positions[0].
  std::vector<int> nodePos;
  s.getNodesPositions(nodePos);
  Rcpp::IntegerVector nodepos(nodePos.begin(), nodePos.end());
  for (size_t ii=0;ii<nodePos.size();ii++) 
    nodepos[ii]=nodePos[ii];
  
  return  Rcpp::List::create(
      Rcpp::Named("edge") = edge,
      Rcpp::Named("labels") = comblabels,
      Rcpp::Named("nleaves") = len,
      Rcpp::Named("nodepos") = nodepos,
      Rcpp::Named("leafcount") = leafcount
    );
}
  
