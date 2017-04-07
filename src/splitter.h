/**                         @file                                           */
#ifndef SPLITTER_H__
#define SPLITTER_H__

#include <list>
#include <algorithm>
#include <stdexcept>
#include "rcpp_binode.h"
#include "testStats.h"
#include "util.h"

using std::list;
using std::ostream;
using std::vector;
using std::pair;
using std::stack;


/** A class to create a splitting haplotype tree *
 * The constructor just creates the root.  
 * The split function splits the leaves depending on the 
 * haplotypes that you have.    */
class splitter {
public:
  /** Constructor for the splitter class                 */
  splitter(const Rcpp::IntegerMatrix haplotypes) 
    :haps(haplotypes),
     samples(haplotypes.nrow()),
     nSNP(haplotypes.ncol()) {
    vector<int> root_labels(haplotypes.nrow());
    for (int i=0; i<haplotypes.nrow(); i++) root_labels[i]=i;
    leaves.push_back(new binode(root_labels));
    root_=leaves.front();
  }
  /** Constructor for the splitter class                 */
  splitter(const splitter &samedata) 
    :haps(samedata.haps),
     samples(samedata.samples),
     nSNP(samedata.nSNP) {
    vector<int> root_labels(samples);
    for (int i=0; i<samples; i++) root_labels[i]=i;
    leaves.push_back(new binode(root_labels));
    root_=leaves.front();
  }
  /** Return a pointer to the root                        */
  binode *root() const {
    return root_;
  }
  /** Return the number of leaves in a tree               */
  int nleaves() const {
    return leaves.size();
  }
  /** index tree                                          */
  void index_tree() {
    int count=0;
    
    std::list<binode *>::iterator bi=begin_leaf();
    while (bi!=end_leaf()) {
      (*bi)->index=count++;
    }
    bi=begin_internal();
    while (bi!=end_internal()) {
      (*bi)->index=count++;
    }
  }
  /** LRN index tree                                          */
  void LRN_index_tree() {
    int count=0;
    LRNIterator<binode> ii(root_); 
    while (!ii.isend()) {
      (*ii) -> LRN_index = count++;
      ++ii;
    }
  }
  /** calculate the top and bottom positions of each node 
   * for drawing.  gap is in proportion to unity.  
  */
  void calculate_top_bottom(double gap);
  /** calculate the top and bottom positions of each node 
   * for a subset of the data marked by individual.
   * gap is in proportion to unity.  
   * These values should be in the centre for joins, at the 
   * top fpr left branches and at the top for right branches
   */
  void calculate_id_top_bottom(const Rcpp::IntegerVector &ind);
  /** Split the haplotypes at a position                  */
  bool split(int position);
  /** Do some splits                                      */
  bool fullsplit(); 

  /** Print the labels at all the nodes                    */
  ostream & printLeaves(ostream &o) {
      list<binode *>::iterator ii=leaves.begin();
      while (ii!=leaves.end()) {
        (*ii++)->printlabels(o);
        o << std::endl; 
      }
      return o;
  }
  /** get the coordinates for a bifurcation tree */
  void get_coordinates(Rcpp::NumericMatrix coords) {
   // Rprintf(" root position = %d\n", root()->position);
  //  Rprintf(" root range is (%g, %g)\n", root()->range.first, root()->range.second);
    coords(0, 0) = root_->position;
    coords(0, 1) = root_->range.first;
    coords(0, 2) = 0;  // squared off at the root
    int index=1;
    root_->left->recursively_get_coords(coords, index);
    coords(index, 0) = root_->position;
    coords(index, 1) = root_->range.first+(root_->range.second-root_->range.first)*root_->left->labels.size()/static_cast<double>(root_->labels.size());
    coords(index++, 2) = 0.0;//0.5;
    root()->right->recursively_get_coords(coords, index);
    coords(index, 0) = root()->position;
    coords(index, 1) = root()->range.second;
    coords(index, 2) = 0;
  }
  /** The first test statistic                            */
  double testStat1(vector<int> &cases) {
    return testStat1(&cases[0],static_cast<int>(cases.size()));
  }
  int testStat1(int *cases, int nc) {
    int count=0;
    list<binode *>::iterator ii=leaves.begin();
    while (ii!=leaves.end()) {
      if ((*ii)->labels.size()>1) {
        count += mismatch2((*ii)->labels,cases,nc);
      }
      ii++;
    }
    ii=internal.begin();
    while (ii!=internal.end()) {
      // look for mismatch nodes
      count += mismatch2((*ii)->labels,cases,nc);
      ii++;
    }
    return count;
  }
  
  std::vector<double> qtlStat( const vector<double> &qvals, int maxk=5, const char *pick="Z") {
    double meanvar[2];
    avevar(qvals, meanvar[0], meanvar[1]);
    double sumx;
    switch(pick[0]) {
      case 'Z': return root_->RecurseZStat2(qvals, sumx, maxk, meanvar);
      case 'A': return root_->RecurseZStatA(qvals, sumx, maxk, meanvar);
      case 'P': return root_->RecurseZStatP(qvals, sumx, maxk, meanvar);
      case 'N': return root_->RecurseZStatN(qvals, sumx, maxk, meanvar);

    default: 
      std::ostringstream oss;
      oss << "The statistic given by '" << pick << "' is not supported yet\n";
      throw std::runtime_error(oss.str());
    }
  }

  std::vector<double> getStat( const vector<int> &cases,int maxk=5,const char *pick="S") {
    double pp=cases.size()/static_cast<double>(samples);
    int CC[2],d[2];
    switch(pick[0]) {
    case 'S': return root_->SevonTestStatisticsFaster(cases,pp,maxk,CC);
    case 'A': return root_->RecurseTestStatistics(cases,pp,maxk,CC,AbsSevonStat);
    case 'Q': return root_->RecurseTestStatistics(cases,pp,maxk,CC,SqSevonStat);
    case 'G': return root_->RecurseTestStatistics(cases,pp,maxk,CC,GStat);
    case 'N': return root_->RecurseTestStatistics(cases,pp,maxk,CC,LogPNorm);   
    case 'P': return root_->RecurseTestStatistics(cases,pp,maxk,CC,LogP);
    case 'C':  {
      std::vector<int> tmp = root_->RecurseCherries(cases);
      std::vector<double> res(8);
      res[0] = tmp[0]-tmp[1];
      res[1] = tmp[0]+tmp[2] - tmp[1];
      res[2] = static_cast<double>(tmp[0])/(tmp[0]+tmp[1]+tmp[2]);
      res[3]=tmp[0];
      res[4]=tmp[1];
      res[5]=tmp[2];
      res[6]=static_cast<double>(root_->labels.size());      // total
      res[7]=static_cast<double>(cases.size());             // ncases
      return res;
    }
    case 'H':  {
      std::vector<double> tmp = root_->RecurseSumHeight(cases,CC);
      tmp[0] /= CC[0];
      tmp[1] /= CC[1]=CC[0];
      std::vector<double> res(4);
      res[0] = tmp[0]-tmp[1];
      res[1] = tmp[0]/(tmp[0]+tmp[1]);
      res[2] = tmp[0];
      res[3] = tmp[1];
      return res;
    }
    case 'T': {
      std::vector<double> tmp = root_->TreeDistanceStatistic(cases,CC,d);
      //      std::cerr << tmp[0] << ' " << tmp[1] << ' " << tmp[2] << std::endl;
      std::vector<double> res(4);
      res[0] = tmp[2]/static_cast<double>(CC[1]) - tmp[0]/static_cast<double>(CC[0]);
      res[1] = tmp[1]/static_cast<double>(CC[1]-CC[0]) -  tmp[0]/static_cast<double>(CC[0]);
      res[2] = (tmp[2]-tmp[1]-tmp[0])/static_cast<double>(CC[1]);
        //     res[2] = (tmp[0]/CC[0])/(tmp[1]/CC[1]);
      res[3] = tmp[2]/(1.0+tmp[1]+tmp[0]);
      //res[2] = 
      // res[2] = tmp[2]/static_cast<double>(CC[1]) - tmp[1]/static_cast<double>(CC[1]-CC[0]);
      //res[3] = tmp[1];//  tmp[2]/static_cast<double>(CC[1]) - (tmp[1]+tmp[0])/static_cast<double>(CC[1]);
      return res;
    }
    case 't': {
      root_->RecurseLeftRightDistances(cases);
      root_->distance(0, 0) = 9999;
      root_->distance(0, 1) = 9999;
      root_->RecurseUpDistances();
      std::vector<double> tmp(6,0.0);
      LRNIterator<binode> i(root_);  

      while (!i.isend()) {

        if ((*i)->localCC[0]==1)                   // Cases to Cases
          tmp[0] += (*i)->distance(0, 0);
     
        if ((*i)->localCC[1]-(*i)->localCC[0]==1)  // Controls to Controls
          tmp[1] += (*i)->distance(1, 0);
       
        if ((*i)->localCC[1]==(*i)->localCC[0])    // Cases to Controls
          tmp[2] +=  (*i)->localCC[0]*(*i)->distance(1, 0);
        
        if ((*i)->localCC[0]==0)                   // Controls to Cases
          tmp[3] += ((*i)->localCC[1]-(*i)->localCC[0])*(*i)->distance(0, 0);
       
        i.nextLeaf();
      }
      
      std::vector<double> res(6);
      CC[0] = root_->localCC[0];
      CC[1] = root_->localCC[1];
 
      res[0] =  (tmp[2]-tmp[0])/static_cast<double>(CC[0]);                  // clustering  of cases
      res[1] =  (tmp[3]-tmp[1])/static_cast<double>(CC[1]-CC[0]);            // clustering of controls
      res[2] =  (tmp[3]+tmp[2]-tmp[1]-tmp[0])/static_cast<double>(CC[1]);    // Clustering of cases and controls
      res[3] =  (tmp[2])/(1.0+tmp[0]);                                       // ratio Case Clustering
      res[4] =  (tmp[3])/(1.0+tmp[1]);                     
      res[5] =  (tmp[2]+tmp[3])/(1.0+tmp[0]+tmp[1]);
       
      return res;
      
    }
    default: 
      std::ostringstream oss;
      oss << "The statistic given by '" << pick << "' is not supported yet\n";
      throw std::runtime_error(oss.str());
    }
  }

  std::vector<std::pair<double, std::vector< binode *> > > 
        getNodes( const vector<int> &cases,int maxk=5, char pick='P') {
    double pp=cases.size()/static_cast<double>(samples);
    int CC[2];
    switch(pick) {
      //   case "S": return root_->SevonTestStatisticsFaster(cases,pp,maxk,CC);
    case 'A': return root_->RecurseTestNodes(cases,pp,maxk,CC,AbsSevonStat);
    case 'Q': return root_->RecurseTestNodes(cases,pp,maxk,CC,SqSevonStat);
    case 'G': return root_->RecurseTestNodes(cases,pp,maxk,CC,GStat);
    case 'P': return root_->RecurseTestNodes(cases,pp,maxk,CC,LogP);
    default: 
      std::ostringstream oss;
      oss << "The statistic given by "" << pick << "" is not supported yet in getNodes\n";
      throw std::runtime_error(oss.str());
    }
  }

  /** The Sevond et al 2006 test statistics - the set of k trees with highest 
   * sum of Z scores                                   */
  std::vector<double> recurseSevonTestStatistic( const vector<int> &cases,int maxk=5) {
    double pp=cases.size()/static_cast<double>(samples);
    int CC[2];
    return root_->SevonTestStatisticsFaster(cases,pp,maxk,CC);
  }
 
  /** A function to get the edges for the tree in a format that 
   * is understandable by ape.  The general method should perhaps
   * be reworked using my iterator for stepping through trees!
   */
  void apesplit(vector<pair<int,int> > &e, vector<vector<int> > &lab) {
    root_->recurseApeSplit(e, lab, nleaves()+1);
  }

  /** An alternative apesplit that gets the numbers of cases and controls
   * for all terminal nodes                                               */
  void apesplit(int *edges, int *ncc, int *cases,int nc);
  /** Alternative to apesplti for building a bifurcating tree             */
  void edges_positions_counts(int *edge1, int *edge2, int *count, int *position, int end_pos);
    
  /** Get the numbers of cases and controls at the internal nodes and the 
      leaves - in lexical order */
  void getCaseControlNodes(int *ncc, int *cases, int nc,int ninternal);
  void getCaseControlLeaves(Rcpp::IntegerMatrix &ncc, const Rcpp::IntegerVector &cases);
  /** Get the index of splits for each node - that is which of the SNPs the node was split on  */
  void getNodesPositions(std::vector<int> &pos);
  /** Get the labels at internal nodes in ape format */
  void getNodesLabels(std::vector<std::vector<int> > &labs);
  /** First need to think about the lengths and check them   */
  void getLengths(int CentrePos,int *posLRnodes, int *posLRtip);
  void getTipLengths(int StartingPoint, Rcpp::IntegerMatrix &posLRtip);
  
  /** Get a vector of pointers to leaves (it doesn't matter is this is NLR or LRN   */
  std::vector<binode *> GetLeaves() {
    std::vector<binode *> xx;
    xx.reserve(nleaves());                 // reserve the correct length
    NLRIterator<binode> ii(root_); 
    ii.nextLeaf();                         // get the first leaf
    while (!ii.isend()) {
      xx.push_back(*ii);
      ii.nextLeaf();
    }
    return xx;
  }
  /** Get a vector of pointers to internal nodes in NLR order               */
  std::vector<binode *> GetLRNInternal() {
    std::vector<binode *> xx;
    xx.reserve(nleaves()-1);                 // reserve the correct length
    LRNIterator<binode> ii(root_); 
    ii.nextInternal();
    while (!ii.isend()) {
      xx.push_back(*ii);
      ii.nextInternal();
    }
    return xx;
  }
  /** Get a vector of internal nodes in NLR order                  */
  std::vector<binode *> GetNLRInternal() {
    std::vector<binode *> xx;
    xx.reserve(nleaves()-1);                 // reserve the correct length
    NLRIterator<binode> ii(root_); 
    while (!ii.isend()) {
      xx.push_back(*ii);
      ii.nextInternal();
    }
    return xx;
  }
  /** Get the set of maximum average lengths (for nodes with > k individuals */
  std::vector<std::pair<int,double> >  maxLengths(int k, int nmax) {
      std::vector<std::pair<int,double> > a;
      return a;	
    }
  std::list<binode *>::iterator begin_leaf() {
    return leaves.begin();
  }
  std::list<binode *>::iterator end_leaf() {
    return leaves.end();
  }
  std::list<binode *>::const_iterator begin_leaf() const {
    return leaves.begin();
  }
  std::list<binode *>::const_iterator end_leaf() const {
    return leaves.end();
  }
  std::list<binode *>::iterator begin_internal() {
    return internal.begin();
  }
  std::list<binode *>::iterator end_internal() {
    return internal.end();
  }
  std::list<binode *>::const_iterator begin_internal() const {
    return internal.begin();
  }
  std::list<binode *>::const_iterator end_internal() const {
    return internal.end();
  }
  binode *first_leaf_matching(const Rcpp::IntegerVector targetPositions, 
                                    const    Rcpp::IntegerVector &target);
private:
  binode *root_;                            /// The root of the tree
  std::list<binode *> leaves;               /// A list of leaves
  std::list<binode *> internal;             /// A list of internal nodes
  Rcpp::IntegerMatrix haps;                 /// A pointer to a matrix of the data
  int samples, nSNP;                        /// the number of samples and SNPs
};




#endif
