/**                         @file                                           */
#ifndef SPLITTER_H__
#define SPLITTER_H__

#include <list>
#include <algorithm>
#include <stdexcept>
#include "rcpp_binode.h"
#include "testStats.h"
//#include "util.h"

using std::list;
using std::ostream;
using std::vector;
using std::pair;
using std::stack;

bool SNPmatch( const IntegerMatrix &hap, const vector<int> &rows, int col);

template<typename T>
bool mismatch(std::vector<T> &labels, T *cases, int nc);

template<typename T>
bool mismatch2(vector<T> &labels, T *cases, int nc);

int count_intersection(vector<int> &a, IntegerVector &b);
/** Note that both of these should be sorted                               */
int count_intersection(vector<int> &a, int *cases, int nc);

/** A class to create a splitting haplotype tree *
 * The constructor just creates the root.  
 * The split function splits the leaves depending on the 
 * haplotypes that you have.    */
class splitter {
public:
  /** Constructor for the splitter class                 */
  splitter( IntegerMatrix haplotypes) 
    :haps(haplotypes)    ,samples(haplotypes.nrow())  ,nSNP(haplotypes.ncol()) {
    vector<int> a(haplotypes.nrow());
    for (int i=0; i<haplotypes.nrow(); i++) a[i]=i;
    leaves.push_back(new binode(a));
    root_=leaves.front();
  }
  binode *root() const {
    return root_;
  }
  /** Return the number of leaves in a tree               */
  int nleaves() const {
    return leaves.size();
  }
  /** calculate the top and bottom positions of each node 
   * for drawing.  gap is in proportion to unity.  
  */
  void calculate_top_bottom(double gap);
  /** Split the haplotypes at a position                  */
  bool split(int position);
  /** Do some splits                                      */
  bool fullsplit(); 
  /** Print the labels at all the nodes                   */
  ostream & printLeaves(ostream &o) {
      list<binode *>::iterator ii=leaves.begin();
      while (ii!=leaves.end()) {
        (*ii++)->printlabels(o);
        o << std::endl; 
      }
      return o;
  }
  /** get the coordinates for a bifurcation tree */
  void get_coordinates(NumericMatrix coords) {
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
  void apesplit(vector<pair<int,int> >  &e, vector<vector<int> > &lab) {
    root_->recurseApeSplit(e,lab,nleaves()+1);
  }
  /** An alternative apesplit that gets the numbers of cases and controls
   * for all terminal nodes                                               */
  void apesplit(int *edges, int *ncc, int *cases,int nc);
  /** Alternative to apesplti for building a bifurcating tree             */
  void edges_positions_counts(int *edge1, int *edge2, int *count, int *position, int end_pos);
    
  /** Get the numbers of cases and controls at the internal nodes and the 
      leaves - in lexical order */
  void getCaseControlNodes(int *ncc, int *cases, int nc,int ninternal);
  void getCaseControlLeaves(IntegerMatrix &ncc,  IntegerVector &cases);
  /** Get the index of splits for each node - that is which of the SNPs the node was split on  */
  void getNodesPositions(std::vector<int> &pos);
  /** Get the labels at internal nodes in ape format */
  void getNodesLabels(std::vector<std::vector<int> > &labs);
  /** First need to think about the lengths and check them   */
  void getLengths(int CentrePos,int *posLRnodes, int *posLRtip);
  void getTipLengths(int StartingPoint, IntegerMatrix &posLRtip);
  
  /** Get the leaves (it doesn'r matter is this is NLR or LRN   */
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
  /** Get a vector of internal nodes in NLR order                  */
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
  std::vector<std::pair<int,double> > maxLengths(int k, int nmax); 
  
private:
  binode *root_;                            /// The root of the tree
  std::list<binode *> leaves;               /// A list of leaves
  std::list<binode *> internal;             /// A list of internal nodes
  IntegerMatrix haps;                       /// A pointer to a matrix of the data
  int samples, nSNP;                        /// the number of samples and SNPs
};

void splitter::apesplit(int *edges, int *ncc, int *cases,int nc)
  {
    int pos=nleaves();
    vector<pair<int,int> > e;
    vector<vector<int> > lab; 
    root_->recurseApeSplit(e,lab,pos+1);

    int start=2*(pos-1);
    assert(start==e.size());
    for (size_t i=0;i<e.size();i++) {
      edges[i]=e[i].first;
      edges[start+i] = e[i].second;
    }
    start=pos;
    assert(pos==lab.size());
    for (size_t i=0;i<lab.size();i++) {
      int n=lab[i].size();
      ncc[i]=count_intersection(lab[i], cases, nc);
      ncc[start+i] = ncc[i];
    }
  }
/** get edges counts and positions for a split tree                        */

void splitter::edges_positions_counts(int *edge1, int *edge2, int *count, int *position, int end_pos)
{
  int pos = nleaves();
  
  std::vector<std::vector<int> > yyy;
  int next_leaf=1;
 
  root_->recurse_edge_count_position(yyy, pos+1, next_leaf, end_pos);
  std::cerr << "size of yyy: " << yyy.size() << std::endl;
  assert(yy.size()==2*(pos-1));
  for (size_t ii=0; ii<yyy.size();ii++) {
    edge1[ii]    = yyy[ii][0];
    edge2[ii]    = yyy[ii][1];
    count[ii]    = yyy[ii][2];
    position[ii] = yyy[ii][3];
  }
}

/** Get the number of cases and controls at the nodes - this is written as 
 * to be used by the R interface   */
void  splitter::getCaseControlNodes(int *ncc, int *cases, int nc,int ninternal) 
  {
    NLRIterator<binode> ii(root_);
    int count=0;
    while (!ii.isend()) {
      ncc[count]=count_intersection((*ii)->labels, cases, nc);
      ncc[ninternal+count] = (*ii)->labels.size()-ncc[count];
      count++;
      ii.nextInternal();
    }
  }
/** Geta vector of node labels - this is written as 
 * to be used by the R interface so is in ape order  */
void splitter::getNodesLabels(std::vector<std::vector<int> > &labs) 
  {
    NLRIterator<binode> ii(root_);
    while (!ii.isend()) {
      labs.push_back((*ii)->labels);
      ii.nextInternal();
    }
  }
  
  /** Get a vector of positions - this is written as 
 * to be used by the R interface so is in ape order  */
void splitter::getNodesPositions(std::vector<int> &pos) 
  {
    NLRIterator<binode> ii(root_);
    while (!ii.isend()) {
      pos.push_back((*ii)->position);
      ii.nextInternal();
    }
  }




/** Get the number of cases and controls at the leaves - this should be in 
 * lexical search order */
void splitter::getCaseControlLeaves(IntegerMatrix &ncc, IntegerVector &cases) 
  {
    NLRIterator<binode> ii(root_);
    ii.nextLeaf();   // the root can never be a leaf
    int count=0;
    while (!ii.isend()) {
      ncc(count, 0) = count_intersection((*ii)->labels, cases);
      ncc(count, 1) = (*ii)->labels.size()-ncc(count, 0);
      count++;
      ii.nextLeaf();
    }
  }

/** Get the length of shared haplotypes at a set of nodes.  Need to search 
 * below the node to work out the minimum to the left and to the right 
 * and return special values if these do not exist (i.e.  the haplotypes
 * stretch to the end of the section of DNA that is being examined.
 *
 * CentrePos is the position of the node above that causes this split
 */
void splitter::getLengths(int CentrePos,int *posLRnode, int *posLRtip) 
  {
    NLRIterator<binode> ii(root_);
    int countNodes=0,countTips=0;
    while (!ii.isend()) {
      if ((*ii)->labels.size()==1) {  
        // a single node at this split.  shared length is 0 or the entire choice of 
        // positions - leave this for another function to decide, and return -1 for 
        // both left and right
        posLRtip[countTips++]=-1;
        posLRtip[countTips++]=-1;
      } else {
        int Left,Right;
        size_t jj;
        std::vector<int> &labs=(*ii)->labels;
        // find the maximum value to the left that is shared
        for (Left=CentrePos;Left>=0;Left--) {
          int MatchSNP=haps(labs[0], Left);
          for (jj=1;jj<labs.size();jj++) {
            if (haps(labs[jj], Left)!=MatchSNP) break;
          }
          // we we have reached the end of labels and there was no split
          // then continue - otherwise break
          if (jj!=labs.size()) break;
        }
        // Left should hold the maximum split to the left
        for (  Right=CentrePos;Right<nSNP;Right++) {
          int MatchSNP=haps(labs[0], Right);
          for (jj=1;jj<labs.size();jj++) {
            if (haps(labs[jj], Right) != MatchSNP) break;
          }
          if (jj!=labs.size()) break;
        }
        // Right should hold the maximum split to the right
        // or nSNP if it has reached the end
        if ((*ii)->isleaf()) {
          posLRtip[countTips++]=Left;
          posLRtip[countTips++]=Right;
        } else {
          posLRnode[countNodes++]=Left;
          posLRnode[countNodes++]=Right;
        }
      }
      ++ii;
    }
  }

void splitter::getTipLengths(int StartingPoint, IntegerMatrix &posLRtip) 
{
  NLRIterator<binode> ii(root_);
  ii.nextLeaf();
  int countTip=0;
  while (!ii.isend()) {
    if ((*ii)->labels.size()==1) {  
      // a single node at this split.  shared length is 0 or the entire choice of 
      // positions - leave this for another function to decide, and return -1 for 
      // both left and right
      posLRtip(countTip, 0) = -1;
      posLRtip(countTip++, 1) = -1;
    } else {
      int Left,Right;
      size_t jj;
      std::vector<int> &labs=(*ii)->labels;
      // find the maximum value to the left that is shared
      for (Left=StartingPoint;Left>=0;Left--) {
        int MatchSNP=haps(labs[0], Left);
        for (jj=1;jj<labs.size();jj++) {
          if (haps(labs[jj], Left) != MatchSNP) break;
        }
        // we we have reached the end of labels and there was no split
        // then continue - otherwise break
        if (jj!=labs.size()) break;
      }
      // Left should hold the maximum split to the left
      for (  Right=StartingPoint;Right<nSNP;Right++) {
        int MatchSNP=haps(labs[0], Right);
        for (jj=1;jj<labs.size();jj++) {
          if (haps(labs[jj], Right)!=MatchSNP) break;
        }
        if (jj!=labs.size()) break;
      }
      // Right should hold the maximum split to the right
      // or nSNP if it has reached the end
      posLRtip(countTip, 0) = Left;
      posLRtip(countTip++, 1) = Right;
    }
    ii.nextLeaf();
  }
}

/** Split a binode based on the SNP at position 
 * 
 * Note that the Left branch is always SNP[position] = 0
 */
bool splitter::split(int position) {
  bool change=false;
  list<binode *>::iterator ii=leaves.begin();
  std::list<binode *> toadd;
  while (ii!=leaves.end()) {
    if ((*ii)->labels.size()>1) {
      if (!SNPmatch(haps, (*ii)->labels, position)) {
        vector< vector<int> > lab(2);
        for (size_t jj=0;jj<(*ii)->labels.size();jj++) {
          lab[haps((*ii)->labels[jj], position)].push_back((*ii)->labels[jj]);
        }
        (*ii)->left=new binode(lab[0], *ii);
        (*ii)->right=new binode(lab[1], *ii);
        (*ii)->position=position;
        toadd.push_back((*ii)->left);
        toadd.push_back((*ii)->right);
        internal.push_back(*ii);
        ii=leaves.erase(ii);
        change=true;
      } else ii++;
    } else ii++;
  }
  leaves.splice(leaves.begin(),toadd);
  return change;
}
/** Function to split the all nodes with >1  
  repeated use of this will produce full splits   */
bool splitter::fullsplit() {
  bool change=false;
  list<binode *>::iterator ii=leaves.begin();
  std::list<binode *> toadd;
  while (ii!=leaves.end()) {
    if ((*ii)->labels.size()>1) {
      vector< vector<int> > lab(2);
      lab[0].push_back((*ii)->labels[0]);
      for (size_t jj=1;jj<(*ii)->labels.size();jj++) {
	lab[1].push_back((*ii)->labels[jj]);
      }
      (*ii)->left=new binode(lab[0]);
      (*ii)->right=new binode(lab[1]);
      (*ii)->position=-99;
      toadd.push_back((*ii)->left);
      toadd.push_back((*ii)->right);
      internal.push_back(*ii);
      ii=leaves.erase(ii);  // remove 
      change=true;
    } else ii++;
  }
  leaves.splice(leaves.begin(),toadd);
  return change;
}
/** Do all the SNPs at position col for rows rows match?                         */
bool SNPmatch(const IntegerMatrix &hap, const vector<int> &rows, int col) 
{
  int SNP = hap(rows[0], col);
  for (size_t j=0;j<rows.size();j++) 
    if (hap(rows[j], col) != SNP) return false;
  return true;
} 
/** Are any of the labels not cases?                                             */
template<typename T>
bool mismatch(vector<T> &labels, T *cases, int nc)
{
  typename vector<T>::iterator curr=labels.begin();
  typename vector<T>::iterator end=labels.end();
  typename vector<T>::iterator fo=std::find_first_of(curr,end,cases,cases+nc);
  if (fo==end) return false; // no labels in cases - all controls
  if (fo!=curr) return true; // first label not a case
  curr++;

  while (curr != end) {
    // find the first label (from curr) that is in the cases
    vector<int>::iterator fo=std::find_first_of(curr,end,cases,cases+nc);
    if (fo!=curr)  //the first label is not  a case
      return true;
    curr++;
  }
  //curr==end so all labels must be in cases!
  return false;
}

/** Are any of the labels not cases?  - an improved and quicker version  */
template<typename T>
bool mismatch2(vector<T> &labels, T *cases, int nc)
{
  int count=count_intersection(labels, cases, nc);
  if (count==static_cast<int>(labels.size())) return false;
  else if (count==0) return false;
  return true;
}
/** Note that both of these should be sorted                               */
int count_intersection(vector<int> &a, IntegerVector &b)
{
  std::vector<int>::iterator al=a.end(),ii=a.begin();
  IntegerVector::iterator bl=b.end(),jj=b.begin();
  int count=0;
  while (ii!=al && jj!=bl) {
    if (*ii<*jj) ++ii;
    else if (*jj<*ii) ++jj;
    else {
      count++;
      ii++;
      jj++;
    }
  }
  return count;
} 
/** Note that both of these should be sorted                               */
int count_intersection(vector<int> &a, int *cases, int nc)
{
  std::vector<int>::iterator al=a.end(),ii=a.begin();
  int *jj=cases,*bl=cases+nc;
  int count=0;
  while (ii!=al && jj!=bl) {
    if (*ii<*jj) ++ii;
    else if (*jj<*ii) ++jj;
    else {
      count++;
      ii++;
      jj++;
    }
  }
  return count;
}

std::vector<std::pair<int,double> >
  splitter::maxLengths(int k, int nmax) {
     std::vector<std::pair<int,double> > a;
     return a;	
  }
  
  /** calculate the top and bottom positions of each node 
   * for drawing.  gap is in proportion to unity.  
   */
  void splitter::calculate_top_bottom(double gap) {
    NLRIterator<binode> ii(root());
    ii.nextLeaf();   // the root can never be a leaf
    std::pair<double, double> last(0.0, (*ii)->labels.size());
    (*ii)->range = last;
    ii.nextLeaf();
    while (!ii.isend()) {
      (*ii)->range.first = gap + last.second;
      (*ii)->range.second = (*ii)->range.first + (*ii)->labels.size();
      last = (*ii)->range;
      ii.nextLeaf();
    }
    // now recurse to get the rest of the nodes
    root_->recurse_calculate_top_bottom();
  }


  
#endif
