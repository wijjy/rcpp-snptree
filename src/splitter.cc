#include "splitter.h"

/** Do all the SNPs at position col for rows rows match?                         */
bool SNPmatch(const Rcpp::IntegerMatrix &hap, const vector<int> &rows, int col) 
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

/** Note that both of these should be sorted                               */
int count_intersection(vector<int> &a, Rcpp::IntegerVector &b)
{
  std::vector<int>::iterator al=a.end(),ii=a.begin();
  Rcpp::IntegerVector::iterator bl=b.end(),jj=b.begin();
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
        (*ii)->left =new binode(lab[0], *ii);
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
  leaves.splice(leaves.begin(), toadd);
  return change;
}


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
    ncc[i]=count_intersection(lab[i], cases, nc);
    ncc[start+i] = ncc[i];
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
void splitter::getCaseControlLeaves(Rcpp::IntegerMatrix &ncc, Rcpp::IntegerVector &cases) 
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

void splitter::getTipLengths(int StartingPoint, Rcpp::IntegerMatrix &posLRtip) 
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



/** Are any of the labels not cases?  - an improved and quicker version  */
template<typename T>
bool mismatch2(vector<T> &labels, T *cases, int nc)
{
  int count=count_intersection(labels, cases, nc);
  if (count==static_cast<int>(labels.size())) return false;
  else if (count==0) return false;
  return true;
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
