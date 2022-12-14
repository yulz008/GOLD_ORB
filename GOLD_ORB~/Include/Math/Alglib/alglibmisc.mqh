//+------------------------------------------------------------------+
//|                                                   alglibmisc.mqh |
//|            Copyright 2003-2012 Sergey Bochkanov (ALGLIB project) |
//|                   Copyright 2012-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
//| Implementation of ALGLIB library in MetaQuotes Language 5        |
//|                                                                  |
//| The features of the library include:                             |
//| - Linear algebra (direct algorithms, EVD, SVD)                   |
//| - Solving systems of linear and non-linear equations             |
//| - Interpolation                                                  |
//| - Optimization                                                   |
//| - FFT (Fast Fourier Transform)                                   |
//| - Numerical integration                                          |
//| - Linear and nonlinear least-squares fitting                     |
//| - Ordinary differential equations                                |
//| - Computation of special functions                               |
//| - Descriptive statistics and hypothesis testing                  |
//| - Data analysis - classification, regression                     |
//| - Implementing linear algebra algorithms, interpolation, etc.    |
//|   in high-precision arithmetic (using MPFR)                      |
//|                                                                  |
//| This file is free software; you can redistribute it and/or       |
//| modify it under the terms of the GNU General Public License as   |
//| published by the Free Software Foundation (www.fsf.org); either  |
//| version 2 of the License, or (at your option) any later version. |
//|                                                                  |
//| This program is distributed in the hope that it will be useful,  |
//| but WITHOUT ANY WARRANTY; without even the implied warranty of   |
//| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the     |
//| GNU General Public License for more details.                     |
//+------------------------------------------------------------------+
#include "ap.mqh"
#include "alglibinternal.mqh"
//+------------------------------------------------------------------+
//| KD-trees                                                         |
//+------------------------------------------------------------------+
class CKDTree
  {
public:
   int               m_n;
   int               m_nx;
   int               m_ny;
   int               m_normtype;
   int               m_kneeded;
   double            m_rneeded;
   bool              m_selfmatch;
   double            m_approxf;
   int               m_kcur;
   double            m_curdist;
   int               m_debugcounter;
   //--- arrays
   int               m_tags[];
   double            m_boxmin[];
   double            m_boxmax[];
   int               m_nodes[];
   double            m_splits[];
   double            m_x[];
   int               m_idx[];
   double            m_r[];
   double            m_buf[];
   double            m_curboxmin[];
   double            m_curboxmax[];
   //--- matrix
   CMatrixDouble     m_xy;
   //--- constructor, destructor
                     CKDTree(void);
                    ~CKDTree(void);
   //--- copy
   void              Copy(CKDTree &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CKDTree::CKDTree(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CKDTree::~CKDTree(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CKDTree::Copy(CKDTree &obj)
  {
//--- copy variables
   m_n=obj.m_n;
   m_nx=obj.m_nx;
   m_ny=obj.m_ny;
   m_normtype=obj.m_normtype;
   m_kneeded=obj.m_kneeded;
   m_rneeded=obj.m_rneeded;
   m_selfmatch=obj.m_selfmatch;
   m_approxf=obj.m_approxf;
   m_kcur=obj.m_kcur;
   m_curdist=obj.m_curdist;
   m_debugcounter=obj.m_debugcounter;
//--- copy arrays
   ArrayCopy(m_tags,obj.m_tags);
   ArrayCopy(m_boxmin,obj.m_boxmin);
   ArrayCopy(m_boxmax,obj.m_boxmax);
   ArrayCopy(m_nodes,obj.m_nodes);
   ArrayCopy(m_splits,obj.m_splits);
   ArrayCopy(m_x,obj.m_x);
   ArrayCopy(m_idx,obj.m_idx);
   ArrayCopy(m_r,obj.m_r);
   ArrayCopy(m_buf,obj.m_buf);
   ArrayCopy(m_curboxmin,obj.m_curboxmin);
   ArrayCopy(m_curboxmax,obj.m_curboxmax);
//--- copy matrix
   m_xy=obj.m_xy;
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CKDTree                          |
//+------------------------------------------------------------------+
class CKDTreeShell
  {
private:
   CKDTree           m_innerobj;
public:
   //--- constructors, destructor
                     CKDTreeShell(void);
                     CKDTreeShell(CKDTree &obj);
                    ~CKDTreeShell(void);
   //--- method                
   CKDTree          *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CKDTreeShell::CKDTreeShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
CKDTreeShell::CKDTreeShell(CKDTree &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CKDTreeShell::~CKDTreeShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Return object of CKDTree             	                           |
//+------------------------------------------------------------------+
CKDTree *CKDTreeShell::GetInnerObj(void)
  {
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Build KD-trees                                                   |
//+------------------------------------------------------------------+
class CNearestNeighbor
  {
public:
   //--- class constants
   static const int  m_splitnodesize;
   static const int  m_kdtreefirstversion;
   //--- constructor, destructor
                     CNearestNeighbor(void);
                    ~CNearestNeighbor(void);
   //--- build
   static void       KDTreeBuild(CMatrixDouble &xy,const int n,const int nx,const int ny,const int normtype,CKDTree &kdt);
   static void       KDTreeBuildTagged(CMatrixDouble &xy,int &tags[],const int n,const int nx,const int ny,const int normtype,CKDTree &kdt);
   static int        KDTreeQueryKNN(CKDTree &kdt,double &x[],const int k,const bool selfmatch);
   static int        KDTreeQueryRNN(CKDTree &kdt,double &x[],const double r,const bool selfmatch);
   static int        KDTreeQueryAKNN(CKDTree &kdt,double &x[],int k,const bool selfmatch,const double eps);
   static void       KDTreeQueryResultsX(CKDTree &kdt,CMatrixDouble &x);
   static void       KDTreeQueryResultsXY(CKDTree &kdt,CMatrixDouble &xy);
   static void       KDTreeQueryResultsTags(CKDTree &kdt,int &tags[]);
   static void       KDTreeQueryResultsDistances(CKDTree &kdt,double &r[]);
   static void       KDTreeQueryResultsXI(CKDTree &kdt,CMatrixDouble &x);
   static void       KDTreeQueryResultsXYI(CKDTree &kdt,CMatrixDouble &xy);
   static void       KDTreeQueryResultsTagsI(CKDTree &kdt,int &tags[]);
   static void       KDTreeQueryResultsDistancesI(CKDTree &kdt,double &r[]);
   //--- serialize
   static void       KDTreeAlloc(CSerializer &s,CKDTree &tree);
   static void       KDTreeSerialize(CSerializer &s,CKDTree &tree);
   static void       KDTreeUnserialize(CSerializer &s,CKDTree &tree);
private:
   //--- private methods
   static void       KDTreeSplit(CKDTree &kdt,const int i1,const int i2,const int d,const double s,int &i3);
   static void       KDTreeGenerateTreeRec(CKDTree &kdt,int &nodesoffs,int &splitsoffs,const int i1,const int i2,const int maxleafsize);
   static void       KDTreeQueryNNRec(CKDTree &kdt,const int offs);
   static void       KDTreeInitBox(CKDTree &kdt,double &x[]);
   static void       KDTreeAllocDataSetIndependent(CKDTree &kdt,const int nx,const int ny);
   static void       KDTreeAllocDataSetDependent(CKDTree &kdt,const int n,const int nx,const int ny);
   static void       KDTreeAllocTemporaries(CKDTree &kdt,const int n,const int nx,const int ny);
  };
//+------------------------------------------------------------------+
//| Initialize constants                                             |
//+------------------------------------------------------------------+
const int CNearestNeighbor::m_splitnodesize=6;
const int CNearestNeighbor::m_kdtreefirstversion=0;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNearestNeighbor::CNearestNeighbor(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNearestNeighbor::~CNearestNeighbor(void)
  {

  }
//+------------------------------------------------------------------+
//| KD-tree creation                                                 |
//| This subroutine creates KD-tree from set of X-values and optional| 
//| Y-values                                                         |
//| INPUT PARAMETERS                                                 |
//|     XY      -   dataset, array[0..N-1,0..NX+NY-1].               |
//|                 one row corresponds to one point.                |
//|                 first NX columns contain X-values, next NY (NY   |
//|                 may be zero)                                     |
//|                 columns may contain associated Y-values          |
//|     N       -   number of points, N>=1                           |
//|     NX      -   space dimension, NX>=1.                          |
//|     NY      -   number of optional Y-values, NY>=0.              |
//|     NormType-   norm type:                                       |
//|                 * 0 denotes infinity-norm                        |
//|                 * 1 denotes 1-norm                               |
//|                 * 2 denotes 2-norm (Euclidean norm)              |
//| OUTPUT PARAMETERS                                                |
//|     KDT     -   KD-tree                                          |
//| NOTES                                                            |
//| 1. KD-tree  creation  have O(N*logN) complexity and              |
//|    O(N*(2*NX+NY)) memory requirements.                           |
//| 2. Although KD-trees may be used with any combination of N and   |
//|    NX, they are more efficient than brute-force search only when |
//|    N >> 4^NX. So they are most useful in low-dimensional tasks   |
//|    (NX=2, NX=3). NX=1  is another inefficient case, because      |
//|    simple  binary  search  (without  additional structures) is   |
//|    much more efficient in such tasks than KD-trees.              |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeBuild(CMatrixDouble &xy,const int n,
                                          const int nx,const int ny,
                                          const int normtype,CKDTree &kdt)
  {
//--- create a variable
   int i=0;
//--- creating array
   int tags[];
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(nx>=1,__FUNCTION__+": NX<1!"))
      return;
//--- check
   if(!CAp::Assert(ny>=0,__FUNCTION__+": NY<0!"))
      return;
//--- check
   if(!CAp::Assert(normtype>=0&&normtype<=2,__FUNCTION__+": incorrect NormType!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(xy)>=n,__FUNCTION__+": rows(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(xy)>=nx+ny,__FUNCTION__+": cols(X)<NX+NY!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(xy,n,nx+ny),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- allocation
   ArrayResizeAL(tags,n);
   for(i=0;i<=n-1;i++)
      tags[i]=0;
//--- function call
   KDTreeBuildTagged(xy,tags,n,nx,ny,normtype,kdt);
  }
//+------------------------------------------------------------------+
//| KD-tree creation                                                 |
//| This subroutine creates KD-tree from set of X-values, integer    |
//| tags and optional Y-values                                       |
//| INPUT PARAMETERS                                                 |
//|     XY      -   dataset, array[0..N-1,0..NX+NY-1].               |
//|                 one row corresponds to one point.                |
//|                 first NX columns contain X-values, next NY (NY   |
//|                 may be zero)                                     |
//|                 columns may contain associated Y-values          |
//|     Tags    -   tags, array[0..N-1], contains integer tags       |
//|                 associated with points.                          |
//|     N       -   number of points, N>=1                           |
//|     NX      -   space dimension, NX>=1.                          |
//|     NY      -   number of optional Y-values, NY>=0.              |
//|     NormType-   norm type:                                       |
//|                 * 0 denotes infinity-norm                        |
//|                 * 1 denotes 1-norm                               |
//|                 * 2 denotes 2-norm (Euclidean norm)              |
//| OUTPUT PARAMETERS                                                |
//|     KDT     -   KD-tree                                          |
//| NOTES                                                            |
//| 1. KD-tree  creation  have O(N*logN) complexity and              |
//|    O(N*(2*NX+NY)) memory requirements.                           |
//| 2. Although KD-trees may be used with any combination of N and   |
//|    NX, they are more efficient than brute-force search only when |
//|    N >> 4^NX. So they are most useful in low-dimensional tasks   |
//|    (NX=2, NX=3). NX=1 is another inefficient case, because simple|
//|    binary search (without additional structures) is much more    |
//|    efficient in such tasks than KD-trees.                        |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeBuildTagged(CMatrixDouble &xy,int &tags[],
                                                const int n,const int nx,
                                                const int ny,
                                                const int normtype,CKDTree &kdt)
  {
//--- create variables
   int i=0;
   int j=0;
   int maxnodes=0;
   int nodesoffs=0;
   int splitsoffs=0;
   int i_=0;
   int i1_=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(nx>=1,__FUNCTION__+": NX<1!"))
      return;
//--- check
   if(!CAp::Assert(ny>=0,__FUNCTION__+": NY<0!"))
      return;
//--- check
   if(!CAp::Assert(normtype>=0&&normtype<=2,__FUNCTION__+": incorrect NormType!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(xy)>=n,__FUNCTION__+": rows(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(xy)>=nx+ny,__FUNCTION__+": cols(X)<NX+NY!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(xy,n,nx+ny),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- initialize
   kdt.m_n=n;
   kdt.m_nx=nx;
   kdt.m_ny=ny;
   kdt.m_normtype=normtype;
//--- Allocate
   KDTreeAllocDataSetIndependent(kdt,nx,ny);
   KDTreeAllocDataSetDependent(kdt,n,nx,ny);
//--- Initial fill
   for(i=0;i<=n-1;i++)
     {
      for(i_=0;i_<=nx-1;i_++)
         kdt.m_xy[i].Set(i_,xy[i][i_]);
      i1_=-nx;
      for(i_=nx;i_<=2*nx+ny-1;i_++)
         kdt.m_xy[i].Set(i_,xy[i][i_+i1_]);
      kdt.m_tags[i]=tags[i];
     }
//--- Determine bounding box
   for(i_=0;i_<=nx-1;i_++)
      kdt.m_boxmin[i_]=kdt.m_xy[0][i_];
   for(i_=0;i_<=nx-1;i_++)
      kdt.m_boxmax[i_]=kdt.m_xy[0][i_];
   for(i=1;i<=n-1;i++)
     {
      for(j=0;j<=nx-1;j++)
        {
         kdt.m_boxmin[j]=MathMin(kdt.m_boxmin[j],kdt.m_xy[i][j]);
         kdt.m_boxmax[j]=MathMax(kdt.m_boxmax[j],kdt.m_xy[i][j]);
        }
     }
//--- prepare tree structure
//--- * MaxNodes=N because we guarantee no trivial splits,i.e.
//---   every split will generate two non-empty boxes
   maxnodes=n;
//--- allocation
   ArrayResizeAL(kdt.m_nodes,m_splitnodesize*2*maxnodes);
   ArrayResizeAL(kdt.m_splits,2*maxnodes);
   nodesoffs=0;
   splitsoffs=0;
   for(i_=0;i_<=nx-1;i_++)
      kdt.m_curboxmin[i_]=kdt.m_boxmin[i_];
   for(i_=0;i_<=nx-1;i_++)
      kdt.m_curboxmax[i_]=kdt.m_boxmax[i_];
//--- function call
   KDTreeGenerateTreeRec(kdt,nodesoffs,splitsoffs,0,n,8);
//--- Set current query size to 0
   kdt.m_kcur=0;
  }
//+------------------------------------------------------------------+
//| K-NN query: K nearest neighbors                                  |
//| INPUT PARAMETERS                                                 |
//|     KDT         -   KD-tree                                      |
//|     X           -   point, array[0..NX-1].                       |
//|     K           -   number of neighbors to return, K>=1          |
//|     SelfMatch   -   whether self-matches are allowed:            |
//|                     * if True, nearest neighbor may be the point |
//|                       itself (if it exists in original dataset)  |
//|                     * if False, then only points with non-zero   |
//|                       distance are returned                      |
//|                     * if not given, considered True              |
//| RESULT                                                           |
//|     number of actual neighbors found (either K or N, if K>N).    |
//| This  subroutine performs query and stores its result in the     |
//| internal structures of the KD-tree. You can use following        |
//| subroutines to obtain these results:                             |
//| * KDTreeQueryResultsX() to get X-values                          |
//| * KDTreeQueryResultsXY() to get X- and Y-values                  |
//| * KDTreeQueryResultsTags() to get tag values                     |
//| * KDTreeQueryResultsDistances() to get distances                 |
//+------------------------------------------------------------------+
static int CNearestNeighbor::KDTreeQueryKNN(CKDTree &kdt,double &x[],
                                            const int k,const bool selfmatch)
  {
//--- check
   if(!CAp::Assert(k>=1,__FUNCTION__+": K<1!"))
      return(-1);
//--- check
   if(!CAp::Assert(CAp::Len(x)>=kdt.m_nx,__FUNCTION__+": Length(X)<NX!"))
      return(-1);
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,kdt.m_nx),__FUNCTION__+": X contains infinite or NaN values!"))
      return(-1);
//--- return result
   return(KDTreeQueryAKNN(kdt,x,k,selfmatch,0.0));
  }
//+------------------------------------------------------------------+
//| R-NN query: all points within R-sphere centered at X             |
//| INPUT PARAMETERS                                                 |
//|     KDT         -   KD-tree                                      |
//|     X           -   point, array[0..NX-1].                       |
//|     R           -   radius of sphere (in corresponding norm), R>0|
//|     SelfMatch   -   whether self-matches are allowed:            |
//|                     * if True, nearest neighbor may be the point |
//|                       itself (if it exists in original dataset)  |
//|                     * if False, then only points with non-zero   |
//|                       distance are returned                      |
//|                     * if not given, considered True              |
//| RESULT                                                           |
//|     number of neighbors found, >=0                               |
//| This subroutine performs query and stores its result in the      | 
//| internal structures of the KD-tree. You can use following        |
//| subroutines to obtain actual results:                            |
//| * KDTreeQueryResultsX() to get X-values                          |
//| * KDTreeQueryResultsXY() to get X- and Y-values                  |
//| * KDTreeQueryResultsTags() to get tag values                     |
//| * KDTreeQueryResultsDistances() to get distances                 |
//+------------------------------------------------------------------+
static int CNearestNeighbor::KDTreeQueryRNN(CKDTree &kdt,double &x[],
                                            const double r,const bool selfmatch)
  {
//--- create variables
   int result=0;
   int i=0;
   int j=0;
//--- check
   if(!CAp::Assert((double)(r)>0.0,__FUNCTION__+": incorrect R!"))
      return(-1);
//--- check
   if(!CAp::Assert(CAp::Len(x)>=kdt.m_nx,__FUNCTION__+": Length(X)<NX!"))
      return(-1);
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,kdt.m_nx),__FUNCTION__+": X contains infinite or NaN values!"))
      return(-1);
//--- Prepare parameters
   kdt.m_kneeded=0;
//--- check
   if(kdt.m_normtype!=2)
      kdt.m_rneeded=r;
   else
      kdt.m_rneeded=CMath::Sqr(r);
//--- change values
   kdt.m_selfmatch=selfmatch;
   kdt.m_approxf=1;
   kdt.m_kcur=0;
//--- calculate distance from point to current bounding box
   KDTreeInitBox(kdt,x);
//--- call recursive search
//--- results are returned as heap
   KDTreeQueryNNRec(kdt,0);
//--- pop from heap to generate ordered representation
//--- last element is not pop'ed because it is already in
//--- its place
   result=kdt.m_kcur;
   j=kdt.m_kcur;
   for(i=kdt.m_kcur;i>=2;i--)
      CTSort::TagHeapPopI(kdt.m_r,kdt.m_idx,j);
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| K-NN query: approximate K nearest neighbors                      |
//| INPUT PARAMETERS                                                 |
//|     KDT         -   KD-tree                                      |
//|     X           -   point, array[0..NX-1].                       |
//|     K           -   number of neighbors to return, K>=1          |
//|     SelfMatch   -   whether self-matches are allowed:            |
//|                     * if True, nearest neighbor may be the point | 
//|                       itself (if it exists in original dataset)  |
//|                     * if False, then only points with non-zero   |
//|                       distance are returned                      |
//|                     * if not given, considered True              |
//|     Eps         -   approximation factor, Eps>=0. eps-approximate|
//|                     nearest neighbor is a neighbor whose distance|
//|                     from X is at most (1+eps) times distance of  |
//|                     true nearest neighbor.                       |
//| RESULT                                                           |
//|     number of actual neighbors found (either K or N, if K>N).    |
//| NOTES                                                            |
//|     significant performance gain may be achieved only when Eps is|
//|     on the order of magnitude of 1 or larger.                    |
//| This subroutine performs query and stores its result in the      |
//| internal structures of the KD-tree. You can use following        |
//| these subroutines to  obtain results:                            |
//| * KDTreeQueryResultsX() to get X-values                          |
//| * KDTreeQueryResultsXY() to get X- and Y-values                  |
//| * KDTreeQueryResultsTags() to get tag values                     |
//| * KDTreeQueryResultsDistances() to get distances                 |
//+------------------------------------------------------------------+
static int CNearestNeighbor::KDTreeQueryAKNN(CKDTree &kdt,double &x[],
                                             int k,const bool selfmatch,
                                             const double eps)
  {
//--- create variables
   int result=0;
   int i=0;
   int j=0;
//--- check
   if(!CAp::Assert(k>0,__FUNCTION__+": incorrect K!"))
      return(-1);
//--- check
   if(!CAp::Assert(eps>=0.0,__FUNCTION__+": incorrect Eps!"))
      return(-1);
//--- check
   if(!CAp::Assert(CAp::Len(x)>=kdt.m_nx,__FUNCTION__+": Length(X)<NX!"))
      return(-1);
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,kdt.m_nx),__FUNCTION__+": X contains infinite or NaN values!"))
      return(-1);
//--- Prepare parameters
   k=MathMin(k,kdt.m_n);
   kdt.m_kneeded=k;
   kdt.m_rneeded=0;
   kdt.m_selfmatch=selfmatch;
//--- check
   if(kdt.m_normtype==2)
      kdt.m_approxf=1/CMath::Sqr(1+eps);
   else
      kdt.m_approxf=1/(1+eps);
   kdt.m_kcur=0;
//--- calculate distance from point to current bounding box
   KDTreeInitBox(kdt,x);
//--- call recursive search
//--- results are returned as heap
   KDTreeQueryNNRec(kdt,0);
//--- pop from heap to generate ordered representation
//--- last element is non pop'ed because it is already in
//--- its place
   result=kdt.m_kcur;
   j=kdt.m_kcur;
   for(i=kdt.m_kcur;i>=2;i--)
      CTSort::TagHeapPopI(kdt.m_r,kdt.m_idx,j);
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| X-values from last query                                         |
//| INPUT PARAMETERS                                                 |
//|     KDT     -   KD-tree                                          |
//|     X       -   possibly pre-allocated buffer. If X is too small |
//|                 to store result, it is resized. If size(X) is    |
//|                 enough to store result, it is left unchanged.    |
//| OUTPUT PARAMETERS                                                |
//|     X       -   rows are filled with X-values                    |
//| NOTES                                                            |
//| 1. points are ordered by distance from the query point (first =  |
//|    closest)                                                      |
//| 2. if  XY is larger than required to store result, only leading  |
//|    part will be overwritten; trailing part will be left          |
//|    unchanged. So if on input XY = [[A,B],[C,D]], and result is   |
//|    [1,2], then on exit we will get XY = [[1,2],[C,D]]. This is   |
//|    done purposely to increase performance; if you want function  | 
//|    to resize array according to result size, use function with   | 
//| same name and suffix 'I'.                                        |
//| SEE ALSO                                                         |
//| * KDTreeQueryResultsXY()            X- and Y-values              |
//| * KDTreeQueryResultsTags()          tag values                   |
//| * KDTreeQueryResultsDistances()     distances                    |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeQueryResultsX(CKDTree &kdt,CMatrixDouble &x)
  {
//--- create variables
   int i=0;
   int k=0;
   int i_=0;
   int i1_=0;
//--- check
   if(kdt.m_kcur==0)
      return;
//--- check
   if(CAp::Rows(x)<kdt.m_kcur||CAp::Cols(x)<kdt.m_nx)
      x.Resize(kdt.m_kcur,kdt.m_nx);
//--- change value
   k=kdt.m_kcur;
   for(i=0;i<=k-1;i++)
     {
      //--- calculation
      i1_=kdt.m_nx;
      for(i_=0;i_<=kdt.m_nx-1;i_++)
         x[i].Set(i_,kdt.m_xy[kdt.m_idx[i]][i_+i1_]);
     }
  }
//+------------------------------------------------------------------+
//| X- and Y-values from last query                                  |
//| INPUT PARAMETERS                                                 |
//|     KDT     -   KD-tree                                          |
//|     XY      -   possibly pre-allocated buffer. If XY is too small|
//|                 to store result, it is resized. If size(XY) is   |
//|                 enough to store result, it is left unchanged.    |
//| OUTPUT PARAMETERS                                                |
//|     XY      -   rows are filled with points: first NX columns    |
//|                 with X-values, next NY columns - with Y-values.  |
//| NOTES                                                            |
//| 1. points are ordered by distance from the query point (first =  |
//|    closest)                                                      |
//| 2. if  XY is larger than required to store result, only leading  |
//|    part will be overwritten; trailing part will be left          |
//|    unchanged. So if on input XY = [[A,B],[C,D]], and result is   |
//|    [1,2], then on exit we will get XY = [[1,2],[C,D]]. This is   |
//|    done purposely to increase performance; if you want function  |
//|    to resize array according to result size, use function with   |
//|    same name and suffix 'I'.                                     |
//| SEE ALSO                                                         |
//| * KDTreeQueryResultsX()             X-values                     |
//| * KDTreeQueryResultsTags()          tag values                   |
//| * KDTreeQueryResultsDistances()     distances                    |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeQueryResultsXY(CKDTree &kdt,CMatrixDouble &xy)
  {
//--- create variables
   int i=0;
   int k=0;
   int i_=0;
   int i1_=0;
//--- check
   if(kdt.m_kcur==0)
      return;
//--- check
   if(CAp::Rows(xy)<kdt.m_kcur || CAp::Cols(xy)<kdt.m_nx+kdt.m_ny)
      xy.Resize(kdt.m_kcur,kdt.m_nx+kdt.m_ny);
//--- change value
   k=kdt.m_kcur;
   for(i=0;i<=k-1;i++)
     {
      //--- calculation
      i1_=kdt.m_nx;
      for(i_=0;i_<=kdt.m_nx+kdt.m_ny-1;i_++)
         xy[i].Set(i_,kdt.m_xy[kdt.m_idx[i]][i_+i1_]);
     }
  }
//+------------------------------------------------------------------+
//| Tags from last query                                             |
//| INPUT PARAMETERS                                                 |
//|     KDT     -   KD-tree                                          |
//|     Tags    -   possibly pre-allocated buffer. If X is too small |
//|                 to store result, it is resized. If size(X) is    |
//|                 enough to store result, it is left unchanged.    |
//| OUTPUT PARAMETERS                                                |
//|     Tags    -   filled with tags associated with points,         |
//|                 or, when no tags were supplied, with zeros       |
//| NOTES                                                            |
//| 1. points are ordered by distance from the query point (first    |
//|    = closest)                                                    |
//| 2. if  XY is larger than required to store result, only leading  |
//|    part will be overwritten; trailing part will be left          |
//|    unchanged. So if on input XY = [[A,B],[C,D]], and result is   |
//|    [1,2], then on exit we will get XY = [[1,2],[C,D]]. This is   |
//|    done purposely to increase performance; if you want function  |
//|    to resize array according to result size, use function with   |
//|    same name and suffix 'I'.                                     |
//| SEE ALSO                                                         |
//| * KDTreeQueryResultsX()             X-values                     |
//| * KDTreeQueryResultsXY()            X- and Y-values              |
//| * KDTreeQueryResultsDistances()     distances                    |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeQueryResultsTags(CKDTree &kdt,int &tags[])
  {
//--- create variables
   int i=0;
   int k=0;
//--- check
   if(kdt.m_kcur==0)
      return;
//--- check
   if(CAp::Len(tags)<kdt.m_kcur)
     {
      //--- allocation
      ArrayResizeAL(tags,kdt.m_kcur);
     }
//--- calculation
   k=kdt.m_kcur;
   for(i=0;i<=k-1;i++)
      tags[i]=kdt.m_tags[kdt.m_idx[i]];
  }
//+------------------------------------------------------------------+
//| Distances from last query                                        |
//| INPUT PARAMETERS                                                 |
//|     KDT     -   KD-tree                                          |
//|     R       -   possibly pre-allocated buffer. If X is too small |
//|                 to store result, it is resized. If size(X) is    |
//|                 enough to store result, it is left unchanged.    |
//| OUTPUT PARAMETERS                                                |
//|     R       -   filled with distances (in corresponding norm)    |
//| NOTES                                                            |
//| 1. points are ordered by distance from the query point (first    |
//|    = closest)                                                    |
//| 2. if  XY is larger than required to store result, only leading  |
//|    part will be overwritten; trailing part will be left          |
//|    unchanged. So if on input XY = [[A,B],[C,D]], and result is   |
//|    [1,2], then on exit we will get XY = [[1,2],[C,D]]. This is   |
//|    done purposely to increase performance; if you want function  |
//|    to resize array according to result size, use function with   |
//|    same name and suffix 'I'.                                     |
//| SEE ALSO                                                         |
//| * KDTreeQueryResultsX()             X-values                     |
//| * KDTreeQueryResultsXY()            X- and Y-values              |
//| * KDTreeQueryResultsTags()          tag values                   |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeQueryResultsDistances(CKDTree &kdt,
                                                          double &r[])
  {
//--- create variables
   int i=0;
   int k=0;
//--- check
   if(kdt.m_kcur==0)
      return;
//--- check
   if(CAp::Len(r)<kdt.m_kcur)
     {
      //--- allocation
      ArrayResizeAL(r,kdt.m_kcur);
     }
   k=kdt.m_kcur;
//--- unload norms
//--- Abs() call is used to handle cases with negative norms
//--- (generated during KFN requests)
   if(kdt.m_normtype==0)
     {
      for(i=0;i<=k-1;i++)
         r[i]=MathAbs(kdt.m_r[i]);
     }
//--- check
   if(kdt.m_normtype==1)
     {
      for(i=0;i<=k-1;i++)
         r[i]=MathAbs(kdt.m_r[i]);
     }
//--- check
   if(kdt.m_normtype==2)
     {
      for(i=0;i<=k-1;i++)
         r[i]=MathSqrt(MathAbs(kdt.m_r[i]));
     }
  }
//+------------------------------------------------------------------+
//| X-values from last query; 'interactive' variant for languages    |
//| like Python which support constructs like "X =                   |
//| KDTreeQueryResultsXI(KDT)" and interactive mode of interpreter.  |
//| This function allocates new array on each call, so it is         |
//| significantly slower than its 'non-interactive' counterpart, but |
//| it is more convenient when you call it from command line.        |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeQueryResultsXI(CKDTree &kdt,CMatrixDouble &x)
  {
//--- memory reset
   x.Resize(0,0);
//--- function call
   KDTreeQueryResultsX(kdt,x);
  }
//+------------------------------------------------------------------+
//| XY-values from last query; 'interactive' variant for languages   |
//| like Python which support constructs like "XY =                  |
//| KDTreeQueryResultsXYI(KDT)" and interactive mode of interpreter. |
//| This function allocates new array on each call, so it is         |
//| significantly slower than its 'non-interactive' counterpart, but | 
//| it is more convenient when you call it from command line.        |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeQueryResultsXYI(CKDTree &kdt,CMatrixDouble &xy)
  {
//--- memory reset
   xy.Resize(0,0);
//--- function call
   KDTreeQueryResultsXY(kdt,xy);
  }
//+------------------------------------------------------------------+
//| Tags from last query; 'interactive' variant for languages like   |
//| Python which  support  constructs  like "Tags =                  |
//| KDTreeQueryResultsTagsI(KDT)" and interactive mode of            |
//| interpreter.                                                     |
//| This function allocates new array on each call, so it is         |
//| significantly slower than its 'non-interactive' counterpart, but |
//| it is more convenient when you call it from command line.        |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeQueryResultsTagsI(CKDTree &kdt,
                                                      int &tags[])
  {
//--- memory reset
   ArrayResizeAL(tags,0);
//--- function call
   KDTreeQueryResultsTags(kdt,tags);
  }
//+------------------------------------------------------------------+
//| Distances from last query; 'interactive' variant for languages   |
//| like Python which support constructs like "R =                   |
//| KDTreeQueryResultsDistancesI(KDT)" and interactive mode of       |
//| interpreter.                                                     |
//| This function allocates new array on each call, so it is         |
//| significantly slower than its 'non-interactive' counterpart, but | 
//| it is more convenient when you call it from command line.        |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeQueryResultsDistancesI(CKDTree &kdt,
                                                           double &r[])
  {
//--- memory reset
   ArrayResizeAL(r,0);
//--- function call
   KDTreeQueryResultsDistances(kdt,r);
  }
//+------------------------------------------------------------------+
//| Serializer: allocation                                           |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeAlloc(CSerializer &s,CKDTree &tree)
  {
//--- Header
   s.Alloc_Entry();
   s.Alloc_Entry();
//--- Data
   s.Alloc_Entry();
   s.Alloc_Entry();
   s.Alloc_Entry();
   s.Alloc_Entry();
//--- allocation
   CApServ::AllocRealMatrix(s,tree.m_xy,-1,-1);
   CApServ::AllocIntegerArray(s,tree.m_tags,-1);
   CApServ::AllocRealArray(s,tree.m_boxmin,-1);
   CApServ::AllocRealArray(s,tree.m_boxmax,-1);
   CApServ::AllocIntegerArray(s,tree.m_nodes,-1);
   CApServ::AllocRealArray(s,tree.m_splits,-1);
  }
//+------------------------------------------------------------------+
//| Serializer: serialization                                        |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeSerialize(CSerializer &s,CKDTree &tree)
  {
//--- Header
   s.Serialize_Int(CSCodes::GetKDTreeSerializationCode());
   s.Serialize_Int(m_kdtreefirstversion);
//--- Data
   s.Serialize_Int(tree.m_n);
   s.Serialize_Int(tree.m_nx);
   s.Serialize_Int(tree.m_ny);
   s.Serialize_Int(tree.m_normtype);
//--- serialization
   CApServ::SerializeRealMatrix(s,tree.m_xy,-1,-1);
   CApServ::SerializeIntegerArray(s,tree.m_tags,-1);
   CApServ::SerializeRealArray(s,tree.m_boxmin,-1);
   CApServ::SerializeRealArray(s,tree.m_boxmax,-1);
   CApServ::SerializeIntegerArray(s,tree.m_nodes,-1);
   CApServ::SerializeRealArray(s,tree.m_splits,-1);
  }
//+------------------------------------------------------------------+
//| Serializer: unserialization                                      |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeUnserialize(CSerializer &s,CKDTree &tree)
  {
//--- create variables
   int i0=0;
   int i1=0;
//--- check correctness of header
   i0=s.Unserialize_Int();
//--- check
   if(!CAp::Assert(i0==CSCodes::GetKDTreeSerializationCode(),__FUNCTION__+": stream header corrupted"))
      return;
   i1=s.Unserialize_Int();
//--- check
   if(!CAp::Assert(i1==m_kdtreefirstversion,__FUNCTION__+": stream header corrupted"))
      return;
//--- Unserialize data
   tree.m_n=s.Unserialize_Int();
   tree.m_nx=s.Unserialize_Int();
   tree.m_ny=s.Unserialize_Int();
   tree.m_normtype=s.Unserialize_Int();
//--- unserializetion
   CApServ::UnserializeRealMatrix(s,tree.m_xy);
   CApServ::UnserializeIntegerArray(s,tree.m_tags);
   CApServ::UnserializeRealArray(s,tree.m_boxmin);
   CApServ::UnserializeRealArray(s,tree.m_boxmax);
   CApServ::UnserializeIntegerArray(s,tree.m_nodes);
   CApServ::UnserializeRealArray(s,tree.m_splits);
//--- function call
   KDTreeAllocTemporaries(tree,tree.m_n,tree.m_nx,tree.m_ny);
  }
//+------------------------------------------------------------------+
//| Rearranges nodes [I1,I2) using partition in D-th dimension with  |
//| S as threshold. Returns split position I3: [I1,I3) and [I3,I2)   |
//| are created as result.                                           |
//| This subroutine doesn't create tree structures, just rearranges  |
//| nodes.                                                           |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeSplit(CKDTree &kdt,const int i1,
                                          const int i2,const int d,
                                          const double s,int &i3)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    ileft=0;
   int    iright=0;
   double v=0;
//--- initialization
   i3=0;
//--- split XY/Tags in two parts:
//--- * [ILeft,IRight] is non-processed part of XY/Tags
//--- After cycle is done, we have Ileft=IRight. We deal with
//--- this element separately.
//--- After this, [I1,ILeft) contains left part, and [ILeft,I2)
//--- contains right part.
   ileft=i1;
   iright=i2-1;
   while(ileft<iright)
     {
      //--- check
      if(kdt.m_xy[ileft][d]<=s)
        {
         //--- XY[ILeft] is on its place.
         //--- Advance ILeft.
         ileft=ileft+1;
        }
      else
        {
         //--- XY[ILeft,..] must be at IRight.
         //--- Swap and advance IRight.
         for(i=0;i<=2*kdt.m_nx+kdt.m_ny-1;i++)
           {
            v=kdt.m_xy[ileft][i];
            kdt.m_xy[ileft].Set(i,kdt.m_xy[iright][i]);
            kdt.m_xy[iright].Set(i,v);
           }
         //--- change values
         j=kdt.m_tags[ileft];
         kdt.m_tags[ileft]=kdt.m_tags[iright];
         kdt.m_tags[iright]=j;
         iright=iright-1;
        }
     }
//--- check
   if(kdt.m_xy[ileft][d]<=s)
      ileft=ileft+1;
   else
      iright=iright-1;
//--- get result
   i3=ileft;
  }
//+------------------------------------------------------------------+
//| Recursive kd-tree generation subroutine.                         |
//| PARAMETERS                                                       |
//|     KDT         tree                                             |
//|     NodesOffs   unused part of Nodes[] which must be filled by   |
//|                 tree                                             |
//|     SplitsOffs  unused part of Splits[]                          |
//|     I1, I2      points from [I1,I2) are processed                |
//| NodesOffs[] and SplitsOffs[] must be large enough.               |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeGenerateTreeRec(CKDTree &kdt,int &nodesoffs,
                                                    int &splitsoffs,const int i1,
                                                    const int i2,const int maxleafsize)
  {
//--- create variables
   int    n=0;
   int    nx=0;
   int    ny=0;
   int    i=0;
   int    j=0;
   int    oldoffs=0;
   int    i3=0;
   int    cntless=0;
   int    cntgreater=0;
   double minv=0;
   double maxv=0;
   int    minidx=0;
   int    maxidx=0;
   int    d=0;
   double ds=0;
   double s=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(!CAp::Assert(i2>i1,__FUNCTION__+": internal error"))
      return;
//--- Generate leaf if needed
   if(i2-i1<=maxleafsize)
     {
      kdt.m_nodes[nodesoffs+0]=i2-i1;
      kdt.m_nodes[nodesoffs+1]=i1;
      nodesoffs=nodesoffs+2;
      //--- exit the function
      return;
     }
//--- Load values for easier access
   nx=kdt.m_nx;
   ny=kdt.m_ny;
//--- select dimension to split:
//--- * D is a dimension number
   d=0;
   ds=kdt.m_curboxmax[0]-kdt.m_curboxmin[0];
   for(i=1;i<=nx-1;i++)
     {
      v=kdt.m_curboxmax[i]-kdt.m_curboxmin[i];
      //--- check
      if(v>ds)
        {
         ds=v;
         d=i;
        }
     }
//--- Select split position S using sliding midpoint rule,
//--- rearrange points into [I1,I3) and [I3,I2)
   s=kdt.m_curboxmin[d]+0.5*ds;
   i1_=(i1) -(0);
   for(i_=0;i_<=i2-i1-1;i_++)
      kdt.m_buf[i_]=kdt.m_xy[i_+i1_][d];
//--- change values
   n=i2-i1;
   cntless=0;
   cntgreater=0;
   minv=kdt.m_buf[0];
   maxv=kdt.m_buf[0];
   minidx=i1;
   maxidx=i1;
   for(i=0;i<=n-1;i++)
     {
      v=kdt.m_buf[i];
      //--- check
      if(v<minv)
        {
         minv=v;
         minidx=i1+i;
        }
      //--- check
      if(v>maxv)
        {
         maxv=v;
         maxidx=i1+i;
        }
      //--- check
      if(v<s)
         cntless=cntless+1;
      //--- check
      if(v>s)
         cntgreater=cntgreater+1;
     }
//--- check
   if(cntless>0&&cntgreater>0)
     {
      //--- normal midpoint split
      KDTreeSplit(kdt,i1,i2,d,s,i3);
     }
   else
     {
      //--- sliding midpoint
      if(cntless==0)
        {
         //--- 1. move split to MinV,
         //--- 2. place one point to the left bin (move to I1),
         //---    others - to the right bin
         s=minv;
         //--- check
         if(minidx!=i1)
           {
            for(i=0;i<=2*kdt.m_nx+kdt.m_ny-1;i++)
              {
               v=kdt.m_xy[minidx][i];
               kdt.m_xy[minidx].Set(i,kdt.m_xy[i1][i]);
               kdt.m_xy[i1].Set(i,v);
              }
            //--- change values
            j=kdt.m_tags[minidx];
            kdt.m_tags[minidx]=kdt.m_tags[i1];
            kdt.m_tags[i1]=j;
           }
         i3=i1+1;
        }
      else
        {
         //--- 1. move split to MaxV,
         //--- 2. place one point to the right bin (move to I2-1),
         //---    others - to the left bin
         s=maxv;
         //--- check
         if(maxidx!=i2-1)
           {
            for(i=0;i<=2*kdt.m_nx+kdt.m_ny-1;i++)
              {
               v=kdt.m_xy[maxidx][i];
               kdt.m_xy[maxidx].Set(i,kdt.m_xy[i2-1][i]);
               kdt.m_xy[i2-1].Set(i,v);
              }
            //--- change values
            j=kdt.m_tags[maxidx];
            kdt.m_tags[maxidx]=kdt.m_tags[i2-1];
            kdt.m_tags[i2-1]=j;
           }
         i3=i2-1;
        }
     }
//--- Generate 'split' node
   kdt.m_nodes[nodesoffs+0]=0;
   kdt.m_nodes[nodesoffs+1]=d;
   kdt.m_nodes[nodesoffs+2]=splitsoffs;
   kdt.m_splits[splitsoffs+0]=s;
   oldoffs=nodesoffs;
   nodesoffs=nodesoffs+m_splitnodesize;
   splitsoffs=splitsoffs+1;
//--- Recirsive generation:
//--- * update CurBox
//--- * call subroutine
//--- * restore CurBox
   kdt.m_nodes[oldoffs+3]=nodesoffs;
   v=kdt.m_curboxmax[d];
   kdt.m_curboxmax[d]=s;
//--- function call
   KDTreeGenerateTreeRec(kdt,nodesoffs,splitsoffs,i1,i3,maxleafsize);
   kdt.m_curboxmax[d]=v;
   kdt.m_nodes[oldoffs+4]=nodesoffs;
   v=kdt.m_curboxmin[d];
   kdt.m_curboxmin[d]=s;
//--- function call
   KDTreeGenerateTreeRec(kdt,nodesoffs,splitsoffs,i3,i2,maxleafsize);
   kdt.m_curboxmin[d]=v;
  }
//+------------------------------------------------------------------+
//| Recursive subroutine for NN queries.                             |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeQueryNNRec(CKDTree &kdt,const int offs)
  {
//--- create variables
   double ptdist=0;
   int    i=0;
   int    j=0;
   int    nx=0;
   int    i1=0;
   int    i2=0;
   int    d=0;
   double s=0;
   double v=0;
   double t1=0;
   int    childbestoffs=0;
   int    childworstoffs=0;
   int    childoffs=0;
   double prevdist=0;
   bool   todive;
   bool   bestisleft;
   bool   updatemin;
//--- Leaf node.
//--- Process points.
   if(kdt.m_nodes[offs]>0)
     {
      i1=kdt.m_nodes[offs+1];
      i2=i1+kdt.m_nodes[offs];
      for(i=i1;i<=i2-1;i++)
        {
         //--- Calculate distance
         ptdist=0;
         nx=kdt.m_nx;
         //--- check
         if(kdt.m_normtype==0)
           {
            for(j=0;j<=nx-1;j++)
               ptdist=MathMax(ptdist,MathAbs(kdt.m_xy[i][j]-kdt.m_x[j]));
           }
         //--- check
         if(kdt.m_normtype==1)
           {
            for(j=0;j<=nx-1;j++)
               ptdist=ptdist+MathAbs(kdt.m_xy[i][j]-kdt.m_x[j]);
           }
         //--- check
         if(kdt.m_normtype==2)
           {
            for(j=0;j<=nx-1;j++)
               ptdist=ptdist+CMath::Sqr(kdt.m_xy[i][j]-kdt.m_x[j]);
           }
         //--- Skip points with zero distance if self-matches are turned off
         if(ptdist==0.0 && !kdt.m_selfmatch)
            continue;
         //--- We CAN'T process point if R-criterion isn't satisfied,
         //--- i.e. (RNeeded<>0) AND (PtDist>R).
         if(kdt.m_rneeded==0.0 || ptdist<=kdt.m_rneeded)
           {
            //--- R-criterion is satisfied, we must either:
            //--- * replace worst point, if (KNeeded<>0) AND (KCur=KNeeded)
            //---   (or skip, if worst point is better)
            //--- * add point without replacement otherwise
            if(kdt.m_kcur<kdt.m_kneeded || kdt.m_kneeded==0)
              {
               //--- add current point to heap without replacement
               CTSort::TagHeapPushI(kdt.m_r,kdt.m_idx,kdt.m_kcur,ptdist,i);
              }
            else
              {
               //--- New points are added or not, depending on their distance.
               //--- If added, they replace element at the top of the heap
               if(ptdist<(double)(kdt.m_r[0]))
                 {
                  //--- check
                  if(kdt.m_kneeded==1)
                    {
                     kdt.m_idx[0]=i;
                     kdt.m_r[0]=ptdist;
                    }
                  else
                     CTSort::TagHeapReplaceTopI(kdt.m_r,kdt.m_idx,kdt.m_kneeded,ptdist,i);
                 }
              }
           }
        }
      //--- exit the function
      return;
     }
//--- Simple split
   if(kdt.m_nodes[offs]==0)
     {
      //--- Load:
      //--- * D  dimension to split
      //--- * S  split position
      d=kdt.m_nodes[offs+1];
      s=kdt.m_splits[kdt.m_nodes[offs+2]];
      //--- Calculate:
      //--- * ChildBestOffs      child box with best chances
      //--- * ChildWorstOffs     child box with worst chances
      if(kdt.m_x[d]<=s)
        {
         childbestoffs=kdt.m_nodes[offs+3];
         childworstoffs=kdt.m_nodes[offs+4];
         bestisleft=true;
        }
      else
        {
         childbestoffs=kdt.m_nodes[offs+4];
         childworstoffs=kdt.m_nodes[offs+3];
         bestisleft=false;
        }
      //--- Navigate through childs
      for(i=0;i<=1;i++)
        {
         //--- Select child to process:
         //--- * ChildOffs      current child offset in Nodes[]
         //--- * UpdateMin      whether minimum or maximum value
         //---                  of bounding box is changed on update
         if(i==0)
           {
            childoffs=childbestoffs;
            updatemin=!bestisleft;
           }
         else
           {
            updatemin=bestisleft;
            childoffs=childworstoffs;
           }
         //--- Update bounding box and current distance
         if(updatemin)
           {
            prevdist=kdt.m_curdist;
            t1=kdt.m_x[d];
            v=kdt.m_curboxmin[d];
            //--- check
            if(t1<=s)
              {
               //--- check
               if(kdt.m_normtype==0)
                  kdt.m_curdist=MathMax(kdt.m_curdist,s-t1);
               //--- check
               if(kdt.m_normtype==1)
                  kdt.m_curdist=kdt.m_curdist-MathMax(v-t1,0)+s-t1;
               //--- check
               if(kdt.m_normtype==2)
                  kdt.m_curdist=kdt.m_curdist-CMath::Sqr(MathMax(v-t1,0))+CMath::Sqr(s-t1);
              }
            kdt.m_curboxmin[d]=s;
           }
         else
           {
            prevdist=kdt.m_curdist;
            t1=kdt.m_x[d];
            v=kdt.m_curboxmax[d];
            //--- check
            if(t1>=s)
              {
               //--- check
               if(kdt.m_normtype==0)
                  kdt.m_curdist=MathMax(kdt.m_curdist,t1-s);
               //--- check
               if(kdt.m_normtype==1)
                  kdt.m_curdist=kdt.m_curdist-MathMax(t1-v,0)+t1-s;
               //--- check
               if(kdt.m_normtype==2)
                  kdt.m_curdist=kdt.m_curdist-CMath::Sqr(MathMax(t1-v,0))+CMath::Sqr(t1-s);
              }
            kdt.m_curboxmax[d]=s;
           }
         //--- Decide: to dive into cell or not to dive
         if(kdt.m_rneeded!=0.0 && kdt.m_curdist>kdt.m_rneeded)
            todive=false;
         else
           {
            //--- check
            if(kdt.m_kcur<kdt.m_kneeded || kdt.m_kneeded==0)
              {
               //--- KCur<KNeeded (i.e. not all points are found)
               todive=true;
              }
            else
              {
               //--- KCur=KNeeded,decide to dive or not to dive
               //--- using point position relative to bounding box.
               todive=kdt.m_curdist<=(double)(kdt.m_r[0]*kdt.m_approxf);
              }
           }
         //--- check
         if(todive)
            KDTreeQueryNNRec(kdt,childoffs);
         //--- Restore bounding box and distance
         if(updatemin)
            kdt.m_curboxmin[d]=v;
         else
            kdt.m_curboxmax[d]=v;
         kdt.m_curdist=prevdist;
        }
      //--- exit the function
      return;
     }
  }
//+------------------------------------------------------------------+
//| Copies X[] to KDT.X[]                                            |
//| Loads distance from X[] to bounding box.                         |
//| Initializes CurBox[].                                            |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeInitBox(CKDTree &kdt,double &x[])
  {
//--- create variables
   int    i=0;
   double vx=0;
   double vmin=0;
   double vmax=0;
//--- calculate distance from point to current bounding box
   kdt.m_curdist=0;
//--- check
   if(kdt.m_normtype==0)
     {
      for(i=0;i<=kdt.m_nx-1;i++)
        {
         vx=x[i];
         vmin=kdt.m_boxmin[i];
         vmax=kdt.m_boxmax[i];
         kdt.m_x[i]=vx;
         kdt.m_curboxmin[i]=vmin;
         kdt.m_curboxmax[i]=vmax;
         //--- check
         if(vx<vmin)
            kdt.m_curdist=MathMax(kdt.m_curdist,vmin-vx);
         else
           {
            //--- check
            if(vx>vmax)
               kdt.m_curdist=MathMax(kdt.m_curdist,vx-vmax);
           }
        }
     }
//--- check
   if(kdt.m_normtype==1)
     {
      for(i=0;i<=kdt.m_nx-1;i++)
        {
         vx=x[i];
         vmin=kdt.m_boxmin[i];
         vmax=kdt.m_boxmax[i];
         kdt.m_x[i]=vx;
         kdt.m_curboxmin[i]=vmin;
         kdt.m_curboxmax[i]=vmax;
         //--- check
         if(vx<vmin)
            kdt.m_curdist=kdt.m_curdist+vmin-vx;
         else
           {
            //--- check
            if(vx>vmax)
               kdt.m_curdist=kdt.m_curdist+vx-vmax;
           }
        }
     }
//--- check
   if(kdt.m_normtype==2)
     {
      for(i=0;i<=kdt.m_nx-1;i++)
        {
         vx=x[i];
         vmin=kdt.m_boxmin[i];
         vmax=kdt.m_boxmax[i];
         kdt.m_x[i]=vx;
         kdt.m_curboxmin[i]=vmin;
         kdt.m_curboxmax[i]=vmax;
         //--- check
         if(vx<vmin)
            kdt.m_curdist=kdt.m_curdist+CMath::Sqr(vmin-vx);
         else
           {
            //--- check
            if(vx>vmax)
               kdt.m_curdist=kdt.m_curdist+CMath::Sqr(vx-vmax);
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| This function allocates all dataset-independent array fields of  | 
//| KDTree, i.e. such array fields that their dimensions do not      |
//| depend on dataset size.                                          |
//| This function do not sets KDT.NX or KDT.NY - it just allocates   |
//| arrays                                                           |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeAllocDataSetIndependent(CKDTree &kdt,
                                                            const int nx,
                                                            const int ny)
  {
//--- allocation
   ArrayResizeAL(kdt.m_x,nx);
   ArrayResizeAL(kdt.m_boxmin,nx);
   ArrayResizeAL(kdt.m_boxmax,nx);
   ArrayResizeAL(kdt.m_curboxmin,nx);
   ArrayResizeAL(kdt.m_curboxmax,nx);
  }
//+------------------------------------------------------------------+
//| This function allocates all dataset-dependent array fields of    |
//| KDTree, i.e. such array fields that their dimensions depend on   |
//| dataset size.                                                    |
//| This function do not sets KDT.N, KDT.NX or KDT.NY -              |
//| it just allocates arrays.                                        |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeAllocDataSetDependent(CKDTree &kdt,
                                                          const int n,
                                                          const int nx,
                                                          const int ny)
  {
//--- allocation
   kdt.m_xy.Resize(n,2*nx+ny);
   ArrayResizeAL(kdt.m_tags,n);
   ArrayResizeAL(kdt.m_idx,n);
   ArrayResizeAL(kdt.m_r,n);
   ArrayResizeAL(kdt.m_x,nx);
   ArrayResizeAL(kdt.m_buf,MathMax(n,nx));
   ArrayResizeAL(kdt.m_nodes,m_splitnodesize*2*n);
   ArrayResizeAL(kdt.m_splits,2*n);
  }
//+------------------------------------------------------------------+
//| This function allocates temporaries.                             |
//| This function do not sets KDT.N,KDT.NX or KDT.NY -               |
//| it just allocates arrays.                                        |
//+------------------------------------------------------------------+
static void CNearestNeighbor::KDTreeAllocTemporaries(CKDTree &kdt,const int n,
                                                     const int nx,const int ny)
  {
//--- allocation
   ArrayResizeAL(kdt.m_x,nx);
   ArrayResizeAL(kdt.m_idx,n);
   ArrayResizeAL(kdt.m_r,n);
   ArrayResizeAL(kdt.m_buf,MathMax(n,nx));
   ArrayResizeAL(kdt.m_curboxmin,nx);
   ArrayResizeAL(kdt.m_curboxmax,nx);
  }
//+------------------------------------------------------------------+
