
// OpenMP will be compiled in if this flag is set to 1 AND the compiler being
// used supports it (i.e. the _OPENMP symbol is defined)
#define USE_OMP 1

#ifdef USE_MPI
#include <mpi.h>

/*
   define one of these three symbols:

   SEDOV_SYNC_POS_VEL_NONE
   SEDOV_SYNC_POS_VEL_EARLY
   SEDOV_SYNC_POS_VEL_LATE
*/

#define SEDOV_SYNC_POS_VEL_EARLY 1
#endif

#include <math.h>
#include <vector>

#include "meshapi/RangeSet.hpp"
#include "meshapi/StaticConstantRelation.hpp"
#include "meshapi/StaticVariableRelation.hpp"
#include "meshapi/DynamicVariableRelation.hpp"
#include "meshapi/Map.hpp"

#include "sidre/sidre.hpp"

//**************************************************
// Allow flexibility for arithmetic representations 
//**************************************************

#define MAX(a, b) ( ((a) > (b)) ? (a) : (b))


// Precision specification
typedef float        real4 ;
typedef double       real8 ;
typedef long double  real10 ;  // 10 bytes on x86

typedef int    Index_t ; // array subscript and loop index
typedef real8  Real_t ;  // floating point representation
typedef int    Int_t ;   // integer representation

enum { VolumeError = -1, QStopError = -2 } ;

inline real4  SQRT(real4  arg) { return sqrtf(arg) ; }
inline real8  SQRT(real8  arg) { return sqrt(arg) ; }
inline real10 SQRT(real10 arg) { return sqrtl(arg) ; }

inline real4  CBRT(real4  arg) { return cbrtf(arg) ; }
inline real8  CBRT(real8  arg) { return cbrt(arg) ; }
inline real10 CBRT(real10 arg) { return cbrtl(arg) ; }

inline real4  FABS(real4  arg) { return fabsf(arg) ; }
inline real8  FABS(real8  arg) { return fabs(arg) ; }
inline real10 FABS(real10 arg) { return fabsl(arg) ; }


// Stuff needed for boundary conditions
// 2 BCs on each of 6 hexahedral faces (12 bits)
#define XI_M        0x00007
#define XI_M_SYMM   0x00001
#define XI_M_FREE   0x00002
#define XI_M_COMM   0x00004

#define XI_P        0x00038
#define XI_P_SYMM   0x00008
#define XI_P_FREE   0x00010
#define XI_P_COMM   0x00020

#define ETA_M       0x001c0
#define ETA_M_SYMM  0x00040
#define ETA_M_FREE  0x00080
#define ETA_M_COMM  0x00100

#define ETA_P       0x00e00
#define ETA_P_SYMM  0x00200
#define ETA_P_FREE  0x00400
#define ETA_P_COMM  0x00800

#define ZETA_M      0x07000
#define ZETA_M_SYMM 0x01000
#define ZETA_M_FREE 0x02000
#define ZETA_M_COMM 0x04000

#define ZETA_P      0x38000
#define ZETA_P_SYMM 0x08000
#define ZETA_P_FREE 0x10000
#define ZETA_P_COMM 0x20000

// MPI Message Tags
#define MSG_COMM_SBN      1024
#define MSG_SYNC_POS_VEL  2048
#define MSG_MONOQ         3072

#define MAX_FIELDS_PER_MPI_COMM 6

// Assume 128 byte coherence
// Assume Real_t is an "integral power of 2" bytes wide
#define CACHE_COHERENCE_PAD_REAL (128 / sizeof(Real_t))

#define CACHE_ALIGN_REAL(n) \
   (((n) + (CACHE_COHERENCE_PAD_REAL - 1)) & ~(CACHE_COHERENCE_PAD_REAL-1))

//////////////////////////////////////////////////////
// Primary data structure
//////////////////////////////////////////////////////

/*
 * The implementation of the data abstraction used for lulesh
 * resides entirely in the Domain class below.  You can change
 * grouping and interleaving of fields here to maximize data layout
 * efficiency for your underlying architecture or compiler.
 *
 * For example, fields can be implemented as STL objects or
 * raw array pointers.  As another example, individual fields
 * m_x, m_y, m_z could be bundled into
 *
 *    struct { Real_t x, y, z ; } *m_coord ;
 *
 * allowing accessor functions such as
 *
 *  "Real_t &x(Index_t idx) { return m_coord[idx].x ; }"
 *  "Real_t &y(Index_t idx) { return m_coord[idx].y ; }"
 *  "Real_t &z(Index_t idx) { return m_coord[idx].z ; }"
 */

#define USE_SIDRE 1

class Domain {


   public:

#if USE_SIDRE==1
  typedef asctoolkit::sidre::DataView * luleshRealData;
  typedef asctoolkit::sidre::DataView * luleshIntData;
  typedef asctoolkit::sidre::DataView * luleshIndexData;
#else
  typedef std::vector<Real_t> luleshRealData;
  typedef std::vector<Int_t> luleshIntData;
  typedef std::vector<Index_t> luleshIndexData;
#endif

   typedef asctoolkit::meshapi::RangeSet               ElemSet;
   typedef asctoolkit::meshapi::RangeSet               NodeSet;
   typedef asctoolkit::meshapi::StaticConstantRelation ElemToNodeRelation;
   typedef asctoolkit::meshapi::StaticConstantRelation ElemFaceAdjacencyRelation;

   typedef asctoolkit::meshapi::RangeSet               RegionSet;
   typedef asctoolkit::meshapi::StaticVariableRelation RegionToElemRelation;


   typedef asctoolkit::meshapi::Map<Index_t>           ElemIndexMap;
   typedef asctoolkit::meshapi::Map<Int_t>             ElemIntMap;
   //typedef asctoolkit::meshapi::Map<Real_t>            ElemRealMap;

   //typedef asctoolkit::meshapi::Map<Index_t>           NodeIndexMap;
   //typedef asctoolkit::meshapi::Map<Int_t>             NodeIntMap;
   //typedef asctoolkit::meshapi::Map<Real_t>            NodeRealMap;

   //typedef asctoolkit::meshapi::Map<Index_t>           RegionIndexMap;
   typedef asctoolkit::meshapi::Map<Int_t>             RegionIntMap;
   //typedef asctoolkit::meshapi::Map<Real_t>            RegionRealMap;



   public:

   // Constructor
   Domain(Int_t numRanks, Index_t colLoc,
          Index_t rowLoc, Index_t planeLoc,
          Index_t nx, Int_t tp, Int_t nr, Int_t balance, Int_t cost);

   //
   // ALLOCATION
   //

   void AllocateNodePersistent(Int_t numNode) // Node-centered
   {

#if USE_SIDRE==1
     m_x = m_DataGroup->createViewAndBuffer("m_x",asctoolkit::sidre::DataType::float64(numNode));
     m_y = m_DataGroup->createViewAndBuffer("m_y",asctoolkit::sidre::DataType::float64(numNode));
     m_z = m_DataGroup->createViewAndBuffer("m_z",asctoolkit::sidre::DataType::float64(numNode));

     m_xd = m_DataGroup->createViewAndBuffer("m_xd",asctoolkit::sidre::DataType::float64(numNode));
     m_yd = m_DataGroup->createViewAndBuffer("m_yd",asctoolkit::sidre::DataType::float64(numNode));
     m_zd = m_DataGroup->createViewAndBuffer("m_zd",asctoolkit::sidre::DataType::float64(numNode));

     m_xdd = m_DataGroup->createViewAndBuffer("m_xdd",asctoolkit::sidre::DataType::float64(numNode));
     m_ydd = m_DataGroup->createViewAndBuffer("m_ydd",asctoolkit::sidre::DataType::float64(numNode));
     m_zdd = m_DataGroup->createViewAndBuffer("m_zdd",asctoolkit::sidre::DataType::float64(numNode));

     m_fx = m_DataGroup->createViewAndBuffer("m_fx",asctoolkit::sidre::DataType::float64(numNode));
     m_fy = m_DataGroup->createViewAndBuffer("m_fy",asctoolkit::sidre::DataType::float64(numNode));
     m_fz = m_DataGroup->createViewAndBuffer("m_fz",asctoolkit::sidre::DataType::float64(numNode));

     m_nodalMass = m_DataGroup->createViewAndBuffer("m_nodalMass",asctoolkit::sidre::DataType::float64(numNode));

#else
      m_x.resize(numNode);  // coordinates
      m_y.resize(numNode);
      m_z.resize(numNode);

      m_xd.resize(numNode); // velocities
      m_yd.resize(numNode);
      m_zd.resize(numNode);

      m_xdd.resize(numNode); // accelerations
      m_ydd.resize(numNode);
      m_zdd.resize(numNode);

      m_fx.resize(numNode);  // forces
      m_fy.resize(numNode);
      m_fz.resize(numNode);

      m_nodalMass.resize(numNode);  // mass
#endif
   }

   void AllocateElemPersistent(Int_t numElem) // Elem-centered
   {
       // elem to node incidence relation
      m_nodelist = ElemToNodeRelation(&m_elemSet, &m_nodeSet);

      // elem adjacencies through face
      m_lxim   = ElemFaceAdjacencyRelation(&m_elemSet,&m_elemSet);
      m_lxip   = ElemFaceAdjacencyRelation(&m_elemSet,&m_elemSet);
      m_letam  = ElemFaceAdjacencyRelation(&m_elemSet,&m_elemSet);
      m_letap  = ElemFaceAdjacencyRelation(&m_elemSet,&m_elemSet);
      m_lzetam = ElemFaceAdjacencyRelation(&m_elemSet,&m_elemSet);
      m_lzetap = ElemFaceAdjacencyRelation(&m_elemSet,&m_elemSet);

#if USE_SIDRE==1
      m_elemBC = m_DataGroup->createViewAndBuffer("m_elemBC",asctoolkit::sidre::DataType::int32(numElem));

      m_e = m_DataGroup->createViewAndBuffer("m_e",asctoolkit::sidre::DataType::float64(numElem));
      m_p = m_DataGroup->createViewAndBuffer("m_p",asctoolkit::sidre::DataType::float64(numElem));

      m_q = m_DataGroup->createViewAndBuffer("m_q",asctoolkit::sidre::DataType::float64(numElem));
      m_ql = m_DataGroup->createViewAndBuffer("m_ql",asctoolkit::sidre::DataType::float64(numElem));
      m_qq = m_DataGroup->createViewAndBuffer("m_qq",asctoolkit::sidre::DataType::float64(numElem));

      m_v = m_DataGroup->createViewAndBuffer("m_v",asctoolkit::sidre::DataType::float64(numElem));

      m_volo = m_DataGroup->createViewAndBuffer("m_volo",asctoolkit::sidre::DataType::float64(numElem));
      m_delv = m_DataGroup->createViewAndBuffer("m_delv",asctoolkit::sidre::DataType::float64(numElem));
      m_vdov = m_DataGroup->createViewAndBuffer("m_vdov",asctoolkit::sidre::DataType::float64(numElem));

      m_arealg = m_DataGroup->createViewAndBuffer("m_arealg",asctoolkit::sidre::DataType::float64(numElem));

      m_ss = m_DataGroup->createViewAndBuffer("m_ss",asctoolkit::sidre::DataType::float64(numElem));

      m_elemMass = m_DataGroup->createViewAndBuffer("m_elemMass",asctoolkit::sidre::DataType::float64(numElem));

#else
      m_elemBC.resize(numElem);

      m_e.resize(numElem);
      m_p.resize(numElem);

      m_q.resize(numElem);
      m_ql.resize(numElem);
      m_qq.resize(numElem);

      m_v.resize(numElem);

      m_volo.resize(numElem);
      m_delv.resize(numElem);
      m_vdov.resize(numElem);

      m_arealg.resize(numElem);

      m_ss.resize(numElem);

      m_elemMass.resize(numElem);
#endif
   }

   void AllocateGradients(Int_t numElem, Int_t allElem)
   {
#if USE_SIDRE==1
     // Position gradients
     m_delx_xi = m_DataGroup->createViewAndBuffer("m_delx_xi",asctoolkit::sidre::DataType::float64(numElem));
     m_delx_eta = m_DataGroup->createViewAndBuffer("m_delx_eta",asctoolkit::sidre::DataType::float64(numElem));
     m_delx_zeta = m_DataGroup->createViewAndBuffer("m_delx_zeta",asctoolkit::sidre::DataType::float64(numElem));

     // Velocity gradients
     m_delv_xi = m_DataGroup->createViewAndBuffer("m_delv_xi",asctoolkit::sidre::DataType::float64(allElem));
     m_delv_eta = m_DataGroup->createViewAndBuffer("m_delv_eta",asctoolkit::sidre::DataType::float64(allElem));
     m_delv_zeta = m_DataGroup->createViewAndBuffer("m_delv_zeta",asctoolkit::sidre::DataType::float64(allElem));
#else
      // Position gradients
      m_delx_xi.resize(numElem);
      m_delx_eta.resize(numElem);
      m_delx_zeta.resize(numElem);

      // Velocity gradients
      m_delv_xi.resize(allElem);
      m_delv_eta.resize(allElem);
      m_delv_zeta.resize(allElem);
#endif
   }

   void DeallocateGradients()
   {
#if USE_SIDRE==1
#else
      m_delx_zeta.clear() ;
      m_delx_eta.clear() ;
      m_delx_xi.clear() ;

      m_delv_zeta.clear() ;
      m_delv_eta.clear() ;
      m_delv_xi.clear() ;
#endif
   }

   void AllocateStrains(Int_t numElem)
   {
#if USE_SIDRE==1
     m_dxx = m_DataGroup->createViewAndBuffer("m_dxx",asctoolkit::sidre::DataType::float64(numElem)) ;
     m_dyy = m_DataGroup->createViewAndBuffer("m_dyy",asctoolkit::sidre::DataType::float64(numElem)) ;
     m_dzz = m_DataGroup->createViewAndBuffer("m_dzz",asctoolkit::sidre::DataType::float64(numElem)) ;
#else
      m_dxx.resize(numElem) ;
      m_dyy.resize(numElem) ;
      m_dzz.resize(numElem) ;
#endif
   }

   void DeallocateStrains()
   {
#if USE_SIDRE==1
#else
      m_dzz.clear() ;
      m_dyy.clear() ;
      m_dxx.clear() ;
#endif
   }
   
   //
   // ACCESSORS
   //

   // Node-centered

#if USE_SIDRE==1
   // Nodal coordinates
   Real_t& x(Index_t idx)    { return static_cast<Real_t*>(m_x->getData())[idx] ; }
   Real_t& y(Index_t idx)    { return static_cast<Real_t*>(m_y->getData())[idx] ; }
   Real_t& z(Index_t idx)    { return static_cast<Real_t*>(m_z->getData())[idx] ; }

   // Nodal velocities
   Real_t& xd(Index_t idx)   { return static_cast<Real_t*>(m_xd->getData())[idx] ; }
   Real_t& yd(Index_t idx)   { return static_cast<Real_t*>(m_yd->getData())[idx] ; }
   Real_t& zd(Index_t idx)   { return static_cast<Real_t*>(m_zd->getData())[idx] ; }

   // Nodal accelerations
   Real_t& xdd(Index_t idx)  { return static_cast<Real_t*>(m_xdd->getData())[idx] ; }
   Real_t& ydd(Index_t idx)  { return static_cast<Real_t*>(m_ydd->getData())[idx] ; }
   Real_t& zdd(Index_t idx)  { return static_cast<Real_t*>(m_zdd->getData())[idx] ; }

   // Nodal forces
   Real_t& fx(Index_t idx)   { return static_cast<Real_t*>(m_fx->getData())[idx] ; }
   Real_t& fy(Index_t idx)   { return static_cast<Real_t*>(m_fy->getData())[idx] ; }
   Real_t& fz(Index_t idx)   { return static_cast<Real_t*>(m_fz->getData())[idx] ; }

   // Nodal mass
   Real_t& nodalMass(Index_t idx) { return static_cast<Real_t*>(m_nodalMass->getData())[idx] ; }

   // Nodes on symmetry planes
   Index_t symmX(Index_t idx) { return static_cast<Index_t*>(m_symmX->getData())[idx] ; }
   Index_t symmY(Index_t idx) { return static_cast<Index_t*>(m_symmY->getData())[idx] ; }
   Index_t symmZ(Index_t idx) { return static_cast<Index_t*>(m_symmZ->getData())[idx] ; }
   bool symmXempty()          { return m_symmX->getNode().dtype().number_of_elements() ; }
   bool symmYempty()          { return m_symmY->getNode().dtype().number_of_elements() ; }
   bool symmZempty()          { return m_symmZ->getNode().dtype().number_of_elements() ; }
#else
   // Nodal coordinates
   Real_t& x(Index_t idx)    { return m_x[idx] ; }
   Real_t& y(Index_t idx)    { return m_y[idx] ; }
   Real_t& z(Index_t idx)    { return m_z[idx] ; }

   // Nodal velocities
   Real_t& xd(Index_t idx)   { return m_xd[idx] ; }
   Real_t& yd(Index_t idx)   { return m_yd[idx] ; }
   Real_t& zd(Index_t idx)   { return m_zd[idx] ; }

   // Nodal accelerations
   Real_t& xdd(Index_t idx)  { return m_xdd[idx] ; }
   Real_t& ydd(Index_t idx)  { return m_ydd[idx] ; }
   Real_t& zdd(Index_t idx)  { return m_zdd[idx] ; }

   // Nodal forces
   Real_t& fx(Index_t idx)   { return m_fx[idx] ; }
   Real_t& fy(Index_t idx)   { return m_fy[idx] ; }
   Real_t& fz(Index_t idx)   { return m_fz[idx] ; }

   // Nodal mass
   Real_t& nodalMass(Index_t idx) { return m_nodalMass[idx] ; }

   // Nodes on symmetry planes
   Index_t symmX(Index_t idx) { return m_symmX[idx] ; }
   Index_t symmY(Index_t idx) { return m_symmY[idx] ; }
   Index_t symmZ(Index_t idx) { return m_symmZ[idx] ; }
   bool symmXempty()          { return m_symmX.empty(); }
   bool symmYempty()          { return m_symmY.empty(); }
   bool symmZempty()          { return m_symmZ.empty(); }
#endif
   //
   // Element-centered
   //
   Index_t  regElemSize(Index_t idx) { return m_regionElementsRel.size(idx) ; }
   Index_t&  regNumList(Index_t idx) { return m_elemRegNum[idx] ; }
   Index_t*  regNumList()            { return &m_elemRegNum[0] ; }
   const Index_t*  regElemlist(Int_t r)    { return &(*m_regionElementsRel.begin(r)); }
   const Index_t&  regElemlist(Int_t r, Index_t idx) { return m_regionElementsRel[r][idx] ; }

   const Index_t*  nodelist(Index_t idx)    { return &(*m_nodelist.begin(idx)) ; }

   // elem adjacencies through face
   const Index_t&  lxim(Index_t idx) { return *m_lxim.begin(idx); }
   const Index_t&  lxip(Index_t idx) { return *m_lxip.begin(idx); }
   const Index_t&  letam(Index_t idx) { return *m_letam.begin(idx); }
   const Index_t&  letap(Index_t idx) { return *m_letap.begin(idx); }
   const Index_t&  lzetam(Index_t idx) { return *m_lzetam.begin(idx); }
   const Index_t&  lzetap(Index_t idx) { return *m_lzetap.begin(idx); }

#if USE_SIDRE==1
   // elem face symm/free-surface flag
   Int_t&  elemBC(Index_t idx) { return static_cast<Int_t*>(m_elemBC->getData())[idx] ; }

   // Principal strains - temporary
   Real_t& dxx(Index_t idx)  { return static_cast<Real_t*>(m_dxx->getData())[idx] ; }
   Real_t& dyy(Index_t idx)  { return static_cast<Real_t*>(m_dyy->getData())[idx] ; }
   Real_t& dzz(Index_t idx)  { return static_cast<Real_t*>(m_dzz->getData())[idx] ; }

   // Velocity gradient - temporary
   Real_t& delv_xi(Index_t idx)    { return static_cast<Real_t*>(m_delv_xi->getData())[idx] ; }
   Real_t& delv_eta(Index_t idx)   { return static_cast<Real_t*>(m_delv_eta->getData())[idx] ; }
   Real_t& delv_zeta(Index_t idx)  { return static_cast<Real_t*>(m_delv_zeta->getData())[idx] ; }

   // Position gradient - temporary
   Real_t& delx_xi(Index_t idx)    { return static_cast<Real_t*>(m_delx_xi->getData())[idx] ; }
   Real_t& delx_eta(Index_t idx)   { return static_cast<Real_t*>(m_delx_eta->getData())[idx] ; }
   Real_t& delx_zeta(Index_t idx)  { return static_cast<Real_t*>(m_delx_zeta->getData())[idx] ; }

   // Energy
   Real_t& e(Index_t idx)          { return static_cast<Real_t*>(m_e->getData())[idx] ; }

   // Pressure
   Real_t& p(Index_t idx)          { return static_cast<Real_t*>(m_p->getData())[idx] ; }

   // Artificial viscosity
   Real_t& q(Index_t idx)          { return static_cast<Real_t*>(m_q->getData())[idx] ; }

   // Linear term for q
   Real_t& ql(Index_t idx)         { return static_cast<Real_t*>(m_ql->getData())[idx] ; }
   // Quadratic term for q
   Real_t& qq(Index_t idx)         { return static_cast<Real_t*>(m_qq->getData())[idx] ; }

   // Relative volume
   Real_t& v(Index_t idx)          { return static_cast<Real_t*>(m_v->getData())[idx] ; }
   Real_t& delv(Index_t idx)       { return static_cast<Real_t*>(m_delv->getData())[idx] ; }

   // Reference volume
   Real_t& volo(Index_t idx)       { return static_cast<Real_t*>(m_volo->getData())[idx] ; }

   // volume derivative over volume
   Real_t& vdov(Index_t idx)       { return static_cast<Real_t*>(m_vdov->getData())[idx] ; }

   // Element characteristic length
   Real_t& arealg(Index_t idx)     { return static_cast<Real_t*>(m_arealg->getData())[idx] ; }

   // Sound speed
   Real_t& ss(Index_t idx)         { return static_cast<Real_t*>(m_ss->getData())[idx] ; }

   // Element mass
   Real_t& elemMass(Index_t idx)  { return static_cast<Real_t*>(m_elemMass->getData())[idx] ; }
#else
   // elem face symm/free-surface flag
   Int_t&  elemBC(Index_t idx) { return m_elemBC[idx] ; }

   // Principal strains - temporary
   Real_t& dxx(Index_t idx)  { return m_dxx[idx] ; }
   Real_t& dyy(Index_t idx)  { return m_dyy[idx] ; }
   Real_t& dzz(Index_t idx)  { return m_dzz[idx] ; }

   // Velocity gradient - temporary
   Real_t& delv_xi(Index_t idx)    { return m_delv_xi[idx] ; }
   Real_t& delv_eta(Index_t idx)   { return m_delv_eta[idx] ; }
   Real_t& delv_zeta(Index_t idx)  { return m_delv_zeta[idx] ; }

   // Position gradient - temporary
   Real_t& delx_xi(Index_t idx)    { return m_delx_xi[idx] ; }
   Real_t& delx_eta(Index_t idx)   { return m_delx_eta[idx] ; }
   Real_t& delx_zeta(Index_t idx)  { return m_delx_zeta[idx] ; }

   // Energy
   Real_t& e(Index_t idx)          { return m_e[idx] ; }

   // Pressure
   Real_t& p(Index_t idx)          { return m_p[idx] ; }

   // Artificial viscosity
   Real_t& q(Index_t idx)          { return m_q[idx] ; }

   // Linear term for q
   Real_t& ql(Index_t idx)         { return m_ql[idx] ; }
   // Quadratic term for q
   Real_t& qq(Index_t idx)         { return m_qq[idx] ; }

   // Relative volume
   Real_t& v(Index_t idx)          { return m_v[idx] ; }
   Real_t& delv(Index_t idx)       { return m_delv[idx] ; }

   // Reference volume
   Real_t& volo(Index_t idx)       { return m_volo[idx] ; }

   // volume derivative over volume
   Real_t& vdov(Index_t idx)       { return m_vdov[idx] ; }

   // Element characteristic length
   Real_t& arealg(Index_t idx)     { return m_arealg[idx] ; }

   // Sound speed
   Real_t& ss(Index_t idx)         { return m_ss[idx] ; }

   // Element mass
   Real_t& elemMass(Index_t idx)  { return m_elemMass[idx] ; }
#endif
   Index_t nodeElemCount(Index_t idx)
   { return m_nodeElemStart[idx+1] - m_nodeElemStart[idx] ; }

   Index_t *nodeElemCornerList(Index_t idx)
   { return &m_nodeElemCornerList[m_nodeElemStart[idx]] ; }

   // Parameters 

   // Cutoffs
   Real_t u_cut() const               { return m_u_cut ; }
   Real_t e_cut() const               { return m_e_cut ; }
   Real_t p_cut() const               { return m_p_cut ; }
   Real_t q_cut() const               { return m_q_cut ; }
   Real_t v_cut() const               { return m_v_cut ; }

   // Other constants (usually are settable via input file in real codes)
   Real_t hgcoef() const              { return m_hgcoef ; }
   Real_t qstop() const               { return m_qstop ; }
   Real_t monoq_max_slope() const     { return m_monoq_max_slope ; }
   Real_t monoq_limiter_mult() const  { return m_monoq_limiter_mult ; }
   Real_t ss4o3() const               { return m_ss4o3 ; }
   Real_t qlc_monoq() const           { return m_qlc_monoq ; }
   Real_t qqc_monoq() const           { return m_qqc_monoq ; }
   Real_t qqc() const                 { return m_qqc ; }

   Real_t eosvmax() const             { return m_eosvmax ; }
   Real_t eosvmin() const             { return m_eosvmin ; }
   Real_t pmin() const                { return m_pmin ; }
   Real_t emin() const                { return m_emin ; }
   Real_t dvovmax() const             { return m_dvovmax ; }
   Real_t refdens() const             { return m_refdens ; }

   // Timestep controls, etc...
   Real_t& time()                 { return m_time ; }
   Real_t& deltatime()            { return m_deltatime ; }
   Real_t& deltatimemultlb()      { return m_deltatimemultlb ; }
   Real_t& deltatimemultub()      { return m_deltatimemultub ; }
   Real_t& stoptime()             { return m_stoptime ; }
   Real_t& dtcourant()            { return m_dtcourant ; }
   Real_t& dthydro()              { return m_dthydro ; }
   Real_t& dtmax()                { return m_dtmax ; }
   Real_t& dtfixed()              { return m_dtfixed ; }

   Int_t&  cycle()                { return m_cycle ; }
   Index_t&  numRanks()           { return m_numRanks ; }

   Index_t&  colLoc()             { return m_colLoc ; }
   Index_t&  rowLoc()             { return m_rowLoc ; }
   Index_t&  planeLoc()           { return m_planeLoc ; }
   Index_t&  tp()                 { return m_tp ; }

   Index_t&  sizeX()              { return m_sizeX ; }
   Index_t&  sizeY()              { return m_sizeY ; }
   Index_t&  sizeZ()              { return m_sizeZ ; }

   Index_t  numReg()              { return m_regionSet.size() ; }
   Int_t&  cost()                 { return m_cost ; }

   Index_t  numElem()  const      { return m_elemSet.size() ; }
   Index_t  numNode()  const      { return m_nodeSet.size() ; }
   
   Index_t&  maxPlaneSize()       { return m_maxPlaneSize ; }
   Index_t&  maxEdgeSize()        { return m_maxEdgeSize ; }
   
   //
   // MPI-Related additional data
   //

#ifdef USE_MPI
   // Communication Work space 
   Real_t *commDataSend ;
   Real_t *commDataRecv ;
   
   // Maximum number of block neighbors 
   MPI_Request recvRequest[26] ; // 6 faces + 12 edges + 8 corners 
   MPI_Request sendRequest[26] ; // 6 faces + 12 edges + 8 corners 
#endif

  private:

   void BuildMesh(Int_t nx, Int_t edgeNodes, Int_t edgeElems);
   void SetupThreadSupportStructures();
   void CreateRegionIndexSets(Int_t nreg, Int_t balance);
   void SetupCommBuffers(Int_t edgeNodes);
   void SetupSymmetryPlanes(Int_t edgeNodes);
   void SetupElementConnectivities(Int_t edgeElems);
   void SetupBoundaryConditions(Int_t edgeElems);

   //
   // IMPLEMENTATION
   //

#if USE_SIDRE==1
   asctoolkit::sidre::DataStore m_DataStore;
   asctoolkit::sidre::DataGroup * m_DataGroup;

#endif
   /* Node-centered */
   luleshRealData m_x ;  /* coordinates */
   luleshRealData m_y ;
   luleshRealData m_z ;

   luleshRealData m_xd ; /* velocities */
   luleshRealData m_yd ;
   luleshRealData m_zd ;

   luleshRealData m_xdd ; /* accelerations */
   luleshRealData m_ydd ;
   luleshRealData m_zdd ;

   luleshRealData m_fx ;  /* forces */
   luleshRealData m_fy ;
   luleshRealData m_fz ;

   luleshRealData m_nodalMass ;  /* mass */

   luleshIndexData m_symmX ;  /* symmetry plane nodesets */
   luleshIndexData m_symmY ;
   luleshIndexData m_symmZ ;

   // Element-centered

   // Region information
   //Int_t    m_numReg ;
   RegionSet m_regionSet;

   Int_t    m_cost; //imbalance cost

   ElemIntMap           m_elemRegNum;
   RegionToElemRelation m_regionElementsRel;

//   Index_t *m_regElemSize ;   // Size of region sets
//   Index_t *m_regNumList ;    // Region number per domain element
//   Index_t **m_regElemlist ;  // region indexset

   ElemToNodeRelation m_nodelist;
   //std::vector<Index_t>  m_nodelist ;     /* elemToNode connectivity */


   ElemFaceAdjacencyRelation m_lxim ;  /* element connectivity across each face */
   ElemFaceAdjacencyRelation m_lxip ;
   ElemFaceAdjacencyRelation m_letam ;
   ElemFaceAdjacencyRelation m_letap ;
   ElemFaceAdjacencyRelation m_lzetam ;
   ElemFaceAdjacencyRelation m_lzetap ;

   luleshIntData    m_elemBC ;  /* symmetry/free-surface flags for each elem face */

   luleshRealData m_dxx ;  /* principal strains -- temporary */
   luleshRealData m_dyy ;
   luleshRealData m_dzz ;

   luleshRealData m_delv_xi ;    /* velocity gradient -- temporary */
   luleshRealData m_delv_eta ;
   luleshRealData m_delv_zeta ;

   luleshRealData m_delx_xi ;    /* coordinate gradient -- temporary */
   luleshRealData m_delx_eta ;
   luleshRealData m_delx_zeta ;
   
   luleshRealData m_e ;   /* energy */

   luleshRealData m_p ;   /* pressure */
   luleshRealData m_q ;   /* q */
   luleshRealData m_ql ;  /* linear term for q */
   luleshRealData m_qq ;  /* quadratic term for q */

   luleshRealData m_v ;     /* relative volume */
   luleshRealData m_volo ;  /* reference volume */
   luleshRealData m_vnew ;  /* new relative volume -- temporary */
   luleshRealData m_delv ;  /* m_vnew - m_v */
   luleshRealData m_vdov ;  /* volume derivative over volume */

   luleshRealData m_arealg ;  /* characteristic length of an element */
   
   luleshRealData m_ss ;      /* "sound speed" */

   luleshRealData m_elemMass ;  /* mass */

   // Cutoffs (treat as constants)
   const Real_t  m_e_cut ;             // energy tolerance 
   const Real_t  m_p_cut ;             // pressure tolerance 
   const Real_t  m_q_cut ;             // q tolerance 
   const Real_t  m_v_cut ;             // relative volume tolerance 
   const Real_t  m_u_cut ;             // velocity tolerance 

   // Other constants (usually setable, but hardcoded in this proxy app)

   const Real_t  m_hgcoef ;            // hourglass control 
   const Real_t  m_ss4o3 ;
   const Real_t  m_qstop ;             // excessive q indicator 
   const Real_t  m_monoq_max_slope ;
   const Real_t  m_monoq_limiter_mult ;
   const Real_t  m_qlc_monoq ;         // linear term coef for q 
   const Real_t  m_qqc_monoq ;         // quadratic term coef for q 
   const Real_t  m_qqc ;
   const Real_t  m_eosvmax ;
   const Real_t  m_eosvmin ;
   const Real_t  m_pmin ;              // pressure floor 
   const Real_t  m_emin ;              // energy floor 
   const Real_t  m_dvovmax ;           // maximum allowable volume change 
   const Real_t  m_refdens ;           // reference density 

   // Variables to keep track of timestep, simulation time, and cycle
   Real_t  m_dtcourant ;         // courant constraint 
   Real_t  m_dthydro ;           // volume change constraint 
   Int_t   m_cycle ;             // iteration count for simulation 
   Real_t  m_dtfixed ;           // fixed time increment 
   Real_t  m_time ;              // current time 
   Real_t  m_deltatime ;         // variable time increment 
   Real_t  m_deltatimemultlb ;
   Real_t  m_deltatimemultub ;
   Real_t  m_dtmax ;             // maximum allowable time increment 
   Real_t  m_stoptime ;          // end time for simulation 


   Int_t   m_numRanks ;

   Index_t m_colLoc ;
   Index_t m_rowLoc ;
   Index_t m_planeLoc ;
   Index_t m_tp ;

   Index_t m_sizeX ;
   Index_t m_sizeY ;
   Index_t m_sizeZ ;

   NodeSet m_nodeSet;
   ElemSet m_elemSet;

   Index_t m_maxPlaneSize ;
   Index_t m_maxEdgeSize ;

   // OMP hack 
   Index_t *m_nodeElemStart ;
   Index_t *m_nodeElemCornerList ;

   // Used in setup
   Index_t m_rowMin, m_rowMax;
   Index_t m_colMin, m_colMax;
   Index_t m_planeMin, m_planeMax ;

} ;

typedef Real_t &(Domain::* Domain_member )(Index_t) ;

struct cmdLineOpts {
   Int_t its; // -i 
   Int_t nx;  // -s 
   Int_t numReg; // -r 
   Int_t numFiles; // -f
   Int_t showProg; // -p
   Int_t quiet; // -q
   Int_t viz; // -v 
   Int_t cost; // -c
   Int_t balance; // -b
};



// Function Prototypes

// lulesh-par
Real_t CalcElemVolume( const Real_t x[8],
                       const Real_t y[8],
                       const Real_t z[8]);

// lulesh-util
void ParseCommandLineOptions(int argc, char *argv[],
                             Int_t myRank, struct cmdLineOpts *opts);
void VerifyAndWriteFinalOutput(Real_t elapsed_time,
                               Domain& locDom,
                               Int_t nx,
                               Int_t numRanks);

// lulesh-viz
void DumpToVisit(Domain& domain, int numFiles, int myRank, int numRanks);

// lulesh-comm
void CommRecv(Domain& domain, Int_t msgType, Index_t xferFields,
              Index_t dx, Index_t dy, Index_t dz,
              bool doRecv, bool planeOnly);
void CommSend(Domain& domain, Int_t msgType,
              Index_t xferFields, Domain_member *fieldData,
              Index_t dx, Index_t dy, Index_t dz,
              bool doSend, bool planeOnly);
void CommSBN(Domain& domain, Int_t xferFields, Domain_member *fieldData);
void CommSyncPosVel(Domain& domain);
void CommMonoQ(Domain& domain);

// lulesh-init
void InitMeshDecomp(Int_t numRanks, Int_t myRank,
                    Int_t *col, Int_t *row, Int_t *plane, Int_t *side);
