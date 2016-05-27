// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_BlockedEpetraLinearObjFactory_hpp__
#define __Panzer_BlockedEpetraLinearObjFactory_hpp__

#include <map>

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_CloneableEvaluator.hpp"

#include "Panzer_GatherOrientation.hpp"
#include "Panzer_GatherSolution_BlockedEpetra.hpp"
#include "Panzer_GatherTangent_BlockedEpetra.hpp"
#include "Panzer_ScatterResidual_BlockedEpetra.hpp"
#include "Panzer_ScatterDirichletResidual_BlockedEpetra.hpp"
#include "Panzer_ThyraObjFactory.hpp"
#include "Panzer_HashUtils.hpp"

#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

namespace panzer {

template <typename Traits,typename LocalOrdinalT>
class BlockedEpetraLinearObjFactory : public LinearObjFactory<Traits>
                                    , public ThyraObjFactory<double> {
public:

   BlockedEpetraLinearObjFactory(const Teuchos::RCP<const Epetra_Comm> & comm,
                                 const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,std::pair<int,int> > > & blkProvider,
                                 const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > > & gidProviders);
   BlockedEpetraLinearObjFactory(const Teuchos::RCP<const Epetra_Comm> & comm,
                                 const Teuchos::RCP<const BlockedDOFManager<LocalOrdinalT,int> > & gidProvider);
   BlockedEpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                                 const Teuchos::RCP<const BlockedDOFManager<LocalOrdinalT,int> > & gidProvider);

   BlockedEpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                                 const Teuchos::RCP<const BlockedDOFManager<LocalOrdinalT,int> > & gidProvider,
                                 const Teuchos::RCP<const BlockedDOFManager<LocalOrdinalT,int> > & colGidProvider);

   virtual ~BlockedEpetraLinearObjFactory();

/*************** Linear object factory methods *******************/

   virtual void readVector(const std::string & identifier,LinearObjContainer & loc,int id) const;

   virtual void writeVector(const std::string & identifier,const LinearObjContainer & loc,int id) const;

   virtual Teuchos::RCP<LinearObjContainer> buildLinearObjContainer() const;

   virtual Teuchos::RCP<LinearObjContainer> buildPrimitiveLinearObjContainer() const 
   { return buildLinearObjContainer(); }

   virtual Teuchos::RCP<LinearObjContainer> buildGhostedLinearObjContainer() const;

   virtual Teuchos::RCP<LinearObjContainer> buildPrimitiveGhostedLinearObjContainer() const 
   { return buildGhostedLinearObjContainer(); }

   virtual void globalToGhostContainer(const LinearObjContainer & container,
                                       LinearObjContainer & ghostContainer,int) const;
   virtual void ghostToGlobalContainer(const LinearObjContainer & ghostContainer,
                                       LinearObjContainer & container,int) const;

   /** Adjust the residual vector and Jacobian matrix (if they exist) for applied
     * dirichlet conditions. The adjustment considers if a boundary condition was
     * set globally and locally and based on that result adjust the ghosted matrix
     * and residual vector so that when they are summed across processors they resulting
     * Dirichlet condition is correct.
     */
   virtual void adjustForDirichletConditions(const LinearObjContainer & localBCRows,
                                             const LinearObjContainer & globalBCRows,
                                             LinearObjContainer & ghostedObjs,
                                             bool zeroVectorRows=false, bool adjustX=false) const;

   /** Adjust a vector by replacing selected rows with the value of the evaluated
     * dirichlet conditions. This is handled through the standard container mechanism.
     */
   virtual void applyDirichletBCs(const LinearObjContainer & counter,
                                  LinearObjContainer & result) const;

   /** Build a GlobalEvaluationDataContainer that handles all domain communication.
     * This is used primarily for gather operations and hides the allocation and usage
     * of the ghosted vector from the user.
     */
   virtual Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData> buildDomainContainer() const;

   virtual Teuchos::MpiComm<int> getComm() const;

   //! Use preconstructed scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatter() const
   { return Teuchos::rcp(new ScatterResidual_BlockedEpetra<EvalT,Traits,LocalOrdinalT,int>(blockedDOFManager_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGather() const
   { return Teuchos::rcp(new GatherSolution_BlockedEpetra<EvalT,Traits,LocalOrdinalT,int>(blockedDOFManager_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherTangent() const
   { return Teuchos::rcp(new GatherTangent_BlockedEpetra<EvalT,Traits,LocalOrdinalT,int>(blockedDOFManager_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherDomain() const
   { return Teuchos::rcp(new GatherSolution_BlockedEpetra<EvalT,Traits,LocalOrdinalT,int>(blockedDOFManager_)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator > buildGatherOrientation() const
   { return Teuchos::rcp(new GatherOrientation<EvalT,Traits,LocalOrdinalT,std::pair<int,int> >(blockProvider_)); }

   //! Use preconstructed dirichlet scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const
   { return Teuchos::rcp(new ScatterDirichletResidual_BlockedEpetra<EvalT,Traits,LocalOrdinalT,int>(blockedDOFManager_)); }

/*************** Generic helper functions for container setup *******************/
   
   /** Initialize container with a specific set of member values.
     *
     * \note This will overwrite everything in the container and zero out values
     *       not requested.
     */
   void initializeContainer(int,LinearObjContainer & loc) const;

   /** Initialize container with a specific set of member values.
     *
     * \note This will overwrite everything in the container and zero out values
     *       not requested.
     */
   void initializeContainer(int mem,BlockedEpetraLinearObjContainer & loc) const;

   /** Initialize container with a specific set of member values.
     *
     * \note This will overwrite everything in the container and zero out values
     *       not requested.
     */
   void initializeGhostedContainer(int,LinearObjContainer & loc) const;

   /** Initialize container with a specific set of member values.
     *
     * \note This will overwrite everything in the container and zero out values
     *       not requested.
     */
   void initializeGhostedContainer(int mem,BlockedEpetraLinearObjContainer & loc) const;

/*************** Thyra based methods *******************/

   //! Get the domain vector space (x and dxdt)
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getThyraDomainSpace() const;

   //! Get the range vector space (f)
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getThyraRangeSpace() const;

   //! Get a domain vector
   Teuchos::RCP<Thyra::VectorBase<double> > getThyraDomainVector() const;

   //! Get a range vector
   Teuchos::RCP<Thyra::VectorBase<double> > getThyraRangeVector() const;

   //! Get a Thyra operator
   Teuchos::RCP<Thyra::LinearOpBase<double> > getThyraMatrix() const;

   // and now the ghosted versions

   //! Get the domain vector space (x and dxdt)
   Teuchos::RCP<Thyra::VectorSpaceBase<double> > getGhostedThyraDomainSpace() const;

   //! Get the range vector space (f)
   Teuchos::RCP<Thyra::VectorSpaceBase<double> > getGhostedThyraRangeSpace() const;

   //! Get a domain vector
   Teuchos::RCP<Thyra::VectorBase<double> > getGhostedThyraDomainVector() const;

   //! Get a range vector
   Teuchos::RCP<Thyra::VectorBase<double> > getGhostedThyraRangeVector() const;

   //! Get a Thyra operator
   Teuchos::RCP<Thyra::BlockedLinearOpBase<double> > getGhostedThyraMatrix() const;

/*************** Epetra based methods *******************/

   //! get the map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getMap(int i) const;

   //! get the map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getColMap(int i) const;

   //! get the ghosted map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getGhostedMap(int i) const;

   //! get the ghosted map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getColGhostedMap(int i) const;

   //! get the graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> getGraph(int i,int j) const;

   //! get the ghosted graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> getGhostedGraph(int i,int j) const;

   //! get importer for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Import> getGhostedImport(int i) const;

   //! get importer for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Import> getColGhostedImport(int i) const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Export> getGhostedExport(int j) const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Export> getColGhostedExport(int j) const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<const Epetra_Comm> getEpetraComm() const;

   Teuchos::RCP<Epetra_CrsMatrix> getEpetraMatrix(int i,int j) const;
   Teuchos::RCP<Epetra_CrsMatrix> getGhostedEpetraMatrix(int i,int j) const;

   //! how many block rows
   int getBlockRowCount() const;

   //! how many block columns
   int getBlockColCount() const;

   Teuchos::RCP<const BlockedDOFManager<LocalOrdinalT,int> > getGlobalIndexer() const
   { return blockedDOFManager_; }

   Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> getRangeGlobalIndexer() const
   { return blockedDOFManager_; }

   Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> getDomainGlobalIndexer() const
   { return blockedDOFManager_; }

   //! exclude a block pair from the matrix
   void addExcludedPair(int rowBlock,int colBlock);

   //! exclude a vector of pairs from the matrix
   void addExcludedPairs(const std::vector<std::pair<int,int> > & exPairs);

   //! Get the unique global indexer this factory was created with.
   Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> getUniqueGlobalIndexerBase() const
   { return blockProvider_; }

protected:
/*************** Generic methods/members *******************/

   // Get the global indexer associated with a particular block
   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > getGlobalIndexer(int i) const;

   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > getColGlobalIndexer(int i) const;

   //! How many gid providers are there?
   std::size_t getGidProviderCount() const;

   //! How many column gid providers are there?
   std::size_t getColGidProviderCount() const;

   //! Allocate the space in the std::vector objects so we can fill with appropriate Epetra data
   void makeRoomForBlocks(std::size_t blockCnt,std::size_t colBlockCnt=0);

   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,std::pair<int,int> > > blockProvider_;
   Teuchos::RCP<const BlockedDOFManager<LocalOrdinalT,int> > blockedDOFManager_;
   std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > > gidProviders_;
   std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > > colGidProviders_;

   bool useColGidProviders_;

   // which block entries are ignored
   std::unordered_set<std::pair<int,int>,panzer::pair_hash> excludedPairs_;
  
/*************** Thyra based methods/members *******************/

   void ghostToGlobalThyraVector(const Teuchos::RCP<const Thyra::VectorBase<double> > & in,
                                 const Teuchos::RCP<Thyra::VectorBase<double> > & out,bool col) const;
   void ghostToGlobalThyraMatrix(const Thyra::LinearOpBase<double> & in,Thyra::LinearOpBase<double> & out) const;
   void globalToGhostThyraVector(const Teuchos::RCP<const Thyra::VectorBase<double> > & in,
                                 const Teuchos::RCP<Thyra::VectorBase<double> > & out,bool col) const;

   mutable Teuchos::RCP<Thyra::ProductVectorSpaceBase<double> > rangeSpace_;
   mutable Teuchos::RCP<Thyra::ProductVectorSpaceBase<double> > domainSpace_;

   mutable Teuchos::RCP<Thyra::ProductVectorSpaceBase<double> > ghostedRangeSpace_;
   mutable Teuchos::RCP<Thyra::ProductVectorSpaceBase<double> > ghostedDomainSpace_;

/*************** Epetra based methods/members *******************/

   void adjustForDirichletConditions(const Epetra_Vector & local_bcs,
                                     const Epetra_Vector & global_bcs,
                                     const Teuchos::Ptr<Epetra_Vector> & f,
                                     const Teuchos::Ptr<Epetra_CrsMatrix> & A,
                                     bool zeroVectorRows) const;

   void ghostToGlobalEpetraVector(int i,const Epetra_Vector & in,Epetra_Vector & out,bool col) const;
   void ghostToGlobalEpetraMatrix(int blockRow,const Epetra_CrsMatrix & in,Epetra_CrsMatrix & out) const;
   void globalToGhostEpetraVector(int i,const Epetra_Vector & in,Epetra_Vector & out,bool col) const;

   // get the map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> buildEpetraMap(int i) const;
   virtual const Teuchos::RCP<Epetra_Map> buildEpetraGhostedMap(int i) const;

   // get the map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> buildEpetraColMap(int i) const;
   virtual const Teuchos::RCP<Epetra_Map> buildEpetraColGhostedMap(int i) const;

   // get the graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> buildEpetraGraph(int i,int j) const;
   virtual const Teuchos::RCP<Epetra_CrsGraph> buildEpetraGhostedGraph(int i,int j) const;

   // storage for Epetra graphs and maps
   Teuchos::RCP<const Epetra_Comm> comm_;
   Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm_;
   Teuchos::RCP<Teuchos::MpiComm<int> > tComm_;

   mutable std::vector<Teuchos::RCP<Epetra_Map> > maps_;
   mutable std::vector<Teuchos::RCP<Epetra_Map> > ghostedMaps_;
   mutable std::vector<Teuchos::RCP<Epetra_Import> > importers_;
   mutable std::vector<Teuchos::RCP<Epetra_Export> > exporters_;

   mutable std::vector<Teuchos::RCP<Epetra_Map> > colMaps_;
   mutable std::vector<Teuchos::RCP<Epetra_Map> > colGhostedMaps_;
   mutable std::vector<Teuchos::RCP<Epetra_Import> > colImporters_;
   mutable std::vector<Teuchos::RCP<Epetra_Export> > colExporters_;

   mutable std::unordered_map<std::pair<int,int>,Teuchos::RCP<Epetra_CrsGraph>,panzer::pair_hash> graphs_ ;
   mutable std::unordered_map<std::pair<int,int>,Teuchos::RCP<Epetra_CrsGraph>,panzer::pair_hash> ghostedGraphs_;
};

}

#endif
