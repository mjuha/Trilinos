/*
// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef SUN_CXX

#include "Thyra_Amesos2LinearOpWithSolveFactory_decl.hpp"

#include "Thyra_Amesos2LinearOpWithSolve.hpp"
//#include "Thyra_EpetraOperatorViewExtractorStd.hpp"
#include "Amesos2.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_TypeTraits.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "Amesos2_KLU2.hpp"

// #ifdef HAVE_AMESOS_PASTIX
// #include "Amesos_Pastix.h"
// #endif
// #ifdef HAVE_AMESOS_LAPACK
// #include "Amesos_Lapack.h"
// #endif
// #ifdef HAVE_AMESOS_MUMPS
// #include "Amesos_Mumps.h"
// #endif
// #ifdef HAVE_AMESOS_SCALAPACK
// #include "Amesos_Scalapack.h"
// #endif
// #ifdef HAVE_AMESOS_UMFPACK
// #include "Amesos_Umfpack.h"
// #endif
// #ifdef HAVE_AMESOS_SUPERLUDIST
// #include "Amesos_Superludist.h"
// #endif
// #ifdef HAVE_AMESOS_SUPERLU
// #include "Amesos_Superlu.h"
// #endif
// #ifdef HAVE_AMESOS_DSCPACK
// #include "Amesos_Dscpack.h"
// #endif
// #ifdef HAVE_AMESOS_PARDISO
// #include "Amesos_Pardiso.h"
// #endif
// #ifdef HAVE_AMESOS_TAUCS
// #include "Amesos_Taucs.h"
// #endif
// #ifdef HAVE_AMESOS_PARAKLETE
// #include "Amesos_Paraklete.h"
// #endif

namespace {

  const std::string epetraFwdOp_str = "epetraFwdOp";

} // namespace

namespace Thyra {


// Parameter names for Paramter List

template<typename Scalar>
const std::string Amesos2LinearOpWithSolveFactory<Scalar>::SolverType_name = "Solver Type";
  
template<typename Scalar>
const std::string Amesos2LinearOpWithSolveFactory<Scalar>::RefactorizationPolicy_name = "Refactorization Policy";

template<typename Scalar>
const std::string Amesos2LinearOpWithSolveFactory<Scalar>::ThrowOnPreconditionerInput_name = "Throw on Preconditioner Input";

template<typename Scalar>
const std::string Amesos2LinearOpWithSolveFactory<Scalar>::Amesos_Settings_name = "Amesos Settings";

// Constructors/initializers/accessors

template<typename Scalar>
Amesos2LinearOpWithSolveFactory<Scalar>::~Amesos2LinearOpWithSolveFactory()
{
#ifdef TEUCHOS_DEBUG
  if(paramList_.get())
    paramList_->validateParameters(
      *this->getValidParameters(),0  // Only validate this level for now!
      );
#endif
}

template<typename Scalar>
Amesos2LinearOpWithSolveFactory<Scalar>::Amesos2LinearOpWithSolveFactory(
  const Amesos2::ESolverType                            solverType
  ,const Amesos2::ERefactorizationPolicy                refactorizationPolicy
  ,const bool                                          throwOnPrecInput
    )
  ://epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd()))
  //,solverType_(solverType)
  solverType_(solverType)
  ,refactorizationPolicy_(refactorizationPolicy)
  ,throwOnPrecInput_(throwOnPrecInput)
{}

// Overridden from LinearOpWithSolveFactoryBase

template<typename Scalar>
bool Amesos2LinearOpWithSolveFactory<Scalar>::isCompatible(
  const LinearOpSourceBase<Scalar> &fwdOpSrc
  ) const
{
  // using Teuchos::outArg;
  // RCP<const LinearOpBase<double> >
  //   fwdOp = fwdOpSrc.getOp();
  // RCP<const Epetra_Operator> epetraFwdOp;
  // EOpTransp epetraFwdOpTransp;
  // EApplyEpetraOpAs epetraFwdOpApplyAs;
  // EAdjointEpetraOp epetraFwdOpAdjointSupport;
  // double epetraFwdOpScalar;
  // epetraFwdOpViewExtractor_->getEpetraOpView(
  //   fwdOp,
  //   outArg(epetraFwdOp), outArg(epetraFwdOpTransp),
  //   outArg(epetraFwdOpApplyAs), outArg(epetraFwdOpAdjointSupport),
  //   outArg(epetraFwdOpScalar)
  //   );
  // if( !dynamic_cast<const Epetra_RowMatrix*>(&*epetraFwdOp) )
  //   return false;
  // return true;
  return true;
}

template<typename Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
Amesos2LinearOpWithSolveFactory<Scalar>::createOp() const
{
  return Teuchos::rcp(new Amesos2LinearOpWithSolve<Scalar>());
}

template<typename Scalar>
void Amesos2LinearOpWithSolveFactory<Scalar>::initializeOp(
  const RCP<const LinearOpSourceBase<Scalar> >    &fwdOpSrc
  ,LinearOpWithSolveBase<Scalar>                                   *Op
  ,const ESupportSolveUse                                          supportSolveUse
  ) const
{
  // using Teuchos::outArg;
//   THYRA_FUNC_TIME_MONITOR("Stratimikos: Amesos2LOWSF");
// #ifdef TEUCHOS_DEBUG
//   TEUCHOS_TEST_FOR_EXCEPT(Op==NULL);
// #endif
  const RCP<const LinearOpBase<Scalar> > 
    fwdOp = fwdOpSrc->getOp();
//   //
//   // Unwrap and get the forward Epetra_Operator object
//   //
//   RCP<const Epetra_Operator> epetraFwdOp;
//   EOpTransp epetraFwdOpTransp;
//   EApplyEpetraOpAs epetraFwdOpApplyAs;
//   EAdjointEpetraOp epetraFwdOpAdjointSupport;
//   double epetraFwdOpScalar;
//   epetraFwdOpViewExtractor_->getEpetraOpView(
//     fwdOp,
//     outArg(epetraFwdOp), outArg(epetraFwdOpTransp),
//     outArg(epetraFwdOpApplyAs), outArg(epetraFwdOpAdjointSupport),
//     outArg(epetraFwdOpScalar)
//     );
//   // Get the Amesos2LinearOpWithSolve object
//   Amesos2LinearOpWithSolve
//     *amesosOp = &Teuchos::dyn_cast<Amesos2LinearOpWithSolve>(*Op);
//   //
//   // Determine if we must start over or not
//   //
//   bool startOver = ( amesosOp->get_amesosSolver()==Teuchos::null );
//   if(!startOver) {
//     startOver =
//       (
//         epetraFwdOpTransp != amesosOp->get_amesosSolverTransp() ||
//         epetraFwdOp.get() != amesosOp->get_epetraLP()->GetOperator()
//         // We must start over if the matrix object changes.  This is a
//         // weakness of Amesos but there is nothing I can do about this right
//         // now!
//         );
//   }
//   //
//   // Update the amesos solver
//   //
//   if(startOver) {
//     //
//     // This LOWS object has not be initialized yet or is not compatible with the existing
//     // 
//     // so this is where we setup everything from the ground up.
//     //
//     // Create the linear problem and set the operator with memory of RCP to Epetra_Operator view!
//     RCP<Epetra_LinearProblem>
//       epetraLP = Teuchos::rcp(new Epetra_LinearProblem());
//     epetraLP->SetOperator(const_cast<Epetra_Operator*>(&*epetraFwdOp));
//     Teuchos::set_extra_data< RCP<const Epetra_Operator> >( epetraFwdOp, epetraFwdOp_str,
//      Teuchos::inOutArg(epetraLP) );
//     // Create the concrete solver
//     RCP<Amesos_BaseSolver>
//       amesosSolver;
//     {
//       THYRA_FUNC_TIME_MONITOR_DIFF("Stratimikos: AmesosLOWSF:InitConstruct",
//         InitConstruct);
//       switch(solverType_) {
//         case Thyra::Amesos2::LAPACK :
//           amesosSolver = Teuchos::rcp(new Amesos_Lapack(*epetraLP));
//           break;
// #ifdef HAVE_AMESOS_KLU
//         case Thyra::Amesos2::KLU :
//           amesosSolver = Teuchos::rcp(new Amesos_Klu(*epetraLP));
//           break;
// #endif
// #ifdef HAVE_AMESOS_PASTIX
//         case Thyra::Amesos2::PASTIX :
//           amesosSolver = Teuchos::rcp(new Amesos_Pastix(*epetraLP));
//           break;
// #endif
// #ifdef HAVE_AMESOS_MUMPS
//         case Thyra::Amesos2::MUMPS :
//           amesosSolver = Teuchos::rcp(new Amesos_Mumps(*epetraLP));
//           break;
// #endif
// #ifdef HAVE_AMESOS_SCALAPACK
//         case Thyra::Amesos2::SCALAPACK :
//           amesosSolver = Teuchos::rcp(new Amesos_Scalapack(*epetraLP));
//           break;
// #endif
// #ifdef HAVE_AMESOS_UMFPACK
//         case Thyra::Amesos2::UMFPACK :
//           amesosSolver = Teuchos::rcp(new Amesos_Umfpack(*epetraLP));
//           break;
// #endif
// #ifdef HAVE_AMESOS_SUPERLUDIST
//         case Thyra::Amesos2::SUPERLUDIST :
//           amesosSolver = Teuchos::rcp(new Amesos_Superludist(*epetraLP));
//           break;
// #endif
// #ifdef HAVE_AMESOS_SUPERLU
//         case Thyra::Amesos2::SUPERLU :
//           amesosSolver = Teuchos::rcp(new Amesos_Superlu(*epetraLP));
//           break;
// #endif
// #ifdef HAVE_AMESOS_DSCPACK
//         case Thyra::Amesos2::DSCPACK :
//           amesosSolver = Teuchos::rcp(new Amesos_Dscpack(*epetraLP));
//           break;
// #endif
// #ifdef HAVE_AMESOS_PARDISO
//         case Thyra::Amesos2::PARDISO :
//           amesosSolver = Teuchos::rcp(new Amesos_Pardiso(*epetraLP));
//           break;
// #endif
// #ifdef HAVE_AMESOS_TAUCS
//         case Thyra::Amesos2::TAUCS :
//           amesosSolver = Teuchos::rcp(new Amesos_Taucs(*epetraLP));
//           break;
// #endif
// #ifdef HAVE_AMESOS_PARAKLETE
//         case Thyra::Amesos2::PARAKLETE :
//           amesosSolver = Teuchos::rcp(new Amesos_Paraklete(*epetraLP));
//           break;
// #endif
//         default:
//           TEUCHOS_TEST_FOR_EXCEPTION(
//             true, std::logic_error
//             ,"Error, the solver type ID = " << solverType_ << " is invalid!"
//             );
//       }
//     }
//     // Set the parameters
//     if(paramList_.get()) amesosSolver->SetParameters(paramList_->sublist("Amesos2 Settings"));
//     // Do the initial factorization
//     {
//       THYRA_FUNC_TIME_MONITOR_DIFF("Stratimikos: AmesosLOWSF:Symbolic", Symbolic);
//       const int err = amesosSolver->SymbolicFactorization();
//       TEUCHOS_TEST_FOR_EXCEPTION( 0!=err, CatastrophicSolveFailure,
//         "Error, SymbolicFactorization() on amesos solver of type \'"<<Teuchos::typeName(*amesosSolver)<<"\'\n"
//         "returned error code "<<err<<"!" );
//     }
//     {
//       THYRA_FUNC_TIME_MONITOR_DIFF("Stratimikos: AmesosLOWSF:Factor", Factor);
//       const int err = amesosSolver->NumericFactorization();
//       TEUCHOS_TEST_FOR_EXCEPTION( 0!=err, CatastrophicSolveFailure,
//         "Error, NumericFactorization() on amesos solver of type \'"<<Teuchos::typeName(*amesosSolver)<<"\'\n"
//         "returned error code "<<err<<"!" );
//     }
//     // Initialize the LOWS object and we are done!
//     amesosOp->initialize(fwdOp,fwdOpSrc,epetraLP,amesosSolver,epetraFwdOpTransp,epetraFwdOpScalar);
//   }
//   else {
//     //
//     // This LOWS object has already be initialized once so we must just reset
//     // the matrix and refactor it.
//     //
//     // Get non-const pointers to the linear problem and the amesos solver.
//     // These const-casts are just fine since the amesosOp in non-const.
//     RCP<Epetra_LinearProblem>
//       epetraLP = Teuchos::rcp_const_cast<Epetra_LinearProblem>(amesosOp->get_epetraLP());
//     RCP<Amesos_BaseSolver>
//       amesosSolver = amesosOp->get_amesosSolver();
//     // Reset the forward operator with memory of RCP to Epetra_Operator view!
//     epetraLP->SetOperator(const_cast<Epetra_Operator*>(&*epetraFwdOp));
//     Teuchos::get_nonconst_extra_data<RCP<const Epetra_Operator> >(epetraLP,epetraFwdOp_str) = epetraFwdOp;
//     // Reset the parameters
//     if(paramList_.get()) amesosSolver->SetParameters(paramList_->sublist(Amesos_Settings_name));
//     // Repivot if asked
//     if(refactorizationPolicy_==Amesos2::REPIVOT_ON_REFACTORIZATION) {
//       THYRA_FUNC_TIME_MONITOR_DIFF("Stratimikos: AmesosLOWSF:Symbolic", Symbolic);
//       const int err = amesosSolver->SymbolicFactorization();
//       TEUCHOS_TEST_FOR_EXCEPTION( 0!=err, CatastrophicSolveFailure,
//         "Error, SymbolicFactorization() on amesos solver of type \'"<<Teuchos::typeName(*amesosSolver)<<"\'\n"
//         "returned error code "<<err<<"!" );
//     }
//     {
//       THYRA_FUNC_TIME_MONITOR_DIFF("Stratimikos: AmesosLOWSF::Factor", Factor);
//       const int err = amesosSolver->NumericFactorization();
//       TEUCHOS_TEST_FOR_EXCEPTION( 0!=err, CatastrophicSolveFailure,
//         "Error, NumericFactorization() on amesos solver of type \'"<<Teuchos::typeName(*amesosSolver)<<"\'\n"
//         "returned error code "<<err<<"!" );
//     }
//     /* ToDo: Put this back in once PrintStatus accepts an std::ostream!
//     OsTab tab(out);
//     amesosSolver->PrintStatus()
//     */
//     // Reinitialize the LOWS object and we are done! (we must do this to get the
//     // possibly new transpose and scaling factors back in)
//     amesosOp->initialize(fwdOp,fwdOpSrc,epetraLP,amesosSolver,epetraFwdOpTransp,epetraFwdOpScalar);
//   }
//   amesosOp->setOStream(this->getOStream());
//   amesosOp->setVerbLevel(this->getVerbLevel());
}

template<typename Scalar>
bool Amesos2LinearOpWithSolveFactory<Scalar>::supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const
{
  return false;
}

template<typename Scalar>
void Amesos2LinearOpWithSolveFactory<Scalar>::initializePreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> >       &fwdOpSrc
  ,const RCP<const PreconditionerBase<Scalar> >      &prec
  ,LinearOpWithSolveBase<Scalar>                                      *Op
  ,const ESupportSolveUse                                             supportSolveUse
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->throwOnPrecInput_, std::logic_error
    ,"Error, the concrete implementation described as \'"<<this->description()<<"\' does not support preconditioners "
    "and has been configured to throw this exception when the  initializePreconditionedOp(...) function is called!"
    );
  // this->initializeOp(fwdOpSrc,Op,supportSolveUse); // Ignore the preconditioner!
}

template<typename Scalar>
void Amesos2LinearOpWithSolveFactory<Scalar>::initializePreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> >       &fwdOpSrc
  ,const RCP<const LinearOpSourceBase<Scalar> >      &approxFwdOpSrc
  ,LinearOpWithSolveBase<Scalar>                                      *Op
  ,const ESupportSolveUse                                             supportSolveUse
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->throwOnPrecInput_, std::logic_error
    ,"Error, the concrete implementation described as \'"<<this->description()<<"\' does not support preconditioners "
    "and has been configured to throw this exception when the  initializePreconditionedOp(...) function is called!"
    );
  // this->initializeOp(fwdOpSrc,Op,supportSolveUse); // Ignore the preconditioner!
}

template<typename Scalar>
void Amesos2LinearOpWithSolveFactory<Scalar>::uninitializeOp(
  LinearOpWithSolveBase<Scalar>                               *Op
  ,RCP<const LinearOpSourceBase<Scalar> >    *fwdOpSrc
  ,RCP<const PreconditionerBase<Scalar> >    *prec
  ,RCP<const LinearOpSourceBase<Scalar> >    *approxFwdOpSrc
  ,ESupportSolveUse                                           *supportSolveUse
  ) const
{
// #ifdef TEUCHOS_DEBUG
//   TEUCHOS_TEST_FOR_EXCEPT(Op==NULL);
// #endif
//   Amesos2LinearOpWithSolve
//     *amesosOp = &Teuchos::dyn_cast<Amesos2LinearOpWithSolve>(*Op);
//   RCP<const LinearOpSourceBase<double> >
//     _fwdOpSrc = amesosOp->extract_fwdOpSrc(); // Will be null if uninitialized!
//   if(_fwdOpSrc.get()) {
//     // Erase the Epetra_Operator view of the forward operator!
//     RCP<Epetra_LinearProblem> epetraLP = amesosOp->get_epetraLP();
//     Teuchos::get_nonconst_extra_data< RCP<const Epetra_Operator> >(
//       epetraLP,epetraFwdOp_str
//       )
//       = Teuchos::null;
//     // Note, we did not erase the address of the operator in
//     // epetraLP->GetOperator() since it seems that the amesos solvers do not
//     // recheck the value of GetProblem()->GetOperator() so you had better not
//     // rest this!
//   }
//   if(fwdOpSrc) *fwdOpSrc = _fwdOpSrc; // It is fine if the client does not want this object back!
//   if(prec) *prec = Teuchos::null; // We never keep a preconditioner!
//   if(approxFwdOpSrc) *approxFwdOpSrc = Teuchos::null; // We never keep an approximate fwd operator!
}

// Overridden from ParameterListAcceptor

template<typename Scalar>
void Amesos2LinearOpWithSolveFactory<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  // TEUCHOS_TEST_FOR_EXCEPT(paramList.get()==NULL);
  // paramList->validateParameters(*this->getValidParameters(),0); // Only validate this level for now!
  // paramList_ = paramList;
  // solverType_ =
  //   Amesos2::solverTypeNameToEnumMap.get<Amesos2::ESolverType>(
  //     paramList_->get(
  //       SolverType_name
  //       ,Amesos2::toString(solverType_)
  //       )
  //     ,paramList_->name()+"->"+SolverType_name
  //     );
  // refactorizationPolicy_ = 
  //   Amesos2::refactorizationPolicyNameToEnumMap.get<Amesos2::ERefactorizationPolicy>(
  //     paramList_->get(
  //       RefactorizationPolicy_name
  //       ,Amesos2::toString(refactorizationPolicy_)
  //       )
  //     ,paramList_->name()+"->"+RefactorizationPolicy_name
  //     );
  // throwOnPrecInput_ = paramList_->get(ThrowOnPreconditionerInput_name,throwOnPrecInput_);
  // Teuchos::readVerboseObjectSublist(&*paramList_,this);
}

template<typename Scalar>
RCP<Teuchos::ParameterList>
Amesos2LinearOpWithSolveFactory<Scalar>::getNonconstParameterList()
{
  return paramList_;
}

template<typename Scalar>
RCP<Teuchos::ParameterList>
Amesos2LinearOpWithSolveFactory<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

template<typename Scalar>
RCP<const Teuchos::ParameterList>
Amesos2LinearOpWithSolveFactory<Scalar>::getParameterList() const
{
  return paramList_;
}

template<typename Scalar>
RCP<const Teuchos::ParameterList>
Amesos2LinearOpWithSolveFactory<Scalar>::getValidParameters() const
{
  return generateAndGetValidParameters();
}

// Public functions overridden from Teuchos::Describable

template<typename Scalar>
std::string Amesos2LinearOpWithSolveFactory<Scalar>::description() const
{
  std::ostringstream oss;
  oss << "Thyra::Amesos2LinearOpWithSolveFactory{";
  oss << "solverType=" << toString(solverType_);
  oss << "}";
  return oss.str();
}

// private

template<typename Scalar>
RCP<const Teuchos::ParameterList>
Amesos2LinearOpWithSolveFactory<Scalar>::generateAndGetValidParameters()
{
  static RCP<Teuchos::ParameterList> validParamList;
  if(validParamList.get()==NULL) {
    validParamList = Teuchos::rcp(new Teuchos::ParameterList("Amesos2"));
    validParamList->set(
      SolverType_name
      ,Amesos2::toString(Amesos2::KLU)
      );
    validParamList->set(RefactorizationPolicy_name,Amesos2::toString(Amesos2::REPIVOT_ON_REFACTORIZATION));
    validParamList->set(ThrowOnPreconditionerInput_name,bool(true));
    // validParamList->sublist(Amesos_Settings_name).setParameters(::Amesos::GetValidParameters());
    Teuchos::setupVerboseObjectSublist(&*validParamList);
  }
  return validParamList;
}

} // namespace Thyra

#endif // SUN_CXX
