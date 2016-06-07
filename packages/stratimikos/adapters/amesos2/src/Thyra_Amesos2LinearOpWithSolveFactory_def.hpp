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
#include "Amesos2.hpp"
#include "Amesos2_Details_LinearSolverFactory.hpp"
#include "Amesos2_Version.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_TypeTraits.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "Amesos2_Superludist.hpp"

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
  //std::cout << "Amesos2LinearOpWithSolveFactory<Scalar>::~Amesos2LinearOpWithSolveFactory()" << std::endl;

}

template<typename Scalar>
Amesos2LinearOpWithSolveFactory<Scalar>::Amesos2LinearOpWithSolveFactory(
  const Amesos2Stratimikos::ESolverType                            solverType
  ,const Amesos2Stratimikos::ERefactorizationPolicy                refactorizationPolicy
  ,const bool                                          throwOnPrecInput
    )
  ://epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd()))
  //,solverType_(solverType)
  solverType_(solverType)
  ,refactorizationPolicy_(refactorizationPolicy)
  ,throwOnPrecInput_(throwOnPrecInput)
{
  //std::cout << "Amesos2LinearOpWithSolveFactory<Scalar>::Amesos2LinearOpWithSolveFactory(...)" << std::endl;
}

// Overridden from LinearOpWithSolveFactoryBase

template<typename Scalar>
bool Amesos2LinearOpWithSolveFactory<Scalar>::isCompatible(
  const LinearOpSourceBase<Scalar> &fwdOpSrc
  ) const
{
  //std::cout << " Amesos2LinearOpWithSolveFactory<Scalar>::isCompatible" << std::endl;

  Teuchos::RCP<const LinearOpBase<Scalar> >
    fwdOp = fwdOpSrc.getOp();
  Teuchos::RCP< const Tpetra_Operator > tpetraFwdOp;
  tpetraFwdOp = ConverterT::getConstTpetraOperator(fwdOp);

  if ( ! dynamic_cast<const Tpetra_CrsMatrix * >(&*tpetraFwdOp) )
    return false;
  return true;
}

template<typename Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
Amesos2LinearOpWithSolveFactory<Scalar>::createOp() const
{
  //std::cout << "Amesos2LinearOpWithSolveFactory<Scalar>::createOp()" << std::endl;
  return Teuchos::rcp(new Amesos2LinearOpWithSolve<Scalar>());
}

template<typename Scalar>
void Amesos2LinearOpWithSolveFactory<Scalar>::initializeOp(
  const RCP<const LinearOpSourceBase<Scalar> >    &fwdOpSrc
  ,LinearOpWithSolveBase<Scalar>                                   *Op
  ,const ESupportSolveUse                                          supportSolveUse
  ) const
{
  //std::cout << "Amesos2LinearOpWithSolveFactory<Scalar>::initializeOp" << std::endl;
   THYRA_FUNC_TIME_MONITOR("Stratimikos: Amesos2LOWSF");
#ifdef TEUCHOS_DEBUG
   TEUCHOS_TEST_FOR_EXCEPT(Op==NULL);
#endif
    
  TEUCHOS_TEST_FOR_EXCEPT(Op==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(fwdOpSrc.get()==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(fwdOpSrc->getOp().get()==NULL);
  RCP<const LinearOpBase<Scalar> > fwdOp = fwdOpSrc->getOp();

  Teuchos::RCP<Teuchos::FancyOStream> 
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  //std::cout << "Unwrap TpetraOperator" << std::endl;
  //
  // Unwrap and get the forward Tpetra::Operator object
  //
  RCP< const Tpetra_Operator > tpetraFwdOp;
  //std::cout << "Tpetra operator" << std::endl;
  tpetraFwdOp = ConverterT::getConstTpetraOperator(fwdOp);
  //std::cout << "amesos operator" << std::endl;
  // Get the Amesos2LinearOpWithSolve object
  Amesos2LinearOpWithSolve<Scalar>
    *amesosOp = &Teuchos::dyn_cast<Amesos2LinearOpWithSolve<Scalar>>(*Op);

  //
  // Determine if we must start over or not
  //
  bool startOver = ( amesosOp->get_amesosSolver()==Teuchos::null );
  if(!startOver) 
    {
      //std::cout << "Start over 2" << std::endl;
      RCP< const Tpetra_Operator > tpetraOp = amesosOp->get_amesosSolver()->getMatrix();
      startOver =
      	(
      	  tpetraFwdOp.get() != tpetraOp.get()
	  // We must start over if the matrix object changes.  This is a
	  // weakness of Amesos but there is nothing I can do about this right
	  // now!
	);
    }
  //
  // Update the amesos solver
  //
  if(startOver) {
    //std::cout << "Start over" << std::endl;
    //
    // This LOWS object has not be initialized yet or is not compatible with the existing
    // 
    // so this is where we setup everything from the ground up.
    //
    // Create the linear problem factory
    Amesos2::Details::LinearSolverFactory<Tpetra_MultiVector,Tpetra_Operator,Scalar> linearsolverfactory;
      
    // Create the concrete solver
    Teuchos::RCP< Trilinos::Details::LinearSolver<Tpetra_MultiVector,Tpetra_Operator,Scalar> > 
      amesosSolver = linearsolverfactory.getLinearSolver("superludist");

    // set 
    amesosSolver->setMatrix(tpetraFwdOp);
    
    // Do the initial factorization
    *out << "Performing factorization ...\n";
    amesosSolver->symbolic();
    amesosSolver->numeric();
    
    
    // Initialize the LOWS object and we are done!
    amesosOp->initialize(fwdOp,fwdOpSrc,amesosSolver);
  }
  else {
    //std::cout << "Start over 3" << std::endl;
    //
    // This LOWS object has already be initialized once so we must just reset
    // the matrix and refactor it.
    //
    // Get non-const pointers to the linear problem and the amesos solver.
    // These const-casts are just fine since the amesosOp in non-const.
    Teuchos::RCP< Trilinos::Details::LinearSolver<Tpetra_MultiVector,Tpetra_Operator,Scalar> >
      amesosSolver = amesosOp->get_amesosSolver();
    
    // set 
    amesosSolver->setMatrix(tpetraFwdOp);

    // Do the initial factorization
    *out << "Performing factorization ...\n";
    amesosSolver->symbolic();
    amesosSolver->numeric();

    // Initialize the LOWS object and we are done!
    amesosOp->initialize(fwdOp,fwdOpSrc,amesosSolver);

  }
  amesosOp->setOStream(this->getOStream());
  amesosOp->setVerbLevel(this->getVerbLevel());
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
  this->initializeOp(fwdOpSrc,Op,supportSolveUse); // Ignore the preconditioner!
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
  this->initializeOp(fwdOpSrc,Op,supportSolveUse); // Ignore the preconditioner!
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
  //std::cout << "Amesos2LinearOpWithSolveFactory<Scalar>::uninitializeOp" << std::endl;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(Op==NULL);
#endif
  Amesos2LinearOpWithSolve<Scalar>
    *amesosOp = &Teuchos::dyn_cast<Amesos2LinearOpWithSolve<Scalar>>(*Op);
  RCP<const LinearOpSourceBase<double> >
    _fwdOpSrc = amesosOp->extract_fwdOpSrc(); // Will be null if uninitialized!
  if(fwdOpSrc) *fwdOpSrc = _fwdOpSrc; // It is fine if the client does not want this object back!
  if(prec) *prec = Teuchos::null; // We never keep a preconditioner!
  if(approxFwdOpSrc) *approxFwdOpSrc = Teuchos::null; // We never keep an approximate fwd operator!
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
      ,Amesos2Stratimikos::toString(Amesos2Stratimikos::SuperLU)
      );
    validParamList->set(RefactorizationPolicy_name,Amesos2Stratimikos::toString(Amesos2Stratimikos::REPIVOT_ON_REFACTORIZATION));
    validParamList->set(ThrowOnPreconditionerInput_name,bool(true));
    // validParamList->sublist(Amesos_Settings_name).setParameters(::Amesos2Stratimikos::GetValidParameters());
    Teuchos::setupVerboseObjectSublist(&*validParamList);
  }
  return validParamList;
}

} // namespace Thyra

#endif // SUN_CXX
