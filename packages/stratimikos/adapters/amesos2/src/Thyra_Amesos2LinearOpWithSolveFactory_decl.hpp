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

#ifndef THYRA_AMESOS2_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
#define THYRA_AMESOS2_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_Amesos2Types.hpp"
#include "Amesos2_DataTypes.hpp"
//#include "Amesos_BaseSolver.h"
#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
//#include "Thyra_EpetraOperatorViewExtractorBase.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Trilinos_Details_LinearSolver.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"

namespace Thyra {

/** \brief Concrete <tt>LinearOpWithSolveFactoryBase</tt> adapter subclass that uses
 * Amesos direct solvers.
 *
 * This class creates objects of type <tt>AmesosLinearOpWithSolve</tt>
 * (through the <tt>LinearOpWithSolveBase</tt> interface) which can then be
 * used to solve for linear systems.  The <tt>%AmesosLinearOpWithSolve</tt>
 * objects created and initialized by this object are completely independent
 * from <tt>*this</tt>.  This allows for multiple
 * <tt>%AmesosLinearOpWithSolve</tt> objects to be created and maintained
 * simultaneously and for <tt>*this</tt> factory object to be destroyed
 * without affecting the created <tt>%AmesosLinearOpWithSolve</tt> objects.
 *
 * ToDo: Mention parameter list usage.
 *
 * <b>Development notes:</b> This class has been designed to allow for "smart"
 * <tt>EpetraLinearOpBase</tt> subclasses that can create an
 * <tt>Epetra_Operator</tt> view on command where the underlying storage may
 * not be an <tt>Epetra</tt> object.  However, the current implementation of
 * at least some of the <tt>Amesos_BaseSolver</tt> subclasses do not allow the
 * <tt>%Epetra_Operator</tt> object to change after construction.  Therefore,
 * this current flawed implementation requires that every call to the
 * <tt>EpetraLinearOpBase::getEpetraOpView()</tt> function return the same
 * <tt>%Epetra_Operator</tt> object.
 *
 * \ingroup Amesos_Thyra_adapters_grp
 */
  template<typename Scalar>
  class Amesos2LinearOpWithSolveFactory : public LinearOpWithSolveFactoryBase<Scalar> {
  public:

  /** \name Parameter names for Paramter List */
  //@{

  /** \brief . */
  static const std::string SolverType_name;
  /** \brief . */
  static const std::string RefactorizationPolicy_name;
  /** \brief . */
  static const std::string ThrowOnPreconditionerInput_name;
  /** \brief . */
  static const std::string Amesos_Settings_name;

  //@}

  // typedef Tpetra::CrsMatrix<Scalar> MAT;
  // typedef Tpetra::MultiVector<Scalar> MV;
  // typedef Tpetra::Operator<Scalar> OP;
  // typedef long long int GO;
  // typedef int LO;

  /** \name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  ~Amesos2LinearOpWithSolveFactory();

  /** \brief Constructor which sets the defaults.
   */
  Amesos2LinearOpWithSolveFactory(
    const Amesos2Stratimikos::ESolverType                 solverType  = Amesos2Stratimikos::KLU
    ,const Amesos2Stratimikos::ERefactorizationPolicy     refactorizationPolicy  = Amesos2Stratimikos::REPIVOT_ON_REFACTORIZATION
    ,const bool                               throwOnPrecInput       = true
    );
    
  /** \brief Set the strategy object used to extract an
   * <tt>Epetra_Operator</tt> view of an input forward operator.
   *
   * This view will then be dynamically casted to <tt>Epetra_RowMatrix</tt>
   * before it is used.
   *
   * The default implementation used is <tt>EpetraOperatorViewExtractorBase</tt>.
   */
  //STANDARD_COMPOSITION_MEMBERS( EpetraOperatorViewExtractorBase, epetraFwdOpViewExtractor );

  //@}

  /** @name Overridden public functions from LinearOpWithSolveFactoryBase */
  //@{

  // This will be used to extract the Tpetra operator that later will be cast into a Tpetra_CrsMatrix
  //typedef Thyra::TpetraOperatorVectorExtraction<Scalar> converterT;

  /** \brief Returns true if <tt>dynamic_cast<const
   * EpetraLinearOpBase*>(fwdOpSrc)!=NULL</tt> .
   */
  bool isCompatible( const LinearOpSourceBase<Scalar> &fwdOpSrc ) const;

  /** \brief . */
  Teuchos::RCP<LinearOpWithSolveBase<Scalar> > createOp() const;

  /** \brief . */
  void initializeOp(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> >    &fwdOpSrc
    ,LinearOpWithSolveBase<Scalar>                                   *Op
    ,const ESupportSolveUse                                          supportSolveUse
    ) const;

  /** \brief Returns <tt>false</tt> . */
  bool supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const;

  /** \brief Throws exception if <tt>this->throwOnPrecInput()==true</tt> and
   * calls <tt>this->initializeOp(fwdOpSrc,Op)</tt> otherwise
   */
  void initializePreconditionedOp(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> >       &fwdOpSrc
    ,const Teuchos::RCP<const PreconditionerBase<Scalar> >      &prec
    ,LinearOpWithSolveBase<Scalar>                                      *Op
    ,const ESupportSolveUse                                             supportSolveUse
    ) const;

  /** \brief Throws exception if <tt>this->throwOnPrecInput()==true</tt> and
   * calls <tt>this->initializeOp(fwdOpSrc,Op)</tt> otherwise
   */
  void initializePreconditionedOp(
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> >       &fwdOpSrc
    ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >      &approxFwdOpSrc
    ,LinearOpWithSolveBase<Scalar>                                      *Op
    ,const ESupportSolveUse                                             supportSolveUse
    ) const;

  /** \brief . */
  void uninitializeOp(
    LinearOpWithSolveBase<Scalar>                               *Op
    ,Teuchos::RCP<const LinearOpSourceBase<Scalar> >    *fwdOpSrc
    ,Teuchos::RCP<const PreconditionerBase<Scalar> >    *prec
    ,Teuchos::RCP<const LinearOpSourceBase<Scalar> >    *approxFwdOpSrc
    ,ESupportSolveUse                                           *supportSolveUse
    ) const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  // /////////////////////////
  // Private data members

  Amesos2Stratimikos::ESolverType                             solverType_;
  Amesos2Stratimikos::ERefactorizationPolicy                  refactorizationPolicy_;
  bool                                            throwOnPrecInput_;
  Teuchos::RCP<Teuchos::ParameterList>    paramList_;

  // /////////////////////////
  // Private member functions

  static Teuchos::RCP<const Teuchos::ParameterList>
  generateAndGetValidParameters();

};

} // namespace Thyra

#endif // THYRA_AMESOS2_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
