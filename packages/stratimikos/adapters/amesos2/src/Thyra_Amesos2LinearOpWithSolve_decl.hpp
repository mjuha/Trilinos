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

#ifndef THYRA_AMESOS2_LINEAR_OP_WITH_SOLVE_DECL_HPP
#define THYRA_AMESOS2_LINEAR_OP_WITH_SOLVE_DECL_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

#include "Amesos2.hpp"
#include "Amesos2_Details_LinearSolverFactory.hpp"
#include "Amesos2_Version.hpp"
#include "Trilinos_Details_LinearSolver.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"
//#include "Thyra_EpetraLinearOpBase.hpp"
//#include "Epetra_LinearProblem.h"
//#include "Amesos_BaseSolver.h"


namespace Thyra {


/** \brief Concrete <tt>LinearOpWithSolveBase</tt> subclass that adapts any
 * <tt>Amesos_BaseSolver</tt> object.
 *
 * See the <tt>LinearOpWithSolveBase</tt> interface for a description of how
 * to use objects of this type.
 *
 * <b>Note:</b> Clients should not generally directly create objects of this
 * type but instead should use <tt>AmesosLinearOpWithSolveFactory</tt>.  Only
 * very sophisticated users should ever directly interact with an object
 * through this subclass interface.
 *
 * \ingroup Amesos_Thyra_adapters_grp
 */

template<typename Scalar>
class Amesos2LinearOpWithSolve: virtual public LinearOpWithSolveBase<Scalar>
{
public:

  typedef int LO;
  typedef int GO;
  typedef typename Tpetra::Vector<Scalar>::node_type NT;
  
  typedef Tpetra::CrsMatrix<Scalar> MAT;
  typedef Tpetra::MultiVector<Scalar> MV;
  typedef Tpetra::Operator<Scalar> OP;

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  Amesos2LinearOpWithSolve();

  /** \brief Calls <tt>this->initialize()</tt>. */
  Amesos2LinearOpWithSolve(
    const Teuchos::RCP<const LinearOpBase<Scalar> > &fwdOp,
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const Teuchos::RCP<Amesos2::Details::LinearSolverFactory<MV,OP,Scalar>> &tpetraLP,
    const Teuchos::RCP< Trilinos::Details::LinearSolver< MV, OP, Scalar > > &amesosSolver,
    // const Teuchos::RCP<Epetra_LinearProblem> &epetraLP,
    // const Teuchos::RCP<Amesos_BaseSolver> &amesosSolver,
    const EOpTransp amesosSolverTransp,
    const Scalar amesosSolverScalar
    );

  /** \brief First initialization.
   *
   * \param  fwdOp
   *           [in] The forward operator for which the factorization exists.
   * \param  epetraLP
   *           [in] The <tt>Epetra_LinearProblem</tt> object that was used to
   *           create the <tt>Amesos_BaseSolver</tt> object
   *           <tt>*amesosSolver</tt>.  Note that the RHS and the LHS
   *           multi-vector pointers in this object will be set and unset
   *           here.
   * \param  amesosSolver
   *           [in] Contains the factored, and ready to go,
   *           <tt>Amesos_BaseSolver</tt> object ready to solve linear system.
   * \param  amesosSolverTransp
   *           [in] Determines if the %Amesos solver should be used as its
   *           transpose or not.
   * \param  amesosSolverScalar
   *           [in] Determines the scaling factor associated with the %Amesos
   *           solver.  The solution to the linear solve is scaled by
   *           <tt>1/amesosSolverScalar</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>fwdOp.get()!=NULL</tt>
   * <li><tt>epetraLP.get()!=NULL</tt>
   * <li><tt>amesosSolver.get()!=NULL</tt>
   * <li><tt>*epetraLP->GetOperator()</tt> is compatible with <tt>*fwdOp</tt>
   * <li><tt>epetraLP->GetLHS()==NULL</tt>
   * <li><tt>epetraLP->GetRHS()==NULL</tt>
   * <li><tt>*amesosSolver</tt> contains the factorization of <tt>*fwdOp</tt> and is
   *     ready to solve linear systems!
   * </ul>
   * 
   * <b>Postconditions:</b><ul>
   * <li><tt>this->get_fwdOp().get() == fwdOp.get()</tt>
   * <li><tt>this->get_epetraLP().get() == epetraLP.get()</tt>
   * <li><tt>this->get_amesosSolver().get() == amesosSolver.get()</tt>
   * <li><tt>this->get_amesosSolverTransp() == amesosSolverTransp</tt>
   * <li><tt>this->get_amesosSolverScalar() == amesosSolverScalar</tt>
   * </ul>
   */
  void initialize(
    const Teuchos::RCP<const LinearOpBase<Scalar> > &fwdOp,
    const Teuchos::RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const Teuchos::RCP<Amesos2::Details::LinearSolverFactory<MV,OP,Scalar>> &tpetraLP,
    const Teuchos::RCP< Trilinos::Details::LinearSolver< MV, OP, Scalar > > &amesosSolver,
    // const Teuchos::RCP<Epetra_LinearProblem> &epetraLP,
    // const Teuchos::RCP<Amesos_BaseSolver> &amesosSolver,
    const EOpTransp amesosSolverTransp,
    const Scalar amesosSolverScalar
    );

  /** \brief Extract the <tt>LinearOpSourceBase<double></tt> object so that it can be modified.
   * 
   * <b>Postconditions:</b><ul>
   * <li><tt>return.get()</tt> is the same as <tt>this->get_fwdOpSrc().get()</tt> before call.
   * <li><tt><tt>this->get_fwdOpSrc().get()==NULL</tt>
   * </ul>
   */
  Teuchos::RCP<const LinearOpSourceBase<Scalar> > extract_fwdOpSrc();

  /** \brief . */
  Teuchos::RCP<const LinearOpBase<Scalar> > get_fwdOp() const;

  /** \brief . */
  Teuchos::RCP<const LinearOpSourceBase<Scalar> > get_fwdOpSrc() const;

  // /** \brief . */
  Teuchos::RCP< Amesos2::Details::LinearSolverFactory<MV,OP,Scalar> > 
  get_linearSolverFactory() const;
  // Teuchos::RCP<Epetra_LinearProblem> get_epetraLP() const;

  // /** \brief . */
  Teuchos::RCP< Trilinos::Details::LinearSolver< MV, OP, Scalar > >
  get_amesosSolver() const;
  // Teuchos::RCP<Amesos_BaseSolver> get_amesosSolver() const;

  /** \brief . */
  EOpTransp get_amesosSolverTransp() const;

  /** \brief . */
  Scalar get_amesosSolverScalar() const;

  /** \brief Uninitialize.
   */
  void uninitialize(
    Teuchos::RCP<const LinearOpBase<Scalar> > *fwdOp = NULL,
    Teuchos::RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc = NULL,
    Teuchos::RCP<Amesos2::Details::LinearSolverFactory<MV,OP,Scalar>> *tpetraLP = NULL,
    Teuchos::RCP< Trilinos::Details::LinearSolver< MV, OP, Scalar > > *amesosSolver = NULL,
    // Teuchos::RCP<Epetra_LinearProblem> *epetraLP = NULL,
    // Teuchos::RCP<Amesos_BaseSolver> *amesosSolver = NULL,
    EOpTransp *amesosSolverTransp = NULL,
    Scalar *amesosSolverScalar = NULL
    );
  
  //@}

  /** @name Overridden public functions from LinearOpBase */
  //@{
  /** \brief. */
  Teuchos::RCP< const VectorSpaceBase<Scalar> > range() const;
  /** \brief. */
  Teuchos::RCP< const VectorSpaceBase<Scalar> > domain() const;
  /** \brief. */
  Teuchos::RCP<const LinearOpBase<Scalar> > clone() const;
  //@}

  /** @name Overridden public functions from Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;
  //@}

protected:

  /** @name Overridden from LinearOpBase  */
  //@{
  /** \brief . */
  virtual bool opSupportedImpl(EOpTransp M_trans) const;
  /** \brief . */
  virtual void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const;
  //@}

  /** @name Overridden from LinearOpWithSolveBase. */
  //@{
  /** \brief . */
  virtual bool solveSupportsImpl(EOpTransp M_trans) const;
  /** \brief . */
  virtual bool solveSupportsSolveMeasureTypeImpl(
    EOpTransp M_trans, const SolveMeasureType& solveMeasureType
    ) const;
  /** \brief . */
  SolveStatus<Scalar> solveImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &B,
    const Ptr<MultiVectorBase<Scalar> > &X,
    const Ptr<const SolveCriteria<Scalar> > solveCriteria
    ) const;
  //@}

private:

  Teuchos::RCP<const LinearOpBase<Scalar> > fwdOp_;
  Teuchos::RCP<const LinearOpSourceBase<Scalar> > fwdOpSrc_;
  Teuchos::RCP< Amesos2::Details::LinearSolverFactory<MV,OP,Scalar> > linearsolverfactory_;
  Teuchos::RCP< Trilinos::Details::LinearSolver<MV,OP,Scalar> > solver_;
  // Teuchos::RCP<Epetra_LinearProblem> epetraLP_;
  // Teuchos::RCP<Amesos_BaseSolver> amesosSolver_;
  EOpTransp amesosSolverTransp_;
  Scalar amesosSolverScalar_;

  void assertInitialized() const;

};

// ///////////////////////////
// Inline members

template<typename Scalar>
inline
Teuchos::RCP<const LinearOpBase<Scalar> >
Amesos2LinearOpWithSolve<Scalar>::get_fwdOp() const
{
  return fwdOp_;
}

template<typename Scalar>
inline
Teuchos::RCP<const LinearOpSourceBase<Scalar> >
Amesos2LinearOpWithSolve<Scalar>::get_fwdOpSrc() const
{
  return fwdOpSrc_;
}

template<typename Scalar>
inline
Teuchos::RCP< Amesos2::Details::LinearSolverFactory< Tpetra::MultiVector<Scalar>, Tpetra::Operator<Scalar>, Scalar > >
Amesos2LinearOpWithSolve<Scalar>::get_linearSolverFactory() const
{
  return linearsolverfactory_;
}

template<typename Scalar>
inline
Teuchos::RCP< Trilinos::Details::LinearSolver< Tpetra::MultiVector<Scalar>, Tpetra::Operator<Scalar>, Scalar > >
Amesos2LinearOpWithSolve<Scalar>::get_amesosSolver() const
{
  return solver_;
}

// template<typename Scalar>
// inline
// Teuchos::RCP< Amesos2::Details::LinearSolverFactory<MV,OP,MT> > 
// Amesos2LinearOpWithSolve<Scalar>::get_linearSolverFactory() const
// {
//   return linearsolverfactory;
//}

// template<typename Scalar>
// inline
// Teuchos::RCP<Amesos_BaseSolver>
// Amesos2LinearOpWithSolve<Scalar>::get_amesosSolver() const
// {
//   return amesosSolver_;
// }

template<typename Scalar>
inline
EOpTransp Amesos2LinearOpWithSolve<Scalar>::get_amesosSolverTransp() const
{
  return amesosSolverTransp_;
}

template<typename Scalar>
inline
Scalar Amesos2LinearOpWithSolve<Scalar>::get_amesosSolverScalar() const
{
  return amesosSolverScalar_;
}

} // namespace Thyra

#endif	// THYRA_AMESOS_LINEAR_OP_WITH_SOLVE_HPP
