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

#ifndef THYRA_AMESOS2_TYPES_HPP
#define THYRA_AMESOS2_TYPES_HPP

//#include "Amesos_ConfigDefs.h"
#include "Amesos2_config.h"
#include "Teuchos_StringToIntMap.hpp"

namespace Thyra {

namespace Amesos2Stratimikos {

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
  enum ESolverType {SuperLU};

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
const int numSolverTypes = 1;

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
extern const ESolverType solverTypeValues[numSolverTypes];

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
extern const char* solverTypeNames[numSolverTypes];

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
extern const bool supportsUnsymmetric[numSolverTypes];

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
inline const char* toString(const ESolverType solverType)
{ return solverTypeNames[solverType]; }

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
extern Teuchos::StringToIntMap solverTypeNameToEnumMap;

/** \brief The policy used on refactoring a matrix.
\ingroup Amesos_Thyra_adapters_grp
*/
enum ERefactorizationPolicy {
  REPIVOT_ON_REFACTORIZATION     ///< Completely new pivoting will be used on refactorizations!
  ,NO_PIVOT_ON_REFACTORIZATION   ///< No piviting, or only minor repivoting, will be used on refactorizations!
};

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
const int numRefactorizationPolices = 2;

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
extern const ERefactorizationPolicy refactorizationPolicyValues[numRefactorizationPolices];

/** \brief . */
extern const char* refactorizationPolicyNames[numRefactorizationPolices];

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
inline const char* toString(const ERefactorizationPolicy refactorizationPolicy)
{ return refactorizationPolicyNames[refactorizationPolicy]; }

/** \brief .
\ingroup Amesos_Thyra_adapters_grp
*/
extern Teuchos::StringToIntMap refactorizationPolicyNameToEnumMap;

}// namespace Amesos2

} // namespace Thyra

#endif // THYRA_AMESOS2_TYPES_HPP
