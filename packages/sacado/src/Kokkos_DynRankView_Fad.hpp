// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_DYN_RANK_VIEW_SACADO_FAD_HPP
#define KOKKOS_DYN_RANK_VIEW_SACADO_FAD_HPP

#include "Sacado_ConfigDefs.h"

// This file is setup to always work even when KokkosContainers (which contains
// Kokkos::DynRankView) isn't enabled.
//
// We also include Kokkos_DynRankView after our specializations to ensure any
// overloads we provide here are in scope before they are used inside the
// DynRankView

#if defined(HAVE_SACADO_KOKKOSCONTAINERS)

#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "Kokkos_View_Fad.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

// Forward declaration of the trait we are specialization
// (see comments above about order of includes)
template <typename Spec> struct DynRankDimTraits;

template <>
struct DynRankDimTraits<ViewSpecializeSacadoFad> {

  enum : size_t{unspecified = ~size_t(0)};

  // Compute the rank of the view from the nonzero dimension arguments.
  // For views of Fad, the rank is one less than the rank determined by the nonzero dimension args
  static size_t computeRank( const size_t N0
                           , const size_t N1
                           , const size_t N2
                           , const size_t N3
                           , const size_t N4
                           , const size_t N5
                           , const size_t N6
                           , const size_t N7 )
  {
    return
      (   (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified && N3 == unspecified && N2 == unspecified && N1 == unspecified && N0 == unspecified) ? 0
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified && N3 == unspecified && N2 == unspecified && N1 == unspecified) ? 0
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified && N3 == unspecified && N2 == unspecified) ? 1
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified && N3 == unspecified) ? 2
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified && N4 == unspecified) ? 3
      : ( (N7 == unspecified && N6 == unspecified && N5 == unspecified) ? 4
      : ( (N7 == unspecified && N6 == unspecified) ? 5
      : ( (N7 == unspecified) ? 6
      : 7 ) ) ) ) ) ) ) );
  }

  // Compute the rank of the view from the nonzero layout arguments.
  template <typename Layout>
  static size_t computeRank( const Layout& layout )
  {
    return computeRank( layout.dimension[0]
                      , layout.dimension[1]
                      , layout.dimension[2]
                      , layout.dimension[3]
                      , layout.dimension[4]
                      , layout.dimension[5]
                      , layout.dimension[6]
                      , layout.dimension[7] );
  }

  // Create the layout for the rank-7 view.
  // For Fad we have to move the fad dimension to the last (rank 8 since the DynRankView is rank-7)
  // LayoutLeft or LayoutRight
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if< (std::is_same<Layout , Kokkos::LayoutRight>::value || std::is_same<Layout , Kokkos::LayoutLeft>::value) , Layout >::type createLayout( const Layout& layout )
  {
    Layout l( layout.dimension[0] != unspecified ? layout.dimension[0] : 1
            , layout.dimension[1] != unspecified ? layout.dimension[1] : 1
            , layout.dimension[2] != unspecified ? layout.dimension[2] : 1
            , layout.dimension[3] != unspecified ? layout.dimension[3] : 1
            , layout.dimension[4] != unspecified ? layout.dimension[4] : 1
            , layout.dimension[5] != unspecified ? layout.dimension[5] : 1
            , layout.dimension[6] != unspecified ? layout.dimension[6] : 1
            , layout.dimension[7] != unspecified ? layout.dimension[7] : 1 );
    const unsigned fad_dim = computeRank(layout);
    const size_t fad_size = layout.dimension[fad_dim];
    l.dimension[fad_dim] = 1;
    l.dimension[7] = fad_size;

    return l;
  }

  //LayoutStride
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static typename std::enable_if< (std::is_same<Layout , Kokkos::LayoutStride>::value) , Layout>::type createLayout( const Layout& layout )
  {
    Layout      l( layout.dimension[0] != unspecified ? layout.dimension[0] : 1
                 , layout.stride[0]
                 , layout.dimension[1] != unspecified ? layout.dimension[1] : 1
                 , layout.stride[1]
                 , layout.dimension[2] != unspecified ? layout.dimension[2] : 1
                 , layout.stride[2]
                 , layout.dimension[3] != unspecified ? layout.dimension[3] : 1
                 , layout.stride[3]
                 , layout.dimension[4] != unspecified ? layout.dimension[4] : 1
                 , layout.stride[4]
                 , layout.dimension[5] != unspecified ? layout.dimension[5] : 1
                 , layout.stride[5]
                 , layout.dimension[6] != unspecified ? layout.dimension[6] : 1
                 , layout.stride[6]
                 , layout.dimension[7] != unspecified ? layout.dimension[7] : 1
                 , layout.stride[7]
                 );
    const unsigned fad_dim = computeRank(layout);
    const size_t fad_size = layout.dimension[fad_dim];
    l.dimension[fad_dim] = 1;
    l.dimension[7] = fad_size;

    return l;
  }

  // Create a view from the given dimension arguments.
  // This is only necessary because the shmem constructor doesn't take a layout.
  template <typename ViewType, typename ViewArg>
  static ViewType createView( const ViewArg& arg
                            , const size_t N0
                            , const size_t N1
                            , const size_t N2
                            , const size_t N3
                            , const size_t N4
                            , const size_t N5
                            , const size_t N6
                            , const size_t N7 )
  {
    typename ViewType::array_layout l( N0, N1, N2, N3, N4, N5, N6, N7 );
    typename ViewType::array_layout l_fad = createLayout(l);
    return ViewType( arg
                   , l_fad.dimension[0]
                   , l_fad.dimension[1]
                   , l_fad.dimension[2]
                   , l_fad.dimension[3]
                   , l_fad.dimension[4]
                   , l_fad.dimension[5]
                   , l_fad.dimension[6]
                   , l_fad.dimension[7] );
  }

};

}
}
}

#endif //defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "Kokkos_DynRankView.hpp"

namespace Kokkos {

// Overload of dimension_scalar() for all dynamic-rank views
template <typename T, typename ... P>
KOKKOS_INLINE_FUNCTION
constexpr unsigned
dimension_scalar(const Experimental::DynRankView<T,P...>& view) {
  return dimension_scalar(view.ConstDownCast());
}

}

#endif // defined(HAVE_SACADO_KOKKOSCONTAINERS)

#endif /* #ifndef KOKKOS_DYN_RANK_VIEW_SACADO_FAD_HPP */
