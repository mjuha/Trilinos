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

#ifndef KOKKOS_VIEW_FACTORY_HPP
#define KOKKOS_VIEW_FACTORY_HPP

// This only works with the experimental view enabled
#include "Sacado_ConfigDefs.h"
#if defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

#include "Sacado_Traits.hpp"
#include "KokkosExp_View_Fad.hpp"
#include "Kokkos_DynRankView_Fad.hpp"

namespace Kokkos {

namespace Impl {

// Class to determine the value_type for a view as a function of one or more
// input views
template <class ... ViewPack>
struct ViewFactoryType {};

template <class View>
struct ViewFactoryType<View> {
  typedef typename View::value_type type;
};

template <class View, class ... ViewPack>
struct ViewFactoryType<View,ViewPack...> {
  typedef typename Sacado::Promote<
    typename View::value_type,
    typename ViewFactoryType<ViewPack...>::type
    >::type type;
};

}

// Function to compute the scalar dimension (e.g., Fad dimesion) from one or
// more views.  It relies on the overload for a single view provided by Sacado
template <class View, class ... ViewPack>
unsigned dimension_scalar(const View& v, const ViewPack&... views) {
  const unsigned dim0 = dimension_scalar(v);
  const unsigned dim1 = dimension_scalar(views...);
  return dim0 >= dim1 ? dim0 : dim1 ;
}

template<typename I, class ... DimPack>
struct dim_count {
  enum {count = dim_count<DimPack...>::count+1};
};

template<typename I>
struct dim_count<I> {
  enum {count=1};
};
// Traits class used to create a view for a given rank and dimension as a
// function of one or more views.  The value_type for the view is determined
// by value_type, and the view is created through the create_view() function.
// The calling code must determine the rank and dimensions of the view to create
// however internal Sacado dimension will be determined automatically.
template <class ... ViewPack>
struct ViewFactory {

  typedef typename Impl::ViewFactoryType<ViewPack...>::type value_type;

  template <class ResultView, class CtorProp, class ... Dims>
  static ResultView
  create_view(const ViewPack& ... views,
              const CtorProp& prop,
              const Dims ... dims) {

    constexpr bool is_dyn_rank = is_dyn_rank_view<ResultView>::value;
    constexpr bool is_fad = Sacado::IsScalarType<typename ResultView::non_const_value_type>::value?false:true; 
 
    typename ResultView::array_layout layout_extern(dims...);
    int rank_extern = dim_count<Dims...>::count;

    if(is_fad)
        layout_extern.dimension[rank_extern] = dimension_scalar(views...); 
 
   if(is_dyn_rank) {
    
      typename ResultView::array_layout layout = Experimental::Impl::reconstructLayout(layout_extern, rank_extern);
    
      const unsigned rank = computeRank(layout);
      if(is_fad) 
      layout.dimension[rank] = dimension_scalar(views...);
      return ResultView(prop, layout);
    } else {
      return ResultView(prop, layout_extern);   
    }
  }

  // Compute the view rank from the dimension arguments
  // This allows the code to work for both static and dynamic-rank views
  template <typename Layout>
  static size_t
  computeRank(const Layout& layout) {
    return computeRank( layout.dimension[0], layout.dimension[1],
                        layout.dimension[2], layout.dimension[3],
                        layout.dimension[4], layout.dimension[5],
                        layout.dimension[6], layout.dimension[7] );
  }

  // Compute the view rank from the dimension arguments
  // This allows the code to work for both static and dynamic-rank views
  static size_t
  computeRank(
    const size_t N0, const size_t N1, const size_t N2, const size_t N3,
    const size_t N4, const size_t N5, const size_t N6, const size_t N7 ) {
    return  ( (N7 == ~size_t(0)) ?
            ( (N6 == ~size_t(0)) ?
            ( (N5 == ~size_t(0)) ?
            ( (N4 == ~size_t(0)) ?
            ( (N3 == ~size_t(0)) ?
            ( (N2 == ~size_t(0)) ?
            ( (N1 == ~size_t(0)) ?
            ( (N0 == ~size_t(0)) ? 0 : 1 ) : 2 ) : 3 ) : 4 ) : 5 ) : 6 ) : 7 ) : 8 );
  }

};

//! Wrapper to simplify use of Sacado ViewFactory
template <typename ResultViewType, typename InputViewType, typename CtorProp,
          typename ... Dims>
typename std::enable_if<
  is_view<InputViewType>::value || is_dyn_rank_view<InputViewType>::value,
  ResultViewType>::type
createDynRankViewWithType(const InputViewType& a,
                          const CtorProp& prop,
                          const Dims... dims)
{
  using view_factory = Kokkos::ViewFactory<InputViewType>;
  return view_factory::template create_view<ResultViewType>(a,prop,dims...);
}

//! Wrapper to simplify use of Sacado ViewFactory
template <typename InputViewType, typename CtorProp, typename ... Dims >
typename std::enable_if<
  is_view<InputViewType>::value || is_dyn_rank_view<InputViewType>::value,
  Kokkos::DynRankView<typename InputViewType::non_const_value_type,
                      typename InputViewType::array_layout,
                      typename InputViewType::device_type> >::type
createDynRankView(const InputViewType& a,
                  const CtorProp& prop,
                  const Dims... dims)
{
  using ResultViewType =
    Kokkos::DynRankView<typename InputViewType::non_const_value_type,
                        typename InputViewType::array_layout,
                        typename InputViewType::device_type>;
  return createDynRankViewWithType<ResultViewType>(a, prop, dims...);
}

//! Wrapper to simplify use of Sacado ViewFactory
template <typename ResultViewType, typename InputViewType, typename CtorProp,
          typename ... Dims>
typename std::enable_if<
  is_view<InputViewType>::value || is_dyn_rank_view<InputViewType>::value,
  ResultViewType>::type
createViewWithType(const InputViewType& a,
                   const CtorProp& prop,
                   const Dims... dims)
{
  using view_factory = Kokkos::ViewFactory<InputViewType>;
  return view_factory::template create_view<ResultViewType>(a,prop,dims...);
}

}

#endif

#endif /* #ifndef KOKKOS_VIEW_FACTORY_HPP */
