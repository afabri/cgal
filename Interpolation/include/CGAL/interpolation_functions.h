// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Julia Floetotto

#ifndef CGAL_INTERPOLATION_FUNCTIONS_H
#define CGAL_INTERPOLATION_FUNCTIONS_H

#include <CGAL/license/Interpolation.h>

#include <CGAL/double.h>
#include <CGAL/use.h>

#include <iterator>
#include <utility>
#include <vector>

namespace CGAL {

//Functor class for accessing the function values/gradients
template< class Map >
struct Data_access
{
  typedef typename Map::mapped_type Data_type;
  typedef typename Map::key_type  Key_type;
  typedef CGAL::cpp11::tuple<Key_type,Data_type,Data_type,bool>  result_type;
  
  Data_access< Map >(const Map& m): map(m){}

  std::pair< Data_type, bool>
  operator()(const Key_type& p) const {
    typename Map::const_iterator mit = map.find(p);
    if(mit!= map.end())
      return std::make_pair(mit->second, true);
    return std::make_pair(Data_type(), false);
  }


  template <typename T1, typename T2>
  CGAL::cpp11::tuple<Key_type,Data_type,Data_type,bool>
  operator()(const std::pair<T1,T2>& t)const
  {
    std::pair< Data_type, bool> oldres = operator()(t.first);
    return CGAL::cpp11::tuple<Key_type,Data_type,Data_type,bool>(t.first, oldres.first, t.second, oldres.second);
  }

  
  const Map& map;
};

//the interpolation functions:
template < class ForwardIterator, class Functor>
typename CGAL::cpp11::tuple_element<1,typename Functor::result_type>::type
linear_interpolation(ForwardIterator first, ForwardIterator beyond,
                     const typename
                       std::iterator_traits<ForwardIterator>::value_type::
                       second_type& norm,
                     Functor function_value)
{
  CGAL_precondition(norm>0);

  typedef typename CGAL::cpp11::tuple_element<1,typename Functor::result_type>::type Value_type;
  Value_type result(0);
  typename Functor::result_type frt;
  for(; first !=beyond; ++first){
    frt = function_value(*first); // tuple<Point,value,bary,bool>
    CGAL_assertion(get<3>(frt));
    result += (get<2>(frt)/norm) * get<1>(frt);
  }

  return result;
}


template < class ForwardIterator, class Functor, class GradFunctor, class Traits>
std::pair< typename CGAL::cpp11::tuple_element<1,typename Functor::result_type>::type, bool>
quadratic_interpolation(ForwardIterator first, ForwardIterator beyond,
                        const typename
                          std::iterator_traits<ForwardIterator>::
                          value_type::second_type& norm,
                        const typename Traits::Point_d& p,
                        Functor function_value,
                        GradFunctor function_gradient,
                        const Traits& traits)
{
  CGAL_precondition(norm >0);
  typedef typename CGAL::cpp11::tuple_element<1,typename Functor::result_type>::type Value_type;
  Value_type result(0);
  typename Functor::result_type frt;
  typename GradFunctor::result_type grad;
  for(; first !=beyond; ++first){
    frt = function_value(*first); // tuple<Point,value,bary,bool>
    grad = function_gradient(*first);
    //test if value and gradient are correctly retrieved:
    CGAL_assertion(get<3>(frt));
    if(!grad.second)
      return std::make_pair(Value_type(0), false);
    result += (get<2>(frt)/norm)
      *( get<1>(frt) + grad.first*
                 traits.construct_scaled_vector_d_object()
         (traits.construct_vector_d_object()(get<0>(frt), p),0.5));
  }
  return std::make_pair(result, true);
}


template < class ForwardIterator, class Functor, class GradFunctor, class Traits>
std::pair< typename CGAL::cpp11::tuple_element<1,typename Functor::result_type>::type, bool>
sibson_c1_interpolation(ForwardIterator first, ForwardIterator beyond,
                        const typename
                          std::iterator_traits<ForwardIterator>::
                          value_type::second_type& norm,
                        const typename Traits::Point_d& p,
                        Functor function_value,
                        GradFunctor function_gradient,
                        const Traits& traits)
{
  CGAL_precondition(norm >0);
  typedef typename CGAL::cpp11::tuple_element<1,typename Functor::result_type>::type Value_type;
  typedef typename Traits::FT                       Coord_type;

  Coord_type term1(0), term2(term1), term3(term1), term4(term1);
  Value_type linear_int(0),gradient_int(0);
  typename Functor::result_type frt;
  typename GradFunctor::result_type grad;

  for(; first !=beyond; ++first){
    frt = function_value(*first);
    grad = function_gradient(*first);
    CGAL_assertion(get<3>(frt));
    if(!grad.second)
      //the values are not correct:
      return std::make_pair(Value_type(0), false);

    Coord_type coeff = get<2>(frt)/norm;
    Coord_type squared_dist = traits.
      compute_squared_distance_d_object()(get<0>(frt), p);
    Coord_type dist = CGAL_NTS sqrt(squared_dist);

    if(squared_dist ==0){
      ForwardIterator it = first;
      CGAL_USE(it);
      CGAL_assertion(++it==beyond);
      return std::make_pair(get<1>(frt), true);
    }
    //three different terms to mix linear and gradient
    //interpolation
    term1 +=  coeff / dist;
    term2 +=  coeff * squared_dist;
    term3 +=  coeff * dist;

    linear_int += coeff * get<1>(frt);

    gradient_int += (coeff/dist)
      * (get<1>(frt) + grad.first *
         traits.construct_vector_d_object()(get<0>(frt), p));
  }

  term4 = term3/ term1;
  gradient_int = gradient_int / term1;

  return std::make_pair((term4* linear_int + term2 * gradient_int)/
                        (term4 + term2), true);
}

//this method works with rational number types:
//modification of Sibson's interpolant without sqrt
//following a proposition by Gunther Rote:
//
// the general scheme:
//  Coord_type inv_weight = f(dist); //i.e. dist^2
//   	term1 +=  coeff/inv_weight;
//	term2 +=  coeff * squared_dist;
//	term3 +=  coeff*(squared_dist/inv_weight);
// 	gradient_int +=  (coeff/inv_weight)*
// 	  (vh->get_value()+  vh->get_gradient()
// 	   *(p - vh->point()));

template < class ForwardIterator, class Functor, class GradFunctor, class Traits>
std::pair< typename CGAL::cpp11::tuple_element<1,typename Functor::result_type>::type, bool>
sibson_c1_interpolation_square(ForwardIterator first, ForwardIterator beyond,
                               const typename
                                 std::iterator_traits<ForwardIterator>::
                                 value_type::second_type& norm,
                               const typename Traits::Point_d& p,
                               Functor function_value,
                               GradFunctor function_gradient,
                               const Traits& traits)
{
  CGAL_precondition(norm >0);
  typedef typename CGAL::cpp11::tuple_element<1,typename Functor::result_type>::type Value_type;
  typedef typename Traits::FT                       Coord_type;

  Coord_type term1(0), term2(term1), term3(term1), term4(term1);
  Value_type linear_int(0),gradient_int(0);
  typename Functor::result_type frt;
  typename GradFunctor::result_type grad;

  for(; first!=beyond; ++first){
    frt = function_value(*first); // tuple<Point,value,bary,bool>
    grad = function_gradient(*first);
    CGAL_assertion(get<3>(frt));
    if(!grad.second)
      //the gradient is not known
      return std::make_pair(Value_type(0), false);

    Coord_type coeff = get<2>(frt)/norm;
    Coord_type squared_dist = traits.
      compute_squared_distance_d_object()(get<0>(frt), p);

    if(squared_dist ==0){
      ForwardIterator it = first;
      CGAL_USE(it);
      CGAL_assertion(++it==beyond);
      return std::make_pair(get<1>(frt),true);
    }
    //three different terms to mix linear and gradient
    //interpolation
    term1 +=  coeff / squared_dist;
    term2 +=  coeff * squared_dist;
    term3 +=  coeff;

    linear_int += coeff * get<1>(frt);

    gradient_int += (coeff/squared_dist) * (get<1>(frt) + grad.first *
                                                           traits.construct_vector_d_object()(get<0>(frt), p));
  }

  term4 = term3/ term1;
  gradient_int = gradient_int / term1;

  return std::make_pair((term4 * linear_int + term2 * gradient_int)/
                        (term4 + term2), true);
}


template < class RandomAccessIterator, class Functor, class
           GradFunctor, class Traits>
std::pair<typename CGAL::cpp11::tuple_element<1,typename Functor::result_type>::type, bool>
farin_c1_interpolation(RandomAccessIterator first,
                       RandomAccessIterator beyond,
                       const typename
                         std::iterator_traits<RandomAccessIterator>::
                         value_type::second_type& norm,
                       const typename Traits::Point_d& /*p*/,
                       Functor function_value,
                       GradFunctor function_gradient,
                       const Traits& traits)
{
  CGAL_precondition(norm >0);
  //the function value is available for all points
  //if a gradient value is not availble: function returns false
  typedef typename CGAL::cpp11::tuple_element<1,typename Functor::result_type>::type Value_type;
  typedef typename Traits::FT                        Coord_type;

  typename Functor::result_type frt;
  typename GradFunctor::result_type grad;
  int n= static_cast<int>(beyond - first);
  if( n==1){
    frt = function_value(*first); // tuple<Point,value,bary,bool>
    CGAL_assertion(get<3>(frt));
    return std::make_pair(get<1>(frt), true);
  }

  //there must be one or at least three NN-neighbors:
  CGAL_assertion(n > 2);

  RandomAccessIterator it2, it;
  typename Functor::result_type frt2;
  Value_type result(0);
  const Coord_type fac3(3);

  std::vector< std::vector<Value_type> >
      ordinates(n,std::vector<Value_type>(n, Value_type(0)));
  // *it was std::pair of Point and barycenter
  for(int i =0; i<n; ++i){
    it = first+i;
    frt = function_value(*it); // tuple<Point,value,bary,bool>
    Coord_type coord_i_square = CGAL_NTS square(get<2>(frt));

    //for later: the function value of it->first:

    CGAL_assertion(get<3>(frt));
    ordinates[i][i] = get<1>(frt);

    //control point = data point
    result += coord_i_square * get<2>(frt)* ordinates[i][i];

    //compute tangent plane control point (one 2, one 1 entry)
    Value_type res_i(0);
    for(int j =0; j<n; ++j){
      if(i!=j){
        it2 = first+j;
        frt2 =  function_value(*it2);

        grad = function_gradient(*it);
        if(!grad.second)
          //the gradient is not known
          return std::make_pair(Value_type(0), false);

        //ordinates[i][j] = (p_j - p_i) * g_i
        ordinates[i][j] = grad.first *
          traits.construct_vector_d_object()(get<0>(frt),get<0>(frt2));

        // a point in the tangent plane:
        // 3( f(p_i) + (1/3)(p_j - p_i) * g_i)
        // => 3*f(p_i) + (p_j - p_i) * g_i
        res_i += (fac3 * ordinates[i][i] + ordinates[i][j])* get<2>(frt2);
      }
    }
    //res_i already multiplied by three
    result += coord_i_square *res_i;
  }

  //the third type of control points: three 1 entries i,j,k
  for(int i=0; i< n; ++i)
    for(int j=i+1; j< n; ++j)
      for(int k=j+1; k<n; ++k){
        // add 6* (u_i*u_j*u_k) * b_ijk
        //  b_ijk = 1.5 * a - 0.5*c
        //where
        //c : average of the three data control points
        //a : 1.5*a = 1/12 * (ord[i][j] + ord[i][k] + ord[j][i] +
        //            ord[j][k] + ord[k][i]+ ord[k][j])
        // =>  6 * b_ijk = 3*(f_i + f_j + f_k) + 0.5*a
        result += (Coord_type(2.0)*( ordinates[i][i]+ ordinates[j][j]+
                                     ordinates[k][k])
                   + Coord_type(0.5)*(ordinates[i][j] + ordinates[i][k]
                                      + ordinates[j][i] +
                                      ordinates[j][k] + ordinates[k][i]+
                                      ordinates[k][j]))
                  *(first+i)->second *(first+j)->second *(first+k)->second ;
      }

  return std::make_pair(result/(CGAL_NTS square(norm)*norm), true);

  
}

} //namespace CGAL

#endif // CGAL_INTERPOLATION_FUNCTIONS_H
