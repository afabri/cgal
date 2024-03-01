// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
//
//
// Author(s) :
//

#ifndef AABB_GEOMTRAITS_H_
#define AABB_GEOMTRAITS_H_


#include <CGAL/Bbox_3.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/AABB_intersections.h>
#include <CGAL/internal/AABB_tree/Has_nested_type_Shared_data.h>
#include <CGAL/internal/AABB_tree/Primitive_helper.h>
#include <boost/optional.hpp>
#include <boost/bind.hpp>
#include <CGAL/internal/AABB_tree/AABB_do_contain_traits.h>
/// \file AABB_traits_d.h

namespace CGAL {

namespace internal{  namespace AABB_tree {

template <class T>
struct Remove_optional  { typedef T type; };

template <class T>
struct Remove_optional< ::boost::optional<T> >  { typedef T type; };

//helper controlling whether extra data should be stored in the AABB_tree traits class
template <class Primitive, bool has_shared_data=Has_nested_type_Shared_data<Primitive>::value>
struct AABB_traits_base;

template <class Primitive>
struct AABB_traits_base<Primitive,false>{};

template <class Primitive>
struct AABB_traits_base<Primitive,true>{
  typename  Primitive::Shared_data m_primitive_data;

  #ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template <class PrimitiveType, typename ... T>
  void set_shared_data(T ... t){
    m_primitive_data=PrimitiveType::construct_shared_data(t...);
  }
  #else
  template <class PrimitiveType>
  void set_shared_data(){
    m_primitive_data=PrimitiveType::construct_shared_data();
  }

  template <class PrimitiveType, class T1>
  void set_shared_data(T1 t1){
    m_primitive_data=PrimitiveType::construct_shared_data(t1);
  }

  template <class PrimitiveType, class T1,class T2,class T3>
  void set_shared_data(T1 t1,T2 t2,T3 t3){
    m_primitive_data=PrimitiveType::construct_shared_data(t1,t2,t3);
  }

  template <class PrimitiveType, class T1,class T2,class T3,class T4>
  void set_shared_data(T1 t1,T2 t2,T3 t3,T4 t4){
    m_primitive_data=PrimitiveType::construct_shared_data(t1,t2,t3,t4);
  }

  template <class PrimitiveType, class T1,class T2,class T3,class T4,class T5>
  void set_shared_data(T1 t1,T2 t2,T3 t3,T4 t4,T5 t5){
    m_primitive_data=PrimitiveType::construct_shared_data(t1,t2,t3,t4,t5);
  }
  #endif
  const typename Primitive::Shared_data& shared_data() const {return m_primitive_data;}
};

} }


template<typename GeomTraits, typename AABBPrimitive, int Dimension>
struct AABB_traits_d:public internal::AABB_tree::AABB_traits_base<AABBPrimitive>
{ 
	typedef AABB_traits_d <GeomTraits, AABBPrimitive,Dimension>  AT;
};

//Template specification for 2D tree
template<typename GeomTraits, typename AABBPrimitive>
struct AABB_traits_d <GeomTraits , AABBPrimitive ,2>:public internal::AABB_tree::AABB_traits_base<AABBPrimitive>
{
  typedef typename CGAL::Object Object;
  typedef typename GeomTraits::FT FT;
  typedef AABBPrimitive Primitive;

  typedef typename std::pair<Object,typename Primitive::Id> Object_and_primitive_id;

  typedef AABB_traits_d <GeomTraits, AABBPrimitive,2>  AT;

  typedef typename std::pair<typename GeomTraits::Point_2, typename Primitive::Id> Point_and_primitive_id;


  template<typename Query>
   struct Intersection_and_primitive_id {
     typedef typename cpp11::result_of<
       typename GeomTraits::Intersect_2(Query, typename Primitive::Datum)
     >::type Intersection_type;

     typedef std::pair<
       typename internal::AABB_tree::Remove_optional<Intersection_type>::type,
       typename Primitive::Id > Type;
   };

  //Dimension dependent definitions

   typedef enum { CGAL_AXIS_X = 0,
                 CGAL_AXIS_Y = 1} Axis;

  typedef typename GeomTraits::Point_2 Point;

  typedef typename GeomTraits::Iso_rectangle_2 Iso_box_d;


  typedef typename CGAL::Bbox_2 Bounding_box;

  /// @}

  typedef typename GeomTraits::Circle_2 Sphere_d;

  typedef typename GeomTraits::Cartesian_const_iterator_2 Cartesian_const_iterator_d;
  typedef typename GeomTraits::Construct_cartesian_const_iterator_2 Construct_cartesian_const_iterator_d;
  typedef typename GeomTraits::Construct_vertex_2 Construct_vertex_d;
  typedef typename GeomTraits::Construct_min_vertex_2 Construct_min_vertex_d;
  typedef typename GeomTraits::Construct_max_vertex_2 Construct_max_vertex_d;
  typedef typename GeomTraits::Construct_iso_rectangle_2 Construct_iso_box_d;

  
  typedef Do_contain_test_traits<AT,GeomTraits,typename Primitive::Datum> Do_primitive_contain_test;
  typedef Do_contain_test_traits<AT,GeomTraits, Iso_box_d> Do_bbox_contain_test;


  AABB_traits_d(){ };

  typedef typename GeomTraits::Compute_squared_distance_2 Compute_squared_distance;
   Compute_squared_distance compute_squared_distance_object() const { return GeomTraits().compute_squared_distance_2_object(); }


   class Sort_primitives
    {
      const AABB_traits_d<GeomTraits,AABBPrimitive,2>& m_traits;
    public:
      Sort_primitives(const AABB_traits_d<GeomTraits,AABBPrimitive,2>& traits)
        : m_traits(traits) {}

      template<typename PrimitiveIterator>
      void operator()(PrimitiveIterator first,
                      PrimitiveIterator beyond,
                      const typename AT::Bounding_box& bbox) const
        {

          PrimitiveIterator middle = first + (beyond - first)/2;
          switch(longest_axis(bbox))
          {
          case AT::CGAL_AXIS_X: // sort along x
            std::nth_element(first, middle, beyond, boost::bind(less_x,_1,_2,m_traits));
            break;
          case AT::CGAL_AXIS_Y: // sort along y
            std::nth_element(first, middle, beyond, boost::bind(less_y,_1,_2,m_traits));
            break;
          default:
            CGAL_error();
          }
        }
    };

   Sort_primitives sort_primitives_object() const {return Sort_primitives(*this);}

  //As the previous proposal

   class Do_intersect {
       const AABB_traits_d<GeomTraits,AABBPrimitive,2>& m_traits;

     public:
   	typedef typename GeomTraits::Iso_rectangle_2 IsoRect;
       Do_intersect(const AABB_traits_d<GeomTraits,AABBPrimitive,2>& traits)
         :m_traits(traits) {}

       template<typename Query>
       bool operator()(const Query& q, const Bounding_box& bbox) const
       {

   		return CGAL::do_intersect(q, bbox);
       }

       template<typename Query>
       bool operator()(const Query& q, const Primitive& pr) const
       {
    	 return true;//GeomTraits().do_intersect_2_object()(q, internal::Primitive_helper<AT>::get_datum(pr,m_traits));
       }
     };

     Do_intersect do_intersect_object() const {return Do_intersect(*this);}


     class Do_contain {

           const AABB_traits_d<GeomTraits,AABBPrimitive,2>& m_traits;

           public:
             Do_contain(const AABB_traits_d<GeomTraits,AABBPrimitive,2>& traits)
               :m_traits(traits) {}
             
             template<typename Query>
             bool operator()(const Query& q, const Bounding_box& bbox) const
             {
				return Do_bbox_contain_test()(q,(Iso_box_d)bbox);
             }

             template<typename Query>
             bool operator()(const Query& q, const Primitive& pr) const
             {
				return Do_primitive_contain_test()(q,internal::Primitive_helper<AT>::get_datum(pr,m_traits));
             }
           };

    Do_contain do_contain_object() const {return Do_contain(*this);}
     
     class Intersection {
       const AABB_traits_d<GeomTraits,AABBPrimitive,2>& m_traits;
     public:
       Intersection(const AABB_traits_d<GeomTraits,AABBPrimitive,2>& traits)
         :m_traits(traits) {}
         #if CGAL_INTERSECTION_VERSION < 2
     template<typename Query>
     boost::optional<typename AT::Object_and_primitive_id>
     operator()(const Query& query, const typename AT::Primitive& primitive) const
     {
       typedef boost::optional<Object_and_primitive_id> Intersection;

       CGAL::Object object = GeomTraits().intersect_2_object()(internal::Primitive_helper<AT>::get_datum(primitive,m_traits),query);
       if ( object.empty() )
         return Intersection();
       else
         return Intersection(Object_and_primitive_id(object,primitive.id()));
     }
         #else
         template<typename Query>
         boost::optional< typename Intersection_and_primitive_id<Query>::Type >
         operator()(const Query& query, const typename AT::Primitive& primitive) const {
           typename cpp11::result_of<typename GeomTraits::Intersect_2(Query, typename Primitive::Datum) >::type
             inter_res = GeomTraits().intersect_2_object()(internal::Primitive_helper<AT>::get_datum(primitive,m_traits),query);
           if (!inter_res)
               return boost::optional<typename Intersection_and_primitive_id<Query>::Type>();
           return boost::make_optional( std::make_pair(*inter_res, primitive.id()) );
         }
         #endif
     };

   Intersection intersection_object() const {return Intersection(*this);}



   // This should go down to the GeomTraits, i.e. the kernel
     class Closest_point {
         typedef typename AT::Point Point;
         typedef typename AT::Primitive Primitive;
       const AABB_traits_d<GeomTraits,AABBPrimitive,2>& m_traits;
     public:
       Closest_point(const AABB_traits_d<GeomTraits,AABBPrimitive,2>& traits)
         : m_traits(traits) {}


       Point operator()(const Point& p, const Primitive& pr, const Point& bound) const
       {
           return CGAL::nearest_point_2(p, internal::Primitive_helper<AT>::get_datum(pr,m_traits), bound);
       }
     };

     // This should go down to the GeomTraits, i.e. the kernel
     // and the internal implementation should change its name from
     // do_intersect to something like does_contain (this is what we compute,
     // this is not the same do_intersect as the spherical kernel)
     class Compare_distance {
         typedef typename AT::Point Point;
         typedef typename AT::FT FT;
         typedef typename AT::Primitive Primitive;
     public:
         template <class Solid>
         CGAL::Comparison_result operator()(const Point& p, const Solid& pr, const Point& bound) const
         {
             return GeomTraits().do_intersect_2_object()
             (GeomTraits().construct_circle_2_object()
             (p, GeomTraits().compute_squared_distance_2_object()(p, bound)), pr)?
             CGAL::SMALLER : CGAL::LARGER;
         }

         template <class Solid>
         CGAL::Comparison_result operator()(const Point& p, const Solid& pr, const FT& sq_distance) const
         {
           return GeomTraits().do_intersect_2_object()
             (GeomTraits().construct_circle_2_object()(p, sq_distance), pr) ?
             CGAL::SMALLER :
             CGAL::LARGER;
         }
     };

     Closest_point closest_point_object() const {return Closest_point(*this);}
     Compare_distance compare_distance_object() const {return Compare_distance();}



   static Axis longest_axis(const Bounding_box& bbox);

     /// Comparison functions
     static bool less_x(const Primitive& pr1, const Primitive& pr2,const AABB_traits_d<GeomTraits,AABBPrimitive,2>& traits)
     { return internal::Primitive_helper<AT>::get_reference_point(pr1,traits).x() < internal::Primitive_helper<AT>::get_reference_point(pr2,traits).x(); }
     static bool less_y(const Primitive& pr1, const Primitive& pr2,const AABB_traits_d<GeomTraits,AABBPrimitive,2>& traits)
     { return internal::Primitive_helper<AT>::get_reference_point(pr1,traits).y() < internal::Primitive_helper<AT>::get_reference_point(pr2,traits).y(); }
};

template<typename GT, typename P>
typename AABB_traits_d<GT,P,2>::Axis
AABB_traits_d<GT,P,2>::longest_axis(const Bounding_box& bbox)
{
  const double dx = bbox.xmax() - bbox.xmin();
  const double dy = bbox.ymax() - bbox.ymin();
  if(dx>=dy)
  {
    return CGAL_AXIS_X;

  }
  else // dy>dx
  {
    return CGAL_AXIS_Y;
  }
}


//Template specification for 3D AABB_traits_d
template<typename GeomTraits, typename AABBPrimitive>
struct AABB_traits_d < GeomTraits , AABBPrimitive ,3>:public internal::AABB_tree::AABB_traits_base<AABBPrimitive>
{
  typedef typename CGAL::Object Object;
  typedef typename GeomTraits::FT FT;
  typedef AABBPrimitive Primitive;
  typedef typename CGAL::Bbox_3 Bounding_box;

  typedef typename std::pair<Object,typename Primitive::Id> Object_and_primitive_id;

  typedef AABB_traits_d <GeomTraits, AABBPrimitive,3>  AT;

  typedef typename std::pair<typename GeomTraits::Point_3, typename Primitive::Id> Point_and_primitive_id;


  template<typename Query>
   struct Intersection_and_primitive_id {
     typedef typename cpp11::result_of<
       typename GeomTraits::Intersect_3(Query, typename Primitive::Datum)
     >::type Intersection_type;

     typedef std::pair<
       typename internal::AABB_tree::Remove_optional<Intersection_type>::type,
       typename Primitive::Id > Type;
   };

  //Dimension dependent definitions
   typedef enum { CGAL_AXIS_X = 0,
                 CGAL_AXIS_Y = 1,
                 CGAL_AXIS_Z = 2} Axis;

  typedef typename GeomTraits::Point_3 Point;

  typedef typename GeomTraits::Iso_cuboid_3 Iso_box_d;


  typedef typename GeomTraits::Sphere_3 Sphere_d;

  typedef typename GeomTraits::Cartesian_const_iterator_3 Cartesian_const_iterator_d;
  typedef typename GeomTraits::Construct_cartesian_const_iterator_3 Construct_cartesian_const_iterator_d;
  typedef typename GeomTraits::Construct_vertex_3 Construct_vertex_d;
  typedef typename GeomTraits::Construct_min_vertex_3 Construct_min_vertex_d;
  typedef typename GeomTraits::Construct_max_vertex_3 Construct_max_vertex_d;
  typedef typename GeomTraits::Construct_iso_cuboid_3 Construct_iso_box_d;

  
  typedef Do_contain_test_traits<AT,GeomTraits,typename Primitive::Datum> Do_primitive_contain_test;
  typedef Do_contain_test_traits<AT,GeomTraits, Iso_box_d> Do_bbox_contain_test;

  
  
  AABB_traits_d(){ };

  typedef typename GeomTraits::Compute_squared_distance_3 compute_squared_distance;
  compute_squared_distance compute_squared_distance_object() const { return GeomTraits().compute_squared_distance_3_object(); }


     class Sort_primitives
      {
        const AABB_traits_d<GeomTraits,AABBPrimitive,3>& m_traits;
      public:
        Sort_primitives(const AABB_traits_d<GeomTraits,AABBPrimitive,3>& traits)
          : m_traits(traits) {}

        template<typename PrimitiveIterator>
        void operator()(PrimitiveIterator first,
                        PrimitiveIterator beyond,
                        const typename AT::Bounding_box& bbox) const
          {

            PrimitiveIterator middle = first + (beyond - first)/2;
            switch(longest_axis(bbox))
            {
            case AT::CGAL_AXIS_X: // sort along x
              std::nth_element(first, middle, beyond, boost::bind(less_x,_1,_2,m_traits));
              break;
            case AT::CGAL_AXIS_Y: // sort along y
              std::nth_element(first, middle, beyond, boost::bind(less_y,_1,_2,m_traits));
              break;
            case AT::CGAL_AXIS_Z: // sort along y
			  std::nth_element(first, middle, beyond, boost::bind(less_z,_1,_2,m_traits));
			  break;
            default:
              CGAL_error();
            }
          }
      };

     Sort_primitives sort_primitives_object() const {return Sort_primitives(*this);}

    //As the previous proposal

     class Do_intersect {

         const AABB_traits_d<GeomTraits,AABBPrimitive,3>& m_traits;

       public:
         Do_intersect(const AABB_traits_d<GeomTraits,AABBPrimitive,3>& traits)
           :m_traits(traits) {}

         template<typename Query>
         bool operator()(const Query& q, const Bounding_box& bbox) const
         {

     		return CGAL::do_intersect(q, bbox);
         }

         template<typename Query>
         bool operator()(const Query& q, const Primitive& pr) const
         {
      	 return GeomTraits().do_intersect_3_object()(q, internal::Primitive_helper<AT>::get_datum(pr,m_traits));
         }
       };

       Do_intersect do_intersect_object() const {return Do_intersect(*this);}


       class Do_contain {

              const AABB_traits_d<GeomTraits,AABBPrimitive,3>& m_traits;

              public:
                Do_contain(const AABB_traits_d<GeomTraits,AABBPrimitive,3>& traits)
                  :m_traits(traits) {}
                
                template<typename Query>
                bool operator()(const Query& q, const Bounding_box& bbox) const
                {
					return Do_bbox_contain_test()(q,(Iso_box_d)bbox);
                }

                template<typename Query>
                bool operator()(const Query& q, const Primitive& pr) const
                {
					return Do_primitive_contain_test()(q,internal::Primitive_helper<AT>::get_datum(pr,m_traits));
                }
              };

       Do_contain do_contain_object() const {return Do_contain(*this);}



       class Intersection {
         const AABB_traits_d<GeomTraits,AABBPrimitive,3>& m_traits;
       public:
         Intersection(const AABB_traits_d<GeomTraits,AABBPrimitive,3>& traits)
           :m_traits(traits) {}
           #if CGAL_INTERSECTION_VERSION < 2
       template<typename Query>
       boost::optional<typename AT::Object_and_primitive_id>
       operator()(const Query& query, const typename AT::Primitive& primitive) const
       {
         typedef boost::optional<Object_and_primitive_id> Intersection;

         CGAL::Object object = GeomTraits().intersect_3_object()(internal::Primitive_helper<AT>::get_datum(primitive,m_traits),query);
         if ( object.empty() )
           return Intersection();
         else
           return Intersection(Object_and_primitive_id(object,primitive.id()));
       }
           #else
           template<typename Query>
           boost::optional< typename Intersection_and_primitive_id<Query>::Type >
           operator()(const Query& query, const typename AT::Primitive& primitive) const {
             typename cpp11::result_of<typename GeomTraits::Intersect_3(Query, typename Primitive::Datum) >::type
               inter_res = GeomTraits().intersect_3_object()(internal::Primitive_helper<AT>::get_datum(primitive,m_traits),query);
             if (!inter_res)
                 return boost::optional<typename Intersection_and_primitive_id<Query>::Type>();
             return boost::make_optional( std::make_pair(*inter_res, primitive.id()) );
           }
           #endif
       };

     Intersection intersection_object() const {return Intersection(*this);}



     // This should go down to the GeomTraits, i.e. the kernel
       class Closest_point {
           typedef typename AT::Point Point;
           typedef typename AT::Primitive Primitive;
         const AABB_traits_d<GeomTraits,AABBPrimitive,3>& m_traits;
       public:
         Closest_point(const AABB_traits_d<GeomTraits,AABBPrimitive,3>& traits)
           : m_traits(traits) {}


         Point operator()(const Point& p, const Primitive& pr, const Point& bound) const
         {
             return CGAL::nearest_point_3(p, internal::Primitive_helper<AT>::get_datum(pr,m_traits), bound);
         }
       };

       // This should go down to the GeomTraits, i.e. the kernel
       // and the internal implementation should change its name from
       // do_intersect to something like does_contain (this is what we compute,
       // this is not the same do_intersect as the spherical kernel)
       class Compare_distance {
           typedef typename AT::Point Point;
           typedef typename AT::FT FT;
           typedef typename AT::Primitive Primitive;
       public:
           template <class Solid>
           CGAL::Comparison_result operator()(const Point& p, const Solid& pr, const Point& bound) const
           {
               return GeomTraits().do_intersect_3_object()
               (GeomTraits().construct_sphere_3_object()
               (p, GeomTraits().compute_squared_distance_3_object()(p, bound)), pr)?
               CGAL::SMALLER : CGAL::LARGER;
           }

           template <class Solid>
           CGAL::Comparison_result operator()(const Point& p, const Solid& pr, const FT& sq_distance) const
           {
             return GeomTraits().do_intersect_3_object()
               (GeomTraits().construct_sphere_3_object()(p, sq_distance),pr) ?
               CGAL::SMALLER :
               CGAL::LARGER;
           }
       };

       Closest_point closest_point_object() const {return Closest_point(*this);}
       Compare_distance compare_distance_object() const {return Compare_distance();}



     static Axis longest_axis(const Bounding_box& bbox);

       static bool less_x(const Primitive& pr1, const Primitive& pr2,const AABB_traits_d<GeomTraits,AABBPrimitive,3>& traits)
       { return internal::Primitive_helper<AT>::get_reference_point(pr1,traits).x() < internal::Primitive_helper<AT>::get_reference_point(pr2,traits).x(); }
       static bool less_y(const Primitive& pr1, const Primitive& pr2,const AABB_traits_d<GeomTraits,AABBPrimitive,3>& traits)
       { return internal::Primitive_helper<AT>::get_reference_point(pr1,traits).y() < internal::Primitive_helper<AT>::get_reference_point(pr2,traits).y(); }
       static bool less_z(const Primitive& pr1, const Primitive& pr2,const AABB_traits_d<GeomTraits,AABBPrimitive,3>& traits)
       { return internal::Primitive_helper<AT>::get_reference_point(pr1,traits).z() < internal::Primitive_helper<AT>::get_reference_point(pr2,traits).z(); }



  //As the previous proposal
};

template<typename GT, typename P>
typename AABB_traits_d<GT,P,3>::Axis
AABB_traits_d<GT,P,3>::longest_axis(const Bounding_box& bbox)
{
  const double dx = bbox.xmax() - bbox.xmin();
  const double dy = bbox.ymax() - bbox.ymin();
  const double dz = bbox.zmax() - bbox.zmin();
  if(dx>=dy)
  {
	  if(dx>=dz)
	  {
		  return CGAL_AXIS_X;
	  }

	  else
	  {
		  return CGAL_AXIS_Z;
	  }

  }
  else // dy>dx
  {
	  if(dy>=dz)
	  {
		  return CGAL_AXIS_Y;
	  }

	  else
	  {
		  return CGAL_AXIS_Z;
	  }
  }
}

}

#endif /* AABB_GEOMTRAITS_H_ */
