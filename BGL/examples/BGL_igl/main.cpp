#include <Eigen/Core>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <igl/triangle_triangle_adjacency.h>
#include <igl/readOFF.h>


namespace CGAL {
struct iglMesh;

template<typename T>
class IGL_Index
{
public:
  typedef boost::uint32_t size_type;
  /// Constructor. %Default construction creates an invalid index.
  /// We write -1, which is <a href="https://en.cppreference.com/w/cpp/types/numeric_limits">
  /// <tt>(std::numeric_limits<size_type>::max)()</tt></a>
  /// as `size_type` is an unsigned type.
  explicit IGL_Index(size_type _idx=(std::numeric_limits<size_type>::max)()) : idx_(_idx) {}

  /// Get the underlying index of this index
  operator size_type() const { return idx_; }

  /// reset index to be invalid (index=(std::numeric_limits<size_type>::max)())
  void reset() { idx_=(std::numeric_limits<size_type>::max)(); }

  /// return whether the index is valid, i.e., the index is not equal to `%std::numeric_limits<size_type>::max()`.
  bool is_valid() const {
    size_type inf = (std::numeric_limits<size_type>::max)();
    return idx_ != inf;
  }

  // Compatibility with OpenMesh handle
  size_type idx() const {
    return idx_;
  }
  // For convenience
  size_type id() const {
    return idx_;
  }

  /// increments the internal index. This operation does not
  /// guarantee that the index is valid or undeleted after the
  /// increment.
  IGL_Index& operator++() { ++idx_; return *this; }
  /// decrements the internal index. This operation does not
  /// guarantee that the index is valid or undeleted after the
  /// decrement.
  IGL_Index& operator--() { --idx_; return *this; }

  /// increments the internal index. This operation does not
  /// guarantee that the index is valid or undeleted after the
  /// increment.
  IGL_Index operator++(int) { IGL_Index tmp(*this); ++idx_; return tmp; }
  /// decrements the internal index. This operation does not
  /// guarantee that the index is valid or undeleted after the
  /// decrement.
  IGL_Index operator--(int) { IGL_Index tmp(*this); --idx_; return tmp; }

  IGL_Index operator+=(std::ptrdiff_t n) { idx_ = size_type(std::ptrdiff_t(idx_) + n); return *this; }

protected:
  size_type idx_;
};

template <class T>
std::size_t hash_value(const IGL_Index<T>&  i)
{
  std::size_t ret = i;
  return ret;
}



 class IGL_Vertex_index
   : public IGL_Index<IGL_Vertex_index>
 {
 public:

   IGL_Vertex_index() : IGL_Index<IGL_Vertex_index>((std::numeric_limits<size_type>::max)()) {}

   explicit IGL_Vertex_index(size_type _idx) : IGL_Index<IGL_Vertex_index>(_idx) {}

   template<class T> bool operator==(const T&) const = delete;
   template<class T> bool operator!=(const T&) const = delete;
   template<class T> bool operator<(const T&) const = delete;

   /// are two indices equal?
   bool operator==(const IGL_Vertex_index& _rhs) const {
     return this->idx_ == _rhs.idx_;
   }

   /// are two indices different?
   bool operator!=(const IGL_Vertex_index& _rhs) const {
     return this->idx_ != _rhs.idx_;
   }

   /// Comparison by index.
   bool operator<(const IGL_Vertex_index& _rhs) const {
     return this->idx_ < _rhs.idx_;
   }


   friend std::ostream& operator<<(std::ostream& os, IGL_Vertex_index const& v)
   {
     return (os << 'v' << (size_type)v );
   }
 };

 class IGL_Face_index
   : public IGL_Index<IGL_Face_index>
 {
 public:

   IGL_Face_index() : IGL_Index<IGL_Face_index>((std::numeric_limits<size_type>::max)()) {}

   explicit IGL_Face_index(size_type _idx) : IGL_Index<IGL_Face_index>(_idx) {}

   template<class T> bool operator==(const T&) const = delete;
   template<class T> bool operator!=(const T&) const = delete;
   template<class T> bool operator<(const T&) const = delete;

   /// are two indices equal?
   bool operator==(const IGL_Face_index& _rhs) const {
     return this->idx_ == _rhs.idx_;
   }

   /// are two indices different?
   bool operator!=(const IGL_Face_index& _rhs) const {
     return this->idx_ != _rhs.idx_;
   }

   /// Comparison by index.
   bool operator<(const IGL_Face_index& _rhs) const {
     return this->idx_ < _rhs.idx_;
   }


   friend std::ostream& operator<<(std::ostream& os, IGL_Face_index const& v)
   {
     return (os << 'f' << (size_type)v );
   }
 };


class IGL_Halfedge_iterator;

struct IGL_Halfedge_index
{
  typedef boost::uint32_t size_type;
  IGL_Face_index fi;
  int i;  //  0, 1, 2
  // bool is_border; // then fi and i are the opposite halfedge

  explicit IGL_Halfedge_index(size_type fi=(std::numeric_limits<size_type>::max)(), size_type i=0)
    : fi(fi), i(i)
  {}

  size_type idx() const
  {
    return fi.idx() + i;
  }

  bool operator==(const IGL_Halfedge_index h) const
  {
      return fi == h.fi && i == h.i;
  }

  bool operator!=(const IGL_Halfedge_index h) const
  {
      return fi != h.fi || i != h.i;
  }

  friend std::ostream& operator<<(std::ostream& os, IGL_Halfedge_index const& h)
  {
    return (os << 'h' << h.fi << " " << h.i);
  }
};


// First version only for meshes without border edges

class IGL_Halfedge_iterator
  : public boost::iterator_facade< IGL_Halfedge_iterator,
                                   IGL_Halfedge_index,
                                   std::random_access_iterator_tag,
                                   IGL_Halfedge_index
                                   >
{
  typedef boost::iterator_facade< IGL_Halfedge_iterator,
                                  IGL_Halfedge_index,
                                  std::random_access_iterator_tag,
                                  IGL_Halfedge_index> Facade;
public:
  IGL_Halfedge_iterator() : hnd_() {}
  IGL_Halfedge_iterator(const IGL_Halfedge_index& h)
    : hnd_(h)
  {}
private:
  friend class boost::iterator_core_access;
  void increment()
  {
    if(hnd_.i < 2){
      ++hnd_.i;
    }else{
      ++hnd_.fi;
      hnd_.i = 0;
    }
  }

  void decrement()
  {
    if(hnd_.i > 0){
      --hnd_.i;
    }else{
      --hnd_.fi;
      hnd_.i = 2;
    }
  }

  void advance(std::ptrdiff_t n)
  {
    // hnd_ += n;
    assert(false); // TBD
  }

  std::ptrdiff_t distance_to(const IGL_Halfedge_iterator& other) const
  {
    assert(false); // TBD
    return 0; // std::ptrdiff_t(other.hnd_) - std::ptrdiff_t(this->hnd_);
  }

  bool equal(const IGL_Halfedge_iterator& other) const
  {
    return this->hnd_ == other.hnd_;
  }

  IGL_Halfedge_index dereference() const { return hnd_; }

  IGL_Halfedge_index hnd_;

};

template<typename Index_>
class IGL_Index_iterator
  : public boost::iterator_facade< IGL_Index_iterator<Index_>,
                                   Index_,
                                   std::random_access_iterator_tag,
                                   Index_
                                   >
{
  typedef boost::iterator_facade< IGL_Index_iterator<Index_>,
                                  Index_,
                                  std::random_access_iterator_tag,
                                  Index_> Facade;
public:
  IGL_Index_iterator() : hnd_() {}
  IGL_Index_iterator(const Index_& h)
    : hnd_(h)
  {}
private:
  friend class boost::iterator_core_access;
  void increment()
  {
    ++hnd_;
  }

  void decrement()
  {
    --hnd_;
  }

  void advance(std::ptrdiff_t n)
  {
    hnd_ += n;
  }

  std::ptrdiff_t distance_to(const IGL_Index_iterator& other) const
  {
    return std::ptrdiff_t(other.hnd_) - std::ptrdiff_t(this->hnd_);
  }

  bool equal(const IGL_Index_iterator& other) const
  {
    return this->hnd_ == other.hnd_;
  }

  Index_ dereference() const { return hnd_; }

  Index_ hnd_;

};



struct iglMesh
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  Eigen::MatrixXi TT, TTi;

  typedef IGL_Vertex_index Vertex_index;
  typedef IGL_Index_iterator<Vertex_index> Vertex_iterator;
  typedef CGAL::Iterator_range<Vertex_iterator> Vertex_range;

  typedef IGL_Face_index Face_index;
  typedef IGL_Index_iterator<Face_index> Face_iterator;
  typedef CGAL::Iterator_range<Face_iterator> Face_range;

  typedef IGL_Halfedge_index Halfedge_index;
  typedef IGL_Halfedge_iterator Halfedge_iterator;
  typedef CGAL::Iterator_range<Halfedge_iterator> Halfedge_range;
};

} // namespace CGAL

namespace boost {
template <>
struct graph_traits<CGAL::iglMesh> {
  typedef typename CGAL::iglMesh::Vertex_index vertex_descriptor;
  typedef typename CGAL::iglMesh::Vertex_iterator vertex_iterator;
  typedef typename CGAL::iglMesh::Face_index face_descriptor;
  typedef typename CGAL::iglMesh::Face_iterator face_iterator;
  typedef typename CGAL::iglMesh::Halfedge_index halfedge_descriptor;
  typedef typename CGAL::iglMesh::Halfedge_iterator halfedge_iterator;

  static vertex_descriptor   null_vertex() { return vertex_descriptor(); }
  static face_descriptor     null_face()   { return face_descriptor(); }
  static halfedge_descriptor null_halfedge()   { return halfedge_descriptor(); }
};

} // namespace boost


namespace CGAL {

int num_vertices(const iglMesh& m)
{
  return m.V.size()/3;
}

int num_faces(const iglMesh& m)
{
  return m.F.size()/3;
}

int num_halfedges(const iglMesh& m)
{
  return m.F.size();
}

typename iglMesh::Vertex_range
vertices(const iglMesh& m)
{
  typedef typename boost::graph_traits<iglMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<iglMesh>::vertex_iterator vertex_iterator;
  return std::make_pair(vertex_iterator(vertex_descriptor(0)), vertex_iterator(vertex_descriptor(num_vertices(m))));
}

typename iglMesh::Face_range
faces(const iglMesh& m)
{
  typedef typename boost::graph_traits<iglMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<iglMesh>::face_iterator face_iterator;
  return std::make_pair(face_iterator(face_descriptor(0)), face_iterator(face_descriptor(num_faces(m))));
}

typename iglMesh::Halfedge_range
halfedges(const iglMesh& m)
{
  typedef typename boost::graph_traits<iglMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<iglMesh>::halfedge_iterator halfedge_iterator;
  return std::make_pair(halfedge_iterator(halfedge_descriptor(0,0)), halfedge_iterator(halfedge_descriptor(num_faces(m),0)));
}


typename boost::graph_traits<iglMesh >::vertex_descriptor
source(typename boost::graph_traits<iglMesh >::halfedge_descriptor h,
       const iglMesh& m)
{
    typedef typename boost::graph_traits<iglMesh >::vertex_descriptor vertex_descriptor;
    std::size_t fi = h.fi.idx();
    std::size_t i = h.i;
  return vertex_descriptor(m.F(fi,i));
}

typename boost::graph_traits<iglMesh >::vertex_descriptor
target(typename boost::graph_traits<iglMesh >::halfedge_descriptor h,
       const iglMesh& m)
{
    typedef typename boost::graph_traits<iglMesh >::vertex_descriptor vertex_descriptor;
    std::size_t fi = h.fi.idx();
    std::size_t i = (h.i == 2)?0:h.i+1;
  return vertex_descriptor(m.F(fi,i));
}

typename boost::graph_traits<iglMesh >::halfedge_descriptor
next(typename boost::graph_traits<iglMesh >::halfedge_descriptor h,
       const iglMesh& m)
{
    typedef typename boost::graph_traits<iglMesh >::halfedge_descriptor halfedge_descriptor;
    std::size_t fi = h.fi.idx();
    std::size_t i = (h.i == 2)?0:h.i+1;
    return halfedge_descriptor(fi,i);
}

typename boost::graph_traits<iglMesh >::halfedge_descriptor
prev(typename boost::graph_traits<iglMesh >::halfedge_descriptor h,
       const iglMesh& m)
{
    typedef typename boost::graph_traits<iglMesh >::halfedge_descriptor halfedge_descriptor;
    std::size_t fi = h.fi.idx();
    std::size_t i = (h.i == 0)?2:h.i-1;
    return halfedge_descriptor(fi,i);
}

typename boost::graph_traits<iglMesh >::face_descriptor
face(typename boost::graph_traits<iglMesh >::halfedge_descriptor h,
       const iglMesh& m)
{
  typedef typename boost::graph_traits<iglMesh >::face_descriptor face_descriptor;
  return face_descriptor(h.fi);
}

typename boost::graph_traits<iglMesh >::halfedge_descriptor
halfedge(typename boost::graph_traits<iglMesh >::face_descriptor f,
       const iglMesh& m)
{
  typedef typename boost::graph_traits<iglMesh >::halfedge_descriptor halfedge_descriptor;
  return halfedge_descriptor(f, 0);
}

typename boost::graph_traits<iglMesh >::halfedge_descriptor
opposite(typename boost::graph_traits<iglMesh >::halfedge_descriptor h,
       const iglMesh& m)
{
  typedef typename boost::graph_traits<iglMesh >::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<iglMesh >::face_descriptor face_descriptor;

  std::size_t fi = h.fi.idx();
  std::size_t i = h.i;
  return halfedge_descriptor(face_descriptor(m.TT.row(fi)[i]), m.TTi.row(fi)[i]);
}



template <typename VEF>
class IGL_Index_pmap
{
public:
  typedef boost::readable_property_map_tag category;
  typedef boost::uint32_t                  value_type;
  typedef boost::uint32_t                  reference;
  typedef VEF                              key_type;

  value_type operator[](const key_type& vd) const
  {
    return vd;
  }

  friend inline
  value_type get(const IGL_Index_pmap& m, const key_type& k)
  {
    return m[k];
  }
};




template<typename Mesh, typename P>
class IGL_point_pmap
{
public:
  typedef boost::read_write_property_map_tag category;
  typedef P value_type;
  typedef P reference;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor key_type;

  IGL_point_pmap()
    : sm_(nullptr)
  {}

  IGL_point_pmap(const Mesh& sm)
    : sm_(&sm)
    {}

  IGL_point_pmap(const IGL_point_pmap& pm)
    : sm_(pm.sm_)
    {}

  reference operator[](key_type v) const
  {
    CGAL_assertion(sm_!=nullptr);
    const auto& xyz = sm_->V.row(v.idx());
    return value_type(xyz[0], xyz[1], xyz[2]);
  }

  inline friend reference get(const IGL_point_pmap<Mesh,P>& pm, key_type v)
  {
    CGAL_precondition(pm.sm_!=nullptr);
    const auto& xyz = pm.sm_->V.row(v.idx());
    return value_type(xyz[0], xyz[1], xyz[2]);
  }

  inline friend void put(const IGL_point_pmap<Mesh,P>& pm, key_type v, const value_type& p)
  {
    CGAL_precondition(pm.sm_!=nullptr);

    //typedef typename OpenMesh::vector_traits<typename OM_Mesh::Point>::value_type Scalar;
    //const_cast<OM_Mesh&>(*pm.sm_).set_point
    //  (v, typename OM_Mesh::Point(Scalar(p[0]), Scalar(p[1]), Scalar(p[2])));
  }

  private:
  const Mesh* sm_;
};

template <>
struct graph_has_property<iglMesh, boost::vertex_index_t>
  : CGAL::Tag_true{};

template <>
struct graph_has_property<iglMesh, boost::face_index_t>
  : CGAL::Tag_true{};

template <>
struct graph_has_property<iglMesh, boost::halfedge_index_t>
  : CGAL::Tag_true{};

template <>
struct graph_has_property<iglMesh, boost::vertex_point_t>
  : CGAL::Tag_true{};

} // namespace CGAL



namespace boost{

template <>
struct property_map<CGAL::iglMesh, boost::vertex_point_t >
{
  typedef typename CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 P;
  typedef CGAL::IGL_point_pmap  <CGAL::iglMesh,P> type;
  typedef CGAL::IGL_point_pmap<CGAL::iglMesh,P> const_type;
};

template <>
struct property_map<CGAL::iglMesh, boost::vertex_index_t >
{
  typedef CGAL::IGL_Index_pmap<typename boost::graph_traits<CGAL::iglMesh >::vertex_descriptor> type;
  typedef CGAL::IGL_Index_pmap<typename boost::graph_traits<CGAL::iglMesh >::vertex_descriptor> const_type;
};

template <>
struct property_map<CGAL::iglMesh, boost::halfedge_index_t >
{
  typedef CGAL::IGL_Index_pmap<typename boost::graph_traits<CGAL::iglMesh >::halfedge_descriptor> type;
  typedef CGAL::IGL_Index_pmap<typename boost::graph_traits<CGAL::iglMesh >::halfedge_descriptor> const_type;
};

template <>
struct property_map<CGAL::iglMesh, boost::face_index_t >
{
  typedef CGAL::IGL_Index_pmap<typename boost::graph_traits<CGAL::iglMesh >::face_descriptor> type;
  typedef CGAL::IGL_Index_pmap<typename boost::graph_traits<CGAL::iglMesh >::face_descriptor> const_type;
};
} // namespace boost


namespace CGAL{

  //template <typename Point>
IGL_Index_pmap<IGL_Vertex_index>
get(const boost::vertex_index_t&, const iglMesh&)
{
  return IGL_Index_pmap<typename boost::graph_traits<iglMesh >::vertex_descriptor>();
}

IGL_Index_pmap<IGL_Face_index>
get(const boost::face_index_t&, const iglMesh&)
{
  return IGL_Index_pmap<typename boost::graph_traits<iglMesh >::face_descriptor>();
}

IGL_Index_pmap<typename boost::graph_traits<iglMesh >::halfedge_descriptor>
get(const boost::halfedge_index_t&, const iglMesh&)
{
  return IGL_Index_pmap<typename boost::graph_traits<iglMesh >::halfedge_descriptor>();
}

// template<typename K>
IGL_point_pmap<iglMesh,
                    typename CGAL::Exact_predicates_inexact_constructions_kernel::Point_3>
get(boost::vertex_point_t, const iglMesh& g)
{
  typedef typename CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 P;
  return IGL_point_pmap<iglMesh, P>(g);
}

} // namespace boost


int main()
{
  CGAL::iglMesh mesh;

  igl::readOFF(  "tet.off", mesh.V, mesh.F);
  igl:: triangle_triangle_adjacency(mesh.F,mesh.TT,mesh.TTi);

  auto pmap = get(CGAL::vertex_point, mesh);
  auto ipmap = get(boost::vertex_index, mesh);
  for(auto v : vertices(mesh)){
    std::cout << v << std::endl;
    std::cout << get(ipmap, v) << "  " << get(pmap, v) << std::endl;
  }

  for(auto f : faces(mesh)){
    std::cout << f << std::endl;
  }

  for(auto h : halfedges(mesh)){
    std::cout << "h = " << h << std::endl;
    std::cout << "opposite(h) = " << opposite(h,mesh) << std::endl;
    std::cout << "o(o(h)) = " << opposite(opposite(h,mesh),mesh) << std::endl;
    std::cout << source(h,mesh) << std::endl;
  }
  std::cout << "#V = " << num_vertices(mesh) << std::endl;




  // std::cout << mesh.TTi.row(0) << std::endl;

  std::cout << "area = " << CGAL::Polygon_mesh_processing::area(mesh, CGAL::parameters::geom_traits(CGAL::Exact_predicates_inexact_constructions_kernel())
  //    .vertex_point_map(pmap)
  ) << std::endl;
  return 0;
}
