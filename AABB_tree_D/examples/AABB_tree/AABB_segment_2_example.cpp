// Author : 

#include <iostream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>


typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Segment_2 Segment;
typedef K::Point_2 Point;

typedef std::list<Segment>::iterator Iterator;
typedef CGAL::AABB_segment_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;


int main()
{
	Point a(0.0, 0.0);
	Point b(2.0, 1.0);
	Point c(3.0, 4.0);
	Point d(1.0, 6.0);
	Point e(-1.0, 3.0);

	std::list<Segment> segments;
	segments.push_back(Segment(a, b));
	segments.push_back(Segment(b, c));
	segments.push_back(Segment(c, d));
	segments.push_back(Segment(d, e));
	segments.push_back(Segment(e, a));

	//   // constructs the AABB tree and the internal search tree for 
	//   // efficient distance computations.
	Tree tree(segments.begin(), segments.end());
	tree.accelerate_distance_queries();

	//   // counts #intersections with a segment query
	Segment segment_query(Point(1.0, 0.0), Point(0.0, 7.0));
	std::cout << tree.number_of_intersected_primitives(segment_query)
		<< " intersections(s) with segment" << std::endl;


	// computes the closest point from a point query 
	Point point_query(1.5, 3.0);
	Point closest = tree.closest_point(point_query);

	std::cerr << "closest point is: " << closest << std::endl;
	return EXIT_SUCCESS;
}
