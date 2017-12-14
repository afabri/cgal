// Copyright (c) 2017  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Maxime Gimeno

// List of named parameters that we use in CGAL
CGAL_add_named_parameter(vertex_point_t, vertex_point, vertex_point_map)
CGAL_add_named_parameter(halfedge_index_t, halfedge_index, halfedge_index_map)
CGAL_add_named_parameter(edge_index_t, edge_index, edge_index_map)
CGAL_add_named_parameter(face_index_t, face_index, face_index_map)

CGAL_add_named_parameter(edge_is_constrained_t, edge_is_constrained, edge_is_constrained_map)
CGAL_add_named_parameter(first_index_t, first_index, first_index)
CGAL_add_named_parameter(number_of_iterations_t, number_of_iterations, number_of_iterations)

// List of named parameters that we use in the package 'Mesh_3'
CGAL_add_named_parameter(vertex_feature_degree_t, vertex_feature_degree, vertex_feature_degree_map)

// List of named parameters used in the package 'Polygon Mesh Processing'
CGAL_add_named_parameter(geom_traits_t, geom_traits, geom_traits)
CGAL_add_named_parameter(vertex_incident_patches_t, vertex_incident_patches, vertex_incident_patches_map)
CGAL_add_named_parameter(density_control_factor_t, density_control_factor, density_control_factor)
CGAL_add_named_parameter(use_delaunay_triangulation_t, use_delaunay_triangulation, use_delaunay_triangulation)
CGAL_add_named_parameter(fairing_continuity_t, fairing_continuity, fairing_continuity)
CGAL_add_named_parameter(sparse_linear_solver_t, sparse_linear_solver, sparse_linear_solver)
CGAL_add_named_parameter(number_of_relaxation_steps_t, number_of_relaxation_steps, number_of_relaxation_steps)
CGAL_add_named_parameter(protect_constraints_t, protect_constraints, protect_constraints)
CGAL_add_named_parameter(relax_constraints_t, relax_constraints, relax_constraints)
CGAL_add_named_parameter(vertex_is_constrained_t, vertex_is_constrained, vertex_is_constrained_map)
CGAL_add_named_parameter(face_patch_t, face_patch, face_patch_map)
CGAL_add_named_parameter(random_uniform_sampling_t, random_uniform_sampling, use_random_uniform_sampling)
CGAL_add_named_parameter(grid_sampling_t, grid_sampling, use_grid_sampling)
CGAL_add_named_parameter(monte_carlo_sampling_t, monte_carlo_sampling, use_monte_carlo_sampling)
CGAL_add_named_parameter(do_sample_edges_t, do_sample_edges, do_sample_edges)
CGAL_add_named_parameter(do_sample_vertices_t, do_sample_vertices, do_sample_vertices)
CGAL_add_named_parameter(do_sample_faces_t, do_sample_faces, do_sample_faces)
CGAL_add_named_parameter(number_of_points_on_faces_t, number_of_points_on_faces, number_of_points_on_faces)
CGAL_add_named_parameter(number_of_points_per_face_t, number_of_points_per_face, number_of_points_per_face)
CGAL_add_named_parameter(grid_spacing_t, grid_spacing, grid_spacing)
CGAL_add_named_parameter(number_of_points_per_edge_t, number_of_points_per_edge, number_of_points_per_edge)
CGAL_add_named_parameter(number_of_points_on_edges_t, number_of_points_on_edges, number_of_points_on_edges)
CGAL_add_named_parameter(nb_points_per_area_unit_t, nb_points_per_area_unit, number_of_points_per_area_unit)
CGAL_add_named_parameter(nb_points_per_distance_unit_t, nb_points_per_distance_unit, number_of_points_per_distance_unit)
CGAL_add_named_parameter(face_partition_t, face_partition, face_partition_map)
CGAL_add_named_parameter(partition_size_t, partition_size, partition_size)

// List of named parameters that we use in the package 'Surface Mesh Simplification'
CGAL_add_named_parameter(get_cost_policy_t, get_cost_policy, get_cost)
CGAL_add_named_parameter(get_placement_policy_t, get_placement_policy, get_placement)
CGAL_add_named_parameter(current_num_edges_t, current_num_edges, current_num_edges)

//to be documented
CGAL_add_named_parameter(face_normal_t, face_normal, face_normal_map)
CGAL_add_named_parameter(random_seed_t, random_seed, random_seed)
CGAL_add_named_parameter(do_project_t, do_project, do_project)

//internal
CGAL_add_named_parameter(weight_calculator_t, weight_calculator, weight_calculator)
