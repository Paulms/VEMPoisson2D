# VEM computes the virtual element solution of a Poisson problem on a polygonal mesh
#
using Tensors
using SparseArrays
using LinearAlgebra
import Base:@propagate_inbounds
# Read mesh
include("mesh.jl")
include("generate_mesh.jl")
include("assembler.jl")

@time mesh = rectangle_mesh(RectangleCell, (10,10), Vec{2}((0.0,0.0)), Vec{2}((1.0,1.0)));

# forcing function
rhs(x::Vec{2}) = 15 * sin(π * x[1]) * sin(π * x[2]);
# Boundary condition
g(x::Vec{2}) = (1 - x[1])*x[2]*sin(π * x[1]);
n_dofs = getnnodes(mesh)
n_polys = 3; # Method has 1 degree of freedom per vertex
assembler = start_assemble(n_dofs)
F = zeros(n_dofs); # Forcing vector
u = zeros(n_dofs); # Degrees of freedom of the virtual element solution
linear_polynomials = ((0,0), (1,0), (0,1)); # Impose an ordering on the linear polynomials
for el_id = 1:getncells(mesh)
	vert_ids = mesh.cells[el_id].nodes # Global IDs of the vertices of this element
	verts = get_coordinates(mesh.cells[el_id],  mesh) # Coordinates of the vertices of this element
	n_sides = getnfaces(mesh.cells[el_id]) # Start computing the geometric information
	area = cell_volume(el_id, mesh)
	centr = cell_centroid(el_id, mesh)
	diameter = cell_diameter(el_id, mesh)
	D = zeros(n_sides, n_polys); D[:, 1] .= 1;
	B = zeros(n_polys, n_sides); B[1, :] .= 1/n_sides;
	for vertex_id = 1:n_sides
		vert = verts[vertex_id]; #This vertex and its neighbours
		prev = verts[mod1(vertex_id - 1, n_sides)];
		nextv = verts[mod1(vertex_id + 1, n_sides)];
		vertex_normal = [nextv[2] - prev[2], prev[1] - nextv[1]]; # Average of edge normals
		for poly_id in 2:n_polys # Only need to loop over non-constant polynomials
			poly_degree = linear_polynomials[poly_id];
			monomial_grad = poly_degree ./ diameter; # Gradient of a linear is constant
			D[vertex_id, poly_id] = dot(vert - centr, poly_degree) / diameter;
			B[poly_id, vertex_id] = 0.5 * dot(monomial_grad, vertex_normal);
		end
	end
	projector = (B*D) \ B; # Compute the local Ritz projector to polynomials
	stabilising_term = (I - D * projector)' * (I - D * projector);
	G = B*D; G[1, :] .= 0;
	local_stiffness = projector' * G * projector + stabilising_term;
	assemble!(assembler,vert_ids,local_stiffness)
	assemble!(F, [dof for dof in vert_ids], rhs(centr)*ones(length(vert_ids)) * area / n_sides)
end
K = end_assemble(assembler)
boundary_nodes = getnodeset(mesh,"boundary")
boundary_vals = Vector{Float64}(undef, length(boundary_nodes))
for (i,node) in enumerate(boundary_nodes)
	boundary_vals[i] = g(mesh.nodes[node].x)
end

internal_dofs = setdiff(1:n_dofs, boundary_nodes)
F = F - K[:, collect(boundary_nodes)] * boundary_vals; # Apply the boundary condition
u[internal_dofs] = K[collect(internal_dofs), collect(internal_dofs)] \ F[collect(internal_dofs)]; # Solve
u[collect(boundary_nodes)] = boundary_vals; # Set the boundary values

#plot solution
# coordinates = get_vertices_matrix(mesh);
# connectivity = get_cells_matrix(mesh);
# using Makie
# poly(coordinates, connectivity, color = u, strokecolor = (:black, 0.6), strokewidth = 4)
