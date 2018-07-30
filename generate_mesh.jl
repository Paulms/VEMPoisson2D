function _generate_2d_nodes!(nodes, nx, ny, LL, LR, UR, UL)
      for i in 0:ny-1
        ratio_bounds = i / (ny-1)

        x0 = LL[1] * (1 - ratio_bounds) + ratio_bounds * UL[1]
        x1 = LR[1] * (1 - ratio_bounds) + ratio_bounds * UR[1]

        y0 = LL[2] * (1 - ratio_bounds) + ratio_bounds * UL[2]
        y1 = LR[2] * (1 - ratio_bounds) + ratio_bounds * UR[2]

        for j in 0:nx-1
            ratio = j / (nx-1)
            x = x0 * (1 - ratio) + ratio * x1
            y = y0 * (1 - ratio) + ratio * y1
            push!(nodes, Node(Vec{2}((x, y))))
        end
    end
end

function _build_cells(face_idx, cells::Vector{Cell{2,N,M}}, el_nodes, el_faces,n_el, faces, facesdict,nodes) where {N,M}
    fill!(el_faces,0)
    #build faces
    for (i,fn_id) in enumerate(reference_edge_nodes(eltype(cells)))
        v1 = el_nodes[fn_id[1]]; v2 = el_nodes[fn_id[2]]
        element = minmax(v1,v2)
        token = Base.ht_keyindex2!(facesdict, element)
        if token > 0
                el_faces[i] = facesdict.vals[token]
                if n_el != faces[facesdict.vals[token],4]
                    faces[facesdict.vals[token],4] = n_el
                end
        else
            face_idx = face_idx + 1
            Base._setindex!(facesdict, face_idx, element, -token)
            faces[face_idx,1] = v1
            faces[face_idx,2] = v2
            faces[face_idx,3] = n_el
            faces[face_idx,4] = 0
            el_faces[i] = face_idx
        end
    end
    #save cell
    cell = _build_cell(eltype(cells),el_nodes,el_faces)
    push!(cells, cell)
    return face_idx
end

#Boundary sets for rectangular mesh
function _get_rectangular_boundary_sets(nodes::Vector{Node{2,T}}, faces::Matrix{Int}, face_idx::Int,
    LL::Vec{2,T}, UR::Vec{2,T}, LR::Vec{2,T}, UL::Vec{2,T}) where {T}
    # Add faces and nodes sets
    # Add faces and nodes sets
    bottomSet = Set{Int}(); bottomNSet = Set{Int}()
    rightSet = Set{Int}(); rightNSet = Set{Int}()
    topSet = Set{Int}(); topNSet = Set{Int}()
    leftSet = Set{Int}(); leftNSet = Set{Int}()
    for k in 1:face_idx
        ref = 0
        if faces[k,4] == 0
            if (nodes[faces[k,1]].x[2] == LL[2] && nodes[faces[k,2]].x[2] == LL[2])
                push!(bottomSet, k); push!(bottomNSet, faces[k,1]); push!(bottomNSet, faces[k,2])
            elseif (nodes[faces[k,1]].x[1] == UR[1] && nodes[faces[k,2]].x[1] == UR[1])
                push!(rightSet, k); push!(rightNSet, faces[k,1]); push!(rightNSet, faces[k,2])
            elseif (nodes[faces[k,1]].x[2] == UR[2] && nodes[faces[k,2]].x[2] == UR[2])
                push!(topSet, k); push!(topNSet, faces[k,1]); push!(topNSet, faces[k,2])
            elseif (nodes[faces[k,1]].x[1] == LL[1] && nodes[faces[k,2]].x[1] == LL[1])
                push!(leftSet, k); push!(leftNSet, faces[k,1]); push!(leftNSet, faces[k,2])
            else
                throw("Face $k belongs to one cell but is not in boundary")
            end
        end
    end
    facesets = Dict("bottom"=>bottomSet, "right"=>rightSet,"left"=>leftSet,
                    "top"=>topSet,"boundary"=>union(bottomSet,rightSet,leftSet,topSet))
    nodesets = Dict("bottom"=>bottomNSet, "right"=>rightNSet,"left"=>leftNSet,
                    "top"=>topNSet,"boundary"=>union(bottomNSet,rightNSet,leftNSet,topNSet))
    return facesets, nodesets
end

#########################
# Triangle Cells 2D   #
#########################
@inline _build_cell(::Type{TriangleCell}, el_nodes, el_faces) = TriangleCell(el_nodes,(el_faces[1],el_faces[2],el_faces[3]))
"""
rectangle_mesh(::Type{TriangleCell}, nel::NTuple{2,Int}, LL::Vec{2,T}, UR::Vec{2,T}) where {T}
Generate a rectangular mesh with triangular cells, where `LL` is the low left vertex
and `UR` is the upper right one. `nel` is a tuple with the number of partions to be
used in each dimension.
"""
function rectangle_mesh(::Type{TriangleCell}, nel::NTuple{2,Int}, LL::Vec{2,T}, UR::Vec{2,T}) where {T}
    LR = Vec{2}((UR[1],LL[2]))
    UL = Vec{2}((LL[1],UR[2]))
    nel_x = nel[1]; nel_y = nel[2]; nel_tot = 2*nel_x*nel_y
    n_nodes_x = nel_x + 1; n_nodes_y = nel_y + 1
    n_nodes = n_nodes_x * n_nodes_y

    # Generate nodes
    nodes = Node{2,T}[]
    _generate_2d_nodes!(nodes, n_nodes_x, n_nodes_y, LL, LR, UR, UL)

    # Generate cells
    faces = Matrix{Int}(undef,nel[1]+nel[2]+3*nel[1]*nel[2],4)
    node_array = reshape(collect(1:n_nodes), (n_nodes_x, n_nodes_y))
    cells = TriangleCell[]
    facesdict = Dict{NTuple{2,Int},Int}()
    el_faces = [0,0,0]
    n_el = 0
    face_idx = 0
    for j in 1:nel_y, i in 1:nel_x
        n_el = n_el + 1
        el_nodes = (node_array[i,j], node_array[i+1,j], node_array[i,j+1])
        face_idx = _build_cells(face_idx, cells, el_nodes, el_faces,n_el, faces, facesdict,nodes) # ◺

        n_el = n_el + 1
        el_nodes = (node_array[i+1,j], node_array[i+1,j+1], node_array[i,j+1])
        face_idx = _build_cells(face_idx,cells, el_nodes, el_faces,n_el, faces, facesdict,nodes) # ◹
    end

    facesets, nodesets = _get_rectangular_boundary_sets(nodes, faces, face_idx, LL, UR, LR, UL)
    return PolygonalMesh{2,3,3,2,T}(cells, nodes, faces, facesets, nodesets)
end

#########################
# Rectangle Cells 2D   #
#########################
@inline _build_cell(::Type{RectangleCell}, el_nodes, el_faces) = RectangleCell(el_nodes,(el_faces[1],el_faces[2],el_faces[3],el_faces[4]))
"""
rectangle_mesh(::Type{RectangleCell}, nel::NTuple{2,Int}, LL::Vec{2,T}, UR::Vec{2,T}) where {T}
Generate a rectangular mesh with triangular cells, where `LL` is the low left vertex
and `UR` is the upper right one. `nel` is a tuple with the number of partions to be
used in each dimension.
"""
function rectangle_mesh(::Type{RectangleCell}, nel::NTuple{2,Int}, LL::Vec{2,T}, UR::Vec{2,T}) where {T}
    LR = Vec{2}((UR[1],LL[2]))
    UL = Vec{2}((LL[1],UR[2]))
    nel_x = nel[1]; nel_y = nel[2]; nel_tot = 2*nel_x*nel_y
    n_nodes_x = nel_x + 1; n_nodes_y = nel_y + 1
    n_nodes = n_nodes_x * n_nodes_y

    # Generate nodes
    nodes = Node{2,T}[]
    _generate_2d_nodes!(nodes, n_nodes_x, n_nodes_y, LL, LR, UR, UL)

    faces = Matrix{Int}(undef,(nel[1]+1)*nel[2]+nel[1]*(nel[2]+1),4)

    # Generate cells
    node_array = reshape(collect(1:n_nodes), (n_nodes_x, n_nodes_y))
    cells = RectangleCell[]
    facesdict = Dict{NTuple{2,Int},Int}()
    el_faces = [0,0,0,0]
    n_el = 0
    face_idx = 0
    for j in 1:nel_y, i in 1:nel_x
        n_el = n_el + 1
        el_nodes = (node_array[i,j], node_array[i+1,j], node_array[i+1,j+1], node_array[i,j+1])
        face_idx = _build_cells(face_idx, cells, el_nodes, el_faces,n_el, faces, facesdict,nodes)
    end

    # Add faces and nodes sets
    facesets, nodesets = _get_rectangular_boundary_sets(nodes, faces, face_idx, LL, UR, LR, UL)
    return PolygonalMesh{2,4,4,2,T}(cells, nodes, faces, facesets, nodesets)
end
