module ModMesh
# Creating the mesh
#
# Hongyang Zhou, hyzhou@umich.edu 05/12/2019

export Mesh, createMesh, createTriMesh

using Triangle

const ϵ = 1e-5
const East_,West_,North_,South_ = 1,2,3,4
const INF = 10000.0 # For determining point inside polygon or not

struct Mesh
   Connectivity::Array{Int,2}
   BClines::Array{Float64,2}
   BCDir::Vector{Int}
   Coords::Array{Float64,2}
   BCnode::Vector{Vector{Int}}
end

"""
Given three colinear points p, q, r, the function checks if
point q lies on line segment 'pr'
"""
function onSegment(p::Vector{Float64}, q::Vector{Float64}, r::Vector{Float64})

   if (q[1] ≤ max(p[1], r[1]) && q[1] ≥ min(p[1], r[1]) &&
      q[2] ≤ max(p[2], r[2]) && q[2] ≥ min(p[2], r[2]))
      return true
   else
      return false
   end

end

"""
To find orientation of ordered triplet (p, q, r).
The function returns following values
0 --> p, q and r are colinear
1 --> clockwise
2 --> counterclockwise
"""
function orientation(p::Vector{Float64}, q::Vector{Float64}, r::Vector{Float64})

   val = (q[2] - p[2])*(r[1] - q[1]) - (q[1] - p[1])*(r[2] - q[2])

   if val == 0
      return 0 # colinear
   else
      return val > 0 ? 1 : 2 # clockwise or counter-clockwise
   end
end

function doIntersect(p1::Vector{Float64},q1::Vector{Float64},
   p2::Vector{Float64}, q2::Vector{Float64})
   # Find the four orientations needed for general and special cases
   o1 = orientation(p1, q1, p2)
   o2 = orientation(p1, q1, q2)
   o3 = orientation(p2, q2, p1)
   o4 = orientation(p2, q2, q1)

   # General case
   if o1 != o2 && o3 != o4 return true end

   # Special Cases
   # p1, q1 and p2 are colinear and p2 lies on segment p1q1
   if o1 == 0 && onSegment(p1, p2, q1) return true end

   # p1, q1 and p2 are colinear and q2 lies on segment p1q1
   if o2 == 0 && onSegment(p1, q2, q1) return true end

   # p2, q2 and p1 are colinear and p1 lies on segment p2q2
   if o3 == 0 && onSegment(p2, p1, q2) return true end

   # p2, q2 and q1 are colinear and q1 lies on segment p2q2
   if o4 == 0 && onSegment(p2, q1, q2) return true end

   return false # Doesn't fall in any of the above cases
end

"Returns true if the point p lies inside the polygon[] with n vertices."
function isInside(point::Vector{Float64},BClines::Array{Float64,2})

   count = 0
   nlines = size(BClines)[2]
   # loop over all vertical boundary lines
   for i=1:nlines
      line = BClines[:,i]
      if doIntersect(line[1:2],line[3:4],point,[INF, point[2]]) count += 1 end
   end

   return isodd(count)
end


# Rules:
# 1. boundaries are ordered counter-clockwisely.
function createMesh(BClines::Array{Float64,2},nI::Int, nJ::Int,
   LatticeType::Int)

   nlines = size(BClines,2)
   BCDir = zeros(Int,nlines)
   Coords = Array{Float64}(undef, 0, 0)

   # Check the orientation of linestyle
   for i = 1:nlines
      if abs(BClines[1,i] - BClines[3,i]) < ϵ
         BClines[2,i] > BClines[4,i] ? BCDir[i] = West_ : BCDir[i] = East_
      else
         BClines[1,i] > BClines[3,i] ? BCDir[i] = North_ : BCDir[i] = South_
      end
   end

   # Find the outermost locations
   xmin,xmax = extrema(BClines[:,1:2:3])
   ymin,ymax = extrema(BClines[:,2:2:4])

   # Shift a little bit to avoid issue at the boundary
   xmesh = range(xmin+ϵ,stop=xmax-ϵ,length=nI)
   ymesh = range(ymin+ϵ,stop=ymax-ϵ,length=nJ)

   dx = xmesh[2] - xmesh[1]
   dy = ymesh[2] - ymesh[1]

   # Check if the point lies inside the region
   for x in xmesh
      for y in ymesh
         if isInside([x,y],BClines)
            # Add to list
            try
               Coords = [Coords [x,y]]
            catch
               Coords = [x,y]
            end
         end
      end
   end

   # Create Connectivity list
   nNode = size(Coords,2)
   if LatticeType == 1 || LatticeType == 2
      Connectivity = fill(-1, (4,nNode))
   else
      Connectivity = fill(-1, (8,nNode))
   end
   for i = 1:nNode
      # East
      iNeiEast =
      findfirst([ abs(Coords[1,i]+dx-Coords[1,j]) < ϵ &&
      abs(Coords[2,i]-Coords[2,j]) < ϵ for j=1:nNode])

      # West
      iNeiWest =
      findfirst([ abs(Coords[1,i]-dx-Coords[1,j]) < ϵ &&
      abs(Coords[2,i]-Coords[2,j]) < ϵ for j=1:nNode])

      # South
      iNeiSouth =
      findfirst([ abs(Coords[1,i]-Coords[1,j]) < ϵ &&
      abs(Coords[2,i]-Coords[2,j]-dy) < ϵ for j=1:nNode])

      # North
      iNeiNorth =
      findfirst([ abs(Coords[1,i]-Coords[1,j]) < ϵ &&
      abs(Coords[2,i]-Coords[2,j]+dy) < ϵ for j=1:nNode])

      # SouthEast
      iNeiSouthEast =
      findfirst([ abs(Coords[1,i]+dx-Coords[1,j]) < ϵ &&
      abs(Coords[2,i]-Coords[2,j]-dy) < ϵ for j=1:nNode])

      # SouthWest
      iNeiSouthWest =
      findfirst([ abs(Coords[1,i]-dx-Coords[1,j]) < ϵ &&
      abs(Coords[2,i]-Coords[2,j]-dy) < ϵ for j=1:nNode])

      # NorthEast
      iNeiNorthEast =
      findfirst([ abs(Coords[1,i]+dx-Coords[1,j]) < ϵ &&
      abs(Coords[2,i]-Coords[2,j]+dy) < ϵ for j=1:nNode])

      # NorthWest
      iNeiNorthWest =
      findfirst([ abs(Coords[1,i]-dx-Coords[1,j]) < ϵ &&
      abs(Coords[2,i]-Coords[2,j]+dy) < ϵ for j=1:nNode])

      if iNeiWest  != nothing Connectivity[1,i] = iNeiWest end
      if iNeiEast  != nothing Connectivity[2,i] = iNeiEast end
      if iNeiSouth != nothing Connectivity[3,i] = iNeiSouth end
      if iNeiNorth != nothing Connectivity[4,i] = iNeiNorth end

      #
      if iNeiSouthEast == nothing &&
         iNeiSouth != nothing && iNeiEast != nothing
         Connectivity[[2,3],i] .= -1
      end
      if iNeiSouthWest == nothing &&
         iNeiSouth != nothing && iNeiWest != nothing
         Connectivity[[1,3],i] .= -1
      end
      if iNeiNorthEast == nothing &&
         iNeiNorth != nothing && iNeiEast != nothing
         Connectivity[[2,4],i] .= -1
      end
      if iNeiNorthWest == nothing &&
         iNeiNorth != nothing && iNeiWest != nothing
         Connectivity[[1,4],i] .= -1
      end
      #
   end

   BCnode = Vector{Vector{Int}}(undef,nlines)
   for i=1:nlines BCnode[i] = [] end

   for i = 1:nNode
      if Connectivity[1,i] == -1
         Coord = [Coords[:,i],[-INF,Coords[2,i]]]
         for iBC = 1:nlines
            if BCDir[iBC]<=2 && doIntersect(
               BClines[1:2,iBC].+[0.0,dy],BClines[3:4,iBC].+[0.0,-dy],
               Coord[1],Coord[2])
               push!(BCnode[iBC],i)
               break
            end
         end
      end
      if Connectivity[2,i] == -1
         Coord = [Coords[:,i],[INF,Coords[2,i]]]
         for iBC = 1:nlines
            if BCDir[iBC]<=2 && doIntersect(
               BClines[1:2,iBC].+[0.0,-dy],BClines[3:4,iBC].+[0.0,dy],
               Coord[1],Coord[2])
               push!(BCnode[iBC],i)
               break
            end
         end
      end
      if Connectivity[3,i] == -1
         Coord = [Coords[:,i],[Coords[1,i],-INF]]
         for iBC = 1:nlines
            if BCDir[iBC] > 2 && doIntersect(
               BClines[1:2,iBC].+[-dx,0.0],BClines[3:4,iBC].+[dx,0.0],
               Coord[1],Coord[2])
               push!(BCnode[iBC],i)
               break
            end
         end
      end
      if Connectivity[4,i] == -1
         Coord = [Coords[:,i],[Coords[1,i],INF]]
         for iBC = nlines:-1:1 # Ugly trick
            if BCDir[iBC] > 2 && doIntersect(
               BClines[1:2,iBC].+[dx,0.0],BClines[3:4,iBC].+[-dx,0.0],
               Coord[1],Coord[2])
               push!(BCnode[iBC],i)
               break
            end
         end
      end
   end

   mesh = Mesh(Connectivity,BClines,BCDir,Coords,BCnode)

   return mesh
end

"""
Create triangle mesh for visualization.
"""
function createTriMesh(mesh::Mesh)
   points = transpose(mesh.Coords)
   # The functions in Triangle does not accept abstract array.
   points = copy(points)
   points_map = Array{Int64,1}(collect(1:1:size(mesh.Coords)[2]))
   # Find edge list: actually the size is known beforehand
   edges_list_1 = Array{Int64}(undef,0)
   edges_list_2 = Array{Int64}(undef,0)

   BCnode = copy(mesh.BCnode)
   BCDir = mesh.BCDir

   nBClines = length(mesh.BCnode)
   push!(BCnode,BCnode[1]) # For easy implementation of the final line

   for i = 1:nBClines
      if BCDir[i] == East_ || BCDir[i] == South_
         for j = 1:size(BCnode[i])[1]
            if j == size(BCnode[i])[1]
               push!(edges_list_1, BCnode[i][j])
               push!(edges_list_2, BCnode[i+1][end])
               continue
            end
            push!(edges_list_1, BCnode[i][j])
            push!(edges_list_2, BCnode[i][j+1])
         end
      else
         for j = size(BCnode[i])[1]:-1:1
            if j == 1
               push!(edges_list_1, BCnode[i][j])
               push!(edges_list_2, BCnode[i+1][end])
               continue
            end
            push!(edges_list_1, BCnode[i][j])
            push!(edges_list_2, BCnode[i][j-1])
         end
      end
   end

   edges_list = [edges_list_1 edges_list_2]

   edge_boundary = [true for i=1:size(edges_list)[1]]

   triMesh = Triangle.constrained_triangulation(points,points_map,edges_list,
   edge_boundary)

   return triMesh
end

end
