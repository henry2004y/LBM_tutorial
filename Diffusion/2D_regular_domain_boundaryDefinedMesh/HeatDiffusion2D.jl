# 2D Heat Diffusion
# LBM
#
# Generalize version on complicated mesh, with 1D indexing of lattices.
# Boundary is defined by edges.
#
# There are plenty of room to improve code performance, if needed.
#
# Hongyang Zhou, hyzhou@umich.edu 05/01/2019

using PyPlot
include("Parameters.jl")
include("Mesh.jl")
include("lattice.jl")
include("BC.jl")

using .ModParameters, .ModMesh, .ModBC, .ModLattice

"""
   solveHeatDiffusion2D!(param,lat,mesh)

Kernel function for 2D LBM.
"""
function solveHeatDiffusion2D!(param::Param,lat::Lattice,mesh::Mesh)

# Set aliases
ω,Δt,TotalTime = param.ω,param.Δt,param.TotalTime
BCtype = param.BCtype
BCDirichlet,BCRiemann  = ModParameters.BCDirichlet, ModParameters.BCRiemann
Temp = param.Temp

nw,w,nN = lat.nw, lat.w, lat.nN
T       = zeros(1,nN)
feq     = Array{Float64,2}(undef,nw,nN)
ftemp   = Array{Float64,2}(undef,nw,nN)
t       = 0.0     # Current simulation time

# Advance
while t < TotalTime

   # Calculate equilibrium distribution
   for iw = 1:nw
      @. feq[iw,:] = w[iw]*T[1,:]
   end

   # Collision
   for iw = 1:nw
      @. lat.f[iw,:] = (1.0-ω)*lat.f[iw,:] + ω.*feq[iw,:]
   end

   # Streaming
   for i=1:nN, iw = 1:nw
      iN = mesh.Connectivity[iw,i]
      (iN > nN || iN < 1) && continue
      ftemp[iw,iN] = lat.f[iw,i]
   end
   #copy!(lat.f,ftemp) # This does not work in Julia 1.0 (but 1.1+)!
   lat.f .= ftemp

   # Remaining BCs
   setBCs!(lat,param,mesh)

   # Calculate macroscopic T
   T .= sum(lat.f, dims=1)

   t += Δt
end

end

"""
   init()

Set parameters, create mesh and lattice structure.
"""
function init()

param = setParameters()
mesh  = createMesh(param.BClines, param.nI, param.nJ, param.LatticeType)

Δx,Δt,LT = param.Δx,param.Δt,param.LatticeType

if LT == 1
   lat = LatticeD2Q4(Δx,Δt,mesh)
elseif LT == 2
   lat = LatticeD2Q5(Δx,Δt,mesh)
elseif LT == 3
   lat = LatticeD2Q9(Δx,Δt,mesh)
end

return param,lat,mesh

end

"""
   plotSol(lat,mesh)

Visualize the solution.
"""
function plotSol(lat::Lattice,mesh::Mesh)
# Visualization
#x = range(0,stop=lat.nI-1,length=lat.nI)
#y = range(0,stop=lat.nJ-1,length=lat.nJ)
T = sum(lat.f, dims=1)

triangles = createTriMesh(mesh)
triangles = transpose(hcat(triangles...))
triangles .-= 1

# This must be fixed for constrained delaunay triangulation!
# Currently there are issues for using triangles: when doing interpolation, it
# destroys the symmetric solution. This is kind of unavoidable, but it is only a
# visual effect.
tricontourf(mesh.Coords[1,:],mesh.Coords[2,:],triangles,T[1,:]); colorbar()
title("2D Heat Diffusion in a Plate")
xlabel("x")
ylabel("y")
grid(true)

# 1D cut plot
#plot(0:nx-1,T[1,:,50],linestyle="-",marker="o",label="LBM")
#xlim(0,30)
end

"""
   main()

Main control function for the 2D diffusion problem.
"""
function main()
   param,lat,mesh = init()

   solveHeatDiffusion2D!(param,lat,mesh)

   if param.DoPlot
      plotSol(lat,mesh)
      #gcf() # This is needed for Windows using Juno
   end

end
