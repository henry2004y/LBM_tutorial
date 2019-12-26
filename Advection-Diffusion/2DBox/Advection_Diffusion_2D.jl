# 2D Advection-Diffusion problem
# LBM
#
# Hongyang Zhou, hyzhou@umich.edu 05/27/2019

# Let us consider flow in a square porous domain.
# Initially the domain is at zero temperature. The upper boundary of the domain
# is kept at zero temperature, while the bottom and vertical right boundaries
# are assumed to be adiabatic. The left vertical wall is subjected to unity
# temperature at time = 0. Find the temperature distribution in the domain at
# time = 200, for material of diffusivity of 1.0. The velocity field of u = 0.1
# and v = 0.2 is kept constant in the domain.

using PyPlot
include("Parameters.jl")
include("lattice.jl")
include("BC.jl")

using .Parameters, .BC, .lattice
using LinearAlgebra

function solveAdvectDiffuse2D!(param::Param,lat::Lattice)

   # Set aliases
   nI,nJ,ω,Δt,TotalTime = param.nI,param.nJ,param.ω,param.Δt,param.TotalTime
   u = param.u

   nw,w,cssq,c = lat.nw, lat.w, lat.cssq, lat.c
   T     = zeros(1,nI,nJ)
   feq   = Array{Float64,3}(undef,nw,nI,nJ)
   ftemp = Array{Float64,3}(undef,nw,nI,nJ)
   t     = 0.0     # Current simulation time

   # Advance
   while t < TotalTime

      # Calculate equilibrium distribution
      for iw = 1:nw
         feq[iw,:,:] .= w[iw].*T[1,:,:].*(1 + dot(c[:,iw],u)/cssq)
      end

      # Collision
      for iw = 1:nw
         lat.f[iw,:,:] .= (1.0-ω).*lat.f[iw,:,:] .+ ω.*feq[iw,:,:]
      end

      # Streaming
      for j=1:nJ, i=1:nI, iw = 1:nw
         iNew = i + c[1,iw]
         jNew = j + c[2,iw]
         (iNew > nI || iNew < 1 || jNew > nJ || jNew < 1) && continue
         ftemp[iw,iNew,jNew] = lat.f[iw,i,j]
      end
      copy!(lat.f,ftemp)

      # Remaining BCs
      setBCs!(lat,param)

      # Calculate macroscopic T
      T .= sum(lat.f, dims=1)

      t += Δt
   end

end


function init_lattice()

   param = setParameters()

   if param.DoEcho printParam(param) end

   Δx,Δt,nI,nJ,LT = param.Δx,param.Δt,param.nI,param.nJ,param.LatticeType

   if LT == 1
      lat = LatticeD2Q4(Δx,Δt,nI,nJ)
   elseif LT == 2
      lat = LatticeD2Q5(Δx,Δt,nI,nJ)
   elseif LT == 3
      lat = LatticeD2Q9(Δx,Δt,nI,nJ)
   end

   return param,lat

end

function printParam(param::Param)
   println("Input parameters for LBM advection-diffusion:")
   println("Lattice type: ",param.LatticeType)
   println("nI, nJ = $(param.nI), $(param.nJ)")
   println("BC type = $(param.BCtype)")
   # may be more to display?
end

function plotSol(lat::Lattice)
   # Visualization
   x = range(0,stop=lat.nI-1,length=lat.nI)
   y = range(0,stop=lat.nJ-1,length=lat.nJ)
   T = sum(lat.f, dims=1)

   # Python follows row major order
   contourf(x,y,T[1,:,:]'); colorbar()
   title("2D Heat Advection-Diffusion in a Plain")
   xlabel("x")
   ylabel("y")
   grid(true)

   # 1D cut plot
   plot(0:lat.nI-1,T[1,:,floor(Int,lat.nJ/2)],
      linestyle="-",marker="o",label="LBM")
end

function main()
   param,lat = init_lattice()

   solveAdvectDiffuse2D!(param,lat)

   plotSol(lat)

end
