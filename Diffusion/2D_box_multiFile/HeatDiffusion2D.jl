# 2D Heat Diffusion
# LBM
#
# This could only work on a rectangular region! To generalize to more
# complicated mesh, I need 1D indexing of lattices.
#
# Hongyang Zhou, hyzhou@umich.edu 05/01/2019

using PyPlot
include("Parameters.jl")
include("lattice.jl")
include("BC.jl")

using .Parameters, .BC, .lattice

function solveHeatDiffusion2D!(param::Param,lat::Lattice)

   # Set aliases
   nI,nJ,ω,Δt,TotalTime = param.nI,param.nJ,param.ω,param.Δt,param.TotalTime

   nw,w  = lat.nw, lat.w
   T     = zeros(1,nI,nJ)
   feq   = Array{Float64,3}(undef,nw,nI,nJ)
   ftemp = Array{Float64,3}(undef,nw,nI,nJ)
   t     = 0.0     # Current simulation time

   # Advance
   while t < TotalTime

      # Calculate equilibrium distribution
      for iw = 1:nw
         @. feq[iw,:,:] = w[iw]*T[1,:,:]
      end

      # Collision
      for iw = 1:nw
         @. lat.f[iw,:,:] = (1.0-ω)*lat.f[iw,:,:] + ω.*feq[iw,:,:]
      end

      # Streaming
      for j=1:nJ, i=1:nI, iw = 1:nw
         iNew = i + lat.c[1,iw]
         jNew = j + lat.c[2,iw]
         (iNew > nI || iNew < 1 || jNew > nJ || jNew < 1) && continue
         ftemp[iw,iNew,jNew] = lat.f[iw,i,j]
      end
      #copy!(lat.f,ftemp) # This does not work in Julia 1.0 (but 1.1+)!
      lat.f .= ftemp

      # Remaining BCs
      setBCs!(lat,param)

      # Calculate macroscopic T
      T = sum(lat.f, dims=1)

      t += Δt
   end

end


function init_lattice()

   param = setParameters()

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


function plotSol(lat::Lattice)
   # Visualization
   x = range(0,stop=lat.nI-1,length=lat.nI)
   y = range(0,stop=lat.nJ-1,length=lat.nJ)
   T = sum(lat.f, dims=1)

   # Python follows row major order
   contourf(x,y,T[1,:,:]'); colorbar()
   title("2D Heat Diffusion in a Plate")
   xlabel("x")
   ylabel("y")
   grid(true)

   # 1D cut plot
   #plot(0:nx-1,T[1,:,50],linestyle="-",marker="o",label="LBM")
   #xlim(0,30)
end

function main()
   param,lat = init_lattice()

   solveHeatDiffusion2D!(param,lat)

   plotSol(lat)

end
