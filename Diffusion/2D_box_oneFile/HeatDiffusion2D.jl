# 2D Heat Diffusion
# LBM
#
# Hongyang Zhou, hyzhou@umich.edu 05/01/2019

# Problem Description
# A two dimensional square slab is subjected to the given boundary conditions.
# Initially the slab was at zero temperature. For time>0, the boundary at x = 0
# is subjected to a high temperature of value 1.0 and other boundaries are kept
# as before. The length of the domain is 100 units. Determine the temperature
# distribution in the slab at time = 400 units (s).
# Thermal diffusivity is 0.25 and xl = yl = 1.0.

using PyPlot

function solveHeatDiffusion2D()
# Parameters
nD = 2             # Dimension of the system
Δx = 1.0           # Discretized space interval
Δt = 1.0           # Discretized time interval
α  = 0.25          # Thermal diffusivity
ω  = inv(nD*Δt*α/Δx^2 + 0.5)  # Relaxation coefficient
TotalTime = 400.   # Total simulation time
nT = Int(TotalTime/Δt)        # Total timesteps
nx = 100           # Number of total space intervals
ny = 100
w  = [0.25,0.25,0.25,0.25]    # Lattice weights
nw = size(w,1)

Tleft  = 1.0       # Left boundary temperature
Tright = 0.0       # Right boundary temperature
Ttop   = 0.0       # Top boundary temperature
Tbot   = 0.0       # Bottom boundary temperature

heatflux  = 100.   # Constant heat flux, [W/mK]
k_thermal = 20.    # Thermal conductivity, [W/mK]

BCDirichlet,BCRiemann = 1,2 # Constants for different boundary types
BCtype = 1         # Type of boundary conditions

UseSource = false  # Logical for calculating source terms
if UseSource
   qg = 1.0        # Rate of heat generation per volume
   ρ  = 1.0        # Density
   C  = 1.0        # Heat capacity
   source = qg/(ρ*C) # Source term. Note that α=k/(ρC)
end

# Initializations
f   = zeros(nw,nx,ny)
feq = zeros(nw,nx,ny) # This can be moved inside the main loop!
T   = zeros(1,nx,ny)
t   = 0.0     # Current simulation time

# Advance
while t < TotalTime

   # Calculate equilibrium distribution
   for iw = 1:nw
      feq[iw,:,:] .= w[iw].*T[1,:,:]
   end

   # Collision
   for iw = 1:nw
      f[iw,:,:] .= (1.0-ω).*f[iw,:,:] .+ ω.*feq[iw,:,:]
   end
   if UseSource
      for iw = 1:nw
         f[iw,:,:] .+= Δt*w[iw]*source
      end
   end

   # Streaming
   f[1,2:nx,:]   .= f[1,1:nx-1,:]
   f[2,1:nx-1,:] .= f[2,2:nx,:]
   f[3,:,2:ny]   .= f[3,:,1:ny-1]
   f[4,:,1:ny-1] .= f[4,:,2:ny]

   # Remaining BCs
   if BCtype == BCDirichlet
      f[1,1,:] .= (w[1]+w[2])*Tleft  .- f[2,1,:]
      f[2,nx,:].= (w[1]+w[2])*Tright .- f[1,nx,:]
      f[3,:,1] .= (w[3]+w[4])*Ttop   .- f[4,:,1]
      f[4,:,ny].= (w[3]+w[4])*Tbot   .- f[3,:,ny]
   elseif BCtype == BCRiemann
      # This is not complete, and not flexible enough!
      f[1,1,:] .= f[1,2,:] .- heatflux*Δx/k_thermal
      #f[2,nx,:] = f[2,nx-1,:] .- heatflux*Δx/k_thermal
   end

   # Calculate macroscopic T
   T .= sum(f, dims=1)

   t += Δt
end


# Visualization
x = range(0,stop=nx-1,length=nx)
y = range(0,stop=ny-1,length=ny)

contourf(x,y,T[1,:,:]); colorbar()
title("2D Heat Diffusion in a Plate")
xlabel("x")
ylabel("y")
grid("on")

# 1D cut plot
#plot(0:nx-1,T[1,:,50],linestyle="-",marker="o",label="LBM")
#xlim(0,30)

end
