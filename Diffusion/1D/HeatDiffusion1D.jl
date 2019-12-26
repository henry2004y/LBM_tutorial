# 1D Heat Diffusion in an Infinite Slab
# LBM
#
# Hongyang Zhou, hyzhou@umich.edu 04/30/2019

# Problem Description
# A slab initially at temperature equal to zero, T=0. For time t>0, the left
# surface of the slab is subjected to a high temperature and equal to unity,
# T=1.0: The slab length is 100 units.
# ∂T/∂t = α ∂^2 T / ∂x^2
# Calculate the temperature distribution in the slab for t = 200, α=0.25.

# Solution
# Let us divide the domain of integration into Δx=1.0 and Δt=1.0.
# We are applying D1Q2 LBM with Dirichlet and Riemann BC.
# Note that in cylindrical and spherical coordinate systems, an axis-symmetric
# radial diffusion is similar to the diffusion in 1D Cartesian coordinates with
# an extra source term.

using PyPlot

function solveHeatDiffusion1D()
# Parameters
Δx = 1.0           # Discretized space interval
Δt = 1.0           # Discretized time interval
ω  = 4/3           # Relaxation coefficient
TotalTime = 200.   # Total simulation time
nT = 200           # Total timesteps
nx = 100           # Number of total space intervals
w  = [0.5,0.5]     # Lattice weights

Tleft  = 1.0       # Left boundary temperature
Tright = 0.0       # Right boundary temperature

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
f   = zeros(2,nx)
feq = zeros(2,nx)
T   = zeros(nx)
t   = 0.0     # Current simulation time

# Advance
while t < TotalTime

   # Calculate equilibrium distribution
   feq[1,:] .= w[1].*T[:]
   feq[2,:] .= w[2].*T[:]

   # Collision
   f[1,:] .= (1.0-ω).*f[1,:] .+ ω.*feq[1,:]
   f[2,:] .= (1.0-ω).*f[2,:] .+ ω.*feq[2,:]
   if UseSource
      f[1,:] .+= Δt*w[1]*source
      f[2,:] .+= Δt*w[2]*source
   end

   # Streaming
   f[1,2:nx] .= f[1,1:nx-1]
   f[2,1:nx-1] .= f[2,2:nx]

   # Remaining BCs
   if BCtype == BCDirichlet
      f[1,1] = Tleft - f[2,1]
      f[2,nx]= Tright- f[1,nx]
   elseif BCtype == BCRiemann
      f[1,1] = f[1,2] + f[2,2] - f[2,1] + heatflux*Δx/k_thermal
      f[2,nx]= Tright - f[1,nx]
   end

   # Calculate macroscopic T
   T .= f[1,:] .+ f[2,:]

   t += Δt
end


# Visualization
plot(0:nx-1,T,linestyle="-",marker="o",label="LBM")
xlim(0,30)
legend()
title("1D Heat Diffusion in an Infinite Slab")
xlabel("x")
ylabel("T")
grid(true)

end
