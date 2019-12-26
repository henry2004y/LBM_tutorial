# 1D advection-diffusion problem.
# LBM
#
# Hongyang Zhou, hyzhou@umich.edu 05/26/2019

# Problem Description
# A slab initially at temperature equal to zero, T=0. For time t>0, the left
# surface of the slab is subjected to a high temperature and equal to unity,
# T=1.0: The slab length is 100 units.
# ∂Φ/∂t + u ∂Φ/∂x = α ∂^2 Φ / ∂x^2
# Calculate the temperature distribution in the slab for
# t = 200, u = 0.1, α=0.25.

# Solution
# Let us divide the domain of integration into Δx=1.0 and Δt=1.0.
# We are applying D1Q2 LBM with Dirichlet and Riemann BC.

using PyPlot

function solveAdvectDiffuse1D()
# Parameters
Δx = 1.0           # Discretized space interval
Δt = 1.0           # Discretized time interval
ω  = 4/3           # Relaxation coefficient
TotalTime = 200.   # Total simulation time
nT = 200           # Total timesteps
nx = 100           # Number of total space intervals
w  = [0.5,0.5]     # Lattice weights
cs = 1/sqrt(2)     # sound speed

u  = 0.1           # constant advection speed
c  = [Δx/Δt,-Δx/Δt]# speed along the streaming direction
ϕᴸ = 1.0           # Left boundary value
ϕᴿ = 0.0           # Right boundary value

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
ϕ   = zeros(nx)
t   = 0.0     # Current simulation time

# Advance
while t < TotalTime

   # Calculate equilibrium distribution
   feq[1,:] .= w[1].*ϕ[:].*(1 + c[1]*u/cs^2)
   feq[2,:] .= w[2].*ϕ[:].*(1 + c[2]*u/cs^2)

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
      f[1,1] = ϕᴸ - f[2,1]
      f[2,nx]= ϕᴿ - f[1,nx]
   elseif BCtype == BCRiemann
      f[1,1] = f[1,2] + f[2,2] - f[2,1] + heatflux*Δx/k_thermal
      f[2,nx]= ϕᴿ - f[1,nx]
   end

   # Calculate macroscopic ϕ
   ϕ .= f[1,:] .+ f[2,:]

   t += Δt
end


# Visualization
plot(0:nx-1,ϕ,linestyle="-",marker="o",label="LBM")
xlim(0,100)
legend()
title("1D Advection-Diffusion Equation")
xlabel("x")
ylabel("T")
grid(true)

end
