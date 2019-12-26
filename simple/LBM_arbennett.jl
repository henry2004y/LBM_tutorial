#!/bin/julia
#
# An implementation of the lattice boltzmann method
#
# @author arbennett
# Collected by henry2004y as a reference.
##


# Constant parameters
const T   = 9000          # Number of timesteps to run
const Re  = 20.0           # Reynolds number
const nx  = 150              # Size of the domain in x direction
const ny  = 50             # Size of the domain in y direction
const q   = 9               # number of flow directions
const v   = 0.05            # velocity
const nu  = (v*ny)/(10*Re)  # viscosity
const tau = (3.0*nu+0.5)    # 1/relaxation
const plot_frequency = T/5

## Set up the flow grids, weights, and indices
#
#  Explain here
#
##
dirs =     [  0   1   0  -1   0   1   -1   -1    1 ;
0   0   1   0  -1   1    1   -1   -1 ]
weights =  [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]
indices =  [  4   5   6;  1   2   3;   7    8    9 ]
no_slip =  [  1   4   5   2   3   8    9    6    7 ]


## Create the obstacle array
#
#  The obst array is a boolean array that uses true to denote where
#  obstacles to flow are.  When the "fluid" encounters one of these
#  cells it will reverse the flow direction using the `no_slip`
#  condition detailed above
#
##
obst_coords = [nx/2, ny/2, 3]  # [x, y, radius]
obst_coords2 = [nx/3,ny/3,4]
y = collect(1:ny)
x = collect(1:nx)'
obst = (x.-obst_coords[1]).^2 .+ (y.-obst_coords[2]).^2 .<= obst_coords[3].^2
obst = (obst + ((x.-obst_coords2[1]).^2 .+ (y.-obst_coords2[2]).^2 .<= obst_coords2[3].^2) .> 0)
#obst[1,:] = true
#obst[end,:] = true
#obst[:,:] = false

##  Set up initial flow profile
#
#  Allocate the velocity array, u, which has 2 points for each grid
#  point, corresponding to the x and y velocity
#
#  Allocate the distribution array, f, which has 9 points for each grid
#  point, corresponding to the directions in dirs
##
u = zeros(ny,nx,2)
vel = zeros(ny,nx,2)
fEq = zeros(ny,nx,9)
rho = ones(ny,nx)

for j=1:ny, i=1:nx
   vel[j,i,2] = v*(1+0.1*sin((j/ny)/(2*pi)))
end

## Calculate equilibrium for some distribution
#
##
function equilib(u, rho)
   dist = zeros(ny,nx,9)
   for i = 1:9
      for x = 1:nx
         for y = 1:ny
            tmp = u[y,x,1] * dirs[1,i] + u[y,x,2] * dirs[2,i]
            dist[y,x,i] = rho[y,x] * weights[i] * (1 + 3*tmp + 4.5*tmp*tmp - 1.5*(u[y,x,1]*u[y,x,1] + u[y,x,2]*u[y,x,2]))
         end
      end
   end
   return dist
end
fIn = equilib(u, rho)
fLeft = fIn[:,1,:]

## Propogate time
#
##
for t=1:T
   if mod(t, plot_frequency) == 0
      println(t)
      writedlm("out/solution.dat."*lpad(string(t),length(string(int(T))),"0"),sqrt(u[:,:,1].*u[:,:,1] + u[:,:,2].*u[:,:,2]))
   end


   # Right wall conditions
   fIn[:,end,:] = fIn[:,end-1,:]
   rho = fIn[:,:,1] + fIn[:,:,2] + fIn[:,:,3] + fIn[:,:,4] + fIn[:,:,5] + fIn[:,:,6] + fIn[:,:,7] + fIn[:,:,8] + fIn[:,:,9]

   u[:,:,1] = (fIn[:,:,2] - fIn[:,:,4] + fIn[:,:,6] - fIn[:,:,7] - fIn[:,:,8] + fIn[:,:,9]) ./ rho[:,:]
   u[:,:,2] = (fIn[:,:,3] - fIn[:,:,5] + fIn[:,:,6] + fIn[:,:,7] - fIn[:,:,8] - fIn[:,:,9]) ./ rho[:,:]

   # Left wall conditions
   u[:,1,:] = vel[:,1,:]
   #rho[:,1] = 1./(1-u[:,1,2]) .* (fIn[:,1,3] + fIn[:,1,6] + fIn[:,1,7] + 2*(fIn[:,1,5] + fIn[:,1,8] + fIn[:,1,9]))
   #println(rho[:,1])

   fEq = fIn - (1/tau) * (fIn - equilib(u, rho))

   for i=1:q
      for x=1:nx
         for y=1:ny
            if obst[y,x]
               fEq[y,x,i] = fIn[y, x, no_slip[i]]
               u[y,x,:] = 0
            end
         end
      end
   end

   for i=1:q
      fIn[:,:,i] = circshift(fEq[:,:,i],dirs[:,i])
   end
end
println("done!")
