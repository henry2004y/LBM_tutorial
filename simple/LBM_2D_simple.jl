# 2D Lattice Boltzmann (BGK) model of a fluid.
#  c4  c3   c2  D2Q9 model. At each timestep, particle densities propagate
#    \  |  /    outwards in the directions indicated in the figure. An
#  c5 -c9 - c1  equivalent 'equilibrium' density is found, and the densities
#    /  |  \    relax towards that state, in a proportion governed by omega.
#  c6  c7   c8
#
# The first Julia version has many traces of Matlab, and it only works before
# version 0.6. The second version looks more like Julia.
#
# There is a lazy version of shifted array package. You can try and replace the
# circshift function with that.
#
# Orignal version by Iain Haslam, March 2006
# First Julia version by Jeroen Wouters, Jun 2015
# Second Julia version by Hongyang Zhou, Mar 2019

using LinearAlgebra

##################
# Configuration
##################
Ω        = 1.0    # Relaxation factor
ρ0       = 1.0    # Initial density
fraction = 0.1    # Average fraction of occupied nodes in the domain
rbody    = 5.0    # Radius of the solid cylinder

const nx,ny = 64,32
const nStep = 4000

ΔU = 1e-7

####################
# Setup of variables
####################
avu, avuOld = 1, 0
ts = 0
#avus = zeros(Float64, 10nStep)

# Single-particle distribution function
# Particle densities conditioned on one of the 9 possible velocities
F = fill(ρ0/9,nx,ny,9)
FEQ = F
nNode = nx*ny
# Offsets of different directions e_a in the F matrix
CI = [0:nNode:nNode*7...]

function setboundary()
   #BOUND = rand(nx,ny) .> 1 - fraction; # random domain
   center = [nx/2, ny/2]
   global BOUND = [norm([i, j]-center) < rbody for i=1:nx, j=1:ny]

   # Matrix offset of each Occupied Node (not physical domain)
   global ON = findall(x->x==true, vec(BOUND))

   global nActiveNodes = sum(1 .- BOUND)

   # Linear indices in F of occupied nodes
   global TO_REFLECT = [ON.+CI[1] ON.+CI[2] ON.+CI[3] ON.+CI[4] ON.+CI[5] ON.+CI[6] ON.+CI[7] ON.+CI[8]]

   # Right <-> Left: 1 <-> 5; Up <-> Down: 3 <-> 7
   #(1,1) <-> (-1,-1): 2 <-> 6; (1,-1) <-> (-1,1): 4 <-> 8
   global REFLECTED = [ON.+CI[5] ON.+CI[6] ON.+CI[7] ON.+CI[8] ON.+CI[1] ON.+CI[2] ON.+CI[3] ON.+CI[4]]
end

setboundary()

###############
# Constants
###############
const w1,w2,w3 = 4/9, 1/9, 1/36 # Weights
const c_squ = 1/3

Ux, Uy = Array{Float64,3}(undef,nx,ny,1), Array{Float64,3}(undef,nx,ny,1)

###############
# main loop
###############
while (ts < nStep) & (abs((avu-avuOld)/avu) > 1e-10)
   global F,avu,avuOld,ts,Ux,Uy # This is only needed outside function call
   ##############
   # Streaming
   ##############
   # particles at (x,y)=[2,1] were at [1,1] before: 1 points right (1,0)
   F[:,:,1] = F[circshift([1:nx...],1),:,1]
   # particles at [1,2] were at [1,1] before: 3 points up (0,1)
   F[:,:,3] = F[:,circshift([1:ny...],1),3]
   # particles at [1,1] were at [2,1] before: 5 points left (-1,0)
   F[:,:,5] = F[circshift([1:nx...],-1),:,5]
   # particles at [1,1] were at [1,2] before: 7 points down (0,-1)
   F[:,:,7] = F[:,circshift([1:ny...],-1),7]

   # particles at [2,2] were at [1, 1] before: 2 points to (1,1)
   F[:,:,2] = F[circshift([1:nx...],1), circshift([1:ny...],1), 2]
   # particles at [1,2] were at [2, 1] before: 4 points to (-1, 1)
   F[:,:,4] = F[circshift([1:nx...],-1),circshift([1:ny...],1), 4]
   # particles at [1,1] were at [2, 2] before: 6 points to (-1, -1)
   F[:,:,6] = F[circshift([1:nx...],-1),circshift([1:ny...],-1),6]
   # particles at [2,1] were at [1, 2] before: 8 points down (1,-1)
   F[:,:,8] = F[circshift([1:nx...],1), circshift([1:ny...],-1),8]

   ρ = sum(F,dims=3)

   # 1,2,8 are moving to the right, 4,5,6 to the left
   # 3, 7 and 9 don't move in the x direction
   Ux = (sum(F[:,:,[1, 2, 8]],dims=3) - sum(F[:,:,[4, 5, 6]],dims=3))./ρ
   # 2,3,4 are moving up, 6,7,8 down
   # 1, 5 and 9 don't move in the y direction
   Uy = (sum(F[:,:,[2, 3, 4]],dims=3) - sum(F[:,:,[6, 7, 8]],dims=3))./ρ

   Ux[1,1:ny] .+= ΔU # Increase inlet pressure

   Ux[ON] .= 0; Uy[ON] .= 0; ρ[ON] .= 0 # Set false cells to 0

   U2 = Ux.^2 + Uy.^2
   U_C2 = Ux + Uy
   U_C4 = -Ux + Uy

   # Calculate equilibrium distribution: stationary (a = 0)
   FEQ[:,:,9] = w1*ρ.*(1 .- U2/(2*c_squ))

   # nearest-neighbours
   FEQ[:,:,1] = w2*ρ.*(1 .+ Ux/c_squ + 0.5*(Ux/c_squ).^2 - 0.5*U2/c_squ)
   FEQ[:,:,3] = w2*ρ.*(1 .+ Uy/c_squ + 0.5*(Uy/c_squ).^2 - 0.5*U2/c_squ)
   FEQ[:,:,5] = w2*ρ.*(1 .- Ux/c_squ + 0.5*(Ux/c_squ).^2 - 0.5*U2/c_squ)
   FEQ[:,:,7] = w2*ρ.*(1 .- Uy/c_squ + 0.5*(Uy/c_squ).^2 - 0.5*U2/c_squ)

   # next-nearest neighbours
   FEQ[:,:,2] = w3*ρ.*(1 .+ U_C2/c_squ + 0.5*(U_C2/c_squ).^2 - 0.5*U2/c_squ)
   FEQ[:,:,4] = w3*ρ.*(1 .+ U_C4/c_squ + 0.5*(U_C4/c_squ).^2 - 0.5*U2/c_squ)
   FEQ[:,:,6] = w3*ρ.*(1 .- U_C2/c_squ + 0.5*(U_C2/c_squ).^2 - 0.5*U2/c_squ)
   FEQ[:,:,8] = w3*ρ.*(1 .- U_C4/c_squ + 0.5*(U_C4/c_squ).^2 - 0.5*U2/c_squ)

   BOUNCEDBACK = F[TO_REFLECT] # Densities bouncing back at next timestep

   F = Ω*FEQ + (1-Ω)*F

   F[REFLECTED] = BOUNCEDBACK

   avuOld = avu
   avu = sum(Ux)/nActiveNodes
   #avus[ts+1] = avu
   ts += 1
end

###############
# Visualization
###############
using PyPlot

figure()
imshow(1 .- BOUND', cmap="hot", interpolation="None", vmin=0., vmax=1., origin="lower");
quiver(1:nx-1, 0:ny-1, Ux[2:nx,:]', Uy[2:nx,:]');
title(string("Flow field after \$ ", string(ts), " \\delta t\$"));
xlabel("x")
ylabel("y")
