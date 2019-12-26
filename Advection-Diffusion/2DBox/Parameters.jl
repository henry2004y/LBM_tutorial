module Parameters

export Param, setParameters

const BCDirichlet,BCRiemann = 1,2 # Constants for different boundary types

using TOML

# Parameters
struct Param
   nD::Int             # Dimension of the system
   LatticeType::Int    # DmQn lattice model
   Δx::Float64         # Discretized space interval
   Δy::Float64         # Discretized space interval (Δy = Δx)
   Δt::Float64         # Discretized time interval
   u ::Vector{Float64} # Advection velocity
   α ::Float64         # Thermal diffusivity
   ω ::Float64         # Relaxation coefficient
   TotalTime::Float64  # Total simulation time
   nT::Int             # Total timesteps
   nI::Int             # Number of space intervals in x
   nJ::Int             # Number of space intervals in y
   T::Vector{Float64}  #  boundary temperatures
   BCtype::Array{Int,1}# Type of boundary conditions
   DoEcho::Bool        # Display parameters on screen

   function Param(nD::Int=2,LatticeType::Int=3,dx::Float64=1.0,dy::Float64=1.0,
      dt::Float64=1.0,
      u::Vector{Float64}=[0.1,0.0],α::Float64=0.25,TotalTime::Float64=400.0,
      nI::Int=100,nJ::Int=100,
      Tᴸ::Float64=1.0,Tᴿ::Float64=0.0,
      Tᴮ::Float64=0.0,Tᵀ::Float64=0.0,
      BCtype::Vector{Int}=[1,1,1,1],
      DoEcho::Bool=true)
      @assert(0 ≤ TotalTime, "Simulation time must be positive!")

      new(nD,LatticeType,dx,dy,dt,u,α,inv(nD*dt*α/dx^2 + 0.5),
         TotalTime,Int(TotalTime/dt),nI,nJ,[Tᴸ,Tᴿ,Tᴮ,Tᵀ],BCtype,DoEcho)
   end
end

function setParameters()
   paramIn = TOML.parsefile("PARAM.toml")
   DoEcho = paramIn["Parameters"]["DoEcho"]
   nD = paramIn["Parameters"]["nD"]
   LatticeType = paramIn["Parameters"]["LatticeType"]
   dx = paramIn["Parameters"]["dx"]
   dy = paramIn["Parameters"]["dy"]
   dt = paramIn["Parameters"]["dt"]
   u = paramIn["Parameters"]["u"]
   α = paramIn["Parameters"]["α"]
   TotalTime = paramIn["Parameters"]["TotalTime"]
   nI = paramIn["Parameters"]["nI"]
   nJ = paramIn["Parameters"]["nJ"]
   Tᴸ = paramIn["Parameters"]["Tleft"]
   Tᴿ = paramIn["Parameters"]["Tright"]
   Tᴮ = paramIn["Parameters"]["Tbot"]
   Tᵀ = paramIn["Parameters"]["Ttop"]

   if LatticeType == "D2Q4"
      LT = 1
   elseif LatticeType == "D2Q5"
      LT = 2
   elseif LatticeType == "D2Q9"
      LT = 3
   end

   nBC = length(paramIn["Parameters"]["BCtype"])
   BCtype = Vector{Int}(undef,nBC)
   for i = 1:nBC
       if paramIn["Parameters"]["BCtype"][i] == "Dirichlet"
           BCtype[i] = 1
       elseif paramIn["Parameters"]["BCtype"][i] == "Riemann"
           BCtype[i] = 2
       end
   end

   param = Param(nD,LT,dx,dy,dt,u,α,TotalTime,nI,nJ,Tᴸ,Tᴿ,Tᴮ,Tᵀ,BCtype,DoEcho)

end

end
