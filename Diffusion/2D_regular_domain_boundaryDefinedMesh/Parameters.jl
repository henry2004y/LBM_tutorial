module ModParameters

export Param, setParameters

const BCDirichlet,BCRiemann = 1,2 # Constants for different boundary types

using TOML

# Parameters
struct Param
   nD::Int             # Dimension of the system
   LatticeType::Int    # DmQn lattice model
   Δx::Float64         # Discretized space interval
   Δt::Float64         # Discretized time interval
   α ::Float64         # Thermal diffusivity
   ω ::Float64         # Relaxation coefficient
   TotalTime::Float64  # Total simulation time
   nT::Int             # Total timesteps

   Temp::Vector{Float64} # Left boundary temperature

   nI::Int
   nJ::Int
   BCtype::Array{Int,1}         # Type of boundary conditions
   BClines::Array{Float64,2}


   DoPlot::Bool

   function Param(BClines::Array{Float64,2},
      nD::Int=2,LatticeType::Int=3,dx::Float64=1.0,dt::Float64=1.0,
      α::Float64=0.25,TotalTime::Float64=400.0,
      Temp::Vector{Float64}=[1.0, 0.0, 0.0, 0.0],
      BCtype::Vector{Int}=[1,1,1,1], nI::Int=5, nJ::Int=5,
      DoPlot::Bool=false)
      @assert(0 ≤ TotalTime, "Simulation time must be positive!")

      new(nD,LatticeType,dx,dt,α,inv(nD*dt*α/dx^2 + 0.5),TotalTime,
         Int(TotalTime/dt),Temp,nI,nJ,BCtype,BClines,DoPlot)
   end
end

"""
   setParameters()

Read parameters from PARAM.toml and construct the parameter list.
"""
function setParameters()

   paramIn = TOML.parsefile("PARAM.toml")
   nD = paramIn["Parameters"]["nD"]
   LatticeType = paramIn["Parameters"]["LatticeType"]
   dx = paramIn["Parameters"]["dx"]
   dt = paramIn["Parameters"]["dt"]
   α = paramIn["Parameters"]["α"]
   TotalTime = paramIn["Parameters"]["TotalTime"]
   Temp = paramIn["Parameters"]["Temp"]
   DoPlot = paramIn["Plots"]["DoPlot"]

   if LatticeType == "D2Q4"
      LT = 1
   elseif LatticeType == "D2Q5"
      LT = 2
   elseif LatticeType == "D2Q9"
      LT = 3
   end

   nI = paramIn["Mesh"]["nI"]
   nJ = paramIn["Mesh"]["nJ"]
   nBCtype = length(paramIn["Mesh"]["BCtype"])
   if nBCtype != size(Temp,1) # This needs to be modified!
      error("Number of boundary values does not match! Check Temp in PARAM!")
   end
   BCtype = Vector{Int}(undef,nBCtype)
   for i = 1:nBCtype
      if paramIn["Mesh"]["BCtype"][i] == "Dirichlet"
         BCtype[i] = 1
      elseif paramIn["Mesh"]["BCtype"][i] == "Riemann"
         BCtype[i] = 2
      end
   end

   BClines = hcat(paramIn["Mesh"]["BClines"]...)

   param = Param(BClines,nD,LT,dx,dt,α,TotalTime,Temp,BCtype,nI,nJ,DoPlot)

   return param
end

end
