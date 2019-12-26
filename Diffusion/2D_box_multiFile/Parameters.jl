module Parameters

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
   nI::Int             # Number of space intervals in x
   nJ::Int             # Number of space intervals in y

   Tᴸ::Float64     # Left boundary temperature
   Tᴿ::Float64     # Right boundary temperature
   Tᵀ::Float64     # Top boundary temperature
   Tᴮ::Float64     # Bottom boundary temperature

   BCtype::Int         # Type of boundary conditions

   function Param(nD::Int=2,LatticeType::Int=3,dx::Float64=1.0,dt::Float64=1.0,
      α::Float64=0.25,TotalTime::Float64=400.0,nI::Int=100,nJ::Int=100,
      Tleft::Float64=1.0,Tright::Float64=0.0,
      Ttop::Float64=0.0,Tbot::Float64=0.0,
      BCtype::Int=1)
      @assert(0 ≤ TotalTime, "Simulation time must be positive!")

      new(nD,LatticeType,dx,dt,α,inv(nD*dt*α/dx^2 + 0.5),TotalTime,
         Int(TotalTime/dt),nI,nJ,Tleft,Tright,Ttop,Tbot,BCtype)
   end
end

function setParameters()
   paramIn = TOML.parsefile("PARAM.toml")
   nD = paramIn["Parameters"]["nD"]
   LatticeType = paramIn["Parameters"]["LatticeType"]
   dx = paramIn["Parameters"]["dx"]
   dt = paramIn["Parameters"]["dt"]
   α = paramIn["Parameters"]["α"]
   TotalTime = paramIn["Parameters"]["TotalTime"]
   nI = paramIn["Parameters"]["nI"]
   nJ = paramIn["Parameters"]["nJ"]
   Tleft = paramIn["Parameters"]["Tleft"]
   Tright = paramIn["Parameters"]["Tright"]
   Ttop = paramIn["Parameters"]["Ttop"]
   Tbot = paramIn["Parameters"]["Tbot"]
   BCtype = paramIn["Parameters"]["BCtype"]

   if LatticeType == "D2Q4"
      LT = 1
   elseif LatticeType == "D2Q5"
      LT = 2
   elseif LatticeType == "D2Q9"
      LT = 3
   end

   if BCtype == "Dirichlet"
      BC = 1
   elseif BCtype == "Riemann"
      BC = 2
   end

   param = Param(nD,LT,dx,dt,α,TotalTime,nI,nJ,Tleft,Tright,Ttop,Tbot,BC)

end

end
