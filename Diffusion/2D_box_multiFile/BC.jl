module BC

using ..lattice: LatticeD2Q4,LatticeD2Q5,LatticeD2Q9
using ..Parameters: Param,BCDirichlet,BCRiemann

export setBCs!

function setBCs!(lat::LatticeD2Q4,param::Param)
   Tᴸ,Tᴿ,Tᵀ,Tᴮ = param.Tᴸ, param.Tᴿ, param.Tᵀ, param.Tᴮ
   nI,nJ = param.nI,param.nJ
   BCtype= param.BCtype
   w = lat.w

   if BCtype == BCDirichlet
      lat.f[1,1,:] .= (w[1]+w[2])*Tᴸ .- lat.f[2,1,:]
      lat.f[2,nI,:].= (w[1]+w[2])*Tᴿ .- lat.f[1,nI,:]
      lat.f[3,:,1] .= (w[3]+w[4])*Tᵀ .- lat.f[4,:,1]
      lat.f[4,:,nJ].= (w[3]+w[4])*Tᴮ .- lat.f[3,:,nJ]
   elseif BCtype == BCRiemann

   end
end

function setBCs!(lat::LatticeD2Q5,param::Param)
   Tᴸ,Tᴿ,Tᵀ,Tᴮ = param.Tᴸ, param.Tᴿ, param.Tᵀ, param.Tᴮ
   nI,nJ = param.nI,param.nJ
   BCtype= param.BCtype
   w = lat.w

   if BCtype == BCDirichlet
      lat.f[1,1,:] .= (w[1]+w[2])*Tᴸ .- lat.f[2,1,:]
      lat.f[5,1,:] .= w[5]*Tᴸ
      lat.f[2,nI,:].= (w[1]+w[2])*Tᴿ .- lat.f[1,nI,:]
      lat.f[5,nI,:].= w[5]*Tᴿ
      lat.f[3,:,1] .= (w[3]+w[4])*Tᵀ .- lat.f[4,:,1]
      lat.f[5,:,1] .= w[5]*Tᵀ
      lat.f[4,:,nJ].= (w[3]+w[4])*Tᴮ .- lat.f[3,:,nJ]
      lat.f[5,:,nJ].= w[5]*Tᴮ
   elseif BCtype == BCRiemann
      #lat.f[1,1,:] .= lat.f[1,2,:]
      #...
   end
end


function setBCs!(lat::LatticeD2Q9,param::Param)
   Tᴸ,Tᴿ,Tᵀ,Tᴮ = param.Tᴸ, param.Tᴿ, param.Tᵀ, param.Tᴮ
   nI,nJ = param.nI,param.nJ
   BCtype= param.BCtype
   w = lat.w

   if BCtype == BCDirichlet
      lat.f[1,1,:] .= (w[1]+w[3])*Tᴸ .- lat.f[3,1,:]
      lat.f[5,1,:] .= (w[5]+w[7])*Tᴸ .- lat.f[7,1,:]
      lat.f[8,1,:] .= (w[8]+w[6])*Tᴸ .- lat.f[6,1,:]
      lat.f[9,1,:] .= w[9]*Tᴸ
      lat.f[3,nI,:].= (w[3]+w[1])*Tᴿ .- lat.f[1,nI,:]
      lat.f[7,nI,:].= (w[7]+w[5])*Tᴿ .- lat.f[5,nI,:]
      lat.f[6,nI,:].= (w[6]+w[8])*Tᴿ .- lat.f[8,nI,:]
      lat.f[9,nI,:].= w[9]*Tᴿ
      lat.f[4,:,1] .= (w[4]+w[2])*Tᵀ .- lat.f[2,:,1]
      lat.f[8,:,1] .= (w[8]+w[6])*Tᵀ .- lat.f[6,:,1]
      lat.f[7,:,1] .= (w[7]+w[5])*Tᵀ .- lat.f[5,:,1]
      lat.f[9,:,1] .= w[9]*Tᵀ
      lat.f[2,:,nJ].= (w[2]+w[4])*Tᴮ .- lat.f[4,:,nJ]
      lat.f[5,:,nJ].= (w[5]+w[8])*Tᴮ .- lat.f[8,:,nJ]
      lat.f[6,:,nJ].= (w[6]+w[7])*Tᴮ .- lat.f[7,:,nJ]
      lat.f[9,:,nJ].= w[9]*Tᴮ
   elseif BCtype == BCRiemann

   end
end

end
