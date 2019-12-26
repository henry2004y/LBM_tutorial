module BC

using ..lattice: LatticeD2Q4,LatticeD2Q5,LatticeD2Q9
using ..Parameters: Param,BCDirichlet,BCRiemann

export setBCs!

function setBCs!(lat::LatticeD2Q4,param::Param)
   T = param.T
   nI,nJ = param.nI,param.nJ
   BCtype= param.BCtype
   w = lat.w

   # Left
   if BCtype[1] == BCDirichlet
      lat.f[1,1,:] .= (w[1]+w[2])*T[1] .- lat.f[2,1,:]
   else
      lat.f[1,1,:] .= lat.f[1,2,:]
   end

   # Right
   if BCtype[2] == BCDirichlet
      lat.f[2,nI,:] .= (w[1]+w[2])*T[2] .- lat.f[1,nI,:]
   else
      lat.f[2,nI,:] .= lat.f[2,nI-1,:]
   end

   # Bottom
   if BCtype[3] == BCDirichlet
      lat.f[3,:,1] .= (w[3]+w[4])*T[3] .- lat.f[4,:,1]
   else
      lat.f[3,:,1] .= lat.f[3,:,2]
   end

   # Top
   if BCtype[4] == BCDirichlet
      lat.f[4,:,nJ] .= (w[3]+w[4])*T[4] .- lat.f[3,:,nJ]
   else
      lat.f[4,:,nJ] .= lat.f[2,nJ-1,:]
   end

end

function setBCs!(lat::LatticeD2Q5,param::Param)
   T = param.T
   nI,nJ = param.nI,param.nJ
   BCtype= param.BCtype
   w = lat.w

   # Left
   if BCtype[1] == BCDirichlet
      lat.f[1,1,:] .= (w[1]+w[2])*T[1] .- lat.f[2,1,:]
      lat.f[5,1,:] .= w[5]*T[1]
   else
      lat.f[1,1,:] .= lat.f[1,2,:]
   end

   # Right
   if BCtype[2] == BCDirichlet
      lat.f[2,nI,:] .= (w[1]+w[2])*T[2] .- lat.f[1,nI,:]
      lat.f[5,nI,:] .= w[5]*T[2]
   else
      lat.f[2,nI,:] .= lat.f[2,nI-1,:]
   end

   # Bottom
   if BCtype[3] == BCDirichlet
      lat.f[3,:,1] .= (w[3]+w[4])*T[3] .- lat.f[4,:,1]
      lat.f[5,:,1] .= w[5]*T[3]
   else
      lat.f[3,:,1] .= lat.f[3,:,2]
   end

   # Top
   if BCtype[4] == BCDirichlet
      lat.f[4,:,nJ] .= (w[3]+w[4])*T[4] .- lat.f[3,:,nJ]
      lat.f[5,:,nJ] .= w[5]*T[4]
   else
      lat.f[4,:,nJ] .= lat.f[4,:,nJ-1]
   end
end


function setBCs!(lat::LatticeD2Q9,param::Param)
   T = param.T
   nI,nJ = param.nI,param.nJ
   BCtype= param.BCtype
   w = lat.w

   if BCtype[1] == BCDirichlet
      lat.f[1,1,:] .= (w[1]+w[3])*T[1] .- lat.f[3,1,:]
      lat.f[5,1,:] .= (w[5]+w[7])*T[1] .- lat.f[7,1,:]
      lat.f[8,1,:] .= (w[8]+w[6])*T[1] .- lat.f[6,1,:]
      lat.f[9,1,:] .= w[9]*T[1]
   elseif BCtype[1] == BCRiemann
      lat.f[1,1,:] .= lat.f[8,2,:]
      lat.f[5,1,:] .= lat.f[1,2,:]
      lat.f[8,1,:] .= lat.f[5,2,:]
   end

   if BCtype[2] == BCDirichlet
      lat.f[3,nI,:] .= (w[3]+w[1])*T[2] .- lat.f[1,nI,:]
      lat.f[7,nI,:] .= (w[7]+w[5])*T[2] .- lat.f[5,nI,:]
      lat.f[6,nI,:] .= (w[6]+w[8])*T[2] .- lat.f[8,nI,:]
      lat.f[9,nI,:] .= w[9]*T[2]
   elseif BCtype[2] == BCRiemann
      lat.f[6,nI,:] .= lat.f[8,nI,:]
      lat.f[3,nI,:] .= lat.f[1,nI,:]
      lat.f[7,nI,:] .= lat.f[5,nI,:]
   end

   if BCtype[3] == BCDirichlet
      lat.f[2,:,1] .= (w[2]+w[4])*T[3] .- lat.f[4,:,1]
      lat.f[5,:,1] .= (w[5]+w[8])*T[3] .- lat.f[8,:,1]
      lat.f[6,:,1] .= (w[6]+w[7])*T[3] .- lat.f[7,:,1]
      lat.f[9,:,1] .= w[9]*T[3]
   elseif BCtype[3] == BCRiemann
      lat.f[2,:,1] .= lat.f[2,:,1]
      lat.f[5,:,1] .= lat.f[5,:,1]
      lat.f[6,:,1] .= lat.f[6,:,1]
   end

   if BCtype[4] == BCDirichlet
      lat.f[4,:,nJ] .= (w[4]+w[2])*T[4] .- lat.f[2,:,nJ]
      lat.f[8,:,nJ] .= (w[8]+w[6])*T[4] .- lat.f[6,:,nJ]
      lat.f[7,:,nJ] .= (w[7]+w[5])*T[4] .- lat.f[5,:,nJ]
      lat.f[9,:,nJ] .= w[9]*T[4]
   elseif BCtype[4] == BCRiemann
      lat.f[4,:,nJ] .= lat.f[4,:,nJ-1]
      lat.f[8,:,nJ] .= lat.f[8,:,nJ-1]
      lat.f[7,:,nJ] .= lat.f[7,:,nJ-1]
   end

end

end
