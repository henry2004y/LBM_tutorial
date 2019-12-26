module BC

using ..lattice: LatticeD2Q4,LatticeD2Q5,LatticeD2Q9
using ..Parameters: Param, BCDirichlet, BCRiemann, Mesh

export setBCs!

function setBCs!(lat::LatticeD2Q4,param::Param,mesh::Mesh)
   Tleft,Tright,Ttop,Tbot = param.Tleft, param.Tright, param.Ttop, param.Tbot
   BCtype= param.BCtype
   BcNodeI1,BcNodeI2,BcNodeJ1,BcNodeJ2 = mesh.BcNodeI1,mesh.BcNodeI2,
      mesh.BcNodeJ1,mesh.BcNodeJ2
   w = lat.w

   if BCtype == BCDirichlet
      for i = eachindex(BcNodeI1) # it was 1:length(BcNodeI1)
         iL,iR = BcNodeI1[i],BcNodeI2[i]
         lat.f[1,iL] = (w[1]+w[2])*Tleft  .- lat.f[2,iL]
         lat.f[2,iR] = (w[1]+w[2])*Tright .- lat.f[1,iR]
      end

      for i = eachindex(BcNodeJ1)
         iB,iT = BcNodeJ1[i],BcNodeJ2[i]
         lat.f[3,iT] = (w[3]+w[4])*Ttop   .- lat.f[4,iT]
         lat.f[4,iB] = (w[3]+w[4])*Tbot   .- lat.f[3,iB]
      end

   elseif BCtype == BCRiemann

   end
end

function setBCs!(lat::LatticeD2Q5,param::Param)
   Tleft,Tright,Ttop,Tbot = param.Tleft, param.Tright, param.Ttop, param.Tbot
   nI,nJ = param.nI,param.nJ
   BCtype= param.BCtype
   w = lat.w

   if BCtype == BCDirichlet
      lat.f[1,1,:] = (w[1]+w[2])*Tleft  .- lat.f[2,1,:]
      lat.f[2,nI,:]= (w[1]+w[2])*Tright .- lat.f[1,nI,:]
      lat.f[3,:,1] = (w[3]+w[4])*Ttop   .- lat.f[4,:,1]
      lat.f[4,:,nJ]= (w[3]+w[4])*Tbot   .- lat.f[3,:,nJ]
   elseif BCtype == BCRiemann

   end
end


function setBCs!(lat::LatticeD2Q9,param::Param)
   Tleft,Tright,Ttop,Tbot = param.Tleft, param.Tright, param.Ttop, param.Tbot
   nI,nJ = param.nI,param.nJ
   BCtype= param.BCtype
   w = lat.w

   if BCtype == BCDirichlet
      lat.f[1,1,:] = (w[1]+w[3])*Tleft  .- lat.f[3,1,:]
      lat.f[5,1,:] = (w[5]+w[7])*Tleft  .- lat.f[7,1,:]
      lat.f[8,1,:] = (w[8]+w[6])*Tleft  .- lat.f[6,1,:]
      lat.f[3,nI,:]= (w[3]+w[1])*Tright .- lat.f[1,nI,:]
      lat.f[7,nI,:]= (w[7]+w[5])*Tright .- lat.f[5,nI,:]
      lat.f[6,nI,:]= (w[6]+w[8])*Tright .- lat.f[8,nI,:]
      lat.f[4,:,1] = (w[4]+w[2])*Ttop   .- lat.f[2,:,1]
      lat.f[8,:,1] = (w[8]+w[6])*Ttop   .- lat.f[6,:,1]
      lat.f[7,:,1] = (w[7]+w[5])*Ttop   .- lat.f[5,:,1]
      lat.f[2,:,nJ]= (w[2]+w[4])*Tbot   .- lat.f[4,:,nJ]
      lat.f[5,:,nJ]= (w[5]+w[8])*Tbot   .- lat.f[8,:,nJ]
      lat.f[6,:,nJ]= (w[6]+w[7])*Tbot   .- lat.f[7,:,nJ]
   elseif BCtype == BCRiemann

   end
end

end
