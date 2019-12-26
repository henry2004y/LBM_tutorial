module ModBC

using ..ModLattice: LatticeD2Q4,LatticeD2Q5,LatticeD2Q9, opp_lat_vec
using ..ModParameters: Param, BCDirichlet, BCRiemann
using ..ModMesh: Mesh, East_,West_,North_,South_

export setBCs!

"""
    setBCs!(lat,param,mesh)

Set the boundary conditions. Note that BCDir: right 1, left 2, top 3, bottom 4
This is consistent with the lattice direction, but different from the
connectivity orientation.
"""
function setBCs!(lat::LatticeD2Q4,param::Param,mesh::Mesh)
   Temp = param.Temp
   BCtype = param.BCtype
   BCDir  = mesh.BCDir
   BCnode = mesh.BCnode
   Conn   = mesh.Connectivity
   w,nw = lat.w, lat.nw

   for iBC = 1:length(BCtype)
      if BCtype[iBC] == BCDirichlet
         i = BCDir[iBC]
         j = opp_lat_vec(lat,i)
         lat.f[j,BCnode[iBC]] = (w[i]+w[j])*Temp[iBC] .- lat.f[i,BCnode[iBC]]

      elseif BCtype[iBC] == BCRiemann
         i = BCDir[iBC]
         for iBCnode in BCnode[iBC]
            lat.f[:,iBCnode] = lat.f[:,Conn[i,iBCnode]]
         end
      end
   end

end

function setBCs!(lat::LatticeD2Q5,param::Param)
   T = param.T
   nI,nJ = param.nI,param.nJ
   BCtype= param.BCtype
   w = lat.w

   if BCtype == BCDirichlet

   elseif BCtype == BCRiemann

   end
end

function setBCs!(lat::LatticeD2Q9,param::Param)
   Tleft,Tright,Ttop,Tbot = param.Tleft, param.Tright, param.Ttop, param.Tbot
   nI,nJ = param.nI,param.nJ
   BCtype= param.BCtype
   w = lat.w

   if BCtype == BCDirichlet

   elseif BCtype == BCRiemann

   end
end

end
