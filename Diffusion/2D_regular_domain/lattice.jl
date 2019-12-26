module lattice

export Lattice,LatticeD2Q9,LatticeD2Q4,LatticeD2Q5

using ..Parameters: Mesh

abstract type Lattice end

# D2Q9 lattice
struct LatticeD2Q9 <: Lattice
  Δx::AbstractFloat
  Δt::AbstractFloat
  cc::AbstractFloat
  cs::AbstractFloat
  cs²::AbstractFloat
  f::Array{Float64, 3}
  c::Matrix{Int64}
  w::Vector{Float64}
  nw::Int

  function LatticeD2Q9(dx::AbstractFloat, dt::AbstractFloat, nI::Int, nJ::Int,
                       rho::AbstractFloat=0.0)
    # Default lattice speed vectors and associated weights
    cdef = [1 0; 0 1; -1 0; 0 -1; 1 1; -1 1; -1 -1; 1 -1; 0 0]'
    wdef = [1.0/9.0; 1.0/9.0; 1.0/9.0; 1.0/9.0; 1.0/36.0; 1.0/36.0;
            1.0/36.0; 1.0/36.0; 4.0/9.0]
    nw = length(wdef)
    # Lattice functions for computing `c`, `cs`, and `cs²`
    cf(dx, dt)    = dx/dt
    csf(dx, dt)   = cf(dx, dt) / sqrt(3)
    cssqf(dx, dt) = csf(dx, dt) * csf(dx, dt)
    f = zeros(Float64, (nw, nI, nJ))
    for k=1:nw
      f[k,:,:] = fill(rho * wdef[k], (nI, nJ))
    end
    new(dx, dt, cf(dx, dt), csf(dx, dt), cssqf(dx, dt), f, cdef, wdef,
      nw)
  end

end

# D2Q4 lattice
struct LatticeD2Q4 <: Lattice
  Δx::AbstractFloat     # Spatial interval
  Δt::AbstractFloat     # Time interval
  cc::AbstractFloat     # Lattice speed
  cs::AbstractFloat
  cs²::AbstractFloat
  f::Array{Float64, 2}  # Distribution
  IJ::Array{Int,2}
  c::Matrix{Int64}      # Speed vectors
  w::Vector{Float64}    # Weights
  nw::Int               # Number of velocity directions
  nN::Int               # Number of lattice nodes

  function LatticeD2Q4(dx::AbstractFloat, dt::AbstractFloat,
                       mesh::Mesh, rho::Float64=0.0)
    # Default lattice speed vectors and associated weights
    cdef = [1 0; -1 0; 0 1; 0 -1]'
    wdef = [1.0/4.0; 1.0/4.0; 1.0/4.0; 1.0/4.0]
    nw = size(cdef, 2)

    # Lattice functions for computing `c`, `cs`, and `cs²`
    cf(dx, dt)    = dx/dt
    csf(dx, dt)   = cf(dx, dt) / sqrt(2)
    cssqf(dx, dt) = csf(dx, dt) * csf(dx, dt)

    nN = size(mesh.connectivity,2)

    f = zeros(Float64, (nw, nN))
    for k=1:nw
      f[k,:] = fill(rho * wdef[k], nN)
    end

    new(dx, dt, cf(dx, dt), csf(dx, dt), cssqf(dx, dt), f, mesh.connectivity,
      cdef, wdef,nw,nN)
  end

#  LatticeD2Q4(dx::AbstractFloat, dt::AbstractFloat, f::Array{Float64, 3}) =
#    new(dx, dt, cf(dx, dt), csf(dx, dt), cssqf(dx, dt), f, cdef, wdef, n)
end

# D2Q4 lattice
struct LatticeD2Q5 <: Lattice
  Δx::AbstractFloat     # Spatial interval
  Δt::AbstractFloat     # Time interval
  cc::AbstractFloat     # Lattice speed
  cs::AbstractFloat
  cs²::AbstractFloat
  f::Array{Float64, 3}  # Distribution
  c::Matrix{Int64}      # Speed vectors
  w::Vector{Float64}    # Weights
  nw::Int               # Number of velocity directions
  nI::Int               # This should not belong to lattice?
  nJ::Int

  function LatticeD2Q5(dx::AbstractFloat, dt::AbstractFloat, nI::Int, nJ::Int,
                       rho::Float64=0.0)
    # Default lattice speed vectors and associated weights
    cdef = [1 0; -1 0; 0 1; 0 -1; 0 0]'
    wdef = [2.0/12.0; 2.0/12.0; 1.0/12.0; 1.0/12.0; 1.0/12.0]
    nw = size(cdef, 2)

    # Lattice functions for computing `c`, `cs`, and `cs²`
    cf(dx, dt)    = dx/dt
    csf(dx, dt)   = cf(dx, dt) / sqrt(3)
    cssqf(dx, dt) = csf(dx, dt) * csf(dx, dt)

    f = zeros(Float64, (nw, nI, nJ))
    for k=1:nw
      f[k,:,:] = fill(rho * wdef[k], (nI, nJ))
    end
    new(dx, dt, cf(dx, dt), csf(dx, dt), cssqf(dx, dt), f, cdef, wdef,
      nw, nI, nJ)
  end

#  LatticeD2Q4(dx::AbstractFloat, dt::AbstractFloat, f::Array{Float64, 3}) =
#    new(dx, dt, cf(dx, dt), csf(dx, dt), cssqf(dx, dt), f, cdef, wdef, n)
end

# Opposite direction lattice vector index
#
# \param lat Lattice
# \param k Index of lattice vector
# \return Index of opposite direction
function opp_lat_vec(lat::LatticeD2Q9, k::Int)
  if k == 1
    return 3
  elseif k == 2
    return 4
  elseif k == 3
    return 1
  elseif k == 4
    return 2
  elseif k == 5
    return 7
  elseif k == 6
    return 8
  elseif k == 7
    return 5
  elseif k == 8
    return 6
  elseif k == 9
    return 9
  else
    error("$k > 9 for a D2Q9 lattice. Only nine vectors (1-9) possible");
  end
end

#=
# Opposite direction lattice vector index
#
# \param lat Lattice
# \param k Index of lattice vector
# \return Index of opposite direction
function opp_lat_vec(lat::LatticeD2Q4, k::Int)
  error("not yet implemented");
end

# Helper function for initializing lattice
function _fill_lat(lat::Lattice, i_range::UnitRange{Int}, j_range::UnitRange{Int},
                   rho::Real)
  for k=1:lat.n, j=j_range, i=i_range
    lat.f[k, i, j] = rho * lat.w[k];
  end
end
=#

end
