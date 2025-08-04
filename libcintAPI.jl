module libcintAPI

# Set LIBCINT to point to your conda-installed libcint library
const LIBCINT = joinpath(@__DIR__, "/opt/libcint/lib/libcint.dylib")  # Adjust for your system  

function CINTcgtos_spheric(bas_id, bas)
    @static if VERSION >= v"1.5.0" # clearer syntax
       @ccall LIBCINT.CINTcgtos_spheric(bas_id::Cint, bas::Ptr{Cint})::Cint
    else # works anyway
       ccall( (:CINTcgtos_spheric, LIBCINT), Cint, (Cint, Ptr{Cint},), bas_id, bas)
    end
end
 
function cint1e_ipnuc_sph(buf, shls, atm, natm, bas, nbas, env)
    @static if VERSION >= v"1.5.0" # clearer syntax
       @ccall LIBCINT.cint1e_ipnuc_sph(
                                       buf  :: Ptr{Cdouble},
                                       shls :: Ptr{Cint},
                                       atm  :: Ptr{Cint},
                                       natm :: Cint,
                                       bas  :: Ptr{Cint},
                                       nbas :: Cint,
                                       env  :: Ptr{Cdouble}
                                      )::Cvoid
    else # works anyway
       ccall( (:cint1e_ipnuc_sph, LIBCINT), Cvoid, (
                                       Ptr{Cdouble},
                                       Ptr{Cint},
                                       Ptr{Cint},
                                       Cint,
                                       Ptr{Cint},
                                       Cint,
                                       Ptr{Cdouble}, ),
                                       buf  ,
                                       shls ,
                                       atm  ,
                                       natm ,
                                       bas  ,
                                       nbas ,
                                       env  ,
                                      )
    end
end

function cint2e_ip_sph(buf, shls, atm, natm, bas, nbas, env)
   @static if VERSION >= v"1.5.0" # clearer syntax
      @ccall LIBCINT.cint2e_ip_sph(
                                      buf  :: Ptr{Cdouble},
                                      shls :: Ptr{Cint},
                                      atm  :: Ptr{Cint},
                                      natm :: Cint,
                                      bas  :: Ptr{Cint},
                                      nbas :: Cint,
                                      env  :: Ptr{Cdouble}
                                     )::Cvoid
   else # works anyway
      ccall( (:cint2e_ip_sph, LIBCINT), Cvoid, (
                                      Ptr{Cdouble},
                                      Ptr{Cint},
                                      Ptr{Cint},
                                      Cint,
                                      Ptr{Cint},
                                      Cint,
                                      Ptr{Cdouble}, ),
                                      buf  ,
                                      shls ,
                                      atm  ,
                                      natm ,
                                      bas  ,
                                      nbas ,
                                      env  ,
                                     )
   end
end

# Additional helper functions (e.g. gto_norm, data preparation) go here.

end # module
