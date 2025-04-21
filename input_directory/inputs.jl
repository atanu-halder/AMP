module inputs
include("../source_directory/AMP.jl")
using .AMP
using DelimitedFiles

function define_rotor()
    Nb = 2                      # Number of blade
    R = 1.016                   # Rotor radius in m
    root_cut = 0.12             # Root-cut as % of radius
    N_elm = 100                 # Number of blade elements (<40 is desirable)
    dr = 1/N_elm
    r = dr/2:dr:1        
    c_bar = [ (r[i]>root_cut) ? 0.08/R : 0 for i in 1:N_elm]      
    stwist = 0*ones(N_elm)      # Sectional twist distrbution along span
    airfoil = readdlm("../data_airfoil/vr12.dat")

    # define rotor
    rotor = AMP.Rotor(Nb,R,root_cut,N_elm,r,c_bar,stwist,airfoil)
  
    return rotor
end

function define_op_params()

    mtip = 0.4435721219191       # tip Mach number
    theta0 = deg2rad(10)         # collective pitch
    thetac, thetas = 0.0, 0.0    # cyclic pitch and shaft tilt
    rot_flag = 1                 # rotational direction 1 -> CCW, -1 -> CW
    w = zeros(3)                 # angular speed of cg
    mu = zeros(3)                # non-dimensional speed of cg
    
    op_params = AMP.Op_params(mtip,theta0,thetac,thetas,mu,w,rot_flag)
end

end