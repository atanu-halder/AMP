## AMP stands for Aerodynamic Modeling of Propulsor
# Lead Developer: Dr. Atanu Halder
# Co-developers: Pedram H Dabaghian

module AMP
using Interpolations
using Statistics

struct Rotor
    Nb          :: Int              # Number of blades
    R           :: Real             # Rotor radius in meter
    root_cut    :: Real             # Root-cut as % of radius
    N_elm       :: Int              # Number of Blade elements
    r           :: Vector{Real}     # radial location of each element's center points
    c_bar       :: Vector{Real}     # sectional distribution of chord/Radius
    stwist      :: Vector{Real}     # sectional twist distribution (in radian)
    airfoil     :: Matrix{Real}     # airfoil look-up table
end

mutable struct Op_params    
    mtip        :: Real             # Tip Mach number
    theta0      :: Real             # Collective pitch      
    thetac      :: Real             # cyclic pitch - cosinusoidal
    thetas      :: Real             # cyclic pitch - sinusoidal
    mu          :: Vector{Real}     # body frame trnaslational velocity components (non-dimensionalized by tip-speed)
    w           :: Vector{Real}     # body frame angular velocity
    rot_flag    :: Int              # flag for rotational direc
end

function bemt_axial(rotor,op_params)    
    itermax = 100
    error_cutoff = 1e-5
    beta = 0.1
    tiploss = 1
       
    dr = 1/rotor.N_elm  
    AR = (1-rotor.root_cut)^2 / (sum(rotor.c_bar)*dr)
    rc = ceil(Int,rotor.root_cut*rotor.N_elm)       
    
    lambdat         = zeros(rotor.N_elm)
    dct             = zeros(rotor.N_elm)
    load_sectional  = zeros((rotor.N_elm,6))  
    check           = zeros(rotor.N_elm)
    
    
    iter_count = 0
    ct = 0
    F_loss = 1
    error = 1
    error_arr = Vector{Real}()      

    while iter_count<itermax && error>error_cutoff  
        iter_count = iter_count+1
        for i in rc:rotor.N_elm
            if rotor.r[i]<rotor.root_cut continue end  
            vt =  rotor.r[i]
            vp = lambdat[i]
            v = sqrt(vt^2 + vp^2)
            phi = atan(vp,vt)
            alpha = rotor.stwist[i]+op_params.theta0-phi
            # cl = 2*pi*alpha + 0.16
            # cd = 0.013
            cl,cd,cm = fn_interp(alpha*180/pi,rotor.airfoil)
            cl = cl/(1+1/AR) 
            dcfxr   = -(0.5/pi)*rotor.c_bar[i]*v^2*dr*( cl*sin(phi) + cd*cos(phi) ) 
            dcfzr    = -(0.5/pi)*rotor.c_bar[i]*v^2*dr*( cl*cos(phi) - cd*sin(phi) )      
            dct[i] = -rotor.Nb*dcfzr

            if (iter_count>1) && (tiploss==1)
                # f_root = 0.5*rotor.Nb*rotor.r[i]/((1-rotor.r[i])*phi)
                f_tip = 0.5*rotor.Nb*(1-rotor.r[i])/(rotor.r[i]*phi)
                f_loss = f_tip
                F_loss = (2/pi)*acos( exp(-f_loss) )
            end
            lambda_old = lambdat[i]
            lambda_root = dct[i]/(4*rotor.r[i]*dr*F_loss) + (0.5*op_params.mu[3])^2
            if lambda_root<0 lambda_root = 0 end
            lambda_new = sqrt(lambda_root) - 0.5*op_params.mu[3] 
            lambdat[i] = lambda_old + beta*(lambda_new-lambda_old)

            load_sectional[i,1] = dcfxr
            load_sectional[i,3] = dcfzr
            load_sectional[i,4] =  rotor.r[i] * dcfzr * op_params.rot_flag 
            load_sectional[i,6] = -rotor.r[i] * dcfxr * op_params.rot_flag
        end
        ct_old = ct
        ct = sum(dct)
        error = abs((ct-ct_old)/ct*100)
        append!(error_arr, error)
    end 
    check = error
    load = sum(load_sectional,dims=1)
    load = load[1,:]
    return check,load,load_sectional
end

function bemt_edgewise(rotor,op_params)

    itermax = 100
    error_cutoff = 1e-5
    beta = 0.1
    delpsi = 10

    dr = 1/rotor.N_elm  
    AR = (1-rotor.root_cut)^2 / (sum(rotor.c_bar)*dr)
    rc = ceil(Int,rotor.root_cut*rotor.N_elm) 
    psi = ( 0:delpsi:(360-delpsi) )*(pi/180)

    mux,muy,muz = op_params.mu
    w1,w2,w3 = op_params.w
    mu = sqrt(mux^2+muy^2+muz^2)
    muxy = sqrt(mux^2+muy^2)
    psi_xy = atan(muy,mux)

    load_sectional = zeros(rotor.N_elm,length(psi),6)
    load_azimuth = zeros(length(psi),6)
    load = zeros(6)
    check = zeros(rotor.N_elm,length(psi))

    ct = 0
    lambda0 = 0
    kx = 0
    ky = 0
    iter_count = 0
    error = 1
    error_arr = []
    
    while iter_count<itermax && error>error_cutoff  
        iter_count = iter_count+1
        for j in 1:length(psi)
            for i in rc:rotor.N_elm
                if rotor.r[i]<rotor.root_cut continue end

                lambdai = lambda0 * (1 + kx*rotor.r[i]*cos(psi[j]) + ky*rotor.r[i]*sin(psi[j]) )
                Ut = rotor.r[i] + muxy*sin(psi[j])  - 0 * ( rotor.r[i] * w3 )
                Up = lambdai - muz - 0 *((w1 *sin(psi[j]) + w2* cos(psi[j]))*rotor.r[i])
                V = sqrt(Ut^2 + Up^2)
                phi = atan(Up, Ut)  #in radian (all zero in the first iteration)
                theta = rotor.stwist[i] + op_params.theta0 + op_params.thetac*cos(psi[j]-psi_xy) + op_params.thetas*sin(psi[j]-psi_xy) 
                alpha = theta - phi  
                # cl = 2*pi*alpha + 0.16
                # cd = 0.013
                cl,cd,cm = fn_interp(alpha*180/pi,rotor.airfoil)
                cm = 0  
                cl = cl/(1+1/AR)
                
                dcfxr = -(0.5/pi)*(rotor.c_bar[i])*V^2*dr*( cl*sin(phi) + cd*cos(phi) ) 
                dcfyr = 0
                dcfzr = -(0.5/pi)*(rotor.c_bar[i])*V^2*dr*( cl*cos(phi) - cd*sin(phi) )
                
                dcmxr = rotor.r[i]*dcfzr 
                dcmyr = cm *0
                dcmzr = -rotor.r[i]*dcfxr                
                
                dcfx = dcfxr*sin(psi[j]) - dcfyr*cos(psi[j]) 
                dcfy = dcfxr*cos(psi[j]) + dcfyr*sin(psi[j]) 
                dcfz = dcfzr

                dcmx = dcmxr*sin(psi[j]) - dcmyr*cos(psi[j])
                dcmy = dcmxr*cos(psi[j]) + dcmyr*sin(psi[j])  
                dcmz = dcmzr
                
                load_sectional[i,j,1] = dcfx  
                load_sectional[i,j,2] = dcfy * op_params.rot_flag
                load_sectional[i,j,3] = dcfz
                load_sectional[i,j,4] = dcmx * op_params.rot_flag
                load_sectional[i,j,5] = dcmy  
                load_sectional[i,j,6] = dcmz * op_params.rot_flag
                check[i,j] = cl #phi*180/np.pi
            end            
        end

        load_azimuth = sum(load_sectional,dims=1)
        load_azimuth = load_azimuth[1,:,:]
        load = mean(load_azimuth,dims=1) 
        load = load[1,:]               
                
        ct_old = ct
        ct = -load[3]*rotor.Nb
        lambda0_old = lambda0
        
        if (iter_count==1)
            lambda0_new = sqrt( 0.5* ( sqrt(mu^4+ct^2)-mu^2 ) )
        else
            lambda0_new = 0.5*ct/sqrt(muxy^2+(-muz+lambda0)^2)
        end
        
        lambda0 = lambda0_old + beta* (lambda0_new - lambda0_old)
        
        if (muxy == 0) 
            chi = 0 
        else
            chi = atan(muxy, -muz+lambda0)
        end
        kx = (15*pi/23)*tan(chi/2)

        error = abs((ct-ct_old)/ct*100)
        push!(error_arr,error)
        
    end

    dcm_xy = [ cos(psi_xy) sin(psi_xy) 0; -sin(psi_xy) cos(psi_xy) 0; 0 0 1 ]
    force = [load[1] load[2] load[3] ]
    moment = [load[4] load[5] load[6] ]
    force2 = force*dcm_xy
    moment2 = moment*dcm_xy
    load2 = [force2 moment2]

    check = error_arr

    return check, load2, load_sectional
end

function fn_interp(xeq,data)
    x = data[:,1]
    L = 1
    R = length(x)
    while ((R-L)>1)
        m = floor(Int,(L+R)*0.5)
        if (xeq > x[m])
            L = m
        elseif (xeq < x[m])
            R = m
        else
            L = m
            break
        end
    end

    s = (xeq-x[L])/(x[R]-x[L])
    h1 = 1-s
    h2 = s

    yeq = zeros(size(data,2)-1)
    
    for j = 2:size(data,2)
        yeq[j-1] = h1*data[L,j]+h2*data[R,j]
    end

    return yeq

end



end