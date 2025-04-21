## Modules helpful to define rotor parameters
# Lead Developer: Dr. Atanu Halder
# Co-developers:

module gv

ft_to_m = 0.3048
inch_to_m = 0.0254
lbs_to_N = 4.44822
hp_to_Watts = 745.7
knots_to_mps = 0.514444

end

module User_functions
using LinearAlgebra
using StaticArrays
using Interpolations

function airfoil_table(data_af,data_high_aoa)
    len1 = size(data_af,1)
    len2 = size(data_high_aoa,1)

    master_table = zeros(2*len2+len1,4)

    for i in 1: size(master_table,1)
        if i<=len2
            master_table[i,1] = -data_high_aoa[len2-i+1,1]
            master_table[i,2] = -data_high_aoa[len2-i+1,2]
            master_table[i,3] = data_high_aoa[len2-i+1,3]
        elseif (i>len2) && (i<=(len2+len1))
            master_table[i,1] = data_af[i-len2,1]
            master_table[i,2] = data_af[i-len2,2]
            master_table[i,3] = data_af[i-len2,3]
            master_table[i,4] = data_af[i-len2,5]
        elseif (i>(len1+len2))
            master_table[i,1] = data_high_aoa[i-len2-len1,1]
            master_table[i,2] = data_high_aoa[i-len2-len1,2]
            master_table[i,3] = data_high_aoa[i-len2-len1,3]            
        end
    end
    return master_table
end

function find_twist_linear(twist,root_cut,r)
    stwist = zeros(length(r))
    theta_tw = twist/(1-root_cut)
    for i in 1:length(r)
        if r[i]>=root_cut
            stwist[i] = theta_tw*(r[i]-0.75)
        end
    end
    return stwist
end

function find_chord_linear(sigma,Nb,taper,root_cut,r)
    cbar_75 = sigma*pi/Nb
    cbar_tip = cbar_75/(1+0.25*(taper-1)/(1-root_cut))
    cbar_root = cbar_tip*taper
    c_bar = zeros(length(r))    
    for i in 1:length(r)
        if r[i]>=root_cut
            c_bar[i] = cbar_tip*( 1 + (taper-1)* (1-r[i]) / (1-root_cut) )
        end
    end
    AR_exact = 2*(1-root_cut)/(cbar_root+cbar_tip)   
    # print('AR_exact = ', AR_exact)    
    
    return c_bar
end



end