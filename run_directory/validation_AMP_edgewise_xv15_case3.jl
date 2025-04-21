## validation of AMP
# experimental data :'Performance of a Mach-Scale Coaxial Counter-Rotating Rotor in Hover'
## Developed by Dr. Atanu Halder

cd(@__DIR__)
include("../source_directory/AMP.jl")
include("../input_directory/inputs_excel.jl")
using .AMP
using .inputs_excel
using XLSX
using Plots
using LaTeXStrings
using Interpolations

rho = 1.2256
sound_speed = 343

file = XLSX.readxlsx("../input_directory/inputs_xv15.xlsx")
rotor,op_params = inputs_excel.read_inputs(file)

c_bar_interp = linear_interpolation(rotor.r,rotor.c_bar)
c_bar_75 = c_bar_interp(0.75)
sigma75 = rotor.Nb*c_bar_75/pi

np = 7
var1a = LinRange(6.4,12.7,np)
var1b = LinRange(5.3,11.8,np)
var1c = LinRange(4.,10.8,np)
var1d = LinRange(3.,10,np)
var1e = LinRange(2,9,np)
var1f = LinRange(1,8,np)
var1g = LinRange(0,7,np)
theta_mat = hcat(var1a, var1b, var1c, var1d, var1e, var1f, var1g)
var2 = [-15,-10, -5, 0, 5, 10, 15]

mu = zeros(3)
mu_infty = 0.17
ct_arr = zeros(length(var1a),length(var2))
cp_arr = zeros(length(var1a),length(var2))



for j in 1:length(var2)
    mu[1] = mu_infty*cos(deg2rad(-var2[j]))
    mu[3] = -mu_infty*sin(deg2rad(-var2[j]))
    op_params.mu = mu
    for i in 1:length(var1a)
    
        op_params.theta0 = deg2rad(theta_mat[i,j])

        # Call BEMT
        check,load,load_sectional = AMP.bemt_edgewise(rotor,op_params)

        # Post-processing
        vtip = op_params.mtip*sound_speed

        ct = -rotor.Nb*load[3]
        cp =  rotor.Nb*load[6]
        ct_arr[i,j] = ct
        cp_arr[i,j] = cp
    end

end

# # T_lbs = T_arr/gv.lbs_to_N

# # P_hp = P_arr/gv.hp_to_Watts
ctbysigma = ct_arr/sigma75
cpbysigma = cp_arr/sigma75


data_exp = XLSX.readtable("../data_validation/xv15_data/xv15_data.xlsx","edgewise3")
exp = hcat(data_exp.data...)

p1=plot(xlabel=L"\textbf{C_T/\sigma}", ylabel=L"\textbf{C_Q/\sigma}",framestyle=:box)
plot!(ctbysigma[:,1],cpbysigma[:,1] ,lw=2,label=L"tilt = -15 deg",lc=1)
plot!(ctbysigma[:,2],cpbysigma[:,2] ,lw=2,label=L"tilt = -10 deg",lc=2)
plot!(ctbysigma[:,3],cpbysigma[:,3] ,lw=2,label=L"tilt = -5 deg",lc=3)
plot!(ctbysigma[:,4],cpbysigma[:,4] ,lw=2,label=L"tilt = 0 deg",lc=4)
plot!(ctbysigma[:,5],cpbysigma[:,5] ,lw=2,label=L"tilt = 5 deg",lc=5)
plot!(ctbysigma[:,6],cpbysigma[:,6] ,lw=2,label=L"tilt = 10 deg",lc=6)
plot!(ctbysigma[:,7],cpbysigma[:,7] ,lw=2,label=L"tilt = 15 deg",lc=7)
# plot!(ctbysigma2,eta2_arr,lw=2,ls=:dash,label=L"\textbf{AMP Edgewise}",)
scatter!(exp[:,2],exp[:,3], label="",color=1)
scatter!(exp[:,5],exp[:,6], label="",color=2)
scatter!(exp[:,8],exp[:,9], label="",color=3)
scatter!(exp[:,11],exp[:,12], label="",color=4)
scatter!(exp[:,14],exp[:,15], label="",color=5)
scatter!(exp[:,17],exp[:,18], label="",color=6)
scatter!(exp[:,20],exp[:,21], label="",color=7)
# scatter!(exp[:,11],exp[:,12], label=L"tilt = 0 deg",color=:black)




