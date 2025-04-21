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

var1 = 0:15

T_arr = zeros(length(var1))
P_arr = zeros(length(var1))
ct_arr = zeros(length(var1))
cp_arr = zeros(length(var1))
FM_arr = zeros(length(var1))

T2_arr = zeros(length(var1))
P2_arr = zeros(length(var1))
ct2_arr = zeros(length(var1))
cp2_arr = zeros(length(var1))
FM2_arr = zeros(length(var1))

for i in 1:length(var1)
    
    op_params.theta0 = deg2rad(var1[i])

#     #Call BEMT
    check,load,load_sectional = AMP.bemt_axial(rotor,op_params)
    check2,load2,load2_sectional = AMP.bemt_edgewise(rotor,op_params)

#     #Post-processing
    vtip = op_params.mtip*sound_speed
    ct = -rotor.Nb*load[3]
    cp =  rotor.Nb*load[6]
    T = ct*rho*(pi*rotor.R^2)*vtip^2
    P  = cp*rho*(pi*rotor.R^2)*vtip^3
    Q = cp*rho*(pi*rotor.R^2)*vtip^2*rotor.R
    FM = (ct^1.5)/(cp*sqrt(2))
    PL = T/P
    # eta = T*va/P
    T_arr[i] = T
    P_arr[i] = P
    ct_arr[i] = ct
    cp_arr[i] = cp
    FM_arr[i] = FM

    ct2 = -rotor.Nb*load2[3]
    cp2 = rotor.Nb*load2[6]
    T2 = ct2*rho*(pi*rotor.R^2)*vtip^2
    P2  = cp2*rho*(pi*rotor.R^2)*vtip^3
    Q2 = cp2*rho*(pi*rotor.R^2)*vtip^2*rotor.R
    FM2 = (ct2^1.5)/(cp2*sqrt(2))
    PL2 = T2/P2
    # eta2 = T2*va/P2
    T2_arr[i] = T2
    P2_arr[i] = P2
    ct2_arr[i] = ct2
    cp2_arr[i] = cp2
    FM2_arr[i] = FM2

end

# # T_lbs = T_arr/gv.lbs_to_N
# # P_hp = P_arr/gv.hp_to_Watts
ctbysigma = ct_arr/sigma75
cpbysigma = cp_arr/sigma75
ctbysigma2 = ct2_arr/sigma75   
cpbysigma2 = cp2_arr/sigma75


data_exp = XLSX.readtable("../data_validation/xv15_data/xv15_data.xlsx","FM")
exp = hcat(data_exp.data...)

p1=plot(xlabel=L"\textbf{C_T}", ylabel=L"\textbf{FM}",framestyle=:box)
plot!(ct_arr,FM_arr,lw=2,label=L"\textbf{AMP Axial}",legend=:bottomright)
plot!(ct2_arr,FM2_arr,lw=2,ls=:dash,label=L"\textbf{AMP Edgewise}",)
scatter!(exp[:,1],exp[:,2], label = L"\textbf{Exp}")

# # lw1 = 2     # linewidth
# # ms1 = 5     # markersize
# # gfs = 16    # guidefontsize, fontsize of xlabel, ylabel
# # tfs = 10    # tickfontsize
# # lfs = 10    # lengend fontsize

data_exp2 = XLSX.readtable("../data_validation/xv15_data/xv15_data.xlsx","blade_loading_hover_v2")
exp2 = hcat(data_exp2.data...)

p2=plot(xlabel=L"\textbf{C_T/\sigma}", ylabel=L"\textbf{C_Q/\sigma}",framestyle=:box)
plot!(ctbysigma,cpbysigma,lw=2,label=L"\textbf{AMP Axial}")
plot!(ctbysigma2,cpbysigma2,lw=2,ls=:dash,label=L"\textbf{AMP Edgewise}",)
scatter!(exp2[:,1],exp2[:,3], label = L"\textbf{Exp}")

plot(p1,p2,layout=(2,1),size=(360,500))


# # p2 = plot(framestyle=:box,guidefontsize=gfs, tickfontsize = tfs, legendfont = lfs )
# # plot!(xlabel=xl, ylabel="Power (hp)")
# # plot!(var1,P_hp, lw=lw1, lc=:red, label="AMP",)
# # scatter!(exp[:,1],exp[:,7], ms = ms1, mc=:black, label = "Exp")
# # scatter!([cfd[1,1]],[cfd[1,4]], ms = ms1, mc=:green, shape=:dtriangle, label = "CFD")

# # plot(p1,p2,layout=(2,1),size=(600,800))



