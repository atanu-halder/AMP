## validation of AMP

## Developed by Dr. Atanu Halder

cd(@__DIR__)
include("../source_directory/AMP.jl")
include("../input_directory/inputs_excel.jl")
include("../input_directory/User_module.jl")
using .AMP
using .inputs_excel
using .gv
using XLSX
using Plots
using LaTeXStrings


rho = 1.2256
sound_speed = 343

file = XLSX.readxlsx("../input_directory/inputs_gofly.xlsx")
rotor,op_params = inputs_excel.read_inputs(file)

var1 = 200:200:1600

T_arr = zeros(length(var1))
P_arr = zeros(length(var1))
T2_arr = zeros(length(var1))
P2_arr = zeros(length(var1))

for i in 1:length(var1)

    rpm = var1[i]
    vtip = rpm*2*pi/60*rotor.R
    op_params.mtip = vtip/sound_speed

    #Call BEMT
    check,load,load_sectional = AMP.bemt_axial(rotor,op_params)
    check2,load2,load2_sectional = AMP.bemt_edgewise(rotor,op_params)

    #Post-processing

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

end

T_lbs = T_arr/gv.lbs_to_N
P_hp = P_arr/gv.hp_to_Watts

T2_lbs = T2_arr/gv.lbs_to_N
P2_hp = P2_arr/gv.hp_to_Watts

data_exp = XLSX.readtable("../data_validation/exp_gofly_subscale_isolated.xlsx","Sheet1")
exp = hcat(data_exp.data...)
data_cfd = XLSX.readtable("../data_validation/cfd_gofly_subscale_isolated.xlsx","Sheet1")
cfd = hcat(data_cfd.data...)

xl = L"\textbf{RPM}"

p1=plot(xlabel=xl, ylabel=L"\textbf{Thrust (lbs)}",framestyle=:box)
plot!(var1,T_lbs,lw=2,label=L"\textbf{AMP Axial}")
plot!(var1,T2_lbs,lw=2,ls=:dash,label=L"\textbf{AMP Edgewise}")
scatter!(exp[:,1],exp[:,3], label = L"\textbf{EXP}")
scatter!([cfd[1,1]],[cfd[1,3]], label = L"\textbf{CFD}")

p2=plot(xlabel=xl, ylabel=L"\textbf{Power (hp)}",framestyle=:box)
plot!(var1,P_hp,lw=2,label=L"\textbf{AMP Axial}",)
plot!(var1,P2_hp,lw=2,s=:dash, label=L"\textbf{AMP Edgewise}",)
scatter!(exp[:,1],exp[:,7],label = L"\textbf{EXP}")
scatter!([cfd[1,1]],[cfd[1,4]], label = L"\textbf{CFD}")

plot(p1,p2,layout=(2,1),size=(360,500))


# lw1 = 2     # linewidth
# ms1 = 5     # markersize
# gfs = 16    # guidefontsize, fontsize of xlabel, ylabel
# tfs = 10    # tickfontsize
# lfs = 10    # lengend fontsize

# p1 = plot(framestyle=:box,guidefontsize=gfs, tickfontsize = tfs, legendfont = lfs )
# plot!(xlabel=xl, ylabel="Thrust (lbs)")
# plot!(var1,T_lbs, lw=lw1, lc=:red, label="AMP",)
# scatter!(exp[:,1],exp[:,3], ms = ms1, mc=:black, label = "Exp")
# scatter!([cfd[1,1]],[cfd[1,3]], ms = ms1, mc=:green, shape=:dtriangle, label = "CFD")

# p2 = plot(framestyle=:box,guidefontsize=gfs, tickfontsize = tfs, legendfont = lfs )
# plot!(xlabel=xl, ylabel="Power (hp)")
# plot!(var1,P_hp, lw=lw1, lc=:red, label="AMP",)
# scatter!(exp[:,1],exp[:,7], ms = ms1, mc=:black, label = "Exp")
# scatter!([cfd[1,1]],[cfd[1,4]], ms = ms1, mc=:green, shape=:dtriangle, label = "CFD")

# plot(p1,p2,layout=(2,1),size=(600,800))



