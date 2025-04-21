## Developed by Dr. Atanu Halder

cd(@__DIR__)
include("../source_directory/AMP.jl")
include("../input_directory/inputs.jl")
include("../input_directory/inputs_excel.jl")
using .AMP
using .inputs
using .inputs_excel
using XLSX

rho = 1.2256
sound_speed = 343

# rotor = inputs.define_rotor()                           # define rotor parameters
# op_params = inputs.define_op_params()                   # define operational parameters     
file = XLSX.readxlsx("../input_directory/inputs.xlsx")
rotor,op_params = inputs_excel.read_inputs(file)
check,load,load_sectional = AMP.bemt_axial(rotor,op_params)  # call bemt
check2,load2,load2_sectional = AMP.bemt_edgewise(rotor,op_params)

#Post-processing
vtip = op_params.mtip*sound_speed
va = op_params.mu[3]*sound_speed

ct = -rotor.Nb*load[3]
cp =  rotor.Nb*load[6]
T = ct*rho*(pi*rotor.R^2)*vtip^2
P  = cp*rho*(pi*rotor.R^2)*vtip^3
Q = cp*rho*(pi*rotor.R^2)*vtip^2*rotor.R
FM = (ct^1.5)/(cp*sqrt(2))
PL = T/P
eta = T*va/P

ct2 = -rotor.Nb*load2[3]
cp2 = rotor.Nb*load2[6]
T2 = ct2*rho*(pi*rotor.R^2)*vtip^2
P2  = cp2*rho*(pi*rotor.R^2)*vtip^3
Q2 = cp2*rho*(pi*rotor.R^2)*vtip^2*rotor.R
FM2 = (ct2^1.5)/(cp2*sqrt(2))
PL2 = T2/P2
eta2 = T2*va/P2

# println(check)
println("FM = ",FM)
println("FM2 = ",FM2)
println("Thrust (N) = ", T)
println("Thrust-2 (N) = ", T2)
println("Power (W) = ", P)
println("Power-2 (W) = ", P2)
println("Torque (Nm) = ", Q)
println("Torque-2 (Nm) = ", Q2)
println("Propulsive Efficiency = ", eta)
println("Propulsive Efficiency 2 = ", eta2)
# println(load_edgewise)






