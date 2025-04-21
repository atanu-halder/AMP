module inputs_excel
include("../source_directory/AMP.jl")
include("../input_directory/User_module.jl")
using .AMP
using .User_functions
using XLSX
using DelimitedFiles
using Interpolations

function read_inputs(file)
    params = file["Sheet1"]
    Nb = params[1,2]
    R = params[2,2]
    root_cut = params[3,2]
    N_elm = params[4,2]
    dr = 1/N_elm
    r = dr/2:dr:1 
    mtip = params[5,2]
    theta0 = deg2rad(params[6,2])
    thetac = deg2rad(params[7,2])
    thetas = deg2rad(params[8,2])
    rot_flag = params[9,2]
    mu=[params[10,2],params[11,2],params[12,2]]
    w=[params[13,2],params[14,2],params[15,2]]
    naf = params[16,2] # number of airfoil section
    raf = params[17,2]
    haoad = params[18,2] # flag for high angle of data

    af_path = "../data_airfoil"
    file1 = params[19,2]
    full_path = joinpath(af_path, file1)
    data_af = readdlm(full_path)

    if haoad == 0
        airfoil = data_af
    else
        data_high_aoa = readdlm("../data_airfoil/high_aoa.dat")
        airfoil = User_functions.airfoil_table(data_af,data_high_aoa)
    end

    sheet2 = file["chord"]
    row2 = sheet2.dimension.stop.row_number
    chd_data = zeros(row2-1,2)

    for i in 2:row2
        chd_data[i-1,1] = sheet2[i,1]
        chd_data[i-1,2] = sheet2[i,2]
    end

    sheet3 = file["twist"]
    row3 = sheet3.dimension.stop.row_number
    twist_data = zeros(row3-1,2)

    for i in 2:row3
        twist_data[i-1,1] = sheet3[i,1]
        twist_data[i-1,2] = sheet3[i,2]
    end

    interp_chd = linear_interpolation(chd_data[:,1], chd_data[:,2])
    interp_twist = linear_interpolation(twist_data[:,1], twist_data[:,2])

    c_bar = zeros(length(r))
    stwist = zeros(length(r))
    for i in 1:length(r)
        if r[i]<root_cut continue end
        c_bar[i] = interp_chd(r[i])
        stwist[i] = deg2rad(interp_twist(r[i]))
    end

    rotor = AMP.Rotor(Nb,R,root_cut,N_elm,r,c_bar,stwist,airfoil)
    op_params = AMP.Op_params(mtip,theta0,thetac,thetas,mu,w,rot_flag)
    return rotor,op_params

end

end