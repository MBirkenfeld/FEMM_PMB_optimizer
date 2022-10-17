function s = make_sim_example(R, geometry_params, mat_Stator, mat_Rotor)
    if ~(length(geometry_params) == 2)
        error('length of geometry_params must be 2')
    end

    h = geometry_params(1);
    b = geometry_params(2);
    
    mag_1 = Magnet(b, h, 0, 0, mat_Stator, -90, 1);
    mag_2 = Magnet(b, h, 0, 0, mat_Rotor, 90, 2);

    mag_halbach_3_1 = mag_1.make_halbach(90, -90, 3);          
    mag_halbach_3_2 = mag_2.make_halbach(-90, 90, 3);             

    mag_halbach_3 = [mag_halbach_3_1; mag_halbach_3_2];  

    s = Simulation_FEMM('planar', 2*pi*R, mag_halbach_3);
end