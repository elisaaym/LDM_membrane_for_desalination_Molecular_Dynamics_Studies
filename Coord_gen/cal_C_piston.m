function [ num_C ] = cal_C_piston( n_y, n_z, size_H2O_y,size_H2O_z, bond_length_C)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
num_C = 0;
y_loc_C_i = 0;
z_loc_C_i = 0;
for z_cnt = 1: n_z
    for y_cnt = 1: n_y
        y_loc_C = y_loc_C_i ;
        z_loc_C = z_loc_C_i + sqrt(3)/2*bond_length_C;
        if (y_loc_C <= size_H2O_y && z_loc_C <= size_H2O_z)
            num_C = num_C+1;
        end        
        y_loc_C = y_loc_C_i + bond_length_C/2;
        z_loc_C = z_loc_C_i;
        if (y_loc_C <= size_H2O_y && z_loc_C <= size_H2O_z)
            num_C = num_C+1;
        end     
        y_loc_C = y_loc_C_i + bond_length_C*3/2;
        z_loc_C = z_loc_C_i;  
        if (y_loc_C <= size_H2O_y && z_loc_C <= size_H2O_z)
            num_C = num_C+1;
        end     
        y_loc_C = y_loc_C_i + bond_length_C*2;
        z_loc_C = z_loc_C_i + sqrt(3)/2*bond_length_C;
        if (y_loc_C <= size_H2O_y && z_loc_C <= size_H2O_z)
            num_C = num_C+1;
        end     
        y_loc_C_i = y_loc_C_i+bond_length_C*3;
    end 
    %Move to next z line
    z_loc_C_i = z_loc_C_i + sqrt(3)*bond_length_C;
    y_loc_C_i = 0;
end
end

