%This script writes a data file as input to read_data in LAMMPS
fid = fopen('data1.system', 'w');

formatspec = 'Start File for LAMMPS\n';
fprintf(fid, formatspec);

fprintf(fid, '\n');

%Write header of the data file

%WATER
%Define parameters for water
density_H2O = 1; %g/cm^3
molar_mass_O = 15.9994;
molar_mass_H = 1.008;
molar_mass_H20 = molar_mass_O + 2*molar_mass_H; %g/mol
av_k = 6.02214086*10^23; %avogadro constant mol^-1

%n = 9; %Multiple of 3*1.42 in the y direction, n must be divisible by 3
bondlength = 3.17;
change_y = bondlength*sin(60/180*pi);
size_H2O_y =44*change_y; %Must be a multiple of a2 to ensure good lattice with periodic BC. Even multiple so we can define z properly
size_H2O_x = 38*bondlength; %Multiple of a1
size_H2O_z = 120;
n_y = ceil(size_H2O_y/change_y/2);
n_z = ceil(size_H2O_z/bondlength);
size_H2O_y_w =40; %Must be a multiple of 3*1.42 to ensure good graphene lattice with periodic BC. Even multiple so we can define z properly
size_H2O_x_w = 40;
size_H2O_z_w = 30;
vol_H2O = size_H2O_x_w*size_H2O_y_w*size_H2O_z_w; %in angstroms
const_a = size_H2O_z_w/size_H2O_y_w;
const_b = size_H2O_x_w/size_H2O_y_w;
mass_H2O = density_H2O*(vol_H2O/10^24);
num_mol_H2O = mass_H2O/molar_mass_H20*av_k;
num_mol_1D_H2O_y = ceil((num_mol_H2O/const_a/const_b)^(1/3));
num_mol_1D_H2O_z = floor(num_mol_1D_H2O_y*const_a);
num_mol_1D_H2O_x =floor(num_mol_1D_H2O_y*const_b);
num_mol_H2O = num_mol_1D_H2O_y*num_mol_1D_H2O_z*num_mol_1D_H2O_x;
density_check = (molar_mass_H20*num_mol_H2O/av_k)/(vol_H2O*10^-24)

%compute other parameters for water
% mass_H2O = num_mol_H2O/av_k*molar_mass_h20; %in grams
% vol_H2O = mass_H2O/density_H2O*10^24; %in angstroms^3
% size_H2O = ceil(vol_H2O^(1/3)*100)/100;
% num_mol_1D_H2O = int8(num_mol_H2O^(1/3));
num_atoms = num_mol_H2O*3;
num_bonds = num_mol_H2O*2;
num_angles = num_mol_H2O;
space_H2O_y = size_H2O_y_w/double(num_mol_1D_H2O_y)
space_H2O_x = size_H2O_x_w/double(num_mol_1D_H2O_x)
space_H2O_z = size_H2O_z_w/double(num_mol_1D_H2O_z)

%Now specify parameters for salt ions from salinity
Salinity = 0;
molar_mass_NA = 22.9898;
molar_mass_CL = 35.452999;
mass_NACL = Salinity*mass_H2O;
num_NACL = int32(av_k/(molar_mass_NA + molar_mass_CL)*mass_NACL);
str = ['Salinity is ', num2str(Salinity)];
disp(str);

molar_mass_C = 12.0107;
molar_mass_S = 32.065;
molar_mass_Mo = 95.94;

[ zm, ym, xm, atomtype_m, n_ym, n_xm,n_m ] = construct_MoS2ca_Mo(size_H2O_y, size_H2O_x,0, 0, 0 );



%2 times the water molecules + salt ions + carbon
num_atoms = num_atoms + n_m;

formatspec = '%d atoms\n';
fprintf(fid, formatspec, num_atoms);
formatspec = '%d bonds\n';
fprintf(fid, formatspec, num_bonds);
formatspec = '%d angles\n';
fprintf(fid, formatspec, num_angles);

fprintf(fid, '\n');


num_atom_types = 7;
num_bond_types = 1;
num_angle_types = 1;

formatspec = '%d atom types\n';
fprintf(fid, formatspec, num_atom_types);
formatspec = '%d bond types\n';
fprintf(fid, formatspec, num_bond_types);
formatspec = '%d angle types\n';
fprintf(fid, formatspec, num_angle_types);

fprintf(fid, '\n');
% Box from structured configuration
formatspec = '%1.3f %1.3f xlo xhi\n';
fprintf(fid, formatspec, 0,size_H2O_x);
formatspec = '%1.3f %1.3f ylo yhi\n';
fprintf(fid, formatspec, 0, size_H2O_y);
formatspec = '%1.3f %1.3f zlo zhi\n';
fprintf(fid, formatspec, 0, size_H2O_z);

% Box for minimum energy configuration
% formatspec = '%1.1f %1.1f xlo xhi\n';
% fprintf(fid, formatspec, -size/2, size/2);
% formatspec = '%1.1f %1.1f ylo yhi\n';
% fprintf(fid, formatspec, -size/2, size/2);
% formatspec = '%1.1f %1.1f zlo zhi\n';
% fprintf(fid, formatspec, -size/2, size/2);

fprintf(fid, '\n');

%Start to write body
%Mass section
fprintf(fid, 'Masses\n\n');

mass_array = [molar_mass_O, molar_mass_H, molar_mass_NA, molar_mass_CL, molar_mass_C, molar_mass_S, molar_mass_Mo];
for i = 1: num_atom_types
    fprintf(fid, '%d  %2.4f\n', i, mass_array(i));
end


fprintf(fid, '\n');

%Bond coeffs section
fprintf(fid, 'Bond Coeffs\n\n');
bond_coeff = zeros(num_bond_types, 2);
%FOr tip3P
bond_coeff(1,1) = 450;
%For TIP4P
%bond_coeff(1,1) = 0;
bond_coeff(1,2) = 0.9572;
for i = 1:num_bond_types
    fprintf(fid, '%d  %d  %1.4f\n', i, bond_coeff(i,1), bond_coeff(i,2));
end

fprintf(fid, '\n');

%Angle Coeffs section
fprintf(fid, 'Angle Coeffs\n\n');
angle_coeff = zeros(num_angle_types, 2);
%For TIP3P
angle_coeff(1,1) = 55;
%fOR tip4p
%angle_coeff(1,1) = 0;
angle_coeff(1,2) = 104.52;

for i = 1:num_angle_types
    fprintf(fid, '%d  %d  %1.4f\n', i, angle_coeff(i,1), angle_coeff(i,2));
end

fprintf(fid, '\n');

%Atoms section
%Oxygen atom of the first water molecule is at 0, others are spaced 3 A
%apart, first in x dir, then y, then z
fprintf(fid, 'Atoms\n\n');

%Water molecules for salt water
%initialize
x_loc = 20;
y_loc = 20;
z_loc = 3.5;
%FOr TIP4P model:
O_charge = -1.040;
H_charge = 0.52;
%For TIP3P model
% O_charge =  -0.830;
% H_charge = 0.415;
OH_length = 0.9572;

HOH_angle = 104.52*pi/180;
OH_length_sin = OH_length*sin(HOH_angle/2);
OH_length_cos = OH_length*cos(HOH_angle/2);
mol_ID = 1;
atom_ID = 1;
formatspec = '%d  %d  %d  %2.4f  %2.4f  %2.4f  %2.4f\n';

for z_cnt = 1:num_mol_1D_H2O_z
    for y_cnt = 1:num_mol_1D_H2O_y
        for x_cnt = 1:num_mol_1D_H2O_x
            x_loc_o = x_loc + space_H2O_x/2;
            y_loc_o = y_loc + space_H2O_y-0.5;
            z_loc_o = z_loc + space_H2O_z/2;
            fprintf(fid, formatspec, atom_ID, mol_ID, 1, O_charge, x_loc_o, y_loc_o, z_loc_o);
            atom_ID = atom_ID+1;
            x_loc_h1 = x_loc_o - OH_length_sin;
            y_loc_h1 = y_loc_o - OH_length_cos;
            z_loc_h1 = z_loc_o;
            fprintf(fid, formatspec, atom_ID, mol_ID, 2, H_charge, x_loc_h1, y_loc_h1, z_loc_h1);
            atom_ID = atom_ID+1;
            x_loc_h2 = x_loc_o + OH_length_sin;
            y_loc_h2 = y_loc_o - OH_length_cos;
            z_loc_h2 = z_loc_o;
            fprintf(fid, formatspec, atom_ID, mol_ID, 2, H_charge, x_loc_h2, y_loc_h2, z_loc_h2);
            atom_ID = atom_ID+1;            
            mol_ID = mol_ID+1;
            x_loc = x_loc + space_H2O_x;
        end
        %reset x_loc
        x_loc = 20;
        y_loc = y_loc+space_H2O_y;
    end
    %reset x_loc, y_loc
    x_loc = 20;
    y_loc = 20;
    z_loc = z_loc+space_H2O_z;
end

 
%Next, boron atoms:
%Next middle of domain (membrane, need to take care of holes)
for i=1:n_m
    atomtype = atomtype_m(i)+5;
    fprintf(fid, formatspec, atom_ID, mol_ID, atomtype, 0, xm(i), ym(i), zm(i));
    atom_ID = atom_ID + 1;
    mol_ID = mol_ID + 1; 
 end

fprintf(fid, '\n');
 
%Velocities Section - all velocities to 0
fprintf(fid, 'Velocities\n\n');
formatspec = '%d  %1.1f  %1.1f  %1.1f\n';
for i = 1:num_atoms
    fprintf(fid, formatspec, i, 0, 0, 0);
end  

fprintf(fid, '\n');

%Bonds Section
fprintf(fid, 'Bonds\n\n');
formatspec = '%d  %d  %d  %d\n';
bond_ID=1;
for i = 1:num_mol_H2O
    fprintf(fid, formatspec, bond_ID, 1, (i-1)*3+1, (i-1)*3+2);
    bond_ID = bond_ID+1;
    fprintf(fid, formatspec, bond_ID, 1, (i-1)*3+1, (i-1)*3+3);
    bond_ID = bond_ID+1;
end

fprintf(fid, '\n');

%Angles Section
fprintf(fid, 'Angles\n\n');
formatspec = '%d  %d  %d  %d  %d\n';
for i = 1:num_mol_H2O
    fprintf(fid, formatspec, i, 1, (i-1)*3+2, (i-1)*3+1, (i-1)*3+3);
end

fprintf(fid, '\n');

fclose(fid);