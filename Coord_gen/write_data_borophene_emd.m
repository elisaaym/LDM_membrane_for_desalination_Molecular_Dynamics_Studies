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
size_H2O_y =34.325; %Must be a multiple of 3*1.42 to ensure good graphene lattice with periodic BC. Even multiple so we can define z properly
size_H2O_z = 55.887;
size_H2O_x = 56.715;
dia_x =0.911;
n_y=ceil(size_H2O_y/1.42/3);
n_z=ceil(size_H2O_z/(1.42*sqrt(3)));
vol_H2O = size_H2O_y*size_H2O_x*size_H2O_z; %in angstroms
const_a = size_H2O_z/size_H2O_y;
const_b = size_H2O_x/size_H2O_y;
mass_H2O = density_H2O*(vol_H2O/10^24);
num_mol_H2O = mass_H2O/molar_mass_H20*av_k;
num_mol_1D_H2O_y = floor((num_mol_H2O/const_a/const_b)^(1/3));
num_mol_1D_H2O_z = ceil(num_mol_1D_H2O_y*const_a);
num_mol_1D_H2O_x =ceil(num_mol_1D_H2O_y*const_b);
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
space_H2O_y = size_H2O_y/double(num_mol_1D_H2O_y)
space_H2O_x = size_H2O_x/double(num_mol_1D_H2O_x)
space_H2O_z = size_H2O_z/double(num_mol_1D_H2O_z)

%Now specify parameters for salt ions from salinity
Salinity = 0.072;
molar_mass_NA = 22.9898;
molar_mass_CL = 35.452999;
mass_NACL = Salinity*mass_H2O;
num_NACL = int32(av_k/(molar_mass_NA + molar_mass_CL)*mass_NACL);
str = ['Salinity is ', num2str(Salinity)];
disp(str);

%SALT
%OR specify parameters for salt ions from num_NACL
% num_NACL = 0;
% molar_mass_NA = 22.9898;
% molar_mass_CL = 35.452999;
% mass_NACL = (molar_mass_NA + molar_mass_CL)/av_k*num_NACL;
% Salinity = mass_NACL/mass_H2O;
% str = ['Salinity is ', num2str(Salinity)];
% disp(str);

%GRAPHENE
%2*graphene - 1 piston, 1 membrane
bond_length_C = 1.42;
%num_C = cal_C_piston( n_y, n_z, size_H2O_y,size_H2O_z, bond_length_C); %the 2 last term is to add the final bordering case.
%num_C=num_C*2;
slitsize = 6.4;

[ xb, yb, zb, atomtype_B, n_yb, n_zb,num_B ] = construct_borophene( size_H2O_y,size_H2O_z,3+size_H2O_x+3, 0, 0, slitsize );

molar_mass_C = 12.0107;
molar_mass_B = 10.811;
%carbon carbon interactions measured by tersoff potential, not through
%bonds

%Water on the permeate side
layers = 20;
num_mol_H2O_p = (num_mol_1D_H2O_y*num_mol_1D_H2O_z)*layers;

%2 times the water molecules + salt ions + carbon
num_atoms = num_atoms+num_mol_H2O_p*3+num_NACL*2 + num_B;
num_bonds = num_bonds+num_mol_H2O_p*2;
num_angles = num_angles+num_mol_H2O_p;

formatspec = '%d atoms\n';
fprintf(fid, formatspec, num_atoms);
formatspec = '%d bonds\n';
fprintf(fid, formatspec, num_bonds);
formatspec = '%d angles\n';
fprintf(fid, formatspec, num_angles);

fprintf(fid, '\n');

num_atom_types = 13;
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
fprintf(fid, formatspec, 0,size_H2O_x+3+dia_x+3+space_H2O_x*layers+3);
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

mass_array = [molar_mass_O, molar_mass_H, molar_mass_NA, molar_mass_CL, molar_mass_C, molar_mass_B, molar_mass_B, molar_mass_B, molar_mass_B, molar_mass_B, molar_mass_B, molar_mass_B, molar_mass_B];
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
x_loc = 3;
y_loc = 0;
z_loc = 0;
%FOr TIP4P model:
O_charge = -1.0484;
H_charge = 0.5242;
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
        x_loc = 3;
        y_loc = y_loc+space_H2O_y;
    end
    %reset x_loc, y_loc
    x_loc = 3;
    y_loc = 0;
    z_loc = z_loc+space_H2O_z;
end

%Water moleculars for pure water
x_loc = 3+size_H2O_x+3+dia_x+3;
y_loc = 0;
z_loc = 0;

for z_cnt = 1:num_mol_1D_H2O_z
    for y_cnt = 1:num_mol_1D_H2O_y
        for x_cnt = 1:layers
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
        x_loc =  3+size_H2O_x+3+dia_x+3;
        y_loc = y_loc+space_H2O_y;
    end
    %reset x_loc, y_loc
    x_loc =  3+size_H2O_x+3+dia_x+3;
    y_loc = 0;
    z_loc = z_loc+space_H2O_z;
end

for i = 1:num_NACL
    x_dim = ceil(double(i)/double((num_mol_1D_H2O_x-1)^2)); %x_dim cannot exceed num_mol_1D/2 
    y_dim = ceil(double(i)/double((num_mol_1D_H2O_z-1)));
    z_dim =  mod(i, num_mol_1D_H2O_z);
%     if (y_dim==0)
%       y_dim = num_mol_1D_H2O_y-1;
%     end
    %xloc random generate
    x_loc_NA = rand()*(size_H2O_x-10)+10;
    x_loc_CL = rand()*(size_H2O_x-10)+10;
    y_loc_NA = y_dim*space_H2O_y;
    y_loc_CL = y_loc_NA;
    z_loc_NA =  double(z_dim)*space_H2O_z+0.5;
    z_loc_CL = z_loc_NA;
    fprintf(fid, formatspec, atom_ID, mol_ID, 3, 1, x_loc_NA, y_loc_NA, z_loc_NA);
    atom_ID = atom_ID + 1;
    mol_ID = mol_ID + 1;  
    fprintf(fid, formatspec, atom_ID, mol_ID, 4, -1, x_loc_CL, y_loc_CL, z_loc_CL);
    atom_ID = atom_ID + 1;
    mol_ID = mol_ID + 1;  
end
 
%Next, carbon atoms:
%First at the left side of the domain
%initialize
% x_loc_C = 0;
% y_loc_C_i = 0;
% z_loc_C_i = 0;
% disp('left piston start atom id')
% atom_ID
% for z_cnt = 1: n_z
%     for y_cnt = 1: n_y
%         y_loc_C = y_loc_C_i ;
%         z_loc_C = z_loc_C_i + sqrt(3)/2*bond_length_C;
%         if (y_loc_C <= size_H2O_y && z_loc_C <= size_H2O_z)
%             fprintf(fid, formatspec, atom_ID, mol_ID, 5, 0, x_loc_C, y_loc_C, z_loc_C);
%             atom_ID = atom_ID + 1;
%             mol_ID = mol_ID + 1; 
%         end
%         y_loc_C = y_loc_C_i + bond_length_C/2;
%         z_loc_C = z_loc_C_i;
%         if (y_loc_C <= size_H2O_y && z_loc_C <= size_H2O_z)
%             fprintf(fid, formatspec, atom_ID, mol_ID, 5, 0, x_loc_C, y_loc_C, z_loc_C);
%             atom_ID = atom_ID + 1;
%             mol_ID = mol_ID + 1; 
%         end    
%         y_loc_C = y_loc_C_i + bond_length_C*3/2;
%         z_loc_C = z_loc_C_i;        
%         if (y_loc_C <= size_H2O_y && z_loc_C <= size_H2O_z)
%             fprintf(fid, formatspec, atom_ID, mol_ID, 5, 0, x_loc_C, y_loc_C, z_loc_C);
%             atom_ID = atom_ID + 1;
%         end    
%         mol_ID = mol_ID + 1;
%         y_loc_C = y_loc_C_i + bond_length_C*2;
%         z_loc_C = z_loc_C_i + sqrt(3)/2*bond_length_C; 
%         if (y_loc_C <= size_H2O_y && z_loc_C <= size_H2O_z)
%             fprintf(fid, formatspec, atom_ID, mol_ID, 5, 0, x_loc_C, y_loc_C, z_loc_C);
%             atom_ID = atom_ID + 1;
%             mol_ID = mol_ID + 1;
%         end    
%         y_loc_C_i = y_loc_C_i+bond_length_C*3;
%     end  
%     %Move to next z line
%     z_loc_C_i = z_loc_C_i + sqrt(3)*bond_length_C;
%     y_loc_C_i = 0;
% end
% disp('left piston end atom id')
% atom_ID-1

% disp('membrane start atom id')
% atom_ID

%Next middle of domain (membrane, need to take care of holes)
for i=1:num_B
    atomtype = atomtype_B(i)+5;
    fprintf(fid, formatspec, atom_ID, mol_ID, atomtype, 0, xb(i), yb(i), zb(i));
    atom_ID = atom_ID + 1;
    mol_ID = mol_ID + 1; 
 end

% disp('right piston start atom id')
% atom_ID
% %Next end of domain
% x_loc_C = 3+size_H2O_x+3+dia_x+3+space_H2O_x*layers+3;
% y_loc_C_i = 0;
% z_loc_C_i = 0;
% for z_cnt = 1: n_z
%     for y_cnt = 1: n_y
%         y_loc_C = y_loc_C_i ;
%         z_loc_C = z_loc_C_i + sqrt(3)/2*bond_length_C;
%         if (y_loc_C <= size_H2O_y && z_loc_C <= size_H2O_z)        
%             fprintf(fid, formatspec, atom_ID, mol_ID, 5, 0, x_loc_C, y_loc_C, z_loc_C);
%             atom_ID = atom_ID + 1;
%             mol_ID = mol_ID + 1; 
%         end    
%         y_loc_C = y_loc_C_i + bond_length_C/2;
%         z_loc_C = z_loc_C_i;
%         if (y_loc_C <= size_H2O_y && z_loc_C <= size_H2O_z)
%             fprintf(fid, formatspec, atom_ID, mol_ID, 5, 0, x_loc_C, y_loc_C, z_loc_C);
%             atom_ID = atom_ID + 1;
%             mol_ID = mol_ID + 1; 
%         end    
%         y_loc_C = y_loc_C_i + bond_length_C*3/2;
%         z_loc_C = z_loc_C_i;     
%         if (y_loc_C <= size_H2O_y && z_loc_C <= size_H2O_z)
%             fprintf(fid, formatspec, atom_ID, mol_ID, 5, 0, x_loc_C, y_loc_C, z_loc_C);
%             atom_ID = atom_ID + 1;
%             mol_ID = mol_ID + 1;
%         end    
%         y_loc_C = y_loc_C_i + bond_length_C*2;
%         z_loc_C = z_loc_C_i + sqrt(3)/2*bond_length_C; 
%         if (y_loc_C <= size_H2O_y && z_loc_C <= size_H2O_z)
%             fprintf(fid, formatspec, atom_ID, mol_ID, 5, 0, x_loc_C, y_loc_C, z_loc_C);
%             atom_ID = atom_ID + 1;
%             mol_ID = mol_ID + 1;
%         end    
%         y_loc_C_i = y_loc_C_i+bond_length_C*3;
%     end 
%     %Move to next z line
%     z_loc_C_i = z_loc_C_i + sqrt(3)*bond_length_C;
%     y_loc_C_i = 0;
% end
% disp('right piston end atom id')
% atom_ID-1
% 
% fprintf(fid, '\n');


% % % Instead of setting the atoms in fix position, use the minimum energy
% % % configuration available from build_water.f
% % 
% % fprintf(fid, 'Atoms\n\n');
% % %FOr TIP4P model:
% % % O_charge = -1.0484;
% % % H_charge = 0.5242;
% % %For TIP3P model
% % O_charge =  -0.830;
% % H_charge = 0.415;
% % OH_length = 0.9572;
% % atom_ID = 1;
% % formatspec = '%d  %d  %d  %2.4f  %2.4f  %2.4f  %2.4f\n';
% % fid2 = fopen('min_energy_water.xyz', 'r');
% % tline = fgetl(fid2);
% % tline = fgetl(fid2);
% % for i = 1:num_mol
% %   c_temp = textscan(fid2, '%c',1);
% %   x_loc_o = textscan(fid2, '%f',1);
% %   y_loc_o = textscan(fid2, '%f',1);
% %   z_loc_o = textscan(fid2, '%f',1);
% %   f_temp = textscan(fid2, '%f', 1);
% %   fprintf(fid, formatspec, atom_ID, i, 1, O_charge, x_loc_o{1}, y_loc_o{1}, z_loc_o{1});
% %   atom_ID = atom_ID+1;
% %   c_temp = textscan(fid2, '%c',1);
% %   x_loc_h1 = textscan(fid2, '%f',1);
% %   y_loc_h1 = textscan(fid2, '%f',1);
% %   z_loc_h1 = textscan(fid2, '%f',1);
% %   f_temp = textscan(fid2, '%f', 1);
% %   fprintf(fid, formatspec, atom_ID, i, 2, H_charge, x_loc_h1{1}, y_loc_h1{1}, z_loc_h1{1}); 
% %   atom_ID = atom_ID+1;
% %   c_temp = textscan(fid2, '%c',1);
% %   x_loc_h2 = textscan(fid2, '%f',1);
% %   y_loc_h2 = textscan(fid2, '%f',1);
% %   z_loc_h2 = textscan(fid2, '%f',1);
% %   f_temp = textscan(fid2, '%f', 1);
% %   fprintf(fid, formatspec, atom_ID, i, 2, H_charge, x_loc_h2{1}, y_loc_h2{1}, z_loc_h2{1});
% %   atom_ID = atom_ID+1;  
% % end
%
% fclose(fid2);
%   
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
for i = 1:num_mol_H2O+num_mol_H2O_p
    fprintf(fid, formatspec, bond_ID, 1, (i-1)*3+1, (i-1)*3+2);
    bond_ID = bond_ID+1;
    fprintf(fid, formatspec, bond_ID, 1, (i-1)*3+1, (i-1)*3+3);
    bond_ID = bond_ID+1;
end

fprintf(fid, '\n');

%Angles Section
fprintf(fid, 'Angles\n\n');
formatspec = '%d  %d  %d  %d  %d\n';
for i = 1:num_mol_H2O+num_mol_H2O_p
    fprintf(fid, formatspec, i, 1, (i-1)*3+2, (i-1)*3+1, (i-1)*3+3);
end

fprintf(fid, '\n');

fclose(fid);