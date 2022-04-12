function [ x, y, z, atomtype, charge, n_y, n_z,n ] = construct_MXeneTi3C2O2ca_charge( size_y, size_z, x_start, y_start, z_start )
%construct_MXene This function returns the coordinates of a 2D Ti2CF2 MXene sheet
%   Following the crystollogrphic data of Ti2CF2 from Khazaei et al in 2012
%   (Novel electronic and magnetic properties of 2D transition metal
%   carbides and nitrides)

%Crystal information
a = 3.04196;
b = 3.04189;
gamma=120.00138; %degrees
Ti1_x=3.49092;
Ti1_y=1.52166;
Ti1_z=1.75631;
q_Ti1 = 2.64;
Ti2_x=6.063;
Ti2_y=1.52183;
Ti2_z=8E-05;
q_Ti2 = 2.34;
Ti3_x=0.91876;
Ti3_y=0;
Ti3_z=0.87788;
q_Ti3 = 2.34;
C1_x=2.23419;
C1_y=1.52145;
C1_z=0;
q_C1 = -2.58;
C2_x=4.74752;
C2_y=0.00064;
C2_z=0.87816;
q_C2 = -2.58;
O1_x=0;
O1_y=1.52105;
O1_z=1.75589;
q_O1 = -1.08;
O2_x=6.98176;
O2_y=1.52155;
O2_z=1.75611;
q_O2 = -1.08;

n_y = ceil(size_y/a);
n_z = ceil(size_z/(b*sin(gamma/180*pi)));

n = n_y*n_z*7;

x = zeros(n,1);
y = zeros(n,1);
z = zeros(n,1);
atomtype = zeros (n,1);
charge = zeros (n,1);

%   Bottom slab first
cnt = 1;
y_init = y_start;
z_init = z_start;
x_init = x_start;
for k = 1:n_z
    %Check if range of y values for this cell is at least partially > 0
    while(y_init+a <0)
        y_init = y_init+a;
    end
    y_init_first = y_init;
    for j = 1:n_y
        %For Ti1:
        y_coord = y_init+Ti1_y;
        z_coord = z_init+Ti1_z;
        x_coord = x_init+Ti1_x;
        if(y_coord <= size_y && z_coord <=size_z)
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 1;
            charge(cnt) = q_Ti1;
            cnt=cnt+1;
        end
        %For Ti2:
        y_coord = y_init+Ti2_y;
        z_coord = z_init+Ti2_z; 
        x_coord = x_init+Ti2_x;
        if(y_coord <= size_y  && z_coord <=size_z) 
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 1;
            charge(cnt) = q_Ti2;
            cnt=cnt+1;
        end
        %For Ti3:
        y_coord = y_init+Ti3_y;
        z_coord = z_init+Ti3_z; 
        x_coord = x_init+Ti3_x;
        if(y_coord <= size_y  && z_coord <=size_z) 
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 1;
            charge(cnt) = q_Ti3;
            cnt=cnt+1;
        end        
        %For C1:
        y_coord = y_init+C1_y;
        z_coord = z_init+C1_z;  
        x_coord = x_init+C1_x;
        if(y_coord <= size_y  && z_coord <=size_z) 
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 2;
            charge(cnt) = q_C1;
            cnt=cnt+1;
        end   
        %For C2:
        y_coord = y_init+C2_y;
        z_coord = z_init+C2_z;  
        x_coord = x_init+C2_x;
        if(y_coord <= size_y  && z_coord <=size_z) 
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 2;
            charge(cnt) = q_C2;
            cnt=cnt+1;
        end           
        %For O1:
        y_coord = y_init+O1_y;
        z_coord = z_init+O1_z;  
        x_coord = x_init+O1_x;
        if(y_coord <= size_y  && z_coord <=size_z) 
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 3;
            charge(cnt) = q_O1;
            cnt=cnt+1;
        end   
        %For O2:
        y_coord = y_init+O2_y;
        z_coord = z_init+O2_z;  
        x_coord = x_init+O2_x;
        if(y_coord <= size_y  && z_coord <=size_z) 
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 3;
            charge(cnt) = q_O2;
            cnt=cnt+1;
        end         
        y_init=y_init+a;
    end
    %update y_init for the next column    
    y_init = y_init_first+b*cos(gamma/180*pi);
    %update z_init for the next column
    z_init = z_init+b*sin(gamma/180*pi);
end 

n = cnt-1;
end

