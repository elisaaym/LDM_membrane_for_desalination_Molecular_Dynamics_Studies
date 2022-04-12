function [ x, y, z, atomtype, n_y, n_z,n ] = construct_MXeneca( size_y, size_z, x_start, y_start, z_start )
%construct_MXene This function returns the coordinates of a 2D Ti2CF2 MXene sheet
%   Following the crystollogrphic data of Ti2CF2 from Khazaei et al in 2012
%   (Novel electronic and magnetic properties of 2D transition metal
%   carbides and nitrides)

%Crystal information
a = 3.05874;
b = 3.05882;
gamma=119.99394; %degrees
Ti1_x=1.24336;
Ti1_y=0.00072;
Ti1_z=0.00051;
Ti2_x=3.53336;
Ti2_y=1.53009;
Ti2_z=0.00051;
C1_x=2.38842;
C1_y=1.53021;
C1_z=1.76651;
F1_x=0;
F1_y=1.52952;
F1_z=0;
F2_x=4.77672;
F2_y=0;
F2_z=0.88334;

n_y = ceil(size_y/a);
n_z = ceil(size_z/(b*sin(gamma/180*pi)));

n = n_y*n_z*5;

x = zeros(n,1);
y = zeros(n,1);
z = zeros(n,1);
atomtype = zeros (n,1);

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
        if(y_coord <= size_y && z_coord <=size_z && y_coord >= 0)
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 1;
            cnt=cnt+1;
        end
        %For Ti2:
        y_coord = y_init+Ti2_y;
        z_coord = z_init+Ti2_z; 
        x_coord = x_init+Ti2_x;
        if(y_coord <= size_y  && z_coord <=size_z && y_coord >= 0) 
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 1;
            cnt=cnt+1;
        end
        %For C1:
        y_coord = y_init+C1_y;
        z_coord = z_init+C1_z;  
        x_coord = x_init+C1_x;
        if(y_coord <= size_y  && z_coord <=size_z && y_coord >= 0) 
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 2;
            cnt=cnt+1;
        end        
        %For F1:
        y_coord = y_init+F1_y;
        z_coord = z_init+F1_z;  
        x_coord = x_init+F1_x;
        if(y_coord <= size_y  && z_coord <=size_z && y_coord >= 0) 
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 3;
            cnt=cnt+1;
        end   
        %For F2:
        y_coord = y_init+F2_y;
        z_coord = z_init+F2_z;  
        x_coord = x_init+F2_x;
        if(y_coord <= size_y  && z_coord <=size_z && y_coord >= 0) 
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 3;
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

