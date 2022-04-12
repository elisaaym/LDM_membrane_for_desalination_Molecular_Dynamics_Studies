function [ x, y, z, atomtype, n_y, n_z,n ] = construct_MoS2( size_y, size_z, x_start, y_start, z_start, slitsize )
%construct_MXene This function returns the coordinates of a 2D Ti2CF2 MXene sheet
%   Following the crystollogrphic data of Ti2CF2 from Khazaei et al in 2012
%   (Novel electronic and magnetic properties of 2D transition metal
%   carbides and nitrides)

%Crystal information
bondlength = 3.17;
change_y = bondlength*sin(60/180*pi);
h = 3.241/2;
Mo_changey = bondlength*sqrt(3)/3;
Mo_changez = bondlength/2;

n_y = ceil(size_y/change_y/2);
n_z = ceil(size_z/bondlength);

n = n_y*n_z*6;

x = zeros(n,1);
y = zeros(n,1);
z = zeros(n,1);
atomtype = zeros (n,1);

y_slit_bottom = size_y/2-slitsize/2;

%   Bottom slab first
cnt = 1;
y_init = y_start;
z_init = z_start;
x_init = x_start;
for k = 1:n_z
    %Check if range of y values for this cell is at least partially > 0
    for j = 1:n_y/2
        %For first layer of S:
        y_coord = y_init+change_y;
        z_coord = z_init;
        x_coord = x_init;
        if(y_coord < y_slit_bottom-change_y/2 && z_coord <=size_z && y_coord >= 0)
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 1;
            cnt=cnt+1;
        end
        y_coord = y_init;
        z_coord = z_init+bondlength/2;
        x_coord = x_init;
        if(y_coord < y_slit_bottom-change_y/2 && z_coord <=size_z && y_coord >= 0)
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 1;
            cnt=cnt+1;
        end     
        %For middle layer of Mo:
        y_coord = y_init + Mo_changey;
        z_coord = z_init + Mo_changez;
        x_coord = x_init+h;
        if(y_coord < y_slit_bottom && z_coord <=size_z && y_coord >= 0)
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 2;
            cnt=cnt+1;
        end        
        y_coord = y_init + Mo_changey - change_y;
        z_coord = z_init + Mo_changez + bondlength/2;
        x_coord = x_init+h;
        if(y_coord < y_slit_bottom && z_coord <=size_z && y_coord >= 0)
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 2;
            cnt=cnt+1;
        end        
        %For second layer of S:
        y_coord = y_init+change_y;
        z_coord = z_init;
        x_coord = x_init+2*h;
        if(y_coord < y_slit_bottom-change_y/2 && z_coord <=size_z && y_coord >= 0)
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 1;
            cnt=cnt+1;
        end
        y_coord = y_init;
        z_coord = z_init+bondlength/2;
        x_coord = x_init;
        if(y_coord < y_slit_bottom-change_y/2 && z_coord <=size_z && y_coord >= 0)
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord+2*h;
            atomtype(cnt) = 1;
            cnt=cnt+1;
        end     
        y_init=y_init+change_y*2;
    end
    %update y_init for the next column    
    y_init = y_start;
    %update z_init for the next column
    z_init = z_init+bondlength;
end 

y_slit_top = max(y)+slitsize;

y_init = y_slit_top;
z_init = z_start;
x_init = x_start;
for k = 1:n_z
    %Check if range of y values for this cell is at least partially > 0
    for j = 1:n_y/2
        %For first layer of S:
        y_coord = y_init+change_y;
        z_coord = z_init;
        x_coord = x_init;
        if(y_coord <= size_y && z_coord <=size_z && y_coord >= y_slit_top+change_y)
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 1;
            cnt=cnt+1;
        end
        y_coord = y_init;
        z_coord = z_init+bondlength/2;
        x_coord = x_init;
        if(y_coord <= size_y && z_coord <=size_z && y_coord >= y_slit_top+change_y)
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 1;
            cnt=cnt+1;
        end     
        %For middle layer of Mo:
        y_coord = y_init + Mo_changey;
        z_coord = z_init + Mo_changez;
        x_coord = x_init+h;
        if(y_coord <= size_y && z_coord <=size_z && y_coord >= y_slit_top)
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 2;
            cnt=cnt+1;
        end        
        y_coord = y_init + Mo_changey - change_y;
        z_coord = z_init + Mo_changez + bondlength/2;
        x_coord = x_init+h;
        if(y_coord <= size_y && z_coord <=size_z && y_coord >= y_slit_top)
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 2;
            cnt=cnt+1;
        end        
        %For second layer of S:
        y_coord = y_init+change_y;
        z_coord = z_init;
        x_coord = x_init+2*h;
        if(y_coord <= size_y && z_coord <=size_z && y_coord >= y_slit_top+change_y)
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord;
            atomtype(cnt) = 1;
            cnt=cnt+1;
        end
        y_coord = y_init;
        z_coord = z_init+bondlength/2;
        x_coord = x_init;
        if(y_coord <= size_y && z_coord <=size_z && y_coord >= y_slit_top+change_y)
            y(cnt) = y_coord;
            z(cnt) = z_coord;
            x(cnt) = x_coord+2*h;
            atomtype(cnt) = 1;
            cnt=cnt+1;
        end     
        y_init=y_init+change_y*2;
    end
    %update y_init for the next column    
    y_init = y_slit_top;
    %update z_init for the next column
    z_init = z_init+bondlength;
end 

n = cnt-1;
end

