function [ x, y, z, atomtype, n_y, n_z,n ] = construct_borophene( size_y, size_z, x_start, y_start, z_start, slitsize )
%construct_borophene This function returns the coordinates of borophene
%   Following the model by Zhou & Jiang 2017 (MD simulation for mechanical
%   properties of borophene: parameterization of valence force field model
%   and stillinger-weber potential

a2 = 1.614;
a1 = 2.866;
h = 0.911;

n_y = ceil(size_y/a2);
n_z = ceil(size_z/a1);

n = 2*n_y*n_z;

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
    for j = 1:n_y/2
        %For atom number 1:
        y_coord = y_init;
        z_coord = z_init;        
        if(y_coord < y_slit_bottom && z_coord <=size_z) 
            y(cnt) = y_init;
            z(cnt) = z_init;
            x(cnt) = x_init;
            
            if (mod(j,2)==1 && mod(k,2)==1)
                atomtype(cnt) = 1;
            elseif (mod(j,2)==0 && mod(k,2)==1)
                atomtype(cnt) = 5;
            elseif (mod(j,2)==1 && mod(k,2)==0)
                atomtype(cnt) =3;
            elseif (mod(j,2)==0 && mod(k,2)==0)
                atomtype(cnt)=7;
            end
            cnt = cnt+1;    
        end
        %For atom number 2:
        y_coord = y_init+ a2/2;
        z_coord = z_init + a1/2;        
        if(y_coord < y_slit_bottom && z_coord <=size_z)         
            y(cnt) = y_init + a2/2;
            z(cnt) = z_init + a1/2;
            x(cnt) = x_init + h;
            if (mod(j,2)==1 && mod(k,2)==1)
                atomtype(cnt) = 2;
            elseif (mod(j,2)==0 && mod(k,2)==1)
                atomtype(cnt) = 6;
            elseif (mod(j,2)==1 && mod(k,2)==0)
                atomtype(cnt) =4;
            elseif (mod(j,2)==0 && mod(k,2)==0)
                atomtype(cnt)=8;
            end            
            cnt = cnt+1;
        end
        %set new y_init
        y_init = y_init + a2;
    end
    %reset y_init
    y_init = y_start;
    %set new z_init
    z_init = z_init + a1;
end 

y_slit_top = max(y)+slitsize;

%next top slab
y_init = y_slit_top;
z_init = z_start;
x_init = x_start;
for k = 1:n_z
    for j = 1:n_y/2
        %For atom number 1:
        y_coord = y_init;
        z_coord = z_init;
        if(y_coord <= size_y&& z_coord <=size_z) 
            y(cnt) = y_init;
            z(cnt) = z_init;
            x(cnt) = x_init;
            if (mod(j,2)==1 && mod(k,2)==1)
                atomtype(cnt) = 1;
            elseif (mod(j,2)==0 && mod(k,2)==1)
                atomtype(cnt) = 5;
            elseif (mod(j,2)==1 && mod(k,2)==0)
                atomtype(cnt) =3;
            elseif (mod(j,2)==0 && mod(k,2)==0)
                atomtype(cnt)=7;
            end
            cnt = cnt+1;
        end

        %For atom number 2:
        y_coord = y_init+ a2/2;
        z_coord = z_init + a1/2;
        if(y_coord <= size_y&& z_coord <=size_z)         
            y(cnt) = y_init + a2/2;
            z(cnt) = z_init + a1/2;
            x(cnt) = x_init + h;
            if (mod(j,2)==1 && mod(k,2)==1)
                atomtype(cnt) = 2;
            elseif (mod(j,2)==0 && mod(k,2)==1)
                atomtype(cnt) = 6;
            elseif (mod(j,2)==1 && mod(k,2)==0)
                atomtype(cnt) =4;
            elseif (mod(j,2)==0 && mod(k,2)==0)
                atomtype(cnt)=8;
            end      
            cnt = cnt+1;
        end
        %set new y_init
        y_init = y_init + a2;
    end
    %reset y_init
    y_init = y_slit_top;
    %set new z_init
    z_init = z_init + a1;
end 

n = cnt-1;
end

