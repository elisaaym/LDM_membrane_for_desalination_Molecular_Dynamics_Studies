[x, y, z, n_y, n_z, n] = construct_borophene( 34.125, 55.9, 0,0, 0,5.88 );

fid = fopen('borophene.xyz', 'w');

formatspec = '%d\n\n';

fprintf(fid, formatspec,n);

formatspec = 'B\t %1.3f\t %1.3f\t %1.3f\n';

for i = 1:n
    fprintf(fid, formatspec, x(i), y(i), z(i));
end

fclose(fid);