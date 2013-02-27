function [A] = readcrd(filename,natom)
fid = fopen(filename,'r');
line = fgetl(fid);
A = fscanf(fid,'%f');
fclose(fid);
[r,c] = size(A);
A = reshape(A,natom*3,r/(natom*3));
