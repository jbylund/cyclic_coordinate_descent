function writecrd(A,filename)
%
% writecrd(A, filename)
% where A is an n x 3 matrix containing atom coordinates and filename is
% the name to put the output

num_atoms = size(A, 1);
A = A';
outid = fopen(filename,'w');
count = 0;
fprintf(outid,'%g atoms\n', num_atoms);
for i = 1:1:(3 * num_atoms)
	fprintf(outid,'%8.3f', A(i));
	if(mod(i, 10) == 0)
		fprintf(outid, '\n');
	end
end
fclose(outid);