function [new_coords] = break_loop(coords, loop_start, loop_end)
%
% [new_coords] = break_loop(coords, loop_start, loop_end)
%

[lengths, angles, dihedrals] = cartesian2internal(coords);
num_atoms = size(lengths,1);
rotable = ones(num_atoms, 1);

for i = 1:1:num_atoms
	if(mod(i,3) == 2)
		rotable(i) = 0;
	end
end

start_index = loop_start;
end_index = loop_end;
loop_len = end_index - start_index;

for i = start_index:1:end_index
	if(rotable(i)==1)
		dihedrals(i) = -pi + rand() * 2 * pi;
	end
end

new_coords = con_molecule(lengths, angles, dihedrals);

for i = (end_index + 1):1:num_atoms
	new_coords(i, :) = coords(i, :);
end