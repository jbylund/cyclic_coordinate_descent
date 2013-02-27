function [new_coords] = con_molecule(bond_lengths, bond_angles, bond_dihedrals)
%
% coords = con_molecule(bond_lengths, bond_angles, bond_dihedrals)
%

bond_lengths = [[0;1;1];bond_lengths];
bond_angles = [[0;0;pi/2];bond_angles];
bond_dihedrals = [[0;0;0];bond_dihedrals];

num_atoms = size(bond_lengths, 1);
new_coords = zeros(num_atoms, 4);

curr_direction = [1, 0, 0, 0];
curr_perpendicular = [0, 0, -1, 0];

% curr_direction
new_coords(2, :) = bond_lengths(2) .* curr_direction; % coords2 correct

curr_direction = new_coords(2, :) - new_coords(1, :); % resets the current direction
curr_direction = curr_direction ./ norm(curr_direction); % resets the current direction

% rotate current direction by bond angle about the current perpendicular
r_mat = make_rot_mat(bond_angles(3), curr_perpendicular);
curr_direction = (r_mat * curr_direction')';
% curr_direction = curr_direction';

% curr_direction
new_coords(3, :) = new_coords(2,:) + bond_lengths(3) .* curr_direction; % coords 3 correct

for i = 4:1:num_atoms
	r_mat = make_rot_mat(bond_dihedrals(i), curr_direction);
	curr_perpendicular = (r_mat * curr_perpendicular')';
	
	r_mat = make_rot_mat(bond_angles(i), curr_perpendicular);
	curr_direction = (r_mat * curr_direction')';
	
	new_coords(i, :) = new_coords(i - 1, :) + bond_lengths(i) .* curr_direction;	
end

new_coords(:,4) = [];

%to undo the hack
new_coords(1,:) = [];
new_coords(1,:) = [];
new_coords(1,:) = [];
%to undo the hack