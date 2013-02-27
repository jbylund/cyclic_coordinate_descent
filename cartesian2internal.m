function [bond_lengths, bond_angles, bond_dihedrals] = cartesian2internal(coords2)
%
% [bond_lengths, bond_angles, bond_dihedrals] = cartesian2internal(coords)
% This function converts from cartesian coordinates to internal coordinates

num_atoms = size(coords2,1);

r = [0, 0, 0; 1, 0, 0; 1, 1, 0];
coords2 = [r;coords2];
num_atoms = num_atoms + 3; % for the hack

bond_lengths = zeros(num_atoms, 1);
temp = zeros(1,3);

for i = 2:1:num_atoms
	temp = coords2(i-1, :) - coords2(i, :);
	bond_lengths(i) = sqrt(sum(temp.^2));
end

v1 = zeros(1,3);
v2 = zeros(1,3);

bond_angles = zeros(num_atoms, 1);
for i = 3:1:num_atoms
	v1 = coords2(i-2,:) - coords2(i-1,:);
	v2 = coords2(i,:) - coords2(i-1,:); 
	bond_angles(i) = pi - acos(dot(v1, v2)/(norm(v1) * norm(v2)));
end
clear v1;
clear v2;

atom1 = zeros(1,3);
atom2 = zeros(1,3);
atom3 = zeros(1,3);
atom4 = zeros(1,3);
v1 = zeros(1,3);
v2 = zeros(1,3);
v3 = zeros(1,3);
normal1 = zeros(1,3); % normal to first plane
normal2 = zeros(1,3); % normal to second plane
tcos = 0;
tsin = 0;

bond_dihedrals = zeros(num_atoms, 1); % need to use atan 2
for i = 4:1:num_atoms
	atom1 = coords2(i - 3, :);
	atom2 = coords2(i - 2, :);
	atom3 = coords2(i - 1, :);
	atom4 = coords2(i + 0, :);
	v1 = atom1 - atom2;
	v2 = atom2 - atom3; % finding the dihedral about v2
	v3 = atom3 - atom4;
	normal1 = cross(v1, v2);
	normal2 = cross(v2, v3);
	normal1 = normal1/(norm(normal1));
	normal2 = normal2/(norm(normal2));
	tcos = dot(normal1, normal2);
	tsin = dot(v2, cross(normal1, normal2))/norm(v2);
	bond_dihedrals(i) = atan2(tsin, tcos);
end

% undo the hack
for i = 1:1:3
	coords2(1,:) = [];
	bond_lengths(1,:) = [];
	bond_angles(1,:) = [];
	bond_dihedrals(1,:) = [];
end
num_atoms = num_atoms - 3;
% undo the hack