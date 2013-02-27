function coords = read_coords(input_file)
%
% read_coords(input_file) returns an n x 3 matrix of the coords in
% input_file, backbone_native.crd
%

cartesian_coordinates = textread(input_file, '%f');
num_atoms = length(cartesian_coordinates)/3;
coords = zeros(3, num_atoms);

for i = 1:1:length(cartesian_coordinates)
	coords(i) = cartesian_coordinates(i);
end

coords = coords';