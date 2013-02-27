function calc_energy = my_energy(coordinates)
%
% calc_energy = my_energy(coordinates)
% Provides an energy estimate of the conformation in coordinates

num_atoms = size(coordinates,1);
calc_energy = 0;

atom1 = zeros(1,3);
atom2 = zeros(1,3);
dist = 0;

for i = 1:1:(num_atoms-1);
	atom1 = coordinates(i,:);
	for j = (i+1):1:num_atoms;
		atom2 = coordinates(j,:);
		dist = sqrt(sum((atom1-atom2).^2));
		calc_energy = calc_energy + (3.4/dist)^12 - (3.4/dist)^6;
	end
end

calc_energy = calc_energy * 10^-6;