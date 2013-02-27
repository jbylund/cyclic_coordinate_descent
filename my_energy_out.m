function calc_energy = my_energy_out(coordinates)
%
% calc_energy = my_energy_out(coordinates)
% This function computes the total energy of a conformation and displays
% the 10 closest atom pairs

num_atoms = size(coordinates,1);
calc_energy = 0;

atom1 = zeros(1,3);
atom2 = zeros(1,3);
dist = 0;

c = 1;
all_dists = zeros((num_atoms)*(num_atoms - 1)/2, 3);

for i = 1:1:(num_atoms-1);
	atom1 = coordinates(i,:);
	for j = (i+1):1:num_atoms;
		atom2 = coordinates(j,:);
		dist = sqrt(sum((atom1-atom2).^2));
		calc_energy = calc_energy + (3.4/dist)^12 - (3.4/dist)^6;
		all_dists(c,:) = [i, j, dist];
		c = c + 1;
	end
end

all_dists = sortrows(all_dists, 3);
calc_energy = calc_energy * 10^-6;

min_dists = zeros(10,3);

for i = 1:1:10
	min_dists(i,:) = all_dists(i,:);
end
clear i;

min_dists