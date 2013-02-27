rot1;
curr_max = my_energy(coords2)
curr_energy = 0;
max_index = 0;

for i = 1:1:num_atoms
	my_dihedrals = bond_dihedrals;
	my_dihedrals(i) = my_dihedrals(i) + pi/3;
	new_coords = con_molecule(num_atoms, bond_lengths, bond_angles, my_dihedrals);
	curr_energy = my_energy(new_coords);
	diff = curr_energy - curr_max;
	if curr_energy > curr_max
		curr_max = curr_energy;
		max_index = i;
	end
end

my_dihedrals = bond_dihedrals;
my_dihedrals(max_index) = my_dihedrals(max_index) + pi/3;
new_coords = con_molecule(num_atoms, bond_lengths, bond_angles, my_dihedrals);
max_index
curr_max
clear i;
clear max_index;
clear curr_energy;