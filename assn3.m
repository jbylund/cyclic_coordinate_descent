coords = read_coords('backbone_native.crd');

broken_coords = coords;
recon_coords = coords;
[lengths, angles, dihedrals] = cartesian2internal(coords);
start_index = 34;
end_index = 46;

for i = 1:1:10
	filename = sprintf('open_coords%d.crd', i);
	broken_coords = break_loop(coords, start_index * 3 - 2, end_index * 3);
	writecrd(broken_coords, filename);
	[lengths, angles, dihedrals] = ccd(broken_coords, start_index * 3 - 2, end_index * 3, coords);
	recon_coords = con_molecule(lengths, angles, dihedrals);
	filename = sprintf('closed_coords%d.crd', i);
	writecrd(recon_coords, filename);
end

clear start_index;
clear end_index;
clear i;
clear filename;