function [new_lengths, new_angles, new_dihedrals] = ccd(coords, loop_start, loop_end, old_coords)
%
% [new_dihedrals] = ccd(coords, loop_start, loop_end, old_coords)
% 

delta = 10^-6;

[real_lengths, real_angles, real_dihedrals] = cartesian2internal(old_coords);
[broken_lengths, broken_angles, broken_dihedrals] = cartesian2internal(coords);

num_atoms = size(real_lengths, 1);
rotable = ones(num_atoms, 1);

for i = 1:1:num_atoms
	if(mod(i,3) == 2)
		rotable(i) = 0;
	end
end

new_molecule = con_molecule(real_lengths, real_angles, broken_dihedrals);

start_index = loop_start;
end_index = loop_end;
target = old_coords(end_index, :);
loop_len = end_index - start_index + 1;

%rand('state',0);

last_index = -1;
curr_index = floor(rand()*loop_len) + start_index;

v0 = zeros(1,3);
vNormal = zeros(1,3);
vEff = zeros(1,3); % vector from intersection of bond and normal plane to effector
vTarget = zeros(1,3); % vector from intersection of bond and normal plane to effector
vEffProj = zeros(1,3); % projection of vEff into plane normal to bond
vTargetProj = zeros(1,3); % projection of vTarget into plane normal to bond
vCross = zeros(1,3);

cos_theta = 0;
sin_theta = 0;
theta = 0;

j = 0;
max_iterations = 1000;
curr_dist = 0;
vEffProjLen = 0;
vTargetProjLen = 0;

while(j < max_iterations && mydist(new_molecule(loop_end,:), target) > delta)
	% curr_dist = mydist(new_molecule(loop_end,:), target)
	curr_index = floor(rand()* loop_len) + start_index;
	while(curr_index == last_index || rotable(curr_index) == 0)
		curr_index = floor(rand()* (loop_len + 1)) + start_index;
	end
	j = j + 1;
	
	v0 = new_molecule(curr_index - 1, :); % second atom defining dihedral
	vNormal = v0 - new_molecule(curr_index - 2, :); % direction of dihedral
	vNormal = vNormal ./ norm(vNormal);
	vEff = new_molecule(end_index, :) - v0; % position of effector relative v0
	vTarget = target - v0; % position of target relative v0
	vEffProj = vEff - vNormal .* dot(vNormal, vEff); % /norm(vNormal)^2);
	vTargetProj = vTarget - vNormal .* dot(vNormal, vTarget); % /norm(vNormal)^2);	
	
	% find angle between vEffProj and vTargetProj, change dihedral by that
	% amount
	
	vEffProjLen = norm(vEffProj);
	vTargetProjLen = norm(vTargetProj);
	vCross = cross(vEffProj, vTargetProj);
	
	if(vEffProjLen * vTargetProjLen > 0) 
		cos_theta = dot(vEffProj, vTargetProj)/(norm(vEffProj)*norm(vTargetProj));
		sin_theta = sign(dot(vNormal, vCross))*norm(vCross)/(norm(vEffProj)*norm(vTargetProj));
		theta = atan2(sin_theta, cos_theta);
	else
		theta = 0;
	end
	
	broken_dihedrals(curr_index) = broken_dihedrals(curr_index) - theta;
	new_molecule = con_molecule(real_lengths, real_angles, broken_dihedrals);
	last_index = curr_index;
end

j
curr_dist = mydist(new_molecule(loop_end,:), target)
new_molecule = con_molecule(real_lengths, real_angles, broken_dihedrals);

for i = end_index:1:num_atoms
	new_molecule(i,:) = old_coords(i,:);
end

[new_lengths, new_angles, new_dihedrals] = cartesian2internal(new_molecule);