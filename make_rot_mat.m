function [rot_mat] = make_rot_mat(r_angle, r_vector)

% goal: find rotation matrix for r_angle about r_vector
% method: rotate vector about x-axis until lies in x-z plane,
% rotate about y-axis until parallel with lies on z-axis
% perform rotation about z-axis, put everything back 
% where you found it.

% normalize r_vector
r_vector = r_vector ./ norm(r_vector);

% find length of projection onto y-z plane
d = sqrt(r_vector(2)^2 + r_vector(3)^2);

a = r_vector(1);
b = r_vector(2);
c = r_vector(3);

r1 = eye(4);
r1(2,2) = c/d;
r1(3,3) = c/d;
r5 = r1;
r1(3,2) = -b/d;
r1(2,3) = b/d;
r5(3,2) = b/d;
r5(2,3) = -b/d;

r2 = eye(4);
r2(1,1) = d;
r2(3,3) = d;
r4 = r2;
r2(3,1) = -a;
r2(1,3) = a;
r4(3,1) = a;
r4(1,3) = -a;

r3 = eye(4);
r3(1,1) = cos(r_angle);
r3(2,2) = cos(r_angle);
r3(2,1) = -sin(r_angle);
r3(1,2) = sin(r_angle);

rot_mat = r1 * r2 * r3 * r4 * r5;