function distance = mydist(p1, p2)
%
% returns the euclidean distance between two points in cartesian n-space
%
distance = sqrt(sum((p1-p2).^2));