function [Q,B] = rotate_align2(A,B)
% Golub & Van Loan, 12.4.1 Rotation of Subspaces
% each row is one `measurement' or data point
% Q is an orthogonal matrix that minimizes || A - BQ ||
% Here A' and B' are equal to A and B with the means subtracted.
  
  m1 = mean(A);
  m2 = mean(B);
  A = A-m1(ones(size(A,1),1),:);
  B2 = B-m2(ones(size(B,1),1),:);
  C=B2'*A;
  [U,S,V]=svd(C);
  Q=U*V';

  % compute aligned B
  B = B2*Q+m1(ones(size(A,1),1),:);
  
