function rmsd(crd,natom,rmsdfile,alignedfile)
% parameters: 
%   crd: A crd formatfile.  First line must be a comment.
%   natom: Number of atoms in the molecule
%   rmsdfile:  File in which to write RMSD measurements
%   alignedfile:  File in which to write aligned structures (.crd)

   A=readcrd(crd,natom);

   refconfs=A(:,1)';
   confs=A';

   [n,m] = size(confs);
   k=m/3; 
   ind=1:m;

  d = zeros(n,1);

  % compare each conf with first conf
  refconf=reshape(refconfs(ind),3,k)';
  for i=1:n
     conf = reshape(confs(i,ind),3,k)';
     [R,conf] = rotate_align2(refconf,conf);
     c = transpose(conf);
     B(:,i) = c(:);
     dist = refconf-conf;
     d(i) = sqrt((1/k)*sum(dist(:).^2));
  end

  rmsdout = fopen(rmsdfile,'w');
  fprintf(rmsdout,'%f\n', [d]);
  fclose(rmsdout);

   %rmsd(refconf,confs)

  writecrd(B,alignedfile);
