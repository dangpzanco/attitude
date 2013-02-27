function C = cayleytransform(Q)
%
% Compute the cayley transform of the matrix Q.
%
I = eye(size(Q));
C = (I + Q)\(I - Q);

end