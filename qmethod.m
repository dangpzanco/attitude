function q = qmethod(Vb, Vi, w)
% 
% q = qmethod(Vb, Vi)
%
% Description: Given two+ body vectors, and their corresponding inertial
% vectors, calculate the DCM rotation from the inertial frame to the body
% frame, using least squares. This solves Wahba's problem using Davenport's
% q-method.
%
% Inputs:  Vb = body vectors in form [vb1, vb2, vb3, etc..], where vbi is a
%               column vector
%          Vi = inertial vectors in form [vi1, vi2, vi3, etc..], where vii 
%               is a column vector (corresponding to vbi)
%           w = a column vector of weight corresponding to each body vector
%
% Outputs:  q = quaternion representing the rotation from the inertial 
%               frame to the body frame in form [q1i,q2j,q3k,q4], where q4
%               is the scalar portion of the quaternion
% 

if nargin < 3, w = ones(1,size(Vb,2)); end

K = getK(Vb, Vi, w);
q = eigK(K);

end


function q = eigK(K)

[V,D] = eig(K);
[~, index] = max(diag(D));
q = V(:,index);

end

function K = getK(Vb, Vi, w)

B = bsxfun(@times,w,Vb)*Vi';
S = B + B';
Z = [B(2,3)-B(3,2);B(3,1)-B(1,3);B(1,2)-B(2,1)];
sigma = trace(B);
K = [S-sigma*eye(3),Z;Z',sigma];

end
