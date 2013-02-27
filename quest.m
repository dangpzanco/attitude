function q = quest(Vb, Vi, w)
% 
% q = quest(Vb, Vi)
%
% Description: Given two+ body vectors, and their corresponding inertial
% vectors, calculate the DCM rotation from the inertial frame to the body
% frame, using least squares. This solves Wahba's problem using a modified
% Davenport's q-method in which the eigen values is solved for using a 
% faster newton-rhapson iteration.
%
% Inputs:   Vb = body vectors in form [vb1; vb2; vb3; etc..], where vbi is a
%                column vector
%           Vi = inertial vectors in form [vi1; vi2; vi3; etc..], where vii 
%                is a column vector (corresponding to vbi)
%            w = a column vector of weight corresponding to each body vector
%
% Outputs:   q = quaternion representing the rotation from the inertial 
%                frame to the body frame in form [q1i,q2j,q3k,q4], where q4
%                is the scalar portion of the quaternion
% 
% Notes:     Iteration not yet supported.
%

if nargin < 3, w = ones(1,size(Vb,2)); end

[Z,S,sigma] = getK(Vb, Vi, w);
lamda = sum(w);
p = ((lamda + sigma)*eye(3) - S)\Z;
q = 1/sqrt(1 + p'*p)*[p; 1];

end

function [Z,S,sigma] = getK(Vb, Vi, w)

B = bsxfun(@times,w,Vb)*Vi';
S = B + B';
Z = [B(2,3)-B(3,2);B(3,1)-B(1,3);B(1,2)-B(2,1)];
sigma = trace(B);
% K = [S-sigma*eye(3),Z;Z',sigma];

end
