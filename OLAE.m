function [q,W,s,S,d] = OLAE(Vb, Vi, w)
% 
% q = OLAE(Vb, Vi)
%
% Description: Given two+ body vectors, and their corresponding inertial
% vectors, calculate the DCM rotation from the inertial frame to the body
% frame, using least squares. This solves Wahba's problem using a Optimum 
% Linear Attitude Estimation technique.
%
% Inputs:  Vb = body vectors in form [vb1, vb2, vb3, etc..], where vbi is a
%               column vector
%          Vi = inertial vectors in form [vi1, vi2, vi3, etc..], where vii 
%               is a column vector (corresponding to vbi)
%           w = a column vector of weight corresponding to each body vector
%
% Outputs:  q = a classical rodriguez parameter representing the rotation 
%               from the inertial frame to the body frame in the form 
%               [q1, q2, q3]
% 

N = size(Vb,2);
if nargin < 3, w = ones(1,N); end

% Create difference matrix:
d = Vb - Vi;
d = reshape(d,numel(d),1);

% Create summation matrix:
s = Vb + Vi;
S = zeros(3*N,3);
for ii=1:3:3*N-2
   S(ii:ii+2, 1:3) = tilde(s(:, (ii-1)/3 + 1));
end

% Create weights matrix:
W = zeros(3*N);
for ii=1:3:3*N-2
    for jj=1:3:3*N-2
        if(ii == jj)
            temp = eye(3)*w((ii-1)/3 + 1);
        else
            temp = zeros(3);
        end
        
        W(ii:ii+2,jj:jj+2) = temp;
    end
end

% Solve for the optimal CRP solution:
q = (S'*W*S)\(S'*W*d);

end