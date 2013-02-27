function BI = triadmethod(b1, b2, i1, i2)
% 
% BI = triadmethod(b1, b2, i1, i2)
%
% Description: Given two body vectors, and their corresponding inertial
% vectors, calculate the DCM rotation from the inertial frame to the body
% frame.
%
% Inputs:  b1 = body vector 1 (most accurate)
%          b2 = body vector 2 (least accurate)
%          i1 = inertial vector 1 (most accurate)
%          i2 = inertial vector 2 (least accurate)
%
% Outputs: BI = DCM representing the rotation from the inertial frame to the
%               body frame
% 

% Allocate memory:
BT = zeros(3);
IT = zeros(3);

% Solve triad method:
BT(:,1) = b1;
IT(:,1) = i1;
BT(:,2) = cross(b1,b2)/norm(cross(b1,b2));
IT(:,2) = cross(i1,i2)/norm(cross(i1,i2));
BT(:,3) = cross(BT(:,1),BT(:,2));
IT(:,3) = cross(IT(:,1),IT(:,2));

% Calculate rotation matrix:
BI = BT*IT';

end