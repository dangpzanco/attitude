function B = BmatEuler312(q)

% BmatEuler312(Q)
%
%	B = BmatEuler312(Q) returns the 3x3 matrix which relates the 
%	body angular velocity vector w to the derivative of
%	(3-1-2) Euler angle vector Q.  
%	
%		dQ/dt = [B(Q)] w
%	

s2 = sin(q(2));
c2 = cos(q(2));
s3 = sin(q(3));
c3 = cos(q(3));

B(1,1) = -s3;
B(1,2) = 0;
B(1,3) = c3;
B(2,1) = c2*c3;
B(2,2) = 0;
B(2,3) = c2*s3;
B(3,1) = s2*s3;
B(3,2) = c2;
B(3,3) = -s2*c3;
B = B/c2;
