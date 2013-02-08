function e = EP2Euler231(q)

% EP2Euler231
%
%	E = EP2Euler231(Q) translates the Euler parameter vector
%	Q into the corresponding (2-3-1) Euler angle set.
%	

q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

e(1) = atan2(-2*(q1*q3-q0*q2), q0*q0+q1*q1-q2*q2-q3*q3);
e(2) = asin(2*(q1*q2+q0*q3));
e(3)= atan2(-2*(q2*q3-q0*q1),q0*q0-q1*q1+q2*q2-q3*q3);
e=e';
