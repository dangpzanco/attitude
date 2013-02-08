function q = Euler1212EP(e)

% Euler1212EP(E)
%
%	Q = Euler1212EP(E) translates the 121 Euler angle
%	vector E into the Euler parameter vector Q.
%

e1 = e(1)/2;
e2 = e(2)/2;
e3 = e(3)/2;

q(1) = cos(e2)*cos(e1+e3);
q(2) = cos(e2)*sin(e1+e3);
q(3) = sin(e2)*cos(e1-e3);
q(4) = sin(e2)*sin(e1-e3);
q=q';
