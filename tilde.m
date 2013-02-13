function xt = tilde(x)
% Return the tilde matrix for a given 3x1 x vector.
    xt = [      0, -x(3),  x(2); 
             x(3),     0, -x(1); 
            -x(2),  x(1),     0  ];
end