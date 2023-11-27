% Function to input model equations
function dx = sys_eqs_1D(t,x,d_1)

% 1-D simplified system (Exp. decaying solution)
dx = [
    -d_1*x(1);
    ];
end