function out = crank(ic, bc, dx, dt, a, b, T)
% here bc is 2 by 100 
    N = T/dt;
    M = (b-a)/dx;
    out = zeros(N,M-2);
    uc = ic(2:M-1); 
    out(1,:) = uc;
    for n = 2:N
        f = @(z) sys(z, uc, bc(:,n-1:n), dx, dt, a, b);
        uc = fsolve(f, zeros(M-2,1));
        %w1 = [bc(1,n),w1',bc(2,n)];
        out(n,:) = uc;
        
    end
    out = [bc(1,:)',out,bc(2,:)'];
end

function F = sys(u, uc, bc, dx, dt, a, b) % lambda = 1
% bc takes arguments xleft and xright
% input specifically at this time step
% uc is solution at previous time step
    M = (b-a)/dx;
    r = -dt/(1i*dx^2);
    F(1) = -r*u(2)+(1+2*r)*u(1)-r*bc(1,2)+(-r*uc(2)-(1-2*r)*uc(1)-r*bc(1,1)) ...
    -(-1i*dt)/4*(abs(uc(1))^2+abs(u(1))^2)*(u(1)+uc(1));
    % here u also goes from 1 to M-2 (actually representing 2 to M-1)
    for j = 2:M-3
        F(j) = -r*u(j+1)+(1+2*r)*u(j)-r*u(j-1)+(-r*uc(j+1)-(1-2*r)*uc(j)-r*uc(j-1)) ...
            -(-1i*dt)/4*(abs(uc(j))^2+abs(u(j))^2)*(u(j)+uc(j));
    end
    F(M-2) = -r*bc(2,2)+(1+2*r)*u(M-2)-r*u(M-3)+(-r*bc(2,1)-(1-2*r)*uc(M-2)-r*uc(M-3)) ...
        -(-1i*dt)/4*(abs(uc(M-2))^2+abs(u(M-2))^2)*(u(M-2)+uc(M-2));
end