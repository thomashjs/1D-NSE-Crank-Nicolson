addpath('/Users/thomas/Desktop/6035/Project')
% need to modify with the correct local path

% Setting up Global Parameters and Variables
a = -10;
b = 10;
T = 0.1;
N = 1000;
M = 200;
x = linspace(a,b,M);
t = linspace(0,T,N);

%% Analytic Solution
sol = pdepe(0, @NLS,@NLSic, @NLSbc,x,t);

%% Plotting Analytic Solution
figure
u = abs(sol);
surf(x,t,u)
xlabel('x')
ylabel('t')
zlabel('|u(x,t)|')

% Below is more exact solution although there is no difference up to
% machine precision
%[S, X] = meshgrid(t,x);
%surf(X,S,abs(uAnt(S,X)));

%% Boundary Condition and Initial Condition
bc = [2^0.5*exp(1i*t')*sech(a),2^0.5*exp(1i*t')*sech(b)']';
ic = 2^0.5 * exp(0.5i*x).*sech(x);

%% Solve with Crank Nicolson
dt = T/N;
dx = (b-a)/M;
ucrank = crank(ic, bc, dx, dt, a, b, T);

%% Plot Crank Nicolson Solution
figure;
surf(x,t,abs(ucrank))
xlabel('x')
ylabel('t')
zlabel('numerical |u(x,t)|')
title('Crank Nicolson Solution')
%% Plot together
figure;
subplot(1,2,1)
surf(x,t,u)
xlabel('x')
ylabel('t')
zlabel('|u(x,t)|')
title('Analytical Solution')

subplot(1,2,2)
surf(x,t,abs(ucrank))
xlabel('x')
ylabel('t')
zlabel('numerical |u(x,t)|')
title('Crank Nicolson Solution')
%% comparison at specific times
figure;
subplot(2,2,1)
plot(x, abs(ucrank(25,:)),x,u(25,:))
title('t = 0.025')
xlabel('x')
ylabel('|u|')
legend('numerical','analytical')

subplot(2,2,2)
plot(x, abs(ucrank(50,:)),x,u(50,:))
title('t = 0.05')
xlabel('x')
ylabel('|u|')
legend('numerical','analytical')

subplot(2,2,3)
plot(x, abs(ucrank(75,:)),x,u(75,:))
title('t = 0.075')
xlabel('x')
ylabel('|u|')
legend('numerical','analytical')

subplot(2,2,4)
plot(x, abs(ucrank(100,:)),x,u(100,:))
title('t = 0.1')
xlabel('x')
ylabel('|u|')
legend('numerical','analytical')

%% Total Spatial Discretization Error by time
err = sum(abs(u-ucrank),2);

figure;
plot(t, err)
xlabel('time')
ylabel('Sum Absolute Error')
title('Spatial Discretization Error')
%% Mass of Crank Nicolson
% Trapezoidal Rule for integral
M1 = abs(ucrank).^2;
M1(:,M) = 0;
M2 = abs(circshift(ucrank,-1,2)).^2; % maybe should shift -1
M2(:,M) = 0;
Mass = dx/2*sum(M1'+ M2')';

figure;
plot(t,Mass)

%% Use norm of u_i to estimate integral for Mass
% can also use right approximation instead of trapezoidal (uncomment
% below) although there is no significant difference
uMass = dx*sum(abs(u').^2)';
%uMass = dx*sum(abs(u').^2)'; % analytic solution mass

figure;
plot(t,uMass,t,Mass);
xlabel('time')
ylabel('mass')
axis([0 T 3.97 3.99])
legend('analytical','numerical')
title('Conservation of Mass')
%%
figure;
plot(t,uMass)
axis([0 T 3.97 3.99])

%% Energy Conservation
% energy for Crank Nicolson
dudx = (circshift(ucrank, -1, 2)-ucrank)/dx;
temp = circshift(dudx,-1,2);
temp(:,M) = 0;
dudx(:,M) = 0;

H1 = (sum(abs(dudx).^2+abs(temp).^2,2)*dx/2)/2;

tmp1 = circshift(ucrank,-1,2);
tmp1 = tmp1(:,1:M-1);
H2 = 1/4*dx/2*sum(abs(ucrank(:,1:M-1)).^4+abs(tmp1).^4,2);

H = H1+H2;

%% energy for analytic solution
uu = @(t,x) uAnt(t,x);
g = @(t) integral(@(x) abs(uu(t,x)).^4, a,b)/4;
gg = zeros(N,1);

for j = 1:N
    gg(j)=g(j*dt);
end

syms z;
syms s;
ff = int(abs(diff(2^0.5 * exp(1i*(0.5*z+s)).*sech(z),z))^2,z,a,b);

s = t;
af = subs(ff);

Hg = gg+1/2*af';
%%
figure;
plot(t,Hg,t,H)
xlabel('time')
ylabel('Hamiltonian energy')
axis([0 T 2 3])
legend('analytical','numerical')
title('Conservation of Energy')

%% Mass & Energy
figure;

subplot(1,2,1)
plot(t,uMass,t,Mass);
xlabel('time')
ylabel('mass')
axis([0 T 3.97 3.99])
legend('analytical','numerical')
title('Conservation of Mass')

subplot(1,2,2)
plot(t,Hg,t,H)
xlabel('time')
ylabel('Hamiltonian energy')
axis([0 T 2 3])
legend('analytical','numerical')
title('Conservation of Energy')