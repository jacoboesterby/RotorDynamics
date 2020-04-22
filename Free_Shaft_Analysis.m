function Free_Shaft_Analysis(n)
Lambda = [4.73,7.8532,10.9956,14.1372];
E  = 2.0E11;
rho = 7800;
Ro = (5/2)/1000;
I = pi*(Ro^4)/4;
A = pi*(Ro^2);
L = 0.4350;
for i=1:n
    if i>4
        Lambda(i) = (2*i+1)*pi/2;
    end
    omega(i) = Lambda(i)^2*sqrt((E*I)/(rho*A))/L^2;
end
keyboard
phi =@(x,lambda) (cosh(lambda*(L-x)/L)+cos(lambda*(L-x)/L))-(cosh(lambda)-cos(lambda))/(sinh(lambda)-sin(lambda))*(sinh(lambda*(L-x)/L)+sin(lambda*(L-x)/L));
disp('Frequencies in rad/s: ')
omega
disp('Frequencies in Hz: ')
f = omega./(2*pi)
figure
hold on
s = linspace(0,L,100);
for i=1:n
    plot(s,phi(s,Lambda(i))./max(abs(phi(s,Lambda(i)))),'linewidth',3);
    leg{i} = sprintf('$phi_{%d}$',i);
end
legend(leg,'interpreter','latex')
    