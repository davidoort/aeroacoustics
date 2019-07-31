function [] = plotXFOILdata()
%PLOTXFOILDATA Temporary helper

close all

data = readmatrix('xf-n0012-il-1000000');

alpha = data(:,1);
cl = data(:,2);
cd = data(:,3);

cl_line = deg2rad(2*pi)*alpha;
cd_quadratic = 0.007+ deg2rad(deg2rad(1))*alpha.^2; %(maybe find best fit instead up to Cdmax)
cd_fourth = 0.007+ deg2rad(deg2rad(0.08))*alpha.^2+(0.0175)^4*7.5*alpha.^4; %(maybe find best fit instead up to Cdmax)


figure(1)
hold on
scatter(alpha,cd)
plot(alpha,cd_quadratic)
xlabel('$\alpha$ [rad]', 'Interpreter', 'latex')
ylabel('$C_d [-]$', 'Interpreter', 'latex')
legend('NACA 0012 XFOIL data','Quadratic fit')
set(gca,'Fontsize',16)
figure(2)
hold on
scatter(alpha,cl)
plot(alpha,cl_line)
xlabel('$\alpha$ [deg]', 'Interpreter', 'latex')
ylabel('$C_l [-]$', 'Interpreter', 'latex')
legend('NACA 0012 XFOIL data','Linear fit')
set(gca,'Fontsize',16)
figure(3)
hold on
scatter(alpha,cd)
plot(alpha,cd_fourth)
xlabel('$\alpha$ [deg]', 'Interpreter', 'latex')
ylabel('$C_d [-]$', 'Interpreter', 'latex')
legend('NACA 0012 XFOIL data','Fourth order polynomial fit')
set(gca,'Fontsize',16)



end

