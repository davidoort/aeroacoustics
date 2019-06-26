function [] = plotXFOILdata()
%PLOTXFOILDATA Temporary helper

data = readmatrix('xf-n0012-il-1000000');

alpha = deg2rad(data(:,1));
cl = data(:,2);
cd = data(:,3);

cl_line = 2*pi*alpha;
cd_quadratic = 0.007+ 1*alpha.^2; %(maybe find best fit instead up to Cdmax)


figure(1)
hold on
scatter(alpha,cd)
plot(alpha,cd_quadratic)
xlabel('$\alpha$ [rad]', 'Interpreter', 'latex')
ylabel('$C_d [-]$', 'Interpreter', 'latex')
legend('XFOIL data','Quadratic fit')
figure(2)
hold on
scatter(alpha,cl)
plot(alpha,cl_line)
xlabel('$\alpha$ [rad]', 'Interpreter', 'latex')
ylabel('$C_l [-]$', 'Interpreter', 'latex')
legend('XFOIL data','Linear fit')

end

