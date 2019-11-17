function [] = plotXFOILdata(airfoil)
%PLOTXFOILDATA Temporary helper

close all

if strcmpi(airfoil,'NACA0012')
    
    data = readmatrix('xf-n0012-il-1000000'); 
    alpha = data(:,1);
    cl_line = deg2rad(2*pi)*alpha;
    cd_quadratic = 0.011+ deg2rad(deg2rad(1))*alpha.^2; %(maybe find best fit instead up to Cdmax)
    cd_fourth = 0.007+ deg2rad(deg2rad(0.08))*alpha.^2+(0.0175)^4*7.5*alpha.^4; %(maybe find best fit instead up to Cdmax)

elseif strcmpi(airfoil,'NACA16006')   
    
    data = readmatrix('xf-naca16006-il-1000000'); 
    alpha = data(:,1);
    cl_line = deg2rad(6.7)*alpha;
    cd_quadratic = 0.003+ deg2rad(deg2rad(3))*alpha.^2; %(maybe find best fit instead up to Cdmax)
    %cd_fourth = 0.007+ deg2rad(deg2rad(0.08))*alpha.^2+(0.0175)^4*7.5*alpha.^4; %(maybe find best fit instead up to Cdmax)
elseif strcmpi(airfoil,'ClarkY')
    data = readmatrix('xf-clarky-il-1000000-n5'); 
    alpha = data(:,1);
    cl_line = deg2rad(6.1)*alpha+0.38;
    cd_quadratic = 0.007+ deg2rad(deg2rad(0.8))*alpha.^2; %(maybe find best fit instead up to Cdmax)
    %cd_fourth = 0.007+ deg2rad(deg2rad(0.08))*alpha.^2+(0.0175)^4*7.5*alpha.^4; %(maybe find best fit instead up to Cdmax)

elseif strcmpi(airfoil,'NACA23015')
    data = readmatrix('xf-naca23015-il-1000000'); 
    alpha = data(:,1);
    cl_line = deg2rad(6.1)*alpha+0.1;
    cd_quadratic = 0.007+ deg2rad(deg2rad(0.8))*alpha.^2; %(maybe find best fit instead up to Cdmax)
    %cd_fourth = 0.007+ deg2rad(deg2rad(0.08))*alpha.^2+(0.0175)^4*7.5*alpha.^4; %(maybe find best fit instead up to Cdmax)
    
end

cl = data(:,2);
cd = data(:,3);

try
    figure(1)
    hold on
    scatter(alpha,cd)
    plot(alpha,cd_quadratic)
    xlabel('$\alpha$ [rad]', 'Interpreter', 'latex')
    ylabel('$C_d [-]$', 'Interpreter', 'latex')
    legend(strcat(num2str(airfoil),' XFOIL data'),'Quadratic fit')
    set(gca,'Fontsize',16)
catch
    close;
    disp('Unable to plot drag polar quadratic')
end

try
    figure(2)
    hold on
    scatter(alpha,cl)
    plot(alpha,cl_line)
    xlabel('$\alpha$ [deg]', 'Interpreter', 'latex')
    ylabel('$C_l [-]$', 'Interpreter', 'latex')
    legend(strcat(num2str(airfoil),' XFOIL data'),'Linear fit')
    set(gca,'Fontsize',16)
catch
    close;
    disp('Unable to plot lift polar')
end

try
    figure(3)
    hold on
    scatter(alpha,cd)
    plot(alpha,cd_fourth)
    xlabel('$\alpha$ [deg]', 'Interpreter', 'latex')
    ylabel('$C_d [-]$', 'Interpreter', 'latex')
    legend(strcat(num2str(airfoil),' XFOIL data'),'Fourth order polynomial fit')
    set(gca,'Fontsize',16)
catch
    close;
    disp('Unable to plot drag polar fourth order')
end


end

