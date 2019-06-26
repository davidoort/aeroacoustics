function [] = diskPlot(r,psi,matrix,state,varargin)
%{
diskPlot finds the induced velocity of a (or multiple) disk element(s)
when given an external axial vel (which could also include the downwash of 
another rotor), a thrust coefficient and a prandtl tip loss function

Inputs:
    r - (matrix [-]) non-dimensional r matrix (similar to the one created using meshgrid)

    psi - (matrix [rad]) dimensional azimuth angle matrix

    matrix - (matrix [?]) matrix with the same size as r and psi to be
    plotted with color coding

Outputs:

    plot - no automated axes or title
    
Other m-files required: none

Other files required: none

Literature referenced: none

Assumptions: none
    
Ideas: none
    

Author: David Oort Alonso, B.Sc, Aerospace Engineering
TU Delft, Faculty of Aerospace Engineering
email address: d.oortalonso@student.tudelft.nl  
Website: https://github.com/davidoort/aeroacoustics
June 2019; Last revision: 12-June-2019
%}

%------------- BEGIN CODE --------------

res = length(r(1,:)); %higher than 1000 will take quite long

if all(all(isnan(matrix)))
    levels = 0;
else
    levels = linspace(min(min(matrix)),max(max(matrix)),res);
end
hold on
contourf(r.*cos(psi),r.*sin(psi),matrix,levels,'LineStyle','none')

if ~isempty(varargin)
    title(strcat(varargin{1}),'Interpreter','latex')
end

xlabel('r'); ylabel('r')
axis equal
set(gca,'FontSize',16)
colorbar

%% Add hub contour

th = 0:pi/50:2*pi;
xunit = r(1,1) * cos(th);
yunit = r(1,1) * sin(th);
plot(xunit,yunit)

dim1 = [.2 .5 .3 .3];
str1 = strcat('Axial vel: ', num2str(state.axial_vel),' m/s');

dim2 = [.2 .55 .3 .3];
str2 = strcat('Tangent vel: ', num2str(state.tangent_vel),' m/s');

dim3 = [.2 .60 .3 .3];
str3 = strcat('Collective: ', num2str(state.collective),' [deg]');

% annotation('textbox',dim1,'String',str1,'FitBoxToText','on','BackgroundColor','white');
% annotation('textbox',dim2,'String',str2,'FitBoxToText','on','BackgroundColor','white');
% annotation('textbox',dim3,'String',str3,'FitBoxToText','on','BackgroundColor','white');

%{
STUPID METHOD IF YOU ARE GIVEN r
try
    r=varargin{2};
    th = 0:pi/50:2*pi;
    xunit = r * cos(th);
    yunit = r * sin(th);
    plot(xunit,yunit)
catch
    warning('Hub size not provided. diskPlot(r,psi,matrix,title,hub_radial_fraction)')
end
%}
end

