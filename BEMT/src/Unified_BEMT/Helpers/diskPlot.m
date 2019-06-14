function [] = diskPlot(r,psi,matrix,varargin)
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

res = 100; %higher than 1000 will take quite long

if all(all(isnan(matrix)))
    levels = 0;
else
    levels = linspace(min(min(matrix)),max(max(matrix)),res);
end
hold on
contourf(r.*cos(psi),r.*sin(psi),matrix,levels,'LineStyle','none')

if ~isempty(varargin)
    title(strcat(varargin{1}))
end

xlabel('r'); ylabel('r')
colorbar

%Add hub contour
try
    r=varargin{2};
    th = 0:pi/50:2*pi;
    xunit = r * cos(th);
    yunit = r * sin(th);
    plot(xunit,yunit)
catch
    warning('Hub size not provided. diskPlot(r,psi,matrix,title,hub_radial_fraction)')
end

end

