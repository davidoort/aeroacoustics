function [] = plotXFOILdata()
%PLOTXFOILDATA Temporary helper

data = readmatrix('xf-n0012-il-1000000');

alpha = data(:,1);
cl = data(:,2);
cd = data(:,3);

cd_quadratic = 0.011+ 0.000175*alpha.^2; %(maybe find best fit instead up to Cdmax)

hold on
scatter(alpha,cd)
plot(alpha,cd_quadratic)

end

