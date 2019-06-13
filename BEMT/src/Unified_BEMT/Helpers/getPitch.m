function pitch = getPitch(r_vec,type,twist,root_pitch)                       
%GETPITCH Summary of this function goes here
%   Detailed explanation goes here


if strcmpi(type,'linear')
    pitch = linspace(root_pitch,root_pitch-twist,length(r_vec))-root_pitch;
elseif strcmpi(type,'ideal')
    kappa = r_vec(1)*root_pitch;
    pitch = ones(size(r_vec))*kappa./r_vec-root_pitch;
end



end

