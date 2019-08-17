function Cd = getCd(method,rotor,atm,phi,chord,lambda,lambda_T,r,psi)

alpha = rotor.pitch-phi; %rad
alpha_0 = rotor(1).aero.alpha_0; %rad

if strcmpi(method,'leishman')
    
    Cd = rotor.aero.Cd0+rotor.aero.D1*(alpha-alpha_0)+rotor.aero.D2*(alpha-alpha_0).^2; %quadratic equation
    
elseif strcmpi(method,'vitleish')
    
    [~,Cd] = get2Dcoeffs(rotor,atm,phi,chord,lambda,lambda_T,r,psi);
    %Cd = rotor.aero.Cd0+rotor.aero.D1*(alpha-alpha_0)+rotor.aero.D2*(alpha-alpha_0).^2; %quadratic equation

end

end