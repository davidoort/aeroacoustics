function cl_alpha = getCl_alpha(rotor,r)

%[Cl,Cd] = get2Dcoeffs(rotor,atm,phi,chord,lambda,lambda_T,r,psi);

cl_alpha = rotor.aero.cl_alpha*ones(size(r));

end