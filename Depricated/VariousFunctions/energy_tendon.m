function energy = energy_tendon(Atendon,FMo,lTr,lTs,shift)
%ENERGY_TENDON
%    ENERGY = ENERGY_TENDON(ATENDON,FMO,LTR,LTS,SHIFT)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    16-Mar-2020 15:58:49

energy = FMo.*(lTr.*5.0-lTr.*shift.*2.0e1).*(-1.0./2.0e1)+(FMo.*lTs.*exp(Atendon.*((lTr+lTs)./lTs-1.99e2./2.0e2)))./(Atendon.*5.0);