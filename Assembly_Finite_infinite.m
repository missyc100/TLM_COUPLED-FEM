function[Keq]=Assembly_Finite_infinite(Fre,MM,CC,KK,SBL,SBR,plus_BLL,plus_BRR)
  Keq=-Fre^2*MM+1i*Fre*CC+KK;
  Keq1=zeros(size(Keq,1),size(Keq,1));
  Keq1(plus_BLL,plus_BLL)=SBL;
  Keq1(plus_BRR,plus_BRR)=SBR;
  Keq=Keq1+Keq;
end

