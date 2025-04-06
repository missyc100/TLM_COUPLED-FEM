function [AAL,ABL,AGL,AML,ADL] = Integrate_TLM(BNL,NOC_BL,Xof_B1,materils_boundary,Material_Elasstic)
%%%%%%≥ı ºæÿ’Û÷µ
AA1L= zeros(2*BNL,2*BNL);  AA2L= zeros(2*BNL,2*BNL); AA3L= zeros(2*BNL,2*BNL);  
AB1L= zeros(2*BNL,2*BNL);  AB2L= zeros(2*BNL,2*BNL); 
AGL= zeros(2*BNL,2*BNL);
AML= zeros(2*BNL,2*BNL); 
AD1L =zeros(2*BNL,2*BNL); AD2L =zeros(2*BNL,2*BNL);
BE=size(NOC_BL,1);

for n = 1:BE
    A1L=F_A1(n,Xof_B1,materils_boundary,Material_Elasstic);     
    % A2L=F_A2(n,BNL,NOC_BL,Xof_B1,materils_boundary,Material_Elasstic);    
    % A3L=F_A3(n,BNL,NOC_BL,Xof_B1,materils_boundary,Material_Elasstic);    
    % 
    B1L=F_B1(n,Xof_B1,materils_boundary,Material_Elasstic);     
    % B2L= F_B2(n,BNL,NOC_BL,Xof_B1,materils_boundary,Material_Elasstic);    
   
    GL=F_Ge(n,Xof_B1,materils_boundary,Material_Elasstic);         
    ML=F_Me(n,Xof_B1,materils_boundary,Material_Elasstic);    
    D1L=F_D1(n,Xof_B1,materils_boundary,Material_Elasstic);     
    % D2L=F_D2(n,BNL,NOC_BL,Xof_B1,materils_boundary,Material_Elasstic);

    node=zeros(1,2);
    degree = zeros(1,4);
    for j1 = 1 : 2
        node(j1) = NOC_BL(n,j1);
        degree(2*j1-1)=2*node(j1)-1;
        degree(2*j1)=2*node(j1);
    end
    AA1L(degree,degree) =  AA1L(degree,degree) + A1L;
    % AA2L(degree,degree) =  AA2L(degree,degree) + A2L;
    % AA3L(degree,degree) =  AA3L(degree,degree) + A3L;
    
    AB1L(degree,degree) =  AB1L(degree,degree) + B1L;
    % AB2L(degree,degree) =  AB2L(degree,degree) + B2L;
    
    AGL(degree,degree) =  AGL(degree,degree) + GL;
    
    AML(degree,degree) =  AML(degree,degree) + ML;
    
    AD1L(degree,degree) =  AD1L(degree,degree) + D1L;
    % AD2L(degree,degree) =  AD2L(degree,degree) + D2L;   
end
AAL=AA1L+AA2L+AA3L;    
ABL=AB1L+AB2L;    
ADL=AD1L+AD2L;          
end