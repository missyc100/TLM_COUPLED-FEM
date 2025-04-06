function [Material_Elasstic] = Materials_matrix(rs1,vt1,cs1,Type_material)
Material_Elasstic=cell(Type_material,2);


for ii=1:Type_material
    cs=cs1(ii,1);vt=vt1(ii,1);rs=rs1(ii,1);
    G=cs^2*rs;
    G=G+0.008*1i*G; % 粘滞阻尼 (1+2i*ebsal)*(G or lamd)
    Et=G*2*(1+vt);
    Lamd=vt*Et/(1+vt)/(1-2*vt);
    cp=sqrt((Lamd+2*G)/rs);
    a1=Lamd+2*G;

    Elastic=[a1,Lamd,0;Lamd,a1,0;0,0,G];
    Material_Elasstic{ii,1}=rs;
    Material_Elasstic{ii,2}=Lamd;
    Material_Elasstic{ii,3}=G;
    Material_Elasstic{ii,4}=cs;
    Material_Elasstic{ii,5}=cp;
    Material_Elasstic{ii,6}=Elastic;
end


end

