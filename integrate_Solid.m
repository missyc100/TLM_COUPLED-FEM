function [Msparse,Ksparse] = integrate_Solid(NOCs,XofNs,Material_Elasstic,materils)
NEs=size(NOCs,1);
mk=8*8;
mm=8*8;
Ani=zeros(mk*NEs,1);Anj=zeros(mk*NEs,1);
Se=zeros(mk*NEs,1);

Me=zeros(mm*NEs,1);

%%%%%%初始矩阵值
for n = 1 :NEs
    SE=Ks(n,NOCs,XofNs,Material_Elasstic,materils) ;
    ME=Ms_J(n,NOCs,XofNs,Material_Elasstic,materils);
    % ME=Ms(n);
    %------形成总刚矩阵 K  M
    node = zeros(1,8);
    node(1:2:end) = 2*NOCs(n,:)-1;
    node(2:2:end) = 2*NOCs(n,:);
    [ni,nj]=meshgrid(node);
    nii=ni(:);
    njj=nj(:);

    Co_se=SE(:);

    Ani((n-1)*8*8+1:(n)*8*8,:)=nii;
    Anj((n-1)*8*8+1:(n)*8*8,:)=njj;
    Se((n-1)*8*8+1:(n)*8*8,:)=(Co_se);

    %     [im,jm,Co_me]=find(ME);
    %     Ami((n-1)*24*24+1:(n)*24*24,:) =node(im).';
    %     Amj((n-1)*24*24+1:(n)*24*24,:) =node(jm).';
    Me((n-1)*8*8+1:(n)*8*8,:)=(ME(:));

end
Ksparse=sparse(Anj,Ani,Se);
Msparse=sparse(Anj,Ani,Me);
end
