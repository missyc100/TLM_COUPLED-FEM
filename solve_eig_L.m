function  [LaL,FbL,jj4]=solve_eig_L(FeiL,ladL,re_degree,esp2)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%%%%%%%%%%%%%ÿһƵ���µ������������󡢺�����ֵ����
LaL=zeros(re_degree,re_degree);
Fb_L=zeros(2*re_degree,re_degree);
jj4=0;
for jj2=1:2*re_degree%%%%%%%�ҳ�ģ̬��������ֵС��0��Ӧ�ĵ���������
    if  ((abs(imag(ladL(jj2,jj2)))<esp2) && (real((ladL(jj2,jj2)))<-esp2 ))...
                     || (imag(ladL(jj2,jj2))>esp2)
        jj4=1+jj4;
        LaL(jj4,jj4)=(ladL(jj2,jj2));
        Fb_L(:,jj4)=FeiL(:,jj2);%%%%%%%%%%�ҳ�����ֵ�鲿С�������������
    end
end
FbL=Fb_L(1:re_degree,:);%
end

