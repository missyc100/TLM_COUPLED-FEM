function  [LaR,FbR,jj4]=solve_eig_R(FeiR,ladR,re_degree,esp2)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%%%%%%%%%%%%%ÿһƵ���µ������������󡢺�����ֵ����
LaR=zeros(re_degree,re_degree);
Fb_R=zeros(2*re_degree,re_degree);
jj4=0;
for jj2=1:2*re_degree%%%%%%%�ҳ�ģ̬��������ֵС��0��Ӧ�ĵ���������
    if    ((abs(imag(ladR(jj2,jj2)))<esp2) && (real((ladR(jj2,jj2)))>esp2 ))...
                     || (imag(ladR(jj2,jj2))<-esp2)
        jj4=1+jj4;
        LaR(jj4,jj4)=(ladR(jj2,jj2));
        Fb_R(:,jj4)=FeiR(:,jj2);%%%%%%%%%%�ҳ�����ֵ�鲿С�������������
    end
end
FbR=Fb_R(1:re_degree,:);%
end

