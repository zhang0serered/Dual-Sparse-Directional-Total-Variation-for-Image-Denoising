function y = shrink_Lp(x,mu,p)
    y = sign(x).*Lp(x,mu,p);     %�����sign������������ʵ���������ţ���������������λ�����Ӽ��
end
function y=Lp(x,mu,p)
y=max(abs(x)-(1/mu).^(p-2).*abs(x).^(p-1),0); 
end
