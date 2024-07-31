function y = shrink_Lp(x,mu,p)
    y = sign(x).*Lp(x,mu,p);     %这里的sign不仅仅可以求实数的正负号，还可以求复数的相位，更加简洁
end
function y=Lp(x,mu,p)
y=max(abs(x)-(1/mu).^(p-2).*abs(x).^(p-1),0); 
end
