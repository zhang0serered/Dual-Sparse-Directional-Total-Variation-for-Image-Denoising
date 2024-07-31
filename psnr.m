

function psnr = psnr(x,y)
% psnr - compute the Peack Signal to Noise Ratio, defined by :
%       PSNR(x,y) = 10*log10( max(max(x),max(y))^2 / mean(|x-y|.^2) ).
%   p = psnr(x,y);
x=double(x);
y=double(y);
if(max(x(:))>10)
    x=x/255;
end
if(max(y(:))>10)
    y=y/255;
end
psnr=10*log10(1/mean((x(:)-y(:)).^2)) ;
end

