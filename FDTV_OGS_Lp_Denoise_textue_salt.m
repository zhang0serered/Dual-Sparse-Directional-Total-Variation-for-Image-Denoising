function [f,k,PSNRp,SSIMp,gmsdp]= FDTV_OGS_Lp_Denoise_textue_salt(g,I,mu1,mu2,beta0,beta1,beta2,beta3,beta4,p,K1, K2)
gama=0.618;
g=double(g);
[M,N]=size(g);
f =zeros(M,N);
z0=zeros(M,N);
z1=zeros(M,N);
z2=zeros(M,N);
z3=zeros(M,N);
z4=zeros(M,N);

Dh=zeros(M,N);
Dv=zeros(M,N);
D45=zeros(M,N);
D135=zeros(M,N);


Dualz0=zeros(M,N);
Dualz1=zeros(M,N);
Dualz2=zeros(M,N);
Dualz3=zeros(M,N);
Dualz4=zeros(M,N);


err = 1;
tol =1e-4;
D1= (psf2otf([1,-1],[M,N])); %频域横向差分
D2 = (psf2otf([1;-1],[M,N])); %频域纵向差分
D3 = (psf2otf([0,1;-1 0],[M,N]));
D4 = (psf2otf([1 0;0 -1],[M,N]));
err=1;
lhs=beta0+beta1*abs(D1).^2+beta2*abs(D2).^2+beta3*abs(D3).^2+beta4*abs(D4).^2;
beta0g=beta0*g;
k=1;
while err>tol
 k
up = f;
Tmp0=z0-Dualz0;
Tmp1=z1-Dualz1;
Tmp2=z2-Dualz2;
Tmp3=z3-Dualz3;
Tmp4=z4-Dualz4;
rhs= fft2(beta0*(Tmp0)+beta0g)+beta1*((fft2(Tmp1)).*conj((D1))+(fft2(Tmp2)).*conj((D2))+(fft2(Tmp3)).*conj((D3))+(fft2(Tmp4)).*conj((D4))); 
F=rhs./(lhs);
f=real(ifft2(F));
PSNRp(k)=psnr(f,I);
SSIMp(k)=ssim(f,I);
gmsdp(k)=GMSD(f,I);
% PSNR(k)=psnr(f,I);
% SSIM(k)=ssim(f,I);
% gmsd(k)=GMSD(f,I);
%[Dh,Dv,D45,D135]=ForwardD(f);
Dh=real(ifft2(D1.*F));
Dv=real(ifft2(D2.*F));
D45=real(ifft2(D3.*F));
D135=real(ifft2(D4.*F));
%%=========================================================
z0=shrink_Lp(f+Dualz0-g,1/beta0,p);
z1=Mogshrink2(Dh+Dualz1,K1,K2, mu1/beta1,5);   
z2=Mogshrink2(Dv+Dualz2,K1,K2, mu2/beta2,5); 
z3=Mogshrink2(D45+Dualz3,K1,K2, mu1/beta3,5);   
z4=Mogshrink2(D135+Dualz4,K1,K2, mu2/beta4,5);


Dualz0=Dualz0+gama*beta0*(f-g-z0);
Dualz1=Dualz1+gama*beta1*(Dh-z1);
Dualz2=Dualz2+gama*beta2*(Dv-z2);
Dualz3=Dualz1+gama*beta3*(D45-z3);
Dualz4=Dualz2+gama*beta4*(D135-z4);
err = norm(up(:)-f(:))/norm(f(:));
%================超过最大迭代次数200，就停止=================
k=k+1;
if k>60
    break;    
end

end
end


function [a] = Mogshrink2(y, K1, K2, lam, Nit)
% [a, cost] = ogshrink2(y, K1, K2, lam, Nit);
% 2D overlapping group shrinkage (OGS)
% Author: Yingpin Chen
y=y(:);
n=length(y);
N=sqrt(n);
a = reshape(y,N,N); 
h1 = ones(K1,K2); % for convolution
for it = 1:Nit
r = sqrt(conv2(abs(a).^2,h1,'same'));
la=conv2(1./(r.*(1+0.001*r+0.001^2*r.*r)),h1 , 'same');
v = 1 + lam*la(:);
a = y./v;
end
a= reshape(a,N,N); 
end

function [Dux,Duy] = ForwardD(U)
Dux = [diff(U,1,2), U(:,1) - U(:,end)];
Duy = [diff(U,1,1); U(1,:) - U(end,:)];
end

function DtXY = Dive(X,Y)
DtXY = [X(:,end) - X(:, 1), -diff(X,1,2)];
DtXY = DtXY + [Y(end,:) - Y(1, :); -diff(Y,1,1)];
end


