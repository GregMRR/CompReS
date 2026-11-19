# Does Sum|w_hat truncated|^2/2 equal the grid mean 1/2<rho|u|^2> to 10 eps? (WP-11 forcing-only test set-up)
import numpy as np
N=32; x=2*np.pi*np.arange(N)/N; X,Y=np.meshgrid(x,x,indexing='ij')
rho=1+0.1*np.cos(X+Y)
k=np.fft.fftfreq(N,1/N); KX,KY=np.meshgrid(k,k,indexing='ij'); K2=KX**2+KY**2
keep=9*(KX**2+KY**2)<N**2
rng=np.random.default_rng(1)
uh=(rng.normal(size=(2,N,N))+1j*rng.normal(size=(2,N,N)))*((K2>=16)&(K2<=36))
uh[0]*=1;
u=np.real(np.fft.ifft2(uh,axes=(1,2)))   # real field in the band
uh=np.fft.fft2(u,axes=(1,2))/N**2
Egrid=0.5*np.mean(rho*(u**2).sum(0))
w=np.sqrt(rho)*u
wh=np.fft.fft2(w,axes=(1,2))/N**2
Espec=0.5*np.sum(np.abs(wh*keep)**2)
Efull=0.5*np.sum(np.abs(wh)**2)
print("grid",Egrid,"spec(full)",Efull,"spec(trunc)",Espec,"rel diff trunc",abs(Espec-Egrid)/Egrid,"eps",np.finfo(float).eps)
om=np.real(np.fft.ifft2(1j*(KX*np.fft.fft2(u[1])-KY*np.fft.fft2(u[0]))))
Og=0.5*np.mean(om**2/rho); s=om/np.sqrt(rho); sh=np.fft.fft2(s)/N**2
print("Omega rel diff trunc",abs(0.5*np.sum(np.abs(sh*keep)**2)-Og)/Og)
