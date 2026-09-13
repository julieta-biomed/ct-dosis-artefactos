import numpy as np
from skimage.transform import radon, iradon
from fantoma import abdomen, a_mu, ventana

N=256; C=N/2
ang_full=np.linspace(0,180,720,endpoint=False)
LES_R=0.06   # radio normalizado -> 7.7 px a N=256
img=abdomen(N, lesion_hu=95, lesion_r=LES_R); mu=a_mu(img)
sino=radon(mu, theta=ang_full, circle=True)

y,x=np.indices((N,N))
cuerpo = a_mu(img)>0
roi_les  = ((y-C*0.86)**2+(x-C*0.62)**2) < (LES_R*C*0.65)**2     # interior de la lesion
roi_fondo= ((y-C*0.92)**2+(x-C*0.80)**2) < (0.05*C)**2          # higado adyacente

def con_ruido(s, I0, semilla=0):
    rng=np.random.default_rng(semilla)
    I=rng.poisson(np.maximum(I0*np.exp(-s),1e-9))
    return -np.log(np.maximum(I,1)/I0)

def reconstruir(s, ang, filtro='hann'):
    r=iradon(s, theta=ang, filter_name=filtro, circle=True, output_size=N)
    a,b=np.polyfit(r[cuerpo].ravel(), img[cuerpo].ravel(),1)
    return a*r+b

def cnr(rec):
    """Contraste-ruido: (media lesion - media fondo) / sigma del fondo."""
    return abs(rec[roi_les].mean()-rec[roi_fondo].mean())/max(rec[roi_fondo].std(),1e-9)

print('Ruido frente a dosis (720 proyecciones, filtro Hann)\n')
print(f'{"I0 (fotones/rayo)":>18} | {"dosis rel.":>10} | {"sigma (HU)":>11} | {"CNR":>7}')
print('-'*56)
base=None; res=[]
for I0 in (3e5,1e5,3e4,1e4,3e3,1e3,3e2):
    rec=reconstruir(con_ruido(sino,I0), ang_full)
    s=rec[roi_fondo].std(); c=cnr(rec)
    if base is None: base=(I0,s)
    print(f'{I0:>18.0e} | {I0/3e5:>10.3f} | {s:>11.1f} | {c:>7.2f}')
    res.append([I0,s,c])
res=np.array(res)
# verificar sigma ~ 1/sqrt(I0)
A=np.polyfit(np.log(res[:,0]), np.log(res[:,1]),1)
print(f'\nAjuste: sigma ~ I0^({A[0]:.3f})   |   teoria: I0^(-0.5)   |   error {abs(A[0]+0.5):.3f}')
np.save('dosis.npy',res); np.save('sino.npy',sino); np.save('img.npy',img)
