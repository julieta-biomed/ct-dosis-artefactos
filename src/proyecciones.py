import numpy as np
from skimage.transform import radon, iradon
from fantoma import abdomen, a_mu

N=256; C=N/2; LES_R=0.06
img=abdomen(N, lesion_hu=95, lesion_r=LES_R); mu=a_mu(img)
y,x=np.indices((N,N)); cuerpo=a_mu(img)>0
roi_les  = ((y-C*0.86)**2+(x-C*0.62)**2) < (LES_R*C*0.65)**2
roi_fondo= ((y-C*0.92)**2+(x-C*0.80)**2) < (0.05*C)**2

def sim(n_proy, I0, semilla=0, filtro='hann'):
    ang=np.linspace(0,180,n_proy,endpoint=False)
    s=radon(mu, theta=ang, circle=True)
    rng=np.random.default_rng(semilla)
    I=rng.poisson(np.maximum(I0*np.exp(-s),1e-9))
    s_r=-np.log(np.maximum(I,1)/I0)
    r=iradon(s_r, theta=ang, filter_name=filtro, circle=True, output_size=N)
    a,b=np.polyfit(r[cuerpo].ravel(), img[cuerpo].ravel(),1)
    rec=a*r+b
    cnr=abs(rec[roi_les].mean()-rec[roi_fondo].mean())/max(rec[roi_fondo].std(),1e-9)
    rmse=np.sqrt(((rec[cuerpo]-img[cuerpo])**2).mean())
    return rec, cnr, rmse, rec[roi_fondo].std()

print('=== A dosis TOTAL constante: reparto entre fotones y proyecciones ===')
print('   dosis total proporcional a  I0 x n_proyecciones = 2.16e8\n')
print(f'{"proyecciones":>13} {"I0/rayo":>10} | {"sigma (HU)":>11} {"RMSE":>8} {"CNR":>7}')
print('-'*56)
TOTAL=720*3e5
filas=[]
for n in (45,90,180,360,720):
    I0=TOTAL/n
    rec,c,rm,s = sim(n, I0)
    print(f'{n:>13} {I0:>10.1e} | {s:>11.1f} {rm:>8.1f} {c:>7.2f}')
    filas.append([n,I0,s,rm,c])
    np.save(f'rec_{n}.npy', rec)
np.save('reparto.npy', np.array(filas))

print('\n=== A igual numero de fotones por rayo: solo varian las vistas ===\n')
print(f'{"proyecciones":>13} | {"sigma (HU)":>11} {"RMSE":>8} {"CNR":>7}')
print('-'*46)
filas2=[]
for n in (30,45,90,180,360,720):
    rec,c,rm,s = sim(n, 3e4)
    print(f'{n:>13} | {s:>11.1f} {rm:>8.1f} {c:>7.2f}')
    filas2.append([n,s,rm,c])
    np.save(f'vistas_{n}.npy', rec)
np.save('vistas.npy', np.array(filas2))
