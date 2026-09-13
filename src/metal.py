import numpy as np
from skimage.transform import radon, iradon
from fantoma import abdomen, a_mu

N=256; C=N/2; LES_R=0.06
ang=np.linspace(0,180,720,endpoint=False)
y,x=np.indices((N,N))

img_sin=abdomen(N, lesion_hu=95, lesion_r=LES_R)
img_con=abdomen(N, lesion_hu=95, lesion_r=LES_R, metal=True)
cuerpo=a_mu(img_sin)>0
metal = np.abs(img_con-15000)<1
lejos = cuerpo & ~metal & (np.hypot(y-C*1.30, x-C*0.78) > 0.35*C)   # lejos del implante

def recon(im, I0=1e5, semilla=0, saturar=None):
    mu=a_mu(im); s=radon(mu, theta=ang, circle=True)
    if saturar is not None:                       # el detector satura: fotones faltantes
        s=np.minimum(s, saturar)
    rng=np.random.default_rng(semilla)
    I=rng.poisson(np.maximum(I0*np.exp(-s),1e-9))
    s_r=-np.log(np.maximum(I,1)/I0)
    r=iradon(s_r, theta=ang, filter_name='hann', circle=True, output_size=N)
    # escalar excluyendo el metal: su valor extremo distorsionaria el ajuste
    m = cuerpo & ~metal
    a,b=np.polyfit(r[m].ravel(), img_sin[m].ravel(),1)
    return a*r+b, s

ref,_   = recon(img_sin)
con,s_m = recon(img_con)
print(f'Atenuacion maxima del sinograma sin metal: {radon(a_mu(img_sin),theta=ang,circle=True).max():.1f}')
print(f'Atenuacion maxima con protesis metalica  : {s_m.max():.1f}')
print(f'  -> I/I0 detras del metal = {np.exp(-s_m.max()):.2e}')
print(f'  -> con 1e5 fotones incidentes llegan {1e5*np.exp(-s_m.max()):.1f} fotones\n')

print(f'{"region":>34} | {"error medio (HU)":>17} | {"desv. (HU)":>11}')
print('-'*70)
for nom, m in [('todo el cuerpo (sin metal)', cuerpo),
               ('todo el cuerpo (con metal)', cuerpo & ~metal),
               ('lejos del implante (con metal)', lejos)]:
    r = ref if 'sin metal' in nom else con
    e = r[m]-img_sin[m]
    print(f'{nom:>34} | {np.abs(e).mean():>17.1f} | {e.std():>11.1f}')
np.save('ref.npy',ref); np.save('con.npy',con); np.save('img_con.npy',img_con)
np.save('img_sin.npy',img_sin)
