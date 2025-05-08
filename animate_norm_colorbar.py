import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- Общие
def r_b(z): return np.tan(np.pi/12) * np.sqrt(2*z - 1)
def load_all():
    zs = np.loadtxt('z_out.txt')
    K=zs.size
    N,M=200,600
    rho = np.loadtxt('rho_out.txt').reshape(K,N,M)
    p,u,v,w = [np.loadtxt(f).reshape(K,N,M) for f in 
        ('p_out.txt','u_out.txt','v_out.txt','w_out.txt')]
    rs = np.loadtxt('r_s_out.txt')
    return zs, rho,p,u,v,w, rs

def make_full_circle(data, odd):
    mid = data[:,1:-1]
    left = -mid[:,::-1] if odd else mid[:,::-1]
    return np.hstack([left,data])

def prepare(zs,rs,xi,z_idx,field,var):
    M = field.shape[2]
    phi = np.linspace(0,np.pi,M)
    phi_full = np.hstack([phi[1:-1]+np.pi,phi])
    z = zs[z_idx]
    rb=r_b(z)
    rsz=rs[z_idx]
    rsf = np.hstack([rsz[1:-1],rsz])
    xi_m,phi_m = np.meshgrid(xi,phi_full, indexing='ij')
    r = rb + xi_m*(rsf-rb)
    phi_rot = np.pi/2-phi_m
    x=r*np.cos(phi_rot)
    y=r*np.sin(phi_rot)
    half=field[z_idx]
    return x,y,make_full_circle(half,var=='v')

def main():
    zs, rho, p, u, v, w, rs = load_all()
    tab = {'rho':rho, 'p':p, 'u':u, 'v':v, 'w':w}
    var = input("Введите (rho,p,u,v,w): ").strip()
    field = tab[var]
    K,N,M = field.shape
    xi = np.linspace(0,1,N)

    fig, ax = plt.subplots(figsize=(7,7))
    x0, y0, d0 = prepare(zs,rs,xi,0,field,var)
    # нормировка для первого кадра
    d0_min, d0_max = d0.min(), d0.max()
    d0n = (d0 - d0_min) / (d0_max - d0_min if d0_max!=d0_min else 1)
    pcm = ax.pcolormesh(x0, y0, d0n, shading='auto', cmap='viridis', vmin=0, vmax=1)
    plt.colorbar(pcm, ax=ax, label=f'{var} (0–1)')
    ax.set_aspect('equal')

    def animate(i):
        x, y, d = prepare(zs,rs,xi,i,field,var)
        # raw min/max
        dmin, dmax = d.min(), d.max()
        # normalize
        dn = (d - dmin) / (dmax - dmin if dmax!=dmin else 1)
        pcm.set_array(dn.ravel())
        # обновляем титул
        ax.set_title(
            fr"z={zs[i]:.3f}, ${var}_{{\min}}$={dmin:.3e}, ${var}_{{\max}}$={dmax:.3e}"
        )
        return pcm,

    FuncAnimation(fig, animate, frames=K, interval=100, blit=False)
    plt.show()

if __name__=='__main__':
    main()

