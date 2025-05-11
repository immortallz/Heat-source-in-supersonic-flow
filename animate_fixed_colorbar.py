import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# --- Общие
def r_b(z): return np.tan(np.pi/12) * np.sqrt(2*z - 1)
def load_all():
    zs = np.loadtxt('z_out.txt')
    K=zs.size
    N,M=500,1500
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
    zs,rho,p,u,v,w,rs = load_all()
    tab={'rho':rho,'p':p,'u':u,'v':v,'w':w}
    var=input("Введите (rho,p,u,v,w): ").strip()
    field=tab[var]
    K,N,M = field.shape
    xi=np.linspace(0,1,N)

    # фикс колорбар
    # vmin=np.nanmin(field)
    # vmax=np.nanmax(field)

    fig, ax = plt.subplots(figsize=(7,7))
    x,y,d = prepare(zs,rs,xi,0,field,var)
    pcm = ax.pcolormesh(x,y,d,shading='auto',cmap='viridis',vmin=1.6,vmax=2.1)
    plt.colorbar(pcm, ax=ax, label=var)
    ax.set_aspect('equal')

    def animate(i):
        x,y,d = prepare(zs,rs,xi,i,field,var)
        pcm.set_array(d.ravel())
        ax.set_title(f'z={zs[i]:.3f}')
        return pcm,

    anim = FuncAnimation(fig, animate, frames=K, interval=100, blit=False)

    # --- сохраняем в GIF ---
    writer = PillowWriter(fps=10)              # или другой FPS на ваше усмотрение
    anim.save("animation.gif", writer=writer)  # имя файла и писатель
    # plt.show()

if __name__=='__main__':
    main()
