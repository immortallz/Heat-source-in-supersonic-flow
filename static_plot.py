import numpy as np
import matplotlib.pyplot as plt

# --- Общие функции и загрузка
def r_b(z):
    return np.tan(np.pi/12) * np.sqrt(2*z - 1)

def load_all():
    zs = np.loadtxt('z_out.txt')
    K = zs.size
    N, M = 200, 600
    rho = np.loadtxt('rho_out.txt').reshape(K,N,M)
    p   = np.loadtxt('p_out.txt').reshape(K,N,M)
    u   = np.loadtxt('u_out.txt').reshape(K,N,M)
    v   = np.loadtxt('v_out.txt').reshape(K,N,M)
    w   = np.loadtxt('w_out.txt').reshape(K,N,M)
    rs  = np.loadtxt('r_s_out.txt')   # shape (K,M)
    return zs, rho, p, u, v, w, rs

def main():
    zs, rho, p, u, v, w, rs = load_all()
    tab = {'rho': rho, 'p': p, 'u': u, 'v': v, 'w': w}
    var = input("Введите (rho,p,u,v,w): ").strip()
    field = tab[var]
    K,N,M = field.shape
    xi = np.linspace(0,1,N)

    # вырезы θ=0 и θ=π
    f0 = field[:,:,0]
    fp = field[:,:,-1]
    r0 = np.zeros_like(f0)
    rp = np.zeros_like(fp)
    for i,z in enumerate(zs):
        rb = r_b(z)
        r0[i] = rb + xi*(rs[i,0]-rb)
        rp[i] = rb + xi*(rs[i,-1]-rb)

    # 2×2
    fig, axs = plt.subplots(2,2, figsize=(12,8))
    pcm = axs[0,0].pcolormesh(xi, zs, f0,   shading='auto', cmap='viridis')
    axs[0,0].set_title(f"{var} (θ=0) ξ–z")
    axs[0,0].set(xlabel='ξ', ylabel='z')
    fig.colorbar(pcm, ax=axs[0,0])

    Zg = np.tile(zs[:,None],(1,N))
    pcm = axs[0,1].pcolormesh(r0, Zg, f0, shading='auto', cmap='viridis')
    axs[0,1].set_title(f"{var} (θ=0) r–z")
    axs[0,1].set(xlabel='r', ylabel='z')
    fig.colorbar(pcm, ax=axs[0,1])

    pcm = axs[1,0].pcolormesh(xi, zs, fp,   shading='auto', cmap='viridis')
    axs[1,0].set_title(f"{var} (θ=π) ξ–z")
    axs[1,0].set(xlabel='ξ', ylabel='z')
    fig.colorbar(pcm, ax=axs[1,0])

    pcm = axs[1,1].pcolormesh(rp, Zg, fp, shading='auto', cmap='viridis')
    axs[1,1].set_title(f"{var} (θ=π) r–z")
    axs[1,1].set(xlabel='r', ylabel='z')
    fig.colorbar(pcm, ax=axs[1,1])

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
