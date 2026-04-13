import numpy as np
import matplotlib.pyplot as plt

def r_b(z): return np.tan(np.pi/12) * np.sqrt(2*z - 1)

def make_full_circle(data, odd):
    N,M = data.shape
    inner = data[:,1:-1]
    left = (-inner[:,::-1] if odd else inner[:,::-1])
    return np.hstack([left, data])

def prepare(zs, rs, xi, z_idx, field, var):
    M = field.shape[2]
    phi = np.linspace(0,np.pi, M)
    phi_full = np.hstack([phi[1:-1]+np.pi, phi])
    z = zs[z_idx]
    rb = r_b(z)
    rsz = rs[z_idx]
    rsf = np.hstack([rsz[1:-1], rsz])
    xi_m, phi_m = np.meshgrid(xi, phi_full, indexing='ij')
    r = rb + xi_m*(rsf - rb)
    phi_rot = np.pi/2 - phi_m
    x = r*np.cos(phi_rot)
    y = r*np.sin(phi_rot)
    half = field[z_idx]
    return x,y,make_full_circle(half, var=='v')

def create_frames_for_rho(z_values):
    data_root_path = '../output_results/'
    zs = np.loadtxt(data_root_path + 'z_out.txt')
    K = zs.size
    N,M = 500,1500
    rho = np.loadtxt(data_root_path + 'rho_out.txt').reshape(K,N,M)
    rs = np.loadtxt(data_root_path + 'r_s_out.txt')

    field = rho
    K, N, M = field.shape
    xi = np.linspace(0, 1, N)

    # --- Выбор конкретного z ---
    # Вариант 1: по индексу

    # Вариант 2: по значению (ближайший слой)
    # target_z = 2.5
    # z_idx = np.argmin(np.abs(zs - target_z))

    # Построение данных для выбранного слоя
        
    vmin_value = rho[:np.argmin(np.abs(zs - 1.8)),:,:].min()
    vmax_value = rho[:np.argmin(np.abs(zs - 1.8)),:,:].max()
    # vmin_value = 1.5

    for z in z_values:
        z_idx = np.argmin(np.abs(zs - z))
        x, y, d = prepare(zs, rs, xi, z_idx, field, 'rho')

        fig, ax = plt.subplots(figsize=(14, 6))

        pcm = ax.pcolormesh(x, y, d, shading='auto', cmap='plasma', vmin=vmin_value, vmax=vmax_value)
        cbar = fig.colorbar(pcm, ax=ax)
        cbar.set_label(r'$\mathrm{кг}/\mathrm{м}^3$', fontsize=20)
        cbar.ax.tick_params(labelsize=15)
        ax.set_aspect('equal')
        ax.set_xlabel('x', fontsize=20)
        ax.set_ylabel('y', fontsize=20)
        ax.tick_params(axis='both', labelsize=15)
        ax.set_title(f'Распределение плотности в сечении z = {zs[z_idx]:.2f}', fontsize=20)
        ax.set_ylim(bottom=0)

        plt.savefig(f"../pictures/for_diploma/wave_frames/z_{z}.png")

create_frames_for_rho([1.0, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.8])