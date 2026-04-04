import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider

# --- Helper to load a file like rho_out.txt
def load_field(filename):
    data, zs = [], []
    with open(filename) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        z = float(lines[i].strip())
        zs.append(z)
        i += 1
        block = []
        while i < len(lines) and lines[i].strip():
            block.append(list(map(float, lines[i].split())))
            i += 1
        data.append(block)
    return np.array(zs), np.array(data)  # zs: (K,), data: (K, N, M)

# --- Body radius function
def r_b(z):
    # return np.tan(np.pi / 12.0) * z # коническое тело
    return np.tan(np.pi / 12.0) * np.sqrt(2*z - 1) # параболическое тело
    return np.tan(np.pi / 12.0) * (z / 20.0 + 0.95) # усеченный конус (конус + цилиндр)

# --- Load all fields
# zs, rho = load_field('rho_out.txt')
# _,   p   = load_field('p_out.txt')
# _,   u   = load_field('u_out.txt')
# _,   v   = load_field('v_out.txt')
# _,   w   = load_field('w_out.txt')

# # --- Load shock radius rs(z, theta)
# rs = []
# with open('r_s_out.txt') as f:
#     for line in f:
#         parts = list(map(float, line.split()))
#         rs.append(parts[1:])
# rs = np.array(rs)  # shape (K, M)

zs = np.loadtxt('z_out.txt')
K = zs.shape[0]
N = 200
M = 600
rho = np.loadtxt('rho_out.txt').reshape(K, N, M)
p = np.loadtxt('p_out.txt').reshape(K, N, M)
u = np.loadtxt('u_out.txt').reshape(K, N, M)
v = np.loadtxt('v_out.txt').reshape(K, N, M)
w = np.loadtxt('w_out.txt').reshape(K, N, M)
rs = np.loadtxt('r_s_out.txt')

# --- Dimensions and normalized xi-grid
K, N, M = rho.shape
xi = np.linspace(0, 1, N)

# --- Map user input to field arrays
tab = {'rho': rho, 'p': p, 'u': u, 'v': v, 'w': w}
while True:
    var = input("Введите, что рисовать (одно из ['rho','p','u','v','w']): ").strip()
    if var in tab:
        field = tab[var]
        break
    print("Неверный ввод.")

# --- Prepare half-circle slices for theta=0 and theta=pi
field0 = field[:, :,  0]
fieldp = field[:, :, -1]
r_phys0 = np.zeros_like(field0)
r_physp = np.zeros_like(fieldp)
for k, z in enumerate(zs):
    rb = r_b(z)
    r_phys0[k] = rb + xi * (rs[k, 0]  - rb)
    r_physp[k] = rb + xi * (rs[k, -1] - rb)

# --- Static 2x2 panels
fig, axs = plt.subplots(2,2, figsize=(12,8))
# xi-z, theta=0
pcm = axs[0,0].pcolormesh(xi, zs, field0, shading='auto', cmap='viridis')
axs[0,0].set_title(f"{var} (θ=0) в ξ–z")
axs[0,0].set_xlabel('ξ')
axs[0,0].set_ylabel('z')
fig.colorbar(pcm, ax=axs[0,0])
# r-z, theta=0
Zg = np.tile(zs[:,None], (1,N))
pcm = axs[0,1].pcolormesh(r_phys0, Zg, field0, shading='auto', cmap='viridis')
axs[0,1].set_title(f"{var} (θ=0) в r–z")
axs[0,1].set_xlabel('r')
axs[0,1].set_ylabel('z')
fig.colorbar(pcm, ax=axs[0,1])
# xi-z, theta=pi
pcm = axs[1,0].pcolormesh(xi, zs, fieldp, shading='auto', cmap='viridis')
axs[1,0].set_title(f"{var} (θ=π) в ξ–z")
axs[1,0].set_xlabel('ξ')
axs[1,0].set_ylabel('z')
fig.colorbar(pcm, ax=axs[1,0])
# r-z, theta=pi
pcm = axs[1,1].pcolormesh(r_physp, Zg, fieldp, shading='auto', cmap='viridis')
axs[1,1].set_title(f"{var} (θ=π) в r–z")
axs[1,1].set_xlabel('r')
axs[1,1].set_ylabel('z')
fig.colorbar(pcm, ax=axs[1,1])
plt.tight_layout()
plt.show()

# --- Reflect half-circle data into full [0,2π]
def make_full_circle(data_half, is_odd):
    # data_half shape (N, M)
    # original angles φ_i = i*(π/(M-1)), i=0..M-1
    # inner half without endpoints: data_half[:,1:-1]
    mid = data_half[:,1:-1]
    if is_odd:
        left = -mid[:, ::-1]
    else:
        left = mid[:, ::-1]
    return np.hstack([left, data_half])  # shape (N, 2M-2)

# --- Cross-section prep
def prepare_cross_section(z_idx, field, var):
    # half-circle phi
    phi = np.linspace(0, np.pi, M)
    # full φ from 0 to 2π without duplicate endpoints
    phi_full = np.hstack([phi[1:-1] + np.pi, phi])
    # compute r mesh
    z = zs[z_idx]
    rb = r_b(z)
    rsz = rs[z_idx]
    rs_full = np.hstack([rsz[1:-1], rsz])
    xi_m, phi_m = np.meshgrid(xi, phi_full, indexing='ij')
    r = rb + xi_m * (rs_full[np.newaxis,:] - rb)
    # rotate zero angle to top: θ=0 maps to φ=π/2
    phi_rot = np.pi/2 - phi_m
    x = r * np.cos(phi_rot)
    y = r * np.sin(phi_rot)
    # reflect field
    half = field[z_idx]  # (N, M)
    full = make_full_circle(half, var=='v')
    return x, y, full

# --- Interactive slider full circle
fig, ax = plt.subplots(figsize=(7,7))
plt.subplots_adjust(bottom=0.25)
x0, y0, d0 = prepare_cross_section(0, field, var)
cax = ax.pcolormesh(x0, y0, d0, shading='auto', cmap='plasma')
cbar = plt.colorbar(cax, ax=ax, label=var)
ax.set_aspect('equal')
ax.set_title(f'z={zs[0]:.3f}')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.axhline(0, color='k', linewidth=0.5)
ax.axvline(0, color='k', linewidth=0.5)
# slider
axsl = plt.axes([0.25,0.1,0.65,0.03])
sl = Slider(axsl, 'Z idx', 0, K-1, valinit=0, valstep=1)
def upd(val):
    i = int(sl.val)
    x,y,dd = prepare_cross_section(i, field, var)
    cax.set_array(dd.ravel())
    ax.set_title(f'z={zs[i]:.3f}')
    fig.canvas.draw_idle()
sl.on_changed(upd)
plt.tight_layout()
plt.show()

# --- Animation fixed colorbar
vmin,vmax = np.nanmin(field), np.nanmax(field)
fig, ax = plt.subplots(figsize=(7,7))
x0,y0,d0 = prepare_cross_section(0, field, var)
cax = ax.pcolormesh(x0,y0,d0,shading='auto',cmap='viridis',vmin=vmin,vmax=vmax)
plt.colorbar(cax, ax=ax, label=var)
ax.set_aspect('equal')
def anim(i):
    x,y,dd = prepare_cross_section(i, field, var)
    cax.set_array(dd.ravel())
    ax.set_title(f'z={zs[i]:.3f}')
    return cax,
ani = FuncAnimation(fig, anim, frames=K, interval=100, blit=False)
plt.show()

# --- Animation normalized [0,1]
fig, ax = plt.subplots(figsize=(7,7))
x0,y0,d0 = prepare_cross_section(0, field, var)
d0n = (d0-d0.min())/(d0.max()-d0.min() if d0.max()!=d0.min() else 1)
caxn = ax.pcolormesh(x0,y0,d0n,shading='auto',cmap='viridis',vmin=0,vmax=1)
plt.colorbar(caxn, ax=ax, label=f'{var} (0–1)')
ax.set_aspect('equal')
def animn(i):
    x,y,dd = prepare_cross_section(i, field, var)
    dnn = (dd-dd.min())/(dd.max()-dd.min() if dd.max()!=dd.min() else 1)
    caxn.set_array(dnn.ravel())
    ax.set_title(f'z={zs[i]:.3f}')
    return caxn,
ani2 = FuncAnimation(fig, animn, frames=K, interval=100, blit=False)
plt.tight_layout()
plt.show()