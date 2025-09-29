# @title
import numpy as np
from math import *
import scipy.optimize as opt
import matplotlib.pyplot as plt
import scipy.integrate as int

# =====================================================
# Étude d’un puits quantique fini : valeurs propres et fonctions d’onde
# =====================================================

# -----------------------------
# Paramètres physiques du problème
# -----------------------------
N = 100               # Nombre d’états considérés (taille de la base)
a = 1e-9              # Demi-largeur de la zone du puits [m]
b = 3e-9              # Demi-largeur totale du domaine d’étude [m]
V0 = 0.5              # Profondeur du puits [eV]
me = 9.31e-31         # Masse de l’électron [kg] (approximation du problème)
hbar_SI = 1.054571628e-34 # Constante de Planck réduite [J·s]
eV_to_J = 1.602176634e-19 # Conversion eV → Joule

# -----------------------------
# 1) Fonction d’état de base (états propres du puits infini)
# -----------------------------
def phi(n, z, b):
    """
    Fonction d’onde de l’état n dans un puits infini centré en 0 et de demi-largeur b.
    Normalisée sur [-b, b].
    """
    return np.sqrt(1.0/b) * np.sin((n * np.pi / (2.0*b)) * (z + b))

# -----------------------------
# 2) Potentiel du puits fini
# -----------------------------
def Vz(z, a, b, V0):
    """
    Définit le profil du potentiel :
      -V0  dans la zone |z| < a      (fond du puits)
       0   dans la barrière a < |z| < b
     ~∞   en dehors |z| ≥ b (approximation numérique par une grande valeur)
    """
    V = np.full_like(z, 1e9)  # Très grande valeur pour simuler les bords infinis
    mask_inside_well = np.abs(z) < a
    mask_inside_barrier = (np.abs(z) >= a) & (np.abs(z) < b)
    V[mask_inside_well] = -V0
    V[mask_inside_barrier] = 0.0
    return V

# -----------------------------
# 3) Élément de matrice du potentiel
# -----------------------------
def pot_2(z, n, m):
    """
    Élément de matrice <n|V|m> avec les fonctions de base phi.
    """
    return phi(n, z, b) * Vz(z, a, b, V0) * phi(m, z, b)

def pot_matrice_element(n, m):
    """
    Calcule l’intégrale <n|V|m> par quadrature de Gauss (quad).
    """
    val, _ = int.quad(pot_2, -b, b, args=(n, m), epsabs=1e-9, epsrel=1e-9)
    return val

# -----------------------------
# 4) Énergie cinétique dans la base infinie
# -----------------------------
def energies_p(n):
    """
    Énergie cinétique (en eV) de l’état n d’un puits infini de largeur 2b.
    """
    T_joule = ((n * hbar_SI * np.pi)**2) / (8.0 * me * (b**2))
    return T_joule / eV_to_J

# -----------------------------
# 5) Construction de l’Hamiltonien
# -----------------------------
def H2(N):
    """
    Construit la matrice Hamiltonienne de taille NxN :
      H_nm = <n|T|m> + <n|V|m>
    où T est diagonal dans la base choisie.
    """
    H = np.zeros((N, N))
    for n in range(1, N+1):
        for m in range(1, N+1):
            T = energies_p(n) if n == m else 0.0
            V_nm = pot_matrice_element(n, m)
            H[n-1, m-1] = T + V_nm
    return H

# -----------------------------
# 6) Diagonalisation de H
# -----------------------------
H = H2(N)
vals, vecs = np.linalg.eig(H)

# Tri des valeurs propres et des vecteurs associés
idx_sort = np.argsort(vals)
vals_c = vals[idx_sort]
vecs_c = vecs[:, idx_sort]

# -----------------------------
# Résultats numériques : valeurs propres
# -----------------------------
print("Les 15 plus basses valeurs propres (en eV) :")
print(vals_c[:15])

# Extraction des états liés (énergies négatives)
bound_idx = np.where(vals_c < 0)[0]
bound_energies = vals_c[bound_idx]
print("\nÉnergies propres < 0 :")
for i, E in enumerate(bound_energies):
    print(f" État {i+1}: E = {E:.6f} eV")
print(f"\nNombre d'états liés (E < 0) = {len(bound_energies)}")

# -----------------------------
# 7) Fonctions d’onde associées aux états propres
# -----------------------------
def valeurs_propres(N):
    """
    Calcule et retourne les valeurs propres et vecteurs propres pour H2(N).
    """
    H = H2(N)
    valeurs, vecteurs = np.linalg.eig(H)
    valeurs = np.round(valeurs, 3)
    return valeurs, vecteurs

valeurs, vecteurs = valeurs_propres(N)
valeurs_norm = np.sort(valeurs)

print(f"\nMatrice H pour N={N} :\n", H2(N))
print(f"\nValeurs propres de la matrice H pour N={N} :\n", valeurs)
print(f"\nVecteurs propres associés :\n", vecteurs)

# Définition de l’espace réel
x_vals = np.linspace(-b, b, 500)

def psi4(c, x):
    """
    Reconstruit la fonction d’onde ψ(x) = Σ c_n φ_n(x) à partir des coefficients c.
    """
    psi_val = np.zeros_like(x)
    for n in range(1, N+1):
        psi_val += c[n-1] * phi(n, x, b)
    return psi_val

# -----------------------------
# Tracé des trois premiers états liés
# -----------------------------
max_to_plot = min(3, len(bound_energies))

plt.figure(figsize=(15, 5))
for i in range(max_to_plot):
    plt.subplot(1, 3, i+1)
    c = vecs[:, bound_idx[i]]
    psi_vals = psi4(c, x_vals)

    # Normalisation
    norm = np.sqrt(np.trapz(np.abs(psi_vals)**2, x_vals))
    psi_vals /= norm

    # Ajustement de signe si nécessaire
    if i == 2:
        psi_vals *= -1

    plt.plot(x_vals * 1e9, psi_vals, lw=2,
             label=f"État {i+1}: E={bound_energies[i]:.4f} eV")
    plt.xlabel("z (nm)")
    plt.ylabel("ψ(z)")
    plt.title(f"État {i+1}")
    plt.grid(True)
    plt.legend()

plt.tight_layout()
plt.show()

# -----------------------------
# 8) Potentiel et niveaux d’énergie dans le puits
# -----------------------------
pot_vals = Vz(x_vals, a, b, V0)

plt.figure(figsize=(8,5))
# Limitation de l’affichage du potentiel
pot_vals_display = np.minimum(pot_vals, 1.0)
plt.plot(x_vals*1e9, pot_vals_display, 'k-', lw=2, label="Potentiel (limité)")

# Superposition des niveaux et fonctions d’onde
echelle_facteur = 0.05
for i in range(max_to_plot):
    E = bound_energies[i]
    c = vecs_c[:, bound_idx[i]]

    plt.hlines(E, -b*1e9, b*1e9, linestyles='--', color='gray')

    if i == 1:
        psi_vals = (-1)*psi4(c, x_vals)
    else:
        psi_vals = psi4(c, x_vals)

    # Normalisation
    norm = np.sqrt(np.trapz(np.abs(psi_vals)**2, x_vals))
    psi_vals /= norm

    # Mise à l’échelle pour affichage
    ondes_vals = E + (echelle_facteur * psi_vals)/30e3
    plt.plot(x_vals*1e9, ondes_vals, lw=2,
             label=f"État {i+1}: E={E:.4f} eV")

plt.xlabel("z (nm)")
plt.ylabel("Énergie (eV)")
plt.title("Fonctions d’onde et niveaux d’énergie")
plt.grid(True)
plt.legend()
plt.ylim(-0.6, 0.2)
plt.xlim(-b*1e9, b*1e9)
plt.show()
