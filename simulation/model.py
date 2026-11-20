"""
Biophysical Continuum Model of Cylindrical Microvascular Pharmacological Sanctuaries
and Inflammation-Driven Clonal Recurrence
Author: Luca Celebrano
Repository: Lukecele/inflammatory-onco-physics

Solves the cylindrical modified Bessel K0 reaction-diffusion transport equation
coupled to the hereditary Volterra immune exhaustion functional.
"""

from dataclasses import dataclass
from typing import Dict, Tuple
import numpy as np
import scipy.special as sp

@dataclass
class ModelParameters:
    # Cylindrical microvascular domain (Krogh tissue cylinder)
    r_v: float = 10.0          # Capillary radius (microns)
    r_ox: float = 200.0        # Viable tumor cord outer boundary / oxygen diffusion limit (microns)
    Nr: int = 60               # Radial grid points

    # Tumor proliferation & carrying capacity
    r_S: float = 0.08          # Sensitive clone proliferation rate (month^-1)
    r_R: float = 0.015         # Persister (DTP) clone residual proliferation rate (month^-1)
    K: float = 1.0             # Normalized local carrying capacity
    
    # Plastic DTP transitions
    alpha: float = 0.05        # Susceptibility of sensitive cells to DTP entry via inflammation (Phi^-1 month^-1)
    beta: float = 0.01         # Spontaneous awakening rate from DTP to sensitive (month^-1)

    # Biophysical Transport & Pharmacodynamics (Bessel K0 decay)
    # Osimertinib verified human parameters: D_eff = 59.29 um^2/s, k = 1.51e-4 s^-1
    # lambda = sqrt(D_eff / k) = 627.6 um
    lambda_C: float = 627.6    # Effective drug diffusion length (microns)
    C_plasma_on: float = 28.6  # Steady state free plasma concentration (nM)
    IC_50: float = 10.0        # Half-maximal inhibitory concentration (nM)
    E_max: float = 0.25        # Maximum drug-induced kill rate (month^-1)
    
    # Inflammatory Field (Phi)
    D_Phi: float = 80.0        # Cytokine diffusion coefficient (um^2 / month)
    k_Phi: float = 0.4         # Inflammatory clearance rate (month^-1)
    sigma_base: float = 0.1    # Baseline chronic systemic/tissue stress (NLR, CRP proxy)
    eta: float = 0.3           # Debris inflammatory release per lysed cell

    # Immunosurveillance & Volterra History-Dependent Exhaustion
    E_0: float = 0.04          # Baseline cytotoxic immune efficacy (month^-1)
    gamma: float = 0.012       # Exhaustion coefficient per cumulative inflammatory unit

class InflammatoryOncoContinuum:
    def __init__(self, params: ModelParameters):
        self.p = params
        self.r = np.linspace(self.p.r_v, self.p.r_ox, self.p.Nr)
        self.dr = self.r[1] - self.r[0]
        
        # Calculate critical sanctuary radius analytically via Bessel K0
        # C(r_crit) = C_0 * K0(r_crit / lambda) / K0(r_v / lambda) = IC_50
        # For lambda = 627.6, r_v = 10.0, C_0 = 28.6, IC_50 = 10.0 -> r_crit = 165.89 um
        self.r_crit = self._compute_r_crit()
        self.sanctuary_volume_fraction = (self.p.r_ox**2 - self.r_crit**2) / (self.p.r_ox**2 - self.p.r_v**2)

    def _compute_r_crit(self) -> float:
        """Find radial distance where C(r) = IC_50."""
        target_ratio = self.p.IC_50 / self.p.C_plasma_on
        k0_rv = sp.kn(0, self.p.r_v / self.p.lambda_C)
        
        # Binary search for r_crit in [r_v, r_ox]
        low, high = self.p.r_v, self.p.r_ox
        for _ in range(50):
            mid = (low + high) / 2.0
            ratio = sp.kn(0, mid / self.p.lambda_C) / k0_rv
            if ratio > target_ratio:
                low = mid
            else:
                high = mid
        return (low + high) / 2.0

    def radial_laplacian(self, u: np.ndarray) -> np.ndarray:
        """
        Cylindrical radial Laplacian:
        Laplacian_r(u) = d^2u/dr^2 + (1/r) * du/dr
        """
        d2u = np.zeros_like(u)
        du_dr = np.zeros_like(u)
        
        # Interior points central differences
        d2u[1:-1] = (u[2:] - 2*u[1:-1] + u[:-2]) / (self.dr**2)
        du_dr[1:-1] = (u[2:] - u[:-2]) / (2 * self.dr)
        
        # Outer boundary r = r_ox: zero flux Neumann (du/dr = 0)
        d2u[-1] = 2 * (u[-2] - u[-1]) / (self.dr**2)
        du_dr[-1] = 0.0
        
        # Inner boundary r = r_v: zero flux for cytokines
        d2u[0] = 2 * (u[1] - u[0]) / (self.dr**2)
        du_dr[0] = 0.0
        
        return d2u + (1.0 / self.r) * du_dr

    def steady_state_drug_profile(self, C_plasma: float) -> np.ndarray:
        """
        Analytical cylindrical Bessel K0 decay from microvessel wall:
        C(r) = C_plasma * K0(r / lambda_C) / K0(r_v / lambda_C)
        """
        if C_plasma <= 1e-9:
            return np.zeros_like(self.r)
        k0_rv = sp.kn(0, self.p.r_v / self.p.lambda_C)
        return C_plasma * sp.kn(0, self.r / self.p.lambda_C) / k0_rv

    def simulate(
        self,
        t_span_months: float = 60.0,
        t_therapy_stop: float = 36.0,
        dt_months: float = 0.01,
        initial_tumor_burden: float = 0.35
    ) -> Dict[str, np.ndarray]:
        p = self.p
        nt = int(t_span_months / dt_months) + 1
        time_points = np.linspace(0, t_span_months, nt)

        # Initialize spatial fields
        N_S = np.ones(p.Nr) * initial_tumor_burden * 0.95
        N_R = np.ones(p.Nr) * initial_tumor_burden * 0.05
        Phi = np.ones(p.Nr) * p.sigma_base / p.k_Phi
        Psi = np.zeros(p.Nr)  # Cumulative inflammatory exposure integral (Volterra memory)
        
        save_interval = int(0.5 / dt_months)

        history_time = []
        history_N_S_mean = []
        history_N_R_mean = []
        history_N_tot_mean = []
        history_Phi_mean = []
        history_Imm_mean = []
        spatial_snapshots = {}

        for step in range(nt):
            t = time_points[step]

            # Drug delivery regime
            is_active_therapy = (t <= t_therapy_stop)
            C_plasma = p.C_plasma_on if is_active_therapy else 0.0

            # Analytical Cylindrical Bessel K0 decay
            C = self.steady_state_drug_profile(C_plasma)

            # Local pharmacodynamics
            H_C = p.E_max * C / (C + p.IC_50)
            lysis_rate = H_C * N_S

            # Local immune surveillance capacity (Volterra memory decay)
            # E_imm(t) = E_0 * exp(-gamma * int_0^t Phi(r, tau) dtau)
            E_imm = p.E_0 * np.exp(-p.gamma * Psi)

            # Inflammatory field cylindrical PDE update
            lap_Phi = self.radial_laplacian(Phi)
            dPhi_dt = p.D_Phi * lap_Phi - p.k_Phi * Phi + p.sigma_base + p.eta * lysis_rate
            Phi += dPhi_dt * dt_months
            Phi = np.maximum(Phi, 0.0)

            # Update cumulative exposure integral
            Psi += Phi * dt_months

            # Tumor population update
            total_N = N_S + N_R
            saturation = np.maximum(1.0 - total_N / p.K, 0.0)

            # DTP transition fluxes
            flux_to_persister = p.alpha * Phi * N_S
            flux_to_sensitive = p.beta * N_R

            dN_S_dt = (p.r_S * N_S * saturation 
                       - lysis_rate 
                       - E_imm * N_S 
                       - flux_to_persister 
                       + flux_to_sensitive)

            dN_R_dt = (p.r_R * N_R * saturation 
                       - E_imm * N_R 
                       + flux_to_persister 
                       - flux_to_sensitive)

            N_S += dN_S_dt * dt_months
            N_R += dN_R_dt * dt_months

            N_S = np.maximum(N_S, 0.0)
            N_R = np.maximum(N_R, 0.0)

            # Record metrics
            if step % save_interval == 0 or step == nt - 1:
                history_time.append(t)
                history_N_S_mean.append(np.average(N_S, weights=self.r))
                history_N_R_mean.append(np.average(N_R, weights=self.r))
                history_N_tot_mean.append(np.average(N_S + N_R, weights=self.r))
                history_Phi_mean.append(np.average(Phi, weights=self.r))
                history_Imm_mean.append(np.average(E_imm, weights=self.r))

            # Store snapshots
            for target_t in [0.0, 12.0, 36.0, 48.0]:
                if abs(t - target_t) < (dt_months / 2):
                    spatial_snapshots[f"month_{int(target_t)}"] = {
                        "r": self.r.copy(),
                        "N_S": N_S.copy(),
                        "N_R": N_R.copy(),
                        "C": C.copy(),
                        "Phi": Phi.copy(),
                        "E_imm": E_imm.copy()
                    }

        return {
            "time": np.array(history_time),
            "N_S": np.array(history_N_S_mean),
            "N_R": np.array(history_N_R_mean),
            "N_tot": np.array(history_N_tot_mean),
            "Phi": np.array(history_Phi_mean),
            "E_imm": np.array(history_Imm_mean),
            "snapshots": spatial_snapshots,
            "r_crit": self.r_crit,
            "sanctuary_fraction": self.sanctuary_volume_fraction
        }
