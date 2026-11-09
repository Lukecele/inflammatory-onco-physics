"""
Biophysical Continuum Model of Inflammation-Driven Clonal Persistence and Recurrence
Authors: Luca Celebrano, Antigravity AI
Repository: /home/luca/inflammatory-onco-physics
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, Tuple

@dataclass
class ModelParameters:
    # Spatial domain
    L: float = 100.0          # Tissue depth from capillary to sanctuary (microns)
    Nx: int = 51              # Spatial grid points

    # Tumor proliferation & carrying capacity
    r_S: float = 0.08         # Sensitive clone proliferation rate (month^-1)
    r_R: float = 0.015        # Persister (DTP) clone residual proliferation rate (month^-1)
    K: float = 1.0            # Normalized local carrying capacity
    
    # Plastic DTP transitions
    alpha: float = 0.05       # Susceptibility of sensitive cells to DTP entry via inflammation (Phi^-1 month^-1)
    beta: float = 0.01        # Spontaneous awakening rate from DTP to sensitive (month^-1)

    # Pharmacodynamics & Transport (Bessel/diffusion)
    D_C: float = 150.0        # Drug diffusion coefficient in stroma (um^2 / month)
    k_C: float = 0.15         # Drug clearance/degradation rate (month^-1)
    E_max: float = 0.25       # Maximum drug-induced kill rate (month^-1)
    IC_50: float = 0.2        # Half-maximal drug inhibitory concentration (normalized)
    C_plasma_on: float = 1.0  # Steady plasma concentration during active therapy
    
    # Inflammatory Field (Phi)
    D_Phi: float = 80.0       # Cytokine diffusion coefficient (um^2 / month)
    k_Phi: float = 0.4        # Inflammatory clearance rate (month^-1)
    sigma_base: float = 0.1   # Baseline chronic systemic/tissue stress (NLR, CRP proxy)
    eta: float = 0.3          # Debris/damage inflammatory release per lysed cell (Phi / lysis)

    # Immunosurveillance & History-Dependent Exhaustion
    E_0: float = 0.04         # Baseline cytotoxic immune efficacy (month^-1)
    gamma: float = 0.012      # Epigenetic exhaustion rate per cumulative inflammatory unit

class InflammatoryOncoContinuum:
    def __init__(self, params: ModelParameters):
        self.p = params
        self.x = np.linspace(0, self.p.L, self.p.Nx)
        self.dx = self.x[1] - self.x[0]

    def laplacian(self, u: np.ndarray) -> np.ndarray:
        """1D Laplacian with Neumann boundary conditions at x=L and Dirichlet/Robin at x=0."""
        d2u = np.zeros_like(u)
        d2u[1:-1] = (u[2:] - 2*u[1:-1] + u[:-2]) / (self.dx**2)
        # Boundaries: zero flux at x=L (tissue interior)
        d2u[-1] = 2 * (u[-2] - u[-1]) / (self.dx**2)
        # Boundary at x=0 (capillary): zero flux for cells/cytokines or Dirichlet
        d2u[0] = 2 * (u[1] - u[0]) / (self.dx**2)
        return d2u

    def simulate(
        self,
        t_span_months: float = 60.0,
        t_therapy_stop: float = 36.0,
        dt_months: float = 0.01,
        initial_tumor_burden: float = 0.35
    ) -> Dict[str, np.ndarray]:
        """
        Simulate tumor dynamics through therapy (0 -> t_therapy_stop) 
        and post-cessation observation (t_therapy_stop -> t_span_months).
        """
        p = self.p
        nt = int(t_span_months / dt_months) + 1
        time_points = np.linspace(0, t_span_months, nt)

        # Initialize spatial fields
        N_S = np.ones(p.Nx) * initial_tumor_burden * 0.95
        N_R = np.ones(p.Nx) * initial_tumor_burden * 0.05
        Phi = np.ones(p.Nx) * p.sigma_base / p.k_Phi
        Psi = np.zeros(p.Nx)  # Cumulative inflammatory exposure integral
        
        # Pre-allocate trajectory storage
        history_time = []
        history_N_S_mean = []
        history_N_R_mean = []
        history_N_tot_mean = []
        history_Phi_mean = []
        history_Imm_mean = []
        
        # Spatial snapshots at key milestones
        spatial_snapshots = {}

        # Save interval (every 0.5 months)
        save_interval = int(0.5 / dt_months)

        for step in range(nt):
            t = time_points[step]

            # Drug boundary condition at capillary (x = 0)
            is_therapy_active = (t <= t_therapy_stop)
            C_plasma = p.C_plasma_on if is_therapy_active else 0.0

            # Quasi-steady state analytical drug penetration across stroma (decaying with distance)
            # lambda_C = sqrt(D_C / k_C)
            lam_C = np.sqrt(p.D_C / p.k_C)
            C = C_plasma * np.cosh((p.L - self.x) / lam_C) / np.cosh(p.L / lam_C)

            # Local pharmacodynamics
            H_C = p.E_max * C / (C + p.IC_50)
            lysis_rate = H_C * N_S

            # Local immune capacity (history-dependent exhaustion)
            E_imm = p.E_0 * np.exp(-p.gamma * Psi)

            # Inflammatory field PDE update
            d2Phi = self.laplacian(Phi)
            dPhi_dt = p.D_Phi * d2Phi - p.k_Phi * Phi + p.sigma_base + p.eta * lysis_rate
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
                history_N_S_mean.append(np.mean(N_S))
                history_N_R_mean.append(np.mean(N_R))
                history_N_tot_mean.append(np.mean(N_S + N_R))
                history_Phi_mean.append(np.mean(Phi))
                history_Imm_mean.append(np.mean(E_imm))

            # Store snapshots at t=0, t=12, t=36 (cessation), t=48 (relapse check)
            for target_t in [0.0, 12.0, 36.0, 48.0]:
                if abs(t - target_t) < (dt_months / 2):
                    spatial_snapshots[f"month_{int(target_t)}"] = {
                        "x": self.x.copy(),
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
            "snapshots": spatial_snapshots
        }
