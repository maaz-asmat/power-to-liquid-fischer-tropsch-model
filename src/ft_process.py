"""Steady-state model of a CO2-based Power-to-Liquid process: RWGS + Fischer-Tropsch with tail-gas recycle.

    fresh H2 + CO2 --> [RWGS 900 C] --> water knock-out --> [FT 225 C, Co catalyst] --> C5+ liquid fuel
                         ^                                        |
                         +---------- tail gas recycle <-----------+--> purge (fuel gas for the RWGS burner)

Units: kmol/h and kg/h. All parameters with literature sources are in `FTParams`.
"""
from dataclasses import dataclass, field
import math
import numpy as np

M = {"H2": 2.016, "CO": 28.010, "CO2": 44.009, "H2O": 18.015, "CH2": 14.027}
LHV_H2 = 33.33                 # kWh/kg (119.96 MJ/kg)
LHV_FUEL = 12.2                # kWh/kg (~44 MJ/kg) for FT C5+ hydrocarbons (diesel/kerosene range)
NM3_H2 = 0.08988               # kg per Nm3 of H2


@dataclass
class FTParams:
    # Fischer-Tropsch synthesis (low-temperature, cobalt catalyst)
    alpha: float = 0.90         # ASF chain-growth probability; cobalt LTFT range 0.70-0.95 (Koenig et al. 2015)
    x_co_pass: float = 0.50     # CO conversion per pass; assumed, typical slurry/fixed-bed value 0.4-0.6
    h2_co_ft: float = 2.05      # H2/CO ratio at the FT inlet, set via the fresh CO2 feed (as in Koenig et al.)
    T_ft_C: float = 225.0       # Koenig et al. 2015
    # Reverse water-gas shift
    T_rwgs_C: float = 900.0     # Koenig et al. 2015
    # recycle
    purge: float = 0.05         # fraction of the tail gas sent to the burner (prevents build-up of light gases)
    # electrolysis (PEM)
    el_kwh_per_nm3: float = 4.3     # system specific energy; Koenig et al. 2015 (69.7 % LHV)
    rectifier_eff: float = 0.96     # AC -> DC; Koenig et al. 2015
    aux_kwh_per_kg_fuel: float = 0.9  # compression, cooling, pumps (assumed; Koenig: ~2 % of input)


def k_rwgs(T_K):
    """Equilibrium constant of CO2 + H2 <-> CO + H2O (inverse of Moe's 1962 WGS correlation)."""
    return 1.0 / math.exp(4577.8 / T_K - 4.33)


def asf_mass_fractions(alpha, n_max=80):
    """Anderson-Schulz-Flory mass (= carbon) fraction of chain length n."""
    n = np.arange(1, n_max + 1)
    w = n * (1 - alpha) ** 2 * alpha ** (n - 1)
    return n, w / w.sum()


def product_cuts(alpha):
    """Mass fractions of the FT product classes (straight-run, before upgrading)."""
    n, w = asf_mass_fractions(alpha)
    cut = lambda lo, hi: float(w[(n >= lo) & (n <= hi)].sum())
    return {"C1-C4 gas": cut(1, 4), "C5-C9 naphtha": cut(5, 9), "C10-C14 kerosene": cut(10, 14),
            "C15-C20 diesel": cut(15, 20), "C21+ wax": cut(21, 999)}


def rwgs_equilibrium(f, T_K):
    """Shift a stream (dict of kmol/h) to RWGS equilibrium; returns the new stream."""
    K = k_rwgs(T_K)
    co2, h2, co, h2o = f["CO2"], f["H2"], f["CO"], f["H2O"]
    lo, hi = -min(co, h2o), min(co2, h2)
    g = lambda x: (co + x) * (h2o + x) - K * (co2 - x) * (h2 - x)
    for _ in range(200):                                    # bisection on the reaction extent
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if g(mid) < 0 else (lo, mid)
    x = 0.5 * (lo + hi)
    out = dict(f)
    out.update(CO2=co2 - x, H2=h2 - x, CO=co + x, H2O=h2o + x)
    return out


def ft_reactor(f, p: FTParams):
    """CO + H2 -> paraffins with ASF distribution. H2 use per C is 3 - alpha (from the ASF average)."""
    usage = 3.0 - p.alpha                                   # mol H2 per mol CO converted (paraffins + water)
    x = min(p.x_co_pass * f["CO"], f["H2"] / usage)         # CO converted, kmol/h
    n, w = asf_mass_fractions(p.alpha)
    gas_c = float(w[n <= 4].sum())                          # carbon fraction ending up in C1-C4
    out = dict(f)
    out["CO"] -= x
    out["H2"] -= usage * x
    out["H2O"] += x
    out["C1-C4"] = f.get("C1-C4", 0.0) + gas_c * x          # kmol C/h (light gases, recycled with the tail gas)
    c5 = (1 - gas_c) * x                                    # kmol C/h in liquid product
    # liquid product mass: C5+ paraffins CnH2n+2 per carbon ~ CH2 plus end-group H
    nn, ww = n[n >= 5], w[n >= 5] / w[n >= 5].sum()
    kg_per_kmolC = float(np.sum(ww * (14.027 * nn + 2.016) / nn))
    return out, c5 * kg_per_kmolC, x


def flowsheet(h2_fresh_kg_h, p: FTParams = FTParams(), tol=1e-9):
    """Solve the recycle loop for a given fresh H2 feed. Fresh CO2 is adjusted so that the
    H2/CO ratio at the FT inlet equals p.h2_co_ft. Returns a dict of results."""
    h2_in = h2_fresh_kg_h / M["H2"]

    def solve(co2_in):
        rec = {"H2": 0.0, "CO": 0.0, "CO2": 0.0, "H2O": 0.0, "C1-C4": 0.0}
        for _ in range(2000):
            feed = {k: rec[k] for k in rec}
            feed["H2"] += h2_in; feed["CO2"] += co2_in
            eq = rwgs_equilibrium(feed, p.T_rwgs_C + 273.15)
            ft_in = dict(eq); ft_in["H2O"] = 0.0             # water knock-out after RWGS
            ft_out, fuel_kg, x = ft_reactor(ft_in, p)
            tail = {k: ft_out[k] for k in ("H2", "CO", "CO2", "C1-C4")}; tail["H2O"] = 0.0
            new = {k: (1 - p.purge) * tail[k] for k in tail}
            if max(abs(new[k] - rec[k]) for k in rec) < tol:
                break
            rec = new
        return ft_in["H2"] / max(ft_in["CO"], 1e-12), dict(ft_in=ft_in, tail=tail, fuel_kg=fuel_kg, x=x,
                                                           rwgs_in=feed, rwgs_out=eq)

    lo, hi = 0.05 * h2_in, 0.6 * h2_in                      # bracket on fresh CO2 (kmol/h)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        ratio, s = solve(mid)
        lo, hi = (lo, mid) if ratio < p.h2_co_ft else (mid, hi)
    co2_in = 0.5 * (lo + hi)
    ratio, s = solve(co2_in)

    purge = {k: p.purge * s["tail"][k] for k in s["tail"]}
    fuel = s["fuel_kg"]
    c_in = co2_in
    c_fuel = c_in - purge["CO"] - purge["CO2"] - purge["C1-C4"]
    el_dc = h2_fresh_kg_h / NM3_H2 * p.el_kwh_per_nm3 / 1000.0          # MW
    el_ac = el_dc / p.rectifier_eff
    aux = fuel * p.aux_kwh_per_kg_fuel / 1000.0
    # heat for the RWGS: reaction enthalpy (+41 kJ/mol) plus the part of the feed pre-heating to 900 C that a
    # feed/effluent heat exchanger cannot recover (assumed 20 %, mean cp ~ 31 J/mol K)
    rwgs_extent = s["rwgs_out"]["CO"] - s["rwgs_in"]["CO"]
    n_feed = sum(s["rwgs_in"].values())
    q_rwgs = (rwgs_extent * 41.2e3 + 0.2 * n_feed * 31.0 * (p.T_rwgs_C - 25)) * 1000 / 3.6e9   # MW
    lhv_purge = (purge["H2"] * M["H2"] * LHV_H2 + purge["CO"] * 28.01 * 2.80
                 + purge["C1-C4"] * 14.027 * 13.3) / 1000.0              # MW (CO 10.1 MJ/kg, light HC ~48 MJ/kg)
    return {
        "h2_fresh_t_h": h2_fresh_kg_h / 1000, "co2_fresh_t_h": co2_in * M["CO2"] / 1000,
        "fuel_t_h": fuel / 1000, "fuel_MW_LHV": fuel * LHV_FUEL / 1000,
        "h2_MW_LHV": h2_fresh_kg_h * LHV_H2 / 1000,
        "electrolyser_MW_AC": el_ac, "aux_MW": aux,
        "carbon_efficiency": c_fuel / c_in,
        "h2_to_fuel_efficiency": fuel * LHV_FUEL / (h2_fresh_kg_h * LHV_H2),
        "ptl_efficiency": fuel * LHV_FUEL / 1000 / (el_ac + aux),
        "recycle_ratio": sum(s["tail"].values()) * (1 - p.purge) / (h2_in + co2_in),
        "co_conversion_overall": s["x"] / (s["x"] + purge["CO"]),
        "rwgs_heat_MW": q_rwgs, "purge_LHV_MW": lhv_purge,
        "h2_co_ft_inlet": ratio, "co2_conversion_rwgs_pass": rwgs_extent / s["rwgs_in"]["CO2"],
    }


def calibrate_purge(h2_fresh_kg_h, target_carbon_eff, p: FTParams = FTParams()):
    """Purge fraction that reproduces a given carbon efficiency (used to match a literature reference plant)."""
    from dataclasses import replace
    lo, hi = 0.001, 0.6
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        ce = flowsheet(h2_fresh_kg_h, replace(p, purge=mid))["carbon_efficiency"]
        lo, hi = (mid, hi) if ce > target_carbon_eff else (lo, mid)
    return 0.5 * (lo + hi)
