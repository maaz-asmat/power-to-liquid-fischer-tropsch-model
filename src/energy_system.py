"""Hourly simulation of an off-grid Power-to-Liquid site: wind + PV -> electrolyser -> H2 storage -> RWGS/FT.

The FT plant cannot follow the weather like the electrolyser can. It needs a minimum load, limited ramp
rates and a long restart after a shutdown. The hydrogen storage decouples the two. The FT plant is
described by the steady-state process model (ft_process.py), evaluated at its design point and assumed
to scale linearly with load (constant efficiency between minimum and full load).
"""
from dataclasses import dataclass, replace
import pathlib
import numpy as np
import pandas as pd
from ft_process import FTParams, flowsheet, NM3_H2

ROOT = pathlib.Path(__file__).resolve().parents[1]


@dataclass
class Site:
    wind_MW: float = 300.0
    pv_MW: float = 200.0
    wind_flh: float = 2500.0          # full-load hours of a modern onshore turbine in coastal Lower Saxony (assumed)
    pv_flh: float = 950.0             # North German PV (assumed)


@dataclass
class Plant:
    electrolyser_MW: float = 250.0    # AC
    el_module_MW: float = 20.0        # the electrolyser is built from stacks/modules of this size
    el_min_load: float = 0.10         # minimum load of one PEM module (~10 %)
    ft_h2_t_h: float = 3.01           # FT/RWGS design H2 intake (Koenig et al. 2015 reference size)
    storage_t: float = 150.0          # usable H2 storage (compressed gas)
    ft_min_load: float = 0.40         # minimum stable load of the FT/RWGS section (assumed)
    ft_ramp_per_h: float = 0.20       # max load change per hour (assumed)
    restart_h: int = 24               # hours offline after a shutdown before it can restart (assumed)
    soc_restart: float = 0.30         # storage level needed for a restart


def load_profiles(csv=ROOT / "data/smard_2024_DE_wind_pv.csv"):
    """Normalised hourly shapes (sum = 1) of German onshore wind and PV generation from SMARD."""
    d = pd.read_csv(csv, parse_dates=["time_utc"])
    return d["time_utc"], d["wind_onshore_MWh"].to_numpy() / d["wind_onshore_MWh"].sum(), \
        d["pv_MWh"].to_numpy() / d["pv_MWh"].sum()


def simulate(site: Site, plant: Plant, ftp: FTParams = FTParams(), profiles=None, keep_series=False):
    t, w, s = profiles if profiles is not None else load_profiles()
    n = len(w)
    p_ren = site.wind_MW * site.wind_flh * w * n / 8760 + site.pv_MW * site.pv_flh * s * n / 8760   # MW
    design = flowsheet(plant.ft_h2_t_h * 1000, ftp)
    kg_h2_per_MWh = 1000.0 / (ftp.el_kwh_per_nm3 / NM3_H2) * ftp.rectifier_eff   # kg H2 per MWh AC
    h2_nom = plant.ft_h2_t_h                       # t/h
    fuel_per_t_h2 = design["fuel_t_h"] / design["h2_fresh_t_h"]
    aux_per_t_fuel = design["aux_MW"] / design["fuel_t_h"]   # MWh per t fuel

    soc, load, off_timer = 0.5 * plant.storage_t, 0.0, 0
    out = np.zeros((n, 6))
    starts = 0
    for i in range(n):
        # FT load decision based on the storage level (hysteresis + ramp limits)
        f = soc / plant.storage_t if plant.storage_t > 0 else 0.0
        if load == 0.0:
            if off_timer > 0:
                off_timer -= 1
            target = 0.0
            if off_timer == 0 and f >= plant.soc_restart:
                target = plant.ft_min_load; starts += 1
        else:
            target = np.clip(plant.ft_min_load + (1 - plant.ft_min_load) * (f - 0.2) / 0.5, plant.ft_min_load, 1.0)
            target = np.clip(target, load - plant.ft_ramp_per_h, load + plant.ft_ramp_per_h)
        # power available for the electrolyser after FT auxiliaries
        aux = target * h2_nom * fuel_per_t_h2 * aux_per_t_fuel
        p_el = min(max(p_ren[i] - aux, 0.0), plant.electrolyser_MW)
        if p_el < plant.el_min_load * plant.el_module_MW:      # below one module at minimum load
            p_el = 0.0
        h2_prod = p_el * kg_h2_per_MWh / 1000.0            # t
        h2_need = target * h2_nom
        avail = soc + h2_prod
        if avail < h2_need:                                # not enough hydrogen: emergency shutdown
            target, h2_need, aux = 0.0, 0.0, 0.0
        if target == 0.0 and load > 0.0:
            off_timer = plant.restart_h
        soc = avail - h2_need
        if soc > plant.storage_t:                          # storage full -> curtail the electrolyser
            excess = soc - plant.storage_t
            h2_prod -= excess
            p_el -= excess * 1000.0 / kg_h2_per_MWh
            soc = plant.storage_t
        load = target
        out[i] = (p_ren[i], p_el, h2_prod, load, soc, max(p_ren[i] - p_el - aux, 0.0))

    df = pd.DataFrame(out, columns=["P_renewable_MW", "P_electrolyser_MW", "H2_produced_t", "FT_load",
                                    "H2_storage_t", "P_unused_MW"])
    fuel = df.FT_load.sum() * h2_nom * fuel_per_t_h2
    kpi = {
        "renewable_GWh": df.P_renewable_MW.sum() / 1000,
        "electrolyser_flh": df.P_electrolyser_MW.sum() / plant.electrolyser_MW,
        "unused_energy_%": 100 * df.P_unused_MW.sum() / df.P_renewable_MW.sum(),
        "H2_t": df.H2_produced_t.sum(),
        "FT_flh": df.FT_load.sum(),
        "FT_hours_online": int((df.FT_load > 0).sum()),
        "FT_starts": starts,
        "fuel_t": fuel,
        "fuel_GWh": fuel * 12.2 / 1000,
        "system_efficiency_%": 100 * fuel * 12.2 / 1000 / (df.P_renewable_MW.sum() / 1000),
        "fuel_per_MW_renewable_t": fuel / (site.wind_MW + site.pv_MW),
    }
    return (kpi, df.assign(time_utc=t.values)) if keep_series else kpi
