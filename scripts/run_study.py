"""Process validation, parameter studies and energy-system scenarios for the Power-to-Liquid model."""
import sys, pathlib, csv
from dataclasses import replace
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from ft_process import FTParams, flowsheet, calibrate_purge, product_cuts
from energy_system import Site, Plant, simulate, load_profiles

OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
C = ["#1f5f8b", "#0e7c7b", "#d97706", "#7c3aed", "#b91c1c"]


def write(name, header, rows):
    with open(OUT / name, "w", newline="") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)


# ===================================================================== 1. process model vs Koenig et al. (2015)
REF = {"h2_fresh_t_h": 3.01, "electrolyser_MW_AC": 150.1, "fuel_t_h": 5.5, "fuel_MW_LHV": 67.1,
       "carbon_efficiency": 0.737, "ptl_efficiency": 0.433, "rwgs_heat_MW": 11.5}
purge = calibrate_purge(3010, REF["carbon_efficiency"])
ftp = FTParams(purge=purge)
r = flowsheet(3010, ftp)
print(f"purge fraction calibrated to the reference carbon efficiency: {purge:.3f}")
rows = []
for k, ref in REF.items():
    rows.append([k, round(r[k], 3), ref, round(100 * (r[k] / ref - 1), 1)])
    print(rows[-1])
write("01_process_validation.csv", ["quantity", "model", "Koenig_et_al_2015", "difference_%"], rows)
write("01_design_point.csv", ["quantity", "value"], [[k, round(v, 4)] for k, v in r.items()] + [["purge_fraction", round(purge, 4)]])

# ===================================================================== 2. FT parameter studies
alphas = np.round(np.arange(0.80, 0.951, 0.025), 3)
rows, cuts_all = [], []
for a in alphas:
    rr = flowsheet(3010, replace(ftp, alpha=a)); cu = product_cuts(a); cuts_all.append(cu)
    rows.append([a, round(rr["fuel_t_h"], 3), round(rr["carbon_efficiency"], 3), round(rr["ptl_efficiency"], 3)]
                + [round(100 * v, 1) for v in cu.values()])
write("02_alpha_study.csv", ["alpha", "fuel_t_h", "carbon_efficiency", "ptl_efficiency"] + [f"{k}_mass_%" for k in cuts_all[0]], rows)

xs = np.round(np.arange(0.3, 0.81, 0.1), 2)
purges = [0.05, 0.10, purge, 0.30]
rows = []
for x in xs:
    for pg in purges:
        rr = flowsheet(3010, replace(ftp, x_co_pass=x, purge=pg))
        rows.append([x, round(pg, 3), round(rr["carbon_efficiency"], 3), round(rr["ptl_efficiency"], 3),
                     round(rr["recycle_ratio"], 2), round(rr["purge_LHV_MW"], 1), round(rr["rwgs_heat_MW"], 1)])
write("02_conversion_purge_study.csv", ["x_co_per_pass", "purge", "carbon_efficiency", "ptl_efficiency",
                                        "recycle_ratio", "purge_LHV_MW", "rwgs_heat_MW"], rows)

fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
bottom = np.zeros(len(alphas))
for i, k in enumerate(cuts_all[0]):
    v = np.array([100 * c[k] for c in cuts_all])
    ax[0].bar(alphas, v, bottom=bottom, width=0.02, label=k, color=(["0.6"] + C)[i]); bottom += v
ax[0].set(xlabel="Chain-growth probability α", ylabel="Mass fraction [%]", title="(a) FT product distribution (ASF)")
ax[0].legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3)
for j, pg in enumerate(purges):
    d = [rr for rr in rows if rr[1] == round(pg, 3)]
    ax[1].plot([q[0] for q in d], [100 * q[2] for q in d], "o-", color=C[j], label=f"purge {pg:.0%}")
    ax[2].plot([q[0] for q in d], [100 * q[3] for q in d], "o-", color=C[j], label=f"purge {pg:.0%}")
ax[1].axhline(73.7, color="k", ls=":", lw=1); ax[1].text(0.31, 74.5, "Koenig et al. (2015)", fontsize=7)
ax[2].axhline(43.3, color="k", ls=":", lw=1)
ax[1].set(xlabel="CO conversion per FT pass", ylabel="Carbon efficiency [%]", title="(b) CO2 → fuel carbon efficiency")
ax[2].set(xlabel="CO conversion per FT pass", ylabel="Power-to-Liquid efficiency [%]", title="(c) Electricity → fuel efficiency")
for a_ in ax[1:]:
    a_.grid(alpha=.3); a_.legend(fontsize=7)
fig.tight_layout(); fig.savefig(OUT / "02_ft_parameter_study.png", dpi=150)

# ===================================================================== 3. energy-system scenarios
prof = load_profiles()
base_site, base_plant = Site(), Plant()
kpi, ts = simulate(base_site, base_plant, ftp, prof, keep_series=True)
print("baseline", {k: round(float(v), 1) for k, v in kpi.items()})
write("03_baseline_kpis.csv", ["kpi", "value"], [[k, round(float(v), 2)] for k, v in kpi.items()])

# (a) storage size x FT minimum load
storages = [10, 25, 50, 100, 150, 250, 400, 600]
rows = []
for ml in (0.2, 0.4, 0.6):
    for st in storages:
        k = simulate(base_site, replace(base_plant, storage_t=st, ft_min_load=ml), ftp, prof)
        rows.append([ml, st, round(st / base_plant.ft_h2_t_h, 1), round(k["FT_flh"]), k["FT_starts"],
                     round(k["fuel_t"]), round(k["unused_energy_%"], 1)])
write("03_storage_minload_study.csv", ["ft_min_load", "storage_t", "storage_hours_of_ft_demand", "FT_full_load_hours",
                                       "FT_starts", "fuel_t_per_year", "unused_renewable_%"], rows)
st_rows = rows

# (b) wind/PV mix at the same annual renewable energy
E = base_site.wind_MW * base_site.wind_flh + base_site.pv_MW * base_site.pv_flh       # MWh
rows = []
for share in (0.0, 0.25, 0.5, 0.75, 0.85, 1.0):
    site = replace(base_site, wind_MW=share * E / base_site.wind_flh, pv_MW=(1 - share) * E / base_site.pv_flh)
    k = simulate(site, base_plant, ftp, prof)
    rows.append([share, round(site.wind_MW), round(site.pv_MW), round(k["FT_flh"]), k["FT_starts"],
                 round(k["fuel_t"]), round(k["unused_energy_%"], 1), round(k["electrolyser_flh"])])
write("03_wind_pv_mix.csv", ["wind_energy_share", "wind_MW", "pv_MW", "FT_full_load_hours", "FT_starts",
                             "fuel_t_per_year", "unused_renewable_%", "electrolyser_full_load_hours"], rows)
mix_rows = rows

# (c) electrolyser size
rows = []
for el in (150, 200, 250, 300, 350, 400):
    k = simulate(base_site, replace(base_plant, electrolyser_MW=el), ftp, prof)
    rows.append([el, round(el / (base_site.wind_MW + base_site.pv_MW), 2), round(k["electrolyser_flh"]),
                 round(k["FT_flh"]), round(k["fuel_t"]), round(k["unused_energy_%"], 1)])
write("03_electrolyser_size.csv", ["electrolyser_MW", "electrolyser_to_renewable_ratio", "electrolyser_full_load_hours",
                                   "FT_full_load_hours", "fuel_t_per_year", "unused_renewable_%"], rows)
el_rows = rows

# ---- figures
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
for j, ml in enumerate((0.2, 0.4, 0.6)):
    d = [q for q in st_rows if q[0] == ml]
    ax[0].plot([q[1] for q in d], [q[3] for q in d], "o-", color=C[j], label=f"FT min. load {ml:.0%}")
    ax[1].plot([q[1] for q in d], [q[4] for q in d], "o-", color=C[j], label=f"FT min. load {ml:.0%}")
ax[0].set(xscale="log", xlabel="Usable H2 storage [t]", ylabel="FT full-load hours [h/a]", title="(a) FT utilisation")
ax[1].set(xscale="log", yscale="log", xlabel="Usable H2 storage [t]", ylabel="FT starts per year", title="(b) FT shutdowns / restarts")
ax[2].plot([q[0] * 100 for q in mix_rows], [q[5] / 1000 for q in mix_rows], "o-", color=C[0], label="fuel [kt/a]")
ax2 = ax[2].twinx()
ax2.plot([q[0] * 100 for q in mix_rows], [q[4] for q in mix_rows], "s--", color=C[2], label="FT starts")
ax[2].set(xlabel="Wind share of renewable energy [%]", ylabel="Fuel [kt per year]", title="(c) Wind/PV mix (same annual energy)")
ax2.set_ylabel("FT starts per year", color=C[2])
for a_ in ax[:2]:
    a_.grid(alpha=.3, which="both"); a_.legend(fontsize=8)
ax[2].grid(alpha=.3)
fig.tight_layout(); fig.savefig(OUT / "03_system_scenarios.png", dpi=150)

# time series: one winter and one summer fortnight
fig, ax = plt.subplots(2, 1, figsize=(13, 6.5), sharey=False)
for a_, (start, title) in zip(ax, [("2024-01-08", "Winter"), ("2024-06-10", "Summer")]):
    s = ts[(ts.time_utc >= start) & (ts.time_utc < np.datetime64(start) + np.timedelta64(14, "D"))]
    a_.fill_between(s.time_utc, s.P_renewable_MW, color="#9ecae1", label="wind + PV [MW]")
    a_.plot(s.time_utc, s.P_electrolyser_MW, color=C[0], lw=1, label="electrolyser [MW]")
    a_.plot(s.time_utc, s.FT_load * 100, color=C[2], lw=2, label="FT load [%]")
    a_.plot(s.time_utc, s.H2_storage_t, color=C[3], lw=1.5, label="H2 storage [t]")
    a_.set(title=f"{title}: two weeks of operation (baseline)", ylabel="MW / % / t"); a_.grid(alpha=.3)
ax[0].legend(fontsize=8, ncol=4, loc="upper right")
fig.tight_layout(); fig.savefig(OUT / "03_operation_timeseries.png", dpi=150)

# monthly production
ts["month"] = ts.time_utc.dt.month
fuel_per_load_h = kpi["fuel_t"] / ts.FT_load.sum()
mon = ts.groupby("month").agg(ren=("P_renewable_MW", "sum"), ft=("FT_load", "sum"))
write("03_monthly.csv", ["month", "renewable_GWh", "fuel_t"],
      [[m, round(r_.ren / 1000, 1), round(r_.ft * fuel_per_load_h)] for m, r_ in mon.iterrows()])
print("done")
