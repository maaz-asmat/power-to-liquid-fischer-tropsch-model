"""Download hourly German wind and solar generation for one year from SMARD (Bundesnetzagentur).

Source: www.smard.de, licence CC BY 4.0. Output: data/smard_<year>_DE_wind_pv.csv (MWh per hour).
Usage:  python scripts/fetch_smard.py 2024
"""
import sys, json, csv, pathlib, urllib.request, datetime as dt

ROOT = pathlib.Path(__file__).resolve().parents[1]
FILTERS = {"wind_onshore_MWh": 4067, "wind_offshore_MWh": 1225, "pv_MWh": 4068}
BASE = "https://www.smard.de/app/chart_data/{f}/DE/{name}"


def get(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.load(r)


def series(fid, year):
    start = dt.datetime(year, 1, 1, tzinfo=dt.timezone.utc).timestamp() * 1000
    end = dt.datetime(year + 1, 1, 1, tzinfo=dt.timezone.utc).timestamp() * 1000
    weeks = get(BASE.format(f=fid, name="index_hour.json"))["timestamps"]
    out = {}
    for w in weeks:
        if w < start - 8 * 86400e3 or w > end:
            continue
        for t, v in get(BASE.format(f=fid, name=f"{fid}_DE_hour_{w}.json"))["series"]:
            if start <= t < end and v is not None:
                out[t] = v
    return out


if __name__ == "__main__":
    year = int(sys.argv[1]) if len(sys.argv) > 1 else 2024
    data = {k: series(f, year) for k, f in FILTERS.items()}
    times = sorted(set.intersection(*(set(d) for d in data.values())))
    path = ROOT / "data" / f"smard_{year}_DE_wind_pv.csv"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time_utc"] + list(FILTERS))
        for t in times:
            w.writerow([dt.datetime.fromtimestamp(t / 1000, dt.timezone.utc).strftime("%Y-%m-%d %H:%M")]
                       + [data[k][t] for k in FILTERS])
    print(f"{len(times)} hours written to {path}")
