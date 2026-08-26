"""Soil evaporation capacitance (SEC) model utilities."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

# Stage-2 evaporation constants (Or & Lehman 2019 appendix workflow)
VAP_JUMP = 10
E2 = 1

# Tropical-cyclone infiltration years in the shipped shapefiles
DEFAULT_YEARS = [2011, 2018, 2020]

SHAPEFILE_TYPES = [
    "Dunes.shp",
    "Alluvium.shp",
    "Hard Limestone.shp",
    "Karstified Limestone.shp",
]

SOIL_TYPE_PARAM_MAP = {
    "dunes": "USMo1",
    "karstified_limestone": "USMo1",
    "alluvium": "USMo7",
    "hard_limestone": "USDk1",
}


def repo_root(start: Path | None = None) -> Path:
    """Return repository root whether cwd is root or notebooks/."""
    base = Path.cwd() if start is None else start
    if not (base / "data").exists() and (base.parent / "data").exists():
        base = base.parent
    return base


def load_soil_params(base_dir: Path | None = None) -> pd.DataFrame:
    root = repo_root(base_dir)
    return pd.read_excel(root / "data" / "soil_param.xlsx", sheet_name="code")


def load_infiltration_gdfs(
    years: list[int] | None = None,
    base_dir: Path | None = None,
) -> dict[int, dict[str, gpd.GeoDataFrame]]:
    root = repo_root(base_dir)
    years = years or DEFAULT_YEARS
    shapefile_base = root / "data" / "Infiltration_output"
    gdfs: dict[int, dict[str, gpd.GeoDataFrame]] = {}
    for year in years:
        year_path = shapefile_base / str(year)
        gdfs[year] = {}
        for shapefile in SHAPEFILE_TYPES:
            key = shapefile.split(".")[0].replace(" ", "_").lower()
            gdfs[year][key] = gpd.read_file(year_path / shapefile)
    return gdfs


def calculate_k_eff(ks: float, alpha: float, h_crit: float, n: float, m: float) -> float:
    """Effective hydraulic conductivity (mm/day)."""
    term = (1 + (alpha * h_crit) ** n) ** (-m)
    inner = 1 - (1 - (1 / (1 + (alpha * h_crit) ** n))) ** m
    k_eff = 4 * ks * np.sqrt(term) * inner**2
    return k_eff * 1000


def simulate_soil_evaporation_volume(
    gdf: gpd.GeoDataFrame,
    param_row: pd.Series,
) -> pd.DataFrame:
    """
    Simulate cumulative evaporation volume (km³) for one soil parameter set.

    Parameters
    ----------
    gdf : GeoDataFrame
        Must contain ``s_mm`` (initial moisture, mm) and ``area_m2``.
    param_row : Series
        One row from the soil parameters table (``site`` column names output).
    """
    delta_theta = param_row["theta_crit"] - param_row["theta_res"] / 2
    dfs: list[pd.DataFrame] = []

    for cell_num, (s, area) in enumerate(
        tqdm(
            zip(gdf["s_mm"], gdf["area_m2"]),
            desc=f"Processing {param_row['site']}",
            total=len(gdf),
        ),
        start=1,
    ):
        theta = s / (param_row["lc"] * 1000)
        if theta >= param_row["theta_sat"]:
            theta = param_row["theta_sat"]
            s = theta * (param_row["lc"] * 1000)
        elif theta <= param_row["theta_res"]:
            continue
        else:
            s = s

        day = 0
        days2 = 0
        t2 = 0
        cumulative_vol = 0.0
        cum_s = s
        es = 0.0
        q = 0.0
        q_current = 0.0
        es_current = 0.0
        es_list: list[float] = []

        while theta > param_row["theta_res"]:
            if theta > param_row["theta_crit"]:
                day += 1
                t2 = day
                theta -= 0.01
                q_pre = q
                theta_eff = (theta - param_row["theta_res"]) / (
                    param_row["theta_sat"] - param_row["theta_res"]
                )
                k_eff = calculate_k_eff(
                    param_row["ks"],
                    param_row["alpha"],
                    param_row["h_crit"],
                    param_row["n"],
                    param_row["m"],
                )
                k_theta = (
                    param_row["ks"]
                    * 1000
                    * (theta_eff) ** 0.5
                    * (
                        1
                        - (1 - theta_eff ** (1 / (1 - (1 / param_row["n"]))))
                        ** (1 - (1 / param_row["n"]))
                    )
                    ** 2
                )
                es_current = (param_row["et"] * k_theta * (1 + (param_row["et"] / k_eff))) / (
                    param_row["et"] + (k_theta * (1 + (param_row["et"] / k_eff)))
                )
                if (cum_s - es_current) / (param_row["lc"] * 1000) > param_row["theta_crit"]:
                    q = k_theta * np.exp(
                        -(day * k_theta)
                        / ((param_row["lc"] * 1000) * (theta - param_row["theta_crit"]))
                    )
                    q_current = q - q_pre
            else:
                q_current = 0.0
                es_previous = es
                days2 += 1
                t = t2 + days2
                es = (
                    np.sqrt(VAP_JUMP)
                    / np.sqrt(VAP_JUMP + (2 * E2 * (t - t2)) / delta_theta)
                    * (2 * E2 * (t - t2) + (delta_theta * VAP_JUMP))
                    - (VAP_JUMP * delta_theta)
                )
                es_current = es - es_previous

            cum_s -= es_current + q_current
            es_vol = es_current * (10**-6) * area * (10**-6)
            cumulative_vol += es_vol
            theta = cum_s / (param_row["lc"] * 1000)
            es_list.append(cumulative_vol)

        dfs.append(pd.DataFrame({cell_num: es_list}))

    if not dfs:
        return pd.DataFrame({param_row["site"]: []})

    soil_evaporation_df = pd.concat(dfs, axis=1)
    soil_evaporation_df_filled = soil_evaporation_df.ffill(axis=0)
    sum_columns = soil_evaporation_df_filled.sum(axis=1)
    return pd.DataFrame({param_row["site"]: sum_columns})


def make_result_key(year: int, soil_type: str) -> str:
    return f"{year}|{soil_type}"


def parse_result_key(key: str) -> tuple[int, str]:
    year_str, soil_type = key.split("|", 1)
    return int(year_str), soil_type


def param_row_for_soil_type(soil_params: pd.DataFrame, soil_type: str) -> pd.Series:
    site = SOIL_TYPE_PARAM_MAP[soil_type.lower()]
    match = soil_params.loc[soil_params["site"] == site]
    if match.empty:
        raise KeyError(f"No parameter row for site {site!r} (soil type {soil_type!r})")
    return match.iloc[0]


def run_all_simulations(
    gdfs: dict[int, dict[str, gpd.GeoDataFrame]],
    soil_params: pd.DataFrame,
    years: list[int] | None = None,
) -> dict[str, pd.DataFrame]:
    years = years or sorted(gdfs.keys())
    results: dict[str, pd.DataFrame] = {}
    for year in years:
        for soil_type, gdf in gdfs[year].items():
            param_row = param_row_for_soil_type(soil_params, soil_type)
            key = make_result_key(year, soil_type)
            results[key] = simulate_soil_evaporation_volume(gdf, param_row)
    return results


def prepare_evaporation_data(results: dict[str, pd.DataFrame]) -> dict[int, pd.DataFrame]:
    """
    Aggregate soil-type series by year.

    Columns with different simulation lengths are forward-filled to the longest
    series before summing (same behavior as the original notebook).
    """
    processed: dict[int, pd.DataFrame] = {}
    years = sorted({parse_result_key(k)[0] for k in results})

    for year in years:
        year_data = {
            parse_result_key(k)[1]: df
            for k, df in results.items()
            if parse_result_key(k)[0] == year
        }
        if not year_data:
            continue

        max_length = max(df.shape[0] for df in year_data.values())
        total_evap_df = pd.DataFrame(index=range(max_length))

        for soil_type, df in year_data.items():
            if df.empty:
                continue
            df_reindexed = df.reindex(range(max_length)).ffill()
            total_evap_df[soil_type] = df_reindexed.iloc[:, 0]

        total_evap_df["Total"] = total_evap_df.sum(axis=1)
        processed[year] = total_evap_df

    return processed


def plot_soil_evaporation(
    processed_data: dict[int, pd.DataFrame],
    save_path: Path | str | None = None,
    show: bool = True,
) -> plt.Figure:
    years = sorted(processed_data.keys())
    labels = ["(a)", "(b)", "(c)"]
    fig, axes = plt.subplots(1, len(years), figsize=(16, 5), gridspec_kw={"wspace": 0.2})
    if len(years) == 1:
        axes = [axes]

    for idx, year in enumerate(years):
        ax = axes[idx]
        df = processed_data[year]
        if df is None or df.empty:
            continue

        soil_cols = [c for c in df.columns if c != "Total"]
        for col in soil_cols:
            ax.plot(df.index, df[col], label=col.replace("_", " ").title())

        ax.plot(df.index, df["Total"], label="Total", linestyle="--", linewidth=2, color="black")
        ax.set_xlabel("Days", fontsize=15)
        ax.set_xlim(0, len(df))
        max_y_value = df[soil_cols].max().max()
        ax.set_ylim(0, max_y_value + 1)

        if idx == 0:
            ax.set_ylabel("Accumulated Soil Evaporation (km$^3$)", fontsize=15)
            ax.legend(loc="best", fontsize=12)

        ax.tick_params(axis="x", rotation=45, labelsize=14)
        ax.tick_params(axis="y", labelsize=14)
        ax.text(
            0.98,
            0.9,
            labels[idx],
            transform=ax.transAxes,
            fontsize=16,
            fontweight="bold",
            va="center",
            ha="right",
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
        )

        if len(df) > 100:
            diff = df["Total"].iloc[-1] - df["Total"].iloc[99]
            print(f"Difference between last value and day 100 for {year}: {diff:.3f} km³")

    if len(years) > 1:
        fig.tight_layout()
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
        print(f"Saved figure: {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig
