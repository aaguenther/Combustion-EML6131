from pathlib import Path
import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "thermodynamics"

_ALIASES = {
    "C(S)": "C_graphite",
    "C(s)": "C_graphite",
    "C": "C_graphite",
    "graphite": "C_graphite",
    "C8H18": "C8H18",
}

_PROPERTY_COLUMNS = {
    "cp": "cp_J_per_mol_K",
    "entropy": "entropy_J_per_mol_K",
    "s": "entropy_J_per_mol_K",
    "minus_g_minus_href_over_t": "minus_G_minus_Href_over_T_J_per_mol_K",
    "delta_h_sensible": "delta_H_sensible_kJ_per_mol",
    "delta_hf": "delta_Hf_kJ_per_mol",
    "delta_gf": "delta_Gf_kJ_per_mol",
    "log_kf": "log10_Kf",
}

def _key(species: str) -> str:
    return _ALIASES.get(species, species)

def load_species(species: str) -> pd.DataFrame:
    key = _key(species)
    path = DATA_DIR / f"{key}.csv"
    if not path.exists():
        raise FileNotFoundError(f"No thermodynamic table found for {species!r}: {path}")
    return pd.read_csv(path)

def get_property(species: str, temperature_K: float, property_name: str) -> float:
    """Linearly interpolate a tabulated property at temperature_K."""
    key = _key(species)
    df = load_species(key)
    if property_name not in _PROPERTY_COLUMNS:
        raise ValueError(f"Unknown property {property_name!r}. Choose from {sorted(_PROPERTY_COLUMNS)}")
    col = _PROPERTY_COLUMNS[property_name]
    x = df["temperature_K"].to_numpy(dtype=float)
    y = df[col].to_numpy(dtype=float)
    if np.isnan(y).all():
        raise ValueError(f"Property {property_name!r} is not available for {species} in the supplied tables.")
    if temperature_K < x.min() or temperature_K > x.max():
        raise ValueError(f"Temperature {temperature_K} K is outside the tabulated range {x.min()}-{x.max()} K for {species}.")
    mask = np.isfinite(y)
    x_valid = x[mask]
    y_valid = y[mask]
    if temperature_K < x_valid.min() or temperature_K > x_valid.max():
        raise ValueError(f"Temperature {temperature_K} K is outside the finite data range for {species} and {property_name}.")
    return float(np.interp(temperature_K, x_valid, y_valid))

def h_standard(species: str, temperature_K: float) -> float:
    """Return H°(T) using Delta Hf°(298.15 K) + [H°(T)-H°(298.15 K)].

    Result is in kJ/mol. This is a derived quantity, not a separate textbook column.
    """
    df = load_species(species)
    row = df.loc[np.isclose(df["temperature_K"], 298.15)]
    if row.empty:
        raise ValueError(f"298.15 K reference row is unavailable for {species}.")
    hf_298 = float(row.iloc[0]["delta_Hf_kJ_per_mol"])
    dh = get_property(species, temperature_K, "delta_h_sensible")
    return hf_298 + dh
