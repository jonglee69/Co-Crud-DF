#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
실험 원자료(시편 치수·질량·ICP 농도·용액 부피·희석배수) → df_calc.py 입력 CSV 변환기

사용법:
    python prepare_input.py examples/raw_lab_data_example.csv -o my_input.csv
    python prepare_input.py raw.csv -o input.csv --area-mode all     # 측면까지 면적에 포함

변환식:
    area_cm2          = 2·L·W / 100                         (--area-mode faces, 기본: 양면)
                      = 2·(L·W + L·T + W·T) / 100           (--area-mode all: 측면 포함)
    dissolved_mass_mg = mass_before_mg − mass_after_mg
    Co_ug             = (Co_mgL − Co_blank_mgL) × dilution × volume_mL      (mg/L × mL = µg)
    u_Co_ug           = Co_ug × Co_RSD_pct / 100            (RSD 입력 시)
    crud_mass_mg      = crud_after_mg − crud_before_mg      (A 시편)
    spike_Co_ug       = spike_conc_mgL × spike_vol_uL / 1000 (D 시편; mg/L × µL = ng)
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RAW_REQUIRED = ["sample_id", "group", "Co_mgL", "Ni_mgL", "volume_mL"]
RAW_OPTIONAL = ["length_mm", "width_mm", "thickness_mm", "area_cm2",
                "mass_before_mg", "mass_after_mg", "dilution",
                "Co_blank_mgL", "Ni_blank_mgL", "Co_RSD_pct", "Ni_RSD_pct",
                "crud_before_mg", "crud_after_mg", "spike_conc_mgL", "spike_vol_uL", "note"]
OUT_COLS = ["sample_id", "group", "area_cm2", "dissolved_mass_mg", "Co_ug", "Ni_ug",
            "u_Co_ug", "u_Ni_ug", "crud_mass_mg", "spike_Co_ug", "note"]


def main(argv=None):
    p = argparse.ArgumentParser(description="실험 원자료 → df_calc.py 입력 CSV")
    p.add_argument("raw", type=Path)
    p.add_argument("-o", "--out", type=Path, default=Path("input.csv"))
    p.add_argument("--area-mode", choices=["faces", "all"], default="faces",
                   help="faces: 앞뒤 양면만(측면 마스킹), all: 측면 포함")
    a = p.parse_args(argv)

    df = pd.read_csv(a.raw, comment="#", skipinitialspace=True, encoding="utf-8-sig")
    df.columns = [c.strip() for c in df.columns]
    miss = [c for c in RAW_REQUIRED if c not in df.columns]
    if miss:
        sys.exit(f"[오류] 필수 열 없음: {miss}")
    for c in RAW_OPTIONAL:
        if c not in df.columns:
            df[c] = np.nan
    num = [c for c in RAW_REQUIRED + RAW_OPTIONAL if c not in ("sample_id", "group", "note")]
    for c in num:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["group"] = df["group"].astype(str).str.strip().str.upper()
    df["dilution"] = df["dilution"].fillna(1.0)
    df["Co_blank_mgL"] = df["Co_blank_mgL"].fillna(0.0)
    df["Ni_blank_mgL"] = df["Ni_blank_mgL"].fillna(0.0)

    # 면적: area_cm2 를 직접 넣었으면 그 값을 우선
    L, W, T = df["length_mm"], df["width_mm"], df["thickness_mm"].fillna(0.0)
    if a.area_mode == "faces":
        calc_area = 2 * L * W / 100.0
    else:
        calc_area = 2 * (L * W + L * T + W * T) / 100.0
    area = df["area_cm2"].fillna(calc_area)
    if area.isna().any():
        sys.exit(f"[오류] 면적을 구할 수 없는 행(치수 또는 area_cm2 필요): {list(df.sample_id[area.isna()])}")

    out = pd.DataFrame({
        "sample_id": df["sample_id"],
        "group": df["group"],
        "area_cm2": area,
        "dissolved_mass_mg": df["mass_before_mg"] - df["mass_after_mg"],
        "Co_ug": (df["Co_mgL"] - df["Co_blank_mgL"]) * df["dilution"] * df["volume_mL"],
        "Ni_ug": (df["Ni_mgL"] - df["Ni_blank_mgL"]) * df["dilution"] * df["volume_mL"],
        "crud_mass_mg": df["crud_after_mg"] - df["crud_before_mg"],
        "spike_Co_ug": df["spike_conc_mgL"] * df["spike_vol_uL"] / 1000.0,
        "note": df["note"],
    })
    out["u_Co_ug"] = out["Co_ug"] * df["Co_RSD_pct"] / 100.0
    out["u_Ni_ug"] = out["Ni_ug"] * df["Ni_RSD_pct"] / 100.0

    # 점검
    warn = []
    neg = out[(out.Co_ug <= 0) | (out.Ni_ug <= 0)]
    if len(neg):
        warn.append(f"Co 또는 Ni 질량이 0 이하(시약 블랭크 과대?): {list(neg.sample_id)}")
    dm = out["dissolved_mass_mg"]
    if (dm < 0).any():
        warn.append(f"용해 후 질량이 더 큼(질량 순서 확인): {list(out.sample_id[dm < 0])}")
    nod = out[out.group.isin(["B", "C", "D", "E"]) & dm.isna()]
    if len(nod):
        warn.append(f"B·C·D·E 인데 용해 전후 질량이 없음(층 두께 계산 불가): {list(nod.sample_id)}")
    if out[out.group == "D"]["spike_Co_ug"].isna().any():
        warn.append("D(스파이크) 시편에 spike_conc_mgL·spike_vol_uL 이 비어 있음")
    bad = set(out.group) - set("ALBCDEF")
    if bad:
        warn.append(f"알 수 없는 group: {sorted(bad)} (A, L, B, C, D, E, F 중 하나)")

    out = out[OUT_COLS]
    out.to_csv(a.out, index=False, encoding="utf-8-sig", float_format="%.6g")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(out.to_string(index=False))
    for w in warn:
        print("[주의]", w)
    print(f"\n[저장] {a.out}  →  python df_calc.py {a.out}")


if __name__ == "__main__":
    main()
