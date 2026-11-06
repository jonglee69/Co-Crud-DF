#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Co 도핑 모의 crud 표면제염 — DF 정량평가 계산기
(Currie 1968 / ISO 11929 검출한계 이론 + GUM 불확도 전파)

사용법:
    python df_calc.py examples/example_data.csv
    python df_calc.py my_data.csv --out results --rel-u-co 0.02 --rel-u-ni 0.02

입력 CSV 형식은 README.md 와 templates/template.csv 참조.

시편 그룹:
    A  제염 전 기준 시편 — 산화막만 선택 용해 (A0 방법 a, 방법 b)
    L  제염액 (A0 방법 c, 확인용 물질수지)
    B  제염 시편 — 표면 얇은 층(~1 µm) 용해
    C  공정 블랭크 — crud 없는 모재, B와 동일 제염·동일 분석 (r, σ0 결정)
    D  스파이크 — B와 동일 시편에 Co 표준 첨가 (회수율)
    E  양성대조 — 불완전 제염 (검출 능력 입증)
    F  모재 원소재 전체 용해 (참고값, 계산에는 사용 안 함)
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REQUIRED = ["sample_id", "group", "area_cm2", "Co_ug", "Ni_ug"]
OPTIONAL = ["u_Co_ug", "u_Ni_ug", "dissolved_mass_mg", "crud_mass_mg", "spike_Co_ug", "note"]
LAYER_GROUPS = ["B", "C", "D", "E"]          # 표면 얇은 층 분석 그룹
MAX_LAYER_UM = 10.0                          # 이보다 두꺼우면 경고


# ----------------------------------------------------------------------------- 입력
def load_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, comment="#", skipinitialspace=True, encoding="utf-8-sig")
    df.columns = [c.strip() for c in df.columns]
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        sys.exit(f"[오류] 필수 열이 없습니다: {missing}\n  필수 열: {REQUIRED}")
    for c in OPTIONAL:
        if c not in df.columns:
            df[c] = np.nan
    df["group"] = df["group"].astype(str).str.strip().str.upper()
    for c in ["area_cm2", "Co_ug", "Ni_ug", "u_Co_ug", "u_Ni_ug",
              "dissolved_mass_mg", "crud_mass_mg", "spike_Co_ug"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    bad = df[df[["area_cm2", "Co_ug", "Ni_ug"]].isna().any(axis=1)]
    if len(bad):
        sys.exit(f"[오류] area_cm2 / Co_ug / Ni_ug 값이 비어 있는 행: {list(bad.sample_id)}")
    if (df["area_cm2"] <= 0).any():
        sys.exit("[오류] area_cm2 는 0보다 커야 합니다.")
    return df


def tq(p: float, n: int, use_z: bool) -> float:
    """단측 분위수: n개 반복이면 자유도 n-1 의 t, 아니면 정규분포 z."""
    if use_z or n < 2:
        return float(stats.norm.ppf(p))
    return float(stats.t.ppf(p, n - 1))


def summarize(values, uncs, conf: float, use_z: bool):
    """평균, 표준편차, n, 평균의 표준불확도, 단측 하한/상한."""
    v = np.asarray(values, float)
    n = len(v)
    if n == 0:
        return None
    mean = float(v.mean())
    if n >= 2:
        sd = float(v.std(ddof=1))
        sem = sd / math.sqrt(n)
        # 반복 산포와 측정 불확도 중 큰 쪽 (보수적)
        sem = max(sem, float(np.sqrt(np.sum(np.asarray(uncs) ** 2))) / n)
    else:
        sd = float("nan")
        sem = float(uncs[0])
    k = tq(conf, n, use_z)
    return dict(n=n, mean=mean, sd=sd, sem=sem, low=mean - k * sem, high=mean + k * sem, k=k)


# ----------------------------------------------------------------------------- 계산
def compute(df: pd.DataFrame, a):
    log = []
    df = df.copy()

    # 측정 불확도 (입력 없으면 상대 정밀도 기본값)
    df["uCo"] = df["u_Co_ug"].fillna(a.rel_u_co * df["Co_ug"].abs())
    df["uNi"] = df["u_Ni_ug"].fillna(a.rel_u_ni * df["Ni_ug"].abs())

    # 용해 층 두께 [µm] = Δm / (ρ·A)
    df["layer_um"] = df["dissolved_mass_mg"] / 1000.0 / (a.rho * df["area_cm2"]) * 1e4

    # ---- 1. 공정 블랭크(C) → 모재 Co/Ni 비 r
    C = df[df.group == "C"]
    nC = len(C)
    if nC < 3:
        sys.exit(f"[오류] 공정 블랭크(C) 시편이 {nC}개입니다. 최소 3개(권장 5개 이상) 필요합니다.")
    if (C["Ni_ug"] <= 0).any():
        sys.exit("[오류] 블랭크의 Ni_ug 가 0 이하입니다.")
    ratios = (C["Co_ug"] / C["Ni_ug"]).to_numpy()
    r_bar = float(ratios.mean())
    s_r = float(ratios.std(ddof=1))
    if nC < 5:
        log.append(f"[주의] 블랭크 n={nC} (<5). 검출한계가 t분포 때문에 커집니다. 5개 이상 권장.")

    # ---- 2. 과잉 Co x = (Co − r̄·Ni)/A,  u(x) (GUM)
    df["x"] = (df["Co_ug"] - r_bar * df["Ni_ug"]) / df["area_cm2"]
    df["u_x"] = np.sqrt(df["uCo"] ** 2 + (r_bar * df["uNi"]) ** 2
                        + (df["Ni_ug"] * s_r) ** 2) / df["area_cm2"]

    # 블랭크는 leave-one-out r 로 다시 계산 (자기 자신으로 만든 r̄ 로 0이 되는 것 방지)
    loo = []
    for i, (_, row) in enumerate(C.iterrows()):
        r_i = float(np.delete(ratios, i).mean())
        loo.append((row["Co_ug"] - r_i * row["Ni_ug"]) / row["area_cm2"])
    loo = np.array(loo)
    df.loc[C.index, "x_loo"] = loo

    # ---- 3. σ0: (i) 모형 — 순신호 0일 때 전파 불확도, (ii) 경험 — 블랭크 LOO 산포
    #          경험값은 배경 Co(=r̄·Ni/A)에 비례하도록 상대값으로 환산해 각 시편에 적용
    bgC = r_bar * C["Ni_ug"] / C["area_cm2"]
    sigma_emp_rel = float(np.std(loo, ddof=1) / bgC.mean())
    bg = r_bar * df["Ni_ug"] / df["area_cm2"]
    uCo0 = df["u_Co_ug"].fillna(a.rel_u_co * r_bar * df["Ni_ug"])   # 순신호 0 이면 Co ≈ r̄·Ni
    df["sigma0_model"] = np.sqrt(uCo0 ** 2 + (r_bar * df["uNi"]) ** 2
                                 + (df["Ni_ug"] * s_r) ** 2) / df["area_cm2"]
    df["sigma0_emp"] = sigma_emp_rel * bg
    df["sigma0"] = df[["sigma0_model", "sigma0_emp"]].max(axis=1)

    k_a = tq(1 - a.alpha, nC, a.use_z)
    k_b = tq(1 - a.beta, nC, a.use_z)
    df["L_C"] = k_a * df["sigma0"]
    df["L_D"] = (k_a + k_b) * df["sigma0"]
    df["L_Q"] = a.kq * df["sigma0"]

    # ---- 4. 초기 crud Co, A0
    A = df[df.group == "A"]
    Lq = df[df.group == "L"]
    a0 = {}
    if len(A):
        a0["a"] = summarize(A["x"], A["u_x"], 1 - a.alpha, a.use_z)
        Ab = A[A["crud_mass_mg"].notna()]
        if len(Ab):
            vb = Ab["crud_mass_mg"] * a.w_co_ppm * 1e-3 / Ab["area_cm2"]      # mg·ppm·1e-3 = µg
            ub = vb * a.u_w_rel
            a0["b"] = summarize(vb, ub, 1 - a.alpha, a.use_z)
            df.loc[Ab.index, "A0_b"] = vb
    if len(Lq):
        a0["c"] = summarize(Lq["x"], Lq["u_x"], 1 - a.alpha, a.use_z)
    if not a0:
        sys.exit("[오류] A0 를 구할 그룹(A 또는 L)이 없습니다.")
    prim = a.a0_method if a.a0_method in a0 else next(m for m in "abc" if m in a0)
    if a.a0_method not in a0:
        log.append(f"[주의] A0 방법 '{a.a0_method}' 데이터가 없어 방법 '{prim}' 을 사용합니다.")
    A0 = a0[prim]

    # 방법 간 일치 (|Δ| ≤ 2·√(sem1²+sem2²))
    consist = []
    ms = list(a0.keys())
    for i in range(len(ms)):
        for j in range(i + 1, len(ms)):
            m1, m2 = a0[ms[i]], a0[ms[j]]
            z = abs(m1["mean"] - m2["mean"]) / math.hypot(m1["sem"], m2["sem"])
            consist.append((ms[i], ms[j], z, z <= 2.0))

    # ---- 5. 제염 시편(B) 판정과 DF
    B = df[df.group == "B"]
    if len(B) == 0:
        sys.exit("[오류] 제염 시편(B)이 없습니다.")
    if len(B) < 3:
        log.append(f"[주의] 제염 시편 n={len(B)} (<3). 3개 이상 권장.")
    df.loc[B.index, "detected"] = B["x"] > B["L_C"]
    det = df.loc[B.index, "detected"].astype(bool)
    LD_max = float(B["L_D"].max())
    xs = summarize(B["x"], B["u_x"], 1 - a.alpha, a.use_z)

    if not det.any():
        case = 1
        res_upper = LD_max
        DF_low = A0["low"] / res_upper
        DF_txt = f"DF ≥ {DF_low:.3g}"
        DF_val, DF_u = None, None
    elif (B["x"] >= B["L_Q"]).all():
        case = 3
        DF_val = A0["mean"] / xs["mean"]
        DF_u = DF_val * math.hypot(A0["sem"] / A0["mean"], xs["sem"] / xs["mean"])
        DF_low = A0["low"] / xs["high"]
        res_upper = xs["high"]
        DF_txt = f"DF = {DF_val:.3g} ± {DF_u:.2g} (표준불확도), 단측 하한 {DF_low:.3g}"
    else:
        case = 2
        res_upper = max(xs["high"], LD_max)
        DF_low = A0["low"] / res_upper
        DF_txt = f"DF ≥ {DF_low:.3g}"
        DF_val, DF_u = None, None

    # ---- 6. 검증 항목
    checks = []
    # 6-1 블랭크 LOO 평균 ≈ 0
    m_loo = float(loo.mean())
    lim = 2 * float(df.loc[C.index, "sigma0"].mean())
    checks.append(("공정 블랭크(C) 과잉 Co 평균 (LOO)", f"{m_loo:.3g} µg/cm² (허용 ±{lim:.3g})",
                   abs(m_loo) <= lim))
    # 6-2 스파이크 회수율
    D = df[df.group == "D"]
    rec = []
    for idx, row in D.iterrows():
        if pd.isna(row["spike_Co_ug"]) or row["spike_Co_ug"] <= 0:
            log.append(f"[주의] 스파이크 {row.sample_id}: spike_Co_ug 값이 없어 회수율 생략.")
            continue
        found = row["Co_ug"] - r_bar * row["Ni_ug"]
        R = 100 * found / row["spike_Co_ug"]
        df.loc[idx, "recovery_pct"] = R
        rec.append(R)
        spike_pa = row["spike_Co_ug"] / row["area_cm2"]
        if not (2 * row["L_D"] <= spike_pa <= 5 * row["L_D"]):
            log.append(f"[참고] 스파이크 {row.sample_id}: 첨가량 {spike_pa:.3g} µg/cm² 는 "
                       f"L_D의 2~5배({2*row['L_D']:.3g}~{5*row['L_D']:.3g}) 범위 밖입니다.")
    if rec:
        checks.append(("스파이크 회수율(D)", ", ".join(f"{r:.0f}%" for r in rec) + " (기준 85~115%)",
                       all(85 <= r <= 115 for r in rec)))
    else:
        checks.append(("스파이크 회수율(D)", "데이터 없음", None))
    # 6-3 양성대조
    E = df[df.group == "E"]
    if len(E):
        ok = bool((E["x"] > E["L_Q"]).all())
        checks.append(("양성대조(E) x > L_Q",
                       ", ".join(f"{x:.3g}" for x in E["x"]) + f" µg/cm² (L_Q ≈ {E['L_Q'].mean():.3g})", ok))
    else:
        checks.append(("양성대조(E)", "데이터 없음", None))
    # 6-4 A0 방법 일치
    for m1, m2, z, ok in consist:
        checks.append((f"A0 방법 {m1}·{m2} 일치", f"|Δ|/u = {z:.2f} (기준 ≤ 2)", ok))
    # 6-5 층 두께
    thick = df[df.group.isin(LAYER_GROUPS) & (df["layer_um"] > MAX_LAYER_UM)]
    lay = df[df.group.isin(LAYER_GROUPS)]["layer_um"].dropna()
    if len(lay):
        checks.append(("분석 층 두께 (B,C,D,E)", f"{lay.min():.2f}~{lay.max():.2f} µm (기준 ≤ {MAX_LAYER_UM:g})",
                       len(thick) == 0))
    else:
        checks.append(("분석 층 두께", "dissolved_mass_mg 미입력", None))

    return dict(df=df, r_bar=r_bar, s_r=s_r, nC=nC, k_a=k_a, k_b=k_b,
                sigma_emp_rel=sigma_emp_rel, a0=a0, prim=prim, A0=A0, xs=xs,
                case=case, DF_low=DF_low, DF_val=DF_val, DF_u=DF_u, DF_txt=DF_txt,
                res_upper=res_upper, LD_max=LD_max, checks=checks, log=log)


# ----------------------------------------------------------------------------- 출력
CASE_TXT = {1: "경우 1 — 모든 제염 시편 미검출 (x < L_C): DF ≥ A0,low / L_D",
            2: "경우 2 — 일부 검출, 정량 불확실 (L_C ≤ x < L_Q): DF ≥ A0,low / max(x_UCL, L_D)",
            3: "경우 3 — 정량 가능 (x ≥ L_Q): DF = A0 / x̄"}
A0_NAME = {"a": "방법 a (산화막 선택 용해)", "b": "방법 b (crud 질량 × Co 농도)",
           "c": "방법 c (제염액 물질수지, 확인용)"}


def f(v, d=3):
    return "—" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:.{d}g}"


def report(R, a, src: Path) -> str:
    df = R["df"]
    L = []
    L.append("# Co-crud 표면제염 DF 계산 결과\n")
    L.append(f"- 입력 파일: `{src.name}`")
    dist = "정규분포 z" if a.use_z else "t 분포 (자유도 %d)" % (R["nC"] - 1)
    L.append(f"- 유의수준 α = {a.alpha}, β = {a.beta}, 분위수: "
             f"{dist} → k_α = {R['k_a']:.3f}, k_β = {R['k_b']:.3f}")
    L.append(f"- 기본 상대 정밀도: Co {a.rel_u_co:.1%}, Ni {a.rel_u_ni:.1%} (u_Co_ug, u_Ni_ug 열이 비어 있을 때 사용)\n")

    L.append("## 1. 최종 결과\n")
    L.append(f"**{R['DF_txt']}** (신뢰수준 {1-a.alpha:.0%})\n")
    L.append(f"- 판정: {CASE_TXT[R['case']]}")
    A0 = R["A0"]
    L.append(f"- A0 ({A0_NAME[R['prim']]}): 평균 {A0['mean']:.4g} µg/cm², 표준불확도 {A0['sem']:.2g}, "
             f"n = {A0['n']}, 단측 하한 A0,low = {A0['low']:.4g} µg/cm²")
    L.append(f"- 잔류 Co 상한 (DF 계산에 사용): {R['res_upper']:.3g} µg/cm²")
    L.append(f"- 제염 시편 과잉 Co 평균: {R['xs']['mean']:.3g} ± {R['xs']['sem']:.2g} µg/cm² (n = {R['xs']['n']})\n")

    L.append("## 2. 모재 보정과 특성 한계\n")
    L.append(f"- 모재 Co/Ni 비 r̄ = {R['r_bar']:.5g}, s_r = {R['s_r']:.3g} "
             f"(상대 {R['s_r']/R['r_bar']:.2%}), 블랭크 n = {R['nC']}")
    L.append(f"- 블랭크 LOO 경험 σ0 / 배경 Co = {R['sigma_emp_rel']:.2%}")
    L.append("- σ0 = max(전파 모형값, 경험값). L_C = k_α·σ0, L_D = (k_α+k_β)·σ0, "
             f"L_Q = {a.kq:g}·σ0\n")

    L.append("## 3. 초기 crud Co (A0) 방법 비교\n")
    L.append("| 방법 | n | 평균 (µg/cm²) | 표준편차 | 표준불확도 | 단측 하한 |")
    L.append("|---|---|---|---|---|---|")
    for m, s in R["a0"].items():
        mark = " ← 사용" if m == R["prim"] else ""
        L.append(f"| {A0_NAME[m]}{mark} | {s['n']} | {s['mean']:.4g} | {f(s['sd'])} | {s['sem']:.2g} | {s['low']:.4g} |")
    L.append("")

    L.append("## 4. 시편별 결과\n")
    L.append("| 시편 | 그룹 | 층 두께 (µm) | x (µg/cm²) | u(x) | σ0 | L_C | L_D | L_Q | 판정 / 비고 |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    order = {g: i for i, g in enumerate("ALBCDEF")}
    for _, r in df.sort_values("group", key=lambda s: s.map(order)).iterrows():
        g = r.group
        note = ""
        if g == "B":
            note = "검출" if r.detected else "미검출 (x < L_C)"
            if r.x >= r.L_Q:
                note += ", 정량 가능"
        elif g == "C":
            note = f"LOO x = {r.x_loo:.3g}"
        elif g == "D" and not pd.isna(r.get("recovery_pct", np.nan)):
            note = f"회수율 {r.recovery_pct:.0f}%"
        elif g == "E":
            note = "x > L_Q ✓" if r.x > r.L_Q else "x ≤ L_Q ✗"
        elif g == "A" and not pd.isna(r.get("A0_b", np.nan)):
            note = f"방법 b = {r.A0_b:.3g}"
        show_lim = g in LAYER_GROUPS
        L.append(f"| {r.sample_id} | {g} | {f(r.layer_um)} | {r.x:.3g} | {r.u_x:.2g} | "
                 f"{f(r.sigma0) if show_lim else '—'} | {f(r.L_C) if show_lim else '—'} | "
                 f"{f(r.L_D) if show_lim else '—'} | {f(r.L_Q) if show_lim else '—'} | {note} |")
    L.append("")

    L.append("## 5. 방법 검증\n")
    L.append("| 항목 | 결과 | 판정 |")
    L.append("|---|---|---|")
    for name, val, ok in R["checks"]:
        L.append(f"| {name} | {val} | {'—' if ok is None else ('합격' if ok else '**불합격**')} |")
    L.append("")
    if R["log"]:
        L.append("## 6. 주의 사항\n")
        L += [f"- {m}" for m in R["log"]]
        L.append("")

    L.append("## 보고 문장 (초안 — 수치 확인 후 사용)\n")
    B = df[df.group == "B"]
    lay = B["layer_um"].dropna()
    lay_txt = f"표면 {lay.mean():.1f} µm 층의 " if len(lay) else ""
    if R["case"] == 1:
        L.append(f"> 제염 후 시편(n = {len(B)}) {lay_txt}과잉 Co는 {R['xs']['mean']:.2g} ± {R['xs']['sem']:.1g} µg/cm²로 "
                 f"결정 한계({B['L_C'].max():.2g} µg/cm²) 미만이었다(Currie/ISO 11929 기준 '미검출'). "
                 f"공정 블랭크(n = {R['nC']})로 구한 검출 한계 {R['LD_max']:.2g} µg/cm²와 "
                 f"제염 전 crud Co {A0['mean']:.3g} ± {A0['sem']:.1g} µg/cm²로부터 "
                 f"**DF ≥ {R['DF_low']:.3g} ({1-a.alpha:.0%} 신뢰수준)**이다.")
    else:
        L.append(f"> 제염 후 시편(n = {len(B)}) {lay_txt}과잉 Co는 {R['xs']['mean']:.3g} ± {R['xs']['sem']:.2g} µg/cm²였으며, "
                 f"제염 전 crud Co {A0['mean']:.3g} µg/cm² 대비 **{R['DF_txt']}**이다.")
    L.append("")
    return "\n".join(L)


def plot(R, path: Path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    df = R["df"]
    sub = df[df.group.isin(LAYER_GROUPS)].copy()
    colors = {"B": "#1f77b4", "C": "#7f7f7f", "D": "#2ca02c", "E": "#d62728"}
    labels = {"B": "B decontaminated", "C": "C process blank (LOO)", "D": "D spike (net)", "E": "E positive ctrl"}
    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=150)
    pos = 0
    ticks, tl = [], []
    for g in LAYER_GROUPS:
        s = sub[sub.group == g]
        for _, r in s.iterrows():
            y = r.x_loo if g == "C" else r.x
            ax.errorbar(pos, y, yerr=2 * r.u_x, fmt="o", color=colors[g], capsize=3,
                        label=labels[g] if r.sample_id == s.sample_id.iloc[0] else None)
            ticks.append(pos); tl.append(str(r.sample_id)); pos += 1
        pos += 0.6
    B = df[df.group == "B"]
    ax.axhline(0, color="k", lw=0.6)
    ax.axhline(B["L_C"].max(), color="orange", ls="--", lw=1, label=f"L_C = {B['L_C'].max():.3g}")
    ax.axhline(B["L_D"].max(), color="red", ls="--", lw=1, label=f"L_D = {B['L_D'].max():.3g}")
    ax.axhline(B["L_Q"].max(), color="purple", ls=":", lw=1, label=f"L_Q = {B['L_Q'].max():.3g}")
    ax.set_xticks(ticks); ax.set_xticklabels(tl, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("excess Co, x (µg/cm²)  [error bar: 2u]")
    ax.set_yscale("symlog", linthresh=max(B["L_C"].max(), 1e-6))
    ax.set_title(f"Residual Co after decontamination — {R['DF_txt'].split(',')[0]}", fontsize=9)
    ymax = float(np.nanmax(np.r_[sub["x"].to_numpy() + 2 * sub["u_x"].to_numpy(), B["L_Q"].max()]))
    ax.set_ylim(top=ymax * 2.5)
    ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return True


def main(argv=None):
    p = argparse.ArgumentParser(description="Co-crud 표면제염 DF 계산기 (Currie/ISO 11929)")
    p.add_argument("csv", type=Path, help="입력 CSV (templates/template.csv 형식)")
    p.add_argument("--out", type=Path, default=Path("results"), help="출력 폴더 (기본 results)")
    p.add_argument("--alpha", type=float, default=0.05, help="1종 오류 확률 (기본 0.05)")
    p.add_argument("--beta", type=float, default=0.05, help="2종 오류 확률 (기본 0.05)")
    p.add_argument("--kq", type=float, default=10.0, help="정량한계 계수 L_Q = kq·σ0 (기본 10)")
    p.add_argument("--rel-u-co", type=float, default=0.02, help="Co 상대 표준불확도 기본값 (기본 0.02)")
    p.add_argument("--rel-u-ni", type=float, default=0.02, help="Ni 상대 표준불확도 기본값 (기본 0.02)")
    p.add_argument("--w-co-ppm", type=float, default=5000.0, help="crud 내 Co 농도 ppm (방법 b, 기본 5000)")
    p.add_argument("--u-w-rel", type=float, default=0.05, help="crud Co 농도의 상대 불확도 (방법 b, 기본 0.05)")
    p.add_argument("--rho", type=float, default=7.9, help="모재 밀도 g/cm³ (층 두께 계산, 기본 7.9)")
    p.add_argument("--a0-method", choices=["a", "b", "c"], default="a", help="DF에 쓸 A0 방법 (기본 a)")
    p.add_argument("--use-z", action="store_true", help="t 대신 정규분포 z 사용 (블랭크가 많을 때)")
    a = p.parse_args(argv)

    if not a.csv.exists():
        sys.exit(f"[오류] 파일이 없습니다: {a.csv}")
    df = load_data(a.csv)
    R = compute(df, a)
    a.out.mkdir(parents=True, exist_ok=True)
    md = report(R, a, a.csv)
    (a.out / "report.md").write_text(md, encoding="utf-8")
    cols = ["sample_id", "group", "area_cm2", "layer_um", "Co_ug", "Ni_ug", "uCo", "uNi", "x", "u_x",
            "x_loo", "sigma0", "L_C", "L_D", "L_Q", "detected", "recovery_pct", "A0_b", "note"]
    out = R["df"].reindex(columns=cols)
    out.to_csv(a.out / "results.csv", index=False, encoding="utf-8-sig", float_format="%.6g")
    has_plot = plot(R, a.out / "plot.png")

    print(md)
    print(f"\n[저장] {a.out/'report.md'}, {a.out/'results.csv'}" + (f", {a.out/'plot.png'}" if has_plot else ""))
    return R


if __name__ == "__main__":
    main()
