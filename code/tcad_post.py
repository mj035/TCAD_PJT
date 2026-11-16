"""tcad_post.py - Sentaurus Process 결과(.plx / .csv) 후처리 및 발표용 그래프 생성

사용법 (서버에서 내려받은 파일이 있는 폴더에서 실행):

  python tcad_post.py plx  n12_t0.plx n12_t60.plx n12_t240.plx --field NetActive --labels "spike 직후,60 min,240 min" --out A_profile.png
  python tcad_post.py winA n*_A.csv  --dxj 10 --drs 10 --out A_window.png      (레시피가 여럿이면 레시피별 파일)
  python tcad_post.py remA n*_A.csv  --tpost 800 --out A_remedy.png            (Part 2: 경계 온도에서 레시피 비교)
  python tcad_post.py rsxj n*_B.csv  --xj 20 --rs 600  --out B_rsxj.png
  python tcad_post.py oxC  n*_C.csv  --out C_orient.png

공통 원칙
  - 농도 축은 log, 범례는 그래프 바깥(프로파일을 가리지 않음)
  - 접합 깊이 Xj 를 점으로 찍고 범례에 숫자로 표기
  - 그림과 같은 이름의 .csv 표를 함께 저장 (발표 자료 표로 그대로 사용)
"""
import argparse
import csv
import glob
import math
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker  # noqa: F401
import matplotlib.colors  # noqa: F401
import matplotlib.font_manager  # noqa: F401

# ---------------------------------------------------------------- 스타일
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SEQ = ["#86b6ef", "#5598e7", "#2a78d6", "#256abf", "#184f95", "#0d366b"]  # 순서형(시간·온도)


def setup():
    for font in ("Malgun Gothic", "AppleGothic", "NanumGothic"):
        if any(f.name == font for f in matplotlib.font_manager.fontManager.ttflist):
            plt.rcParams["font.family"] = font
            break
    plt.rcParams.update({
        "axes.unicode_minus": False, "mathtext.fontset": "dejavusans",
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "grid.linestyle": "-",
        "axes.spines.top": False, "axes.spines.right": False,
        "font.size": 11, "axes.titlesize": 13, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "lines.linewidth": 2.0, "legend.frameon": False,
    })


def seq_colors(n):
    if n <= 1:
        return [SEQ[2]]
    idx = np.linspace(0, len(SEQ) - 1, n).round().astype(int)
    return [SEQ[i] for i in idx]


NOTE = ""


def save(fig, out, rows=None, header=None):
    if NOTE:
        fig.text(0.01, -0.02, NOTE, ha="left", va="top", fontsize=8.5, color=MUTED, transform=fig.transFigure)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print("saved", out)
    if rows:
        tab = os.path.splitext(out)[0] + ".csv"
        with open(tab, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.writer(f)
            if header:
                w.writerow(header)
            w.writerows(rows)
        print("saved", tab)


# ---------------------------------------------------------------- PLX 읽기
def read_plx(path):
    """{곡선이름: (x[um], y)} 반환. 따옴표 줄 = 곡선 이름, 숫자 두 개 줄 = 데이터."""
    curves, name, xs, ys = {}, None, [], []

    def flush():
        if name is not None and xs:
            curves[name] = (np.array(xs), np.array(ys))

    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.replace(",", " ").split()
            try:
                if len(parts) < 2:
                    raise ValueError
                x, y = float(parts[0]), float(parts[1])
                xs.append(x)
                ys.append(y)
            except ValueError:
                flush()
                name, xs, ys = s.strip('"').strip(), [], []
    flush()
    return curves


def pick(curves, key):
    for k in curves:
        if k.lower() == key.lower():
            return curves[k]
    for k in curves:
        if k.lower().startswith(key.lower()):
            return curves[k]
    return None


def junction_um(x, net):
    """NetActive 부호가 바뀌는 가장 깊은 위치. 부호가 없으면(절댓값 저장) 가장 깊은 국소 최솟값."""
    sgn = np.sign(net)
    idx = np.where(sgn[:-1] * sgn[1:] < 0)[0]
    if len(idx):
        i = idx[-1]
        return float(x[i] + (x[i + 1] - x[i]) * (0 - net[i]) / (net[i + 1] - net[i]))
    a = np.log10(np.abs(net) + 1.0)
    cand = [i for i in range(1, len(a) - 1) if a[i] < a[i - 1] and a[i] < a[i + 1] and a[i] < np.median(a) - 0.3]
    return float(x[cand[-1]]) if cand else float("nan")


def mob_p(n):   # Masetti (boron), cm2/Vs
    n = np.maximum(n, 1.0)
    return 44.9 * np.exp(-9.23e16 / n) + (470.5 - 44.9) / (1 + (n / 2.23e17) ** 0.719) - 29.0 / (1 + (6.1e20 / n) ** 2.0)


def mob_n(n):   # Masetti (phosphorus), cm2/Vs
    n = np.maximum(n, 1.0)
    return 68.5 + (1414 - 68.5) / (1 + (n / 9.2e16) ** 0.711) - 56.1 / (1 + (3.41e20 / n) ** 1.98)


def sheet_res(x, net, xj, ptype_top=True):
    """Rs = 1 / (q * int(mu * |net|) dx), 표면~Xj 구간. x[um], net[cm-3]."""
    if not np.isfinite(xj):
        return float("nan")
    m = (x >= 0) & (x <= xj)
    if m.sum() < 2:
        return float("nan")
    c = np.abs(net[m])
    mu = mob_p(c) if ptype_top else mob_n(c)
    g = 1.602e-19 * np.trapezoid(mu * c, x[m] * 1e-4)
    return 1.0 / g if g > 0 else float("nan")


# ---------------------------------------------------------------- plx : 프로파일 겹쳐 그리기
def cmd_plx(a):
    files = [f for pat in a.files for f in sorted(glob.glob(pat))] or a.files
    labels = a.labels.split(",") if a.labels else [os.path.splitext(os.path.basename(f))[0] for f in files]
    colors = seq_colors(len(files)) if a.ordered else CAT
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    rows = []
    for i, (f, lab) in enumerate(zip(files, labels)):
        cur = read_plx(f)
        d = pick(cur, a.field)
        if d is None:
            print("!! field", a.field, "not in", f, "->", list(cur))
            continue
        x, y = d
        col = colors[i % len(colors)]
        net = pick(cur, "NetActive")
        xj = junction_um(*net) if net is not None else float("nan")
        rs = sheet_res(net[0], net[1], xj, not a.ntop) if net is not None else float("nan")
        txt = lab if not np.isfinite(xj) else f"{lab}  (Xj {xj*1000:.1f} nm)"
        ax.plot(x * 1000, np.abs(y) + 1.0, color=col, label=txt)
        if np.isfinite(xj):
            # 접합 위치를 x축 위에 점으로 표시 (곡선이 급락하는 지점이라 축 위가 가장 읽기 쉬움)
            ax.plot([xj * 1000], [0.0], "o", ms=9, color=col, mec=SURFACE, mew=2, zorder=6,
                    transform=ax.get_xaxis_transform(), clip_on=False)
        dose = float(np.trapezoid(np.abs(y)[x >= 0], x[x >= 0] * 1e-4))
        rows.append([lab, os.path.basename(f), f"{xj*1000:.2f}", f"{rs:.1f}", f"{dose:.3e}", f"{np.abs(y).max():.3e}"])
    ax.set_yscale("log")
    ax.set_xlabel("깊이 x (nm)")
    ax.set_ylabel(f"{a.field} (cm$^{{-3}}$)")
    if a.xmax:
        ax.set_xlim(0 if not a.xmin else a.xmin, a.xmax)
    if a.ymin or a.ymax:
        ax.set_ylim(a.ymin or None, a.ymax or None)
    ax.set_title(a.title or f"{a.field} 프로파일")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), title="조건 (점 = 접합 위치)", alignment="left")
    save(fig, a.out, rows, ["label", "file", "Xj_nm", "Rs_ohm_sq(Masetti)", "dose_cm2", "peak_cm3"])


# ---------------------------------------------------------------- 공통 CSV 읽기
def read_rows(patterns):
    rows = []
    for pat in patterns:
        for f in sorted(glob.glob(pat)) or [pat]:
            with open(f, encoding="utf-8-sig") as fh:
                rows += list(csv.DictReader(fh))
    return rows


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


# ---------------------------------------------------------------- winA : thermal budget window
def recipe_key(r):
    pai = "PAI" if r.get("pai", "0").strip() == "1" else "noPAI"
    c = "C" + r.get("cdose", "0") if fnum(r.get("cdose", "0")) > 0 else "noC"
    return f'{pai}_{c}_E{r.get("Eimp_keV", "")}_spk{r.get("Tspike_C", "")}'


def cmd_winA(a):
    allrows = [r for r in read_rows(a.files) if fnum(r["t_min"]) > 0]
    recipes = sorted({recipe_key(r) for r in allrows})
    for rk in recipes:
        rows = [r for r in allrows if recipe_key(r) == rk]
        out = a.out if len(recipes) == 1 else os.path.splitext(a.out)[0] + "_" + rk + os.path.splitext(a.out)[1]
        print("recipe", rk, "->", out)
        _winA_one(rows, a, out)


def _winA_one(rows, a, out):
    Ts = sorted({fnum(r["Tpost_C"]) for r in rows})
    ts = sorted({fnum(r["t_min"]) for r in rows})
    Z = np.full((len(Ts), len(ts)), np.nan)
    R = np.full_like(Z, np.nan)
    for r in rows:
        i, j = Ts.index(fnum(r["Tpost_C"])), ts.index(fnum(r["t_min"]))
        Z[i, j], R[i, j] = fnum(r["dXj_pct"]), fnum(r["dRs_pct"])
    ok = (np.abs(Z) <= a.dxj) & (np.abs(R) <= a.drs)

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8), gridspec_kw={"width_ratios": [1, 1, 0.9]},
                             layout="constrained")
    fig.get_layout_engine().set(w_pad=0.25)
    for ax, M, name, lim in ((axes[0], Z, "ΔXj (%)", a.dxj), (axes[1], R, "ΔRs (%)", a.drs)):
        A = np.abs(M)
        vmax = np.nanmax(A) if np.isfinite(A).any() else 1
        im = ax.imshow(A, origin="lower", aspect="auto", cmap=matplotlib.colors.LinearSegmentedColormap.from_list(
            "seq", ["#f0efec", "#9ec5f4", "#2a78d6", "#0d366b"]), vmin=0, vmax=max(vmax, lim * 2))
        for i in range(len(Ts)):
            for j in range(len(ts)):
                if np.isfinite(M[i, j]):
                    dark = A[i, j] > 0.55 * max(vmax, lim * 2)
                    ax.text(j, i, f"{M[i, j]:+.1f}", ha="center", va="center", fontsize=10,
                            color="#ffffff" if dark else INK, fontweight="bold" if A[i, j] > lim else "normal")
        ax.set_xticks(range(len(ts)), [f"{t:g}" for t in ts])
        ax.set_yticks(range(len(Ts)), [f"{T:g}" for T in Ts])
        ax.set_xlabel("누적 열처리 시간 (min)")
        ax.set_ylabel("후속 열공정 온도 (°C)")
        ax.set_title(f"{name}  ·  기준 |Δ| ≤ {lim:g}% (굵은 숫자 = 초과)")
        ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(False)

    # 허용 시간 t_allow(T): |ΔXj| 가 기준에 도달하는 시간 (log-log 보간) -> Arrhenius
    ax = axes[2]
    k = 8.617e-5
    pts = []
    for i, T in enumerate(Ts):
        z = np.abs(Z[i])
        m = np.isfinite(z) & (z > 0)
        if m.sum() < 2 or z[m].max() < a.dxj or z[m].min() > a.dxj:
            continue
        t_allow = float(np.exp(np.interp(math.log(a.dxj), np.log(z[m]), np.log(np.array(ts)[m]))))
        pts.append((T, t_allow))
    table = [[f"{T:g}"] + [f"{Z[i, j]:.2f}" for j in range(len(ts))] + [f"{R[i, j]:.2f}" for j in range(len(ts))]
             + ["".join("O" if ok[i, j] else "X" for j in range(len(ts)))] for i, T in enumerate(Ts)]
    if len(pts) >= 2:
        T_, t_ = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
        inv = 1.0 / (k * (T_ + 273.15))
        slope, icpt = np.polyfit(inv, np.log(t_), 1)
        ax.plot(inv, t_, "o", ms=8, color=CAT[0], mec=SURFACE, mew=2, zorder=5, label="ΔXj 기준 도달 시간")
        xx = np.linspace(inv.min(), inv.max(), 50)
        ax.plot(xx, np.exp(icpt + slope * xx), color=CAT[0], lw=2, label=f"Arrhenius 맞춤  Ea = {slope:.2f} eV")
        ax.plot(xx, np.exp(np.log(t_[0]) + 3.46 * (xx - inv[0])), color=MUTED, lw=1.5,
                label="이론 기울기 (B 고유 확산 3.46 eV)")
        for T, t in pts:
            ax.annotate(f"{T:g}°C\n{t:.0f} min", (1.0 / (k * (T + 273.15)), t), textcoords="offset points",
                        xytext=(8, 6), fontsize=9, color=INK2)
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.margins(x=0.18, y=0.25)
        ax.set_xlabel("1/kT (eV$^{-1}$)   ← 고온   저온 →")
        ax.set_ylabel(f"허용 시간 (min), |ΔXj| = {a.dxj:g}%")
        ax.set_title("허용 thermal budget")
        ax.legend(loc="upper left", fontsize=9)
        print(f"Ea_eff = {slope:.3f} eV ; t_allow:", ", ".join(f"{T:g}C={t:.1f}min" for T, t in pts))
    else:
        ax.text(0.5, 0.5, "허용 시간 점이 2개 미만\n(온도 split 을 늘리거나 기준 % 조정)", ha="center", va="center",
                transform=ax.transAxes, color=INK2)
        ax.set_axis_off()
    save(fig, out, table, ["Tpost_C"] + [f"dXj%_{t:g}min" for t in ts] + [f"dRs%_{t:g}min" for t in ts] + ["window(O/X)"])


# ---------------------------------------------------------------- remA : 처방 비교 (Part 2)
def cmd_remA(a):
    rows = [r for r in read_rows(a.files) if fnum(r["t_min"]) > 0 and fnum(r["Tpost_C"]) == a.tpost]
    if not rows:
        sys.exit(f"Tpost={a.tpost} 인 행이 없습니다. --tpost 값을 CSV의 Tpost_C 값 중 하나로 주세요.")
    recipes = sorted({recipe_key(r) for r in rows})
    base = a.base or [k for k in recipes if k.startswith("noPAI_noC")]
    base = base[0] if isinstance(base, list) and base else (base if isinstance(base, str) else recipes[0])
    recipes = [base] + [k for k in recipes if k != base]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), layout="constrained")
    table = []
    for i, rk in enumerate(recipes):
        rs = sorted([r for r in rows if recipe_key(r) == rk], key=lambda r: fnum(r["t_min"]))
        t = [fnum(r["t_min"]) for r in rs]
        col = CAT[i % 8]
        lab = "기준 (baseline)" if rk == base else rk
        for ax, key in ((axes[0], "dXj_pct"), (axes[1], "dRs_pct")):
            y = [fnum(r[key]) for r in rs]
            ax.plot(t, y, "-o", color=col, ms=7, mec=SURFACE, mew=1.5, label=lab)
        brow = {fnum(r["t_min"]): r for r in rows if recipe_key(r) == base}
        for r in rs:
            b = brow.get(fnum(r["t_min"]))
            def imp(key):
                if not b or fnum(b[key]) == 0:
                    return ""
                return f"{(1 - abs(fnum(r[key])) / abs(fnum(b[key]))) * 100:.0f}"
            table.append([rk, a.tpost, r["t_min"], r["dXj_pct"], imp("dXj_pct"), r["dRs_pct"], imp("dRs_pct")])
    for ax, name in ((axes[0], "ΔXj (%)"), (axes[1], "ΔRs (%)")):
        ax.axhline(a.crit, color=MUTED, lw=1)
        ax.axhline(-a.crit, color=MUTED, lw=1)
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
        ax.set_xlabel("누적 열처리 시간 (min)")
        ax.set_ylabel(name)
        ax.set_title(f"{name} @ {a.tpost:g} °C  (회색선 = ±{a.crit:g} % 기준)")
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), title="레시피", alignment="left")
    save(fig, a.out, table, ["recipe", "Tpost_C", "t_min", "dXj_pct", "Xj_improve_%_vs_base", "dRs_pct", "Rs_improve_%_vs_base"])


# ---------------------------------------------------------------- rsxj : Rs-Xj trade-off
def cmd_rsxj(a):
    rows = read_rows(a.files)
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    groups, table = {}, []
    for r in rows:
        key = ("Ge PAI" if r["pai"].strip() == "1" else "PAI 없음") + (" + C" if fnum(r["cdose"]) > 0 else "")
        groups.setdefault(key, []).append(r)
    for gi, (key, rs) in enumerate(groups.items()):
        x = [fnum(r["Xj_nm"]) for r in rs]
        y = [fnum(r["Rs_ohm"]) for r in rs]
        ax.plot(x, y, "o", ms=9, color=CAT[gi % 8], mec=SURFACE, mew=2, label=key, zorder=5)
        for r, xi, yi in zip(rs, x, y):
            ax.annotate(f'{r["rta"]} {r["Tpk_C"]}°C', (xi, yi), textcoords="offset points", xytext=(8, 5),
                        fontsize=8.5, color=INK2)
            table.append([key, r["rta"], r["Tpk_C"], r["EB_keV"], r["cdose"], r["Xj_nm"], r["Rs_ohm"],
                          r["abrupt_nm_dec"], r["active_pct"],
                          "O" if (xi <= a.xj and yi <= a.rs) else "X"])
    xl, yl = ax.get_xlim(), ax.get_ylim()
    ax.fill_between([xl[0], a.xj], yl[0], a.rs, color="#cde2fb", alpha=0.55, zorder=1, lw=0)
    ax.annotate(f"목표 영역\nXj ≤ {a.xj:g} nm\nRs ≤ {a.rs:g} Ω/sq", (xl[0], yl[0]), textcoords="offset points",
                xytext=(8, 8), va="bottom", ha="left", fontsize=9.5, color=INK2)
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_xlabel("접합 깊이 Xj (nm)")
    ax.set_ylabel("면저항 Rs (Ω/sq)")
    ax.set_title("Rs–Xj trade-off (왼쪽 아래일수록 좋음)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), title="주입 조건", alignment="left")
    save(fig, a.out, table, ["group", "rta", "Tpk_C", "EB_keV", "cdose", "Xj_nm", "Rs_ohm", "abrupt_nm_dec", "active_pct", "target"])


# ---------------------------------------------------------------- oxC : 결정면별 tox
def cmd_oxC(a):
    rows = read_rows(a.files)
    ambs = [x for x in ("O2", "H2O", "ISSG") if any(r["amb"] == x for r in rows)]
    oris = [x for x in ("100", "110", "111") if any(r["orient"] == x for r in rows)]
    name = {"O2": "Dry O$_2$", "H2O": "Wet H$_2$O", "ISSG": "ISSG (라디칼)"}
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    w, table = 0.26, []
    for oi, o in enumerate(oris):
        xs, ys = [], []
        for ai, am in enumerate(ambs):
            sel = [r for r in rows if r["amb"] == am and r["orient"] == o]
            if not sel:
                continue
            r = sel[-1]
            t = fnum(r["tox_nm"])
            base = [q for q in rows if q["amb"] == am and q["orient"] == "100"]
            ratio = t / fnum(base[-1]["tox_nm"]) if base else float("nan")
            xs.append(ai + (oi - (len(oris) - 1) / 2) * w)
            ys.append(t)
            ax.text(xs[-1], t, f"{t:.2f}\n(×{ratio:.2f})", ha="center", va="bottom", fontsize=9, color=INK2)
            table.append([am, o, r["Tox_C"], r["t_min"], f"{t:.3f}", f"{ratio:.3f}",
                          r.get("Cs_cm3", ""), r.get("X17_nm", ""), r.get("Bret_pct", "")])
        ax.bar(xs, ys, width=w - 0.02, color=CAT[oi], label=f"({o})", zorder=3)
    ax.set_xticks(range(len(ambs)), [name[x] for x in ambs])
    ax.set_ylabel("산화막 두께 tox (nm)")
    ax.set_title("같은 레시피를 결정면만 바꿔 적용했을 때 tox (괄호 = (100) 대비 비율)")
    ax.grid(axis="x", visible=False)
    ax.margins(y=0.18)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), title="결정면", alignment="left")
    save(fig, a.out, table, ["amb", "orient", "Tox_C", "t_min", "tox_nm", "ratio_vs_100", "Cs_cm3", "X17_nm", "Bret_pct"])


def main():
    setup()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("plx", help="plx 프로파일 겹쳐 그리기 + Xj/Rs 표")
    s.add_argument("files", nargs="+")
    s.add_argument("--field", default="NetActive")
    s.add_argument("--labels", default="")
    s.add_argument("--title", default="")
    s.add_argument("--xmin", type=float, default=0)
    s.add_argument("--xmax", type=float, default=0, help="nm")
    s.add_argument("--ymin", type=float, default=0)
    s.add_argument("--ymax", type=float, default=0)
    s.add_argument("--categorical", dest="ordered", action="store_false", help="순서 없는 조건 비교(색상 구분)")
    s.add_argument("--ntop", action="store_true", help="표면층이 n형일 때")
    s.add_argument("--out", default="profile.png")
    s.set_defaults(func=cmd_plx)

    s = sub.add_parser("winA", help="주제 A: thermal budget window")
    s.add_argument("files", nargs="+")
    s.add_argument("--dxj", type=float, default=10.0, help="허용 |ΔXj| %")
    s.add_argument("--drs", type=float, default=10.0, help="허용 |ΔRs| %")
    s.add_argument("--out", default="A_window.png")
    s.set_defaults(func=cmd_winA)

    s = sub.add_parser("remA", help="주제 A Part 2: 한 온도에서 레시피별 ΔXj·ΔRs 비교와 개선율")
    s.add_argument("files", nargs="+")
    s.add_argument("--tpost", type=float, required=True, help="비교할 후속 열처리 온도 (window 경계 온도)")
    s.add_argument("--crit", type=float, default=10.0)
    s.add_argument("--base", default="", help="기준 레시피 키 (기본: noPAI_noC_... 자동)")
    s.add_argument("--out", default="A_remedy.png")
    s.set_defaults(func=cmd_remA)

    s = sub.add_parser("rsxj", help="주제 B: Rs-Xj trade-off")
    s.add_argument("files", nargs="+")
    s.add_argument("--xj", type=float, default=20.0)
    s.add_argument("--rs", type=float, default=600.0)
    s.add_argument("--out", default="B_rsxj.png")
    s.set_defaults(func=cmd_rsxj)

    s = sub.add_parser("oxC", help="주제 C: 결정면별 tox")
    s.add_argument("files", nargs="+")
    s.add_argument("--out", default="C_orient.png")
    s.set_defaults(func=cmd_oxC)

    p.add_argument("--note", default="", help="그림 아래에 작은 회색 글씨로 넣을 주석 (예: 예비 모델 표시)")
    a = p.parse_args()
    global NOTE
    NOTE = a.note
    a.func(a)


if __name__ == "__main__":
    sys.exit(main())
