"""Self-contained HTML report for a proteome-wide scan."""
from __future__ import annotations

import base64
import datetime as dt
import html
import json
import re
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from . import plots
from .ranking import TIERS

CSS = """
body{margin:0;background:#fcfcfb;color:#0b0b0b;font:15px/1.6 -apple-system,'Segoe UI',Roboto,Arial,sans-serif}
main{max-width:980px;margin:0 auto;padding:32px 20px 64px}
h1{font-size:1.7rem;margin:0 0 4px;letter-spacing:-.01em}h2{font-size:1.2rem;margin:40px 0 8px;border-bottom:1px solid #e7e6e2;padding-bottom:6px}
.sub{color:#52514e;margin:0 0 24px}.muted{color:#52514e;font-size:.9rem}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:18px 0}
.tile{border:1px solid #e7e6e2;border-radius:10px;padding:12px 14px;background:#fff}.tile b{display:block;font-size:1.6rem;line-height:1.1}
.tile span{color:#52514e;font-size:.82rem}
table{border-collapse:collapse;width:100%;font-size:.84rem;margin:8px 0}th,td{padding:6px 8px;border-bottom:1px solid #eeede9;text-align:right;white-space:nowrap}
th{color:#52514e;font-weight:600;background:#f5f5f2;position:sticky;top:0}td:first-child,th:first-child,td:nth-child(2),th:nth-child(2){text-align:left}
.scroll{overflow-x:auto;max-height:560px;border:1px solid #e7e6e2;border-radius:8px}
img{max-width:100%;border:1px solid #e7e6e2;border-radius:8px;background:#fcfcfb}
code{background:#f1f0ec;padding:1px 5px;border-radius:4px;font-size:.85em}
li{margin:4px 0}
"""


def _img(path):
    return f'<img alt="" src="data:image/png;base64,{base64.b64encode(Path(path).read_bytes()).decode()}">'


def _fmt(v, kind):
    if pd.isna(v):
        return ""
    if kind == "p":
        return f"{v:.2e}" if v < 0.001 else f"{v:.3f}"
    if kind == "f2":
        return f"{v:.2f}"
    return html.escape(str(v))


def _table(d: pd.DataFrame):
    cols = [("tier", "Tier", "s"), ("gene", "Protein", "s"), ("mr_OR", "OR", "f2"), ("mr_OR_lci", "95% CI low", "f2"),
            ("mr_OR_uci", "95% CI high", "f2"), ("mr_p", "MR p", "p"), ("mr_q", "MR q (FDR)", "p"),
            ("n_instruments", "Instr.", "s"), ("coloc_PP.H4", "PP.H4", "f2"), ("coloc_PP.H3", "PP.H3", "f2"),
            ("flags", "Flags", "s")]
    cols = [c for c in cols if c[0] in d]
    head = "".join(f"<th>{h}</th>" for _, h, _ in cols)
    body = "".join("<tr>" + "".join(f"<td>{_fmt(r[c], k)}</td>" for c, _, k in cols) + "</tr>" for _, r in d.iterrows())
    return f'<div class="scroll"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def build_report(results_path, region_dir, out_path, highlight=(), params_path=None, n_regional=6):
    res = pd.read_csv(results_path, sep="\t")
    tmp = Path(tempfile.mkdtemp())
    region_dir = Path(region_dir)
    tested = res.mr_p.notna().sum()
    counts = res.tier.value_counts()

    plots.volcano(res, tmp / "volcano.png", highlight=highlight)
    has_forest = plots.forest(res, tmp / "forest.png")

    show = res[res.tier.isin(TIERS[:3])].head(n_regional)
    show = pd.concat([show, res[res.gene.isin(highlight)]]).drop_duplicates("gene")
    regional_html = []
    for _, r in show.iterrows():
        key = re.sub(r"[^A-Za-z0-9_.-]", "_", r.protein_id)
        reg = region_dir / f"{key}.region.tsv.gz"
        if not reg.exists():
            continue
        region = pd.read_csv(reg, sep="\t")
        inst_f = region_dir / f"{key}.instruments.tsv"
        inst = pd.read_csv(inst_f, sep="\t") if inst_f.exists() and inst_f.stat().st_size > 1 else pd.DataFrame()
        title = f"{r.gene} \u2014 " + (f"Tier {r.tier}" if r.tier[0].isdigit() else r.tier)
        sub = (f"MR OR {r.mr_OR:.2f} (p = {_fmt(r.mr_p, 'p')}) \u00b7 coloc PP.H4 = {r['coloc_PP.H4']:.2f}, "
               f"PP.H3 = {r['coloc_PP.H3']:.2f}")
        f = tmp / f"reg_{key}.png"
        plots.regional(region, inst, title, f, subtitle=sub)
        regional_html.append(_img(f))

    params = json.loads(Path(params_path).read_text()) if params_path else {}
    param_html = "".join(f"<li><code>{html.escape(k)}</code> = {html.escape(str(v))}</li>" for k, v in params.items())
    top = res[res.tier.isin(TIERS[:4]) | res.gene.isin(highlight)]

    tiles = [(len(res), "proteins in scan"), (tested, "with a valid cis instrument"),
             (int(counts.get(TIERS[0], 0)), "Tier 1: MR + coloc"), (int(counts.get(TIERS[1], 0)), "Tier 2: suggestive coloc"),
             (int(counts.get(TIERS[2], 0)), "Tier 3: distinct variants")]
    tiles_html = "".join(f'<div class="tile"><b>{n}</b><span>{t}</span></div>' for n, t in tiles)

    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>PD proteome MR scan</title><style>{CSS}</style></head>
<body><main>
<h1>Proteome-wide cis-MR and colocalisation scan for Parkinson’s disease</h1>
<p class="sub">Generated {dt.date.today().isoformat()} · {len(res)} proteins · full table: <code>scan_results.tsv</code></p>
<div class="tiles">{tiles_html}</div>

<h2>How to read the tiers</h2>
<ul>
<li><b>Tier 1</b> — MR significant after FDR correction <i>and</i> colocalisation PP.H4 ≥ 0.8: protein and PD most likely share one causal variant. Strongest genetic support.</li>
<li><b>Tier 2</b> — MR significant, PP.H4 0.5–0.8: supportive but not conclusive.</li>
<li><b>Tier 3</b> — MR significant but PP.H3 ≥ 0.5: the protein and PD signals come from <i>different</i> variants, so the MR estimate is probably driven by LD rather than the protein.</li>
<li><b>Tier 4</b> — MR significant, colocalisation inconclusive (often an underpowered outcome signal).</li>
</ul>

<h2>MR estimates across the proteome</h2>
{_img(tmp / "volcano.png")}
<p class="muted">Each point is one protein’s primary MR estimate (Wald ratio for one instrument, multiplicative random-effects IVW for several). Colour and marker shape show the evidence tier.</p>

{"<h2>Prioritised proteins</h2>" + _img(tmp / "forest.png") if has_forest else ""}

<h2>Results table</h2>
{_table(top) if len(top) else "<p>No protein reached FDR significance.</p>"}

<h2>Regional association plots</h2>
<p class="muted">Top: cis-pQTL association; bottom: Parkinson’s disease GWAS over the same variants. Blue points are the MR instruments.</p>
{"".join(regional_html) or "<p>None to show.</p>"}

<h2>Methods</h2>
<ul>
<li>cis window: gene body ± window (GRCh38 from the Olink protein map); variants joined to the outcome on GRCh37 position and alleles.</li>
<li>QC: INFO, MAF, allele harmonisation (strand flips resolved; palindromic variants kept only when MAF &lt; 0.42 and allele frequencies agree).</li>
<li>Instruments: p &lt; threshold and F ≥ threshold, LD-clumped against the reference panel; if the panel lacks the variants, the single lead variant is used.</li>
<li>MR: Wald ratio (1 instrument) or multiplicative random-effects IVW (≥2), with MR-Egger and weighted median as sensitivity (≥3). Benjamini–Hochberg FDR across all tested proteins.</li>
<li>Colocalisation: approximate Bayes factor method (single causal variant per trait), priors p1, p2, p12 as below; PP.H4 also reported at p12 = 1e-6 and 1e-4.</li>
</ul>
{"<p class='muted'>Run parameters:</p><ul>" + param_html + "</ul>" if param_html else ""}

<h2>Limitations</h2>
<ul>
<li>Single-causal-variant colocalisation can report H3 at loci with several independent pQTLs even when one of them is shared with PD. Follow up loci flagged <code>multiple_independent_signals</code> or <code>H3_at_PD_locus_check_multi_signal</code> with SuSiE-based colocalisation (coloc.susie) using an LD matrix.</li>
<li>cis-MR assumes the instrument acts on PD only through the protein; protein-altering variants can change assay binding (epitope effects) rather than protein level.</li>
<li>Results are genetic prioritisation evidence, not proof of causality or druggability.</li>
</ul>
</main></body></html>"""
    Path(out_path).write_text(doc)
    return out_path
