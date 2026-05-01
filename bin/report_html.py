#!/usr/bin/env python3
import argparse, html, base64, json
from pathlib import Path
import pandas as pd
import datetime
from bs4 import BeautifulSoup

def read_run_report_meta(run_report_html: Path):
    if not run_report_html or not run_report_html.exists():
        return {}
    soup = BeautifulSoup(run_report_html.read_text(errors="ignore"), "html.parser")
    lines=[ln.strip() for ln in soup.get_text("\n").splitlines()]
    def val_after(key):
        for i,ln in enumerate(lines):
            if ln == key:
                for j in range(i+1, min(i+15, len(lines))):
                    if lines[j]:
                        return lines[j]
        return ""
    meta = {}
    meta["Flow cell type"] = val_after("Flow cell type")
    meta["Flow cell ID"] = val_after("Flow cell ID")
    meta["Kit type"] = val_after("Kit type")
    meta["Expansion kit"] = val_after("Expansion kit")
    meta["Basecalling"] = val_after("Basecalling")
    return {k:v for k,v in meta.items() if v}

def read_key_value_file(path: Path):
    if not path or not path.exists():
        return []
    items = []
    for line in path.read_text(errors="ignore").splitlines():
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        items.append((k.strip(), v.strip()))
    return items

def kv_items_to_html(items):
    if not items:
        return ""
    rows = "".join(
        f"<tr><th>{html.escape(str(k))}</th><td>{html.escape(str(v))}</td></tr>"
        for k, v in items
    )
    return f'<table class="kv-table"><tbody>{rows}</tbody></table>'

def hide_composition_columns(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return df
    def is_hidden(col: str) -> bool:
        key = str(col).strip().lower().replace("-", "").replace("_", "").replace(" ", "")
        return key in {"nopv", "nopv2"}
    keep = [c for c in df.columns if not is_hidden(c)]
    return df[keep].copy()

def is_fasta_path(v: str) -> bool:
    v = (v or "").strip().lower()
    return v.endswith(".fasta") or v.endswith(".fa") or v.endswith(".fna")

def df_to_html_table(df, table_id, base_dir: Path, embed_fastas: bool = True):
    cols = df.columns.tolist()

    # Cache for embedded FASTA content to keep HTML generation fast
    _cache = {}  # rel_or_abs_path -> base64 text

    # columns that represent per-group mapped reads in the composition table
    GROUP_COUNT_COLS = set([
        "NonPolioEV","WPV1","WPV2","WPV3","Sabin1-related","Sabin2-related","Sabin3-related","nOPV"
    ])

    def is_read_count_col(colname: str) -> bool:
        c = colname.strip().lower()
        if c == "unmapped" or "unmapped" in c:
            return False
        if "reads" in c:
            return True
        if colname in GROUP_COUNT_COLS:
            return True
        return False

    def fmt_mut(raw):
        s = str(raw).strip()
        if s == "" or s.lower() == "nan":
            return ""
        try:
            return str(int(float(s)))
        except Exception:
            return s

    def _as_path(v: str) -> Path:
        p = Path(str(v))
        return p if p.is_absolute() else (base_dir / p)

    def _base64_for_file(p: Path) -> str:
        key = str(p)
        if key in _cache:
            return _cache[key]
        try:
            b = p.read_bytes()
        except Exception:
            _cache[key] = ""
            return ""
        b64 = base64.b64encode(b).decode("ascii")
        _cache[key] = b64
        return b64

    def _badge_for_group(g: str) -> str:
        gs = (g or '').strip()
        if not gs:
            return ''
        key = gs.lower()
        cls = 'badge badge-generic'
        if key == 'wpv1':
            cls = 'badge badge-wpv1'
        elif key == 'wpv2':
            cls = 'badge badge-wpv2'
        elif key == 'wpv3':
            cls = 'badge badge-wpv3'
        elif key.startswith('sabin1'):
            cls = 'badge badge-sabin1'
        elif key.startswith('sabin2'):
            cls = 'badge badge-sabin2'
        elif key.startswith('sabin3'):
            cls = 'badge badge-sabin3'
        elif key.startswith('nopv') or key.startswith('nopv'.replace('pv','pv')) or key.startswith('nopv'):
            cls = 'badge badge-nopv'
        elif key.startswith('nonpolioev') or key.startswith('non polio'):
            cls = 'badge badge-npe'
        return "<span class='%s'>%s</span>" % (cls, html.escape(gs))

    def _pill_for_yesno(v) -> str:
        s = str(v).strip().lower()
        if s in ('yes','y','true','1','pass','passed'):
            return "<span class='pill pill-yes'>Yes</span>"
        if s in ('no','n','false','0','fail','failed'):
            return "<span class='pill pill-no'>No</span>"
        if s == '' or s == 'nan':
            return ""
        return "<span class='pill pill-unk'>%s</span>" % html.escape(str(v))



    rows = []
    for _, r in df.iterrows():
        tds = []
        for c in cols:
            v = r.get(c, "")
            if pd.isna(v):
                v = ""
            raw = v
            s = str(v)

            # FASTA links open in a new tab so the sequence can be viewed directly.
            if s.strip() and (c.lower().endswith("fasta") or is_fasta_path(s) or c.strip().lower() == "sequence (vp1)"):
                p = _as_path(s)
                if embed_fastas and p.exists() and p.is_file():
                    b64 = _base64_for_file(p)
                    if b64:
                        cell_html = f"<a href='#' onclick='return openEmbeddedFasta(event, {json.dumps(b64)}, {json.dumps(p.name)});'>FASTA</a>"
                    else:
                        cell_html = f'<a href="{html.escape(s)}" target="_blank" rel="noopener noreferrer">FASTA</a>'
                else:
                    cell_html = f'<a href="{html.escape(s)}" target="_blank" rel="noopener noreferrer">FASTA</a>'
            else:
                # Mutations integer formatting
                if c.strip().lower() == "mutations":
                    cell_html = html.escape(fmt_mut(raw))
                else:
                    cell_html = html.escape(s)

            # Bold sample names everywhere
            if c.strip().lower() == "sample":
                cell_html = f"<b>{cell_html}</b>"

            # Bold read counts > 50 (mapped categories) in all tables, excluding Unmapped
            if is_read_count_col(c):
                try:
                    n = int(float(str(raw).strip())) if str(raw).strip() != "" else 0
                    if n > 50:
                        cell_html = f"<b>{cell_html}</b>"
                except Exception:
                    pass

            # Color badges for key categorical columns
            if c.strip().lower() == "reference group":
                cell_html = _badge_for_group(str(raw))
            if c.strip().lower() == "consensus recovered":
                cell_html = _pill_for_yesno(raw)
            tds.append(f"<td>{cell_html}</td>")
        rows.append("<tr>" + "".join(tds) + "</tr>")

    thead = "<thead><tr>" + "".join([f"<th>{html.escape(c)}</th>" for c in cols]) + "</tr></thead>"
    tbody = "<tbody>" + "".join(rows) + "</tbody>"
    return f'<table class="sortable" id="{table_id}">{thead}{tbody}</table>'

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--run-report-html", default="")
    ap.add_argument("--title", default="MinION DDNS Sequencing Report", help="Report title shown in HTML")
    ap.add_argument("--run-name", default="", help="Run name shown under the title")
    ap.add_argument("--embed-fastas", type=int, default=1, help="1 to embed FASTA downloads inside HTML (self-contained). 0 to keep external links.")
    ap.add_argument("--signal-threshold", type=int, default=50, help="WPV1 signal mapped reads threshold used for display (default 50).")
    args = ap.parse_args()

    title = args.title
    run_name = (args.run_name or '').strip()
    report_date = datetime.date.today().strftime('%d/%m/%Y')
    run_line = f"<div class=\"muted small\"><b>Run</b>: {html.escape(run_name)}</div>" if run_name else ""

    out = Path(args.out)
    meta = read_run_report_meta(Path(args.run_report_html)) if args.run_report_html else {}
    versions_txt = out / "00_logs" / "versions.txt"
    kv_items = read_key_value_file(versions_txt)
    version_keys = {"python", "minimap2", "samtools", "cutadapt", "filtlong", "medaka"}
    param_items = [(k, v) for k, v in kv_items if k not in version_keys]
    tool_items = [(k, v) for k, v in kv_items if k in version_keys]

    info_csv = out / "Table_Sample summary information.csv"
    comp_csv = out / "Table_Composition of samples.csv"
    top_csv  = out / "Table_Top hits.csv"
    sumgrp_csv = out / "Summary_Unique_Samples_By_Group.csv"
    confgrp_csv = out / "Summary_Confirmed_By_Group.csv"
    runsum_csv = out / "Summary_Run_Overview.csv"
    sum_tsv  = out / "summary.tsv"
    wpv1_sig_csv = out / "Table_WPV1_signal.csv"
    wpv1_conf_csv = out / "Table_WPV1_confirmed.csv"
    conf_all_csv = out / "Table_Confirmed_All.csv"
    conf_low_csv = out / "Table_Confirmed_LowReads.csv"


    # Load "confirmed low reads" table safely
    try:
        if conf_low_csv.exists() and conf_low_csv.stat().st_size > 0:
            df_conf_low = pd.read_csv(conf_low_csv)
        else:
            df_conf_low = pd.DataFrame()
    except Exception:
        df_conf_low = pd.DataFrame()

    def safe_read_csv(p, **kw):
        try:
            if not p.exists():
                return pd.DataFrame()
            try:
                if p.stat().st_size == 0:
                    return pd.DataFrame()
            except Exception:
                pass
            return pd.read_csv(p, **kw)
        except pd.errors.EmptyDataError:
            return pd.DataFrame()
        except Exception:
            return pd.DataFrame()

    df_info = safe_read_csv(info_csv)
    df_comp = safe_read_csv(comp_csv)
    df_top  = safe_read_csv(top_csv)
    df_sumgrp = safe_read_csv(sumgrp_csv)
    df_confgrp = safe_read_csv(confgrp_csv)
    df_runsum = safe_read_csv(runsum_csv)
    df_sum  = safe_read_csv(sum_tsv, sep="\t")
    df_wpv1_sig = safe_read_csv(wpv1_sig_csv)
    df_wpv1_conf = safe_read_csv(wpv1_conf_csv)
    df_conf_all = safe_read_csv(conf_all_csv)
    df_comp = hide_composition_columns(df_comp)


    total_samples = len(df_comp) if not df_comp.empty else 0
    pass_count = int((df_sum["status"]=="PASS").sum()) if (not df_sum.empty and "status" in df_sum.columns) else 0

    meta_html = ""
    params_html = ""

    summary_html = ""
    wpv1_thr = int(getattr(args, 'signal_threshold', 50))

    summary_blocks = []
    order = ["WPV1","WPV2","WPV3","Sabin1-related","Sabin2-related","Sabin3-related","nOPV","NonPolioEV"]
    if not df_confgrp.empty and "Group" in df_confgrp.columns and "Confirmed_Barcodes" in df_confgrp.columns:
        cards = []
        m = {str(r["Group"]): int(r["Confirmed_Barcodes"]) for _, r in df_confgrp.iterrows()}
        for g in order:
            cards.append(f"""
              <div class="kpi">
                <div class="kpi-label">{html.escape(g)}</div>
                <div class="kpi-value">{m.get(g,0)}</div>
              </div>
            """)
        summary_blocks.append(
            '<div><h2>Summary</h2><div class="muted small">Confirmed barcode-level detections by reference group. One barcode can count in multiple groups, but only once per group.</div><div class="kpi-grid">' + ''.join(cards) + '</div></div>'
        )
    elif not df_sumgrp.empty and "Group" in df_sumgrp.columns and "Unique_Samples" in df_sumgrp.columns:
        cards = []
        m = {str(r["Group"]): int(r["Unique_Samples"]) for _, r in df_sumgrp.iterrows()}
        for g in order:
            cards.append(f"""
              <div class="kpi">
                <div class="kpi-label">{html.escape(g)}</div>
                <div class="kpi-value">{m.get(g,0)}</div>
              </div>
            """)
        summary_blocks.append('<h2>Summary</h2><div class="kpi-grid">' + ''.join(cards) + '</div>')

    if not df_runsum.empty and "Metric" in df_runsum.columns and "Value" in df_runsum.columns:
        metric_order = ["Total Samples", "Any Confirmed", "Mixed Confirmed", "No Confirmed Consensus"]
        cards = []
        m = {str(r["Metric"]): int(r["Value"]) for _, r in df_runsum.iterrows()}
        for metric in metric_order:
            cards.append(f"""
              <div class="kpi">
                <div class="kpi-label">{html.escape(metric)}</div>
                <div class="kpi-value">{m.get(metric,0)}</div>
              </div>
            """)
        summary_blocks.append(
            '<div><div class="muted small" style="margin:6px 0 8px;">Run summary</div><div class="kpi-grid">' + ''.join(cards) + '</div></div>'
        )

    summary_html = "".join(summary_blocks)
    if meta:
        items = "".join([f"<li><b>{html.escape(k)}</b>: {html.escape(v)}</li>" for k,v in meta.items()])
        meta_html = f"<h2>Run metadata</h2><ul>{items}</ul>"

    if param_items or tool_items:
        blocks = []
        if param_items:
            blocks.append(
                "<div class=\"subcard\"><h3 class=\"subhead\" style=\"margin-top:0;\">Run settings</h3>" +
                kv_items_to_html(param_items) +
                "</div>"
            )
        if tool_items:
            blocks.append(
                "<div class=\"subcard\"><h3 class=\"subhead\" style=\"margin-top:0;\">Tool versions</h3>" +
                kv_items_to_html(tool_items) +
                "</div>"
            )
        params_html = "<h2>Parameters used</h2><div class=\"kv-grid\">" + "".join(blocks) + "</div>"

    style = r'''
    body { font-family: Arial, sans-serif; margin: 18px; color: #111; }
    h1 { margin-bottom: 6px; }
    .muted { color: #555; }
    .subhead { margin: 12px 0 6px; font-size: 14px; }
    .grid { display: grid; grid-template-columns: 1fr; gap: 16px; }
    .card { border: 1px solid #ddd; border-radius: 12px; padding: 14px; background: #fff; }
    table { width: 100%; border-collapse: collapse; font-size: 14px; }
    th, td { border-bottom: 1px solid #eee; padding: 8px; text-align: left; }
    th { position: sticky; top: 0; background: #fafafa; }
    .toolbar { display: flex; gap: 10px; align-items: center; margin: 10px 0 14px; flex-wrap: wrap; }
    .btn { padding: 7px 10px; border: 1px solid #ccc; border-radius: 10px; background: #fafafa; cursor: pointer; }
    .btn:hover { background: #f0f0f0; }
    input { padding: 8px; border: 1px solid #ccc; border-radius: 10px; width: 340px; }
    .small { font-size: 13px; }
    .header { display:flex; align-items:flex-start; justify-content:space-between; gap:12px; }
    .header-left { display:flex; flex-direction:column; }
    .datebox { font-size:13px; color:#555; margin-top: 8px; white-space: nowrap; }
    th { cursor: pointer; }
    th.sorted-asc::after { content: ' ▲'; font-size: 11px; color:#666; }
    th.sorted-desc::after { content: ' ▼'; font-size: 11px; color:#666; }
    /* color theme */

    :root {
      --bg: #f6f8fb;
      --card: #ffffff;
      --text: #0f172a;
      --muted: #64748b;
      --border: #e7eaf0;
      --shadow: 0 10px 26px rgba(15, 23, 42, 0.08);
      --blue: #2563eb;
      --purple: #7c3aed;
      --orange: #f97316;
      --green: #16a34a;
      --teal: #0891b2;
      --slate: #334155;
    }

    body { background: var(--bg); color: var(--text); }
    a { color: var(--blue); text-decoration: none; }
    a:hover { text-decoration: underline; }

    .card {
      border: 1px solid var(--border);
      border-radius: 14px;
      padding: 14px;
      background: var(--card);
      box-shadow: var(--shadow);
    }

    .section-info { border-top: 5px solid var(--blue); }
    .section-comp { border-top: 5px solid var(--purple); }
    .section-top  { border-top: 5px solid var(--orange); }

    .subcard {
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 12px;
      background: #fbfdff;
      margin-top: 12px;
    }
    .kv-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 12px;
      margin: 12px 0 18px;
    }
    .kv-table th, .kv-table td {
      padding: 7px 9px;
      vertical-align: top;
    }
    .kv-table th {
      width: 38%;
      background: #f8fafc;
      font-weight: 700;
    }
    .sub-wpv1sig { border-left: 6px solid var(--orange); }
    .sub-wpv1conf { border-left: 6px solid var(--green); }
    .sub-confall { border-left: 6px solid var(--teal); }

    table {
      width: 100%;
      border-collapse: separate;
      border-spacing: 0;
      font-size: 14px;
      border: 1px solid var(--border);
      border-radius: 12px;
      overflow: hidden;
      background: #fff;
    }
    thead th {
      background: #f1f5f9;
      color: var(--slate);
      font-weight: 700;
      border-bottom: 1px solid var(--border);
    }
    tbody tr:nth-child(even) { background: #f8fafc; }
    tbody tr:hover { background: #eef2ff; }

    .btn { background: #ffffff; border: 1px solid var(--border); border-radius: 10px; }
    .btn:hover { background: #f1f5f9; }
    input { background: #ffffff; border: 1px solid var(--border); }

    .kpi {
      border: 1px solid var(--border);
      border-radius: 16px;
      background: #ffffff;
      box-shadow: 0 8px 18px rgba(15, 23, 42, 0.06);
    }

    .badge {
      display: inline-flex;
      align-items: center;
      padding: 2px 10px;
      border-radius: 999px;
      border: 1px solid var(--border);
      font-size: 12px;
      font-weight: 700;
      line-height: 18px;
      white-space: nowrap;
    }
    .badge-generic { background: #f1f5f9; color: #0f172a; }
    .badge-wpv1 { background: #fee2e2; color: #7f1d1d; border-color: #fecaca; }
    .badge-wpv2 { background: #ffedd5; color: #7c2d12; border-color: #fed7aa; }
    .badge-wpv3 { background: #e0f2fe; color: #0c4a6e; border-color: #bae6fd; }
    .badge-sabin1 { background: #dcfce7; color: #14532d; border-color: #bbf7d0; }
    .badge-sabin2 { background: #dbeafe; color: #1e3a8a; border-color: #bfdbfe; }
    .badge-sabin3 { background: #f3e8ff; color: #581c87; border-color: #e9d5ff; }
    .badge-nopv { background: #cffafe; color: #164e63; border-color: #a5f3fc; }
    .badge-npe { background: #e2e8f0; color: #0f172a; border-color: #cbd5e1; }

    .pill {
      display: inline-flex;
      align-items: center;
      padding: 2px 10px;
      border-radius: 999px;
      font-size: 12px;
      font-weight: 800;
      border: 1px solid var(--border);
      line-height: 18px;
    }
    .pill-yes { background: #dcfce7; color: #14532d; border-color: #bbf7d0; }
    .pill-no  { background: #e2e8f0; color: #334155; border-color: #cbd5e1; }
    .pill-unk { background: #f1f5f9; color: #334155; border-color: #e2e8f0; }

    .dot { display:inline-block; width:10px; height:10px; border-radius:999px; margin-right:10px; vertical-align:middle; }
    .dot-info { background: var(--blue); }
    .dot-comp { background: var(--purple); }
    .dot-top { background: var(--orange); }

    .datebox {
      padding: 8px 10px;
      border: 1px solid var(--border);
      border-radius: 12px;
      background: #ffffff;
      box-shadow: 0 8px 18px rgba(15, 23, 42, 0.06);
    }
    .muted { color: var(--muted); }

'''

    js = r'''
function getCellText(td) {
  if (!td) return "";
  return (td.innerText || td.textContent || "").trim();
}

function isNumericColumn(rows, colIndex) {
  let seen = 0;
  for (const r of rows) {
    const t = getCellText(r.children[colIndex]);
    if (t === "") continue;
    seen += 1;
    const n = Number(t.replace(/,/g, ""));
    if (!Number.isFinite(n)) return false;
  }
  return seen > 0;
}

function sortTableByColumn(table, colIndex, direction) {
  const tbody = table.tBodies[0];
  const rows = Array.from(tbody.rows);
  const numeric = isNumericColumn(rows, colIndex);

  rows.sort((a, b) => {
    const av = getCellText(a.children[colIndex]);
    const bv = getCellText(b.children[colIndex]);
    if (numeric) {
      const an = av === "" ? Number.NEGATIVE_INFINITY : Number(av.replace(/,/g,""));
      const bn = bv === "" ? Number.NEGATIVE_INFINITY : Number(bv.replace(/,/g,""));
      return direction === "asc" ? (an - bn) : (bn - an);
    } else {
      const cmp = av.localeCompare(bv, undefined, {numeric:true, sensitivity:"base"});
      return direction === "asc" ? cmp : -cmp;
    }
  });

  rows.forEach(r => tbody.appendChild(r));
}

function enableSorting() {
  document.querySelectorAll("table.sortable").forEach(table => {
    const ths = table.querySelectorAll("thead th");
    ths.forEach((th, idx) => {
      th.addEventListener("click", () => {
        const allTh = table.querySelectorAll("thead th");
        allTh.forEach(h => { h.classList.remove("sorted-asc"); h.classList.remove("sorted-desc"); });

        const current = th.getAttribute("data-sort") || "none";
        const next = (current === "asc") ? "desc" : "asc";
        th.setAttribute("data-sort", next);
        th.classList.add(next === "asc" ? "sorted-asc" : "sorted-desc");

        sortTableByColumn(table, idx, next);
      });
    });
  });
}

document.addEventListener("DOMContentLoaded", enableSorting);
    function filterTables() {
      const q = document.getElementById('q').value.toLowerCase();
      const tables = document.querySelectorAll('table');
      tables.forEach(t => {
        const rows = t.querySelectorAll('tbody tr');
        rows.forEach(r => {
          const txt = r.innerText.toLowerCase();
          r.style.display = txt.includes(q) ? '' : 'none';
        });
      });
    }

        function tableToCSV(tableId) {
          const table = document.getElementById(tableId);
          const rows = Array.from(table.querySelectorAll('tr'));
          return rows.map(row => {
            const cells = Array.from(row.querySelectorAll('th,td'));
            return cells.map(cell => {
              const txt = cell.innerText.replace(/\r?\n/g, ' ').trim();
              return '"' + txt.replace(/"/g, '""') + '"';
            }).join(',');
          }).join('\n');
        }

        async function copyCSV(tableId) {
          const csv = tableToCSV(tableId);
          try {
            await navigator.clipboard.writeText(csv);
            alert('Copied CSV to clipboard');
          } catch(e) {
            alert('Copy failed (browser permissions).');
          }
        }

        function downloadCSV(tableId, filename) {
          const csv = tableToCSV(tableId);
          const blob = new Blob([csv], {type: 'text/csv;charset=utf-8;'});
          const url = URL.createObjectURL(blob);
          const a = document.createElement('a');
          a.href = url;
          a.download = filename;
          document.body.appendChild(a);
          a.click();
          document.body.removeChild(a);
          URL.revokeObjectURL(url);
        }

        function escapeHtmlText(s) {
          return String(s)
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;');
        }

        function openEmbeddedFasta(event, b64, filename) {
          if (event) event.preventDefault();
          try {
            const text = atob(b64);
            const win = window.open('', '_blank');
            if (!win) {
              alert('Browser blocked the FASTA tab. Please allow pop-ups for this report.');
              return false;
            }
            const title = filename || 'consensus.fasta';
            win.document.open();
            win.document.write(`<!doctype html><html><head><meta charset="utf-8"><title>${escapeHtmlText(title)}</title><style>body{font-family:Consolas,Menlo,monospace;margin:16px;background:#fff;color:#111}pre{white-space:pre-wrap;word-break:break-word}</style></head><body><pre>${escapeHtmlText(text)}</pre></body></html>`);
            win.document.close();
          } catch (e) {
            console.error(e);
            alert('Could not open the FASTA in a new tab.');
          }
          return false;
        }
    '''

    info_table = df_to_html_table(df_info, "info", out, bool(args.embed_fastas)) if not df_info.empty else "<p>No sample summary table found.</p>"
    wpv1_sig_table = df_to_html_table(df_wpv1_sig, "wpv1sig", out, bool(args.embed_fastas)) if not df_wpv1_sig.empty else "<p class=\"muted\">No samples met the WPV1 signal threshold.</p>"
    wpv1_conf_table = df_to_html_table(df_wpv1_conf, "wpv1conf", out, bool(args.embed_fastas)) if not df_wpv1_conf.empty else "<p class=\"muted\">No samples with recovered WPV1 consensus.</p>"
    conf_all_table  = df_to_html_table(df_conf_all,  "confall", out, bool(args.embed_fastas)) if not df_conf_all.empty else "<p class=\"muted\">No confirmed consensus found.</p>"
    conf_low_table  = df_to_html_table(df_conf_low,  "conflow", out, bool(args.embed_fastas)) if not df_conf_low.empty else "<p class=\"muted\">No recovered consensus with mapped reads below 50.</p>"

    comp_table = df_to_html_table(df_comp, "comp", out, bool(args.embed_fastas)) if not df_comp.empty else "<p>No composition table found.</p>"
    top_table  = df_to_html_table(df_top,  "top",  out, bool(args.embed_fastas)) if not df_top.empty else "<p>No top hits table found.</p>"

    html_doc = f'''<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>{style}.kpi-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 12px; margin: 12px 0 18px; }}
.kpi {{ border: 1px solid #eee; border-radius: 16px; padding: 14px 14px; background: #fff; display:flex; align-items:center; justify-content:space-between; gap:12px; }}
.kpi-label {{ font-size: 13px; color: #666; margin: 0; }}
.kpi-value {{ font-size: 34px; font-weight: 700; line-height: 1; margin: 0; }}
</style>
</head>
<body>
<div class="header">
  <div class="header-left">
    <h1 style="margin:0;">{title}</h1>
    {run_line}
    <div class="muted small">Output folder: {html.escape(str(out))} | Samples: {total_samples}</div>
  </div>
  <div class="datebox">{report_date}</div>
</div>

{meta_html}

{summary_html}

<div class="toolbar">
  <input id="q" type="text" placeholder="Search any table" oninput="filterTables()">
  <div class="muted small">Search any value</div>
</div>

<div class="grid">
<div class="card">
  <div style="display:flex; gap:10px; align-items:center; justify-content:space-between; flex-wrap:wrap;">
    <h2 style="margin:0;"><span class="dot dot-info"></span>Sample summary information</h2>
    <div>
      <button class="btn" onclick="copyCSV('info')">Copy CSV</button>
      <button class="btn" onclick="downloadCSV('info','Table_Sample_summary_information.csv')">Save as CSV</button>
    </div>
  </div>
  <div class="muted small">One row per sample per reference group. Each row shows the top finding in that group, ranked by mapped reads and then supporting metrics. Breadth is calculated using the active min_cov threshold.</div>
  {info_table}
  <hr style="border:none; border-top:1px solid #eee; margin:14px 0;">

  <div style="display:flex; gap:10px; align-items:center; justify-content:space-between; flex-wrap:wrap; margin-top:12px;">
    <h3 class="subhead" style="margin:0;">Confirmed (all consensus recovered)</h3>
    <div>
      <button class="btn" onclick="copyCSV('confall')">Copy CSV</button>
      <button class="btn" onclick="downloadCSV('confall','Table_Confirmed_All.csv')">Save as CSV</button>
    </div>
  </div>
  {conf_all_table}

</div>
<div class="card">
  <div style="display:flex; gap:10px; align-items:center; justify-content:space-between; flex-wrap:wrap;">
    <h2 style="margin:0;"><span class="dot dot-comp"></span>Composition of samples</h2>
    <div>
      <button class="btn" onclick="copyCSV('comp')">Copy CSV</button>
      <button class="btn" onclick="downloadCSV('comp','Table_Composition_of_samples.csv')">Save as CSV</button>
    </div>
  </div>
  {comp_table}
</div>
<div class="card">
  <div style="display:flex; gap:10px; align-items:center; justify-content:space-between; flex-wrap:wrap;">
    <h2 style="margin:0;"><span class="dot dot-top"></span>Reported hits</h2>
    <div>
      <button class="btn" onclick="copyCSV('top')">Copy CSV</button>
      <button class="btn" onclick="downloadCSV('top','Table_Top_hits.csv')">Save as CSV</button>
    </div>
  </div>
  <div class="muted small">All reportable hits are listed here. FASTA links are embedded in this HTML when embed-fastas is enabled and open in a new tab.</div>
  {top_table}
</div>
</div>

{params_html}

<script>{js}</script>
</body>
</html>'''

    (out / "report.html").write_text(html_doc)

if __name__ == "__main__":
    main()
