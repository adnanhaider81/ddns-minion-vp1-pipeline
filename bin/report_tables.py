#!/usr/bin/env python3
import argparse, csv, re
from pathlib import Path
import pandas as pd

def fmt_int(v):
    """Format numeric-like values as integer strings (e.g. 13.0 -> 13)."""
    try:
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return ""
    except Exception:
        pass
    s = str(v).strip()
    if s == "" or s.lower() == "nan":
        return ""
    try:
        return str(int(float(s)))
    except Exception:
        return s

COMP_COLUMNS = ["Sample","Barcode","Unmapped","NonPolioEV","WPV1","WPV2","WPV3","Sabin1-related","Sabin2-related","Sabin3-related","nOPV"]
INFO_COLUMNS = ["Sample","Barcode","Hit Rank","Ref Id","Reference Group","Mapped Reads","Breadth","Mean Depth","Consensus Recovered","Mutations","Number of Ns","Consensus FASTA"]

def fmt_pct(v):
    """Format fractional breadth values as percentages (e.g. 1.0 -> 100.00%)."""
    try:
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return ""
    except Exception:
        pass
    s = str(v).strip()
    if s == "" or s.lower() == "nan":
        return ""
    try:
        return f"{float(s) * 100:.2f}%"
    except Exception:
        return s

def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", (s or "").strip())

def normalize_barcode(value: str) -> str:
    s = (value or "").strip()
    if not s:
        return ""
    sl = s.lower()
    if sl in {"unclassified", "uncl"}:
        return "unclassified"
    m = re.search(r"(\d+)", s)
    if m:
        return f"barcode{int(m.group(1)):02d}"
    return s

def barcode_sort_key(value: str):
    m = re.match(r"barcode(\d+)$", value or "")
    if m:
        return (0, int(m.group(1)))
    if value == "unclassified":
        return (1, 0)
    return (2, value or "")

def parse_ref_headers(ref_headers_tsv: Path):
    ref_to_group = {}
    ref_to_cluster = {}
    with ref_headers_tsv.open() as f:
        for line in f:
            line=line.rstrip("\n")
            if not line:
                continue
            ref_id, header = line.split("\t", 1)
            m = re.search(r"ddns_group=([^\s]+)", header)
            ref_to_group[ref_id] = m.group(1) if m else "Other"
            m2 = re.search(r"cluster=([^\s]+)", header)
            ref_to_cluster[ref_id] = m2.group(1) if m2 else ref_id
    return ref_to_group, ref_to_cluster

def read_barcode_map_csv(barcodes_csv: Path):
    out = {}
    with barcodes_csv.open(newline="") as f:
        r=csv.DictReader(f)
        if not r.fieldnames or "barcode" not in r.fieldnames or "sample" not in r.fieldnames:
            raise SystemExit("barcodes.csv must have columns: barcode,sample")
        for row in r:
            b=normalize_barcode(str(row["barcode"]))
            s=str(row["sample"]).strip()
            if b:
                out[b] = s
    return out

def list_barcode_folders(fastq_pass: Path):
    bcs = []
    for p in fastq_pass.glob("barcode*"):
        if p.is_dir() and re.match(r"^barcode\d+$", p.name):
            bcs.append(p.name)
    def key(x):
        m=re.match(r"barcode(\d+)$", x)
        return int(m.group(1)) if m else 10**9
    return sorted(set(bcs), key=key)

def read_idxstats_full(p: Path):
    unmapped = 0
    refs = []
    if not p.exists():
        return unmapped, refs
    with p.open() as f:
        for line in f:
            line=line.rstrip("\n")
            if not line:
                continue
            parts=line.split("\t")
            if len(parts) < 4:
                continue
            ref_id, ref_len, mapped, unm = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
            if ref_id == "*":
                unmapped = unm
            else:
                refs.append((ref_id, ref_len, mapped))
    return unmapped, refs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--fastq-pass", required=True)
    ap.add_argument("--barcode-map", required=True)
    ap.add_argument("--min-breadth", type=float, default=0.85)
    ap.add_argument("--max-hits", type=int, default=3)
    ap.add_argument("--min-cons-len", type=int, default=800)
    ap.add_argument("--max-n-frac", type=float, default=0.30)
    ap.add_argument("--signal-threshold", type=int, default=50, help="Mapped reads threshold for WPV1 signal table (default 50).")
    ap.add_argument("--min-group-reads", type=int, default=50, help="In Sample summary information, keep only the best reference per group with mapped reads >= this threshold")
    args = ap.parse_args()

    out = Path(args.out)
    fastq_pass = Path(args.fastq_pass)
    barcodes_csv = Path(args.barcode_map)

    ref_headers = out / "ref" / "ref_headers.tsv"
    if not ref_headers.exists():
        raise SystemExit(f"Missing {ref_headers}")

    ref_to_group, ref_to_cluster = parse_ref_headers(ref_headers)
    barcode_map = read_barcode_map_csv(barcodes_csv)
    folders = list_barcode_folders(fastq_pass)

    if barcode_map:
        all_barcodes = sorted(set(barcode_map.keys()), key=barcode_sort_key)
    else:
        all_barcodes = sorted(set(folders), key=barcode_sort_key)

    summary_tsv = out / "summary.tsv"
    df_sum = pd.read_csv(summary_tsv, sep="\t") if summary_tsv.exists() else pd.DataFrame()
    sum_map = {str(r["sample_id"]): r for _, r in df_sum.iterrows()} if not df_sum.empty else {}

    hits_tsv = out / "top_hits_report.tsv"
    df_hits = pd.read_csv(hits_tsv, sep="\t") if hits_tsv.exists() else pd.DataFrame()

    cluster_to_seqid = {}
    seq_counter = 1

    comp_rows=[]
    info_rows=[]
    top_rows=[]

    for barcode_folder in all_barcodes:
        sample_name = barcode_map.get(barcode_folder, barcode_folder)
        safe_sample = safe_name(sample_name)
        sample_id = f"{safe_sample}__{barcode_folder}"

        reads_filtered = None
        best_ref = ""
        if sample_id in sum_map:
            r = sum_map[sample_id]
            try:
                reads_filtered = int(r["reads_filtered"])
            except Exception:
                reads_filtered = None
            best_ref = str(r.get("best_ref","") or "")

        sample_dir = out / "per_sample" / sample_id
        idx_full = sample_dir / "idxstats.full.tsv"
        unmapped_reads, refs = read_idxstats_full(idx_full)

        counts = {"NonPolioEV":0, "WPV1":0, "WPV2":0, "WPV3":0, "Sabin1-related":0, "Sabin2-related":0, "Sabin3-related":0, "nOPV":0}
        # Composition table uses top-reference mapped reads per group (max), not sum across all refs.
        total_mapped = 0
        for ref_id, ref_len, mapped_reads in refs:
            if ref_len <= 0 or mapped_reads <= 0:
                continue
            grp = ref_to_group.get(ref_id, "Other")
            if isinstance(grp, str) and grp.startswith("nOPV"):
                grp = "nOPV"
            if grp not in counts:
                grp = "NonPolioEV"
            counts[grp] = max(counts[grp], mapped_reads)  # top ref per group
            total_mapped += mapped_reads

        if reads_filtered is None:
            reads_filtered = total_mapped + unmapped_reads
        elif unmapped_reads == 0 and reads_filtered > total_mapped:
            # Older benchmark outputs wrote idxstats from a BAM that had already
            # dropped unmapped reads, so reconstruct the filtered-read remainder.
            unmapped_reads = max(0, int(reads_filtered) - int(total_mapped))

        comp_rows.append({
            "Sample": sample_name,
            "Barcode": barcode_folder,
            "Unmapped": int(unmapped_reads),
            "NonPolioEV": int(counts["NonPolioEV"]),
            "WPV1": int(counts["WPV1"]),
            "WPV2": int(counts["WPV2"]),
            "WPV3": int(counts["WPV3"]),
            "Sabin1-related": int(counts["Sabin1-related"]),
            "Sabin2-related": int(counts["Sabin2-related"]),
            "Sabin3-related": int(counts["Sabin3-related"]),
            "nOPV": int(counts["nOPV"]),
        })

        best_group = ref_to_group.get(best_ref, "") if best_ref else ""
        if reads_filtered == 0 and total_mapped == 0 and unmapped_reads == 0:
            sample_class = "NoReads"
            ref_group = ""
        elif total_mapped == 0 and unmapped_reads > 0:
            sample_class = "Unmapped"
            ref_group = ""
        elif best_group:
            sample_class = best_group
            ref_group = best_group
        elif total_mapped > 0:
            sample_class = "Mapped"
            ref_group = ""
        else:
            sample_class = "NoReads"
            ref_group = ""

        cons_id = ""
        if best_ref and total_mapped > 0:
            cluster = ref_to_cluster.get(best_ref, best_ref)
            if cluster not in cluster_to_seqid:
                cluster_to_seqid[cluster] = f"SEQ{seq_counter:02d}"
                seq_counter += 1
            cons_id = cluster_to_seqid[cluster]


        # Sample summary information:
        # One row per sample per reference group, selecting the top reported reference in that group.
        # Ranking is mapped reads first, then breadth, then mean depth, then best hit rank.
        if not df_hits.empty:
            sub = df_hits[df_hits['sample_id'] == sample_id].copy()
            if len(sub):
                if 'ddns_group' not in sub.columns:
                    sub['ddns_group'] = sub['ref_id'].map(ref_to_group).fillna('Other')
                sub['mapped_reads_num'] = pd.to_numeric(sub.get('mapped_reads', 0), errors='coerce').fillna(0).astype(int)
                sub['breadth_num'] = pd.to_numeric(sub.get('breadth', 0), errors='coerce').fillna(0.0)
                sub['mean_depth_num'] = pd.to_numeric(sub.get('mean_depth', 0), errors='coerce').fillna(0.0)
                if 'hit_rank' in sub.columns:
                    sub['hit_rank_num'] = pd.to_numeric(sub['hit_rank'], errors='coerce').fillna(9999).astype(int)
                else:
                    sub['hit_rank_num'] = 9999

                sub['consensus_fasta'] = sub.get('consensus_fasta', '').fillna('').astype(str)
                sub['has_cons'] = sub['consensus_fasta'].str.len().gt(0)
                sub['consensus_len'] = pd.to_numeric(sub.get('consensus_len', 0), errors='coerce').fillna(0).astype(int)
                sub['n_count'] = pd.to_numeric(sub.get('n_count', 0), errors='coerce').fillna(0).astype(int)
                sub['n_frac'] = sub.apply(lambda r: (r['n_count']/r['consensus_len']) if r['consensus_len']>0 else 1.0, axis=1)
                sub['final_status_norm'] = sub['final_status'].fillna('').astype(str).str.strip().str.upper() if 'final_status' in sub.columns else ''
                sub['status_norm'] = sub['status'].fillna('').astype(str).str.strip().str.upper() if 'status' in sub.columns else ''
                sub['recovered'] = sub['final_status_norm'].eq('CONSENSUS_PASS') | sub['status_norm'].eq('PASS')

                # Keep the usual mapped-read threshold, but do not hide rows when a masked consensus FASTA exists.
                # This preserves FASTA links for outcomes like HIGH_N_FRAC or NO_CALLABLE_BASES.
                sub = sub[(sub['mapped_reads_num'] >= args.min_group_reads) | (sub['has_cons'])].copy()
                if len(sub):
                    sub = sub.sort_values(
                        ['ddns_group','mapped_reads_num','breadth_num','mean_depth_num','hit_rank_num','ref_id'],
                        ascending=[True, False, False, False, True, True]
                    )
                    best_rows = sub.groupby('ddns_group', as_index=False).head(1)
                    best_rows = best_rows.sort_values(
                        ['mapped_reads_num','breadth_num','mean_depth_num','hit_rank_num','ddns_group','ref_id'],
                        ascending=[False, False, False, True, True, True]
                    )
                    for _, r in best_rows.iterrows():
                        info_fasta = str(r.get('consensus_fasta','') or '')
                        if str(r.get('final_status_norm','')) == 'NO_CALLABLE_BASES' or str(r.get('status_norm','')) == 'NO_CALLABLE_BASES':
                            info_fasta = ''
                        info_rows.append({
                            'Sample': sample_name,
                            'Barcode': barcode_folder,
                            'Hit Rank': int(r.get('hit_rank_num', 9999)) if int(r.get('hit_rank_num', 9999)) != 9999 else '',
                            'Ref Id': str(r.get('ref_id','')),
                            'Reference Group': str(r.get('ddns_group','Other')),
                            'Mapped Reads': int(r.get('mapped_reads_num', 0)),
                            'Breadth': fmt_pct(r.get('breadth_num', 0.0)),
                            'Mean Depth': float(r.get('mean_depth_num', 0.0)),
                            'Consensus Recovered': 'Yes' if bool(r.get('recovered', False)) else 'No',
                            'Mutations': fmt_int(r.get('mutations','')),
                            'Number of Ns': fmt_int(r.get('n_count','')),
                            'Consensus FASTA': info_fasta
                        })

        if not df_hits.empty:
            sub = df_hits[(df_hits["sample_id"]==sample_id)].copy()
            if len(sub):
                sub["has_cons"] = sub["consensus_fasta"].fillna("").astype(str).str.len().gt(0).astype(int)
                sub["consensus_len"] = pd.to_numeric(sub["consensus_len"], errors="coerce").fillna(0).astype(int)
                sub["n_count"] = pd.to_numeric(sub["n_count"], errors="coerce").fillna(0).astype(int)
                sub["n_frac"] = sub.apply(lambda r: (r["n_count"]/r["consensus_len"]) if r["consensus_len"]>0 else 1.0, axis=1)
                sub["final_status_norm"] = sub["final_status"].fillna("").astype(str).str.strip().str.upper() if "final_status" in sub.columns else ""
                sub["status_norm"] = sub["status"].fillna("").astype(str).str.strip().str.upper() if "status" in sub.columns else ""
                sub["recovered"] = sub["final_status_norm"].eq("CONSENSUS_PASS") | sub["status_norm"].eq("PASS")
                sub = sub.sort_values(["recovered","breadth","mapped_reads","hit_rank"], ascending=[False,False,False,True])
                for _, hr in sub.iterrows():
                    top_rows.append({
                        "Sample": sample_name,
                        "Barcode": barcode_folder,
                        "Hit Rank": int(hr.get("hit_rank",0)) if pd.notna(hr.get("hit_rank",0)) else "",
                        "Ref Id": hr.get("ref_id",""),
                        "Reference Group": hr.get("ddns_group",""),
                        "Mapped Reads": int(hr.get("mapped_reads",0)) if pd.notna(hr.get("mapped_reads",0)) else 0,
                        "Breadth": fmt_pct(hr.get("breadth","")),
                        "Mean Depth": hr.get("mean_depth",""),
                        "Consensus Recovered": "Yes" if bool(hr.get("recovered",False)) else "No",
                        "Mutations": fmt_int(hr.get("mutations","")),
                        "Number of Ns": fmt_int(hr.get("n_count","")),
                        "Consensus FASTA": hr.get("consensus_fasta","")
                    })
            else:
                top_rows.append({
                    "Sample": sample_name, "Barcode": barcode_folder, "Hit Rank":"",
                    "Ref Id":"", "Reference Group":"", "Mapped Reads":0, "Breadth":"",
                    "Mean Depth":"", "Consensus Recovered":"No", "Mutations":"", "Number of Ns":"", "Consensus FASTA":""
                })
        else:
            top_rows.append({
                "Sample": sample_name, "Barcode": barcode_folder, "Hit Rank":"",
                "Ref Id":"", "Reference Group":"", "Mapped Reads":0, "Breadth":"",
                "Mean Depth":"", "Consensus Recovered":"No", "Mutations":"", "Number of Ns":"", "Consensus FASTA":""
            })

    pd.DataFrame(comp_rows)[COMP_COLUMNS].to_csv(out / "Table_Composition of samples.csv", index=False)
    df_info = pd.DataFrame(info_rows)
    if df_info.empty:
        df_info = pd.DataFrame(columns=INFO_COLUMNS)
    else:
        for c in INFO_COLUMNS:
            if c not in df_info.columns:
                df_info[c] = ""
        df_info = df_info[INFO_COLUMNS]
    df_info.to_csv(out / "Table_Sample summary information.csv", index=False)

    # Extra summary tables for signal vs confirmed reporting
    thr_signal = int(getattr(args, "signal_threshold", 50))
    def _is_yes(v):
        s = str(v).strip().lower()
        return s in ("yes","y","true","1","pass","passed")

    # Ensure empty outputs keep headers so report_html can read them safely
    empty_cols = list(df_info.columns) if (df_info is not None and hasattr(df_info, "columns") and len(df_info.columns)) else INFO_COLUMNS
    if df_info is None or df_info.empty:
        pd.DataFrame(columns=empty_cols).to_csv(out / "Table_WPV1_signal.csv", index=False)
        pd.DataFrame(columns=empty_cols).to_csv(out / "Table_WPV1_confirmed.csv", index=False)
        pd.DataFrame(columns=empty_cols).to_csv(out / "Table_Confirmed_All.csv", index=False)
        pd.DataFrame(columns=empty_cols).to_csv(out / "Table_Confirmed_LowReads.csv", index=False)
    else:
        df_tmp = df_info.copy()
        if "Mapped Reads" in df_tmp.columns:
            df_tmp["Mapped Reads"] = pd.to_numeric(df_tmp["Mapped Reads"], errors="coerce").fillna(0).astype(int)
        else:
            df_tmp["Mapped Reads"] = 0
        # Confirmed All (any reference group): Consensus Recovered == YES
        if "Consensus Recovered" in df_tmp.columns:
            df_conf_all = df_tmp[df_tmp["Consensus Recovered"].apply(_is_yes)].copy()
        else:
            df_conf_all = df_tmp.iloc[0:0].copy()
        df_conf_all.to_csv(out / "Table_Confirmed_All.csv", index=False)
        # WPV1-only subset
        if "Reference Group" in df_tmp.columns:
            df_wpv1 = df_tmp[df_tmp["Reference Group"].astype(str) == "WPV1"].copy()
        else:
            df_wpv1 = df_tmp.iloc[0:0].copy()
        # WPV1 Signal: mapped reads >= threshold
        df_sig = df_wpv1[df_wpv1["Mapped Reads"] >= thr_signal].copy() if len(df_wpv1) else df_wpv1
        df_sig.to_csv(out / "Table_WPV1_signal.csv", index=False)
        # WPV1 Confirmed: Consensus Recovered == YES (independent of threshold)
        if "Consensus Recovered" in df_wpv1.columns:
            df_conf = df_wpv1[df_wpv1["Consensus Recovered"].apply(_is_yes)].copy()
        else:
            df_conf = df_wpv1.iloc[0:0].copy()
        df_conf.to_csv(out / "Table_WPV1_confirmed.csv", index=False)

        # Additional table: consensus recovered with mapped reads below the main summary threshold
        # This captures recovered consensuses that are excluded from Table_Sample summary information because of --min-group-reads filtering.
        low_thr = int(getattr(args, "min_group_reads", 50))
        # Build from df_hits (not df_info) so we do not lose low-read recovered consensuses
        df_low = pd.DataFrame(columns=empty_cols)
        try:
            if not df_hits.empty:
                sub2 = df_hits.copy()
                # Ensure group column exists
                if "ddns_group" not in sub2.columns:
                    sub2["ddns_group"] = sub2["ref_id"].map(ref_to_group).fillna("Other")
                # Numeric columns
                sub2["mapped_reads_num"] = pd.to_numeric(sub2.get("mapped_reads", 0), errors="coerce").fillna(0).astype(int)
                sub2["breadth_num"] = pd.to_numeric(sub2.get("breadth", 0), errors="coerce").fillna(0.0)
                sub2["mean_depth_num"] = pd.to_numeric(sub2.get("mean_depth", 0), errors="coerce").fillna(0.0)
                sub2["hit_rank_num"] = pd.to_numeric(sub2.get("hit_rank", 9999), errors="coerce").fillna(9999).astype(int)
                sub2["consensus_fasta"] = sub2.get("consensus_fasta", "").fillna("").astype(str)
                sub2["has_cons"] = sub2["consensus_fasta"].str.len().gt(0)
                sub2["consensus_len"] = pd.to_numeric(sub2.get("consensus_len", 0), errors="coerce").fillna(0).astype(int)
                sub2["n_count"] = pd.to_numeric(sub2.get("n_count", 0), errors="coerce").fillna(0).astype(int)
                sub2["n_frac"] = sub2.apply(lambda r: (r["n_count"]/r["consensus_len"]) if r["consensus_len"]>0 else 1.0, axis=1)
                sub2["recovered"] = (sub2["has_cons"]) & (sub2["consensus_len"] >= args.min_cons_len) & (sub2["n_frac"] <= args.max_n_frac)

                # Keep recovered only, mapped reads < low_thr
                sub2 = sub2[(sub2["recovered"] == True) & (sub2["mapped_reads_num"] < low_thr)].copy()

                if len(sub2):
                    # Exclude rows already present in Confirmed_All (which is derived from df_info, typically >= low_thr)
                    try:
                        existing_keys = set(
                            tuple(x)
                            for x in df_conf_all[["Sample","Barcode","Reference Group"]].fillna("").astype(str).values.tolist()
                        )
                    except Exception:
                        existing_keys = set()

                    # Convert df_hits rows into output schema; pick best ref per sample per group
                    sub2 = sub2.sort_values(
                        ["sample_id","ddns_group","mapped_reads_num","breadth_num","mean_depth_num","hit_rank_num"],
                        ascending=[True, True, False, False, False, True]
                    )
                    best2 = sub2.groupby(["sample_id","ddns_group"], as_index=False).head(1)

                    low_rows = []
                    for _, r2 in best2.iterrows():
                        # sample_id was created as safe_sample__barcodeXX; we have sample_name and barcode_folder via parsing
                        sid = str(r2.get("sample_id",""))
                        # Find matching barcode and sample name from earlier loop inputs if possible
                        # Safer: use the barcode map and sample_id suffix
                        b = ""
                        m = re.search(r"__(barcode\d+)$", sid)
                        if m:
                            b = m.group(1)
                        # best-effort sample name: from barcode map if possible, else prefix of sample_id
                        sname = ""
                        if b and b in barcode_map:
                            sname = barcode_map[b]
                        else:
                            sname = sid.split("__")[0] if "__" in sid else sid

                        grp = str(r2.get("ddns_group","Other"))
                        key = (str(sname), str(b), grp)
                        if key in existing_keys:
                            continue

                        low_rows.append({
                            "Sample": sname,
                            "Barcode": b,
                            "Hit Rank": int(r2.get("hit_rank_num", 9999)) if int(r2.get("hit_rank_num", 9999)) != 9999 else "",
                            "Ref Id": str(r2.get("ref_id","")),
                            "Reference Group": grp,
                            "Mapped Reads": int(r2.get("mapped_reads_num", 0)),
                            "Breadth": fmt_pct(r2.get("breadth_num", 0.0)),
                            "Mean Depth": float(r2.get("mean_depth_num", 0.0)),
                            "Consensus Recovered": "Yes",
                            "Mutations": fmt_int(r2.get("mutations","")),
                            "Number of Ns": fmt_int(r2.get("n_count","")),
                            "Consensus FASTA": str(r2.get("consensus_fasta","") or "")
                        })

                    df_low = pd.DataFrame(low_rows)
                    if df_low.empty:
                        df_low = pd.DataFrame(columns=empty_cols)
                    else:
                        for c in empty_cols:
                            if c not in df_low.columns:
                                df_low[c] = ""
                        df_low = df_low[empty_cols]
        except Exception:
            df_low = pd.DataFrame(columns=empty_cols)

        df_low.to_csv(out / "Table_Confirmed_LowReads.csv", index=False)

    # Summary counts (unique samples) by reference group, based on Sample summary information table
    df_info = pd.DataFrame(info_rows)
    groups = ["WPV1","WPV2","WPV3","Sabin1-related","Sabin2-related","Sabin3-related","nOPV","NonPolioEV"]
    rows = []
    if not df_info.empty and "Reference Group" in df_info.columns:
        rgcol = df_info["Reference Group"].fillna("").astype(str)
        for g in groups:
            if g == "nOPV":
                df_g = df_info[rgcol.str.startswith("nOPV")]
            else:
                df_g = df_info[df_info["Reference Group"] == g]
            rows.append({"Group": g, "Unique_Samples": int(df_g[["Sample","Barcode"]].drop_duplicates().shape[0])})
    else:
        for g in groups:
            rows.append({"Group": g, "Unique_Samples": 0})
    pd.DataFrame(rows).to_csv(out / "Summary_Unique_Samples_By_Group.csv", index=False)

    # Confirmed summary cards: count unique sample+barcode detections per reference group,
    # allowing one barcode to contribute to multiple groups but only once per group.
    confirmed_groups = ["WPV1","WPV2","WPV3","Sabin1-related","Sabin2-related","Sabin3-related","nOPV","NonPolioEV"]
    confirmed_rows = []
    run_overview_rows = []
    if df_hits.empty:
        for g in confirmed_groups:
            confirmed_rows.append({"Group": g, "Confirmed_Barcodes": 0})
        total_samples = len(all_barcodes)
        run_overview_rows = [
            {"Metric": "Total Samples", "Value": total_samples},
            {"Metric": "Any Confirmed", "Value": 0},
            {"Metric": "Mixed Confirmed", "Value": 0},
            {"Metric": "No Confirmed Consensus", "Value": total_samples},
        ]
    else:
        df_confsum = df_hits.copy()
        if "ddns_group" not in df_confsum.columns:
            df_confsum["ddns_group"] = df_confsum["ref_id"].map(ref_to_group).fillna("Other")
        df_confsum["final_status_norm"] = df_confsum["final_status"].fillna("").astype(str).str.strip().str.upper() if "final_status" in df_confsum.columns else ""
        df_confsum["status_norm"] = df_confsum["status"].fillna("").astype(str).str.strip().str.upper() if "status" in df_confsum.columns else ""
        df_confsum["recovered"] = df_confsum["final_status_norm"].eq("CONSENSUS_PASS") | df_confsum["status_norm"].eq("PASS")
        df_confsum = df_confsum[df_confsum["recovered"] == True].copy()

        def _norm_group(v):
            s = str(v or "").strip()
            if s.startswith("nOPV"):
                return "nOPV"
            return s

        if not df_confsum.empty:
            df_confsum["summary_group"] = df_confsum["ddns_group"].apply(_norm_group)
            df_confuniq = df_confsum[["sample_id","summary_group"]].drop_duplicates()
        else:
            df_confuniq = pd.DataFrame(columns=["sample_id","summary_group"])

        for g in confirmed_groups:
            n = int(df_confuniq[df_confuniq["summary_group"] == g]["sample_id"].nunique()) if not df_confuniq.empty else 0
            confirmed_rows.append({"Group": g, "Confirmed_Barcodes": n})

        total_samples = len(all_barcodes)
        any_confirmed = int(df_confuniq["sample_id"].nunique()) if not df_confuniq.empty else 0
        if not df_confuniq.empty:
            mixed_confirmed = int((df_confuniq.groupby("sample_id")["summary_group"].nunique() > 1).sum())
        else:
            mixed_confirmed = 0
        no_confirmed = max(total_samples - any_confirmed, 0)
        run_overview_rows = [
            {"Metric": "Total Samples", "Value": total_samples},
            {"Metric": "Any Confirmed", "Value": any_confirmed},
            {"Metric": "Mixed Confirmed", "Value": mixed_confirmed},
            {"Metric": "No Confirmed Consensus", "Value": no_confirmed},
        ]

    pd.DataFrame(confirmed_rows).to_csv(out / "Summary_Confirmed_By_Group.csv", index=False)
    pd.DataFrame(run_overview_rows).to_csv(out / "Summary_Run_Overview.csv", index=False)

    df_top_out = pd.DataFrame(top_rows)
    if not df_top_out.empty and "Mutations" in df_top_out.columns:
        df_top_out["Mutations"] = df_top_out["Mutations"].apply(fmt_int)
    df_top_out.to_csv(out / "Table_Top hits.csv", index=False)
    # Final group counts using recovered hits
    groups = ["WPV1","Sabin1-related","Sabin2-related","Sabin3-related"]
    counts_out = []
    if df_hits.empty:
        for g in groups:
            counts_out.append({"Group": g, "Recovered_Consensus": 0, "Unique_Samples": 0})
    else:
        df_rec = df_hits.copy()
        df_rec["consensus_fasta"] = df_rec["consensus_fasta"].fillna("").astype(str)
        df_rec["has_cons"] = df_rec["consensus_fasta"].str.len().gt(0)
        df_rec["consensus_len"] = pd.to_numeric(df_rec["consensus_len"], errors="coerce").fillna(0).astype(int)
        df_rec["n_count"] = pd.to_numeric(df_rec["n_count"], errors="coerce").fillna(0).astype(int)
        df_rec["n_frac"] = df_rec.apply(lambda r: (r["n_count"]/r["consensus_len"]) if r["consensus_len"]>0 else 1.0, axis=1)
        df_rec["final_status_norm"] = df_rec["final_status"].fillna("").astype(str).str.strip().str.upper() if "final_status" in df_rec.columns else ""
        df_rec["status_norm"] = df_rec["status"].fillna("").astype(str).str.strip().str.upper() if "status" in df_rec.columns else ""
        df_rec["recovered"] = df_rec["final_status_norm"].eq("CONSENSUS_PASS") | df_rec["status_norm"].eq("PASS")
        df_rec = df_rec[df_rec["recovered"] == True].copy()

        for g in groups:
            if "ddns_group" in df_rec.columns:
                df_g = df_rec[df_rec["ddns_group"] == g]
            else:
                df_g = df_rec.iloc[0:0]
            recovered_consensus = int(len(df_g))
            unique_samples = int(df_g["sample_id"].nunique()) if recovered_consensus else 0
            counts_out.append({"Group": g, "Recovered_Consensus": recovered_consensus, "Unique_Samples": unique_samples})

    pd.DataFrame(counts_out).to_csv(out / "Final_Group_Counts.csv", index=False)

if __name__ == "__main__":
    main()
