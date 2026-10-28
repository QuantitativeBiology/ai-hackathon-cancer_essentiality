import json, sys
from pathlib import Path
ROOT=Path(r"C:\Users\ritab\Desktop\project")
f = sys.argv[1] if len(sys.argv)>1 else "repro_out.json"
me=json.load(open(f))
crit=json.loads((ROOT/"configs/transfer_criteria.json").read_text(encoding="utf-8"))
print("selftest", me["selftest_ap_vs_sklearn"])
rows=[]
def add(name, mine, proj):
    d = None if (mine is None or proj is None) else mine-proj
    rows.append((name, mine, proj, d))
for bg in ["RPE1","hPSC"]:
    P=json.load(open(ROOT/f"results/transfer/setB/{bg}.json"))
    S=json.load(open(ROOT/f"results/transfer/setB/{bg}_secondary.json")) if False else None
    m=me["setB"][bg]
    add(f"{bg} plan n_genes", m["plan_n_genes"], P["plan"]["n_genes"]); add(f"{bg} plan n_blocks", m["plan_n_blocks"], P["plan"]["n_blocks"])
    for e,(ss,mn,mc) in m["consensus_qc"].items():
        fq=crit["backgrounds"][bg]["harmonisation_validation"]["consensus_qc"][e]
        add(f"{bg} {e} consensus SSMD", ss, fq["ssmd"]); add(f"{bg} {e} median NEG eval", mn, fq["median_neg_eval"]); add(f"{bg} {e} median CEG eval", mc, fq["median_ceg_eval"])
    for key in ["Q3a_primary","Q3a_effect_label","Q3b_primary","Q3b_vs_B0","Q3b_rare_vs_B1n"]:
        st=P["statistics"][key]
        add(f"{bg} {key} T point", m[key]["point"], st["T"]["point"])
        add(f"{bg} {key} T_ensemble point", m[key]["ensemble_point"], st["T_ensemble"]["point"])
        for e,v in m[key]["per_entry"].items(): add(f"{bg} {key} entry {e}", v, st["T"]["per_entry"][e])
        for e,v in m[key]["absT"].items(): add(f"{bg} {key} absT {e}", v, st["T"]["model_abs_per_entry"][e])
        for e,v in m[key]["absComp"].items(): add(f"{bg} {key} absComp {e}", v, st["T"]["comparator_abs_per_entry"][e])
        rc=crit["backgrounds"][bg]["replicate_ceiling"][key]
        add(f"{bg} {key} ceiling point", m[key]["ceiling_point"], rc["ceiling_skill"])
        for e,v in m[key]["ceiling_per_entry"].items(): add(f"{bg} {key} ceiling entry {e}", v, rc["per_entry"][e])
        if "boot" in m[key]:
            add(f"{bg} {key} ci_lo", m[key]["boot"]["ci_lo"], st["T"]["ci_lo"]); add(f"{bg} {key} ci_hi", m[key]["boot"]["ci_hi"], st["T"]["ci_hi"]); add(f"{bg} {key} se", m[key]["boot"]["se"], st["T"]["se"])
        if "swap" in m[key]:
            for pid,v in m[key]["swap"].items():
                k=[s for s in st["swap"] if s.endswith(pid)][0]
                add(f"{bg} {key} swap {pid} point", v["point"], st["swap"][k]["point"])
                if "ci_lo" in v: add(f"{bg} {key} swap {pid} ci_lo", v["ci_lo"], st["swap"][k]["ci_lo"]); add(f"{bg} {key} swap {pid} ci_hi", v["ci_hi"], st["swap"][k]["ci_hi"])
    if bg=="RPE1":
        for key in ["Q3a_primary","Q3b_primary"]:
            pass
SB=json.load(open(ROOT/"results/transfer/setB/summary.json"))
for g,v in me["H1"]["per_gene"].items():
    add(f"H1 model delta {g}", v["model"], SB["hypotheses"]["H1"]["model_delta"][g]); add(f"H1 obs delta {g}", v["obs"], SB["hypotheses"]["H1"]["observed_delta"][g])
for k,pk in [("pooled",None),("brown","brown2019"),("bowden","bowden2020")]:
    h = SB["hypotheses"]["H2"] if pk is None else SB["hypotheses"]["H2"]["per_study"][pk]
    add(f"H2 {k} n", me["H2"][k]["n"], h["n_genes"]); add(f"H2 {k} spearman", me["H2"][k]["spearman"], h["spearman"]); add(f"H2 {k} sign", me["H2"][k]["sign"], h["sign_agreement"])
H3=SB["hypotheses"]["H3"]
add("H3 n_obs_dep", me["H3"]["n_obs_dep"], H3["n_observed_depleted"]); add("H3 n_pred_T", me["H3"]["n_pred_T"], H3["n_predicted_T"]); add("H3 n_pred_noEC", me["H3"]["n_pred_noEC"], H3["n_predicted_noEC"])
for g,v in me["H3"]["T_pred"].items(): add(f"H3 T pred {g}", v, H3["per_gene"][g]["T"]); add(f"H3 noEC pred {g}", me["H3"]["noEC_pred"][g], H3["per_gene"][g]["noEC"])
A=json.load(open(ROOT/"results/transfer/setA/summary.json"))
for key in ["Q3a_primary","Q3b_primary"]:
    for ln,v in me["setA"][key]["per_line"].items():
        p=A["statistics"][key]["per_line"][ln]
        add(f"setA {key} {ln} point", v["point"], p["point"]); add(f"setA {key} {ln} n_U", v["n_U"], p["n_universe"])
        add(f"setA {key} {ln} pct", v["pct"], A["statistics"][key]["descriptive"]["per_line"][ln]["test_percentile"])
        if "ci_lo" in v: add(f"setA {key} {ln} ci_lo", v["ci_lo"], p["ci_lo"]); add(f"setA {key} {ln} ci_hi", v["ci_hi"], p["ci_hi"])
    add(f"setA {key} median pct", me["setA"][key]["median_pct"], A["statistics"][key]["descriptive"]["median_test_percentile"])
mx=0
for r in rows:
    flag = "" if r[3] is None or abs(r[3])<=1e-4 else "  <-- DIFF"
    if r[3] is not None: mx=max(mx,abs(r[3]))
    print(f"{r[0]:70s} {float(r[1]):>14.6g} {float(r[2]) if r[2] is not None else float('nan'):>14.6g} {(r[3] if r[3] is not None else float('nan')):>11.2e}{flag}")
print("n compared", len(rows), "max abs diff", mx)
print(json.dumps({k:me[k] for k in ["comparator_recompute_maxabs","n_dev_rows"]},indent=0))
print(json.dumps({bg:me["setB"][bg]["diag_Q3a"] for bg in ["RPE1","hPSC"]},indent=0))
print(json.dumps({bg:me["setB"][bg]["final_eq_affine"] for bg in ["RPE1","hPSC"]}))
print(json.dumps(me["H3"]["POU5F1_in_mair_arms"]), me["H3"]["dep"], me["H3"]["result"], me["H1"]["result"], me["H2"]["pooled"]["result"])
for bg in ["RPE1","hPSC"]:
    for key in ["Q3a_primary","Q3b_primary","Q3b_rare_vs_B1n"]:
        print(bg,key,"n_pos",me["setB"][bg][key]["n_pos"],"n_U",me["setB"][bg][key]["n_U"])
