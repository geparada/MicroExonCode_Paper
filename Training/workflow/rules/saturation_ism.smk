import pandas as pd

known_mics_hg38 = pd.read_table("resources/inputs/wt/MicEvents.hg38.tab.gz")
known_mics_hg38 = known_mics_hg38.loc[known_mics_hg38.label == "hi"]["event"]

rule saturation_ism_make_variant_files_hg38:
  input:
    "resources/inputs/wt/MicEvents.hg38.tab.gz",
    "resources/genomes/hg38.2bit"
  output:
    [f"resources/inputs/saturation_mutagenesis/{event}_variants.hg38.tab.gz" for event in known_mics_hg38]
  params:
    window=300,
    label="hi"
  conda:
    "../envs/base_env.yml"
  script:
    "../scripts/make_saturation_ism_variant_files.py"

known_mics_mm10 = pd.read_table("resources/inputs/wt/MicEvents.mm10.tab.gz")
known_mics_mm10 = known_mics_mm10.loc[known_mics_mm10.label == "hi"]["event"]

rule saturation_ism_make_variant_files_mm10:
  input:
    "resources/inputs/wt/MicEvents.mm10.tab.gz",
    "resources/genomes/mm10.2bit"
  output:
    [f"resources/inputs/saturation_mutagenesis/{event}_variants.mm10.tab.gz" for event in known_mics_mm10]
  params:
    window=300,
    label="hi"
  conda:
    "../envs/base_env.yml"
  script:
    "../scripts/make_saturation_ism_variant_files.py"

rule predict_saturation_in_silico_mutagenesis:
  input:
    expand(
      "results/predictions/{model}/saturation_mutagenesis/predictions.{event}_variants.hg38.csv.gz",
      event=known_mics_hg38,
      model=["known_microexons", "novel_microexons", "original_model"],
    ),
    expand("results/predictions/{model}/wt/predictions.MicEvents.hg38.csv.gz",
      model=["known_microexons", "novel_microexons", "original_model"],
    ),
    expand(
      "results/predictions/{model}/saturation_mutagenesis/predictions.{event}_variants.mm10.csv.gz",
      event=known_mics_mm10,
      model=["known_microexons", "novel_microexons", "original_model"],
    ),
    expand("results/predictions/{model}/wt/predictions.MicEvents.mm10.csv.gz",
      model=["known_microexons", "novel_microexons", "original_model"],
    )

putative_detected_hg38, = glob_wildcards(
  "resources/inputs/putative_ISM_input/{event}_variants.hg38.tab.gz"
)
rule predict_saturation_in_silico_mutagenesis_putative_detected:
  input:
    expand(
      "results/predictions/{model}/putative_ISM_input/predictions.{event}_variants.hg38.csv.gz",
      event=putative_detected_hg38,
      model=["known_microexons", "novel_microexons", "original_model"]
    )

control_variants, = glob_wildcards(
  "resources/inputs/Control_variants/{event}.tab.gz"
)
rule predict_control_variants:
  input:
    expand(
      "results/predictions/{model}/Control_variants/predictions.{event}.csv.gz",
      event=control_variants,
      model=["known_microexons", "novel_microexons", "original_model"]
    )

mssng_variants, = glob_wildcards(
  "resources/inputs/MSSNG/{event}.tab.gz"
)
rule predict_mssng_variants:
  input:
    expand(
      "results/predictions/{model}/MSSNG/predictions.{event}.csv.gz",
      event=mssng_variants,
      model=["known_microexons", "novel_microexons", "original_model"]
    )



# rule predict_saturation_in_silico_mutagenesis:
#   input:
#     model="results/models/{model}/xgboost_model.json",
#     exons="resources/saturation_mutagenesis/MicEvents.variants.{species}.neural.tab.gz",
#     ss_strength="resources/saturation_mutagenesis/MicEvents.variants.{species}.neural.features.ss_strength.csv.gz",
#     gc_content="resources/saturation_mutagenesis/MicEvents.variants.{species}.neural.features.gc_content.csv.gz",
#     bp_scores="resources/saturation_mutagenesis/MicEvents.variants.{species}.neural.features.svm_bp_scores.csv.gz",
#     manual_rbp="resources/saturation_mutagenesis/MicEvents.variants.{species}.neural.features.manual_rbp_binding.csv.gz",
#     cis_bp="resources/saturation_mutagenesis/MicEvents.variants.{species}.neural.features.aggregated_matches.csv.gz",
#     cossmo="resources/saturation_mutagenesis/MicEvents.variants.{species}.neural.features.cossmo_features.csv.gz",
#     seqweaver="resources/saturation_mutagenesis/MicEvents.variants.{species}.neural.features.seqweaver.csv.gz",
#     nsr100="resources/saturation_mutagenesis/MicEvents.variants.{species}.neural.features.nsr100_features.csv.gz",
#     kmers="resources/saturation_mutagenesis/MicEvents.variants.{species}.neural.features.kmer_features.csv.gz"
#   output:
#     "results/saturation_ism/{model}/ism_predictions.MicEvents.variants.{species}.csv.gz"
#   conda:
#     "../envs/xgboost.yml"
#   wildcard_constraints:
#     species="\w{2,6}\d{1,2}"
#   script:
#     "../scripts/predict_microexons.py"

# rule variant_compute_dpsi:
#   input:
#     wt_pred_file="results/predictions/{model}/wt/predictions.MicEvents.{species}.csv.gz",
#     mt_pred_file=rules.predict_saturation_in_silico_mutagenesis.output[0]
#   output: "results/saturation_ism/{model}/ism_predictions.MicEvents.variants.{species}.dpsi.csv.gz"
#   script: "../scripts/saturation_mutagenesis_compute_dpsi.py"
#
# rule compute_ism_heatmap:
#   input:
#     events="resources/inputs/wt/MicEvents.hg38.tab.gz",
#     genome="resources/genomes/hg38.2bit",
#     ism_score_diff=rules.variant_compute_dpsi.output[0],
#     hexamer_pvals="resources/saturation_mutagenesis/hexamer_p_vals.pkl"
#   output:
#     motif_pvals="results/saturation_ism/{model}/motif_pvals.csv",
#     neural_high_pred_diff_heatmap="results/saturation_ism/{model}/neural_high_pred_diff_heatmap.csv",
#     neural_high_pred_diff_heatmap_by_pval="results/saturation_ism/{model}/neural_high_pred_diff_heatmap_by_pval.csv"
#   conda: "../envs/base_env.yml"
#   script: "../scripts/compute_ism_heatmap.py"
#
# rule plot_ism_heatmap:
#   input: "results/saturation_ism/{model}/neural_high_pred_diff_heatmap.csv"
#   output: "results/saturation_ism/{model}/ism_heatmap.pdf"
#   conda: "../envs/r.yml"
#   script: "../scripts/plot_ism_heatmap.R"
