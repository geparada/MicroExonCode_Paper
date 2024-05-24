from math import ceil

import yaml
import pathlib

N_TCAG_SUBTASKS = 500

feature_keys = yaml.load(
  open(pathlib.Path(config["model_path"]) / "feature_set.yaml"),
  Loader=yaml.FullLoader
)['feature_set']

def divup(a, b):
  return int(ceil(float(a) / b))


# rule predictions:
#   input:
#     model=pathlib.Path(config["model_path"]) / "xgboost_model.json",
#     exons=config["input_file_schema"],
#     length=rules.features_length.output[0],
#     ss_strength=rules.features_ss_strength.output[0],
#     gc_content=rules.features_gc_content.output[0],
#     bp_scores=rules.features_svm_bp.output[0],
#     manual_rbp=rules.features_manual_rbp_motifs.output[0],
#     cis_bp=rules.features_cisbp.output[0],
#     cossmo=rules.features_cossmo.output[0],
#     seqweaver=rules.features_seqweaver.output[0],
#     nsr100=rules.features_nsr100.output[0],
#     kmers=rules.features_kmers.output[0]
#   output:
#     "results/{task}/predictions.{prefix}.{species}.csv.gz"
#   params:
#     feature_keys=feature_keys
#   conda:
#     "../envs/predictions.yml"
#   wildcard_constraints:
#     species="\w{2,6}\d{1,2}"
#   script:
#     "../scripts/predict_microexons.py"

rule predictions_for_model:
  input:
    model="results/models/{model}/xgboost_model.json",
    exons=config["input_file_schema"],
    length=rules.features_length.output[0],
    ss_strength=rules.features_ss_strength.output[0],
    gc_content=rules.features_gc_content.output[0],
    bp_scores=rules.features_svm_bp.output[0],
    manual_rbp=rules.features_manual_rbp_motifs.output[0],
    cis_bp=rules.features_cisbp.output[0],
    cossmo=rules.features_cossmo.output[0],
    seqweaver=rules.features_seqweaver.output[0],
    nsr100=rules.features_nsr100.output[0],
    kmers=rules.features_kmers.output[0],
    cisbp_features="resources/cisbp_feature_names.txt"
  output:
    "results/predictions/{model}/{task}/predictions.{prefix}.{species}.csv.gz"
  params:
    feature_keys=feature_keys
  conda:
    "../envs/predictions.yml"
  wildcard_constraints:
    # species="\w{2,6}\d{1,2}",
    model="known_microexons|novel_microexons|original_model|known_microexons_hg38"
  # priority: 1000
  script:
    "../scripts/predict_microexons.py"


rule compute_shapley_values:
  input:
    model = pathlib.Path(config["model_path"]) / "xgboost_model.pkl",
    exons = config["input_file_schema"],
    length = rules.features_length.output[0],
    ss_strength = rules.features_ss_strength.output[0],
    gc_content = rules.features_gc_content.output[0],
    bp_scores = rules.features_svm_bp.output[0],
    manual_rbp = rules.features_manual_rbp_motifs.output[0],
    cis_bp = rules.features_cisbp.output[0],
    cossmo = rules.features_cossmo.output[0],
    seqweaver = rules.features_seqweaver.output[0],
    nsr100 = rules.features_nsr100.output[0],
    kmers = rules.features_kmers.output[0]
  output:
    "results/{task}/shapley_values.{prefix}.{species}.csv.gz"
  conda:
    "../envs/predictions.yml"
  wildcard_constraints:
    species="\w{2,6}\d{1,2}"
  script:
    "../scripts/compute_shapley_values.py"


# def get_tcag_full_genotype_inputs(wildcards):
#   files = []
#   for file in pathlib.Path("resources/inputs/tcag_full_genotypes").glob("*.tab.gz"):
#     task_id = file.name[:file.name.index(".tab.gz")]
#     files.append(f"results/predictions/known_microexons/tcag_full_genotypes/predictions.{task_id}.csv.gz")
#
#   return files
#
#
# def get_tcag_full_genotype_inputs_for_subtask(wildcards):
#   subtask_id = int(wildcards["subtask_id"])
#   input_files = list(pathlib.Path("resources/inputs/tcag_full_genotypes").glob("*.tab.gz"))
#   n_files_per_task = divup(len(input_files), N_TCAG_SUBTASKS)
#
#   files = []
#   for i, file in enumerate(input_files):
#     if (i // n_files_per_task) == subtask_id:
#       files.append(str(file))
#
#   return files
#
#
# rule make_tcag_subtasks:
#   input: get_tcag_full_genotype_inputs_for_subtask
#   output: "resources/inputs/tcag_full_genotypes_aggregated/tcag_full_genotypes_part_{subtask_id}.hg38.tab.gz"
#   wildcard_constraints:
#     subtask_id="\d{3}"
#   run:
#     import pandas as pd
#     df = []
#     for file in input:
#       df.append(pd.read_table(file, index_col="event"))
#
#     df = pd.concat(df)
#     df.to_csv(output[0], sep="\t", index_label="event")


# rule predict_tcag_full_genotypes:
#   # input: get_tcag_full_genotype_inputs
#   # input: expand(rules.make_tcag_subtasks.output[0], subtask_id=[f"{i:03d}" for i in range(500)])
#   input:
#     expand("results/predictions/{model}/tcag_full_genotypes_aggregated/"
#       "predictions.tcag_full_genotypes_part_{subtask_id}.hg38.csv.gz",
#       model=["known_microexons", "novel_microexons"],
#       subtask_id=[f"{i:03d}" for i in range(500)])
