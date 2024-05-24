from snakemake.remote.SFTP import RemoteProvider
import pathlib
# SFTP = RemoteProvider(username=config_predict_microexons["tcag_credentials"]["username"], password=config_predict_microexons["tcag_credentials"]["password"])

SUBSET, SAMPLE_ID = glob_wildcards(
  pathlib.Path(config["tcag_input_path"]) / "microexon_overlap.{subset}/{sample_id,[^/]+}.tsv.gz"
)
SUBSET, SAMPLE_ID = zip(*[(s, i) for s, i in zip(SUBSET, SAMPLE_ID) if "/" not in s])

variant_subsets = ["MSSNG_CG", "MSSNG_ILMN", "SSC", "SPARK_WGS_1", "SPARK_WGS_2", "SPARK_WGS_3"]
# variant_subsets = ["SPARK_WGS_1", "SPARK_WGS_2", "SPARK_WGS_3"]


subworkflow microexon_prediction:
  workdir: "../microexon-code-workflow/workflow"

def get_metadata_path(wildcards):
  if wildcards.subset.startswith("MSSNG"):
    return "resources/tcag/metadata/MSSNG_metadata.tsv"
  else:
    return f"resources/tcag/metadata/{wildcards.subset}_metadata.tsv"


def variant_files_for_subset(wildcards):
  return ["resources/tcag/filtered_variants/{subset}/{sample_id}.tsv.gz".format(subset=subset, sample_id=sample_id)
          for subset, sample_id in zip(SUBSET, SAMPLE_ID) if subset == wildcards.subset]


rule filter_tcag_file:
  """This rule filters on only those variants that have the 'high_quality' columns set to 'true'."""
  input:
    pathlib.Path(config["tcag_input_path"]) / "microexon_overlap.{subset}/{sample_id}.tsv.gz"
  output:
    "resources/tcag/filtered_variants/{subset}/{sample_id}.tsv.gz"
  wildcard_constraints:
    subset="[A-Z_123]+"
  params:
    filter_field_name="high_quality",
    filter_value="true"
  script:
    "../scripts/filter_tcag_file.py"


rule filter_all_tcag_files:
  input:
    expand(rules.filter_tcag_file.output[0], zip, subset=SUBSET, sample_id=SAMPLE_ID)


rule fix_tcag_files:
  """Remove dynamically named column that causes problems for hail"""
  input:
    pathlib.Path(config["tcag_input_path"]) / "microexon_overlap.{subset}/{sample_id}.tsv.gz"
  output:
    pathlib.Path(config["tcag_input_path"]) / "microexon_overlap.{subset}/fixed/{sample_id}.tsv.gz"
  wildcard_constraints:
    sample_id="[^/]+"
  shell:
    "gzip -dc {input} | cut -f11 --complement | bgzip > {output}"

rule fix_all_tcag_files:
  input: expand(rules.fix_tcag_files.output[0], zip, subset=SUBSET, sample_id=SAMPLE_ID)


rule collect_variant_statistics:
  input:
    variants=variant_files_for_subset,
    regions=[
        "resources/MSSNG.variants.2020-07-24/mic_events.neural-high.300nt_flanking.hg38.bed.gz",
        "resources/MSSNG.variants.2020-07-24/mic_events.neural-high.flanking_exons.300nt_flanking.hg38.bed.gz",
        "resources/MSSNG.variants.2020-07-24/mic_events.non-neural-high.300nt_flanking.hg38.bed.gz",
        "resources/MSSNG.variants.2020-07-24/mic_events.putative.300nt_flanking.hg38.bed.gz",
        "resources/MSSNG.variants.2020-07-24/non_mic_event.neural-high.300nt_flanking.hg38.bed.gz"
    ],
    metadata=get_metadata_path
  output:
    html_output=report("results/tcag/statistics/variant_counts.{subset}.html", caption="../report/variant_statistics.rst"),
    csv_output="results/tcag/statistics/variant_counts.{subset}.csv"
  conda:
    "../../microexon-code-workflow/workflow/envs/base_env.yml"
  script:
    "../scripts/collect_variant_statistics.py"


rule collect_all_variant_statistics:
  input:
    expand(rules.collect_variant_statistics.output, subset=variant_subsets)


rule sort_variants:
  """Collect variants into files based on ASD affection and type of affected microexons"""
  input:
    variants=variant_files_for_subset,
    mic_events="resources/MSSNG.variants.2020-07-24/{mic_set}.hg38.bed.gz",
    metadata=get_metadata_path
  wildcard_constraints:
    subset="[A-Z_123]+"
  conda:
    "../envs/pybedtools.yaml"
  output:
   "results/tcag/sorted_variants/{subset}.{mic_set}.tsv.gz"
  script:
    "../scripts/sort_variants.py"


rule sort_variants_for_all_mic_sets:
  input:
    expand(
      rules.sort_variants.output[0],
      mic_set=[
        "mic_events.neural-high.300nt_flanking",
        "mic_events.neural-high.flanking_exons.300nt_flanking",
        "mic_events.non-neural-high.300nt_flanking",
        "mic_events.putative.300nt_flanking",
        "non_mic_event.neural-high.300nt_flanking"
      ],
      subset=set(SUBSET)
    )


rule sort_novel_mic_variants:
  """Separate rule to sort variants that apply to novel microexons"""
  input:
    variants=variant_files_for_subset,
    mic_events="resources/MSSNG.variants.2020-07-24/mic_events.putative.300nt_flanking.hg38.bed.gz",
    novel_mic_list="resources/novel_microexon_hg38/hannes_total.2M.hg38.neuronal_expanded.txt",
    metadata=get_metadata_path
  wildcard_constraints:
    subset="[A-Z_123]+"
  conda:
    "../envs/pybedtools.yaml"
  output:
   "results/tcag/sorted_variants/{subset}.mic_events.novel.300nt_flanking.tsv.gz"
  script:
    "../scripts/sort_novel_mic_variants.py"


rule prepare_model_input_files:
  input:
    variant_files=variant_files_for_subset,
    mic_events="resources/inputs/MicEvents.hg38.tab.gz",
    metadata=get_metadata_path,
    novel_mic_list="resources/novel_microexon_hg38/hannes_total.2M.hg38.neuronal_expanded.txt"
  output:
    "resources/inputs/tcag/MicEvents.hg38.tcag_variants.{subset}.hg38.tab.gz"
  wildcard_constraints:
    subset="[A-Z_123]+"
  conda:
    "../envs/base_env.yml"
  script:
    "../scripts/prepare_model_input_files.py"


rule split_model_input_files:
  input:
    rules.prepare_model_input_files.output[0]
  output:
    expand("resources/inputs/tcag/MicEvents.hg38.tcag_variants.split_{split:02d}.{{subset}}.hg38.tab.gz", split=range(50))
  conda:
    "../envs/base_env.yml"
  script:
    "../scripts/split_tsv.py"

rule tcag_predictions:
  input:
    expand("results/predictions/{model}/tcag/predictions.MicEvents.hg38.tcag_variants.split_{split:02d}.{subset}.hg38.csv.gz",
           split=range(50), subset=variant_subsets, model=["known_microexons", "novel_microexons", "original_model"])
