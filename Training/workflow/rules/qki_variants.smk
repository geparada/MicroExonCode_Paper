# subworkflow microexon_prediction:
#   workdir: "../microexon-code-workflow/workflow"

rule generate_qki_variants:
  input:
    "resources/inputs/MicEvents.{genome}.tab.gz",
    "resources/genomes/{genome}.2bit"
  output:
    "resources/inputs/qki/MicEvents.{genome}.qki.{motif}.{genome}.tab.gz"
  params:
    window_size=150,
    event_label_filter="hi"
  wildcard_constraints:
    motif="[ACGTWSMKRYBDHVN]+"
  conda:
    "../envs/base_env.yml"
  script:
    "../scripts/generate_qki_variants.py"

rule split_qki_input_files:
  input:
    rules.generate_qki_variants.output[0]
  output:
    expand("resources/inputs/qki/MicEvents.{{genome}}.qki.split_{split:02d}.{{motif}}.{{genome}}.tab.gz",
      split=range(10))
  conda:
    "../envs/base_env.yml"
  script:
    "../scripts/split_tsv.py"


rule qki_predictions:
  input:
    expand("results/qki/predictions.MicEvents.{genome}.qki.split_{split:02d}.{motif}.{genome}.csv.gz",
      split=range(10), motif=["ACTAA", "ACTAAY"], genome=["hg38", "mm10"])
