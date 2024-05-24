import pandas as pd

double_mutant_motifs = [line.rstrip() for line in open("config/double_mutant_motifs.txt")]

rule double_mutants_generate_inputs:
  input:
    event_file="resources/MicEvents.hg38.neural.tab.gz",
    genome_file="resources/genomes/hg38.2bit"
  output:
    "resources/inputs/double_mutants/MicEvents.neural.double_mutants.{event_id}.hg38.tab.gz"
  params:
    motifs=double_mutant_motifs,
  conda:
    "../envs/base_env.yml"
  script:
    "../scripts/generate_double_mutant_variants.py"

def get_double_mutant_variant_event_ids(wildcards):
  events = pd.read_csv("resources/MicEvents.hg38.neural.tab.gz", sep="\t", index_col="event")
  events = events.loc[events.label == 'hi']
  event_ids = events.index
  return [f"results/double_mutants/predictions.MicEvents.neural.double_mutants.{event_id}.hg38.csv.gz"
          for event_id in event_ids]


rule double_mutant_predictions:
  input:
    get_double_mutant_variant_event_ids
