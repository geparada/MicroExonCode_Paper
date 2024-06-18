storage:
    provider="http",

rule download_genome_2bit:
  """Download the genomes"""
  input:
    lambda wildcards: storage.http(
      f"http://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.species}/bigZips/{wildcards.species}.2bit")
  output:
    "resources/genomes/{species}.2bit"
  wildcard_constraints:
    species="\w{2,6}\d{1,2}"
  shell:
    "mv {input} {output}"
