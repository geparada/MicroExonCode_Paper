from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()
FTP = FTPRemoteProvider()

rule download_genome_2bit:
  """Download the genomes"""
  input:
    lambda wildcards: HTTP.remote(
      f"hgdownload.soe.ucsc.edu/goldenPath/{wildcards.species}/bigZips/{wildcards.species}.2bit", insecure=True)
  output:
    "resources/genomes/{species}.2bit"
  wildcard_constraints:
    species="\w{2,6}\d{1,2}"
  shell:
    "mv {input} {output}"
