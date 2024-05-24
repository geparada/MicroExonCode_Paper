from snakemake.remote.SFTP import RemoteProvider
SFTP = RemoteProvider(username=config["sftp_username"], password=config["sftp_password"])


VCF_FILTER = (
  "'FILTER=\"PASS\" &" 
  "  FORMAT/DP>=10 &" 
  "  ((GT=\"het\" & (TYPE=\"snp\" | TYPE=\"mnp\") & GQ>=99 & AD[:1]/(FORMAT/DP)>=0.3 & AD[:1]/(FORMAT/DP)<0.8) | "
  "   (GT=\"het\" & TYPE=\"indel\" & GQ>=90 & AD[:1]/(FORMAT/DP)>=0.3 & AD[:1]/(FORMAT/DP)<0.8) | "
  "   (GT=\"AA\" & GQ>=25 & AD[:1]/(FORMAT/DP)>=0.8))'"
)

rule download_tcag:
  input:
    vcf=SFTP.remote("tcagfts.ccm.sickkids.ca/inbox/{dataset}/{sample}.vcf.gz"),
    index = SFTP.remote("tcagfts.ccm.sickkids.ca/inbox/{dataset}/{sample}.vcf.gz.tbi")
  output:
    vcf=temp("resources/tcag_full_genotypes/vcf/{dataset}/{sample}.vcf.gz"),
    index=temp("resources/tcag_full_genotypes/vcf/{dataset}/{sample}.vcf.gz.tbi")
  resources:
    parallel_genomes=1
  shell:
    "cp {input.vcf} {output.vcf} && cp {input.index} {output.index}"

rule consensus:
  input:
    g="resources/genomes/hg38.fa.gz",
    vcf=ancient(rules.download_tcag.output.vcf),
    index=ancient(rules.download_tcag.output.index)
  output:
    twobit=temp("resources/genomes/[{dataset}--{sample}].2bit"),
    chain="results/tcag_full_genotypes/consensus/{dataset}/chain/{sample}.chain"
  log:
    "results/tcag_full_genotypes/consensus/{dataset}/logs/{sample}.txt",
  conda: "../envs/bcftools.yml"
  # group: "compute_consensus_genome"
  shell:
    "bcftools consensus "
    "  {input.vcf} "
    "  -f {input.g} "
    "  -H A "
    "  -o $TMPDIR/{wildcards.sample}.fa " 
    "  -c {output.chain} " 
    f"  -i  {VCF_FILTER} "
    "2> {log} && "
    "faToTwoBit $TMPDIR/{wildcards.sample}.fa {output.twobit} && "
    "rm $TMPDIR/{wildcards.sample}.fa"

rule intersect_variants:
  input:
    vcf=ancient(rules.download_tcag.output.vcf),
    index=ancient(rules.download_tcag.output.index),
    # vcf="resources/tcag_full_genotypes/vcf/{sample}.vcf.gz",
    # index="resources/tcag_full_genotypes/vcf/{sample}.vcf.gz.tbi",
    events_bed="resources/tcag_full_genotypes/wt_bed/MicEvents_known_novel.feature_intervals.hg38.bed"
  output:
    "results/tcag_full_genotypes/bed/MicEvents_known_novel.variant_overlaps.[{dataset}--{sample}].bed"
  conda: "../envs/bcftools.yml"
  # group: "compute_consensus_genome"
  shell:
    f"bcftools view {{input.vcf}} -i {VCF_FILTER} | "
     "bedtools intersect -a {input.events_bed} -b stdin -loj > {output}"

rule postprocess_intersect_variants:
  input:
    rules.intersect_variants.output[0]
  output:
    "results/tcag_full_genotypes/variant_intersection/MicEvents_known_novel.variants_overlaps.[{dataset}--{sample}].tsv"
  conda: "../envs/pybedtools.yaml"
  # priority: 1000
  script: "../scripts/postprocess_intersect_variants.py"


rule liftover:
  input:
    chain=rules.consensus.output.chain,
    up="resources/tcag_full_genotypes/wt_bed/MicEvents_known_novel.hg38.upInt.bed",
    dn="resources/tcag_full_genotypes/wt_bed/MicEvents_known_novel.hg38.dnInt.bed"
  output:
    up="results/tcag_full_genotypes/bed/MicEvents_known_novel.[{dataset}--{sample}].upInt.bed",
    dn="results/tcag_full_genotypes/bed/MicEvents_known_novel.[{dataset}--{sample}].dnInt.bed",
    up_un="results/tcag_full_genotypes/bed/MicEvents_known_novel.[{dataset}--{sample}].upInt.unmapped.bed",
    dn_un="results/tcag_full_genotypes/bed/MicEvents_known_novel.[{dataset}--{sample}].dnInt.unmapped.bed"
  conda : "../envs/bcftools.yml"
  # group: "compute_consensus_genome"
  shell:
    """liftOver {input.up} {input.chain} {output.up} {output.up_un}
       liftOver {input.dn} {input.chain} {output.dn} {output.dn_un}"""


rule bed_to_mic_events:
  input:
    wt_events="resources/inputs/tcag_full_genotypes/MicEvents_known_novel.hg38.tab.gz",
    up_bed=rules.liftover.output.up,
    dn_bed=rules.liftover.output.dn,
    up_un_bed=rules.liftover.output.up_un,
    dn_un_bed=rules.liftover.output.dn_un,
    wt_genome="resources/genomes/hg38.2bit",
    mt_genome=rules.consensus.output.twobit
  output:
    mt_events="resources/inputs/tcag_full_genotypes/MicEvents_known_novel.[{dataset}--{sample}].tab.gz",
    skipped_events="resources/inputs/tcag_full_genotypes/MicEvents_known_novel.skipped_events.[{dataset}--{sample}].csv"
  params:
    intron_window=300,
    exon_window=100
  conda: "../envs/base_env.yml"
  # group: "compute_consensus_genome"
  script: "../scripts/bed_to_mic_events.py"


def get_all_samples(wildcards):
  yield "results/predictions/known_microexons/tcag_full_genotypes/predictions.MicEvents_known_novel.hg38.csv.gz"
  # for dataset in ["MSSNG_ILMN"]:
  for dataset in ["MSSNG_CG", "MSSNG_ILMN", "SSC", "SPARK_WGS_1", "SPARK_WGS_2", "SPARK_WGS_3"]:
    with open(f"resources/tcag_full_genotypes/{dataset}_children_ids.txt", "rt") as f:
      for sample in f:
        yield "results/predictions/known_microexons/tcag_full_genotypes/" \
              f"predictions.MicEvents_known_novel.[{dataset}--{sample.rstrip()}].csv.gz"
        yield "results/tcag_full_genotypes/variant_intersection/" \
              f"MicEvents_known_novel.variants_overlaps.[{dataset}--{sample.rstrip()}].tsv"

rule tcag_predictions_for_dataset:
  input: get_all_samples

EVENTS, = glob_wildcards("resources/inputs/tcag_single_variant_analysis/{event}_single_variants.hg38.tab.gz")

rule tcag_single_variant_predictions:
  input:
    expand(
      "results/predictions/known_microexons/tcag_single_variant_analysis/"
      "predictions.{event}_single_variants.hg38.csv.gz",
      event=EVENTS
    )

SSC_EVENTS, = glob_wildcards("resources/inputs/tcag_single_variant_analysis/{event}_single_variants.hg38.tab.gz")

rule ssc_single_variant_predictions:
  input:
    expand(
      "results/predictions/known_microexons/ssc_single_variant_analysis/"
      "predictions.{event}_variants.hg38.csv.gz",
      event=SSC_EVENTS
    )

