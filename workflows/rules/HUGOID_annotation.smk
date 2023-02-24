rule HUGOID_annotation:
    input:
        f"results/{config['run_name']}/DRIMSeq/genes_DTU.tsv",
        f"results/{config['run_name']}/DRIMSeq/transcripts_proportions.tsv",
        f"results/{config['run_name']}/DRIMSeq/StageR.tsv",
    output:
        f"results/{config['run_name']}/DRIMSeq/genes_DTU.csv",
        f"results/{config['run_name']}/DRIMSeq/transcripts_proportions.csv",
        f"results/{config['run_name']}/DRIMSeq/StageR.csv"
    params:
        DTU = f"results/{config['run_name']}/DRIMSeq/genes_DTU.tsv",
        PROP = f"results/{config['run_name']}/DRIMSeq/transcripts_proportions.tsv",
        STAG = f"results/{config['run_name']}/DRIMSeq/StageR.tsv",
        out_path = config["run_name"],
    conda:
        "../envs/drimseq.yml"
    shell:
        r"""
        Rscript scripts/HUGOID_annotation.R --input1 {params.DTU} --input2 {params.PROP} --input3 {params.STAG} --out_path {params.out_path}
        """
