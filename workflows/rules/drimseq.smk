rule drimseq:
    input:
        # f"data/processed/{config['run_name']}/{{sample}}/{{sample}}_salmon/quant.sf"
    output:
        temp(f"results/{config['run_name']}/DRIMSeq/genes_DTU.tsv"),
        temp(f"results/{config['run_name']}/DRIMSeq/transcripts_proportions.tsv"),
        temp(f"results/{config['run_name']}/DRIMSeq/StageR.tsv"),

    params:
        samples = f"data/raw/samples_new.csv",
        out_path = config["run_name"],
        gtf = f"/scratch/project/open-25-50/gencode/gencode.v36.annotation.gtf",
        alpha = 0.05,
        mge = 10,
        mfe = 10,
        mfp = 0.1,
    log:
        #log = "data/logs/drimseq.log"
    conda:
        "../envs/drimseq.yml"
    shell:
        r"""
        Rscript scripts/drimseq.R --input1 {params.samples} --input2 {params.gtf} --out_path {params.out_path} --alpha {params.alpha} --mge {params.mge} --mfe {params.mfe} --mfp {params.mfp}
        """