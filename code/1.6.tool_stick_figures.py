import wolf

def tool_stick_figures(
        maf,
        sig_genes,
        uniprot_id_map="gs://broad-institute-gdac/reference/UniProt/UniProtKB_AC_to_ID_human_Release_2017_10_25.shelf",
        uniprot_swiss_db="gs://broad-institute-gdac/reference/UniProt/Uniprot_SWISS_human_Release_2017_10_25.shelf",
        localize_maf=False):

    #maf can be quite massive with union of WGS samples
    if localize_maf:
        maf_disk=wolf.LocalizeToDisk(files={"maf":maf})
        maf_file=maf_disk["maf"]
    else:
        maf_file=maf

    tool_stick_figures=wolf.Task(
            name="tool_stick_figures",
            inputs={"maf":maf_file, "sig_genes":sig_genes, "uniprot_id_map":uniprot_id_map, "uniprot_swiss_db":uniprot_swiss_db, "package":"true" },
            script="""
            set -euxo pipefail
            python /src/generate_figs.py \
			--jpg \
			--pdf \
			${maf} \
			${sig_genes} \
			${uniprot_swiss_db} \
			${uniprot_id_map} \
			/src/mutfig

		if ${package}; then
			zip -r tool_stick_figures . -x \
				"fc-[a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9]-[a-f0-9][a-f0-9][a-f0-9][a-f0-9]-[a-f0-9][a-f0-9][a-f0-9][a-f0-9]-[a-f0-9][a-f0-9][a-f0-9][a-f0-9]-[a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9][a-f0-9]/*" \
				lost+found/\* \
				broad-institute-gdac/\* \
				"tmp.[a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9]/*" \
				exec.sh
		fi
            """,
            outputs={"tool_stick_figures_pkg":"tool_stick_figures.zip"},
            docker="broadgdac/tool_stick_figures:22" # author: D. Heiman
            )
    return tool_stick_figures

maf="gs://2-large-wgs/mafs/20230503.su2c.mark.jco.rsp.mmrf.hg19.nowm.nomgip.maf"
sig_genes="gs://2-large-wgs/SU2C_MutSig_Saturation_n_1030/outdir/sig_genes.txt"

with wolf.Workflow(
        workflow=tool_stick_figures,
        common_task_opts={
                       "retry": 5,
                       "cleanup_job_workdir": True
                       }) as w:
            w.run(RUN_NAME="SU2C_tool_stick_figures", maf=maf, sig_genes=sig_genes)
