import wolf
import dalmatian
import pandas as pd
from wolf.localization import LocalizeToDisk
from wolf.fc import WorkspaceInputConnector, SyncToWorkspace


#####
#####
# Unfortunately some @PG command tags inherited the sample IDs from terra
# which are not the same as in dbGAP. To homogenize, we perl replace terra sample
# id with the walkup id.
# Thankfully doesn't affect @RG and is easy to fix
#####
#####

WORKSPACE = "broad-getzlab-mm-germline-t/MM_WGS_SU2C_2021_10_15"

_DEBUG = False

WIC = WorkspaceInputConnector(WORKSPACE)

S = WIC.samples

bucket = "gs://2-large-wgs/hg38-dbgap"

docker="gcr.io/broad-getzlab-workflows/samtools@sha256:8074df347e20ca7f39646914eb9fcccde30a9ca85d4dd10b60fe99d5b78a223d"

S = S.loc[~pd.isna(S["dbgap_SAMPLE_ID"])]
S = S.loc[pd.isna(S["hg38_dbgap_ready_bam"])]

if _DEBUG:
    S = S.iloc[[0]]


def reheader_workflow(bam, bai, walkup_id, sample_name, catissue, terra_sample=None, terra=None, bucket=None):
    local_bam = LocalizeToDisk(files={"bam":bam, "bai":bai})
    reheadered_bam = wolf.Task(
            name="reheader",
            inputs={"bam":local_bam["bam"], "bai":local_bam["bai"], "walkup_id":walkup_id, "sample_name":sample_name, "catissue":catissue},
            script="""
            set -euxo pipefail
            out="${sample_name}.bam"
            ln -s $bam bam.bam
            ln -s $bai bam.bai
            samtools reheader -c 'perl -pe "s/${catissue}/${walkup_id}/g"' bam.bam > $out
            samtools index -@7 -o $(basename ${out}).bai $out
            rm bam.bam bam.bai
            samtools idxstats $out > $(basename ${out}).idxstats
            """,
            outputs={"bam":"*.bam","bai":"*.bai", "idxstats":"*.idxstats"},
            docker=docker,
            resources={"mem": "1G", "cpus-per-task": "8"},
            use_scratch_disk=True
            )
    # doesn't make sense not to have a bucket set. should make mandatory?
    if bucket:
        # dividing because not sure how to grab specific [0, 1, ...] index without a name?
        # do they come back in same order as sent?
        bam_cloud_path = wolf.UploadToBucket(files=[reheadered_bam["bam"]], bucket=bucket)
        bai_cloud_path = wolf.UploadToBucket(files=[reheadered_bam["bai"]], bucket=bucket)
        if terra:
            terra_sync = wolf.fc.SyncToWorkspace(nameworkspace=terra,
                                                 entity_type="sample",
                                                 entity_name=terra_sample,
                                                 attr_map={"hg38_dbgap_ready_bam":bam_cloud_path["cloud_path"],
                                                           "hg38_dbgap_ready_bai":bai_cloud_path["cloud_path"],
                                                           "hg38_idxstats":reheadered_bam["idxstats"]})
        if not _DEBUG:
            wolf.localization.DeleteDisk(
                name="DeleteBam",
                inputs={
                    "disk": reheadered_bam["bam"],
                    "upstream": terra_sync if terra is not None else bam_cloud_path
                })
        return bam_cloud_path, bai_cloud_path
    return None

with wolf.Workflow(workflow=reheader_workflow) as w:
    for sample, s in S.iterrows():
        w.run(
                run_name=s["dbgap_SAMPLE_ID"],
                walkup_id=s["walkup_id"],
                sample_name=s["dbgap_SAMPLE_ID"],
                catissue=sample,
                bam=s["hg38_analysis_ready_bam"],
                bai=s["hg38_analysis_ready_bam_index"],
                terra_sample=sample,
                terra=WORKSPACE,
                bucket=bucket)

