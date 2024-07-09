
setwd('~/Dropbox (Partners HealthCare)/2_Projects/ms-wgs/code/')

source("utils.R")

# Load sample details -----------------------------------------------------

# Clinical annotation contains sensitive data - no github upload
clinical.annotation <- read_xlsx("../DuttaAlberge_CleanClinicalRef.xlsx")

# Terra sheets must remain controlled-access - no github upload
# hg38
samples.su2c <- read_xlsx("../annot/AlbergeDutta_TerraSheets.xlsx", sheet = "su2c_sample")
participants.su2c <- read_xlsx("../annot/AlbergeDutta_TerraSheets.xlsx", sheet = "su2c_participant")
pairs.su2c <- read_xlsx("../annot/AlbergeDutta_TerraSheets.xlsx", sheet = "su2c_pair")

# hg38 non-su2c (mark + oben + wtc)
samples.nonsu2c <- read_xlsx("../annot/AlbergeDutta_TerraSheets.xlsx", sheet = "non_su2c_sample")
participants.nonsu2c <- read_xlsx("../annot/AlbergeDutta_TerraSheets.xlsx", sheet = "non_su2c_participant")
pairs.nonsu2c <- read_xlsx("../annot/AlbergeDutta_TerraSheets.xlsx", sheet = "non_su2c_pair")

# hg19 ctc
samples.ctc <- read_xlsx("../annot/AlbergeDutta_TerraSheets.xlsx", sheet = "ctc_sample")
participants.ctc <- read_xlsx("../annot/AlbergeDutta_TerraSheets.xlsx", sheet = "ctc_participant")
pairs.ctc <- read_xlsx("../annot/AlbergeDutta_TerraSheets.xlsx", sheet = "ctc_pair")

# participant ID overlap

sum(participants.su2c$`entity:participant_id` %in% participants.nonsu2c$`entity:participant_id`)
sum(participants.su2c$`entity:participant_id` %in% participants.ctc$`entity:participant_id`)
sum(participants.nonsu2c$`entity:participant_id` %in% participants.ctc$`entity:participant_id`)
# if all zero pursue

participants.su2c <- participants.su2c %>% select(-c(Gender))
participants.nonsu2c <- participants.nonsu2c %>% select(-c(Gender))
participants.ctc <- participants.ctc %>% select(-c(Gender, Stage, Age))


clinical.and.terra <- clinical.annotation %>%
  left_join(participants.su2c, by=c("Participant_ID"="entity:participant_id")) %>%
  left_join(participants.nonsu2c, by=c("Participant_ID"="entity:participant_id")) %>%
  left_join(participants.ctc, by=c("Participant_ID"="entity:participant_id")) %>%
  arrange(Participant_ID)


# check who is missing in su2c terra
clinical.annotation$Participant_ID[!(clinical.annotation$Participant_ID %in% c(participants.su2c$`entity:participant_id`, participants.ctc$`entity:participant_id`, participants.nonsu2c$`entity:participant_id`))]
# 2023 01 23 == only pM11272 missing - since CLIA / delivery we do hg19 here
# conversely

participants.su2c$`entity:participant_id`[!(participants.su2c$`entity:participant_id` %in% clinical.annotation$Participant_ID)]
# looks fine

# overlap between pairs of three cohorts
intersect( pairs.ctc$`entity:pair_id` , pairs.su2c$`entity:pair_id` )
intersect( pairs.nonsu2c$`entity:pair_id` , pairs.su2c$`entity:pair_id` )
intersect( pairs.nonsu2c$`entity:pair_id` , pairs.ctc$`entity:pair_id` )

# overlap between column names
intersect(colnames(clinical.and.terra), colnames(pairs.ctc))
intersect(colnames(clinical.and.terra), colnames(pairs.su2c))
intersect(colnames(clinical.and.terra), colnames(pairs.nonsu2c))
# "ploidy"  "purity"  "TiN_WGS" only in all three

# merge pair sets
three.pair.sets <- rbindlist(list(select(pairs.su2c, c(`entity:pair_id`, case_sample, control_sample, participant)),
                                  select(pairs.nonsu2c, c(`entity:pair_id`, case_sample, control_sample, participant)),
                                  select(pairs.ctc, c(`entity:pair_id`, case_sample, control_sample, participant))))

clinical.and.terra.pairs <- clinical.and.terra %>%
  left_join(three.pair.sets, by=c("Participant_ID"="participant"))

# add sample QC
samples.su2c$quick_cov_estimate
samples.nonsu2c$quick_cov_estimate
samples.ctc$quick_cov_estimate <- samples.ctc$wgs_coverage

# overlap between column names
intersect(colnames(clinical.and.terra), colnames(samples.ctc))
intersect(colnames(clinical.and.terra), colnames(samples.su2c))
intersect(colnames(clinical.and.terra), colnames(samples.nonsu2c))

# overlap between column names
intersect(samples.ctc$`entity:sample_id`, samples.su2c$`entity:sample_id`)
intersect(samples.ctc$`entity:sample_id`, samples.nonsu2c$`entity:sample_id`)
intersect(samples.su2c$`entity:sample_id`, samples.nonsu2c$`entity:sample_id`)

# TODO
# there's room for more QC directly in the TERRA sheets

three.samples.sets <- rbindlist(list(select(samples.ctc, c(`entity:sample_id`, participant, quick_cov_estimate)),
               select(samples.su2c, c(`entity:sample_id`, participant, quick_cov_estimate)),
               select(samples.nonsu2c, c(`entity:sample_id`, participant, quick_cov_estimate))))

clinical.and.terra.pairs.samples <- clinical.and.terra.pairs %>%
  mutate(MuTect_Patient=ifelse(!is.na(CTF_ID), CTF_ID, Participant_ID),
         `Sample ID`=Participant_ID) %>%
  left_join(three.samples.sets, by=c("Participant_ID"="participant", "case_sample"="entity:sample_id"), suffix=c("", "_case")) %>%
  left_join(three.samples.sets, by=c("Participant_ID"="participant", "control_sample"="entity:sample_id"), suffix=c("", "_control"))

# annotate detection power for SNVs clonal and subclonal

clinical.and.terra.pairs.samples <- clinical.and.terra.pairs.samples |>
  rowwise() |>
  mutate(power_ccf_1=calculate_tumor_power(depth=ifelse(is.na(quick_cov_estimate), 60, round(quick_cov_estimate)), # optimal would be to generate 10k random poisson counts from coverage 
                                           tlod = 6.13,
                                           eps = 1E-3,
                                           rho = purity,
                                           psi = 2,
                                           CCF = 1),
         power_ccf_0.5=calculate_tumor_power(depth=ifelse(is.na(quick_cov_estimate), 60, round(quick_cov_estimate)), # optimal would be to generate 10k random poisson counts from coverage 
                                           tlod = 6.13,
                                           eps = 1E-3,
                                           rho = purity,
                                           psi = 2,
                                           CCF = 0.5),
         )

# FIXME add secondary, serial, matched tissues pNULL# FIXME add secondary, serial, matched tissues pair_ids
# currently only reference is annotated
clinical.and.terra.ref.pairs.samples <- clinical.and.terra.pairs.samples %>%
  group_by(Participant_ID) %>%
  filter( n()==1 | ( n()>=2 & `entity:pair_id`==Reference_Pair ) )

# all this will be useful for softwares like MuTect that in the default pipeline name their tumor sample as pair + "_tumor"

# who's missing in reference
clinical.and.terra.ref.pairs.samples %>% filter(Participant_ID %nin% clinical.and.terra.ref.pairs.samples$Participant_ID) %>% pull(`entity:pair_id`)

write_tsv(clinical.and.terra.ref.pairs.samples, "../annot/2023_v0_participant_annotation.tsv")

# TODO list
# check reciprocally missing second / third sample


