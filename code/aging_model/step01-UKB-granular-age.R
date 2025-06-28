source("code/functions/packages.R")
source("code/functions/utils.R")

data_path <- "data/ukb"
out_path <- "data/ukb"

patient_data_01 <- read.csv(
  file.path(data_path, "data_participant.csv"),
  row.names = 1
)
patient_data_02 <- read.csv(
  file.path(data_path, "p53.csv"),
  row.names = 1
)
common_ids <- intersect(
  rownames(patient_data_01),
  rownames(patient_data_02)
)

patient_data_01 <- patient_data_01[common_ids, ]
patient_data_02 <- patient_data_02[common_ids, ]
patient_data <- cbind(patient_data_01, patient_data_02)
patient_data$eid <- rownames(patient_data)

cols <- c(
  "eid",
  "p34",
  "p52",
  "p53_i0",
  "p53_i2",
  "p53_i3",
  "p21022"
)

patient_data$birth_year <- patient_data$p34
patient_data$birth_month <- patient_data$p52
patient_data$recruitment_date <- patient_data$p53_i0
patient_data$recruitment_date_ins2 <- patient_data$p53_i2
patient_data$recruitment_date_ins3 <- patient_data$p53_i3
patient_data$recruitment_age <- patient_data$p21022

patient_data$recruitment_date <- as.Date(
  patient_data$recruitment_date,
  format = "%Y-%m-%d"
)
patient_data$recruitment_date_ins2 <- as.Date(
  patient_data$recruitment_date_ins2,
  format = "%Y-%m-%d"
)
patient_data$recruitment_date_ins3 <- as.Date(
  patient_data$recruitment_date_ins3,
  format = "%Y-%m-%d"
)

# make estimated birth date a first day of birth month and year
patient_data$est_birth_date <- ymd(
  paste0(
    patient_data$birth_year, "-",
    patient_data$birth_month, "-01"
  )
)
patient_data$est_birth_date <- as.Date(
  patient_data$est_birth_date,
  format = "%Y-%m-%d"
)

# make granular age as difference between recruitment date and estimated birth date
patient_data$age_granular <- as.numeric(
  patient_data$recruitment_date - patient_data$est_birth_date
)
patient_data$age_granular <- patient_data$age_granular / 365.25

# code age at 2014+ imaging visit (ins = 2)
patient_data$time2 <- as.numeric(
  patient_data$recruitment_date_ins2 - patient_data$recruitment_date
) / 365.25
patient_data$age_granular_ins2 <- patient_data$age_granular + patient_data$time2

# code age at 2019+ repeat imaging visit (ins = 3)
patient_data$time3 <- as.numeric(
  patient_data$recruitment_date_ins3 - patient_data$recruitment_date
) / 365.25
patient_data$age_granular_ins3 <- patient_data$age_granular + patient_data$time3

cols <- c(
  "eid",
  "birth_year",
  "birth_month",
  "recruitment_date",
  "recruitment_date_ins2",
  "recruitment_date_ins3",
  "recruitment_age",
  "est_birth_date",
  "age_granular",
  "age_granular_ins2",
  "age_granular_ins3"
)

data <- as.data.frame(patient_data)[cols]

data <- data[which(!is.na(data$recruitment_age)), ]

hist(data$recruitment_age)
hist(data$age_granular)

summary(data$recruitment_age)
summary(data$age_granular)

write.csv(
  data,
  row.names = FALSE,
  file.path(out_path, "granular_age_april_07_2025.csv")
)
