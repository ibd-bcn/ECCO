## code to prepare `batch` dataset goes here
library("tidyverse")
library("readxl")

files_raw <- read.table("data-raw/list_fastq.txt", header = FALSE,
                        stringsAsFactors = FALSE)
files_p <- gsub("_read[12].fastq.gz", "", files_raw$V1) %>%
    sort()


dim(files_p) <- c(2, 292)
bat <- files_p %>%
    t() %>%
    as.data.frame(stringAsFactors = FALSE) %>%
    rename(read1 = V1, read2 = V2) %>%
    mutate(file = basename(as.character(read1)),
           read1 = paste0(as.character(read1), "_read1.fastq.gz"),
           read2 = paste0(as.character(read2), "_read2.fastq.gz"),
           folder = dirname(read1),
           sample = str_extract(file, "[:alnum:]*_"),
           folder = gsub("\\./", "", folder),
           f = folder,
           sample = gsub("_", "", sample),
           ) %>%
    pivot_wider(names_from = folder, values_from = folder) %>%
    mutate()


dup_samples <- bat %>%
    group_by(sample) %>%
    count() %>%
    filter(n != 1) %>%
    pull(sample)

bat %>%
    filter(sample %in% dup_samples) %>%
    arrange(sample) %>%
    write_csv(path = "data-raw/duplicated_files.csv")

exc <- read_excel("data-raw/For Juanjo.xlsx", n_max = 162)


missing <- exc$`AGILENT/ARRAY CODE`[!exc$`AGILENT/ARRAY CODE` %in% bat$sample]
length(missing) # 8 samples missing (quality control)

# Check the number of samples per group
exc %>%
    mutate(GROUP = toupper(GROUP)) %>%
    filter(
        !`AGILENT/ARRAY CODE` %in% missing &
            GROUP != "UC"
    ) %>%
    group_by(GROUP, TYPE, LOCATION, `SAMPLE STATUS`) %>%
    count() %>%
    arrange(n)

# At least 6 per group
exc %>%
    mutate(GROUP = toupper(GROUP)) %>%
    filter(
        !`AGILENT/ARRAY CODE` %in% missing &
            GROUP != "UC"
    ) %>%
    group_by(GROUP, TYPE, LOCATION) %>%
    count() %>%
    arrange(n)

exc %>%
    mutate(GROUP = toupper(GROUP)) %>%
    filter(
        !`AGILENT/ARRAY CODE` %in% missing &
            GROUP != "UC"
    ) %>%
    mutate(type = case_when(endsWith(`AGILENT/ARRAY CODE`, "V") ~ "OLD",
                            TRUE ~ "NEW")) %>%
    count(type)

usethis::use_data(batch, overwrite = TRUE)
