# library(rhdf5)
# library(dplyr)
# library(tibble)
# library(ggplot2)
# library(lubridate)
# library(zoo)
# 
# .getSignalMetaData <- function(file) {
#     
#     dat <- h5readAttributes(file, "/UniqueGlobalKey/channel_id/")
#     dat$sampling_rate <- IONiseR:::.getSamplingRate(file)
#     dat$exp_start_time <- as.numeric(h5readAttributes(file, "/UniqueGlobalKey/tracking_id")$exp_start_time)
#     return( dat )
# }
# 
# .getEvents <- function(file) {
#     
#     fid <- H5Fopen(file)
#     on.exit(H5Fclose(fid))
#     
#     exists <- IONiseR:::.groupExistsObj(fid, group = "/Analyses/EventDetection_000/Reads")
#     if(!exists) {
#         start_time <- duration <- num_events <- median_signal <- NA
#     } else {
#         ## get the Read_No., this changes in every file
#         gid <- H5Gopen(fid, "/Analyses/EventDetection_000/Reads")
#         read_number_char <- h5ls(gid)[1,"name"]
#         H5Gclose(gid)
#         
#         ## Open the group
#         gid <- H5Gopen(fid, paste0("/Analyses/EventDetection_000/Reads/", read_number_char)) 
#         did <- H5Dopen(gid, "Events")
#         events <- as_tibble(H5Dread(did, bit64conversion = "int", compoundAsDataFrame = TRUE))
#         H5Dclose(did)
#         
#         H5Gclose(gid)
#     } 
#     
#     return(events)
# }
# 
# .getRaw <- function(file) {
#     
#     fid <- H5Fopen(file)
#     on.exit(H5Fclose(fid))
#     
#     exists <- IONiseR:::.groupExistsObj(fid, group = "/Raw/Reads")
#     if(!exists) {
#         start_time <- duration <- num_events <- median_signal <- NA
#     } else {
#         ## get the Read_No., this changes in every file
#         gid <- H5Gopen(fid, "/Raw/Reads")
#         read_number_char <- h5ls(gid)[1,"name"]
#         H5Gclose(gid)
#         
#         ## Open the group
#         gid <- H5Gopen(fid, paste0("/Raw/Reads/", read_number_char)) 
#         did <- H5Dopen(gid, "Signal")
#         signal = H5Dread(did)
#         H5Dclose(did)
#         
#         ## get the starting time
#         aid <- H5Aopen(gid, name = "start_time")
#         startTime <- H5Aread(aid)
#         H5Aclose(aid)
#         H5Gclose(gid)
#         
#         raw <- tibble(signal, time = seq(startTime, length.out = length(signal)))
#     } 
#     
#     return(raw)
# }
# 
# # [experiment start time]	/UniqueGlobalKey/tracking_id/{exp_start_time}
# # [read start time]	/Raw/Reads/Read_<nnn>/{start_time}
# # [sampling rate]	/UniqueGlobalKey/channel_id/{sampling_rate}
# # [channels 0pA adc]	/UniqueGlobalKey/channel_id/{offset}
# # [digitisable range in pA]	/UniqueGlobalKey/channel_id/{range}
# # [digitisation]	/UniqueGlobalKey/channel_id/{digitisation}
# # x-axis: [signal time in UTC] = [experiment start time] + ([read start time] + <signal_index_0_based>) / [sampling rate]
# # y-axis: [current in pA] = (<signal_value> + [channels 0pA adc] ) * [digitisable range in pA] / [digitisation]
# # <signal_value> is the column value in the dataset at /Raw/Reads/Read_###/Signal
# # <signal_index_0_based> is the row index starting at 0 for the time conversion formula to work
# 
# f1 <- "~/projects/bioconductor/IONiseR/inst/extdata/example_v2.fast5"
# f1 <- "/mnt/data/randomEncounters/MichaelVeronesi/Nanopore Sample 2/graveley_HP_20160627_FN_MN16664_sequencing_run_sample_id_84084_ch276_read1110_strand.fast5"
# s1 <- .getSignalMetaData(f1)
# t1 <- .getRaw(f1) %>%
#     mutate(pA = signal + s1$offset * (s1$range / s1$digitisation),
#            utc = (time / s1$sampling_rate) + as.numeric(s1$exp_start_time))
# 
# slice(t1,500:5000 ) %>%
#     mutate(date = as_datetime(utc)) %>%
# ggplot( aes(x = date, y = pA)) + geom_line()
# 
# t2 <- slice(t1,2500:4000 ) %>%
#     mutate(date = as_datetime(utc)) %>%
#     mutate(mean.back = rollmean(x = pA, 15, align = "right", fill = NA),
#            mean.forward = rollmean(x = pA, 15, align = "left", fill = NA)) %>%
#     mutate(mad.back = 3 * rollapply(data = pA, width = 15, align = "right", fill = NA, FUN = mad)) %>%
#     mutate(mad.forward = 3 * rollapply(data = pA, width = 15, align = "left", fill = NA, FUN = mad)) %>%
#     mutate(jump = ifelse((pA<mean.back-mad.back|pA>mean.back+mad.back) &
#                         (pA>mean.forward-mad.forward&pA<mean.forward+mad.forward), 
#                         TRUE, FALSE)) %>%
#     mutate(jump2 = rollapply(data = jump, width = 15, align = "right", fill = FALSE, FUN = function(x) { ifelse(sum(x, na.rm=TRUE)==1&x[15], TRUE, FALSE)}))
# 
#     
# ggplot(t2, aes(x = date, y = pA)) + geom_point() +
#     geom_point(data = filter(t2, jump2==TRUE), aes(x = date, y = pA), color = "cyan")
