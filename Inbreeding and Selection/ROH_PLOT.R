library(dplyr)
library(gridExtra)
library(stringr)

#Dingoes and Reference
roh <- read.csv('FILENAME', sep="") #Replace FILENAME with path and name of .hom plink output file 
info <- read.csv('PATH/ROH_sample.txt', sep="") #Replace PATH with path to ROH_sample.txt

genome_size <- 2203764842 
roh <- merge(roh, info[c(1,2)], by='IID')

######################
## ROH Count/Length ##
######################

ROH_num <- roh %>% 
  group_by(IID) %>%
  summarise(NROH = n())

ROH_length <- roh %>%
  group_by(IID) %>%
  summarise(ROH_Total=sum(KB))

#Merge data
ROH_final <- merge(merge(ROH_num, ROH_length[c(1,2)],  by='IID'), info[c(1,2)], by='IID')
ROH_final <- ROH_final %>% filter(Group!='Outgroup')

#Estimate FROH for all ROH size classes
ROH_final$FROH <- (ROH_final$ROH_Total*1000)/genome_size #Converts from KBs

##################
## Size Classes ##
##################

ROH_class <- roh %>% 
  group_by(IID)

ROH_class <- ROH_class %>%
  mutate(size_class = case_when(
    KB < 1000 ~ "<1 Kb",
    KB >= 1000 & KB < 2000 ~ "1-2 Kb",
    KB >= 2000 & KB < 5000 ~ "2-5 Kb",
    KB >= 5000 ~ ">5 Kb",
    TRUE ~ "Other"
  ))

ROH_class <- ROH_class %>%
  mutate(size_class = case_when(
    KB < 1000 ~ "<1 Kb",
    KB >= 1000 & KB < 2000 ~ "1-2 Kb",
    KB >= 2000 ~ ">2 Kb",
    TRUE ~ "Other"
  ))

#Calculate mean ROH lengths
ROH_class %>% group_by(FID) %>% mutate(mean_length = mean(KB)) %>%
  distinct(IID, .keep_all = TRUE)

ROH_class_counts_ind <- ROH_class %>%
  group_by(IID, Group, size_class) %>%
  summarize(count = n(), .groups = 'drop')

ROH_class_counts_all <- ROH_class_counts_ind %>%
  group_by(Group, size_class) %>%
  summarize(mean_count = mean(count),
            se_count = sd(count) / sqrt(n()), .groups = 'drop')

desired_order <- c("<1 Kb","1-2 Kb",">2 Kb")
ROH_class_counts_all$size_class <- factor(ROH_class_counts_all$size_class, levels = desired_order)
ROH_class_counts_all <- ROH_class_counts_all %>% arrange(size_class)

ROH_class_counts_dingo <- ROH_class_counts_all %>%
  filter(str_detect(Group, "Dingo|NewGuinea|European") & !str_detect(Group, "Aus_European"))

ROH_class_proportion <- ROH_class %>%
  group_by(Group, size_class) %>%
  summarise(count = n()) %>%
  group_by(Group) %>%
  mutate(proportion = count / sum(count) * 100)

merged_df <- ROH_class_proportion %>%
  left_join(pop_num, by = "Group")

merged_df <- merged_df %>%
  mutate(mean = count.x / count.y)

desired_order <- c("<1 Kb","1-2 Kb","2-5 Kb",">5 Kb")
ROH_class_proportion$size_class <- factor(ROH_class_proportion$size_class, levels = desired_order)
ROH_class_proportion <- ROH_class_proportion %>% arrange(size_class)
