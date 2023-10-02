#set working directory
setwd("")

#in the original spreadsheet, I changed the 0-4 scale codes of 2015-2022 to the following scale:
# 0=0
# 0;=1
# ;=2
# ;C=3
# 1;=4
# 1=5
# 2=6
# 3=7
# 3+=8
# 4=9

#this conversion is based on miller et al 2021

library(tidyverse) #several packages including ggplot2 and dplyr

#read in data
raw_OCR_survey <- read.csv("Summary_1992-2022.csv")

# check data
str(raw_OCR_survey)
table(raw_OCR_survey$HiFi)
table(raw_OCR_survey$TAM.O.405)

# change the integer type to numeric
raw_OCR_survey$TAM.O.405 <- as.numeric(raw_OCR_survey$TAM.O.405)

#clean up the isolate names to the get rid of "-" between state and number for some of the years (95LA-026 -> 95LA026)
raw_OCR_survey <- raw_OCR_survey %>%
  mutate(isolate= gsub("([A-Za-z])-", "\\1", isolate)) #only deleting "-" that follow a letter

# Doing a race clone correction for the data BEFORE changing from more comprehensive scale 0-4 or 0-9 to binary avir (0) or vir (1)
# if two or more single pustule isolates from the same collection (usually submitted as Pca leaf or single plant) ex. 16TX31-1 and 16TX31-2, have the
# exact same virulence profile then that is an artificial oversampling of that race.

# First, make a secondary ID column that only has the initial ID before the "-"
raw_OCR_survey_1 <- raw_OCR_survey %>% separate(isolate, c("new_id", "second"), sep = "-", remove=FALSE) %>% select(-second)

str(raw_OCR_survey_1)

# dropping all duplicate rows based on new_id and the differential line scores
raw_OCR_survey_1 <- raw_OCR_survey_1 %>% distinct(new_id,Pc14,Pc35,Pc36,Pc39,Pc40,Pc45,Pc46,Pc48,Pc50,Pc51,Pc52,Pc53,Pc54,Pc55,Pc56,Pc57,Pc58,Pc59,Pc60,Pc61,Pc62,Pc63,Pc64,Pc67,Pc68,Pc70,Pc71,Pc91,Pc94,Pc96,MARVELOUS,H548,IAB605Xsel.,WIX4361.9,TAM.O.405,Belle,HiFi,Leggett,Stainless, .keep_all = TRUE)

# dropping the new_id column
raw_OCR_survey <- raw_OCR_survey_1[,-2]

#split dataframe into two based on the year, exclude 1992 because the data was collected by a different person and only includes accessions from MN, so not really part of the national survey
OCR_before2015 <- filter(raw_OCR_survey, year < 2015 & year >=1993)
OCR_after2015 <- filter(raw_OCR_survey, year >= 2015)

table(OCR_before2015$Pc14)
table(OCR_after2015$Pc14)

#for before 2015, take columns 8-47 with the rating data and change to 1 if virulent (>=3) or 0 if avirulent
OCR_before2015[8:47] <- lapply(OCR_before2015[8:47], function(x) ifelse(x>=3, 1, 0))

#for 2015 and after, take columns 8-47 with the rating data and change to 1 if virulent (>=7) or 0 if avirulent
OCR_after2015[8:47] <- lapply(OCR_after2015[8:47], function(x) ifelse(x>=7, 1, 0))

# recombine the two datasets
OCR_survey <- rbind(OCR_before2015,OCR_after2015)

str(OCR_survey)

OCR_survey %>% count(region)

#for each isolate, calculate the average amount of virulences and the virulence count, append to the end of the dataframe
OCR_survey$meanvirulence = rowMeans(OCR_survey[,c(8:47)], na.rm = TRUE)
OCR_survey$countvirulence = rowSums(OCR_survey[,c(8:47)], na.rm = TRUE)


#filter out 2013 data because only two samples and add in later
OCR_wout_2013<- OCR_survey %>% filter(year!=2013)

#these two observations must have typos for values Pc91, Pc94,Pc96, IAB605Xsel., WIX4361.9, TAM.O.405, Belle, HiFi, Leggett, and Stainless
# because they were not being tested for at that time.
# replace with NAs to clean up the data set
OCR_wout_2013 %>%
  mutate(across("Pc91":"Pc96", ~ifelse(isolate=="96WI064", NA, .))) %>% 
  mutate(across("IAB605Xsel.":"TAM.O.405", ~ifelse(isolate=="96WI064", NA, .))) %>% 
  mutate(across("Pc91":"Pc96", ~ifelse(isolate=="01TX010", NA, .))) %>% 
  mutate(across("IAB605Xsel.":"TAM.O.405", ~ifelse(isolate=="01TX010", NA, .))) -> OCR_wout_2013_clean

#### Pca race analysis for all study years only using the 30 original lines and excluding the isolates with missing data #######


library(data.table)

OCR_strains_allminusNAs <- OCR_wout_2013_clean %>%
  select(-c(Belle, Stainless,Leggett,HiFi,IAB605Xsel.,WIX4361.9,TAM.O.405,Pc91,Pc94,Pc96)) %>% 
  na.omit()
setDT(OCR_strains_allminusNAs)[,group :=.GRP,by = .(Pc14,Pc35,Pc36,Pc38,Pc39,Pc40,Pc45,Pc46,Pc48,Pc50,Pc51,Pc52,Pc53,Pc54,Pc55,Pc56,Pc57,Pc58,Pc59,Pc60,Pc61,Pc62,Pc63,Pc64,Pc67,Pc68,Pc70,Pc71,H548,MARVELOUS)]

OCR_n_strains_allminusNAs <-OCR_strains_allminusNAs %>% 
  count(group) %>% 
  arrange(-n)

strain_counts_all <-OCR_n_strains_allminusNAs %>% 
  count(n, name= "frequency") %>% 
  # rename(nn=frequency) %>% 
  mutate(totalNofstrains= n*frequency) %>%
  mutate(percentage= totalNofstrains/sum(totalNofstrains)) %>% 
  mutate(percentageofstrains= frequency/sum(frequency))
sum(strain_counts_all$frequency)
sum(strain_counts_all$totalNofstrains)

isolate1688_all <-OCR_strains_allminusNAs %>% 
  filter(group==1688)
isolate1688_all %>% count(region, year)

isolate561_all <-OCR_strains_allminusNAs %>% 
  filter(group==561)
isolate561_all %>% count(region, year)

isolate1676_all <-OCR_strains_allminusNAs %>% 
  filter(group==1676)
isolate1676_all %>% count(region, year)

isolate26_all <-OCR_strains_allminusNAs %>% 
  filter(group==26)
isolate26_all %>% count(region, year)

# Sorting out and printing the isolates with the least amount of virulences
head(OCR_strains_allminusNAs[order(OCR_strains_allminusNAs$countvirulence),], n=20L)

# determining how many races are exclusive to each region when there are at least 5 isolates
# sort out any races with less than 5 isolates
OCR_5isolates_plus <- OCR_n_strains_allminusNAs[OCR_n_strains_allminusNAs$n >= 5]
OCR_5isolates_plus_data <-dplyr::filter(OCR_strains_allminusNAs, group %in% OCR_5isolates_plus$group)
str(OCR_5isolates_plus_data)
OCR_5isolates_plus_data %>% count(group,region)
by_year <-OCR_5isolates_plus_data %>% count(group,year)

#####

# For all races with more than one isolate, determine what percentage were sampled first in the North, South, or both in the same year
# sorting out any races with less than 2 isolates
OCR_2isolates_plus <- OCR_n_strains_allminusNAs[OCR_n_strains_allminusNAs$n >= 2] #this step is technically unnecessary but keep in case want to increase the filter for strains that are present 3 or more times
OCR_2isolates_plus_data <-dplyr::filter(OCR_strains_allminusNAs, group %in% OCR_2isolates_plus$group)
str(OCR_2isolates_plus_data)
# only keep strains that are present in both the north and south
# counting the number of distinct regions for the strains, keep only those with both regions (n_regions=2)
strains_in_both_regions <- OCR_2isolates_plus_data %>% 
  group_by(group) %>% 
  summarise(n_regions=n_distinct(region),.groups = 'drop') %>% 
  as.data.frame() %>% 
  filter(n_regions==2)
OCR_2isolates_plus_data <-dplyr::filter(OCR_2isolates_plus_data, group %in% strains_in_both_regions$group)

# only retain the first instance of a race (aka group) + region combo
NvS_strain_appearance <-OCR_2isolates_plus_data %>% 
  group_by(group, region) %>% 
  filter(year==min(year)) %>% 
  slice(1) %>% #only taking the first observation per region
  select(year, region, group) %>% 
  # spread into wide format
  ungroup() %>% 
  pivot_wider(names_from=region, values_from=year) %>% 
  mutate(NvS_quant=North-South) %>% 
  mutate(NvS_qual= case_when(NvS_quant>0 ~ "S", NvS_quant<0 ~ "N", NvS_quant==0 ~ "same")) %>% 
  as.data.frame() #makes it easier to see the data

NvS_strain_appearance

NvS_strain_appearance %>% 
count(NvS_qual)

#visualize with a bar chart (sqrt scale) with numbers
ggplot(strain_counts_all, aes(x=n, y=frequency)) + 
  geom_bar(stat='identity', fill="light blue")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text(aes(label=frequency), vjust=-0.3, size=2.5)+
  # scale_y_sqrt()+
  coord_trans(y= 'sqrt')+
  scale_y_continuous(limits=c(0,2600), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.5,47), expand=c(0,0))+
  xlab("Number of Isolates per Race")+
  ylab("Number of Races \n(Square root axis transformation)")+
  ggtitle("Oat Crown Rust Races")

ggsave("Figure 2 Histogram races.tiff", scale=1, dpi=600, width=178, height=100, units="mm")

######Now comparing the races of the northern and southern isolates##########
####first looking at the northern strains
OCR_N_strains_allminusNAs <- OCR_wout_2013_clean %>%
  filter(region=="North") %>%
  select(-c(Belle, Stainless,Leggett,HiFi,IAB605Xsel.,WIX4361.9,TAM.O.405,Pc91,Pc94,Pc96)) %>% 
  na.omit()
setDT(OCR_N_strains_allminusNAs)[,group :=.GRP,by = .(Pc14,Pc35,Pc36,Pc38,Pc39,Pc40,Pc45,Pc46,Pc48,Pc50,Pc51,Pc52,Pc53,Pc54,Pc55,Pc56,Pc57,Pc58,Pc59,Pc60,Pc61,Pc62,Pc63,Pc64,Pc67,Pc68,Pc70,Pc71,H548, MARVELOUS)]

str(OCR_wout_2013_clean)
str(OCR_N_strains_allminusNAs)

OCR_N_n_strains_allminusNAs <-OCR_N_strains_allminusNAs %>% 
  count(group) %>% 
  arrange(-n)

strain_N_counts_all <-OCR_N_n_strains_allminusNAs %>% 
  count(n, name= "frequency") %>% 
  # rename(nn=frequency) %>% 
  mutate(totalNofstrains= n*frequency) %>%
  mutate(percentage= totalNofstrains/sum(totalNofstrains)) %>% 
  mutate(percentageofstrains= frequency/sum(frequency))
sum(strain_N_counts_all$frequency)
sum(strain_N_counts_all$totalNofstrains)

ggplot(strain_N_counts_all, aes(x=n, y=frequency)) + 
  geom_bar(stat='identity', fill="light blue")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text(aes(label=frequency), vjust=-0.3, size=3)+
  # scale_y_sqrt()+
  coord_trans(y= 'sqrt')+
  scale_y_continuous(limits=c(0,2600), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.5,47), expand=c(0,0))+
  xlab("Number of Isolates per Race")+
  ylab("Number of Races \n(Square root axis transformation)")+
  ggtitle("Northern Region Oat Crown Rust Races")

####Now the south

OCR_S_strains_allminusNAs <- OCR_wout_2013_clean %>%
  filter(region=="South") %>%
  select(-c(Belle, Stainless,Leggett,HiFi,IAB605Xsel.,WIX4361.9,TAM.O.405,Pc91,Pc94,Pc96)) %>% 
  na.omit()
setDT(OCR_S_strains_allminusNAs)[,group :=.GRP,by = .(Pc14,Pc35,Pc36,Pc38,Pc39,Pc40,Pc45,Pc46,Pc48,Pc50,Pc51,Pc52,Pc53,Pc54,Pc55,Pc56,Pc57,Pc58,Pc59,Pc60,Pc61,Pc62,Pc63,Pc64,Pc67,Pc68,Pc70,Pc71,H548, MARVELOUS)]

str(OCR_wout_2013_clean)
str(OCR_S_strains_allminusNAs)

OCR_S_n_strains_allminusNAs <-OCR_S_strains_allminusNAs %>% 
  count(group) %>% 
  arrange(-n)

strain_S_counts_all <-OCR_S_n_strains_allminusNAs %>% 
  count(n, name= "frequency") %>% 
  # rename(nn=frequency) %>% 
  mutate(totalNofstrains= n*frequency) %>%
  mutate(percentage= totalNofstrains/sum(totalNofstrains)) %>% 
  mutate(percentageofstrains= frequency/sum(frequency))
sum(strain_S_counts_all$frequency)
sum(strain_S_counts_all$totalNofstrains)

ggplot(strain_S_counts_all, aes(x=n, y=frequency)) + 
  geom_bar(stat='identity', fill="light blue")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text(aes(label=frequency), vjust=-0.3, size=3)+
  # scale_y_sqrt()+
  coord_trans(y= 'sqrt')+
  scale_y_continuous(limits=c(0,2600), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.5,47), expand=c(0,0))+
  xlab("Number of Isolates per Race")+
  ylab("Number of Races \n(Square root axis transformation)")+
  ggtitle("Southern Region Oat Crown Rust Races")


####Presenting the data in one bar graph

# Combine the two datasets
strain_N_counts_all <- strain_N_counts_all %>%
  mutate(Region="North")

strain_S_counts_all <- strain_S_counts_all %>%
  mutate(Region="South")

strain_regional_counts_all <- bind_rows(strain_N_counts_all, strain_S_counts_all)
str(strain_regional_counts_all)

#graphing the percentage of races in order to better compare the two, because North has more observations
ggplot(strain_regional_counts_all, aes(x=n, y=percentageofstrains, fill=Region)) + 
  geom_bar(stat='identity', position=position_dodge())+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text(aes(label=round(percentageofstrains*100, digits=1)), position = position_dodge(width = 1), vjust=-0.3, size=1.65)+
  # scale_y_sqrt()+
  coord_trans(y= 'sqrt')+
  # scale_y_continuous(limits=c(0,2600), expand=c(0,0))+
  # scale_x_continuous(limits=c(-0.5,47), expand=c(0,0))+
  xlab("Number of Isolates per Race")+
  ylab("Percentage of Races \n(Square root axis transformation)")+
  # scale_fill_brewer(palette="Paired")+
  ggtitle("Regional Oat Crown Rust Races")

ggsave("Sup Figure S1 Histogram races N vs S.tiff", scale=0.9, dpi=300, width=2000, height=1500,units="px")


#create boxplot showing countvirulence by year with the alpha (transparency) corresponding to the relative amount of samples

# counting up the number of observations per year and dividing by the largest to get values between 0 and 1
samples_per_year <- OCR_wout_2013_clean %>% count(year, sort = FALSE) %>% mutate(n_adjusted=n/max(n))

typeof(samples_per_year)
str(samples_per_year)


boxplot1 <-ggplot(OCR_wout_2013_clean, aes(x=factor(year, levels=1993:2022), y=countvirulence)) + 
  geom_boxplot(fill="dark orange", alpha=samples_per_year$n_adjusted, outlier.alpha = 1, outlier.size=0.5, notch=TRUE)+
  xlab("Year")+
  ylab("Number of virulences per rust isolate")+
  ggtitle("Virulences per rust isolate")+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), breaks=factor(1993:2022), drop=FALSE)+
  scale_y_continuous(limits = c(0,40), expand = c(0, 0))+
  # geom_jitter(color="black", size=0.2, alpha=0.3) + #not as much of a fan of jitter for this plot
  annotate("rect", xmin = 0.4, xmax = 12.5, ymin = 30, ymax = 40,fill = "grey")+ #adding grey rectangles because didn't have all 40 lines in the beginning
  annotate("rect", xmin = 12.5, xmax = 17.5, ymin = 36, ymax = 40,fill = "grey")+
  annotate("point", x = 21, y = 16, size=0.5) + #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot
  annotate("point", x = 21, y = 24, size=0.5) + #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot 
  theme(text=element_text(size=9))
  
boxplot1

#######manually print out a legend color scale for the side because unable to print out a legend using alpha

library(ggpubr)

boxplot_legend <-ggplot(samples_per_year, aes(x=factor(year, levels=1993:2022), y= n_adjusted, fill=n)) + 
  geom_point()+
  scale_fill_gradient(low=rgb(1,0.549,0,0), high=rgb(1,0.549,0,1), limit = c(0,544), name="Number of \nIsolates",
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+ #manually specifying dark orange with alpha 1 and alpha 0
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8.5))

boxplot_legend
boxplot_legend <- get_legend(boxplot_legend) #extracting the legend
boxplot_legend <- as_ggplot(boxplot_legend) #making the legend a ggplot object
boxplot_legend

boxplotwlegend <- ggarrange(boxplot1, boxplot_legend, widths=c(1,0.15),
                          nrow=1)
boxplotwlegend
ggsave("Figure 3 boxplotwlegend.tiff", scale=1, dpi=600, width=178, height=125, units="mm", path=".", bg="white")

#############################

# Regression analysis for the mean isolate virulence per year (using the mean for only the original 30 lines)
OCR_vir_regression<-OCR_wout_2013_clean
OCR_vir_regression$countvirulenceOG30 = rowSums(OCR_wout_2013_clean[,c(8:35,39,40)], na.rm = TRUE)
tail(OCR_vir_regression)

#take the mean for each year
# counting up the number of observations per year and dividing by the largest to get values between 0 and 1
meanvir <-aggregate(OCR_vir_regression$countvirulenceOG30, list(OCR_vir_regression$year), FUN=mean)
# Group.1=year
# x=meanvirulence

lmVir1 =lm(x~Group.1, data=meanvir)
summary(lmVir1)


###################################################################
#building a heatmap of virulences per differential line over time #
###################################################################

library("colorspace")

str(OCR_wout_2013_clean)
heatmap_data <- OCR_wout_2013_clean %>% group_by(year) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100))

### change data into long format
heatmap_data_long <-gather(heatmap_data, differential_line, percent_virulent, Pc14:Stainless)
str(heatmap_data_long)

#createheatmap
heatmap_year <-ggplot(heatmap_data_long, aes(factor(year), factor(differential_line, levels=rev(unique(differential_line))), fill=percent_virulent))+
  geom_tile(show.legend = FALSE)+
  scale_fill_continuous_sequential(palette = "YlOrRd", rev = TRUE)+
  theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_text(aes(label = round(percent_virulent, 1)), color = "black", size = 3)+
  ylab("Differential Line")+
  ggtitle("Percentage of virulent isolates per year")
heatmap_year

#Order the heatmap rows by hierarchical clustering so that patterns are easier to see

heatmap_data_matrix <- as.matrix(heatmap_data[, -1]) # -1 to omit categories from matrix
# Cluster based on euclidean distance
clust <- hclust(dist(t(heatmap_data_matrix)))
rust_clust <- as.dendrogram(hclust(dist(t(heatmap_data_matrix))))

library("ggdendro")

# Create dendrogram
dendro_data_for_plot <- dendro_data(rust_clust, type = "rectangle")
dendro_rust <- ggplot(segment(dendro_data_for_plot)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() +
  theme_dendro()+
  theme(plot.margin = unit(c(5.5,0,5.5,5.5), "points"), text=element_text(size=7.5))+
  scale_y_reverse(expand = c(0, 0))+
  scale_x_continuous(expand=c(0.015, 0.015))+
  annotate("text", x=29.0, y=170, label= "A", color="red", cex=3)+ #manually designate the groups of isolates that I discuss in the article
  annotate("text", x=13.3, y=230, label= "B", color="red", cex=3)+
  annotate("text", x=4.7, y=230, label= "C", color="red", cex=3)+
  annotate("text", x=24.2, y=169, label= "D", color="red", cex=3)

dendro_rust

#ordered heatmap using dendrogram
heatmap_year_ordered <-ggplot(heatmap_data_long, aes(factor(year), factor(differential_line, levels=rev(unique(differential_line))), fill=percent_virulent))+
  geom_tile(show.legend = FALSE)+
  scale_y_discrete(limits = colnames(heatmap_data_matrix)[clust$order])+
  scale_fill_continuous_sequential(palette = "YlOrRd", rev = TRUE)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size=11), plot.margin = unit(c(5.5,5.5,5.5,0), "points"), text=element_text(size=7.5))+ #decrease plot margin on left side so lines up closer to the dendrogram
  geom_text(aes(label = round(percent_virulent, 1)), color = "black", size = 2)+
  # ylab("Differential Line")+
  ggtitle("Percentage of virulent isolates per year \n")+
  annotate("segment", x = 20.5, xend = 20.5, y = 40.5, yend = 0.5, colour = "gray", linewidth=1, alpha=0.4)
heatmap_year_ordered

# create summary figure grouped by decade
# 1993-2002, 2003-2012, 2013-2022
# Use the clean data with the two observations from 2013
OCR_survey %>%
  mutate(across("Pc91":"Pc96", ~ifelse(isolate=="96WI064", NA, .))) %>% 
  mutate(across("IAB605Xsel.":"TAM.O.405", ~ifelse(isolate=="96WI064", NA, .))) %>% 
  mutate(across("Pc91":"Pc96", ~ifelse(isolate=="01TX010", NA, .))) %>% 
  mutate(across("IAB605Xsel.":"TAM.O.405", ~ifelse(isolate=="01TX010", NA, .))) -> OCR_survey_clean

# now subset data into three decades and take the mean
OCR_survey_clean %>%
  filter(year<=2002) %>% 
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_1993_to_2002
OCR_survey_clean %>%
  filter(year>2002, year<2013) %>% 
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_2003_to_2012
OCR_survey_clean %>%
  filter(year>=2013) %>% 
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_2013_to_2022
# bind the means together
Mean_by_decade <- bind_rows("1993-2002"= Mean_1993_to_2002, "2003-2012"= Mean_2003_to_2012, "2013-2022"= Mean_2013_to_2022, .id="Decade")

### change data into long format
Mean_by_decade_long <-gather(Mean_by_decade, differential_line, percent_virulent, Pc14:Stainless)
str(Mean_by_decade_long)

#ordered heatmap by decade
heatmap_decade_ordered <- ggplot(Mean_by_decade_long, aes(factor(Decade), factor(differential_line, levels=rev(unique(differential_line))), fill=percent_virulent))+
  geom_tile(show.legend = FALSE)+
  scale_y_discrete(limits = colnames(heatmap_data_matrix)[clust$order])+
  # scale_fill_gradient2(low = "light blue", mid = "white", high = "red", midpoint = 0.5)+
  # scale_fill_gradient2(low = "white", high = "red")+
  # scale_fill_distiller(palette = "YlGnBu")+
  # scale_fill_viridis_c(option="I")+
  scale_fill_continuous_sequential(palette = "YlOrRd", rev = TRUE)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5), text=element_text(size=7.5))+
  scale_x_discrete(labels=c("1993-2002" = "1993-\n2002", "2003-2012" = "2003-\n2012",
                            "2013-2022" = "2013-\n2022"))+
  geom_text(aes(label = round(percent_virulent, 1)), color = "black", size = 2)+
  ggtitle("Percentage of \nvirulent isolates \nper decade")
heatmap_decade_ordered

#arrange the two ordered plots together with the dendrogram
arrangedplotsordered <- ggarrange(dendro_rust, heatmap_year_ordered, heatmap_decade_ordered,
                                  labels=c("a", "","b"),
                                  nrow=1,
                                  align="h",
                                  widths=c(0.08, 1, 0.24))
arrangedplotsordered

annotate_figure(arrangedplotsordered,
                left = text_grob("Differential Line", rot = 90))


####### Linear regression for each pc gene over time ###########
# https://www.datacamp.com/tutorial/linear-regression-R
# lm([target] ~ [predictor / features], data = [data source])

# doing all the lines together actually ignores the years that don't have observations for all of the lines (only takes last 10 years)
# do separately for the three sets of lines depending on how long they have been used

lmOG =lm(cbind(Pc14,Pc35,Pc36,Pc38,Pc39,Pc40,Pc45,Pc46,Pc48,Pc50,Pc51,Pc52,Pc53,Pc54,Pc55,Pc56,Pc57,Pc58,Pc59,Pc60,Pc61,Pc62,Pc63,Pc64,Pc67,Pc68,Pc70,Pc71,H548, MARVELOUS)~year, data=heatmap_data)
summary(lmOG)

lm2005 =lm(cbind(Pc91,Pc94,Pc96,IAB605Xsel.,WIX4361.9,TAM.O.405)~year, data=heatmap_data)
summary(lm2005)

lm2010 =lm(cbind(Belle,HiFi,Leggett,Stainless)~year, data=heatmap_data)
summary(lm2010)

# Was not able to efficiently copy out the slope, p value, and r2 with code so did manually, results are in lmdata.csv

lmdata <- read.csv("lmdata.csv")
lmdata <- lmdata %>% 
  mutate(slope= case_when(Pvalue <= 0.05 ~ year)) %>%
  mutate(signif= case_when(Pvalue <=0.05 & Pvalue > 0.01 ~ "*",
                           Pvalue <=0.01 & Pvalue > 0.001 ~ "**",
                           Pvalue <=0.001 ~ "***")) %>%
  mutate(slope_sig=paste(as.character(round(slope,1)),signif, sep=" "))%>% #combining the slope and significance into one string
  mutate(slope_sig = na_if(slope_sig,"NA NA")) %>% #turning "Na Na" into a real <NA>
  mutate(y=1) #adding a y=1 for plotting in the heatmap

### creating an ordered heatmap to add to the final heatmap to communicate the slope and significance
heatmap_slope <-ggplot(lmdata, aes(factor(y), factor(gene, levels=rev(unique(gene))), fill=year))+
  geom_tile(show.legend = FALSE)+
  scale_y_discrete(limits = colnames(heatmap_data_matrix)[clust$order])+
  # scale_fill_gradient2(low = "light blue", mid = "white", high = "red", midpoint = 0.5)+
  # scale_fill_gradient2(low = "white", high = "red")+
  # scale_fill_distiller(palette = "YlGnBu")+
  # scale_fill_viridis_c(option="I")+
  scale_fill_continuous_diverging(palette = "Blue-Red 3", rev = FALSE, na.value="white")+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        plot.title = element_blank(), 
        plot.margin = unit(c(5.5,5.5,5.5,0), "points"), #decrease plot margin on left side so lines up closer to the dendrogram
        text=element_text(size=7.5))+ 
  scale_x_discrete(labels =c("1" = "change\nper year"))+
  geom_text(aes(label = slope_sig), color = "black", size = 2)+
  # ylab("Differential Line")+
  ggtitle("Slope")
heatmap_slope

#arrange all the plots together
arrangedplotsordered_slope <- ggarrange(dendro_rust, heatmap_year_ordered, heatmap_slope, heatmap_decade_ordered,
                                        labels=c("a", "","","b"),
                                        nrow=1,
                                        align="h",
                                        widths=c(0.08, 1, 0.06, 0.22))
arrangedplotsordered_slope

annotate_figure(arrangedplotsordered_slope,
                left = text_grob("Differential Line", rot = 90, size=8))


ggsave("Figure 4 heatmap 450dpi.tiff", scale=1, dpi=450, width=246, height=178, units="mm", path=".", bg="white")

#################################################
#Plotting the heat map of observations per state#
#################################################

library("rnaturalearth")
library("sf")

states_us <- ne_states(country="United States of America",returnclass = 'sf')
class(states_us)

#working from the original dataframe, need to rename state column to postal so that it merges correctly with the sf object 
str(OCR_survey)
samples_per_state <- OCR_survey %>% count(state, sort = FALSE, name="rust_count") %>% rename(postal="state")
#add region to the samples_per_state by matching columns from OCR survey
# samples_per_state <-left_join(samples_per_state, OCR_survey, by= "postal")
# merge(samples_per_state, OCR_survey[, c("state", "region")], by="state")
# %>% rename(postal="state")

#merge the sf and the samples_per_state data frame
states_us <- merge(states_us, samples_per_state, all = TRUE)

#read in a simple CSV with the regions N vs S defined and merge with states_us
regions_key <- read.csv("Postal codes and regions.csv")
states_us <-merge(states_us, regions_key, by="postal")

#statistics on the percentages for each state for the results section
sum(samples_per_state$rust_count)
samples_per_state$perc_per_state =(samples_per_state$rust_count/sum(samples_per_state$rust_count)*100)
arrange(samples_per_state, rust_count)
# Stats on North vs South 
str(OCR_survey)
nrow(OCR_survey[OCR_survey$region == 'South', ])
nrow(OCR_survey[OCR_survey$region == 'North', ])

(mainland <- ggplot(data = states_us) +
    geom_sf(data=states_us, aes(fill = rust_count)) +
    scale_fill_gradient(low="light blue", high="dark blue", trans= "sqrt", name="Sample count")+
    geom_sf_label(aes(label = rust_count), size = 2, fontface = "bold", label.padding = unit(0.15, "lines")) +
    geom_sf(fill = "transparent", color = "black", linewidth=0.3, #adding a thick border around the two regions
            data = . %>% group_by(rust_region) %>% summarise()) +
    coord_sf(crs = st_crs(2163), xlim = c(-2500000, 2500000), ylim = c(-2300000, 730000)))

(alaska <- ggplot(data = states_us) +
    geom_sf(data=states_us, aes(fill = rust_count)) +
    scale_fill_gradient(low="light blue", high="dark blue", trans= "sqrt")+
    theme(legend.position = "none")+
    geom_sf_text(aes(label = rust_count), size = 2, fontface = "bold") +
    geom_sf_label(aes(label = rust_count), size = 2, fontface = "bold", label.padding = unit(0.15, "lines")) +
    coord_sf(crs = st_crs(3467), xlim = c(-2400000, 1600000), ylim = c(200000, 2500000), expand = FALSE, datum = NA))+
  xlab("") + ylab("")

(hawaii  <- ggplot(data = states_us) +
    geom_sf(data=states_us, aes(fill = rust_count)) +
    scale_fill_gradient(low="light blue", high="dark blue", trans= "sqrt")+
    theme(legend.position = "none")+
    coord_sf(crs = st_crs(4135), xlim = c(-161, -154), ylim = c(18, 
                                                                23), expand = FALSE, datum = NA))

mainland + xlab("Longitude") + ylab("Latitude") +
  ggtitle("1993-2022 Rust survey samples by state")+
  annotation_custom(
    grob = ggplotGrob(alaska + theme(axis.title.x=element_blank(), axis.title.y=element_blank())),
    xmin = -2750000,
    xmax = -2750000 + (1600000 - (-2400000))/2.5,
    ymin = -2450000,
    ymax = -2450000 + (2500000 - 200000)/2.5
  ) +
  annotation_custom(
    grob = ggplotGrob(hawaii),
    xmin = -1250000,
    xmax = -1250000 + (-154 - (-161))*120000,
    ymin = -2450000,
    ymax = -2450000 + (23 - 18)*120000
  )

ggsave("Figure 1 Rust sample map.tiff", scale=1, dpi=600, width=178, units="mm")

##############################################################
#Comparing the southern to northern accessions with a heatmap#
##############################################################

str(OCR_wout_2013_clean)
heatmap_data_byregion <- OCR_wout_2013_clean %>% group_by(region, year) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100))
str(heatmap_data_byregion)
head(heatmap_data_byregion, n=20L)
tail(heatmap_data_byregion, n=20L)

# Taking a closer look at the differences in numbers between northern and southern isolates in the same year to see if I can compare percentages
samples_per_year_perRegion <- OCR_wout_2013_clean %>% count(year, region, sort = FALSE)
samples_per_year_perRegion


ggplot(samples_per_year_perRegion, aes(x=year, y=n, fill=region)) +
  geom_bar(stat='identity', position='dodge')

heatmap_data_byregion
#arrange by year
heatmap_data_byregion <- arrange(heatmap_data_byregion, year)


# subtract the north by the south for each year by creating northern and southern data tibbles and extracting them from one another
heatmap_data_north <- heatmap_data_byregion %>%
  filter(region=="North")

heatmap_data_south <- heatmap_data_byregion %>%
  filter(region=="South")

#subtract north from south
heatmap_data_diff <-(heatmap_data_north[,c(3:42)])-(heatmap_data_south[,c(3:42)])

#re-add the year to the data set for visualization
year_vec <- pull(heatmap_data_north, year)
heatmap_data_diff<- heatmap_data_diff %>%
  mutate(year=year_vec,
         .before=Pc14)

### change data into long format
heatmap_data_diff_long <-gather(heatmap_data_diff, differential_line, percent_virulent, Pc14:Stainless)

#heatmap
heatmap_NvsS <- ggplot(heatmap_data_diff_long, aes(factor(year), factor(differential_line, levels=rev(unique(differential_line))), fill=percent_virulent))+
  geom_tile()+
  scale_y_discrete(limits = colnames(heatmap_data_matrix)[clust$order])+
  xlab("Year")+
  labs(fill = "North Vir Perc minus \nSouth Vir Perc   \n")+
  # ylab("Differences in virulence percentages from North to South")+
  scale_fill_continuous_diverging(palette = "Purple_Green", rev = TRUE, limits=c(-100,100))+
  theme(plot.title = element_text(hjust = 0.5, size=11), axis.title.y=element_blank(), axis.title.x=element_blank(), legend.position="bottom", text=element_text(size=8))+
  geom_text(aes(label = round(percent_virulent, 1)), color = "black", size = 2.2)+
  ggtitle("Difference in percentage of virulent isolates per year by region")+
  annotate("segment", x = 20.5, xend = 20.5, y = 40.5, yend = 0.5, colour = "gray", linewidth=1, alpha=0.6)
heatmap_NvsS


#######Decade heatmap comparing N vs S########

OCR_survey_clean %>%
  filter(region=="North") %>%
  filter(year<=2002) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_North_1993_to_2002
OCR_survey_clean %>%
  filter(region=="South") %>%
  filter(year<=2002) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_South_1993_to_2002

OCR_survey_clean %>%
  filter(region=="North") %>%
  filter(year>2002, year<2013) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_North_2003_to_2012
OCR_survey_clean %>%
  filter(region=="South") %>%
  filter(year>2002, year<2013) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_South_2003_to_2012

OCR_survey_clean %>%
  filter(region=="North") %>%
  filter(year>=2013) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_North_2013_to_2022
OCR_survey_clean %>%
  filter(region=="South") %>%
  filter(year>=2013) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_South_2013_to_2022

# bind the means together
South_decade_mean <- bind_rows("1993-2002"= Mean_South_1993_to_2002, "2003-2012"= Mean_South_2003_to_2012, "2013-2022"= Mean_South_2013_to_2022, .id="Decade")
North_decade_mean <- bind_rows("1993-2002"= Mean_North_1993_to_2002, "2003-2012"= Mean_North_2003_to_2012, "2013-2022"= Mean_North_2013_to_2022, .id="Decade")
str(North_decade_mean)

# subtract mean of south from mean of north
heatmap_data_mean_diff <- North_decade_mean[,c(2:41)] - South_decade_mean[,c(2:41)]
str(heatmap_data_mean_diff)

# Re-add the decade names to dataset
decade_vec <- pull(South_decade_mean, Decade)
heatmap_data_mean_diff<- heatmap_data_mean_diff %>%
  mutate(decade=decade_vec,
         .before=Pc14)
str(heatmap_data_mean_diff)

#put into long form
heatmap_data_mean_long <-gather(heatmap_data_mean_diff, differential_line, percent_virulent, Pc14:Stainless)
str(heatmap_data_mean_long)

#ordered heatmap
heatmap_mean_diff <- ggplot(heatmap_data_mean_long, aes(factor(decade), factor(differential_line, levels=rev(unique(differential_line))), fill=percent_virulent))+
  geom_tile(show.legend = FALSE)+
  scale_y_discrete(limits = colnames(heatmap_data_matrix)[clust$order])+
  # labs(fill = "North Vir Perc\n      minus \nSouth Vir Perc")+
  scale_fill_continuous_diverging(palette = "Purple_Green", rev = TRUE, limits=c(-100,100))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5), text=element_text(size=8))+
  scale_x_discrete(labels=c("1993-2002" = "1993-\n2002", "2003-2012" = "2003-\n2012",
                            "2013-2022" = "2013-\n2022"))+
  geom_text(aes(label = round(percent_virulent, 1)), color = "black", size = 2.2)+
  ggtitle("Differences \nby decade")
heatmap_mean_diff

# Manually creating dendrogram groups as a graphic
line_coord <- data.frame(x1 = c(0,0,0,0,0,0,0,0,1,1,1,1), x2 = c(1,1,1,1,1,1,1,1,1,1,1,1), 
                         y1 = c(1,8,9,20,22,25,26,37,1,9,22,26), y2 = c(1,8,9,20,22,25,26,37,8,20,25,37))

dendro_rust_groups <- ggplot(line_coord, aes(x = x1, y = y1, xend = x2, yend = y2)) + 
  geom_segment() +
  # coord_flip() +
  theme_dendro()+
  theme(plot.margin = unit(c(5.5,0,5.5,5.5), "points"))+
  scale_y_continuous(expand = c(0, 0), limits=c(0.5,40.5))+ #forsome reason have to manually adjust the y axis limits so that it lines up properly on the combined graphic 
  scale_x_reverse(expand=c(0.015, 0.015), limits=c(3.5,0))+
  annotate("text", x=2.5, y=(37-26)/2+26, label= "A", color="red", cex=3)+
  annotate("text", x=2.5, y=(25-22)/2+22, label= "D", color="red", cex=3)+
  annotate("text", x=2.5, y=(20-9)/2+9, label= "B", color="red", cex=3)+
  annotate("text", x=2.5, y=(8-1)/2+1, label= "C", color="red", cex=3)

dendro_rust_groups

# Arranging all three graphs together with the dendrogram groups clearly marked
arrangeddiffplotswgroups <- ggarrange(dendro_rust_groups, heatmap_NvsS, heatmap_mean_diff,
                               labels=c("", "a", "b"),
                               nrow=1,
                               align="h",
                               widths=c(0.035,1,0.22),
                               common.legend = TRUE,
                               legend="bottom")
arrangeddiffplotswgroups <- annotate_figure(arrangeddiffplotswgroups,
                left = text_grob("Differences in virulence percentages from North to South", rot = 90, size=10))

ggsave("Figure 6 N vs S heatmap 450dpi.tiff", scale=1, dpi=450, width=246, height=178, units="mm", path=".", bg="white")


#################################################################################################
#Boxplot comparing the differences in virulences between northern and southern isolates per year#
#################################################################################################

str(OCR_wout_2013_clean)
NvS_boxplot <-ggplot(OCR_wout_2013_clean, aes(x=factor(year, levels=1993:2022), y=countvirulence, fill=region, group = interaction(year, region))) + 
  # geom_boxplot(alpha=samples_per_year$n_adjusted, outlier.alpha = 1, outlier.size=0.5)+
  geom_boxplot(outlier.size=0.5)+
  xlab("Year")+
  ylab("Number of virulences per rust isolate")+
  ggtitle("Virulences per rust isolate")+
  theme_bw()+
  theme(legend.position="bottom")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), breaks=factor(1993:2022), drop=FALSE)+
  scale_y_continuous(limits = c(0,40), expand = c(0, 0))+
  # geom_jitter(color="black", size=0.2, alpha=0.3) + #not as much of a fan of jitter for this plot
  annotate("rect", xmin = 0.4, xmax = 13.5, ymin = 30, ymax = 40,fill = "grey")+ #adding grey rectangles because didn't have all 40 lines in the beginning
  annotate("rect", xmin = 13.5, xmax = 18.5, ymin = 36, ymax = 40,fill = "grey")+
  annotate("point", x = 21, y = 16, size=0.5) + #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot
  annotate("point", x = 21, y = 24, size=0.5) #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot

NvS_boxplot

###### Wilcox test #######

library(rstatix)

OCR_wout_2013_long_ttest <- OCR_wout_2013_clean %>%
  select(isolate, year, region, countvirulence)

wilcox.stat.test <- OCR_wout_2013_long_ttest %>%
  group_by(year) %>%
  wilcox_test(countvirulence ~ region) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
wilcox.stat.test %>% print(n=30)

wilcox.stat.test <- wilcox.stat.test %>%
  add_xy_position(x = "year")%>%
  select(-x) %>% #remove x because improperly calculated and will not show up in the right place on the graph
  add_column(x=1, .before = 'xmin') %>%
  mutate(x= year-1992)%>%
  mutate(xmin= x-0.2) %>%
  mutate(xmax= x+0.2)
wilcox.stat.test

# Add the significances to the graph
str(OCR_wout_2013_clean)
ggplot(OCR_wout_2013_clean, aes(x=factor(year, levels=1993:2022), y=countvirulence, fill=region, group = interaction(year, region))) +
  geom_boxplot(outlier.size=0.5)+
  xlab("Year")+
  ylab("Number of virulences per rust isolate")+
  ggtitle("Virulences per rust isolate")+
  theme_bw()+
  theme(legend.position="bottom")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), breaks=factor(1993:2022), drop=FALSE)+
  scale_y_continuous(limits = c(0,40), expand = c(0, 0))+
  # geom_jitter(color="black", size=0.2, alpha=0.3) + #not as much of a fan of jitter for this plot
  annotate("rect", xmin = 0.4, xmax = 13.5, ymin = 30, ymax = 40,fill = "grey")+ #adding grey rectangles because didn't have all 40 lines in the beginning
  annotate("rect", xmin = 13.5, xmax = 18.5, ymin = 36, ymax = 40,fill = "grey")+
  annotate("point", x = 21, y = 16, size=0.5) + #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot
  annotate("point", x = 21, y = 24, size=0.5) + #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot
  stat_pvalue_manual(wilcox.stat.test, label = "p.adj.signif", inherit.aes = FALSE, bracket.nudge.y= -0.5, hide.ns=TRUE)

ggsave("Sup Figure S2 NvsS.tiff", scale=1, dpi=300, width=2500, height=1800, units="px", path=".")

# #now make the same graph but for combining with barplot below
NvS_boxplot <-ggplot(OCR_wout_2013_clean, aes(x=factor(year, levels=1993:2022), y=countvirulence, fill=region, group = interaction(year, region))) +
  # geom_boxplot(alpha=samples_per_year$n_adjusted, outlier.alpha = 1, outlier.size=0.5)+
  geom_boxplot(outlier.size=0.5)+
  ylab("Number of virulences per rust isolate")+
  ggtitle("Virulences per rust isolate")+
  theme_bw()+
  theme(legend.position="none", axis.title.x=element_blank())+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), breaks=factor(1993:2022), drop=FALSE)+
  scale_y_continuous(limits = c(0,40), expand = c(0, 0))+
  # geom_jitter(color="black", size=0.2, alpha=0.3) + #not as much of a fan of jitter for this plot
  annotate("rect", xmin = 0.4, xmax = 13.5, ymin = 30, ymax = 40,fill = "grey")+ #adding grey rectangles because didn't have all 40 lines in the beginning
  annotate("rect", xmin = 13.5, xmax = 18.5, ymin = 36, ymax = 40,fill = "grey")+
  annotate("point", x = 21, y = 16, size=0.5) + #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot
  annotate("point", x = 21, y = 24, size=0.5) + #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot
  stat_pvalue_manual(wilcox.stat.test, label = "p.adj.signif", inherit.aes = FALSE, bracket.nudge.y= -0.5, hide.ns=TRUE)
NvS_boxplot

# Graphing the data for north and south in order to compare the sample sizes by year
samples_per_year_perRegion <- OCR_survey_clean %>% count(year, region, sort = FALSE)
samples_per_year_perRegion


Sample_Number_by_region_bargraph <-ggplot(samples_per_year_perRegion, aes(x=factor(year, levels=1993:2022), y=n, fill=region)) +
  geom_bar(stat='identity', position='dodge', width= 0.65)+
  theme_bw()+
  theme(legend.position="bottom")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), breaks=factor(1993:2022), drop=FALSE)+
  xlab("Year")+
  ylab("Number of isolates by region")
Sample_Number_by_region_bargraph


arrangedNvsS <- ggarrange(NvS_boxplot, Sample_Number_by_region_bargraph,
                          labels=c("a", "b"),
                          nrow=2,
                          align="v")
arrangedNvsS
ggsave("NvsSarranged.tiff", scale=1, dpi=300, width=4000, height=2300, units="px", path=".")


#####################################
#Performing polychorric correlations#
#####################################

library("polycor")

# perform the analysis without Marvelous because it is causing singularity errors (due to it always being susceptible)

OCR_for_polychor <-OCR_wout_2013_clean[,8:47]
OCR_for_polychorNoMarv <-OCR_for_polychor[,-32]

#change from 0 and 1 into characters so that it doesn't mistake them for numbers and use the pearson's correlation on it
OCR_for_polychorNoMarv[OCR_for_polychorNoMarv==0] <- "Av"
OCR_for_polychorNoMarv[OCR_for_polychorNoMarv==1] <- "V"

OCR_correlation <- hetcor(OCR_for_polychorNoMarv, ML = TRUE, std.err = TRUE, use="pairwise.complete.obs", 
                          bins=4, pd=TRUE, parallel=FALSE, ncores=detectCores(logical=FALSE), thresholds=FALSE)

tail(OCR_correlation$std.errors)

#checking values of individual correlations
OCR_correlation$correlations["Pc38" ,"Pc63"]
OCR_correlation$correlations["Pc68" ,"HiFi"]
OCR_correlation$correlations["Pc48" ,"Pc52"]
OCR_correlation$correlations["Pc62" ,"Pc64"]

OCR_correlation$correlations["Pc91" ,"HiFi"]
OCR_correlation$correlations["Pc35" ,"Pc58"]

OCR_correlation$correlations
class(OCR_correlation$correlations)

OCR_correlation_matrix <- OCR_correlation$correlations
str(OCR_correlation_matrix)
class(OCR_correlation_matrix)
# reorder the correlation matrix before further melting
# using the rev() to reverse the order of the names because corrplot uses the opposite order of ggplot heatmaps
order <-colnames(heatmap_data_matrix)[rev(clust$order)]
order <- order[-20] #removing Marvelous from the name order vector
OCR_corr_ordered <-OCR_correlation_matrix[,order]
OCR_corr_ordered2 <-OCR_corr_ordered[order, ]
is.vector(clust$order)

library("pgirmess")
get_upper_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

OCR_lower_tri <- get_upper_tri(OCR_corr_ordered2)

library(corrplot)

tiff(file = "Figure 5 correlations.tiff", width = 178, height = 178, units = "mm", res = 600)
corrplot(OCR_lower_tri, type='lower', tl.col = 'black', na.label='X', na.label.col='gray', tl.cex=0.75)
dev.off()


#######################
#additional statistics#
#######################

# calculate the virulence count of the original 30 differentials in 2022
meanvirulence2022 <-OCR_survey %>%
  filter(year==2022)
mean(rowSums(meanvirulence2022[,c(8:35,39,40)], na.rm = TRUE))

meanvirulence1993 <-OCR_survey %>%
  filter(year==1993)
mean(rowSums(meanvirulence1993[,c(8:35,39,40)], na.rm = TRUE))


#plotting all the points is messy but here is the code
dotplot_long <-gather(heatmap_data, differential_line, percent_virulent, Pc14:Stainless)
str(dotplot_long)
ggplot(dotplot_long, aes(x=year, y=percent_virulent, group=differential_line))+
  geom_point()+
  geom_line(aes(col=differential_line))

###checking percent virulence to individual differentials
##plotting one at a time
ggplot(heatmap_data, aes(x=year, y=Pc45))+
  geom_point()+
  ylim(0,100)

# checking geographical distribution for a particular year
OCR_wout_2013_clean %>% 
  filter(year==2018) %>% 
  group_by(state) %>% 
  summarise(total_count=n(),.groups = 'drop') %>% 
  as.data.frame()
