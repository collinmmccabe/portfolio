survey <- read.csv("survey_data.csv") # read in the dataset

# make sure that factors are being interpreted as such
survey$symptomatic <- as.factor(survey$symptomatic)
survey$inf_group <- as.factor(survey$inf_group)
survey$sex <- as.factor(survey$sex)

# subset the data for only those individuals that were symptomatic
symptomatic.Y <- survey[which(survey$symptomatic == "Y"),]