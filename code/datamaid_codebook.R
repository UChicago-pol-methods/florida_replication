library(dataMaid)

## read in data

df_analysis<-readRDS("../data/df_for_analysis_processed.rds")
df_keep<-df_analysis


## Primary outcome:  general trans tolerance attitudes
attr(df_keep$sender_id, 'shortDescription') <- "Unique identifier."
attr(df_keep$therm_trans_t1, 'shortDescription') <- "How do you feel towards transgender people? The higher the number, the warmer or more favorable you feel toward that person, the lower the number, the colder or less favorable you feel. You can pick any number between 0 and 10.
Coded from 1-10."
attr(df_keep$gender_norm_sexchange_t1, 'shortDescription') <- "I would support a friend choosing to have a sex change; 
Coded as: -1 = Disagree, 0 = No opinion/don’t know, 1 = Agree."
attr(df_keep$gender_norm_moral_t1,  'shortDescription')<-"It is morally wrong for a man to present himself as awoman in public; 
Coded as: -1 = Agree, 0 = No opinion/don’t know, 1 = Disagree; reverse coded."
attr(df_keep$gender_norm_abnormal_t1, 'shortDescription')<-"A man who identifies as a woman is psychologically abnormal;
Coded as: -1 = Agree, 0 = No opinion/don’t know, 1 = Disagree; reverse coded."
attr(df_keep$gender_norm_trans_moral_wrong_t1, 'shortDescription')<-"Saying you are a gender that is different than the one you were born with is morally wrong;
Coded as: -1 = Agree, 0 = No opinion/don’t know, 1 = Disagree; reverse coded."
attr(df_keep$trans_teacher_t1, 'shortDescription')<-"Transgender women (people who identify as women but were designated male at birth) should not be allowed to serve as public school teachers;
Coded as: -1 = Agree, 0 = No opinion/don’t know, 1 = Disagree; reverse coded."
attr(df_keep$trans_bathroom_t1, 'shortDescription')<-"It would be wrong to allow a transgender woman (a personwho identifies as a woman but was designated male at birth) to use the woman’s restroom;
Coded as: -1 = Agree, 0 = No opinion/don’t know, 1 = Disagree; reverse coded."
attr(df_keep$trans_bathroom_t1, 'shortDescription')<-"It would be wrong to allow a transgender woman (a personwho identifies as a woman but was designated male at birth) to use the woman’s restroom;
Coded as: -1 = Agree, 0 = No opinion/don’t know, 1 = Disagree; reverse coded."
attr(df_keep$gender_norm_dress_t1, 'shortDescription')<-"Men should dress like men and women should dress likewomen;
Coded as: -1 = Agree, 0 = No opinion/don’t know, 1 = Disagree; reverse coded."

## Secondary outcome:  policy only
attr(df_keep$florida_trans_policy_t1, 'shortDescription')<-"Updates to the the Florida State Medicaid policy willexclude gender-affirming care in state Medicaid coverage. Do you favor or opposethis new policy?;
Coded as: -3 = Strongly favor, -2 = Favor, -1 = Somewhat favor,0 = Neither favor nor oppose, 1 = Somewhat oppose, 2 = Oppose, 3 = Strongly oppose"
attr(df_keep$florida_trans_policy2_t1, 'shortDescription')<-"Some people say it’s important to provide gender-affirming health care to transgender people. Other people have concerns about the risks with this type of health care, and do not want gender-affirming care for transgender people included in our state Medicaid coverage. What do you think? Do you agree or disagree that Florida policy should protect transgender people from discrimination?;
Coded as: -3 = Strongly disagree, -2 = Disagree, -1 = Somewhatdisagree, 0 = Neither agree nor disagree, 1 = Somewhat agree, 2 = Agree, 3 = Stronglyagree"


## Covariates

attr(df_keep$age_t0,'shortDescription')<-"Pre-survey, age. How old are you?"
attr(df_keep$gender_t0,'shortDescription')<-"Pre-survey, gender. Do you describe yourself as a man, a woman, or in some other way? Coding: 1 = Male, 0 = otherwise."
attr(df_keep$ideology_t0,'shortDescription')<-"Pre-survey, ideology. When it comes to your political views, how would you describe yourself? Coding: -3 = Very liberal, -2 = Liberal, -1 = Somewhat liberal, 0 = Middle of the road, 1 = Somewhat conservative, 2 = Conservative, 3 = Very conservative."
attr(df_keep$pid_t0,'shortDescription')<-"Pre-survey, party ID. Generally speaking, do you consider yourself a...; Coding: -3 = Strong Democrat, -2 = Not very strongDemocrat, -1 = Closer to the Democratic Party, 0 = Not closer to either party, 1 = Closer to the Republican Party, 2 = Not very strong Republican, 3 = Strong Republican."
attr(df_keep$pol_interest_t0,'shortDescription')<-"Pre-survey, interest in politics. How interested are you in politics? Coding: -2 = Not much interested, -1 = Somewhat interested, 0 = Not sure, 1 = Very much interested."
attr(df_keep$healthcare_t0,'shortDescription')<-"Pre-survey, views on healthcare. Which comes closest to your view about providing health care in the United States? Coding: Factors: 1. The Govern-ment should provide everyone with health care and pay for it with tax dollars; 2. Companies should be required to provide health insurance for their employees and the government should provide subsidies for those who are not working or retired; 3. Health insurance should be voluntary. Individuals should either buy insurance or obtain it through their employers as they do currently. The elderly and the very poor should be covered by Medicare and Medicaid as they are currently."
attr(df_keep$climate_t0,'shortDescription')<-"Pre-survey, climate change thermomete. Should the federal government take actions to slow the effects of climate change and global warming even if it costs some jobs and makes life inconvenient for Americans or should maintaining jobs and our standard of living be given priority? Between the two positions below, on a scale of 0 to 10, where would you put yourself? You can pick any number between 0 and 10. Coding: 0-10. 0: The federal government should take actions to slow the effects of climate change 10: Maintaining jobs and standard of living should be given priority."
attr(df_keep$religion_t0,'shortDescription')<-"Pre-survey, religion. Would you say your religion provides some guidance in your day-to-day living, quite a bit of guidance, or a great deal of guidance in your day-to-day life? Coding: -1 = Some guidance, 0 = Quite a bit ofguidance, 1 = A great deal of guidance."
attr(df_keep$abortion_t0,'shortDescription')<-"Pre-survey, views on abortion. Under what circumstances should abortion be legal? Coding: Factors: 1 = Abortion shouldalways be legal. There should be no restrictions on abortion. 2 = Abortion shouldbe legal, but with some restrictions (such as for minors or late-term abortions). 3 = Abortion should only be legal in special circumstances, such as when the life of themother is in danger. 4 = Abortion should be illegal. It should never be allowed."
attr(df_keep$immigration_t0,'shortDescription')<-"Pre-survey, views on immigration. Which comes closest to your view about illegal immigration? Coding: -1 = Illegal immigrants should be arrested and deported as quickly as possible regardless of their circumstances. 1 = Illegal immigrants now living in the U.S. should be allowed to become citizens if they pay a fine and meet other requirements"



## data transformation for better visualization 
## temporarily treating 999 as NA 
df_keep[, c(2:11,13:21)][df_keep[,c(2:11,13:21)] == 999] <- NA

df_keep<-df_keep[,-c(22,23)]
## write out
makeCodebook(df_keep, vol = "", reportTitle = NULL, 
             file = '../test-data-codebook.Rmd', 
             replace = TRUE)

