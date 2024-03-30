

setwd("") ### YOUR FILE PATH HERE
library(foreign)
library(sandwich)
library(stargazer)
library(entropy)


# entropy measure
entropy_fun<-function(x){
	freqs<-x[x!=0 & !is.na(x)]
	(-1)* sum(freqs*log2(freqs))
}

# herfindahl index
herfindahl_fun<-function(x){
	freqs<-x[x!=0 & !is.na(x)]
	sum(freqs^2)
}

# function for lagging single variables:
lagger<-function(variable, country, year, laglength){
	
	country<-as.character(country)
	laggedvar<-rep(NA,length(variable))
	
	leadingNAs<-rep(NA,laglength)
	countryshift<-c(leadingNAs, country[1:(length(country)-laglength)])
	
	variableshift<-c(leadingNAs, variable[1:(length(variable)-laglength)])
	
	replacementrefs<-country==countryshift
	replacementrefs[is.na(replacementrefs)==T]<-FALSE
	laggedvar[replacementrefs]<-variableshift[replacementrefs]
	
	laggedvar
	}



eventdata<-read.dta("eventdata_icews_threedigit.dta") 

cameocodes<-c("govcit_CAMEO150", "govcit_CAMEO151", "govcit_CAMEO152", "govcit_CAMEO153", "govcit_CAMEO154", 
	"govcit_CAMEO170", "govcit_CAMEO171", "govcit_CAMEO172", "govcit_CAMEO173", "govcit_CAMEO174","govcit_CAMEO175", 
	"govcit_CAMEO180", "govcit_CAMEO181", "govcit_CAMEO182", "govcit_CAMEO183", "govcit_CAMEO184", "govcit_CAMEO185", "govcit_CAMEO186",
	"govcit_CAMEO190", "govcit_CAMEO191", "govcit_CAMEO192", "govcit_CAMEO193", "govcit_CAMEO194", "govcit_CAMEO195", "govcit_CAMEO196",
	"govcit_CAMEO200", "govcit_CAMEO201", "govcit_CAMEO202", "govcit_CAMEO203", "govcit_CAMEO204")



# summary table of event types
cameocodes_num<-gsub("govcit_CAMEO","",cameocodes)
stargazer(cbind(cameocodes_num, cameocodes_num, as.numeric(apply(eventdata[cameocodes],2, sum)) , as.numeric(round(apply(eventdata[cameocodes],2, sum)/sum(eventdata$govcit_matconf),3))),
	 align=T, type="latex", column.sep.width="0pt", font.size="footnotesize", no.space=T)
 

# creating entropy and other diversity measures

eventdata_props<-apply(eventdata[,cameocodes],2,function(x){x/eventdata$govcit_matconf})
eventdata_props<-as.data.frame(eventdata_props)
names(eventdata_props)<-paste(names(eventdata_props),"prop",sep="_")

eventdata$eventsum<-rowSums(eventdata[, cameocodes])
eventdata$entropy<-apply(eventdata_props, 1, entropy_fun)


# alternative measures
eventdata$herfindahl<-apply(eventdata_props, 1, herfindahl_fun)
eventdata$code173share<-(eventdata$govcit_CAMEO173/eventdata$eventsum)
eventdata$entropy_173vsothers<-apply(cbind(eventdata$code173share, (1-eventdata$code173share)), 1, entropy_fun)

eventdata$entropy_mm<-apply(eventdata[,cameocodes], 1, entropy.MillerMadow,unit="log2")
eventdata$entropy_mm<-ifelse(is.na(eventdata$entropy_mm),0,eventdata$entropy_mm)

eventdata_omit173<-eventdata[cameocodes]
eventdata_omit173<-eventdata_omit173[,names(eventdata_omit173)!="govcit_CAMEO173"]
eventdata_omit173_props<- apply(eventdata_omit173,2,function(x){x/ (rowSums(eventdata_omit173))})
eventdata$entropy_omit173<-apply(eventdata_omit173_props, 1, entropy_fun)


# coarse-grained entropy by 2digit not 3digit
eventdata$cameo15<- eventdata$govcit_CAMEO150 + eventdata$govcit_CAMEO151 + eventdata$govcit_CAMEO152 + eventdata$govcit_CAMEO153 + eventdata$govcit_CAMEO154
eventdata$cameo17<- eventdata$govcit_CAMEO170 + eventdata$govcit_CAMEO171 + eventdata$govcit_CAMEO172 + eventdata$govcit_CAMEO173 + eventdata$govcit_CAMEO174 + eventdata$govcit_CAMEO175
eventdata$cameo18<- eventdata$govcit_CAMEO180 + eventdata$govcit_CAMEO181 + eventdata$govcit_CAMEO182 + eventdata$govcit_CAMEO183 + eventdata$govcit_CAMEO184 + eventdata$govcit_CAMEO185 + eventdata$govcit_CAMEO186
eventdata$cameo19<- eventdata$govcit_CAMEO190 + eventdata$govcit_CAMEO191 + eventdata$govcit_CAMEO192 + eventdata$govcit_CAMEO193 + eventdata$govcit_CAMEO194 + eventdata$govcit_CAMEO195 + eventdata$govcit_CAMEO196
eventdata$cameo20<- eventdata$govcit_CAMEO200 + eventdata$govcit_CAMEO201 + eventdata$govcit_CAMEO202 + eventdata$govcit_CAMEO203 + eventdata$govcit_CAMEO204

eventdata$entropy_2digit<-apply(eventdata[,c("cameo15","cameo17","cameo18","cameo19","cameo20")]/rowSums(eventdata[,c("cameo15","cameo17","cameo18","cameo19","cameo20")]), 1, entropy_fun)


# without one-a-day filtering
eventdataNOFILTER<-read.dta("eventdata_icews_threedigitNOFILTER.dta")
eventdata_propsNOFILTER<-apply(eventdataNOFILTER[,cameocodes],2,function(x){x/eventdataNOFILTER$govcit_matconf})
eventdata_propsNOFILTER<-as.data.frame(eventdata_propsNOFILTER)
names(eventdata_propsNOFILTER)<-paste(names(eventdata_propsNOFILTER),"prop",sep="_")
eventdataNOFILTER$entropyNOFILTER<-apply(eventdata_propsNOFILTER, 1, entropy_fun)
eventdataNOFILTER$eventsumNOFILTER<-rowSums(eventdataNOFILTER[, cameocodes])


# 1995 publisher restricted data
eventdata_pub95<-read.dta("eventdata_icews_threedigit_publisher_1995.dta")
eventdata_pub95_props<-apply(eventdata_pub95[,cameocodes],2,function(x){x/eventdata_pub95$govcit_matconf})
eventdata_pub95_props<-as.data.frame(eventdata_pub95_props)
names(eventdata_pub95_props)<-paste(names(eventdata_pub95_props),"prop",sep="_")
eventdata_pub95$entropy_pub95<-apply(eventdata_pub95_props, 1, entropy_fun)
eventdata_pub95$eventsum_pub95<-rowSums(eventdata_pub95[, cameocodes])


# goldstein scores
data_gold<-read.dta("eventdata_icews_goldstein.dta") #this is the file correcting the 2002 issue



# read in country data and merge 

countrydata<-read.csv("JPR_countrydata.csv")

d<-merge(eventdata, countrydata, by=c("country","year"), all.x=T, all.y=F)
d<-merge(d, eventdataNOFILTER[,c("country","year","entropyNOFILTER","eventsumNOFILTER")],by=c("country","year"), all.x=T, all.y=F)
d<-merge(d, eventdata_pub95[,c("country","year","entropy_pub95","eventsum_pub95")],by=c("country","year"), all.x=T, all.y=F)
d<-merge(d, data_gold[,c("country","year","goldsteinrepress")],by=c("country","year"),all.x=T,all.y=F)


d$logeventcount<-log(d$eventsum+1)
d$logeventcount_lag<-lagger(d$logeventcount, d$country, d$year, 1)

d$entropy_lag<-lagger(d$entropy, d$country, d$year, 1)
d$entropy_omit173_lag<-lagger(d$entropy_omit173, d$country, d$year, 1)
d$code173share_lag<-lagger(d$code173share, d$country, d$year, 1)
d$entropy_2digit_lag<-lagger(d$entropy_2digit, d$country, d$year, 1)
d$entropyNOFILTER_lag<-lagger(d$entropyNOFILTER, d$country, d$year, 1)
d$herfindahl_lag<-lagger(d$herfindahl, d$country, d$year, 1)
d$entropy_mm_lag<-lagger(d$entropy_mm, d$country, d$year, 1)

d$entropy_pub95_lag<-lagger(d$entropy_pub95, d$country, d$year, 1)

d$HRscores_rev_lag<-lagger(d$HRscores_rev, d$country, d$year, 1)
d$PTS_State_lag<-lagger(d$PTS_State, d$country, d$year, 1)
d$CIRI_rev_lag<-lagger(d$CIRI_rev, d$country, d$year, 1)
d$vdem_HR_rev_lag<-lagger(d$vdem_HR_rev, d$country, d$year, 1)

d$goldsteinrepress<-d$goldsteinrepress*-1
d$goldsteinrepress_lag<-lagger(d$goldsteinrepress, d$country, d$year, 1)










##########

# PLOTS


# entropy examples plot

barplotmatrix<-cbind(
c(1,0,0,0,0), c(0.5,0.5,0,0,0), c(1/3,1/3,1/3,0,0), c(1/4,1/4,1/4,1/4,0),c(0.2,0.2,0.2,0.2,0.2),  
c(0.15,0.15,0.15,0.15, 0.4), c(0.075, 0.15, 0.175, 0.25, 0.35), c(0.1,0.1,0.1,0.1, 0.6), c(0.05,0.05,0.05,0.05, 0.8), c(0.05,0,0.05,0, 0.9) )

tiff(file="entropy_example_distributions.tiff",height=7,width=8, units="in",res=300)
barplot(barplotmatrix, names.arg=round(apply(barplotmatrix,2,entropy_fun),3) , xlab="Entropy", main="Entropy for Example Distributions",
	col=c("gray90","gray70","gray50","gray30","gray10"))
dev.off()

pdf(file="entropy_example_distributions.pdf",height=7,width=8)
barplot(barplotmatrix, names.arg=round(apply(barplotmatrix,2,entropy_fun),3) , xlab="Entropy", main="Entropy for Example Distributions",
	col=c("gray90","gray70","gray50","gray30","gray10"))
dev.off()



# plots of average entropy over time

d_subset<-d[d$year %in% 1996:2016,]

# only countries that are consistently observed
country_select_events5<-names(tapply(d_subset$eventsum, d_subset$country, function(x){length(x[x>=5])})[tapply(d_subset$eventsum, d_subset$country, function(x){length(x[x>=5])})==21])
country_select_events10<-names(tapply(d_subset$eventsum, d_subset$country, function(x){length(x[x>=10])})[tapply(d_subset$eventsum, d_subset$country, function(x){length(x[x>=10])})==21])

tiff(file="annualworldentropyaverages1_new.tiff",height=7,width=8, units="in",res=300)

par(mfrow=c(2,2))
plot(1996:2016, tapply(d_subset$entropy[d_subset$eventsum>=5], d_subset$year[d_subset$eventsum>=5], mean, na.rm=T),type="n",
	xlab="Year",ylab="Entropy",ylim=c(1.4,2.0), main="Annual Average Entropy:\nObservations with >=5 Events")
abline(h=seq(1.4,2,0.1), col="gray80")
points(1996:2016, tapply(d_subset$entropy[d_subset$eventsum>=5], d_subset$year[d_subset$eventsum>=5], mean, na.rm=T), pch=19)
segments(x0=1996:2016, x1=1996:2016, 
	y0=tapply(d_subset$entropy[d_subset$eventsum>=5], d_subset$year[d_subset$eventsum>=5], function(x){t.test(x)$conf.int[[1]]}),
	y1=tapply(d_subset$entropy[d_subset$eventsum>=5], d_subset$year[d_subset$eventsum>=5], function(x){t.test(x)$conf.int[[2]]}) )
abline(lm(as.numeric(tapply(d_subset$entropy[d_subset$eventsum>=5], d_subset$year[d_subset$eventsum>=5], mean, na.rm=T)) ~ c(1996:2016)), lty=2, col="gray30")

plot(1996:2016, tapply(d_subset$entropy[d_subset$eventsum>=10], d_subset$year[d_subset$eventsum>=10], mean, na.rm=T),type="n",
	xlab="Year",ylab="Entropy",ylim=c(1.4,2.0), main="Annual Average Entropy:\nObservations with >=10 Events")
points(1996:2016, tapply(d_subset$entropy[d_subset$eventsum>=10], d_subset$year[d_subset$eventsum>=10], mean, na.rm=T),pch=19)
abline(h=seq(1.4,2,0.1), col="gray80")
segments(x0=1996:2016, x1=1996:2016, 
	y0=tapply(d_subset$entropy[d_subset$eventsum>=10], d_subset$year[d_subset$eventsum>=10], function(x){t.test(x)$conf.int[[1]]}),
	y1=tapply(d_subset$entropy[d_subset$eventsum>=10], d_subset$year[d_subset$eventsum>=10], function(x){t.test(x)$conf.int[[2]]}) )
abline(lm(as.numeric(tapply(d_subset$entropy[d_subset$eventsum>=10], d_subset$year[d_subset$eventsum>=10], mean, na.rm=T)) ~ c(1996:2016)), lty=2, col="gray30")


plot(1996:2016, tapply(d_subset$entropy[d_subset$country %in% country_select_events5], d_subset$year[d_subset$country %in% country_select_events5], mean, na.rm=T),type="n",
	xlab="Year",ylab="Entropy",ylim=c(1.4,2.0), main="Annual Average Entropy:\nCountries with Consistently >=5 Events")
points(1996:2016, tapply(d_subset$entropy[d_subset$country %in% country_select_events5], d_subset$year[d_subset$country %in% country_select_events5], mean, na.rm=T),pch=19)
abline(h=seq(1.4,2,0.1), col="gray80")
segments(x0=1996:2016, x1=1996:2016, 
	y0=tapply(d_subset$entropy[d_subset$country %in% country_select_events5], d_subset$year[d_subset$country %in% country_select_events5], function(x){t.test(x)$conf.int[[1]]}),
	y1=tapply(d_subset$entropy[d_subset$country %in% country_select_events5], d_subset$year[d_subset$country %in% country_select_events5], function(x){t.test(x)$conf.int[[2]]}) )
abline(lm(as.numeric(tapply(d_subset$entropy[d_subset$country %in% country_select_events5], d_subset$year[d_subset$country %in% country_select_events5], mean, na.rm=T)) ~ c(1996:2016)), lty=2, col="gray30")


plot(1996:2016, tapply(d_subset$entropy[d_subset$country %in% country_select_events10], d_subset$year[d_subset$country %in% country_select_events10], mean, na.rm=T),type="n",
	xlab="Year",ylab="Entropy",ylim=c(1.4,2.0), main="Annual Average Entropy:\nCountries with Consistently >=10 Events")
points(1996:2016, tapply(d_subset$entropy[d_subset$country %in% country_select_events10], d_subset$year[d_subset$country %in% country_select_events10], mean, na.rm=T),pch=19)
abline(h=seq(1.4,2,0.1), col="gray80")
segments(x0=1996:2016, x1=1996:2016, 
	y0=tapply(d_subset$entropy[d_subset$country %in% country_select_events10], d_subset$year[d_subset$country %in% country_select_events10], function(x){t.test(x)$conf.int[[1]]}),
	y1=tapply(d_subset$entropy[d_subset$country %in% country_select_events10], d_subset$year[d_subset$country %in% country_select_events10], function(x){t.test(x)$conf.int[[2]]}) )
abline(lm(as.numeric(tapply(d_subset$entropy[d_subset$country %in% country_select_events10], d_subset$year[d_subset$country %in% country_select_events10], mean, na.rm=T)) ~ c(1996:2016)), lty=2, col="gray30")
dev.off()




pdf(file="annualworldentropyaverages1_new.pdf",height=7,width=8)

par(mfrow=c(2,2))
plot(1996:2016, tapply(d_subset$entropy[d_subset$eventsum>=5], d_subset$year[d_subset$eventsum>=5], mean, na.rm=T),type="n",
	xlab="Year",ylab="Entropy",ylim=c(1.4,2.0), main="Annual Average Entropy:\nObservations with >=5 Events")
abline(h=seq(1.4,2,0.1), col="gray80")
points(1996:2016, tapply(d_subset$entropy[d_subset$eventsum>=5], d_subset$year[d_subset$eventsum>=5], mean, na.rm=T), pch=19)
segments(x0=1996:2016, x1=1996:2016, 
	y0=tapply(d_subset$entropy[d_subset$eventsum>=5], d_subset$year[d_subset$eventsum>=5], function(x){t.test(x)$conf.int[[1]]}),
	y1=tapply(d_subset$entropy[d_subset$eventsum>=5], d_subset$year[d_subset$eventsum>=5], function(x){t.test(x)$conf.int[[2]]}) )
abline(lm(as.numeric(tapply(d_subset$entropy[d_subset$eventsum>=5], d_subset$year[d_subset$eventsum>=5], mean, na.rm=T)) ~ c(1996:2016)), lty=2, col="gray30")

plot(1996:2016, tapply(d_subset$entropy[d_subset$eventsum>=10], d_subset$year[d_subset$eventsum>=10], mean, na.rm=T),type="n",
	xlab="Year",ylab="Entropy",ylim=c(1.4,2.0), main="Annual Average Entropy:\nObservations with >=10 Events")
points(1996:2016, tapply(d_subset$entropy[d_subset$eventsum>=10], d_subset$year[d_subset$eventsum>=10], mean, na.rm=T),pch=19)
abline(h=seq(1.4,2,0.1), col="gray80")
segments(x0=1996:2016, x1=1996:2016, 
	y0=tapply(d_subset$entropy[d_subset$eventsum>=10], d_subset$year[d_subset$eventsum>=10], function(x){t.test(x)$conf.int[[1]]}),
	y1=tapply(d_subset$entropy[d_subset$eventsum>=10], d_subset$year[d_subset$eventsum>=10], function(x){t.test(x)$conf.int[[2]]}) )
abline(lm(as.numeric(tapply(d_subset$entropy[d_subset$eventsum>=10], d_subset$year[d_subset$eventsum>=10], mean, na.rm=T)) ~ c(1996:2016)), lty=2, col="gray30")


plot(1996:2016, tapply(d_subset$entropy[d_subset$country %in% country_select_events5], d_subset$year[d_subset$country %in% country_select_events5], mean, na.rm=T),type="n",
	xlab="Year",ylab="Entropy",ylim=c(1.4,2.0), main="Annual Average Entropy:\nCountries with Consistently >=5 Events")
points(1996:2016, tapply(d_subset$entropy[d_subset$country %in% country_select_events5], d_subset$year[d_subset$country %in% country_select_events5], mean, na.rm=T),pch=19)
abline(h=seq(1.4,2,0.1), col="gray80")
segments(x0=1996:2016, x1=1996:2016, 
	y0=tapply(d_subset$entropy[d_subset$country %in% country_select_events5], d_subset$year[d_subset$country %in% country_select_events5], function(x){t.test(x)$conf.int[[1]]}),
	y1=tapply(d_subset$entropy[d_subset$country %in% country_select_events5], d_subset$year[d_subset$country %in% country_select_events5], function(x){t.test(x)$conf.int[[2]]}) )
abline(lm(as.numeric(tapply(d_subset$entropy[d_subset$country %in% country_select_events5], d_subset$year[d_subset$country %in% country_select_events5], mean, na.rm=T)) ~ c(1996:2016)), lty=2, col="gray30")


plot(1996:2016, tapply(d_subset$entropy[d_subset$country %in% country_select_events10], d_subset$year[d_subset$country %in% country_select_events10], mean, na.rm=T),type="n",
	xlab="Year",ylab="Entropy",ylim=c(1.4,2.0), main="Annual Average Entropy:\nCountries with Consistently >=10 Events")
points(1996:2016, tapply(d_subset$entropy[d_subset$country %in% country_select_events10], d_subset$year[d_subset$country %in% country_select_events10], mean, na.rm=T),pch=19)
abline(h=seq(1.4,2,0.1), col="gray80")
segments(x0=1996:2016, x1=1996:2016, 
	y0=tapply(d_subset$entropy[d_subset$country %in% country_select_events10], d_subset$year[d_subset$country %in% country_select_events10], function(x){t.test(x)$conf.int[[1]]}),
	y1=tapply(d_subset$entropy[d_subset$country %in% country_select_events10], d_subset$year[d_subset$country %in% country_select_events10], function(x){t.test(x)$conf.int[[2]]}) )
abline(lm(as.numeric(tapply(d_subset$entropy[d_subset$country %in% country_select_events10], d_subset$year[d_subset$country %in% country_select_events10], mean, na.rm=T)) ~ c(1996:2016)), lty=2, col="gray30")
dev.off()





# plot of global entropy over time

# overall global entropy over time -- putting together all events each year, regardless of country
year_global_entropy<-NULL
year_global_props<-NULL
for(i in 1996:2016){
	#colsumsbyyear<-colSums(d[d$year==i,grep("govcit_CAMEO",names(d))])
	colsumsbyyear<-colSums(d[d$year==i,cameocodes])
	propsbyyear<-colsumsbyyear/sum(colsumsbyyear)
	year_global_props<-cbind(year_global_props, propsbyyear)
	year_global_entropy<-c(year_global_entropy, entropy_fun(propsbyyear) )
}


# resample to get CIs
year_sample_lo<-NULL
year_sample_hi<-NULL

for(j in 1996:2016){
yearsums<-colSums(d[d$year==j,cameocodes])
yearvec<-NULL
for(i in 1:length(yearsums)){
	yearvec<-c(yearvec, rep(names(yearsums)[i], yearsums[i]) )
}
if(length(yearvec)!=sum(yearsums)) print("error")

sample_entropy<-NULL
for(i in 1:1000){
year_newsample<-sample(yearvec, size=length(yearvec), replace=T)
sample_entropy<-c(sample_entropy, entropy_fun(table(year_newsample)/length(yearvec)) )
}
year_sample_lo<-c(year_sample_lo, sort(sample_entropy)[25] )
year_sample_hi<-c(year_sample_hi, sort(sample_entropy)[975] )
}


tiff(file="annualglobalentropy_new.tiff",height=7,width=8, units="in",res=300)
plot(1996:2016, year_global_entropy,type="n", 	xlab="Year",ylab="Global Entropy",ylim=c(1.9,2.6), main="Annual Global Entropy:\nAll Repression Events" )
abline(h=seq(1.9,2.6,0.1), col="gray80")
points(1996:2016, year_global_entropy, pch=19)
segments(x0=1996:2016, x1=1996:2016, 
	y0=year_sample_lo,year_sample_hi )
abline(lm(year_global_entropy ~ c(1996:2016)), lty=2, col="gray30")
dev.off()



pdf(file="annualglobalentropy_new.pdf",height=7,width=8)
plot(1996:2016, year_global_entropy,type="n", 	xlab="Year",ylab="Global Entropy",ylim=c(1.9,2.6), main="Annual Global Entropy:\nAll Repression Events" )
abline(h=seq(1.9,2.6,0.1), col="gray80")
points(1996:2016, year_global_entropy, pch=19)
segments(x0=1996:2016, x1=1996:2016, 
	y0=year_sample_lo,year_sample_hi )
abline(lm(year_global_entropy ~ c(1996:2016)), lty=2, col="gray30")
dev.off()



# plots of country averages over time for different measures

d_subset<-d[d$year %in% 1996:2016 & d$eventsum>=10 & !is.na(d$HRscores_rev),]
country_avgs1<-tapply(d_subset$entropy, d_subset$country, mean)

tiff(file="entropy_vs_fariss_entropy_vs_count_allyears_v2_new.tiff",height=7,width=13, units="in",res=300)
par(mfrow=c(1,2))

plot(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), tapply(log(d_subset$eventsum), d_subset$country, mean, na.rm=T), type="n", xlab="HR Scores (Reversed)",ylab="Log Event Count",main="HR Scores and Event Count\nCountry Averages 1996-2016")
text(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), tapply(log(d_subset$eventsum), d_subset$country, mean, na.rm=T), labels=names(country_avgs1), cex=0.7 )
abline(lm(tapply(log(d_subset$eventsum), d_subset$country, mean, na.rm=T) ~tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T)), lty=2)


plot(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_avgs1, type="n", xlab="HR Scores (Reversed)",ylab="Entropy",main="HR Scores and Repression Entropy\nCountry Averages 1996-2016")
text(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_avgs1, labels=names(country_avgs1), cex=0.7 )
abline(lm(country_avgs1~tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T)), lty=2)
dev.off()

pdf(file="entropy_vs_fariss_entropy_vs_count_allyears_v2_new.pdf",height=7,width=13)
par(mfrow=c(1,2))

plot(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), tapply(log(d_subset$eventsum), d_subset$country, mean, na.rm=T), type="n", xlab="HR Scores (Reversed)",ylab="Log Event Count",main="HR Scores and Event Count\nCountry Averages 1996-2016")
text(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), tapply(log(d_subset$eventsum), d_subset$country, mean, na.rm=T), labels=names(country_avgs1), cex=0.7 )
abline(lm(tapply(log(d_subset$eventsum), d_subset$country, mean, na.rm=T) ~tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T)), lty=2)


plot(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_avgs1, type="n", xlab="HR Scores (Reversed)",ylab="Entropy",main="HR Scores and Repression Entropy\nCountry Averages 1996-2016")
text(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_avgs1, labels=names(country_avgs1), cex=0.7 )
abline(lm(country_avgs1~tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T)), lty=2)
dev.off()




# illustrative country plots

d_subset<-d[d$year>1995 & d$country=="India",]
tiff(file="4timeplots_India.tiff",height=7,width=8, units="in",res=300)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="India:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="India:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="India:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="India:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Serbia",]
tiff(file="4timeplots_Serbia.tiff",height=7,width=8, units="in",res=300)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Yugoslavia/Serbia:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Yugoslavia/Serbia:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Yugoslavia/Serbia:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Yugoslavia/Serbia:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Turkey",]
tiff(file="4timeplots_Turkey.tiff",height=7,width=8, units="in",res=300)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Turkey:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Turkey:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Turkey:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Turkey:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Brazil",]
tiff(file="4timeplots_Brazil.tiff",height=7,width=8, units="in",res=300)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Brazil:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Brazil:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Brazil:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Brazil:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Kenya",]
tiff(file="4timeplots_Kenya.tiff",height=7,width=8, units="in",res=300)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Kenya:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Kenya:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Kenya:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Kenya:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Iran",]
tiff(file="4timeplots_Iran.tiff",height=7,width=8, units="in",res=300)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Iran:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Iran:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Iran:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Iran:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Afghanistan",]
tiff(file="4timeplots_Afghanistan.tiff",height=7,width=8, units="in",res=300)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Afghanistan:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Afghanistan:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Afghanistan:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Afghanistan:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Rwanda",]
tiff(file="4timeplots_Rwanda.tiff",height=7,width=8, units="in",res=300)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Rwanda:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Rwanda:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Rwanda:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Rwanda:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Singapore",]
tiff(file="4timeplots_Singapore.tiff",height=7,width=8, units="in",res=300)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Singapore:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Singapore:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Singapore:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Singapore:\nPolitical Terror Scale")
dev.off()




d_subset<-d[d$year>1995 & d$country=="India",]
pdf(file="4timeplots_India.pdf",height=7,width=8)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="India:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="India:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="India:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="India:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Serbia",]
pdf(file="4timeplots_Serbia.pdf",height=7,width=8)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Yugoslavia/Serbia:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Yugoslavia/Serbia:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Yugoslavia/Serbia:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Yugoslavia/Serbia:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Turkey",]
pdf(file="4timeplots_Turkey.pdf",height=7,width=8)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Turkey:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Turkey:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Turkey:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Turkey:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Brazil",]
pdf(file="4timeplots_Brazil.pdf",height=7,width=8)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Brazil:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Brazil:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Brazil:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Brazil:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Kenya",]
pdf(file="4timeplots_Kenya.pdf",height=7,width=8)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Kenya:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Kenya:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Kenya:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Kenya:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Iran",]
pdf(file="4timeplots_Iran.pdf",height=7,width=8)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Iran:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Iran:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Iran:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Iran:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Afghanistan",]
pdf(file="4timeplots_Afghanistan.pdf",height=7,width=8)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Afghanistan:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Afghanistan:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Afghanistan:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Afghanistan:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Rwanda",]
pdf(file="4timeplots_Rwanda.pdf",height=7,width=8)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Rwanda:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Rwanda:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Rwanda:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Rwanda:\nPolitical Terror Scale")
dev.off()

d_subset<-d[d$year>1995 & d$country=="Singapore",]
pdf(file="4timeplots_Singapore.pdf",height=7,width=8)
par(mfrow=c(2,2))
plot(1996:2016,d_subset$entropy,type="l", ylim=c(0,3),
	ylab="Entropy", xlab="",main="Singapore:\nRepression Entropy")
plot(1996:2016,log(d_subset$eventsum),type="l",ylim=c(0,8.5),
	ylab="Log Count", xlab="",main="Singapore:\nLog Repression Event Count")
plot(1996:2016,d_subset$HRscores_rev,type="l",ylim=c(-2.5,2.5),
	ylab="HR Score", xlab="",main="Singapore:\nHuman Rights Score (Reversed)")
plot(1996:2016,d_subset$PTS_State,type="l",ylim=c(1,5),
	ylab="PTS Score", xlab="",main="Singapore:\nPolitical Terror Scale")
dev.off()






# alternate Goldstein repression score version of plots

d_subset<-d[d$year %in% 1996:2016 & d$eventsum>=10 & !is.na(d$HRscores_rev),]

tiff(file="entropy_vs_fariss_entropy_vs_goldsteinrep_allyears_v2_new.tiff",height=7,width=13, units="in",res=300)
par(mfrow=c(1,2))

plot(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), tapply(d_subset$goldsteinrepress, d_subset$country, mean, na.rm=T), type="n", xlab="HR Scores (Reversed)",ylab="Goldstein Scale (Reversed)",main="HR Scores and Goldstein Scale\nCountry Averages 1996-2016")
text(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), tapply(d_subset$goldsteinrepress, d_subset$country, mean, na.rm=T), labels=names(country_avgs1), cex=0.7 )
abline(lm(tapply(d_subset$goldsteinrepress, d_subset$country, mean, na.rm=T) ~tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T)), lty=2)


plot(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_avgs1, type="n", xlab="HR Scores (Reversed)",ylab="Entropy",main="HR Scores and Repression Entropy\nCountry Averages 1996-2016")
text(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_avgs1, labels=names(country_avgs1), cex=0.7 )
abline(lm(country_avgs1~tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T)), lty=2)

dev.off()


pdf(file="entropy_vs_fariss_entropy_vs_goldsteinrep_allyears_v2_new.pdf",height=7,width=13)
par(mfrow=c(1,2))

plot(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), tapply(d_subset$goldsteinrepress, d_subset$country, mean, na.rm=T), type="n", xlab="HR Scores (Reversed)",ylab="Goldstein Scale (Reversed)",main="HR Scores and Goldstein Scale\nCountry Averages 1996-2016")
text(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), tapply(d_subset$goldsteinrepress, d_subset$country, mean, na.rm=T), labels=names(country_avgs1), cex=0.7 )
abline(lm(tapply(d_subset$goldsteinrepress, d_subset$country, mean, na.rm=T) ~tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T)), lty=2)


plot(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_avgs1, type="n", xlab="HR Scores (Reversed)",ylab="Entropy",main="HR Scores and Repression Entropy\nCountry Averages 1996-2016")
text(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_avgs1, labels=names(country_avgs1), cex=0.7 )
abline(lm(country_avgs1~tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T)), lty=2)

dev.off()



# herfindahl comparisons

country_herf_avgs1<-tapply(-1*d_subset$herfindahl, d_subset$country, mean)

tiff(file="entropy_vs_fariss_herfindahl_vs_fariss_allyears_v2_new.tiff",height=7,width=13, units="in",res=300)
par(mfrow=c(1,2))
plot(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_avgs1, type="n", xlab="HR Scores (Reversed)",ylab="Entropy",main="HR Scores and Repression Entropy\nCountry Averages 1996-2016")
text(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_avgs1, labels=names(country_avgs1), cex=0.7 )
abline(lm(country_avgs1~tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T)), lty=2)

plot(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_herf_avgs1, type="n", xlab="HR Scores (Reversed)",ylab="Inverse Herfindahl",main="HR Scores and Repression Inverse Herfindahl\nCountry Averages 1996-2016")
text(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_herf_avgs1, labels=names(country_avgs1), cex=0.7 )
abline(lm(country_herf_avgs1~tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T)), lty=2)
dev.off()

pdf(file="entropy_vs_fariss_herfindahl_vs_fariss_allyears_v2_new.pdf",height=7,width=13)
par(mfrow=c(1,2))
plot(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_avgs1, type="n", xlab="HR Scores (Reversed)",ylab="Entropy",main="HR Scores and Repression Entropy\nCountry Averages 1996-2016")
text(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_avgs1, labels=names(country_avgs1), cex=0.7 )
abline(lm(country_avgs1~tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T)), lty=2)

plot(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_herf_avgs1, type="n", xlab="HR Scores (Reversed)",ylab="Inverse Herfindahl",main="HR Scores and Repression Inverse Herfindahl\nCountry Averages 1996-2016")
text(tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T), country_herf_avgs1, labels=names(country_avgs1), cex=0.7 )
abline(lm(country_herf_avgs1~tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T)), lty=2)
dev.off()


# compare R2
entrop<-lm(country_avgs1~tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T))
summary(entrop)

herf<-lm(country_herf_avgs1~tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T))
summary(herf)

# compare correlation
cor(country_avgs1,tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T))
cor(country_herf_avgs1,tapply(d_subset$HRscores_rev, d_subset$country, mean, na.rm=T))


tiff(file="entropy_vs_herfindahl_allyears_v2_new.tiff",height=7,width=8, units="in",res=300)
par(mfrow=c(1,1))
plot(country_herf_avgs1, country_avgs1, type="n", xlab="Inverse Herfindahl",ylab="Entropy",main="Repression: Inverse Herfindahl Vs. Entropy\nCountry Averages 1996-2016")
text(country_herf_avgs1, country_avgs1, labels=names(country_avgs1), cex=0.7 )
abline(lm(country_avgs1~country_herf_avgs1), lty=2)
dev.off()

tiff(file="entropy_vs_herfindahl_densities.tiff",height=7,width=13, units="in",res=300)
par(mfrow=c(1,2))
entden<-density(d_subset$entropy)
plot(entden, main="Kernel Density of Repression (Entropy)")
polygon(entden, col="red", border="blue") 
herfden<-density(-1*d_subset$herfindahl)
plot(herfden, main="Kernel Density of Repression (Inverse Herfindahl)")
polygon(herfden, col="blue", border="red") 
dev.off()


pdf(file="entropy_vs_herfindahl_allyears_v2_new.pdf",height=7,width=8)
par(mfrow=c(1,1))
plot(country_herf_avgs1, country_avgs1, type="n", xlab="Inverse Herfindahl",ylab="Entropy",main="Repression: Inverse Herfindahl Vs. Entropy\nCountry Averages 1996-2016")
text(country_herf_avgs1, country_avgs1, labels=names(country_avgs1), cex=0.7 )
abline(lm(country_avgs1~country_herf_avgs1), lty=2)
dev.off()

pdf(file="entropy_vs_herfindahl_densities.pdf",height=7,width=13)
par(mfrow=c(1,2))
entden<-density(d_subset$entropy)
plot(entden, main="Kernel Density of Repression (Entropy)")
polygon(entden, col="red", border="blue") 
herfden<-density(-1*d_subset$herfindahl)
plot(herfden, main="Kernel Density of Repression (Inverse Herfindahl)")
polygon(herfden, col="blue", border="red") 
dev.off()


library(e1071)
kurtosis(d_subset$entropy)
kurtosis(d_subset$herfindahl)










##########

# MODELS

# main results table A.II

# restrict to observations with at least 10 events observed
d1<-d[d$year>1995 & d$eventsum>=10, ]


model2a<-lm(entropy~ protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount
	,data=d1)

model2b<-lm(entropy~ protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_lag
	,data=d1)

model2c<-lm(entropy~ protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_lag+time
	,data=d1)

model2d<-lm(entropy~ protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_lag+as.factor(year)
	,data=d1)

model2e<-lm(entropy~ protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_lag+as.factor(year)+as.factor(e_regionpol)
	,data=d1)

model2f<-lm(entropy~ protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_lag+as.factor(year)+as.factor(country)
	,data=d1)

model2g<-lm(entropy~ protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_lag+HRscores_rev
	,data=d1)

model2h<-lm(entropy~ protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_lag+PTS_State
	,data=d1)

model2i<-lm(entropy~ protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_lag+CIRI_rev
	,data=d1)

model2j<-lm(entropy~ protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_lag+vdem_HR_rev
	,data=d1)

stargazer(model2a, model2b, model2c, model2d, model2e, model2f, model2g, model2h, model2i, model2j,
	se=list(
		sqrt(diag(vcovCL(model2a, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model2b, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model2c, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model2d, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model2e, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model2f, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model2g, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model2h, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model2i, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model2j, cluster = d1$country, type = "HC1")))
		),
		add.lines=list(c("Year FE", "", "", "", "X", "X", "X", "", "", "", ""),
		c("Region FE",				"", "", "",  "",  "X", "",  "", "", "", ""),
		 c("Country FE",			"", "", "",  "",  "", "X", "", "", "", "")),
	 align=T, type="latex", omit=c("factor"), omit.stat=c("f", "ser", "rsq"), column.sep.width="0pt", font.size="footnotesize", no.space=T, 
	 star.char = c("\\dagger", "*", "**"), notes="\\dagger p<0.1; * p<0.05; ** p<0.01", notes.append=F)



# results table A.III: alternate events-based DVs

model3a<-lm(entropyNOFILTER~ 
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropyNOFILTER_lag
	,data=d1[d1$eventsumNOFILTER>=10,])

model3b<-lm(entropy_2digit~ 
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_2digit_lag
	,data=d1)

model3c<-lm(entropy_omit173~ 
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_omit173_lag
	,data=d1)

model3d<-lm(code173share~ 
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	code173share_lag
	,data=d1)

model3e<-lm(herfindahl~ 
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	herfindahl_lag
	,data=d1)

model3f<-lm(entropy_mm~ 
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_mm_lag
	,data=d1)

model3g<-lm(entropy_pub95~ 
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_pub95_lag
	,data=d1[d1$eventsum_pub95>=10,])

stargazer(model3a, model3b, model3c, model3d, model3e, model3f, model3g,
	se=list(
		sqrt(diag(vcovCL(model3a, cluster = d1$country[d1$eventsumNOFILTER>=10], type = "HC1"))),
		sqrt(diag(vcovCL(model3b, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model3c, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model3d, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model3e, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model3f, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model3g, cluster = d1$country[d1$eventsum_pub95>=10], type = "HC1")))
		),
	 align=T, type="latex", omit=c("factor"), omit.stat=c("f", "ser", "rsq"), column.sep.width="0pt", font.size="footnotesize", no.space=T,
	 star.char = c("\\dagger", "*", "**"), notes="\\dagger p<0.1; * p<0.05; ** p<0.01", notes.append=F)



# results table A.IV: alternate non-events DVs

model4a<-lm(logeventcount~ 
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+
	logeventcount_lag
	,data=d1)

model4b<-lm(HRscores_rev~
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+
	HRscores_rev_lag
	,data=d1)

model4c<-lm(PTS_State~
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+
	PTS_State_lag
	,data=d1)

model4d<-lm(CIRI_rev~
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+
	CIRI_rev_lag
	,data=d1)

model4e<-lm(goldsteinrepress~
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+
	goldsteinrepress_lag
	,data=d1)

model4f<-lm(vdem_HR_rev~
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+
	vdem_HR_rev_lag
	,data=d1)

stargazer(model4a, model4b, model4c, model4d, model4e, model4f,
	se=list(
		sqrt(diag(vcovCL(model4a, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model4b, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model4c, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model4d, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model4e, cluster = d1$country, type = "HC1"))),
		sqrt(diag(vcovCL(model4f, cluster = d1$country, type = "HC1")))
		),
	 align=T, type="latex", omit=c("factor"), omit.stat=c("f", "ser", "rsq"), column.sep.width="0pt", font.size="footnotesize", no.space=T,
	 star.char = c("\\dagger", "*", "**"), notes="\\dagger p<0.1; * p<0.05; ** p<0.01", notes.append=F)



# results table A.V: varying sample restrictions.

model5a<-lm(entropy~
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_lag
	,data=d[d$year>1995 & d$eventsum>0, ])

model5b<-lm(entropy~
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_lag
	,data=d[d$year>1995 & d$eventsum>=5, ])

model5c<-lm(entropy~
	protest_dummy+civwar+p_polity2+treatytotal1+loggdpcap+logpop+logeventcount+
	entropy_lag
	,data=d[d$year>1995 & d$eventsum>=20, ])

stargazer(model5a, model5b, model5c,
	se=list(
		sqrt(diag(vcovCL(model5a, cluster = d$country[d$year>1995 & d$eventsum>0], type = "HC1"))),
		sqrt(diag(vcovCL(model5b, cluster = d$country[d$year>1995 & d$eventsum>=5], type = "HC1"))),
		sqrt(diag(vcovCL(model5c, cluster = d$country[d$year>1995 & d$eventsum>=20], type = "HC1")))
		),
	 align=T, type="latex", omit=c("factor"), omit.stat=c("f", "ser", "rsq"), column.sep.width="0pt", font.size="footnotesize", no.space=T,
	 star.char = c("\\dagger", "*", "**"), notes="\\dagger p<0.1; * p<0.05; ** p<0.01", notes.append=F)














