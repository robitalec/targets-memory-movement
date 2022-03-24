############################
#####LOAD PACKAGES##########
############################
library(adehabitatHR)
library(adehabitatLT)
library(circular)
library(raster)
library(rgeos)
library(geosphere)
library(rgdal)

###########################
###LOADING DATA############
###########################
#ALL DATA
head(data) 
hog=data[data$id=="150011"|data$id=="150041",] #GETTING ONLY TWO INDIVIDUALS (a short example to increase velocity)
hog$id=as.vector(hog$id)
#MAP
map=raster("Map/map.tif") #HABITAT RASTER

############################
#####DATA PREPAIRING########
############################
ids=unique(hog$id)
memory=random_memory=vector("list",length(ids))
names(memory)=ids

for(j in 1:length(ids)){ #FOR BY INDIVIDUAL
temp=hog[hog$id==ids[j],]

#################################################
#TRAJECTORY - SUBSAMPLING PER HOUR###############
################TURNING ANGLES AND STEP LENGTHS##
#################################################
traj=as.ltraj(temp[,c("x","y")],id=temp$id,date=temp$date_time)
refda <- temp$date_time[1]
traj1=setNA(traj,refda,1/12,units="hour")
traj1=sett0(traj1, refda, 1/12, units="hour")
traj1=redisltraj(na.omit(traj1), 1/12*3600, type = "time")
traj2=redisltraj(na.omit(traj1), 1*60*60, type = "time")
data1=ld(traj2)

##################################################
#########HABITAT USED#############################
#################INTERCECTING RELOCATION WITH MAP#
##################################################
#class=data.frame(map.index=1:length(map@data@attributes[[1]]$CLASS),habitat=map@data@attributes[[1]]$CLASS)
#habitat=data.frame(map.index=extract(map,SpatialPoints(data1[,c("x","y")])));habitat$id=1:nrow(habitat)
#habitat=merge(habitat,class,all.x=T,all.y=F)
#habitat=data.frame(habitat=habitat[order(habitat$id),"habitat"])
#habitat$habitat[habitat$habitat=="Wet High Probability"]="Wet Medium Probability"
#habitat$habitat[habitat$habitat=="Wet Low Probability"]="Dry"
#habitat$habitat=factor(habitat$habitat)

######################################
######################################
######BROWNIAN BRIDGE KERNEL##########
################INTENSIVITY DENSITY###
######################################
######################################
hour=as.POSIXlt(data1$date)$hour
datelt=as.POSIXlt(data1$date)$yday
start=data1[datelt>unique(datelt)[4],]#STARTING INTHE FOUTH DAY
l=as.numeric(rownames(start)[1])

d_rm=d_rmt=d_lm=d_lmt=numeric()
random_end=matrix(,,18);random_end=random_end[-1,]
colnames(random_end)=c("id","step","x","y","date","dist","rel.angle","abs.ablge","R_rm","R_rmt","R_lm","R_lmt","habitat","c1","c2","s1","s2","presabs") 
for(i in 1:(nrow(start)-1)){ #FOR BY STEP
day=as.POSIXlt(start$date[i])

####################################
#CUTTING DATA FOR MEMORY ESTIMATION#
####################################
rm=dl(data1[data1$date<day & data1$date>day-3*24*60*60,]) #HERE IS THE RECENT MEMORY DELAY 3 DAYS
lm=dl(data1[data1$date<day,]) #LONG MEMORY
rmd=dl(rm[[1]][as.POSIXlt(rm[[1]]$date)$hour>=6 & as.POSIXlt(rm[[1]]$date)$hour<=18,]) #RECENT MEMORY DAY
rmn=dl(rm[[1]][as.POSIXlt(rm[[1]]$date)$hour<6 | as.POSIXlt(rm[[1]]$date)$hour>18,]) #RECENT MEMORY NIGHT
lmd=dl(lm[[1]][as.POSIXlt(lm[[1]]$date)$hour>=6 & as.POSIXlt(lm[[1]]$date)$hour<=18,]) #LONG MEMORY DAY
lmn=dl(lm[[1]][as.POSIXlt(lm[[1]]$date)$hour<6 | as.POSIXlt(lm[[1]]$date)$hour>18,])  #LONG MEMORY NIGHT

####################
#MEMORY ESTIMATION##
####################
id_rm=BRB(rm, D=3, Tmax=2*60*60, Lmin=15, hmin=30, type="ID", radius = 30,maxt = 0.5*60*60,grid=80)
id_rm1=getvolumeUD(id_rm)

id_rmd=BRB(rmd, D=3, Tmax=2*60*60, Lmin=15, hmin=30, type="ID", radius = 30,maxt = 0.5*60*60,grid=80)
id_rmd1=getvolumeUD(id_rmd)

id_rmn=BRB(rmn, D=3, Tmax=2*60*60, Lmin=15, hmin=30, type="ID", radius = 30,maxt = 0.5*60*60,grid=80)
id_rmn1=getvolumeUD(id_rmn)

id_lm=BRB(lm, D=3, Tmax=2*60*60, Lmin=15, hmin=30, type="ID", radius = 30,maxt = 0.5*60*60,grid=80)
id_lm1=getvolumeUD(id_lm)

id_lmd=BRB(lmd, D=3, Tmax=2*60*60, Lmin=15, hmin=30, type="ID", radius = 30,maxt = 0.5*60*60,grid=80)
id_lmd1=getvolumeUD(id_lmd)

id_lmn=BRB(lmn, D=3, Tmax=2*60*60, Lmin=15, hmin=30, type="ID", radius = 30,maxt = 0.5*60*60,grid=80)
id_lmn1=getvolumeUD(id_lmn)

##############################
##############################
######SAMPLE RANDOM STEPS#####
##############################
##############################

###########################################
#EXTRACTING TURNING ANGLE AND STEP LENGTHS#
###########################################
tangle=circular(data1$rel.angle,units="rad")
tangle=(conversion.circular(tangle,units="degrees"))[!is.na(tangle)]
dist=data1$dist;dist=dist[!is.na(dist)]

###########################
#PROJECT 50 RANDOMS STEPS##
###########################
current_location=SpatialPoints(start[i,c("x","y")],proj4string=CRS("+proj=utm +datum=WGS84 +zone=21 +south"))
future_location=SpatialPoints(start[i+1,c("x","y")],proj4string=CRS("+proj=utm +datum=WGS84 +zone=21 +south"))
projected_locations=destPoint(p=spTransform(current_location,CRS=CRS("+proj=longlat")),b=sample(tangle,50,T),d=sample(dist,50,T))
projected_sp=SpatialPoints(projected_locations,proj4string=CRS("+proj=longlat +datum=WGS84"))
random=spTransform(projected_sp,CRS=CRS("+proj=utm +datum=WGS84 +zone=21 +south"))
random=SpatialPoints(data.frame(rbind(future_location,random)))

################################################################
#INTERCEPTING OBSERVED AND 50 RANDOM STEPS WITH THE MEMORY MAP##
################################################################
R_rm=adehabitatMA::join(random,id_rm1)###RECENT MEMORY
if(length(is.na(R_rm)[is.na(R_rm)==TRUE])>0){
	for(e in seq(1,50,0.5)){
	id_rm=BRB(rm, D=3, Tmax=2*60*60, Lmin=15, hmin=30, type="ID", radius = 30,maxt = 0.5*60*60,grid=80,extent=e)
	id_rm1=getvolumeUD(id_rm)
	R_rm=adehabitatMA::join(random,id_rm1)
		if(length(is.na(R_rm)[is.na(R_rm)==TRUE])==0){"break"()}
	}
}

if(day$hour>=6 & day$hour<=18){R_rmt=adehabitatMA::join(random,id_rmd1)##RECENT DAY
	if(length(is.na(R_rmt)[is.na(R_rmt)==TRUE])>0){
		for(e in seq(1,50,0.5)){
		id_rmd=BRB(rmd, D=3, Tmax=2*60*60, Lmin=15, hmin=30, type="ID", radius = 30,maxt = 0.5*60*60,grid=80,extent=e)
		id_rmd1=getvolumeUD(id_rmd)
		R_rmt=adehabitatMA::join(random,id_rmd1)
			if(length(is.na(R_rmt)[is.na(R_rmt)==TRUE])==0){"break"()}
		}
	}
}

if(day$hour<6 | day$hour>18){R_rmt=adehabitatMA::join(random,id_rmn1)##RECENT NIGTH
	if(length(is.na(R_rmt)[is.na(R_rmt)==TRUE])>0){
		for(e in seq(1,50,0.5)){
		id_rmn=BRB(rmn, D=3, Tmax=2*60*60, Lmin=15, hmin=30, type="ID", radius = 30,maxt = 0.5*60*60,grid=80,extent=e)
		id_rmn1=getvolumeUD(id_rmn)
		R_rmt=adehabitatMA::join(random,id_rmn1)
			if(length(is.na(R_rmt)[is.na(R_rmt)==TRUE])==0){"break"()}
		}
	}
}


R_lm=adehabitatMA::join(random,id_lm1)##LONG MEMORY
if(length(is.na(R_lm)[is.na(R_lm)==TRUE])>0){
	for(e in seq(1,50,0.5)){
	id_lm=BRB(lm, D=3, Tmax=2*60*60, Lmin=15, hmin=30, type="ID", radius = 30,maxt = 0.5*60*60,grid=80,extent=e)
	id_lm1=getvolumeUD(id_lm)
	R_lm=adehabitatMA::join(random,id_lm1)
		if(length(is.na(R_lm)[is.na(R_lm)==TRUE])==0){"break"()}
	}
}

if(day$hour>=6 & day$hour<=18){R_lmt=adehabitatMA::join(random,id_lmd1)##LONG DAY
	if(length(is.na(R_lmt)[is.na(R_lmt)==TRUE])>0){
		for(e in seq(1,50,0.5)){
		id_lmd=BRB(lmd, D=3, Tmax=2*60*60, Lmin=15, hmin=30, type="ID", radius = 30,maxt = 0.5*60*60,grid=80,extent=e)
		id_lmd1=getvolumeUD(id_lmd)
		R_lmt=adehabitatMA::join(random,id_lmd1)
			if(length(is.na(R_lmt)[is.na(R_lmt)==TRUE])==0){"break"()}
		}
	}
}

if(day$hour<6 | day$hour>18){R_lmt=adehabitatMA::join(random,id_lmn1)##LONG NIGHT
	if(length(is.na(R_lmt)[is.na(R_lmt)==TRUE])>0){
		for(e in seq(1,50,0.5)){
		id_lmn=BRB(lmn, D=3, Tmax=2*60*60, Lmin=15, hmin=30, type="ID", radius = 30,maxt = 0.5*60*60,grid=80,extent=e)
		id_lmn1=getvolumeUD(id_lmn)
		R_lmt=adehabitatMA::join(random,id_lmn1)
			if(length(is.na(R_lmt)[is.na(R_lmt)==TRUE])==0){"break"()}
		}
	}
}

####################################################
#########HABITAT AVAILABLE##########################
#################INTERCECTING RANDOM STEPS WITH MAP#
####################################################
habitat_random=data.frame(map.index=extract(map,random));habitat_random$id=1:nrow(habitat_random)
habitat_random=merge(habitat_random,class,all.x=T,all.y=F)
habitat_random=data.frame(habitat=habitat_random[order(habitat_random$id),"habitat"])
habitat_random$habitat[habitat_random$habitat=="Wet High Probability"]="Wet Medium Probability"
habitat_random$habitat[habitat_random$habitat=="Wet Low Probability"]="Dry"
habitat_random$habitat=factor(habitat_random$habitat)

#############
#HARMONICS###
#############
hours=as.POSIXlt(start[i,]$date)$hour
c1r=cos(hours*2*pi/24)
c2r=cos(hours*4*pi/24)
s1r=sin(hours*2*pi/24)
s2r=sin(hours*4*pi/24)

######################################################
#STORING THE "FOR" WITHIN INDIVIDUALS - STEP BY STEP##
######################################################
random_res=data.frame(id=start[i,"id"],step=i,x=start[i,"x"],y=start[i,"y"],date=start[i,"date"],dist=start[i,"dist"],rel.angle=start[i,"rel.angle"],
abs.ablge=start[i,"abs.angle"],R_rm,R_rmt,R_lm,R_lmt,habitat_random,c1=c1r,c2=c2r,s1=s1r,s2=s2r,presabs=c(1,rep(0,50)))
random_end=rbind(random_end,random_res)
cat(ids[j],i,round(i/nrow(start),3),"\n")
}
##################################
#STORING THE "FOR" OF INDIVIDUALS#
##################################
random_memory[[j]]=random_end
}

####################################################
#THE RESULT IS A LIST (EACH ELEMNT IS A INDIVIDUAL)#
#FROM LIST TO DATAFRAME#############################
####################################################
ssf=do.call("rbind",random_memory)

#########################
#SSF MEMORY_BIASED MODEL#
#########################
#Running an example - Recent Memory SSF
model=coxph(Surv(step,presabs)~habitat*c1+habitat*c2+habitat*s1+habitat*s2+R_rm+strata(paste(id,step)),data=ssf) #RECENT MEMORY
summary(model)
source("ssf_graph.R")
plot.ssf(model,"Forest",col=3)



