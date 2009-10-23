# Copyright (C) 2009 Peter Gething
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

require(deldir)
require(seqinr)
##########################################################################
getdeclusteredsample<-function(tablepath,prop,minSample=c(),DECLUSTER=TRUE,MAKEPLOT=FALSE){
               

 ## handle paramter passes from python and set defaults where necesssary
    if(DECLUSTER=="False") DECLUSTER<-FALSE 
    if(DECLUSTER=="True") DECLUSTER<-TRUE
    if(MAKEPLOT=="False") MAKEPLOT<-FALSE 
    if(MAKEPLOT=="True") MAKEPLOT<-TRUE 
 
 ## check required packages are present 
    ip=installed.packages()
    if(!any(ip[,1]=="seqinr")){
        print("ERROR in get_declustered_sample!!! Package seqinr not loaded - returning error")
        return(-9999)
    }    
    if(!any(ip[,1]=="deldir")){
        print("ERROR in get_declustered_sample!!! Package deldir not loaded - returning error")
        return(-9999)
    }    


 ## read in table (if full path not given, assume we are using the current working duirectry)

  # ensure all slashes are forward 
    tempvec<-s2c(tablepath)
    tempvec[which(tempvec=='\\')]<-"/"
    tablepath<-c2s(tempvec)

  # if there are any slashes then we assume this is the full path and import directly as stated
    if(any(tempvec=='/')){
        fullpath=tablepath
        tableIN<-read.csv(fullpath)
    }
        
  # if not then we assume the given filename is residing in the curent working directory
    if(!any(tempvec=='/')){
        wd<-system("pwd",intern=T)
        fullpath<-paste(wd,"/",tablepath,sep="")
        print(paste("WARNING from get_declustered_sample!! not given path so assuming file",tablepath,"is in working directory:",wd))
        tableIN<-read.csv(fullpath)
    }
    
   # decompose full path into directory and file prefix for later use
   temp<-strsplit(fullpath,"/")[[1]][length((strsplit(fullpath,"/"))[[1]])]
   filesuffix<-strsplit(temp,".csv")[[1]][1]
   directory<-strsplit(fullpath,filesuffix)[[1]][1]

 ## check for and define lat and lon columns
    if( (!any(names(tableIN)=="lon")) | (!any(names(tableIN)=="lat"))){
        print(paste("ERROR in get_declustered_sample!!! table input from path",tablepath,"does not contain lon and lat columns, returning error!!!"))
        return(-9999)
    }
    lon<-tableIN$lon
    lat<-tableIN$lat

 ## check for and define t column
    if( (any(names(tableIN)=="t"))){
        ST<-TRUE
        time<-tableIN$t

        #  order table by time of survey (so most recent surveys come first)
        Order<-order(time,decreasing=TRUE)
        tableIN<-tableIN[Order,]
        
    }else{
        ST<-FALSE
        print("WARNING!! no 't' column founs, so assuming a spatial only database")
    }

  # now define index so can focus only on those records that are not spatial duplciates of preceding records
    SpatUniqueIndex<-!(duplicated(cbind(lon,lat)))
    Ndata<-nrow(tableIN[SpatUniqueIndex,])
    
  # initialise null weight if not declustering  
    weights<-rep(1,Ndata)
    
 ## generate Thiessen polygons around points and define weights using area of each polygon
    if (DECLUSTER){
        ThiessObj<-deldir(lon[SpatUniqueIndex],lat[SpatUniqueIndex],frac=0.0000000000000000001) 
        ThiessSummary<-ThiessObj$summary
        weights<-sqrt(ThiessSummary[,dimnames(ThiessSummary)[[2]]=="del.wts"])
    }

 ## determine size of validation set (defined by prop but optional floor of minSample)
    SampleSize<-ceiling(length(tableIN[,1])*prop)
    
    if(class(minSample)!="NULL"){
        SampleSize<-max(minSample,SampleSize)
        if(minSample>Ndata){
            print(paste("ERROR in get_declustered_sample!!! specified minimum sample size (",minSample,") exceeds available spatially unique locations (",Ndata,")",sep=""))
            return(-9999)
        }
    }

 ## draw sample with weights proportional to area
    listID<-1:Ndata
    print("range(weights):")
    print(range(weights))
    SampleID<-sample(listID, SampleSize, replace = FALSE, prob = weights^100)

 ## define index for full length table defining whether in validation or thinned set
    ValSetIndex<-rep(FALSE,nrow(tableIN))
    ValSetIndex[SpatUniqueIndex][SampleID]<-TRUE

 ## define hold-out tables for non-buffer data set
    tableHOLDOUT<-tableIN[ValSetIndex,]

 ## define thinned table for combined buffer and non-buffer set
    tableTHINNED<-tableIN[!ValSetIndex,]

 ## optionally create and export plot of all points and show those sampled
    if(MAKEPLOT){
        pdf(paste(directory,filesuffix,"_sampleMap.pdf"))
        plot(lon,lat,ylab="",xlab="",type="n")
        points(tableTHINNED$lon,tableTHINNED$lat,col=3,pch=16,cex=0.5)
        points(tableHOLDOUT$lon,tableHOLDOUT$lat,col=2,pch=16,cex=0.5)
        if(DECLUSTER)title(main=paste("Declustered ",prop," sample (in red)\n(",SampleSize," of ",nrow(tableIN)," with ",Ndata," spatially unique)",sep=""))
        if(!DECLUSTER)title(main=paste("Random (non-declustered) ",prop," sample (in red)\n(",SampleSize," of ",nrow(tableIN)," with ",Ndata," spatially unique)",sep=""))
        dev.off()
    }
    
 ## export metadata tables of complete, thinned and holdout sets
    write.table(tableHOLDOUT,paste(directory,filesuffix,"HOLDOUT.csv",sep=""),row.names=F,sep=",")
    write.table(tableTHINNED,paste(directory,filesuffix,"THINNED.csv",sep=""),row.names=F,sep=",")

 ## return data tallies
    return(0)
}


#getdeclusteredsample(tablepath="\\home\\pwg\\tempval\\testtable.csv",prop=0.1,minSample=30,DECLUSTER=TRUE,MAKEPLOT=TRUE)








