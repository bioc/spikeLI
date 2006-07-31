Ivsc<- function(probe_set,probe="NULL",outfile="NULL"){

if (!exists("hgu")) {
    cat("\nhgu does not exist yet please load the value using data(hgu)\n");
    quit("yes");
}
if (!exists("SPIKE_INA")) {
    cat("\nSPIKE_INA does not exits, please use data(SPIKE_INA) to set the object\n");
    quit("yes");
}
##    hgu<-read.table("HGU");
    probe_set<-probe_set[1]
    if(probe=="NULL") probe<-1
    else probe<-probe[1]

    lgth<-length(hgu$Ipm[hgu$Probe.Set.Name==probe_set & hgu$conc==0])
    lgth<-min(lgth,length(hgu$Ipm[hgu$Probe.Set.Name==probe_set & hgu$conc==16]))
    lgth<-min(lgth,length(hgu$Ipm[hgu$Probe.Set.Name==probe_set & hgu$conc==64]))

    if(probe>lgth){cat("Error: probe input number cannot exceed \n")
		   cat(paste("       the probe length (=",lgth,")\n\n",sep=""))}
    concts<-hgu$conc[hgu$Probe.Set.Name==probe_set]
    concts<-levels(factor(concts))
    concts<-as.double(concts)
    
    par(mfrow=c(1, 1));
    Y<-c(100,15000)
    laby<-"I"
    labx<-"c (picoM)"
    tit<-paste(probe_set,probe,sep="")
    plot(range(concts)*1.1,Y,type="n",log="y",main=tit,las=1,xlab=labx,ylab=laby);

    dIpm<-double()
    dImm<-double()
    concs<-double()

    for(i in (1:length(concts)))
         {conc<-concts[i]
         Ipm<-hgu$Ipm[hgu$Probe.Set.Name==probe_set & hgu$conc==conc]
         Imm<-hgu$Imm[hgu$Probe.Set.Name==probe_set & hgu$conc==conc]
         kp<-length(Ipm)/lgth
         cc<-seq(probe,lgth*kp,lgth)
         dIpm<-c(dIpm,Ipm[cc])
         dImm<-c(dImm,Imm[cc])
         concs<-c(concs,rep(conc,times=kp))
        }

      points(concs,dIpm,col="blue",lwd=2.3)
      points(concs,dImm,col="red",lwd=2.3)
      text(range(concts)[2],600,"PM",col="blue")
      text(range(concts)[2],400,"MM",col="red")

## Non-linear fitting with the Langmuir isotherm

##     cat("----\nNon-linear fit I vs. c for the PM\n----\n")
      nonlin.fit<-nls(dIpm ~ I0 + A*concs/(K+concs),start=list(I0=1,A=10000,K=200))

      lines(concs,fitted.values(nonlin.fit))
      fit_PM<-fitted.values(nonlin.fit)

##       Imax_PM<-as.integer(summary(nonlin.fit)$parameters[2])
      Imax_PM<-as.integer(summary(nonlin.fit)$parameters[1]+summary(nonlin.fit)$parameters[2])
      text(range(concts)[2]*0.9,200,paste("Imax_PM=",Imax_PM))
     
##      cat("----\nNon-linear fit I vs. c for the MM\n----\n")
      nonlin.fit<-nls(dImm ~ I0 + A*concs/(K+concs),start=list(I0=1,A=10000,K=200))
      lines(concs,fitted.values(nonlin.fit),lty=2)
      fit_MM<-fitted.values(nonlin.fit)

      Imax_MM<-as.integer(summary(nonlin.fit)$parameters[1]+summary(nonlin.fit)$parameters[2])
      text(range(concts)[2]*0.9,150,paste("Imax_MM=",Imax_MM))

      if(outfile!="NULL")
        {cat("# PM\n",file=outfile)
            for(i in (1:length(concs))) 
                cat(paste(concs[i],dIpm[i],"\n"),file=outfile,append=TRUE)
         cat("&\n",file=outfile,append=TRUE)
         for(i in (1:length(concs))) 
            cat(paste(concs[i],dImm[i],"\n"),file=outfile,append=TRUE)
             cat("&\n",file=outfile,append=TRUE)
             for(i in (1:length(concs))) 
                 cat(paste(concs[i],fit_PM[i],"\n"),file=outfile,append=TRUE)
                 cat("&\n",file=outfile,append=TRUE)
                 for(i in (1:length(concs))) 
                     cat(paste(concs[i],fit_MM[i],"\n"),file=outfile,append=TRUE)
                     cat("&\n",file=outfile,append=TRUE)
## 	          cat(paste(),file=outfile,append=TRUE)
         }
}
