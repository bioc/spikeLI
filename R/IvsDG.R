IvsDG<- function(probe_set,conc,outfile="NULL"){

##    hgu<-read.table("HGU");
    probe_set<-probe_set[1]

## Determine the chipset type HGU95/HGU133

        fic<-hgu$FILE[hgu$Probe.Set.Name == probe_set];
	fichier<- fic

## Calculate I-I0 for the PM

        Ipm <- double();
        I0pm <- double();
        Ipm<-hgu$Ipm[hgu$Probe.Set.Name == probe_set & hgu$conc==conc];
        I0pm<- 0*Ipm+ hgu$Ipm[hgu$Probe.Set.Name == probe_set&hgu$conc == 0];
        difIpm<-Ipm-I0pm[1:length(Ipm)];

        for(i in (1:length(Ipm))) if(difIpm[i]<0.) difIpm[i]=1.e-8

## Calculate I-I0 for the MM

        Imm <- double();
        I0mm <- double();
        Imm<-hgu$Imm[hgu$Probe.Set.Name == probe_set & hgu$conc==conc];
        I0mm<- 0*Imm+hgu$Imm[hgu$Probe.Set.Name == probe_set&hgu$conc == 0];
        difImm<-Imm-I0mm[1:length(Imm)];

	for(i in (1:length(Imm))) if(difImm[i]<0.) difImm[i]=1.e-8

## idem : je fais une boucle qui prend les probe_sets les uns apres les autres

        gPM <- double();
        gMM <- double();
        gRNA <- double();

## Read free energies from the file

        gPM <-hgu$DGpm[hgu$Probe.Set.Name == probe_set & hgu$conc==0];
        gMM <-hgu$DGmm[hgu$Probe.Set.Name == probe_set & hgu$conc==0];
        gRNA<-hgu$DGRNA[hgu$Probe.Set.Name == probe_set & hgu$conc==0];

        gPM <-gPM[1:length(Ipm)]
        gMM <-gMM[1:length(Ipm)]
        gRNA <-gRNA[1:length(Ipm)]

## compte le nombre de probes
        nbProbes<-length(hgu$Ipm[hgu$Probe.Set.Name == probe_set&hgu$conc == 0])
        nbProbes<-min(nbProbes,length(hgu$Ipm[hgu$Probe.Set.Name == probe_set&hgu$conc == 8]))
        nbProbes<-min(nbProbes,length(hgu$Ipm[hgu$Probe.Set.Name == probe_set&hgu$conc == 64]))
        expTmp <- integer();

## Calculate alpha

        alpha <- (-gRNA+46.5)*(-0.686)
        alpha <- exp(alpha);
        alpha <- alpha + 1;
        alpha <- 1 / alpha;

        beta<-0.74;

        A<-10000;

        conc<-conc*1e-12

## Comparison between X and X'

        xPM <- conc*exp(beta*gPM);
        xMM <- conc*exp(beta*gMM);

## Langmuir isotherm

        xpPM <- 10^seq(-5, 4, 0.2);
        LangPM <- (10000*xpPM)/(1+xpPM);
        xpPM <- log(xpPM)/beta;
        dg<-seq(20,40,0.2)
        cexpdg<-conc*exp(beta*dg)
        Langdg<-A*cexpdg/(1+cexpdg)
        cexpdgp4<-cexpdg*4
        Langdgp4<-A*cexpdgp4/(1+cexpdgp4)
        cexpdgd4<-cexpdg/4
        Langdgd4<-A*cexpdgd4/(1+cexpdgd4)

    ##indice du tableau
        kn <- length(gPM)/nbProbes;
        txt <- rep(paste(seq(1, nbProbes, 1)), times=kn);
    ##indice experience ? afficher
        exp<-rep(paste(expTmp), times=kn);

    ##dimensionne le graphique
        X <- c(24, 38);
        Y <- c(100,2*10^4);

## Correction for hybridization in solution
        correct<-log(alpha)/beta
        gPM<-gPM+correct
        gMM<-gMM+correct

## Create the graphic window
        par(mfrow=c(1,2))

        plot(c(0,nbProbes),Y,type="n",log="y",main=probe_set,las=1,xlab="probes",ylab="I-I0");
        lines(seq(1,nbProbes,1),difIpm,col="blue",lwd=2.3)
        points(seq(1,nbProbes,1),difIpm,col="blue",pch=19)
        lines(seq(1,nbProbes,1),difImm,col="red",lwd=2.3)
        points(seq(1,nbProbes,1),difImm,col="red",pch=19)
        text(2,15000,"PM",col="blue",font=2,adj=0,cex=1.1);
        text(6,15000,"MM",col="red",font=2,adj=0,cex=1.1);
        plot(X,Y,log="y",type="n",main=probe_set,las=1,xlab="Delta G + RT log(alpha)",ylab="I-I0");
        text(2,200,paste("c=",conc*1e12," pM"),font=2);
##	text(2,200,paste("c=",conc*1e12," pM"),font=2,cex=1.8);
        text(gPM,difIpm,txt,col="blue",font=2);
##        text(gPM,difIpm,txt,col="blue",cex=1.1,font=2);
        text(gMM,difImm,txt,col="red",font=2);
##        text(gMM,difImm,txt,col="red",cex=1.1,font=2);
        lines(dg,Langdg)
        lines(dg,Langdgp4,lty=4,col="green")
        lines(dg,Langdgd4,lty=4,col="green")

## Fit for the concentration
        gtot<-c(gPM,gMM[gMM>0])
        difItot<-c(difIpm,difImm[gMM>0])
        conctot<-1e12*difItot*exp(-beta*gtot)/(A-difItot)
        difItot<-difItot[conctot>0]
        conctot<-conctot[conctot>0]	
##        print(conctot)

##afficher le nom du probe_set et la legende
        text(24.5,15000,paste("c=",conc*1e12," pM"),font=2,cex=1.1,adj=0);
        text(33,15000,fichier,font=2,adj=0,cex=1.1);
        text(30,150,paste("median=",round(median(conctot),0),"pM\n"),font=2,adj=0,cex=1.0);
        text(30,125,paste("mean=",round(exp(mean(log(conctot))),0),"pM\n"),
                        font=2,adj=0,cex=1.0);
##probe="NULL",outfile="NULL"
## 
        l1<-length(conctot)
        l2<-l1*0.8

        sdv<-(mean(log(conctot))-log(conctot))^2
        for(i in (1:(l1-l2))){conctot<-conctot[sdv!=max(sdv)];sdv<-sdv[sdv!=max(sdv)];}
##        print(conctot)
        text(30,105,paste("median_80=",round(median(conctot),0),"pM\n"),
                        font=2,adj=0,cex=1.0);
        text(30,90,paste("mean_80=",round(exp(mean(log(conctot))),0),"pM\n"),
                        font=2,adj=0,cex=1.0);

        if(outfile!="NULL"){
                cat("# Perfect matches \n",file=outfile)
                for(i in 1:length(gPM)) 
                        cat(paste(gPM[i],difIpm[i],"\n"),file=outfile,append=TRUE)
                        cat("&\n# Mismatches \n",file=outfile,append=TRUE)
                for(i in 1:length(gPM)) 
                        cat(paste(gMM[i],difImm[i],"\n"),file=outfile,append=TRUE)
                        cat("&\n",file=outfile,append=TRUE)
                for(i in 1:length(dg)) 
                        cat(paste(dg[i],Langdg[i],"\n"),file=outfile,append=TRUE)
                        cat("&\n",file=outfile,append=TRUE)
        }
}
