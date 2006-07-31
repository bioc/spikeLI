"collapse" <-
function(probe_set, param="NULL", probes="NULL", output="NULL",filename="NULL"){

##    hgu<-read.table("HGU");

if (!exists("hgu")) {
    cat("\nhgu does not exist yet please load the value using data(hgu)\n");
    quit("yes");
}
if (!exists("SPIKE_INA")) {
    cat("\nSPIKE_INA does not exits, please use data(SPIKE_INA) to set the object\n");
    quit("yes");
}

    nbGenes<-length(probe_set);

    if (output=="PS"){
        path <- system("pwd",TRUE);
        if (filename=="NULL")
                {postscript(paste(path,paste("/",probe_set[1],sep=""),sep=""));}
        else    {postscript(paste(path,paste("/",filename,sep=""),sep=""));}
    }
## check the length of the vector "probe_set" to visualize correctly
    else{
        if (nbGenes>2){
            par(mfrow=c(2, 2));
	    if (nbGenes>4){
	        probe_set<-probe_set[1:4];
                cat("Warning: The program cannot show more than 4 probe sets simultaneously!\n");
                cat("         Only the first 4 will be shown. If you want to see more\n");
		cat("         the option output=\"PS\", filename=\"file.ps\" prints \n");
		cat("         all probe sets into a postscript file \n");
	    }
        }else{par(mfrow=c(nbGenes, 2));}
    }


##  Loop of each probe set in probe_set
    for(probe_setTmp in probe_set){
        continue<-1;
## To which chipset is "probe_set": HGU95/HGU133?
        fic<-hgu$FILE[hgu$Probe.Set.Name == probe_setTmp];
        fichier<- fic[1]
## Calculation of I-I0 for PM
        Ipm <- double();
        I0pm <- double();
        Ipm<-hgu$Ipm[hgu$Probe.Set.Name == probe_setTmp];
## Error output if the input "probe_set" is not in the spike_in experience
        if(length(Ipm)==0)
            {
            cat(paste("Error: ",probe_setTmp,"cannot be found in spike_in sets\n"));
            continue<-0;}

        I0pm<- 0*Ipm+ hgu$Ipm[hgu$Probe.Set.Name == probe_setTmp&hgu$conc == 0];
        difIpm<-Ipm-I0pm;

## Calculation of I-I0 for MM
        Imm <- double();
        I0mm <- double();
        Imm<-hgu$Imm[hgu$Probe.Set.Name == probe_setTmp];
        I0mm<- 0*Imm+ hgu$Imm[hgu$Probe.Set.Name == probe_setTmp&hgu$conc == 0];
        difImm<-Imm-I0mm;

## idem : je fais une boucle qui prend les probe_sets les uns apres les autres
        gPM <- double();
        gMM <- double();
        gRNA <- double();
        conc <- double();
## recupere les donnees dans le fichier hgu133
        gPM <- hgu$DGpm[hgu$Probe.Set.Name == probe_setTmp];
        gMM <- hgu$DGmm[hgu$Probe.Set.Name == probe_setTmp];
        gRNA <- hgu$DGRNA[hgu$Probe.Set.Name == probe_setTmp];
        conc <- hgu$conc[hgu$Probe.Set.Name == probe_setTmp]*10^(-12);

## compte le nombre de probes
        if (continue==1){
            nbProbes<-min(length(hgu$Ipm[hgu$Probe.Set.Name==probe_setTmp&hgu$conc==0]),
                length(hgu$Ipm[hgu$Probe.Set.Name==probe_setTmp&hgu$conc==8]), 
                length(hgu$Ipm[hgu$Probe.Set.Name==probe_setTmp&hgu$conc==16]));
        expTmp <- integer();
            if(probes[1]=="NULL") {
                expTmp<-seq(1, nbProbes, 1);
            }  ##faire un tableau de la m?me taille que les autres avec des 0 pour les exp?riences ? ne pas afficher.
            else{
                taille<-length(probes);
                for (i in seq(1, nbProbes, 1)){
                    ok=0;
                    for (j in seq(1, taille, 1)){
                        if(probes[j]==i){
                            expTmp<-c(expTmp, i);
                            ok=1;
                        }
                    }
                    if(ok==0){expTmp<-c(expTmp, 0);}
                }
            }


## Calculation of alpha (reduction of target due to hybridization in solution)
        alpha <- (-gRNA+46.5)*(-0.686)
        alpha <- exp(alpha);
        alpha <- alpha+1;
        alpha <- 1/alpha;

## Parameters used (beta=inverse effective temperature, A=maximal intensity value)
        beta<-0.74;
        A<-10000;

## Scaling variables X et X'
        xprimePM <- alpha*conc*exp(beta*gPM);
        xprimeMM <- alpha*conc*exp(beta*gMM);

        xPM <- conc*exp(beta*gPM);
        xMM <- conc*exp(beta*gMM);

## Langmuir isotherm (printed in black in the graph)
        xpPM <- 10^seq(-5, 4, 0.2);
        LangPM <- (10000*xpPM)/(1+xpPM);

## Pour le calcul des indices de dispersion
        langTheoPM <-(10000*xprimePM)/(1+xprimePM);
        langTheoMM <-(10000*xprimeMM)/(1+xprimeMM);

        difIpmTmp<-difIpm
        langTheoPMtmp <- langTheoPM;
        langTheoPM <- langTheoPM [langTheoPM>0 & difIpmTmp>0 & !langTheoPM==1 
                & !difIpmTmp==1 & langTheoPM<10^3 & langTheoPM>10^1] ;
        difIpmTmp <- difIpmTmp [langTheoPMtmp>0 & difIpmTmp>0 & !langTheoPMtmp==1 
                & !difIpmTmp==1 & langTheoPMtmp<10^3 & langTheoPMtmp>10^1];

        difImmTmp<-difImm;
        langTheoMMtmp <- langTheoMM;
        langTheoMM <- langTheoMM [langTheoMM>0 & difImmTmp>0 & !langTheoMM==1 
                & !difImmTmp==1 & langTheoMM<10^3 & langTheoMM>10^1];
        difImmTmp <- difImmTmp [langTheoMMtmp>0 & difImmTmp>0 & !langTheoMMtmp==1 
                & !difImmTmp==1 & langTheoMMtmp<10^3 & langTheoMMtmp>10^1];

        indDispPM <-abs(((log(langTheoPM)-log(difIpmTmp))/log(langTheoPM)));
        indicePM <- mean(indDispPM);

        indDispMM <-abs(((log(langTheoMM)-log(difImmTmp))/log(langTheoMM)));
        indiceMM <- mean(indDispMM);

## 	print(indDispMM);

        indGlobal <- c(indDispPM,indDispMM);
        indiceGlobal <- mean(indGlobal);

        taille80<-.8*length(indGlobal);
        indice80<- double();
        indiceTmp<- sort(indGlobal);

        for (i in seq(1:taille80)){
                indice80<-c(indice80,min(indiceTmp));
                indiceTmp<-indiceTmp[-1];
        }
        indice80<-mean(indice80);

        ##indice du tableau
                kn <- length(xprimePM)/nbProbes;
                txt <- rep(paste(seq(1, nbProbes, 1)), times=kn);
        ##indice experience ? afficher
                exp<-rep(paste(expTmp), times=kn);

## Dimension the graph
                X <- c(10^(-5), 10^4);
                Y <- c(0.1, 10^5);

## Create the graphic window
                plot(X,Y,log="xy",type="n",main="Langmuir + hybrid. sol.",xlab="X'",ylab="I-I0");
## Trace the Langmuir isotherm
                lines(xpPM, LangPM);

## Trace the measured intensities for different probes and different concentrations
                if (!param == "pm" && !param == "mm" && !param == "PM" && !param == "MM"){

## print the PM and MM data
                text(xprimePM[txt==exp & xprimePM>0],difIpm[txt==exp & xprimePM>0], 
                        txt[txt == exp & xprimePM>0], col="blue", pch=1, cex=1.0);
                        text(xprimeMM[txt==exp & xprimePM>0],difImm[txt == exp & xprimePM>0], 
                        txt[txt == exp & xprimePM>0], col="red", pch=1, cex=1.0);

## print the probe set name, the set and PM/MM 
                        text(10^-5,10^{4.7},probe_setTmp,font=2,adj=0);
                        text(10^3,10^-1,fichier,font=2);
                        text(10^2.5,10^2,"PM",col="blue",font=2,adj=0);
                        text(10^2.5,10^1.5,"MM",col="red",font=2,adj=0);

## print the dispersion indices
                        shift<-(min(2,nbGenes)-1)
                        text(10^(-.8+shift),10^.9,paste("disp PM: ",round(indicePM,2)), 
                            col="blueviolet", font=1,adj=0);
                        text(10^(-.8+shift),10^0.4,paste("disp MM: ",round(indiceMM,2)), 
                            col="blueviolet",font=1,adj=0);
                        text(10^(-.8+shift),10^-.1,paste("disp: ",round(indiceGlobal,2)), 
                            col="blueviolet",font=1,adj=0);
                        text(10^(-.8+shift),10^-.6,paste("disp_80%: ",round(indice80,2)), 
                            col="blueviolet",font=1,adj=0);

## If there are at most two probe sets compare the two models:
## Left: Langmuir + hybridization in solution / Right: Langmuir
		    if(nbGenes<3){
                        plot(X,Y,log="xy",type="n",main="Langmuir",xlab="X",ylab="I-I0");
                        ##trace la courbe de Langmuir
                        lines(xpPM, LangPM);
                        text(xPM[txt==exp & xprimePM>0],difIpm[txt==exp & xprimePM>0],
                            txt[txt==exp & xprimePM>0], col="blue", pch=1, cex=1.0);
                        text(xMM[txt==exp & xprimePM>0],difImm[txt==exp & xprimePM>0], 
                            txt[txt==exp & xprimePM>0], col="red", pch=1, cex=1.0);
                        text(10^-5, 10^{4.7}, probe_setTmp, font=2, adj=0);
                        text(10^3, 10^-1, fichier, font=2);
                        text(10^2.5, 10^2, "PM", col="blue", font=2,adj=0);
                        text(10^2.5, 10^1.5,"MM", col="red", font=2,adj=0);
                    }
                }

                if (param == "pm"||param == "PM"){

##afficher les differentes probes PM
                text(xprimePM[xprimePM>0 & txt == exp], difIpm[xprimePM>0 & txt == exp], txt[txt == exp], col="blue", pch=1, cex=1.0);

            ##afficher le nom du probe_set et la legende
                text(10^-5, 10^{4.7}, probe_setTmp, font=2, adj=0);
                text(10^3, 10^-1, fichier, font=2);
                text(10^2.5, 10^2, "PM", col="blue", font=2,adj=0);
                text(10^-.5, 10^.6, paste("disp PM: ",round(indicePM,2)), col="blueviolet", font=1,adj=0);
                if(nbGenes<3){
                    plot(X, Y, log="xy", type="n", main="Langmuir X", xlab="X", ylab="I-I0");
                        ##trace la courbe de Langmuir
                        lines(xpPM, LangPM);
                           text(xPM[txt == exp & xprimePM>0], difIpm[txt == exp & xprimePM>0], txt[txt == exp & xprimePM>0], col="blue", pch=1, cex=1.0);

                        text(10^-5, 10^{4.7}, probe_setTmp, font=2, adj=0);
                        text(10^3, 10^-1, fichier, font=2);
                        text(10^2.5, 10^2, "PM", col="blue", font=2,adj=0);
                    }
                }

                if (param == "mm"||param == "MM"){

                ##afficher les differentes probes MM
                text(xprimeMM[xprimePM>0 & txt == exp], difImm[xprimePM>0 & txt == exp], txt[xprimePM>0 & txt==exp], col="red", pch=1, cex=1.0);

                ##afficher le nom du probe_set et la legende
                text(10^-5, 10^{4.7}, probe_setTmp, font=2, adj=0);
                text(10^3, 10^-1, fichier, font=2);
                text(10^2.5, 10^1.5, "MM", col="red", font=2,adj=0);
                text(10^-.5, 10^0.3, paste("disp MM: ",round(indiceMM,2)), col="blueviolet", font=1,adj=0);
                if(nbGenes<3){
                        plot(X, Y, log="xy", type="n", main="Langmuir X", xlab="X", ylab="I-I0");
                ##trace la courbe de Langmuir
                        lines(xpPM, LangPM);
                        text(xMM[txt == exp & xprimePM>0], difImm[txt == exp & xprimePM>0], txt[txt == exp & xprimePM>0], col="red", pch=1, cex=1.0);
                        text(10^-5, 10^{4.7}, probe_setTmp, font=2, adj=0);
               	        text(10^3, 10^-1, fichier, font=2);
               	        text(10^2.5, 10^1.5, "MM", col="red", font=2,adj=0);
                    }
                }
            }##else{dev.off();}
        }
        if (output=="PS"){
            dev.off();
        }
}
