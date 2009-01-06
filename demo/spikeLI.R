## A quick demo to see the package functions
##         -- Paul

cont <- function()
	{ ssw<-readline("\n\nHit return to continue") }

refresh <- function(gene){

conc133<-c(0.,0.125*2^(0:12));

for (k in conc133[2:length(conc133)] ) 
    { 	readline("\n\nPress any key to get IvsDG with higher concentration");
        IvsDG(gene, k); }
}


spike_in <-function()
{
	cat("\n\n\n\n\n\n\t\t-------------------------------------")
	cat("\n\t\t Welcome to the SpikeLI package demo ")
	cat("\n\t\t-------------------------------------")
	cat("\n\nFor convenience the probe sets of the Affymetrix spike-in")
	cat("\n(or Latin square) experiment are stored in some vectors\n")
	cat("\nSPIKE_IN contains all the probe sets of the HGU133 experiment\n") 
	cat("\n--------------------------------------------------------------\n") 
}

spike_inh<-function()
{
	cat("\n\n\n-----------------------------------------------------------------------------") 
	cat("\nSPIKE_INH contains all the *human* probe sets of the HGU133 experiment") 
	cat("\nSPIKE_INB contains all the *bacterial* probe sets of the HGU133 experiment") 
	cat("\nSPIKE_INA contains all the *artificial* probe sets of the HGU133 experiment\n")
	cat("\n-----------------------------------------------------------------------------\n\n\n\n") 
}

spike_intro<-function()
{
	cat("\n\n\n\n\n\n\-----------------------------------------------------------------------------") 
	cat("\nSpikeLI: functions")
	cat("\n-----------------------------------------------------------------------------") 
	cat("\nThe SpikeLI package uses three functions based on the Lamgmuir isotherm: Ivsc, IvsDG, collapse") 
	cat("\nIvsc: Plot Intensity versus concentration") 
	cat("\nIvsDG: Plot Intensity vs. free energy, for all probes at a chosen concentration.") 
	cat("\ncollapse: Plot intensity versus all concentrations for all probes vs. free energy\n") 
	cat("\nthe demo also shows gradual increase of the concentration on the IvsDG function\n") 
	cat("\n-----------------------------------------------------------------------------\n\n\n\n") 
}


spike_inIc<-function()
{
	cat("\n-----------------------------------------------------------------------------") 
	cat("\nIvsc: Intensity vs. concentration")
	cat("\n-----------------------------------------------------------------------------")  
	cat("\nPlot the measured Intensities of the spike-in concentration for a given probe.") 
	cat("\nThis takes into account the Perfect Matched value and the Mis-Matched values.") 
	cat("\n-----------------------------------------------------------------------------\n\n\n\n") 
}


spike_inIDG<-function()
{
	cat("\n\n\n\n\n\n\-----------------------------------------------------------------------------") 
	cat("\nIvsDG: Intensity vs. free energy")
	cat("\n-----------------------------------------------------------------------------") 
	cat("\nThis plots measured intensities at a given spike-in concentration as a function of the hybridization free energies") 
	cat("\nThe first graph is the intensity for each probe. The second will plot this intensity against the free energies.");
	cat("\nThe solid line of the second graph is the Langmuir isotherm. Both graphs take into account the MisMatch and Perfect Match values") 
	cat("\n\n") 
	cat("\n-----------------------------------------------------------------------------\n\n\n\n") 
}


spike_incoll<-function()
{
	cat("\n\n\n\n\n\n\-----------------------------------------------------------------------------") 
	cat("\ncollapse")
	cat("\n-----------------------------------------------------------------------------") 
	cat("\nIt allows to plot on a same graph all spike-in data at different concentrations of a given probe set") 
	cat("\nUp to four probe sets can be visualised simultaneously") 
	cat("\nProbes strongly deviating from the curve should be considered as outliers\n") 
	cat("\n-----------------------------------------------------------------------------\n\n\n\n") 
}

spike_inprog<-function()
{
	cat("\n-----------------------------------------------------------------------------")
	cat("\nGradual increase of concentration on IvsDG function")
	cat("\n-----------------------------------------------------------------------------")  
	cat("\nWhile pressing the space bar, you will gradually increase the intensity of the probe set") 
	cat("\naccording to the value in the latin square matrix of the microArray data of the probe set") 
	cat("\nThis gives us a better idea of the Intensity at different concentrations and helps us spot outliers") 
	cat("\nBy plotting the intensity against the free energy functions we can see which probe sets are iregular at different values") 
	cat("\n-----------------------------------------------------------------------------\n\n\n\n") 
}

spike_in(); SPIKE_IN; cont();
spike_inh(); SPIKE_INH; SPIKE_INB; SPIKE_INA; cont();

# readline("\n\nPress any key to see Ivsc function");

spike_intro()

spike_inIc()
Ivsc(SPIKE_IN[3]);

# readline("\n\nPress any key to see IvsDG function");
cont();


spike_inIDG()
IvsDG(SPIKE_IN[3], 512);

# readline("\n\nPress any key to see collapse function");
cont();


spike_incoll()

collapse(SPIKE_IN[3]);
cont();

spike_inprog()

# readline("\n\nPress any key to see the variation of the concentration on the IvsDG function");
cont();

refresh(SPIKE_IN[3]);
