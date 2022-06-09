#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "structDefinitions.h"
#include "chiralCorrelation.h"

int main(int argc, char const *argv[])
{
	if (argc == 1)
	{
		printf("REQUIRED ARGUMENTS:\n~~~~~~~~~~~~~~~~~~~\n\n[*] argv[0] = ./main\n[*] argv[1] = input dump filename\n[*] argv[2] = input dihedral filename.\n\n");
		exit (1);
	}

	int isXZ_dump = -1, isXZ_dihedral = -1;
	char *pipeString;
	pipeString = (char *) malloc (500 * sizeof (char));

	FILE *inputDump, *inputDihedral, *output;

	// Setting pointer for input dump file
	if (strstr (argv[1], ".xz"))
	{
		snprintf (pipeString, 500, "xzcat %s", argv[1]); inputDump = popen (pipeString, "r"); isXZ_dump = 1;
	}
	else
	{
		inputDump = fopen (argv[1], "r"); isXZ_dump = 0;
	}

	// Setting pointer for input dihedral file
	if (strstr (argv[2], ".xz"))
	{
		snprintf (pipeString, 500 , "xzcat %s", argv[2]); inputDihedral = popen (pipeString, "r"); isXZ_dihedral = 1;
	}
	else
	{
		inputDihedral = fopen (argv[2], "r"); isXZ_dihedral = 0;
	}

	output = fopen ("chiralCorrelation.output", "w");

	int nAtoms = getNAtoms (inputDump), nDihedrals = getNDihedrals (inputDihedral), nTimeframes_dump = countTimeframes_dump (inputDump, nAtoms), nTimeframes_dihedral = countTimeframes_dihedral (inputDihedral, nDihedrals);

	if (strstr (argv[1], ".xz"))
	{
		fclose (inputDump);
		FILE *inputDump;
		snprintf (pipeString, 500, "xzcat %s", argv[1]); 
		inputDump = popen (pipeString, "r"); 
		isXZ_dump = 1;
	}

	if (strstr (argv[2], ".xz"))
	{
		fclose (inputDihedral);
		FILE *inputDihedral;
		snprintf (pipeString, 500 , "xzcat %s", argv[2]); 
		inputDihedral = popen (pipeString, "r"); 
		isXZ_dihedral = 1;
	}

	if (nTimeframes_dump == 0)
	{
		printf("Fault found in input dump file...\n");
		exit (1);
	}
	else printf("nTimeframes_dump: %d\n", nTimeframes_dump);

	if (nTimeframes_dihedral == 0)
	{
		printf("Fault found in input dihedral file...\n");
		exit (1);
	}
	else printf("nTimeframes_dihedral: %d\n", nTimeframes_dihedral);

	// Check the planar density distribution
	computePlanarDensity (inputDump, nTimeframes_dump, nAtoms);

	if (strstr (argv[1], ".xz"))
	{
		fclose (inputDump);
		FILE *inputDump;
		snprintf (pipeString, 500, "xzcat %s", argv[1]); 
		inputDump = popen (pipeString, "r"); 
		isXZ_dump = 1;
	}
	
	if (strstr (argv[2], ".xz"))
	{
		fclose (inputDihedral);
		FILE *inputDihedral;
		snprintf (pipeString, 500 , "xzcat %s", argv[2]); 
		inputDihedral = popen (pipeString, "r"); 
		isXZ_dihedral = 1;
	}

	// Initializing structs to store dump and dihedral values
	// Dihedral values from all timesteps will be stored in a single struct
	DUMP_ENTRIES *dump;
	DIHEDRAL_ENTRIES *dihedral;
	dihedral = (DIHEDRAL_ENTRIES *) malloc (nDihedrals * nTimeframes_dump * sizeof (DIHEDRAL_ENTRIES));

	int currentTimestep_dump, currentTimestep_dihedral;

	// Reading and storing all dihedral data
	rewind (inputDump);
	rewind (inputDihedral);

	int regionNecessary = 1;
	if (nAtoms < 2000) regionNecessary = 0;

	SIM_BOUNDARY regionOfInterest;
	if (regionNecessary == 1)
	{
		printf("\nEnter the region of interest. Only the atoms that fall within this region will be considered for further analysis.\n");
		printf("\nxlo: "); scanf ("%f", &regionOfInterest.xlo); printf("\n");
		printf("xhi: "); scanf ("%f", &regionOfInterest.xhi); printf("\n");
		printf("ylo: "); scanf ("%f", &regionOfInterest.ylo); printf("\n");
		printf("yhi: "); scanf ("%f", &regionOfInterest.yhi); printf("\n");
		printf("zlo: "); scanf ("%f", &regionOfInterest.zlo); printf("\n");
		printf("zhi: "); scanf ("%f", &regionOfInterest.zhi); printf("\n");
	}
	else
	{
		regionOfInterest.xlo = -3000; regionOfInterest.ylo = -3000; regionOfInterest.zlo = -3000;
		regionOfInterest.xhi = 3000; regionOfInterest.yhi = 3000; regionOfInterest.zhi = 3000;
	}

	for (int i = 0; i < nTimeframes_dump; ++i)
	{
		printf("Reading dihedral frame: %d/%d               \r", i, nTimeframes_dump);
		fflush (stdout);
		dump = readDump (inputDump, nAtoms, &currentTimestep_dump);

		readDihedral (&dihedral, inputDihedral, dump, i, nDihedrals, nAtoms, currentTimestep_dump, &currentTimestep_dihedral, regionOfInterest);
	}

	printf("\n");

	CHIRAL_CORRELATION *corr;
	corr = (CHIRAL_CORRELATION *) malloc (nTimeframes_dump * sizeof (CHIRAL_CORRELATION));
	computeCorrelation (dihedral, nDihedrals, nTimeframes_dump, &corr);

	// Printing header information
	fprintf(output, "t gp gm tt tgm tgp gg\n");

	for (int i = 0; i < nTimeframes_dump; ++i)
	{
		fprintf(output, "%f %f %f %f %f %f %f\n", 
			logf (corr[i].correlationTrans), 
			logf (corr[i].correlationGauchePlus), 
			logf (corr[i].correlationGaucheMinus), 
			logf (corr[i].correlationTransTrans), 
			logf (corr[i].correlationTransGaucheMinus), 
			logf (corr[i].correlationTransGauchePlus), 
			logf (corr[i].correlationGaucheGauche));
	}

	free (dump);
	free (dihedral);
	free (corr);
	fclose (inputDump);
	fclose (inputDihedral);
	fclose (output);
	return 0;
}