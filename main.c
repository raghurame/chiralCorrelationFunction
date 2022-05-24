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

	FILE *inputDump, *inputDihedral;
	inputDump = fopen (argv[1], "r");
	inputDihedral = fopen (argv[2], "r");

	int nAtoms = getNAtoms (inputDump), nDihedrals = getNDihedrals (inputDihedral), nTimeframes_dump = countTimeframes_dump (inputDump, nAtoms), nTimeframes_dihedral = countTimeframes_dihedral (inputDihedral, nDihedrals);

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

	// Initializing structs to store dump and dihedral values
	// Dihedral values from all timesteps will be stored in a single struct
	DUMP_ENTRIES *dump;
	DIHEDRAL_ENTRIES *dihedral;
	dihedral = (DIHEDRAL_ENTRIES *) malloc (nDihedrals * nTimeframes_dump * sizeof (DIHEDRAL_ENTRIES));

	int currentTimestep_dump, currentTimestep_dihedral;

	// Reading and storing all dihedral data
	rewind (inputDump);
	rewind (inputDihedral);
	for (int i = 0; i < nTimeframes_dump; ++i)
	{
		printf("Reading dihedral frame: %d/%d               \r", i, nTimeframes_dump);
		fflush (stdout);
		dump = readDump (inputDump, nAtoms, &currentTimestep_dump);
		readDihedral (&dihedral, inputDihedral, dump, i, nDihedrals, nAtoms, currentTimestep_dump, &currentTimestep_dihedral);
	}

	printf("\n");

	CHIRAL_CORRELATION *corr;
	corr = (CHIRAL_CORRELATION *) malloc (nTimeframes_dump * sizeof (CHIRAL_CORRELATION));
	computeCorrelation (dihedral, nDihedrals, nTimeframes_dump, &corr);

	for (int i = 0; i < nTimeframes_dump; ++i)
	{
		printf("%f %f %f %f %f %f %f\n", 
			logf (corr[i].correlationTrans), 
			logf (corr[i].correlationGauchePlus), 
			logf (corr[i].correlationGaucheMinus), 
			logf (corr[i].correlationTransTrans), 
			logf (corr[i].correlationTransGaucheMinus), 
			logf (corr[i].correlationTransGauchePlus), 
			logf (corr[i].correlationGaucheGauche));
		usleep (100000);
	}

	free (dump);
	free (dihedral);
	fclose (inputDump);
	fclose (inputDihedral);
	return 0;
}