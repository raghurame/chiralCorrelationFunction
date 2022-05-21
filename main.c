#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

typedef struct dumpEntries
{
	int id, molType, atomType, ix, iy, iz;
	float x, y, z;
} DUMP_ENTRIES;

typedef struct dihedralEntries
{
	int id, atom1, atom2, atom3, atom4;
	int isTrans, isGauchePlus, isGaucheMinus, isBackbone;
	float angle;
} DIHEDRAL_ENTRIES;

int getNAtoms (FILE *inputDump)
{
	int nAtoms;
	rewind (inputDump);
	char lineString[2000];

	while (fgets (lineString, 2000, inputDump) != NULL)
	{
		if (strstr (lineString, "ITEM: NUMBER OF ATOMS"))
		{
			fgets (lineString, 2000, inputDump);
			sscanf (lineString, "%d\n", &nAtoms);
			break;
		}
	}

	return nAtoms;
}

int getNDihedrals (FILE *inputDihedral)
{
	int nDihedrals;
	rewind (inputDihedral);
	char lineString[2000];

	while (fgets (lineString, 2000, inputDihedral) != NULL)
	{
		if (strstr (lineString, "ITEM: NUMBER OF ENTRIES"))
		{
			fgets (lineString, 2000, inputDihedral);
			sscanf (lineString, "%d\n", &nDihedrals);
			break;
		}
	}

	return nDihedrals;
}

int countTimeframes_dump (FILE *inputDump, int nAtoms)
{
	int nTimeframes = 0, nLines = 0;
	rewind (inputDump);
	char lineString[2000];

	while (fgets (lineString, 2000, inputDump) != NULL)
	{
		nLines++;

		if (strstr (lineString, "ITEM: TIMESTEP"))
			nTimeframes++;
	}

	if ((nLines / (nAtoms + 9)) == nTimeframes)
		return nTimeframes;
	else
		return 0;

}

int countTimeframes_dihedral (FILE *inputDihedral, int nDihedrals)
{
	int nTimeframes = 0, nLines = 0;
	rewind (inputDihedral);
	char lineString[2000];

	while (fgets (lineString, 2000, inputDihedral) != NULL)
	{
		nLines++;

		if (strstr (lineString, "ITEM: TIMESTEP"))
			nTimeframes++;
	}

	if ((nLines / (nDihedrals + 9)) == nTimeframes)
		return nTimeframes;
	else
		return 0;

	return nTimeframes;
}

DUMP_ENTRIES *readDump (FILE *inputDump, int nAtoms, int *currentTimestep_dump)
{
	DUMP_ENTRIES *dump;
	dump = (DUMP_ENTRIES *) malloc (nAtoms * sizeof (DUMP_ENTRIES));

	char lineString[2000];

	while (fgets (lineString, 2000, inputDump) != NULL)
	{
		if (strstr (lineString, "ITEM: TIMESTEP"))
		{
			fgets (lineString, 2000, inputDump);
			sscanf (lineString, "%d\n", &(*currentTimestep_dump));
		}

		if (strstr (lineString, "ITEM: ATOMS"))
		{
			for (int i = 0; i < nAtoms; ++i)
			{
				fgets (lineString, 2000, inputDump);
				sscanf (lineString, "%d %d %d %f %f %f %*f %*f %*f %d %d %d\n", &dump[i].id, &dump[i].molType, &dump[i].atomType, &dump[i].x, &dump[i].y, &dump[i].z, &dump[i].ix, &dump[i].iy, &dump[i].iz);
			}

			return dump;			
		}
	}

}

int getIndex1d (int i, int j, int width)
{
	// width is the max value of 'j'
	int index1d = 0;
	index1d = (width * i) + j;
	return index1d;
}

int readDihedral (DIHEDRAL_ENTRIES **dihedral, FILE *inputDihedral, DUMP_ENTRIES *dump, int currentStep, int nDihedrals, int currentTimestep_dump, int *currentTimestep_dihedral)
{
	char lineString[2000];
	int dihAtom1, dihAtom2, dihAtom3, dihAtom4, arrayIndex;

	while (fgets (lineString, 2000, inputDihedral) != NULL)
	{
		if (strstr (lineString, "ITEM: TIMESTEP"))
		{
			fgets (lineString, 2000, inputDihedral);
			sscanf (lineString, "%d", &(*currentTimestep_dihedral));
		}

		if (((*currentTimestep_dihedral) == currentTimestep_dump) && (strstr (lineString, "ITEM: ENTRIES")))
		{
			for (int i = 0; i < nDihedrals; ++i)
			{
				arrayIndex = getIndex1d (currentStep, i, nDihedrals);

				fgets (lineString, 2000, inputDihedral);
				sscanf (lineString, "%d %f %d %d %d %d\n", &(*dihedral)[i].id, &(*dihedral)[i].angle, &(*dihedral)[i].atom1, &(*dihedral)[i].atom2, &(*dihedral)[i].atom3, &(*dihedral)[i].atom4);

				dihAtom1 = (*dihedral)[i].atom1;
				dihAtom2 = (*dihedral)[i].atom2;
				dihAtom3 = (*dihedral)[i].atom3;
				dihAtom4 = (*dihedral)[i].atom4;

				if ((dihAtom1 != 3) && (dihAtom2 != 3) && (dihAtom3 != 3) && (dihAtom4 != 3)) (*dihedral)[i].isBackbone = 1;
				else (*dihedral)[i].isBackbone = 0;

				// Resetting the values before assigning '1'
				(*dihedral)[i].isGaucheMinus = 0; (*dihedral)[i].isGauchePlus = 0; (*dihedral)[i].isTrans = 0;

				if (((*dihedral)[i].angle < -38) && ((*dihedral)[i].angle > -116)) (*dihedral)[i].isGaucheMinus = 1;
				else if (((*dihedral)[i].angle > 34) && ((*dihedral)[i].angle < 115)) (*dihedral)[i].isGauchePlus = 1;
				else if ((((*dihedral)[i].angle < -117) && ((*dihedral)[i].angle > -180)) || (((*dihedral)[i].angle < 180) && ((*dihedral)[i].angle > 116))) (*dihedral)[i].isTrans = 1;
			}

			return 1;
		}
	}
}

void computeCorrelation (DIHEDRAL_ENTRIES *dihedral, int nDihedrals, int nTimeframes_dump, int **lHandCorrelation, int **rHandCorrelation, int **exCoilCorrelation, int **coilCorrelation)
{
	int arrayIndex1, arrayIndex2;

	for (int currentDihedral = 0; currentDihedral < nDihedrals; ++currentDihedral)
	{
		// Time lag
		for (int lag = 0; lag < nTimeframes_dump; ++lag)
		{
			// Going through all the elements
			for (int currentElement = 0; currentElement < (nTimeframes_dump - lag); ++currentElement)
			{
				arrayIndex1 = getIndex1d (currentElement, currentDihedral, nDihedrals);
				arrayIndex2 = getIndex1d (currentElement + lag, currentDihedral, nDihedrals);
				printf("lag: %d; %d[%d] %d[%d]\n", lag, dihedral[arrayIndex1].id, arrayIndex1, dihedral[arrayIndex2].id, arrayIndex2);
				sleep (1);
			}
		}		
	}
}

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
		dump = readDump (inputDump, nAtoms, &currentTimestep_dump);
		// TO DO:
		// Values are currently stored only for i = 0.
		// Fix the bug and save the values for all timeframes
		readDihedral (&dihedral, inputDihedral, dump, i, nDihedrals, currentTimestep_dump, &currentTimestep_dihedral);
	}

	// Calculating correlation function
	int *lHandCorrelation, *rHandCorrelation, *exCoilCorrelation, *coilCorrelation;
	lHandCorrelation = (int *) malloc (nTimeframes_dump * sizeof (int));
	rHandCorrelation = (int *) malloc (nTimeframes_dump * sizeof (int));
	exCoilCorrelation = (int *) malloc (nTimeframes_dump * sizeof (int));
	coilCorrelation = (int *) malloc (nTimeframes_dump * sizeof (int));

	// TO DO:
	// Once the above todo is complete, then continue with this function
	// Compute autocorrelation of the timeseries for every dihedral
	// Calculate the average correlation function for every conformation
	// Later, implement the same for melt by taking distance from the crystal substrate into consideration
	computeCorrelation (dihedral, nDihedrals, nTimeframes_dump, &lHandCorrelation, &rHandCorrelation, &exCoilCorrelation, &coilCorrelation);

	free (dump);
	free (dihedral);
	fclose (inputDump);
	fclose (inputDihedral);
	return 0;
}