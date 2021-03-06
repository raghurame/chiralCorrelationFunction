#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "structDefinitions.h"
#include "chiralCorrelation.h"

void computePlanarDensity (FILE *inputDump, int nTimeframes_dump, int nAtoms)
{
	rewind (inputDump);

	char lineString[2000];
	SIM_BOUNDARY bounds;
	bounds.xlo = 0; bounds.ylo = 0; bounds.zlo = 0; bounds.xhi = 0; bounds.yhi = 0; bounds.zhi = 0; bounds.xLength = 0; bounds.yLength = 0; bounds.zLength = 0;
	int *xBins, *yBins, *zBins, nXBins = 0, nYBins = 0, nZBins = 0, atomType;
	float xBinLength = 3, yBinLength = 3, zBinLength = 3, xCoord, yCoord, zCoord;

	printf("Entered computePlanarDensity function...\n");

	for (int i = 0; i < nTimeframes_dump; ++i)
	{
		for (int j = 0; j < 9; ++j)
		{
			fgets (lineString, 2000, inputDump);
			if (j == 5 && bounds.xlo == 0 && bounds.xhi == 0) { sscanf (lineString, "%f %f\n", &bounds.xlo, &bounds.xhi); }
			if (j == 6 && bounds.ylo == 0 && bounds.yhi == 0) { sscanf (lineString, "%f %f\n", &bounds.ylo, &bounds.yhi); }
			if (j == 7 && bounds.zlo == 0 && bounds.zhi == 0) { sscanf (lineString, "%f %f\n", &bounds.zlo, &bounds.zhi); }
		}

		for (int j = 0; j < nAtoms; ++j)
		{
			if ((j == 0) && (nXBins == 0) && (nYBins == 0) && (nZBins == 0))
			{
				printf("xlo %f\txhi %f\nylo %f\tyhi %f\nzlo %f\tzhi %f\n", bounds.xlo, bounds.xhi, bounds.ylo, bounds.yhi, bounds.zlo, bounds.zhi);
				bounds.xLength = bounds.xhi - bounds.xlo;
				bounds.yLength = bounds.yhi - bounds.ylo;
				bounds.zLength = bounds.zhi - bounds.zlo;

				nXBins = (int) ceil (bounds.xLength / xBinLength);
				nYBins = (int) ceil (bounds.yLength / yBinLength);
				nZBins = (int) ceil (bounds.zLength / zBinLength);

				printf(" nXBins: %d\n nYBins: %d\n nZBins: %d\n", nXBins, nYBins, nZBins);

				xBins = (int *) calloc (nXBins, sizeof (int));
				yBins = (int *) calloc (nYBins, sizeof (int));
				zBins = (int *) calloc (nZBins, sizeof (int));
			}

			fgets (lineString, 2000, inputDump);
			sscanf (lineString, "%*d %*d %d %f %f %f\n", &atomType, &xCoord, &yCoord, &zCoord);

			if (atomType == 1 || atomType == 2)
			{
				for (int k = 0; k < nXBins; ++k)
				{
					if (xCoord > (bounds.xlo + (k * xBinLength)) && xCoord <= (bounds.xlo + ((k + 1) * xBinLength)))
					{ 
						xBins[k]++; 
						break;
					}
				}

				for (int k = 0; k < nYBins; ++k)
				{
					if (yCoord > (bounds.ylo + (k * yBinLength)) && yCoord <= (bounds.ylo + ((k + 1) * yBinLength)))
					{ 
						yBins[k]++; 
						break; 
					}
				}

				for (int k = 0; k < nZBins; ++k)
				{
					if (zCoord > (bounds.zlo + (k * zBinLength)) && zCoord <= (bounds.zlo + ((k + 1) * zBinLength))) 
					{ 
						zBins[k]++; 
						break; 
					}
				}
			}
		}

		printf("Calculating density distribution: %d/%d           \r", i, nTimeframes_dump);
		fflush (stdout);
	}

	rewind (inputDump);

	FILE *xDistributionOutput, *yDistributionOutput, *zDistributionOutput;
	xDistributionOutput = fopen ("xDistribution.output", "w");
	yDistributionOutput = fopen ("yDistribution.output", "w");
	zDistributionOutput = fopen ("zDistribution.output", "w");

	int xBins_max = 0, yBins_max = 0, zBins_max = 0;

	for (int i = 0; i < nXBins; ++i) { if (xBins[i] > xBins_max) xBins_max = xBins[i]; }
	for (int i = 0; i < nYBins; ++i) { if (yBins[i] > yBins_max) yBins_max = yBins[i]; }
	for (int i = 0; i < nZBins; ++i) { if (zBins[i] > zBins_max) zBins_max = zBins[i]; }

	for (int i = 0; i < nXBins; ++i) fprintf(xDistributionOutput, "%f %d\n", bounds.xlo + (i * xBinLength), xBins[i]);
	for (int i = 0; i < nYBins; ++i) fprintf(yDistributionOutput, "%f %d\n", bounds.ylo + (i * yBinLength), yBins[i]);
	for (int i = 0; i < nZBins; ++i) fprintf(zDistributionOutput, "%f %d\n", bounds.zlo + (i * zBinLength), zBins[i]);

	fclose (xDistributionOutput);
	fclose (yDistributionOutput);
	fclose (zDistributionOutput);
}

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
		{
			nTimeframes++;
			// printf("Counting dump timeframes: %d        \r", nTimeframes);
			// fflush (stdout);			
		}
	}

	printf("\n");

	// if ((nLines / (nAtoms + 9)) == nTimeframes)
		return nTimeframes;
	// else
	// 	return 0;
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
		{
			nTimeframes++;
			// printf("Counting dihedral timeframes: %d         \r", nTimeframes);
			// fflush (stdout);
		}
	}

	// printf("\n");

	// if ((nLines / (nDihedrals + 9)) == nTimeframes)
		return nTimeframes;
	// else
	// 	return 0;
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

int readDihedral (DIHEDRAL_ENTRIES **dihedral, FILE *inputDihedral, DUMP_ENTRIES *dump, int currentStep, int nDihedrals, int nAtoms, int currentTimestep_dump, int *currentTimestep_dihedral, SIM_BOUNDARY regionOfInterest)
{
	char lineString[2000];
	int dihAtom1, dihAtom2, dihAtom3, dihAtom4, dihAtom1Type, dihAtom2Type, dihAtom3Type, dihAtom4Type, arrayIndex, arrayIndex1;
	DIHEDRAL_ENTRIES *dihTemp;
	dihTemp = (DIHEDRAL_ENTRIES *) malloc (nDihedrals * sizeof (DIHEDRAL_ENTRIES));

	COORDINATES atom1coords, atom2coords, atom3coords, atom4coords;

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
				fgets (lineString, 2000, inputDihedral);
				sscanf (lineString, "%d %f %d %d %d %d\n", 
					&dihTemp[i].id, 
					&dihTemp[i].angle, 
					&dihTemp[i].atom1, 
					&dihTemp[i].atom2, 
					&dihTemp[i].atom3, 
					&dihTemp[i].atom4);
			}

			goto processFurther;
		}
	}

	processFurther: ;

	// Re-arranging the dihedral entries in ascending order
	int currentDihID = 0;
	for (int i = 0; i < nAtoms; ++i)
	{
		for (int j = 0; j < nDihedrals; ++j)
		{
			if (dihTemp[j].atom2 == (i + 1))
			{
				arrayIndex = getIndex1d (currentStep, currentDihID, nDihedrals);

				(*dihedral)[arrayIndex].id = dihTemp[j].id;
				(*dihedral)[arrayIndex].angle = dihTemp[j].angle;
				(*dihedral)[arrayIndex].atom1 = dihTemp[j].atom1;
				(*dihedral)[arrayIndex].atom2 = dihTemp[j].atom2;
				(*dihedral)[arrayIndex].atom3 = dihTemp[j].atom3;
				(*dihedral)[arrayIndex].atom4 = dihTemp[j].atom4;
				currentDihID++;
			}
		}
	}

	// Assigning isBackbone, isTrans, isGauchePlus, isGaucheMinus
	for (int i = 0; i < nDihedrals; ++i)
	{
		arrayIndex = getIndex1d (currentStep, i, nDihedrals);

		dihAtom1 = (*dihedral)[arrayIndex].atom1;
		dihAtom2 = (*dihedral)[arrayIndex].atom2;
		dihAtom3 = (*dihedral)[arrayIndex].atom3;
		dihAtom4 = (*dihedral)[arrayIndex].atom4;

		atom1coords.x = dump[dihAtom1 - 1].x;
		atom1coords.y = dump[dihAtom1 - 1].y;
		atom1coords.z = dump[dihAtom1 - 1].z;
		atom2coords.x = dump[dihAtom2 - 1].x;
		atom2coords.y = dump[dihAtom2 - 1].y;
		atom2coords.z = dump[dihAtom2 - 1].z;
		atom3coords.x = dump[dihAtom3 - 1].x;
		atom3coords.y = dump[dihAtom3 - 1].y;
		atom3coords.z = dump[dihAtom3 - 1].z;
		atom4coords.x = dump[dihAtom4 - 1].x;
		atom4coords.y = dump[dihAtom4 - 1].y;
		atom4coords.z = dump[dihAtom4 - 1].z;

		dihAtom1Type = dump[dihAtom1 - 1].atomType;
		dihAtom2Type = dump[dihAtom2 - 1].atomType;
		dihAtom3Type = dump[dihAtom3 - 1].atomType;
		dihAtom4Type = dump[dihAtom4 - 1].atomType;

		if ((dihAtom1Type != 3) && (dihAtom2Type != 3) && (dihAtom3Type != 3) && (dihAtom4Type != 3)) (*dihedral)[arrayIndex].isBackbone = 1;
		else (*dihedral)[arrayIndex].isBackbone = 0;

		// Resetting all of the values before assigning '1'
		(*dihedral)[arrayIndex].isGaucheMinus = 0; 
		(*dihedral)[arrayIndex].isGauchePlus = 0; 
		(*dihedral)[arrayIndex].isTrans = 0; 
		(*dihedral)[arrayIndex].isTransTrans_previous = 0; 
		(*dihedral)[arrayIndex].isTransGauchePlus_previous = 0; 
		(*dihedral)[arrayIndex].isTransGaucheMinus_previous = 0; 
		(*dihedral)[arrayIndex].isGaucheGauche_previous = 0; 
		(*dihedral)[arrayIndex].isTransTrans_next = 0; 
		(*dihedral)[arrayIndex].isTransGauchePlus_next = 0; 
		(*dihedral)[arrayIndex].isTransGaucheMinus_next = 0; 
		(*dihedral)[arrayIndex].isGaucheGauche_next = 0;

		if (
			((*dihedral)[arrayIndex].angle < -38) && 
			((*dihedral)[arrayIndex].angle > -116) && 
			(atom2coords.x > regionOfInterest.xlo) && 
			(atom2coords.x < regionOfInterest.xhi) && 
			(atom2coords.y > regionOfInterest.ylo) && 
			(atom2coords.y < regionOfInterest.yhi) && 
			(atom2coords.z > regionOfInterest.zlo) && 
			(atom2coords.z < regionOfInterest.zhi) && 
			(atom3coords.x > regionOfInterest.xlo) && 
			(atom3coords.x < regionOfInterest.xhi) && 
			(atom3coords.y > regionOfInterest.ylo) && 
			(atom3coords.y < regionOfInterest.yhi) && 
			(atom3coords.z > regionOfInterest.zlo) && 
			(atom3coords.z < regionOfInterest.zhi)
			)
		{
			(*dihedral)[arrayIndex].isGaucheMinus = 1;
		}
		else if (
			((*dihedral)[arrayIndex].angle > 34) && 
			((*dihedral)[arrayIndex].angle < 115) && 
			(atom2coords.x > regionOfInterest.xlo) && 
			(atom2coords.x < regionOfInterest.xhi) && 
			(atom2coords.y > regionOfInterest.ylo) && 
			(atom2coords.y < regionOfInterest.yhi) && 
			(atom2coords.z > regionOfInterest.zlo) && 
			(atom2coords.z < regionOfInterest.zhi) && 
			(atom3coords.x > regionOfInterest.xlo) && 
			(atom3coords.x < regionOfInterest.xhi) && 
			(atom3coords.y > regionOfInterest.ylo) && 
			(atom3coords.y < regionOfInterest.yhi) && 
			(atom3coords.z > regionOfInterest.zlo) && 
			(atom3coords.z < regionOfInterest.zhi)
			)
		{
			(*dihedral)[arrayIndex].isGauchePlus = 1;
		}
		else if (
			((((*dihedral)[arrayIndex].angle < -117) && ((*dihedral)[arrayIndex].angle > -180)) || 
			(((*dihedral)[arrayIndex].angle < 180) && ((*dihedral)[arrayIndex].angle > 116))) && 
			(atom2coords.x > regionOfInterest.xlo) && 
			(atom2coords.x < regionOfInterest.xhi) && 
			(atom2coords.y > regionOfInterest.ylo) && 
			(atom2coords.y < regionOfInterest.yhi) && 
			(atom2coords.z > regionOfInterest.zlo) && 
			(atom2coords.z < regionOfInterest.zhi) && 
			(atom3coords.x > regionOfInterest.xlo) && 
			(atom3coords.x < regionOfInterest.xhi) && 
			(atom3coords.y > regionOfInterest.ylo) && 
			(atom3coords.y < regionOfInterest.yhi) && 
			(atom3coords.z > regionOfInterest.zlo) && 
			(atom3coords.z < regionOfInterest.zhi)
			)
		{
			(*dihedral)[arrayIndex].isTrans = 1;
		}
	}

	// Assigning isTransTrans, isTransGauchePlus, isTransGaucheMinus, isGaucheGauche
	for (int i = 0; i < nDihedrals; ++i)
	{
		arrayIndex = getIndex1d (currentStep, i, nDihedrals);

		// Search for dihedral status of neighboring bond from chain backbone
		if ((*dihedral)[arrayIndex].isBackbone == 0)
		{
			(*dihedral)[arrayIndex].isTransTrans_previous = 0; (*dihedral)[arrayIndex].isTransGauchePlus_previous = 0; (*dihedral)[arrayIndex].isTransGaucheMinus_previous = 0; (*dihedral)[arrayIndex].isGaucheGauche_previous = 0;
			(*dihedral)[arrayIndex].isTransTrans_next = 0; (*dihedral)[arrayIndex].isTransGauchePlus_next = 0; (*dihedral)[arrayIndex].isTransGaucheMinus_next = 0; (*dihedral)[arrayIndex].isGaucheGauche_next = 0;
		}
		else if ((*dihedral)[arrayIndex].isBackbone == 1)
		{
			for (int j = 0; j < nDihedrals; ++j)
			{
				arrayIndex1 = getIndex1d (currentStep, j, nDihedrals);

				// Checking for dihedral of neighboring backbone bond
				if (
					((*dihedral)[arrayIndex].atom2 == (*dihedral)[arrayIndex1].atom1 && (*dihedral)[arrayIndex].atom3 == (*dihedral)[arrayIndex1].atom2) 
					|| 
					((*dihedral)[arrayIndex].atom2 == (*dihedral)[arrayIndex1].atom3 && (*dihedral)[arrayIndex].atom3 == (*dihedral)[arrayIndex1].atom4)
					)
				{
					// Checking for TT conformation
					if ((*dihedral)[arrayIndex].isTrans == 1 && (*dihedral)[arrayIndex1].isTrans == 1 && arrayIndex1 > arrayIndex)
					{
						(*dihedral)[arrayIndex].isTransTrans_next = 1; 
						(*dihedral)[arrayIndex1].isTransTrans_previous = 1;
					}
					if ((*dihedral)[arrayIndex].isTrans == 1 && (*dihedral)[arrayIndex1].isTrans == 1 && arrayIndex1 < arrayIndex)
					{
						(*dihedral)[arrayIndex].isTransTrans_previous = 1; 
						(*dihedral)[arrayIndex1].isTransTrans_next = 1;
					}

					// Checking for TG+ conformation
					if ((*dihedral)[arrayIndex].isTrans == 1 && (*dihedral)[arrayIndex1].isGauchePlus == 1 && arrayIndex1 > arrayIndex)
					{
						(*dihedral)[arrayIndex].isTransGauchePlus_next = 1; 
						(*dihedral)[arrayIndex1].isTransGauchePlus_previous = 1;
					}
					if ((*dihedral)[arrayIndex].isTrans == 1 && (*dihedral)[arrayIndex1].isGauchePlus == 1 && arrayIndex1 < arrayIndex)
					{
						(*dihedral)[arrayIndex].isTransGauchePlus_previous = 1; 
						(*dihedral)[arrayIndex1].isTransGauchePlus_next = 1;
					}

					// Checking for G+T conformation
					if ((*dihedral)[arrayIndex].isGauchePlus == 1 && (*dihedral)[arrayIndex1].isTrans == 1 && arrayIndex1 > arrayIndex)
					{
						(*dihedral)[arrayIndex].isTransGauchePlus_next = 1; 
						(*dihedral)[arrayIndex1].isTransGauchePlus_previous = 1;
					}
					if ((*dihedral)[arrayIndex].isGauchePlus == 1 && (*dihedral)[arrayIndex1].isTrans == 1 && arrayIndex1 < arrayIndex)
					{
						(*dihedral)[arrayIndex].isTransGauchePlus_previous = 1; 
						(*dihedral)[arrayIndex1].isTransGauchePlus_next = 1;
					}

					// Checking for TG- conformation
					if ((*dihedral)[arrayIndex].isTrans == 1 && (*dihedral)[arrayIndex1].isGaucheMinus == 1 && arrayIndex1 > arrayIndex)
					{
						(*dihedral)[arrayIndex].isTransGaucheMinus_next = 1; 
						(*dihedral)[arrayIndex1].isTransGaucheMinus_previous = 1;
					}
					if ((*dihedral)[arrayIndex].isTrans == 1 && (*dihedral)[arrayIndex1].isGaucheMinus == 1 && arrayIndex1 < arrayIndex)
					{
						(*dihedral)[arrayIndex].isTransGaucheMinus_previous = 1; 
						(*dihedral)[arrayIndex1].isTransGaucheMinus_next = 1;
					}

					// Checking for G-T conformation
					if ((*dihedral)[arrayIndex].isGaucheMinus == 1 && (*dihedral)[arrayIndex1].isTrans == 1 && arrayIndex1 > arrayIndex)
					{
						(*dihedral)[arrayIndex].isTransGaucheMinus_next = 1; 
						(*dihedral)[arrayIndex1].isTransGaucheMinus_previous = 1;
					}
					if ((*dihedral)[arrayIndex].isGaucheMinus == 1 && (*dihedral)[arrayIndex1].isTrans == 1 && arrayIndex1 < arrayIndex)
					{
						(*dihedral)[arrayIndex].isTransGaucheMinus_previous = 1; 
						(*dihedral)[arrayIndex1].isTransGaucheMinus_next = 1;
					}

					// Checking for G(+/-)G(+/-) conformation
					if (
						((*dihedral)[arrayIndex].isGauchePlus == 1 && (*dihedral)[arrayIndex1].isGauchePlus == 1) ||
						((*dihedral)[arrayIndex].isGauchePlus == 1 && (*dihedral)[arrayIndex1].isGaucheMinus == 1) ||
						((*dihedral)[arrayIndex].isGaucheMinus == 1 && (*dihedral)[arrayIndex1].isGauchePlus == 1) ||
						((*dihedral)[arrayIndex].isGaucheMinus == 1 && (*dihedral)[arrayIndex1].isGaucheMinus == 1)
						)
					{
						if (arrayIndex1 > arrayIndex)
						{
							(*dihedral)[arrayIndex].isGaucheGauche_next = 1; 
							(*dihedral)[arrayIndex1].isGaucheGauche_previous = 1;
						}
						else if (arrayIndex1 < arrayIndex)
						{
							(*dihedral)[arrayIndex].isGaucheGauche_previous = 1; 
							(*dihedral)[arrayIndex1].isGaucheGauche_next = 1;
						}
					}
				}
			}
		}
	}

	free (dihTemp);
}

void computeCorrelation (DIHEDRAL_ENTRIES *dihedral, int nDihedrals, int nTimeframes_dump, CHIRAL_CORRELATION **corr)
{
	int arrayIndex1, arrayIndex2, *nItems, chirality1, chirality2;
	nItems = (int *) calloc (nTimeframes_dump, sizeof (int));

	for (int i = 0; i < nTimeframes_dump; ++i)
	{
		(*corr)[i].correlationTrans = 0; (*corr)[i].correlationGauchePlus = 0; (*corr)[i].correlationGaucheMinus = 0; (*corr)[i].correlationTransTrans = 0; (*corr)[i].correlationTransGaucheMinus = 0; (*corr)[i].correlationTransGauchePlus = 0; (*corr)[i].correlationGaucheGauche = 0;
	}

	// Autocorrelation will be calculated for all dihedrals
	for (int currentDihedral = 0; currentDihedral < nDihedrals; ++currentDihedral)
	{
		printf("Computing correlation for dihedrals (%d/%d)...          \r", currentDihedral, nDihedrals);
		fflush (stdout);
		// Time lag
		for (int lag = 0; lag < nTimeframes_dump; ++lag)
		{
			// Going through all the elements
			for (int currentElement = 0; currentElement < (nTimeframes_dump - lag); ++currentElement)
			{
				arrayIndex1 = getIndex1d (currentElement, currentDihedral, nDihedrals);
				arrayIndex2 = getIndex1d (currentElement + lag, currentDihedral, nDihedrals);
				(*corr)[lag].correlationTrans += (float)dihedral[arrayIndex1].isTrans * (float)dihedral[arrayIndex2].isTrans;
				(*corr)[lag].correlationGauchePlus += (float)dihedral[arrayIndex1].isGauchePlus * (float)dihedral[arrayIndex2].isGauchePlus;
				(*corr)[lag].correlationGaucheMinus += (float)dihedral[arrayIndex1].isGaucheMinus * (float)dihedral[arrayIndex2].isGaucheMinus;

				chirality1 = 0; chirality2 = 0;
				if (dihedral[arrayIndex1].isTransTrans_previous == 1 || dihedral[arrayIndex1].isTransTrans_next == 1) chirality1 = 1;
				if (dihedral[arrayIndex2].isTransTrans_previous == 1 || dihedral[arrayIndex2].isTransTrans_next == 1) chirality2 = 1;
				(*corr)[lag].correlationTransTrans += (float)chirality1 * (float)chirality2;

				chirality1 = 0; chirality2 = 0;
				if (dihedral[arrayIndex1].isTransGaucheMinus_previous == 1 || dihedral[arrayIndex1].isTransGaucheMinus_next == 1) chirality1 = 1;
				if (dihedral[arrayIndex2].isTransGaucheMinus_previous == 1 || dihedral[arrayIndex2].isTransGaucheMinus_next == 1) chirality2 = 1;
				(*corr)[lag].correlationTransGaucheMinus += (float)chirality1 * (float)chirality2;

				chirality1 = 0; chirality2 = 0;
				if (dihedral[arrayIndex1].isTransGauchePlus_previous == 1 || dihedral[arrayIndex1].isTransGauchePlus_next == 1) chirality1 = 1;
				if (dihedral[arrayIndex2].isTransGauchePlus_previous == 1 || dihedral[arrayIndex2].isTransGauchePlus_next == 1) chirality2 = 1;
				(*corr)[lag].correlationTransGauchePlus += (float)chirality1 * (float)chirality2;

				chirality1 = 0; chirality2 = 0;
				if (dihedral[arrayIndex1].isGaucheGauche_previous == 1 || dihedral[arrayIndex1].isGaucheGauche_next == 1) chirality1 = 1;
				if (dihedral[arrayIndex2].isGaucheGauche_previous == 1 || dihedral[arrayIndex2].isGaucheGauche_next == 1) chirality2 = 1;
				(*corr)[lag].correlationGaucheGauche += (float)chirality1 * (float)chirality2;

				nItems[lag]++;
			}
		}		
	}

	printf("\n");

	for (int i = 0; i < nTimeframes_dump; ++i)
	{
		(*corr)[i].correlationTrans /= nItems[i];
		(*corr)[i].correlationTrans /= (*corr)[0].correlationTrans;

		(*corr)[i].correlationGauchePlus /= nItems[i];
		(*corr)[i].correlationGauchePlus /= (*corr)[0].correlationGauchePlus;

		(*corr)[i].correlationGaucheMinus /= nItems[i];
		(*corr)[i].correlationGaucheMinus /= (*corr)[0].correlationGaucheMinus;

		(*corr)[i].correlationTransTrans /= nItems[i];
		(*corr)[i].correlationTransTrans /= (*corr)[0].correlationTransTrans;

		(*corr)[i].correlationTransGauchePlus /= nItems[i];
		(*corr)[i].correlationTransGauchePlus /= (*corr)[0].correlationTransGauchePlus;

		(*corr)[i].correlationTransGaucheMinus /= nItems[i];
		(*corr)[i].correlationTransGaucheMinus /= (*corr)[0].correlationTransGaucheMinus;

		(*corr)[i].correlationGaucheGauche /= nItems[i];
		(*corr)[i].correlationGaucheGauche /= (*corr)[0].correlationGaucheGauche;

		// printf("%f\n", corr[i].correlationTrans);
		// sleep (1);
	}
}