#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

int getNAtoms (FILE *input)
{
	int nAtoms;
	char lineString[2000];

	while (fgets (lineString, 2000, input) != NULL)
	{
		if (strstr (lineString, "ITEM: NUMBER OF ATOMS"))
		{
			fgets (lineString, 2000, input);
			sscanf (lineString, "%d\n", &nAtoms);
			break;
		}
	}

	return nAtoms;
}

void saveDump (FILE *input, FILE *output, int nAtoms, int beginningID, int endID)
{
	char lineString[2000];
	int isHeaderline = 0, isAtomline = 0, isNAtomsline = 0, nAtoms_new = (endID - beginningID + 1);
	int atoms_id, atoms_mol, atoms_type, nAtoms_old;
	int nTimesteps = 0;

	while (fgets (lineString, 2000, input) != NULL)
	{
		if (isNAtomsline == 1 && isHeaderline == 1 && isAtomline == 0) {
			fgets (lineString, 2000, input);
			sscanf (lineString, "%d\n", &nAtoms_old);
			fprintf(output, "%d\n", nAtoms_new);
			isNAtomsline = 0; }

		if (strstr (lineString, "ITEM: NUMBER OF ATOMS")) {
			isNAtomsline = 1; }

		if (strstr (lineString, "ITEM: TIMESTEP")) {
			nTimesteps++;
			printf("Current timestep: %d                                            \r", nTimesteps);
			fflush (stdout);
			isHeaderline = 1; 
			isAtomline = 0; }

		if (isHeaderline == 1) {
			fprintf(output, "%s", lineString); }

		if ((isAtomline == 1) && (isHeaderline == 0))
		{
			sscanf (lineString, "%d %d %d\n", &atoms_id, &atoms_mol, &atoms_type);

			if ((atoms_id <= endID) && (atoms_id >= beginningID)) {
				fprintf(output, "%s", lineString); }
		}

		if (strstr (lineString, "ITEM: ATOMS")) {
			isHeaderline = 0; 
			isAtomline = 1; }
	}

	printf("\nRead/Write complete !\n");
}

int main(int argc, char const *argv[])
{
	if (argc == 1)
	{
		printf("\nREQUIRED ARGUMENTS:\n~~~~~~~~~~~~~~~~~~~\n\n[~] argv[0] = program\n[~] argv[1] = input dump file name\n[~] argv[2] = output fixed dump file name\n[~] argv[3] = beginning of atom range to consider\n[~] argv[4] = end of atom range to consider.\n\n");
		exit (1);
	}
	if (argc != 5)
	{
		printf("\nERROR: INSUFFICIENT/INCORRECT ARGUMENTS PASSED !\n~~~~~~\n\nREQUIRED ARGUMENTS:\n~~~~~~~~~~~~~~~~~~~\n\n[~] argv[0] = program\n[~] argv[1] = input dump file name\n[~] argv[2] = output fixed dump file name\n[~] argv[3] = beginning of atom range to consider\n[~] argv[4] = end of atom range to consider.\n\n");
		exit (1);
	}

	FILE *input, *output;
	output = fopen (argv[2], "w");

	char *pipeString;
	pipeString = (char *) malloc (500 * sizeof (char));

	if (strstr (argv[1], ".xz")) {
		snprintf (pipeString, 500, "xzcat %s", argv[1]);
		input = popen (pipeString, "r"); }
	else {
		input = fopen (argv[1], "r"); }

	int nAtoms = getNAtoms (input);

	// rewinding the input file
	if (strstr (argv[1], ".xz")) {
		fclose (input);
		FILE *input;
		input = popen (pipeString, "r"); }
	else {
		rewind (input); }

	saveDump (input, output, nAtoms, atoi (argv[3]), atoi (argv[4]));

	char *compressionString;
	compressionString = (char *) malloc (100 * sizeof (char));
	snprintf (compressionString, 100, "xz %s -ve", argv[2]);

	system (compressionString);

	fclose (input);
	fclose (output);
	return 0;
}