#ifndef CHIRALCORRELATION_H
#define CHIRALCORRELATION_H

int getNAtoms (FILE *inputDump);
int getNDihedrals (FILE *inputDihedral);
int countTimeframes_dump (FILE *inputDump, int nAtoms);
int countTimeframes_dihedral (FILE *inputDihedral, int nDihedrals);
DUMP_ENTRIES *readDump (FILE *inputDump, int nAtoms, int *currentTimestep_dump);
int getIndex1d (int i, int j, int width);
int readDihedral (DIHEDRAL_ENTRIES **dihedral, FILE *inputDihedral, DUMP_ENTRIES *dump, int currentStep, int nDihedrals, int nAtoms, int currentTimestep_dump, int *currentTimestep_dihedral);
void computeCorrelation (DIHEDRAL_ENTRIES *dihedral, int nDihedrals, int nTimeframes_dump, int **lHandCorrelation, int **rHandCorrelation, int **exCoilCorrelation, int **coilCorrelation);

#endif