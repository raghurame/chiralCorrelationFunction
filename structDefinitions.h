#ifndef CHIRALCORRELATION_STRUCTS_H
#define CHIRALCORRELATION_STRUCTS_H

typedef struct dumpEntries
{
	int id, molType, atomType, ix, iy, iz;
	float x, y, z;
} DUMP_ENTRIES;

typedef struct dihedralEntries
{
	int id, atom1, atom2, atom3, atom4;
	int isTrans, isGauchePlus, isGaucheMinus, isBackbone;
	int isTransTrans_previous, isTransGauchePlus_previous, isTransGaucheMinus_previous, isGaucheGauche_previous;
	int isTransTrans_next, isTransGauchePlus_next, isTransGaucheMinus_next, isGaucheGauche_next;
	float angle;
} DIHEDRAL_ENTRIES;

typedef struct correlation
{
	float correlationTrans, correlationGauchePlus, correlationGaucheMinus;
	float correlationTransTrans, correlationTransGaucheMinus, correlationTransGauchePlus, correlationGaucheGauche;
} CHIRAL_CORRELATION;

#endif