import time
from sklearn.linear_model import LinearRegression
import numpy as np
import decimal

def extract_numbers (line):
	lineString = line.replace ("\t", ",").replace (" ", ",").replace ("\n", "")
	for item in lineString.split (","):
		try:
			yield decimal.Decimal (item)
		except:
			pass

def main ():
	trans = []
	gauchePlus = []
	gaucheMinus = []
	transTrans = []
	transGaucheMinus = []
	transGauchePlus = []
	gaucheGauche = []

	with open ("chiralCorrelation.output", "r") as inputCorrelationFile:
		for line in inputCorrelationFile:
			lineArray = list (extract_numbers (line))
			if (len (lineArray) > 0):
				trans.append (lineArray[0])
				gauchePlus.append (lineArray[1])
				gaucheMinus.append (lineArray[2])
				transTrans.append (lineArray[3])
				transGaucheMinus.append (lineArray[4])
				transGauchePlus.append (lineArray[5])
				gaucheGauche.append (lineArray[6])
			
	nTrans = np.array (trans [int (len (trans) / 2):])
	nGauchePlus = np.array (gauchePlus [int (len (gauchePlus) / 2):])
	nGaucheMinus = np.array (gaucheMinus [int (len (gaucheMinus) / 2):])
	nTransTrans = np.array (transTrans [int (len (transTrans) / 2):])
	nTransGaucheMinus = np.array (transGaucheMinus [int (len (transGaucheMinus) / 2):])
	nTransGauchePlus = np.array (transGauchePlus [int (len (transGauchePlus) / 2):])
	nGaucheGauche = np.array (gaucheGauche [int (len (gaucheGauche) / 2):])

	nTrans_x = np.array ([x for x in range (0, nTrans.size)]).reshape((-1, 1))
	nGauchePlus_x = np.array ([x for x in range (0, nGauchePlus.size)]).reshape((-1, 1))
	nGaucheMinus_x = np.array ([x for x in range (0, nGaucheMinus.size)]).reshape((-1, 1))
	nTransTrans_x = np.array ([x for x in range (0, nTransTrans.size)]).reshape((-1, 1))
	nTransGaucheMinus_x = np.array ([x for x in range (0, nTransGaucheMinus.size)]).reshape((-1, 1))
	nTransGauchePlus_x = np.array ([x for x in range (0, nTransGauchePlus.size)]).reshape((-1, 1))
	nGaucheGauche_x = np.array ([x for x in range (0, nGaucheGauche.size)]).reshape((-1, 1))

	model = LinearRegression ()

	model.fit (nTrans_x, nTrans)
	trans_rSquare = model.score (nTrans_x, nTrans)
	trans_intercept = model.intercept_
	trans_slope = model.coef_

	model.fit (nGauchePlus_x, nGauchePlus)
	gauchePlus_rSquare = model.score (nGauchePlus_x, nGauchePlus)
	gauchePlus_intercept = model.intercept_
	gauchePlus_slope = model.coef_

	model.fit (nGaucheMinus_x, nGaucheMinus)
	gaucheMinus_rSquare = model.score (nGaucheMinus_x, nGaucheMinus)
	gaucheMinus_intercept = model.intercept_
	gaucheMinus_slope = model.coef_

	model.fit (nTransTrans_x, nTransTrans)
	transTrans_rSquare = model.score (nTransTrans_x, nTransTrans)
	transTrans_intercept = model.intercept_
	transTrans_slope = model.coef_

	model.fit (nTransGaucheMinus_x, nTransGaucheMinus)
	transGaucheMinus_rSquare = model.score (nTransGaucheMinus_x, nTransGaucheMinus)
	transGaucheMinus_intercept = model.intercept_
	transGaucheMinus_slope = model.coef_

	model.fit (nTransGauchePlus_x, nTransGauchePlus)
	transGauchePlus_rSquare = model.score (nTransGauchePlus_x, nTransGauchePlus)
	transGauchePlus_intercept = model.intercept_
	transGauchePlus_slope = model.coef_

	model.fit (nGaucheGauche_x, nGaucheGauche)
	gaucheGauche_rSquare = model.score (nGaucheGauche_x, nGaucheGauche)
	gaucheGauche_intercept = model.intercept_
	gaucheGauche_slope = model.coef_

	with open ("chiralCorrelation.fit", "w") as outputFile:
		outputFile.write ("trans: {}\t{}\t{}\ngauchePlus: {}\t{}\t{}\ngaucheMinus: {}\t{}\t{}\ntransTrans: {}\t{}\t{}\ntransGaucheMinus: {}\t{}\t{}\ntransGauchePlus: {}\t{}\t{}\ngaucheGauche: {}\t{}\t{}\n".format (trans_rSquare, trans_intercept, trans_slope, gauchePlus_rSquare, gauchePlus_intercept, gauchePlus_slope, gaucheMinus_rSquare, gaucheMinus_intercept, gaucheMinus_slope, transTrans_rSquare, transTrans_intercept, transTrans_slope, transGaucheMinus_rSquare, transGaucheMinus_intercept, transGaucheMinus_slope, transGauchePlus_rSquare, transGauchePlus_intercept, transGauchePlus_slope, gaucheGauche_rSquare, gaucheGauche_intercept, gaucheGauche_slope))

if __name__ == '__main__':
	main()