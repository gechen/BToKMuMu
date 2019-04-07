
   
	double matrix[nPar][nPar];  // bin 0,1,9
	gMinuit->mnemat(&matrix[0][0],nPar);
	cout<<"ERROR MATRIX: "<<endl;
	
//	Save Fitting results
	double chi2Val=0;
	double arrPar[nPar], arrParErr[nPar];
	for (int i = 0; i < nPar; i++) {arrPar[i]=0; arrParErr[i]=0;}
	double errPar1[nPar], errPar2[nPar], errPar3[nPar], errPar4[nPar], errPar5[nPar], errPar6[nPar], errPar7[nPar]; 
	for (int iPar = 0; iPar < nPar; iPar++) {
		arrPar[iPar]    = f1_model->GetParameter(iPar);
		arrParErr[iPar] = f1_model->GetParError(iPar);
		chi2Val         = f1_model->GetChisquare();
		errPar1[iPar]   = matrix[0][iPar];
		errPar2[iPar]   = matrix[1][iPar];
		errPar3[iPar]   = matrix[2][iPar];
		errPar4[iPar]   = matrix[3][iPar];
		errPar5[iPar]   = matrix[4][iPar];
		errPar6[iPar]   = matrix[5][iPar];
		errPar7[iPar]   = matrix[6][iPar];
//		if (iPar < 3) arrParErr[iPar] = readParam(iBin,"accErr", iPar+1);
//		else arrParErr[iPar] = readParam(iBin,"recoErr", iPar-3);
	}
	
//	return output;
	writeParam(iBin,"accXrecoEff",   arrPar,   nPar);  // f1_model_format_2
	writeParam(iBin,"accXrecoEffErr",arrParErr,nPar);  // f1_model_format_2
	writeParam(iBin,"EffErr1",       errPar1,  nPar);  
	writeParam(iBin,"EffErr2",       errPar2,  nPar);  
	writeParam(iBin,"EffErr3",       errPar3,  nPar);  
	writeParam(iBin,"EffErr4",       errPar4,  nPar);  
	writeParam(iBin,"EffErr5",       errPar5,  nPar);  
	writeParam(iBin,"EffErr6",       errPar6,  nPar);  
	writeParam(iBin,"EffErr7",       errPar7,  nPar); 
