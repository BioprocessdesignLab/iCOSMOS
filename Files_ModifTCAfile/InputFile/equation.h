void formula(complex<double> *x, complex<double> *fx, complex<double> **v, complex<double> *sv, int nd1, int ni, int cv, complex<double> **stoin, double t)
{
	int i, j;
	//*****************************************************************************
	complex<double> s0;
	complex<double> s1;
	//****************************************************************   
	int cv1;
	cv1 = cv + 1;
	double M = 1.0e-3;
	double M1 = 1.0e+3;
	//****************************************************************
	complex<double> s, sv14, sv15, sv16, sv17, sv18, sv19, sv20, sv21, sv22, sv23, sv24;
	complex<double> sv25, sv26, sv27, sv28, sv29, sv30, sv31, sv32, sv33, sv34, sv35, sv36;
	complex<double> sv37, sv38, sv39, sv40, sv41, sv42, sv43, sv44, sv45, denom1, denom2, denom3, denom;

	sv14 = M * x[14] * M * x[6] * M * x[46] / (0.2 * M * x[6] + 2.0 * M * x[46] * (1.0 + M * x[48] / 0.025) + M * x[6] * M * x[46]);
	denom1 = 0.0700 * M * x[10] * M * x[47] + 0.002 * M * x[10] * M * x[46] + M * x[46] * M * x[47] + M * x[10] * M * x[46] * M * x[47] + 1.5 * M * x[11] * M * x[48];
	denom2 = 0.0700 * M * x[10] * M * x[47] * M * x[48] / 0.018 + 0.002 * M * x[10] * M * x[46] * M * x[11] + 1.5 * M * x[10] * M * x[11] * M * x[48] / 0.75;
	denom = denom1 + denom2;
	sv15 = M * x[15] * M * x[10] * M * x[46] * M * x[47] / denom;
	sv16 = M * x[16] * (M * x[11] - M * x[12] / 10.0) / (0.1 + M * x[11] + M * x[12] / 10.0);
	sv17 = M * x[17] * (M * x[12] - M * x[13] / 10.0) / (0.1 + M * x[12] + M * x[13] / 10.0);
	denom1 = 0.31 * 1.33 + 1.33 * M * x[46] + 0.1 * M * x[13] + M * x[46] * M * x[13] + 0.27 * M * x[48] + 0.04 * M * x[1] + M * x[48] * M * x[1] + 0.27 * M * x[46] * M * x[48] / 0.31 + 0.100 * M * x[13] * M * x[1] / 0.27;
	denom2 = M * x[46] * M * x[13] * M * x[48] / 0.04 + M * x[13] * M * x[48] * M * x[1] / 3.3 + M * x[46] * M * x[13] * M * x[1] / 0.17 + M * x[46] * M * x[48] * M * x[1] / 0.31;
	denom3 = 0.31 * M * x[46] * M * x[13] * M * x[48] * M * x[1] / (0.31 * 0.27 * 0.04);
	denom = denom1 + denom2 + denom3;
	sv18 = M * x[18] * (M * x[13] * M * x[46] - M * x[48] * M * x[1]) / denom;
	sv19 = M * x[19] * M * x[13] * M * x[7] / ((0.37 + M * x[13]) * (0.1 + M * x[7]));
	sv20 = M * x[20] * M * x[8];
	denom1 = 0.11 * M * x[5] * M * x[47] + 0.01 * M * x[5] * M * x[46] + 0.14 * M * x[47] * M * x[46] + M * x[5] * M * x[47] * M * x[46] + 0.14 * 1.70 * M * x[3] * M * x[48] / 0.02 + 0.11 * M * x[5] * M * x[47] * M * x[48] / 0.05;
	denom2 = 0.01 * M * x[5] * M * x[46] * M * x[3] / 0.02 + 0.14 * 1.7 * M * x[5] * M * x[3] * M * x[48] / (0.18 * 0.02);
	denom = denom1 + denom2;
	sv21 = M * x[21] * M * x[5] * M * x[47] * M * x[46] / denom;
	sv22 = M * x[22] * M * x[2];
	sv23 = M * x[23] * M * x[7];
	sv24 = M * x[24] * M * x[3] * M * x[2] / (0.007 * M * x[3] + 0.01 * (1.0 + M * x[47] / 0.11) * M * x[2] + M * x[3] * M * x[2]);
	sv25 = M * x[25] * M * x[9];
	sv26 = M * x[26] * M * x[4] * M * x[46] / (0.34 * M * x[4] + 0.13 * (1.0 + M * x[48] / 0.02) * M * x[46] + M * x[4] * M * x[46]);
	sv27 = M * x[27] * M * x[6];
	denom = 0.33 * M * x[7] + 0.46 * M * x[10] + M * x[7] * M * x[10] + 9.4 * M * x[1] / 9.5 + 0.1 * M * x[6] / 9.5 + M * x[1] * M * x[6] / 9.5 + 9.4 * M * x[7] * M * x[1] / (0.46 * 9.5) + 0.46 * M * x[10] * M * x[6] / 9.4;
	sv28 = M * x[28] * (M * x[7] * M * x[10] - M * x[1] * M * x[6] / 9.5) / denom;
	denom = 0.19 * M * x[8] + 0.43 * M * x[10] + M * x[8] * M * x[10] + 15.0 * M * x[5] / 9.5 + 0.87 * M * x[6] / 9.5 + M * x[5] * M * x[6] / 9.5 + 15.0 * M * x[8] * M * x[5] / (0.43 * 9.5) + 0.43 * M * x[10] * M * x[6] / 15.0;
	sv29 = -M * x[29] * (M * x[8] * M * x[10] - M * x[5] * M * x[6] / 9.5) / denom;
	sv30 = M * x[30] * M * x[1];
	sv31 = M * x[31] * M * x[7];
	sv32 = M * x[32] * M * x[11];
	sv33 = M * x[33] * M * x[1];
	sv34 = M * x[34];
	sv35 = M * x[35];
	sv36 = M * x[36];
	sv37 = M * x[37];
	sv38 = M * x[38];
	sv39 = M * x[39];
	sv40 = M * x[40] * M * x[7];
	sv41 = M * x[41] * M * x[3];
	sv42 = M * x[42] * M * x[11];
	sv43 = M * x[43] * M * x[12];
	sv44 = M * x[44] * M * x[8];
	sv45 = M * x[45] * M * x[6];


	sv[1] = M1 * sv14;
	sv[2] = M1 * sv15;
	sv[3] = M1 * sv16;
	sv[4] = M1 * sv17;
	sv[5] = M1 * sv18;
	sv[6] = M1 * sv19;
	sv[7] = M1 * sv20;
	sv[8] = M1 * sv21;
	sv[9] = M1 * sv22;
	sv[10] = M1 * sv23;
	sv[11] = M1 * sv24;
	sv[12] = M1 * sv25;
	sv[13] = M1 * sv26;
	sv[14] = M1 * sv27;
	sv[15] = M1 * sv28;
	sv[16] = M1 * sv29;
	sv[17] = M1 * sv30;
	sv[18] = M1 * sv31;
	sv[19] = M1 * sv32;
	sv[20] = M1 * sv33;
	sv[21] = M1 * sv34;
	sv[22] = M1 * sv35;
	sv[23] = M1 * sv36;
	sv[24] = M1 * sv37;
	sv[25] = M1 * sv38;
	sv[26] = M1 * sv39;
	sv[27] = M1 * sv40;
	sv[28] = M1 * sv41;
	sv[29] = M1 * sv42;
	sv[30] = M1 * sv43;
	sv[31] = M1 * sv44;
	sv[32] = M1 * sv45;


	//****************************************************************

	for (i = 1; i < nd1; i++) {
		s0 = 0.0;
		s1 = 0.0;
		for (j = 1; j < cv1; j++) {

			if (real(stoin[i][j]) > 0.0) s0 = s0 + stoin[i][j] * sv[j];
			if (real(stoin[i][j]) < 0.0) s1 = s1 + stoin[i][j] * sv[j];
		}

		v[i][0] = s0;
		v[i][1] = -s1;
		fx[i] = s0 + s1;
	}
}