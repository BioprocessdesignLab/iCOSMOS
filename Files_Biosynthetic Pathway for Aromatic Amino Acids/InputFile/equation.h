void formula(complex<double> *x, complex<double> *fx, complex<double> **v, complex<double> *sv, int nd1, int ni, int cv, complex<double> **stoin, double t)
{
	int i, j;
	//*****************************************************************************
	complex<double> s0;
	complex<double> s1;
	//****************************************************************   
	int cv1;
	cv1 = cv + 1;

	//�����ł�sv�́A(NV)�̂���
	//stoin�͕��f��??

	double k021 = 0.00062, k022 = 0.11, k023 = 0.1;
	double k031 = 193.0, k032 = 2.558, k033 = 2.326;
	double k041 = 6253.0; //double R4,L4,T4;
	double k051 = 0.25, k052 = 6.273; //double A5;
	double k061 = 0.3964, k062 = 311.2; //double R6,L6,T6;
	double k101 = 2.0, k102 = 0.008665, k103 = 0.9281;//k103=0.921
	double k111 = 50.0, k112 = 40.0;
	double k121 = 50.0;
	double k131 = 1.0;
	double k141 = 1.269, k142 = 2.0, k143 = 0.9281, k144 = 16.0;
	double k151 = 0.8;

	//�F�����A�~�m�_�������f��
	complex<double> sv1, sv2, sv3, sv4, sv5, sv6, sv7, sv8, sv9, sv10, ADP, R, L, T, A;
	complex<double> sv11, sv12, sv13, sv14, sv15, sv16, sv17, sv18;

	//Original Model(MM�^��)
	sv1 = x[13] - (132.5 * x[2]);//original
	sv2 = x[14] / (k021 / (x[1] * x[5]) + k022 / x[1] + k023 / x[5] + 1.0);
	sv3 = x[15] * pow(x[2], 8.51) / (k031 + pow(x[2], 8.51)) / (k032 / x[2] / x[2] + k033 / x[2] + 1.0); // original
	ADP = (-x[5] + sqrt(12.0*x[5] - 3.0*x[5] * x[5])) / 2.0;
	R = 1.0 + 0.5714 * x[2] + 16.67 * x[5] + 95.24 * x[2] * x[5];
	L = (1.0 + 0.76 *(3.0 - x[5] - ADP)) / (1.0 + 40.0 * (3.0 - x[5] - ADP));
	T = 1.0 + 0.0002857 * x[2] + 16.67 * x[5] + 0.004762 * x[2] * x[5];
	sv4 = x[16] * x[2] * x[5] * R / (R * R + k041 * L * L * T * T);
	A = 1.0 + (3.0 - ADP - x[5]) / 1.1 + ADP / 1.5 + x[5] / 2.5; //Stephanopoulos
	sv5 = x[17] / (1.0 + k051 / x[3] + (0.09375 + k052 / x[3]) * A);
	R = 1.0 + 157.0 * x[4] + 0.2 * ADP + 3.14 * x[4] * ADP;
	L = (1.0 + 0.05 * x[3]) / (1.0 + 5.0 * x[3]);
	T = 1.0 + 0.02 * x[4] + 0.2 * ADP + 0.004 * x[4] * ADP;
	sv6 = x[18] * x[4] * ADP * (R + k061 * L * L * T) / (R * R + k062 * L * L * T * T);
	sv7 = x[19] * sv5;
	sv8 = x[20] * x[5];
	sv9 = x[21] * (0.79 / (1.0 + x[10] / 53.0) + 0.2 / (1.0 + x[11] / 40.0) + 0.01 / (1.0 + x[12] / 16.0)) / ((0.0002 / x[4] / x[2] + 0.006 / x[4]) * (1.0 + 50.0 * x[6]) + 0.1 / x[2] + 1.0);//0.006or0.06
	sv10 = x[22] * x[6] * x[4] * x[5] / (k101 + x[6]) / (k102 + x[4]) / (k103 + x[5]);
	sv11 = x[23] * x[7] / (2.0 + x[7]) / (1.0 + x[10] / k111) / (1.0 + x[11] / k112);
	sv12 = x[24] * x[8] / (1.0 + x[8]) / (1.0 + x[10] / k121);
	sv13 = x[25] * x[8] / (k131 + x[8]);
	sv14 = x[26] * x[9] * x[7] * x[5] / (k141 + x[9]) / (k142 + x[7]) / (k143 + x[5]) / (1.0 + x[12] / k144);
	sv15 = x[27] * ADP * x[9] / (k151 + x[9]);
	sv16 = x[28] * x[10];
	sv17 = x[29] * x[11];
	sv18 = x[30] * x[12];

	sv[1] = sv1;
	sv[2] = sv2;
	sv[3] = sv3;
	sv[4] = sv4;
	sv[5] = sv5;
	sv[6] = sv6;
	sv[7] = sv7;
	sv[8] = sv8;
	sv[9] = sv9;
	sv[10] = sv10;
	sv[11] = sv11;
	sv[12] = sv12;
	sv[13] = sv13;
	sv[14] = sv14;
	sv[15] = sv15;
	sv[16] = sv16;
	sv[17] = sv17;
	sv[18] = sv18;


	//****************************************************************

	for (i = 1; i<nd1; i++) {
		s0 = 0.0;
		s1 = 0.0;
		for (j = 1; j< cv1; j++) {
			if (real(stoin[i][j])>0.0) s0 = s0 + stoin[i][j] * sv[j];
			if (real(stoin[i][j])<0.0) s1 = s1 + stoin[i][j] * sv[j];
//			printf("%23.15e\n", s1);
		}
		v[i][0] = s0;
		v[i][1] = -s1;
//		printf("%23.15e\t%23.15e\n",s0,s1);
		
		fx[i] = s0 + s1;
	}
//	getchar();

}