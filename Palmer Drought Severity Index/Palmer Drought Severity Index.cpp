/*
* The PDSI does not consider the human impacts on the water balance,
* such as the irrigation.
*/

/* Table 
*	4.0	or	more		| Extremely wet
*	3.0	to 3.99		| Very wet
*	2.0	to 2.99		| Moderately wet
*	1.0	to 1.99		| Slightly wet
*	0.5	to 0.99		| incipent wet spell
*	0.49 to -0.49	| Near normal
*	-0.5 to - 0.99	| incipent dry spell
*	-1.0 to -1.99	| Mild drought
*	-2.0 to -2.99	| Moderate drought
*	-3.0 to -3.99	| Severe drought
*	-4.0 or less	| Extreme drought
*/

/*
* Usually computed monthly, but could also be calculated weekly.
*/
double PDSI(const double& precipitation, const double& temperature, const double& AWC)
{
	double alpha, beta, gamma, delta;
	double Ssi, Sui;

	double Rp = AWC - (Ss - +Sui);
	double ROp = Ssi + Sui;
/*
	double PLs = min(ETp, Sui)
	double PL = PLs + (ETp - PLs) * Sui / AWC;
*/
	double Lp = (ETp + (AWC / Sui - 1) * fmin(ETp, Sui)) * Sui / AWC;

	double ET_hat = ETp * alpha;
	double R_hat = Rp * beta;
	double RO_hat = RO * gamma;
	double L_hat = Lp * delta;

	double P_hat = ET_hat + R_hat + RO_hat - L_hat;
	double d = precipitation - P_hat;

	double K = 1.5 * ln(((ETp + R + RO) / (P + L) + 2.8) / d) + 0.5;


	double Z = K * d;


	return 0.0;
}

double PHDI()
{
	return 0.0;
}

double SPI()
{
	return 0.0;
}

double SPEI()
{
	return 0.0;
}


//remote sensing
double NDVI()
{
	return 0.0;
}

double NDWI()
{
	return 0.0;
}

double TDVI(const double& Ts, const double& Tsmin, const double& Tsmax)
{
	return (Ts - Tsmin) / (Tsmax - Tsmin);
}

double VSWI()
{
	return 0.0;
}