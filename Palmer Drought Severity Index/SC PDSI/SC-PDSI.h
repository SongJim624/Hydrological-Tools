#ifndef _SC_PDSI_
#define _SC_PDSI_
class PDSI
{
private:
//actual value
    double E, RE, RO, L; 
//potential value    
//    double Ep, REp, ROp, Lp; 
//coefficients
    double alpha, beta, gamma, delta;
//climatic constant
    double K;
//duration factors
    double a, b;

private:
    double Calibration();

public:
    PDSI();
    double SC_PDSI()
};
#endif //!_SC_PDSI_