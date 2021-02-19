#ifndef _HYDROLOGICAL_TOOLS_PALMER_DROUGHT_SEVERITY INDEX_SC_PDSI_
#define _HYDROLOGICAL_TOOLS_PALMER_DROUGHT_SEVERITY INDEX_SC_PDSI_

#include <list>
class SCPDSI
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
    SCPDSI(const char* records);
    SCPDSI(const std::list<double>& P, const std::list<double>& PET);

    std::list<double> Weekly(const std::list<double>& P, const std::list<double>& PET);
    std::list<double> Monthly(const std::list<double>& AWC);
};
#endif //!_HYDROLOGICAL_TOOLS_PALMER_DROUGHT_SEVERITY INDEX_SC_PDSI_