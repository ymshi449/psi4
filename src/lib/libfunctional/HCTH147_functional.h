#ifndef HCTH147_functional_h
#define HCTH147_functional_h
/**********************************************************
* HCTH147_functional.h: declarations for HCTH147_functional for KS-DFT
* Robert Parrish, robparrish@gmail.com
* Autogenerated by MATLAB Script on 25-May-2011
*
***********************************************************/
#include "functional.h"

namespace psi { namespace functional {

class HCTH147_Functional : public Functional {
public:
    HCTH147_Functional(int npoints, int deriv);
    virtual ~HCTH147_Functional();
    virtual void computeRKSFunctional(boost::shared_ptr<RKSFunctions> prop);
    virtual void computeUKSFunctional(boost::shared_ptr<UKSFunctions> prop);
};
}}
#endif

