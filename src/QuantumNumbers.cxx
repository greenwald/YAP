#include "QuantumNumbers.h"

namespace yap {

//-------------------------
QuantumNumbers::QuantumNumbers(unsigned char twoJ, char P, char C, char I, char G) :
    twoJ_(twoJ),
    P_(P),
    C_(C),
    I_(I),
    G_(G)
{}

//-------------------------
bool operator== (const QuantumNumbers& lhs, const QuantumNumbers& rhs)
{
    return (lhs.twoJ_ == rhs.twoJ_
            && lhs.P_ == rhs.P_
            && lhs.C_ == rhs.C_
            && lhs.I_ == rhs.I_
            && lhs.G_ == rhs.G_);
}

}
