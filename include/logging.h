#ifndef yap_logging_h
#define yap_logging_h

/// \file

#include "easylogging++.h"

#include <string>

namespace yap {

/**
 * \file logging.h
 * \brief Logging system using easylogging++
 *
 * Usage:
 * you have to include this header and call the following macro line in your program:
 * \code
 * #include "logging.h"
 * INITIALIZE_EASYLOGGINGPP
 * \endcode
 *
 * To create a log message, put the following line into the code:
 * \code
 * LOG(level) << "put your log message here";
 * \endcode
 *
 * Available levels:
 * INFO, WARNING, DEBUG, ERROR, FATAL, TRACE, VERBOSE
 *
 * To completely disable DEBUG logging messages on preprocessor level, use the macro DEBUG(...),
 * and define ELPP_DISABLE_DEBUG_LOGS
 *
 * A complete manual for easylogging++ is available at:
 * https://github.com/easylogging/easyloggingpp/blob/master/README.md
 */

#define Max_Log_File_Size 1000

/// disable logging for lvl
/// \param lvl (Global, Trace, Debug, Fatal, Error, Warning, Verbose, Info)
inline void disableLogs(el::Level lvl)
{
    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.set(lvl, el::ConfigurationType::Enabled, "0");
    el::Loggers::reconfigureLogger("default", defaultConf);
}

/// just print the debug message without any additional info
/// \param lvl (Global, Trace, Debug, Fatal, Error, Warning, Verbose, Info)
inline void plainLogs(el::Level lvl)
{
    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.set(lvl, el::ConfigurationType::Format, "%msg");
    el::Loggers::reconfigureLogger("default", defaultConf);
}

/// \def DEBUG(x)
/// Provides way to have debug output ignored by the compiler.
#ifdef ELPP_DISABLE_DEBUG_LOGS
#define DEBUG(x)
#else
#define DEBUG(x) LOG(DEBUG) << x;
#endif

/// \def FDEBUG(x)
/// Provides way to have debug output ignored by the compiler.
/// Pretty logging output: prepends function name to x
#ifdef ELPP_DISABLE_DEBUG_LOGS
#define FDEBUG(x)
#else
#define FDEBUG(x) LOG(DEBUG) << std::string(ELPP_FUNC) + ": " << x;
#endif

/// \def FLOG(x)
/// Pretty logging output: prepends function name to x
#define FLOG(x) LOG( x ) << std::string(ELPP_FUNC) + ": "

/// multiline stream
template <typename T>
T& multiline(T& t, std::string s, std::string end_of_line = "")
{
    for (size_t p = 0; p < s.size() and p != std::string::npos; ++p) {
        size_t p_n = s.find_first_of('\n', p);
        if (p_n != std::string::npos)
            t << s.substr(p, p_n - p) << end_of_line;
        else
            t << s.substr(p);
        p = p_n;
    }
    return t;
}

}

#endif
