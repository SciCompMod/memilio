#include <epidemiology_io/io.h>
#include "epidemiology/utils/compiler_diagnostics.h"

MSVC_WARNING_DISABLE_PUSH(4268)
#include <tixi.h>
#include <hdf5.h>
MSVC_WARNING_POP

void epi::io_cleanup()
{
    tixiCleanup();
    H5close();
}
