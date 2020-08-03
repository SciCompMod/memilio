#include <epidemiology_io/io.h>

#include <tixi.h>
#include <hdf5.h>

void epi::io_cleanup()
{
    tixiCleanup();
    H5close();
}
