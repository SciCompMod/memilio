#ifndef MIO_MPI_H
#define MIO_MPI_H

#include "memilio/config.h"
#ifdef MEMILIO_HAS_MPI
#include "mpi.h"

#include <cstdint>
#include <climits>

// set MPI_Type for size_t
#ifndef MIO_MPI_SIZE_T
#if SIZE_MAX == UCHAR_MAX
#define MIO_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
#define MIO_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
#define MIO_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#define MIO_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#define MIO_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
#error "Cannot determine MPI_Type for size_t, please define MIO_MPI_SIZE_T."
#endif
#endif // MIO_MPI_SIZE_T

// set MPI_Type for the return type of mio::RandomNumberGenerator
#ifndef MIO_MPI_RANDOM_NUMBER_T
#if UINT_FAST64_MAX == UCHAR_MAX
#define MIO_MPI_RANDOM_NUMBER_T MPI_UNSIGNED_CHAR
#elif UINT_FAST64_MAX == USHRT_MAX
#define MIO_MPI_RANDOM_NUMBER_T MPI_UNSIGNED_SHORT
#elif UINT_FAST64_MAX == UINT_MAX
#define MIO_MPI_RANDOM_NUMBER_T MPI_UNSIGNED
#elif UINT_FAST64_MAX == ULONG_MAX
#define MIO_MPI_RANDOM_NUMBER_T MPI_UNSIGNED_LONG
#elif UINT_FAST64_MAX == ULLONG_MAX
#define MIO_MPI_RANDOM_NUMBER_T MPI_UNSIGNED_LONG_LONG
#else
#error "Cannot determine MPI_Type for unit_fast64_t, please define MIO_MPI_RANDOM_NUMBER_T."
#endif
#endif // MIO_MPI_RANDOM_NUMBER_T

#endif // MEMILIO_HAS_MPI

#endif // MIO_MPI_H