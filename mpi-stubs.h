/*
 * MPI stubs for serial runs of QuanSimBench
 *
 * This file defines only the MPI datatypes, functions, and macros
 * required by QuanSimBench.  It is not a general-purpose MPI stub
 * library.
 */

#ifndef _MPI_H_
#define _MPI_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define QSB_MPI_STUBS 1

#ifdef __GNUC__
# define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
# define UNUSED(x) UNUSED_ ## x
#endif

#define MPI_SUCCESS         0
#define MPI_COMM_WORLD     10
#define MPI_SUM            20
#define MPI_DOUBLE         sizeof(double)
#define MPI_COMPLEX        sizeof(float complex)
#define MPI_STATUS_IGNORE  NULL

typedef int MPI_Comm;
typedef int MPI_Op;
typedef int MPI_Request;
typedef int MPI_Status;
typedef size_t MPI_Datatype;

static int
MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
              MPI_Datatype datatype, MPI_Op UNUSED(op), MPI_Comm UNUSED(comm))
{
  memcpy(recvbuf, sendbuf, datatype*count);
  return MPI_SUCCESS;
}

static int
MPI_Comm_rank(MPI_Comm UNUSED(comm), int *rank)
{
  *rank = 0;
  return MPI_SUCCESS;
}

static int
MPI_Comm_size(MPI_Comm UNUSED(comm), int *size)
{
  *size = 1;
  return MPI_SUCCESS;
}

static int
MPI_Finalize()
{
  exit(0);
  return MPI_SUCCESS;
}

static int
MPI_Init(int *UNUSED(argc), char ***UNUSED(argv))
{
  return MPI_SUCCESS;
}

static int
MPI_Irecv(void *UNUSED(buf), int UNUSED(count), MPI_Datatype UNUSED(datatype),
            int UNUSED(source), int UNUSED(tag), MPI_Comm UNUSED(comm),
            MPI_Request *UNUSED(request))
{
  return MPI_SUCCESS;
}

static int
MPI_Isend(void *UNUSED(buf), int UNUSED(count), MPI_Datatype UNUSED(datatype),
            int UNUSED(dest), int UNUSED(tag), MPI_Comm UNUSED(comm),
            MPI_Request *UNUSED(request))
{
  return MPI_SUCCESS;
}

static int
MPI_Wait(MPI_Request *UNUSED(request), MPI_Status *UNUSED(status))
{
  return MPI_SUCCESS;
}

#endif
