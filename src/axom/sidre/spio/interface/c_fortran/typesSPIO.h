// typesSPIO.h
// This is generated code, do not edit
//
// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
//
// Produced at the Lawrence Livermore National Laboratory
//
// LLNL-CODE-741217
//
// All rights reserved.
//
// This file is part of Axom.
//
// For details about use and distribution, please read axom/LICENSE.
//
// For C users and C++ implementation

#ifndef TYPESSPIO_H
#define TYPESSPIO_H


#ifdef __cplusplus
extern "C" {
#endif

struct s_SPIO_iomanager
{
  void* addr;       /* address of C++ memory */
  int idtor;        /* index of destructor */
};
typedef struct s_SPIO_iomanager SPIO_iomanager;

struct s_SPI_SHROUD_capsule_data
{
  void* addr;       /* address of C++ memory */
  int idtor;        /* index of destructor */
};
typedef struct s_SPI_SHROUD_capsule_data SPI_SHROUD_capsule_data;

void SPIO_SHROUD_memory_destructor(SPI_SHROUD_capsule_data* cap);

#ifdef __cplusplus
}
#endif

#endif  // TYPESSPIO_H
