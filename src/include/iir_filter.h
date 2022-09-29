/*****************************************************************************
 * Filename              :   iir_filter.h
 * Author                :   Giulio Dalla Vecchia
 * Origin Date           :   11 dic 2020
 *
 * Copyright (c) 2022 Giulio Dalla Vecchia. All rights reserved.
 *
 ******************************************************************************/

/** @file iir_filter.h
 *  @brief Questo file mostra l'interfaccia per la gestione di filtri di tipo IIR
 *         passa-basso o passa-banda
 */

#ifndef IIR_FILTER_H__
#define IIR_FILTER_H__

/*****************************************************************************
 * Includes
 ******************************************************************************/
#include <stdint.h>
#include "lw_math.h"

#ifdef __cplusplus
extern "C"{
#endif

/**
 * \defgroup        IIRFilter IIR Filter
 * \brief
 * \{
 */

/*****************************************************************************
 * Module Preprocessor Constants
 ******************************************************************************/

/*****************************************************************************
 * Module Preprocessor Macros
 ******************************************************************************/

/*****************************************************************************
 * Module Typedefs
 ******************************************************************************/

/**
 * Struttura che descrive un filtro passa-basso del primo ordine
 */
typedef struct {
  fix_t A;
  fix_t B;
  fix_t LastOut;
  fix_t LastIn;
  fix_t sat;
} IIRFoLpf_t;

/**
 * Struttura che descrive le configurazioni necessarie per
 * un filtro passa-basso del primo ordine
 */
typedef struct {
  fix_t cutOffFreq;
  fix_t sampleFreq;
  fix_t sat;
} IIRFoLpfCfg_t;

/**
 * Struttura che descrive un filtro passa-banda del secondo ordine
 */
typedef struct {
  fix_t A;
  fix_t B;
  fix_t C;
  fix_t D;
  fix_t LastOut[2];
  fix_t LastIn[2];
  fix_t sat;
  uint32_t index;
} IIRSoBpf_t;

/**
 * Struttura che descrive le configurazioni necessarie per
 * un filtro passa-banda del secondo ordine
 */
typedef struct {
  fix_t centralFreq;
  fix_t sampleFreq;
  fix_t qFactor;
  fix_t sat;
} IIRSoBpfCfg_t;

/*****************************************************************************
 * Module Variable Definitions
 ******************************************************************************/

/*****************************************************************************
 * Function Prototypes
 ******************************************************************************/

/** API for low pass filter */

void iir_lpf_fo_init(IIRFoLpf_t *const me, IIRFoLpfCfg_t const *const cfgPtr);

fix_t iir_lpf_fo_process(IIRFoLpf_t *const me, fix_t const *in,
    fix_t *const out, size_t len);

/** API for band pass filter */

void iir_bpf_so_init(IIRSoBpf_t *const me, IIRSoBpfCfg_t const *const cfgPtr);

fix_t iir_bpf_so_process(IIRSoBpf_t *const me, fix_t const *in,
    fix_t *const out, size_t len);

/**
 * \}
 */

#ifdef __cplusplus
} // extern "C"
#endif

#endif /*IIR_FILTER_H__*/

/*** End of File *************************************************************/
