/******************************************************************************
 * Filename              :   iir_filter.c
 * Author                :   Giulio Dalla Vecchia
 * Origin Date           :   11 dic 2020
 *
 * Copyright (c) 2022 Giulio Dalla Vecchia. All rights reserved.
 *
 ******************************************************************************/

/** @file iir_filter.c
 *  @brief Questo file definisce l'interfaccia per il controllo di filtri
 *         digitali
 */

/*****************************************************************************
 * Includes
 ******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "iir_filter.h"

/**
 * \addtogroup        IIRFilter
 * \brief
 * \{
 */

/*****************************************************************************
 * Module Preprocessor Constants
 ******************************************************************************/

#define FIX_PI    (0x19220) /* PI convertito in formato Q15 */

/*****************************************************************************
 * Module Preprocessor Macros
 ******************************************************************************/

#define SAT(x, MAX) ( ((x) > MAX) ? MAX : ( (-(x) > MAX) ? (-MAX) : (x) ) )

/*****************************************************************************
 * Module Typedefs
 ******************************************************************************/

/*****************************************************************************
 * Function Prototypes
 ******************************************************************************/

/*****************************************************************************
 * Module Variable Definitions
 ******************************************************************************/

/*****************************************************************************
 * Function Definitions
 ******************************************************************************/

/**
 * Questa funzione inizializza i coefficienti del filtro passa-basso del primo
 * ordine discretizzato con il metodo di Tustin [S = 2 / T * ((1 - z^-1) / (1 + z^-1))]
 * I coefficienti del filtro sono tutti in formato Q15
 *
 * L'equazione del filtro risulta: \n
 *
 * Y(z) = A * U(z) + A * U(z - 1) + B * Y(z - 1)
 *
 * dove: A = [T / (T + 2*tau)], B = [(T - 2*tau) / (T + 2*tau)], tau = 1 / (2 * PI * cutOffFreq)
 *
 * T = sampleTime
 *
 * @param me: puntatore alla memoria dove è allocata la struttura del filtro
 * @param cfgPtr: puntatore alle configurazioni del filtro
 */
void iir_lpf_fo_init(IIRFoLpf_t *const me, IIRFoLpfCfg_t const *const cfgPtr) {

  fix_t tau;
  fix_t val;
  fix_t twoTau;
  fix_t sampleTime;

  me->sat = cfgPtr->sat;

  me->LastIn = 0;
  me->LastOut = 0;

  sampleTime = FDIV(INT2FIX(1, 15), cfgPtr->sampleFreq, 15);

  val = FMUL(INT2FIX(2, 15), FIX_PI, 15); /* 2 * PI */
  val = FMUL(val, cfgPtr->cutOffFreq, 15); /* (2 * PI) * cutOffFreq */

  tau = FDIV(INT2FIX(1, 15), val, 15); /* 1 / ((2 * PI) * cutOffFreq) */

  twoTau = FMUL(INT2FIX(2, 15), tau, 15); /* 2 * tau */
  val = FADD(twoTau, sampleTime); /* (2 * tau) + sampleTime */
  me->A = FDIV(sampleTime, val, 15); /* sampleTime / ((2 * tau) + sampleTime) */

  /* ((2 * tau) - sampleTime) / ((2 * tau) + sampleTime) */
  me->B = FDIV(FSUB(twoTau, sampleTime), val, 15);
}

/**
 * Questa funzione processa il campione in ingresso tramite il filtro digitale
 * passa-basso del primo ordine inizializzato precedentemente.
 *
 * Y(z) = A * U(z) + A * U(z - 1) + B * Y(z - 1)
 *
 * @param me: puntatore al filtro digitale
 * @param in: buffer di dati in ingresso in formato Q15
 * @param out: buffer di dati in uscita in formato Q15
 * @param len: numero di campioni
 * @return: ultimo valore in uscita dal filtro
 */
fix_t iir_lpf_fo_process(IIRFoLpf_t *const me, fix_t const *in,
    fix_t *const out, size_t len) {

  fix_t x;
  fix_t y;
  fix_t z;
  fix_t val;

  for (uint32_t i = 0U; i < len; ++i) {
    x = FMUL(me->A, in[i], 15);
    y = FMUL(me->A, me->LastIn, 15);
    z = FMUL(me->B, me->LastOut, 15);
    val = FADD(x, y);
    me->LastOut = SAT(FADD(val,z), me->sat);
    if(out != NULL) {
      out[i] = me->LastOut;
    }
    me->LastIn = in[i];
  }

  return me->LastOut;
}

/**
 * Questa funzione inizializza i coefficienti del filtro passa-banda del secondo
 * ordine discretizzato con il metodo di Tustin [S = 2 / T * ((1 - z^-1) / (1 + z^-1))]
 * I coefficienti del filtro sono tutti in formato Q15
 *
 * L'equazione del filtro risulta: \n
 *
 * Y(z) = [D * U(z - 1) + D * U(z) - A * Y(z - 1) - B * Y(z - 2)] / C
 *
 * dove:
 *
 * A = 4 - ((2 * T * wn) / Q) + (wn^2 * T^2)
 * B = [2 * (wn^2 * T^2)] - 8
 * C = 4 + ((2 * T * wn) / Q) + (wn^2 * T^2)
 * D = ((2 * T * wn) / Q)
 *
 * T = sampleTime
 *
 * @param me: puntatore alla memoria dove è allocata la struttura del filtro
 * @param cfgPtr: puntatore alle configurazioni del filtro
 */
void iir_bpf_so_init(IIRSoBpf_t *const me, IIRSoBpfCfg_t const *const cfgPtr) {

  fix_t wn;
  fix_t wn2;
  fix_t x;
  fix_t y;
  fix_t sampleTime;
  fix_t sampleTime2;
  fix_t val;

  sampleTime = FDIV(INT2FIX(1, 15), cfgPtr->sampleFreq, 15);
  sampleTime2 = FMUL(sampleTime, sampleTime, 15);

  /* Calcolo Wn del filtro */
  wn = FMUL(FIX_PI, INT2FIX(2, 15), 15);
  wn = FMUL(wn, cfgPtr->centralFreq, 15);

  /* Calcolo Wn al quadrato */
  wn2 = FMUL(wn, wn, 15);

  /* ((Wn^2) * (sampleTime^2)) */
  x = FMUL(wn2, sampleTime2, 15);

  /* (2 * sampleTime * Wn) */
  y = FMUL(INT2FIX(2, 15), sampleTime, 15);
  y = FMUL(y, wn, 15);

  /* ((2 * sampleTime * Wn) / Q) */
  val = FDIV(y, cfgPtr->qFactor, 15);
  me->D = -val;

  /* ((2 * sampleTime * Wn) / Q) + ((Wn^2) * (sampleTime^2)) */
  val = FADD(x, val);

  /* 4 - ((2 * sampleTime * Wn) / Q) + ((Wn^2) * (sampleTime^2)) */
  me->A = FSUB(INT2FIX(4, 15), val);

  /* 4 + ((2 * sampleTime * Wn) / Q) + ((Wn^2) * (sampleTime^2)) */
  me->C = FADD(INT2FIX(4, 15), val);

  /* (2 * ((Wn^2) * (sampleTime^2))) - 8 */
  val = FMUL(INT2FIX(2, 15), x, 15);
  me->B = FSUB(val, INT2FIX(8, 15));

  me->index = 0U;
  me->sat = cfgPtr->sat;
  me->LastIn[0] = 0U;
  me->LastIn[1] = 0U;
  me->LastOut[0] = 0U;
  me->LastOut[1] = 0U;
}

/**
 * Questa funzione processa il campione in ingresso tramite il filtro digitale
 * passa-banda del secondo ordine inizializzato precedentemente.
 *
 * Y(z) = [D * U(z - 1) + D * U(z) - A * Y(z - 1) - B * Y(z - 2)] / C
 *
 * @param me: puntatore al filtro digitale
 * @param in: buffer di dati in ingresso in formato Q15
 * @param out: buffer di dati in uscita in formato Q15
 * @param len: numero di campioni
 * @return: ultimo valore in uscita dal filtro
 */
fix_t iir_bpf_so_process(IIRSoBpf_t *const me, fix_t const *in,
    fix_t *const out, size_t len) {

  fix_t x, y, z, w;
  fix_t val;
  uint32_t lastIndex;

  for (uint32_t i = 0U; i < len; ++i) {
    lastIndex = (me->index == 0U) ? 1U : 0U;

    x = FMUL(me->D, me->LastIn[me->index], 15);
    y = FMUL(me->D, in[i], 15);
    z = FMUL(me->A, me->LastOut[me->index], 15);
    w = FMUL(me->B, me->LastOut[lastIndex], 15);

    val = FADD(x, y);
    val = FADD(val, z);
    val = FADD(val, w);
    val = FDIV(val, me->C, 15);
    me->LastOut[me->index] = SAT(val, me->sat);

    me->LastIn[me->index] = in[i];

    me->index = (me->index + 1U) % 2U;
  }

  return me->LastOut[(me->index == 1U) ? 0U : 1U];
}

/**
 * \}
 */
