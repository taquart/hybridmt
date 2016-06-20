//-----------------------------------------------------------------------------
// Source: moment_tensor.h
// Module: FOCIMT
// Main routine.
//
// Copyright (c) 2013-2016, Grzegorz Kwiatek.
//
// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:
//
// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//-----------------------------------------------------------------------------
#ifndef MOMENT_TENSOR_H_
#define MOMENT_TENSOR_H_
//---------------------------------------------------------------------------
#define FOCIMT_MAXCHANNEL 128
#define FOCIMT_MIN_ALLOWED_CHANNELS 8
#define FOCIMT_SQ(x) (pow(x,2.0))
#define FOCIMT_SEP "\t"
#define FOCIMT_SEP2 " "
#define FOCIMT_NEWLINE "\n"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <list>
#include <algorithm>
#include <vector>

#include "trinity_library.h"
/*
#include <triexceptions/exceptions.h>
#include <trilib/fortranmath.h>
#include <trilib/georoutines.h>
#include <trilib/string.h>
#include <trilib/tristat.h>
#include <tricairo/tricairo_meca.h>
*/

//---------------------------------------------------------------------------
namespace Taquart {
  //! Norm type used in calculation of the moment tensor solution.
  /*! Defines norm used in seismic moment tensor inversion.
   * \ingroup Foci
   */
  enum NormType {
    ntL1, /*!< L1 norm used. */
    ntL2 /*!< L2 norm used. */
  };

  //! Seismic moment tensor solution type.
  /*! Determines the type of seismic moment tensor.
   *  \ingroup Foci
   */
  enum SolutionType {
    stFullSolution = 1, /*!< Full solution.*/
    stTraceNullSolution = 2, /*!< Trace-null solution.*/
    stDoubleCoupleSolution = 3 /*!< Double-couple solution.*/
  };
}

//---------------------------------------------------------------------------
#endif /* MOMENT_TENSOR_H_ */
