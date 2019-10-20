/*
 Copyright 2016 Emanuele Vespa, Imperial College London

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef VOLUME_H
#define VOLUME_H

// Data types definitions
#include <iostream>
#include <se/voxel_traits.hpp>

/******************************************************************************
 *
 * KFusion Truncated Signed Distance Function voxel traits
 *
****************************************************************************/

typedef struct {
  float x;
  float y;
} SDF;

template<>
struct voxel_traits<SDF> {
  typedef SDF value_type;
  static inline value_type empty(){ return {1.f, -1.f}; }
  static inline value_type initValue(){ return {1.f, 0.f}; }
};

/******************************************************************************
 *
 * Bayesian Fusion voxel traits and algorithm specificic defines
 * state: -1 unknown, 0 free, 1 occupied
 *
******************************************************************************/


typedef enum class voxel_state: std::int8_t
{
  kFree     = 0,
  kFrontier = 1,
  kUnknown  = 2,
  kOccupied = 3
};
typedef struct {
    float x;
    double y;
    voxel_state st;
} OFusion;

template<>
struct voxel_traits<OFusion> {
  typedef struct  {
    float x; // occ prob
    double y; // timestamp
    voxel_state st;
  } value_type;
  static inline value_type empty(){ return {0.f, 0.f, voxel_state::kUnknown}; }
  static inline value_type initValue(){ return {0.f, 0.f, voxel_state::kUnknown}; }
};
static  std::ostream& operator<<(std::ostream& os, const voxel_state & dt)
{
  return os << static_cast<int>(dt);
}

// Windowing parameters
#define DELTA_T   1.f
#define CAPITAL_T 4.f

#define INTERP_THRESH 0.05f
#define SURF_BOUNDARY 0.f
#define TOP_CLAMP     1000.f
#define BOTTOM_CLAMP  (-TOP_CLAMP)
#define THRESH_OCC 0.8f
#define THRESH_FREE 0.2f
#define THRESH_FREE_LOG log2(0.03f / (1.f - 0.03f))
#define THRESH_OCC_LOG log2(0.8f / (1.f - 0.8f))
#endif