/*************************************************************************
 *                                                                       *
 * Open Dynamics Engine, Copyright (C) 2001,2002 Russell L. Smith.       *
 * All rights reserved.  Email: russ@q12.org   Web: www.q12.org          *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of EITHER:                                  *
 *   (1) The GNU Lesser General Public License as published by the Free  *
 *       Software Foundation; either version 2.1 of the License, or (at  *
 *       your option) any later version. The text of the GNU Lesser      *
 *       General Public License is included with this library in the     *
 *       file LICENSE.TXT.                                               *
 *   (2) The BSD-style license that is included with this library in     *
 *       the file LICENSE-BSD.TXT.                                       *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the files    *
 * LICENSE.TXT and LICENSE-BSD.TXT for more details.                     *
 *                                                                       *
 *************************************************************************/

// object, body, and world structs.


#ifndef _ODE__PRIVATE_OBJECTS_H_
#define _ODE__PRIVATE_OBJECTS_H_

#include "ode/common.h"
#include "ode/memory.h"
#include "error.h"


class dxStepWorkingMemory;
class dxWorldProcessContext;

// some body flags

enum {
    dxBodyFlagFiniteRotation =        1,  // use finite rotations
    dxBodyFlagFiniteRotationAxis =    2,  // use finite rotations only along axis
    dxBodyDisabled =                  4,  // body is disabled
    dxBodyNoGravity =                 8,  // body is not influenced by gravity
    dxBodyAutoDisable =               16, // enable auto-disable on body
    dxBodyLinearDamping =             32, // use linear damping
    dxBodyAngularDamping =            64, // use angular damping
    dxBodyMaxAngularSpeed =           128,// use maximum angular speed
    dxBodyGyroscopic =                256 // use gyroscopic term
};


// base class that does correct object allocation / deallocation

struct dBase {
    void *operator new (size_t size) { return dAlloc (size); }
    void *operator new (size_t, void *p) { return p; }
    void operator delete (void *ptr, size_t size) { dFree (ptr,size); }
    void *operator new[] (size_t size) { return dAlloc (size); }
    void operator delete[] (void *ptr, size_t size) { dFree (ptr,size); }
};


// base class for bodies and joints

struct dObject : public dBase {
    dxWorld *world;		// world this object is in
    dObject *next;		// next object of this type in list
    dObject **tome;		// pointer to previous object's next ptr
    int tag;			// used by dynamics algorithms
    void *userdata;		// user settable data

    explicit dObject(dxWorld *w): world(w), next(NULL), tome(NULL), tag(0), userdata(NULL) {}
    virtual ~dObject();
};


// auto disable parameters
struct dxAutoDisable {
    dReal idle_time;		// time the body needs to be idle to auto-disable it
    int idle_steps;		// steps the body needs to be idle to auto-disable it
    unsigned int average_samples;     // size of the average_lvel and average_avel buffers
    dReal linear_average_threshold;   // linear (squared) average velocity threshold
    dReal angular_average_threshold;  // angular (squared) average velocity threshold

    dxAutoDisable() {}
    explicit dxAutoDisable(void *);
};


// damping parameters
struct dxDampingParameters {
    dReal linear_scale;  // multiply the linear velocity by (1 - scale)
    dReal angular_scale; // multiply the angular velocity by (1 - scale)
    dReal linear_threshold;   // linear (squared) average speed threshold
    dReal angular_threshold;  // angular (squared) average speed threshold

    dxDampingParameters() {}
    explicit dxDampingParameters(void *);
};


// quick-step parameters
struct dxQuickStepParameters {
    int num_iterations;		// number of SOR iterations to perform
    dReal w;			// the SOR over-relaxation parameter

    dxQuickStepParameters() {}
    explicit dxQuickStepParameters(void *);
};


// contact generation parameters
struct dxContactParameters {
    dReal max_vel;		// maximum correcting velocity
    dReal min_depth;		// thickness of 'surface layer'

    dxContactParameters() {}
    explicit dxContactParameters(void *);
};

// position vector and rotation matrix for geometry objects that are not
// connected to bodies.
struct dxPosR {
    dVector3 pos;
    dMatrix3 R;
};

#define dWORLDSTEP_RESERVEFACTOR_DEFAULT    1.2f
#define dWORLDSTEP_RESERVESIZE_DEFAULT      65536U


#endif // #ifndef _ODE__PRIVATE_OBJECTS_H_
