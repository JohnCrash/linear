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
#include "linear.h"
#include "config.h"
#include "ode/common.h"
#include "ode/error.h"
#include "ode/memory.h"
#include "ode/util.h"

#include <new>


#define dMIN(A,B)  ((A)>(B) ? (B) : (A))
#define dMAX(A,B)  ((B)>(A) ? (B) : (A))


//****************************************************************************
// World processing context management

dxWorldProcessMemArena *dxWorldProcessMemArena::ReallocateMemArena (
    dxWorldProcessMemArena *oldarena, size_t memreq, 
    const dxWorldProcessMemoryManager *memmgr, float rsrvfactor, unsigned rsrvminimum)
{
    dxWorldProcessMemArena *arena = oldarena;
    bool allocsuccess = false;

    size_t nOldArenaSize; 
    void *pOldArenaBuffer = NULL;

    do {
        size_t oldmemsize = oldarena ? oldarena->GetMemorySize() : 0;
        if (oldarena == NULL || oldmemsize < memreq) {
            nOldArenaSize = oldarena ? dxWorldProcessMemArena::MakeArenaSize(oldmemsize) : 0;
            pOldArenaBuffer = oldarena ? oldarena->m_pArenaBegin : NULL;

            if (!dxWorldProcessMemArena::IsArenaPossible(memreq)) {
                break;
            }

            size_t arenareq = dxWorldProcessMemArena::MakeArenaSize(memreq);
            size_t arenareq_with_reserve = AdjustArenaSizeForReserveRequirements(arenareq, rsrvfactor, rsrvminimum);
            size_t memreq_with_reserve = memreq + (arenareq_with_reserve - arenareq);

            if (oldarena != NULL) {
                oldarena->m_pArenaMemMgr->m_fnFree(pOldArenaBuffer, nOldArenaSize);
                oldarena = NULL;

                // Zero variables to avoid another freeing on exit
                pOldArenaBuffer = NULL;
                nOldArenaSize = 0;
            }

            // Allocate new arena
            void *pNewArenaBuffer = memmgr->m_fnAlloc(arenareq_with_reserve);
            if (pNewArenaBuffer == NULL) {
                break;
            }

            arena = (dxWorldProcessMemArena *)dEFFICIENT_PTR(pNewArenaBuffer);

            void *blockbegin = dEFFICIENT_PTR(arena + 1);
            void *blockend = dOFFSET_EFFICIENTLY(blockbegin, memreq_with_reserve);

            arena->m_pAllocBegin = blockbegin;
            arena->m_pAllocEnd = blockend;
            arena->m_pArenaBegin = pNewArenaBuffer;
            arena->m_pAllocCurrentOrNextArena = NULL;
            arena->m_pArenaMemMgr = memmgr;
        }

        allocsuccess = true;
    }
    while (false);

    if (!allocsuccess) {
        if (pOldArenaBuffer != NULL) {
            dIASSERT(oldarena != NULL);
            oldarena->m_pArenaMemMgr->m_fnFree(pOldArenaBuffer, nOldArenaSize);
        }
        arena = NULL;
    }

    return arena;
}

void dxWorldProcessMemArena::FreeMemArena (dxWorldProcessMemArena *arena)
{
    size_t memsize = arena->GetMemorySize();
    size_t arenasize = dxWorldProcessMemArena::MakeArenaSize(memsize);

    void *pArenaBegin = arena->m_pArenaBegin;
    arena->m_pArenaMemMgr->m_fnFree(pArenaBegin, arenasize);
}


size_t dxWorldProcessMemArena::AdjustArenaSizeForReserveRequirements(size_t arenareq, float rsrvfactor, unsigned rsrvminimum)
{
    float scaledarena = arenareq * rsrvfactor;
    size_t adjustedarena = (scaledarena < SIZE_MAX) ? (size_t)scaledarena : SIZE_MAX;
    size_t boundedarena = (adjustedarena > rsrvminimum) ? adjustedarena : (size_t)rsrvminimum;
    return dEFFICIENT_SIZE(boundedarena);
}


dxWorldProcessMemArena *dxAllocateTemporaryWorldProcessMemArena(
    size_t memreq, const dxWorldProcessMemoryManager *memmgr/*=NULL*/, const dxWorldProcessMemoryReserveInfo *reserveinfo/*=NULL*/)
{
    const dxWorldProcessMemoryManager *surememmgr = memmgr ? memmgr : &g_WorldProcessMallocMemoryManager;
    const dxWorldProcessMemoryReserveInfo *surereserveinfo = reserveinfo ? reserveinfo : &g_WorldProcessDefaultReserveInfo;
    dxWorldProcessMemArena *arena = dxWorldProcessMemArena::ReallocateMemArena(NULL, memreq, surememmgr, surereserveinfo->m_fReserveFactor, surereserveinfo->m_uiReserveMinimum);
    return arena;
}

void dxFreeTemporaryWorldProcessMemArena(dxWorldProcessMemArena *arena)
{
    dxWorldProcessMemArena::FreeMemArena(arena);
}

