/* Copyright (C) 2002-2012 The SSI Project. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   3. Neither the name of the project nor the names of its contributors 
      may be used to endorse or promote products derived from this software 
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/


#ifdef HAVE_CONFIG_H
	#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN32_H
	#include "lis_config_win32.h"
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
        #include <malloc.h>
#endif
#include <string.h>
#include <stdarg.h>
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_matrix_create_f
 * lis_matrix_duplicate_f
 * lis_matrix_destroy_f
 ************************************************/
#ifdef USE_FORTRAN

#undef __FUNC__
#define __FUNC__ "lis_matrix_create_f"
void lis_matrix_create_f(LIS_Comm_f *comm, LIS_MATRIX_F *Amat, LIS_INT *ierr)
{
	LIS_MATRIX			A;
	LIS_Comm			c_comm;

	LIS_DEBUG_FUNC_IN;

	#ifdef USE_MPI
		if( *comm==lis_comm_world_f )
		{
			c_comm = MPI_COMM_WORLD;
		}
		else
		{
			c_comm = MPI_Comm_f2c(*comm);
		}
	#else
		c_comm = *comm;
	#endif

	*ierr = lis_matrix_create(c_comm,&A);
	if( *ierr )	return;

	A->origin = LIS_ORIGIN_1;
	*Amat = LIS_P2V(A);
	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_size_f"
void lis_matrix_set_size_f(LIS_MATRIX_F *A, LIS_INT *local_n, LIS_INT *global_n, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_size((LIS_MATRIX)LIS_V2P(A),*local_n,*global_n);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_duplicate_f"
void lis_matrix_duplicate_f(LIS_MATRIX_F *Ain, LIS_MATRIX_F *Aout, LIS_INT *ierr)
{
	LIS_MATRIX	A;

	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_duplicate((LIS_MATRIX)LIS_V2P(Ain),&A);
	*Aout = LIS_P2V(A);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_destroy_f"
void lis_matrix_destroy_f(LIS_MATRIX_F *Amat, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_destroy((LIS_MATRIX)LIS_V2P(Amat));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_range_f"
void lis_matrix_get_range_f(LIS_MATRIX_F *A, LIS_INT *is, LIS_INT *ie, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_get_range((LIS_MATRIX)LIS_V2P(A),is,ie);
	(*is)++;
	(*ie)++;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_size_f"
void lis_matrix_get_size_f(LIS_MATRIX_F *A, LIS_INT *local_n, LIS_INT *global_n, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_get_size((LIS_MATRIX)LIS_V2P(A),local_n,global_n);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_value_f"
void lis_matrix_set_value_f(LIS_INT *flag, LIS_INT *i, LIS_INT *j, LIS_SCALAR *value, LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_value(*flag,*i,*j,*value,(LIS_MATRIX)LIS_V2P(A));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_type_f"
void lis_matrix_set_type_f(LIS_MATRIX_F *A, LIS_INT *matrix_type, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_type((LIS_MATRIX)LIS_V2P(A),*matrix_type);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_type_f"
void lis_matrix_get_type_f(LIS_MATRIX_F *A, LIS_INT *matrix_type, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_get_type((LIS_MATRIX)LIS_V2P(A),matrix_type);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_f"
void lis_matrix_malloc_f(LIS_MATRIX_F *A, LIS_INT *nnz_row, LIS_INT *nnz, LIS_INT *ierr)
{
	LIS_MATRIX AA;

	LIS_DEBUG_FUNC_IN;

	AA    = (LIS_MATRIX)LIS_V2P(A);
	*ierr = lis_matrix_malloc(AA,*nnz_row,nnz);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_crs_f"
void lis_matrix_malloc_crs_f(LIS_INT *n, LIS_INT *nnz, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_crs(*n, *nnz, &ptr, &index, &value);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_ccs_f"
void lis_matrix_malloc_ccs_f(LIS_INT *n, LIS_INT *nnz, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *ierr)

{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_ccs(*n, *nnz, &ptr, &index, &value);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_msr_f"
void lis_matrix_malloc_msr_f(LIS_INT *n, LIS_INT *nnz, LIS_INT *ndz, LIS_INT *index, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_msr(*n, *nnz, *ndz, &index, &value);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_dia_f"
void lis_matrix_malloc_dia_f(LIS_INT *n, LIS_INT *nnd, LIS_INT *index, LIS_SCALAR *value, LIS_INT *ierr)

{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_dia(*n, *nnd, &index, &value);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_ell_f"
void lis_matrix_malloc_ell_f(LIS_INT *n, LIS_INT *maxnzr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_ell(*n, *maxnzr, &index, &value);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_jds_f"
void lis_matrix_malloc_jds_f(LIS_INT *n, LIS_INT *nnz, LIS_INT *maxnzr, LIS_INT *perm, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_INT *ierr)

{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_jds(*n, *nnz, *maxnzr, &perm, &ptr, &index, &value);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_bsr_f"
void lis_matrix_malloc_bsr_f(LIS_INT *n, LIS_INT *bnr, LIS_INT *bnc, LIS_INT *bnnz, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_bsr(*n, *bnr, *bnc, *bnnz, &bptr, &bindex, &value); 

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_malloc_malloc_bsc_f"
void lis_matrix_malloc_bsc_f(LIS_INT *n, LIS_INT *bnr, LIS_INT *bnc, LIS_INT *bnnz, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_bsc(*n, *bnr, *bnc, *bnnz, &bptr, &bindex, &value);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_coo_f"
void lis_matrix_malloc_coo_f(LIS_INT *nnz, LIS_INT *row, LIS_INT *col, LIS_SCALAR *value, LIS_INT *ierr)

{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_coo(*nnz, &row, &col, &value);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_dns_f"
void lis_matrix_malloc_dns_f(LIS_INT *n, LIS_INT *np, LIS_SCALAR *value, LIS_INT *ierr)

{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_dns(*n, *np, &value);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_vbr_f"
void lis_matrix_malloc_vbr_f(LIS_INT *n, LIS_INT *nnz, LIS_INT *nr, LIS_INT *nc, LIS_INT *bnnz, LIS_INT *row, LIS_INT *col, LIS_INT *ptr, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_malloc_vbr(*n, *nnz, *nr, *nc, *bnnz, &row, &col, &ptr, &bptr, &bindex, &value);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_crs_f"
void lis_matrix_set_crs_f(LIS_INT *nnz, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_crs(*nnz,ptr,index,value,((LIS_MATRIX)LIS_V2P(A)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_ccs_f"
void lis_matrix_set_ccs_f(LIS_INT *nnz, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_ccs(*nnz,ptr,index,value,((LIS_MATRIX)LIS_V2P(A)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_msr_f"
void lis_matrix_set_msr_f(LIS_INT *nnz, LIS_INT *ndz, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_msr(*nnz,*ndz,index,value,((LIS_MATRIX)LIS_V2P(A)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_dia_f"
void lis_matrix_set_dia_f(LIS_INT *nnd, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_dia(*nnd,index,value,((LIS_MATRIX)LIS_V2P(A)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_ell_f"
void lis_matrix_set_ell_f(LIS_INT *maxnzr, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_ell(*maxnzr,index,value,((LIS_MATRIX)LIS_V2P(A)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_jds_f"
void lis_matrix_set_jds_f(LIS_INT *nnz, LIS_INT *maxnzr, LIS_INT *perm, LIS_INT *ptr, LIS_INT *index, LIS_SCALAR *value, LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_jds(*nnz,*maxnzr,perm,ptr,index,value,((LIS_MATRIX)LIS_V2P(A)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_bsr_f"
void lis_matrix_set_bsr_f(LIS_INT *bnr, LIS_INT *bnc, LIS_INT *bnnz, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_bsr(*bnr,*bnc,*bnnz,bptr,bindex,value,((LIS_MATRIX)LIS_V2P(A)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_bsc_f"
void lis_matrix_set_bsc_f(LIS_INT *bnr, LIS_INT *bnc, LIS_INT *bnnz, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_bsc(*bnr,*bnc,*bnnz,bptr,bindex,value,((LIS_MATRIX)LIS_V2P(A)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_coo_f"
void lis_matrix_set_coo_f(LIS_INT *nnz, LIS_INT *row, LIS_INT *col, LIS_SCALAR *value, LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_coo(*nnz,row,col,value,((LIS_MATRIX)LIS_V2P(A)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_dns_f"
void lis_matrix_set_dns_f(LIS_SCALAR *value, LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_dns(value,((LIS_MATRIX)LIS_V2P(A)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_vbr_f"
void lis_matrix_set_vbr_f(LIS_INT *nnz, LIS_INT *nr, LIS_INT *nc, LIS_INT *bnnz, LIS_INT *row, LIS_INT *col, LIS_INT *ptr, LIS_INT *bptr, LIS_INT *bindex, LIS_SCALAR *value, LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matrix_set_vbr(*nnz,*nr,*nc,*bnnz,row,col,ptr,bptr,bindex,value,((LIS_MATRIX)LIS_V2P(A)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_assemble_f"
void lis_matrix_assemble_f(LIS_MATRIX_F *A, LIS_INT *ierr)
{
	LIS_MATRIX AA;

	LIS_DEBUG_FUNC_IN;

	AA    = (LIS_MATRIX)LIS_V2P(A);
	*ierr = lis_matrix_assemble(AA);
	*A    = LIS_P2V(AA);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_f"
void lis_matrix_convert_f(LIS_MATRIX_F *Ain, LIS_MATRIX_F *Aout, LIS_INT *ierr)
{
	LIS_MATRIX AAin,AAout;

	LIS_DEBUG_FUNC_IN;

	AAin    = (LIS_MATRIX)LIS_V2P(Ain);
	AAout   = (LIS_MATRIX)LIS_V2P(Aout);
	*ierr   = lis_matrix_convert(AAin,AAout);
	*Aout   = LIS_P2V(AAout);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copy_f"
void lis_matrix_copy_f(LIS_MATRIX_F *Ain, LIS_MATRIX_F *Aout, LIS_INT *ierr)
{
	LIS_MATRIX AAin,AAout;

	LIS_DEBUG_FUNC_IN;

	AAin    = (LIS_MATRIX)LIS_V2P(Ain);
	AAout   = (LIS_MATRIX)LIS_V2P(Aout);
	*ierr   = lis_matrix_copy(AAin,AAout);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_diagonal_f"
void lis_matrix_get_diagonal_f(LIS_MATRIX_F *A, LIS_VECTOR_F *d, LIS_INT *ierr)
{
	LIS_MATRIX AA;
	LIS_VECTOR dd;

	LIS_DEBUG_FUNC_IN;

	AA    = (LIS_MATRIX)LIS_V2P(A);
	dd    = (LIS_VECTOR)LIS_V2P(d);
	*ierr = lis_matrix_get_diagonal(AA,dd);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scaling_f"
void lis_matrix_scaling_f(LIS_MATRIX_F *A, LIS_VECTOR_F *b, LIS_VECTOR_F *d, LIS_INT *action, LIS_INT *ierr)
{
	LIS_MATRIX AA;
	LIS_VECTOR dd,bb;

	LIS_DEBUG_FUNC_IN;

	AA    = (LIS_MATRIX)LIS_V2P(A);
	bb    = (LIS_VECTOR)LIS_V2P(b);
	dd    = (LIS_VECTOR)LIS_V2P(d);
	*ierr = lis_matrix_scaling(AA,bb,dd,*action);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_blocksize_f"
void lis_matrix_set_blocksize_f(LIS_MATRIX_F *A, LIS_INT *bnr, LIS_INT *bnc, LIS_INT *row, LIS_INT *col, LIS_INT *ierr)
{
	LIS_MATRIX AA;

	LIS_DEBUG_FUNC_IN;

	AA    = (LIS_MATRIX)LIS_V2P(A);
	if( *row==0 )
	{
		*ierr = lis_matrix_set_blocksize(AA,*bnr,*bnc,NULL,NULL);
	}
	else
	{
		*ierr = lis_matrix_set_blocksize(AA,*bnr,*bnc,row,col);
	}

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matvec_f"
void lis_matvec_f(LIS_MATRIX_F *A, LIS_VECTOR_F *x, LIS_VECTOR_F *y, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matvec((LIS_MATRIX)LIS_V2P(A),(LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_matvect_f"
void lis_matvect_f(LIS_MATRIX_F *A, LIS_VECTOR_F *x, LIS_VECTOR_F *y, LIS_INT *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_matvect((LIS_MATRIX)LIS_V2P(A),(LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#endif
