!! @generated from magma2_zfortran.F90, fortran z -> c, Wed Nov  1 17:17:00 2023

module magma2_cfortran

use magma2_common
implicit none

!! =============================================================================
!! Fortran interfaces to C functions
interface

    !! -------------------------------------------------------------------------
    !! CPU interfaces (matrix in CPU memory)
    subroutine magma_cgetrf( m, n, A, lda, ipiv, info ) &
    bind(C, name="magma_cgetrf")
        use iso_c_binding
        integer(c_int),            value  :: m, n, lda
        complex(c_float_complex), target :: A(lda,*)
        integer(c_int),            target :: ipiv(*)
        integer(c_int),            target :: info  !! int*
    end subroutine

    subroutine magma_cpotrf( uplo, n, A, lda, info ) &
    bind(C, name="magma_cpotrf")
        use iso_c_binding
        integer(c_int),            value  :: uplo
        integer(c_int),            value  :: n, lda
        complex(c_float_complex), target :: A(lda,*)
        integer(c_int),            target :: info  !! int*
    end subroutine

    !! -------------------------------------------------------------------------
    !! GPU interfaces (matrix in GPU memory)
    subroutine magma_cgetrf_gpu( m, n, dA, lda, ipiv, info ) &
    bind(C, name="magma_cgetrf_gpu")
        use iso_c_binding
        integer(c_int), value  :: m, n, lda
        type(c_ptr),    value  :: dA
        integer(c_int), target :: ipiv(*)
        integer(c_int), target :: info  !! int*
    end subroutine

    subroutine magma_cpotrf_gpu( uplo, n, dA, lda, info ) &
    bind(C, name="magma_cpotrf_gpu")
        use iso_c_binding
        integer(c_int), value  :: uplo, n, lda
        type(c_ptr),    value  :: dA
        integer(c_int), target :: info  !! int*
    end subroutine

    subroutine magma_cgetri_gpu( n, dA, lda, ipiv, dwork, lwork, info ) &
        bind(C, name="magma_cgetri_gpu")
            use iso_c_binding
            integer(c_int), value  :: n
            type(c_ptr),    value  :: dA
            integer(c_int), value  :: lda
            integer(c_int), target :: ipiv(*)
            type(c_ptr),    value  :: dwork
            integer(c_int), value  :: lwork
            integer(c_int), target :: info  !! int*
    end subroutine

    function magma_get_cgetri_nb( n ) result(size) &
        bind(C, name="magma_get_cgetri_nb")
            use iso_c_binding
            integer(c_int) :: n
            integer(c_int) :: size
    end function 

    !! -------------------------------------------------------------------------
    !! batched GPU interfaces (all arrays in GPU memory)
    subroutine magma_cgetrf_batched( &
        m, n, dA_array, lda, ipiv_array, info_array, batchcount, queue ) &
    bind(C, name="magma_cgetrf_batched")
        use iso_c_binding
        integer(c_int), value  :: m, n, lda, batchcount
        type(c_ptr),    value  :: dA_array    !! double_complex**
        type(c_ptr),    value  :: ipiv_array  !! int**
        type(c_ptr),    value  :: info_array  !! int*
        type(c_ptr),    value  :: queue
    end subroutine

    !! -------------------------------------------------------------------------
    !! BLAS (matrices in GPU memory)
    subroutine magma_caxpy( &
        n, &
        alpha, dx, incx, &
               dy, incy, &
        queue ) &
    bind(C, name="magma_caxpy")
        use iso_c_binding
        integer(c_int),             value :: n, incx, incy
        complex(c_float_complex),  value :: alpha
        type(c_ptr),                value :: dx, dy
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    subroutine magma_cgemv( &
        transA, m, n, &
        alpha, dA, lda, &
               dx, incx, &
        beta,  dy, incy, &
        queue ) &
    bind(C, name="magma_cgemv")
        use iso_c_binding
        integer(c_int),             value :: transA, m, n, lda, incx, incy
        complex(c_float_complex),  value :: alpha, beta
        type(c_ptr),                value :: dA, dx, dy
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

    subroutine magma_cgemm( &
        transA, transB, m, n, k, &
        alpha, dA, lda, &
               dB, ldb, &
        beta,  dC, ldc, &
        queue ) &
    bind(C, name="magma_cgemm")
        use iso_c_binding
        integer(c_int),             value :: transA, transB, m, n, k, lda, ldb, ldc
        complex(c_float_complex),  value :: alpha, beta
        type(c_ptr),                value :: dA, dB, dC
        type(c_ptr),                value :: queue  !! queue_t
    end subroutine

end interface

!! =============================================================================
!! Fortran routines & functions
contains

    !! -------------------------------------------------------------------------
    !! malloc wrappers
    integer(c_int) function magma_cmalloc( ptr, n )
        use iso_c_binding
        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_cmalloc = magma_malloc( ptr, n*sizeof_complex )
    end function

    integer(c_int) function magma_cmalloc_cpu( ptr, n )
        use iso_c_binding
        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_cmalloc_cpu = magma_malloc_cpu( ptr, n*sizeof_complex )
    end function

    integer(c_int) function magma_cmalloc_pinned( ptr, n )
        use iso_c_binding
        type(c_ptr),       target :: ptr  !! void**
        integer(c_size_t), value  :: n
        
        magma_cmalloc_pinned = magma_malloc_pinned( ptr, n*sizeof_complex )
    end function

    !! -------------------------------------------------------------------------
    !! set/get wrappers
    subroutine magma_csetmatrix( &
        m, n, hA_src, lda, dB_dst, ldb, queue )
        use iso_c_binding
        integer(c_int),            value  :: m, n, lda, ldb
        complex(c_float_complex), target :: hA_src(lda,*)
        type(c_ptr),               value  :: dB_dst
        type(c_ptr),               value  :: queue
        
        call magma_setmatrix_internal( &
                m, n, int(sizeof_complex), c_loc(hA_src), lda, dB_dst, ldb, queue, &
                "magma_csetmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

    subroutine magma_cgetmatrix( &
        m, n, dA_src, lda, hB_dst, ldb, queue )
        use iso_c_binding
        integer(c_int),            value  :: m, n, lda, ldb
        type(c_ptr),               value  :: dA_src
        complex(c_float_complex), target :: hB_dst(ldb,*)
        type(c_ptr),               value  :: queue
        
        call magma_getmatrix_internal( &
                m, n, int(sizeof_complex), dA_src, lda, c_loc(hB_dst), ldb, queue, &
                "magma_cgetmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine

end module
