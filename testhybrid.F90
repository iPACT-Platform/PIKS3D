PROGRAM MAIN
    USE MPI
    implicit none
    INTEGER :: MPI_ERR, proc, nproc
    call MPI_Init(MPI_ERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc, MPI_ERR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, MPI_ERR)
    !$OMP parallel
    PRINT*, "THIS IS"
    !$OMP end parallel
    CALL MPI_FINALIZE(MPI_ERR)
END
