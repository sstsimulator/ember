//
// Created by Teranishi, Keita on 10/19/21.
//
#include<stdio.h>
#include<mpi.h>

int main(int argc, char ** argv)
{
    MPI_Init(&argc,&argv);
    MPI_Finalize();
}
