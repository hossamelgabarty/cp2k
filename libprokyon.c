#include <stdio.h>

int prokyon_train_forces(double* coords, double* forces, double energy, int natoms);
int prokyon_atoms(char** labels, int natoms);

int prokyon_atoms(char** labels, int natoms)
{
      int i;
      printf("atoms labels:\n");
      for (i=0; i<natoms; i++)
      {
         printf("%s\n", labels[i]);
      }
      printf("\n");

}

int prokyon_train_forces(double* coords, double* forces, double energy, int natoms)
{
   printf("Hi from C!\n");
   printf("natoms is %d\n", natoms);
   printf("energy is %f\n", energy);
   printf("coord0 is %f\n", coords[0]);
   printf("force4 is %f\n", forces[4]);
   coords[0] = 0.0;

   return 0;
}
