#ifndef TSP_TSP_H
#define TSP_TSP_H
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif


/* dernier minimum trouv� */
extern int minimum;
int getMin();
void setMin(int);

int present (int city, int hops, tsp_path_t path, uint64_t vpres);
void tsp (int hops, int len, uint64_t vpres, tsp_path_t path, long long int *cuts, tsp_path_t sol, int *sol_len, int *min_local);

#ifdef __cplusplus
}
#endif


#endif
