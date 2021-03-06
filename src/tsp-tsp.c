#include <assert.h>
#include <string.h>
#include <stdint.h>
#include "tsp-types.h"
#include "tsp-genmap.h"
#include "tsp-print.h"
#include "tsp-tsp.h"
#include "tsp-lp.h"
#include "tsp-hkbound.h"
#include <pthread.h>

/* dernier minimum trouv� */
int minimum;
#define REFRESH_MIN_RATE 1
/////////////////////////////////////////
//A UTILISER QU'AVEC SET et GET
//ZONE DE PROTECTION DU MINIMUM
pthread_mutex_t mutex_get_min;
int getMin(){
	pthread_mutex_lock(&mutex_get_min);
	int m = minimum;
	pthread_mutex_unlock(&mutex_get_min);
	return m;
}
void setMin(int m){
	pthread_mutex_lock(&mutex_get_min);
	// Afin d'�tre sur que la valeur du min reste min 
	if(m < minimum){
		minimum = m;
	}
	pthread_mutex_unlock(&mutex_get_min);
}
//FIN DE ZONE DE PROTECTION DU MIN
///////////////////////////////////////////

/* r�solution du probl�me du voyageur de commerce */
int present (int city, int hops, tsp_path_t path, uint64_t vpres)
{
	(void) hops;
	(void) path;
	return (vpres & (1<<city)) != 0;
}

pthread_mutex_t mutex_tsp;
pthread_mutex_t mutex_tsp2;
void tsp (int hops, int len, uint64_t vpres, tsp_path_t path, long long int *cuts, tsp_path_t sol, int *sol_len,int* min_loc)
{
	// On fait la cpy locale du minimum 
	// et si on trouve mieux on en fera la cpy en exclu mutuelle 
	if (len + cutprefix[(nb_towns-hops)] >= *min_loc) {
		// S�curisation de cuts
		pthread_mutex_lock(&mutex_tsp);
		(*cuts)++ ;
		pthread_mutex_unlock(&mutex_tsp);
		return;
	}

	/* calcul de l'arbre couvrant comme borne inf�rieure */
	if ((nb_towns - hops) > 6 &&
			lower_bound_using_hk(path, hops, len, vpres) >= *min_loc) {
		// S�curisation de cuts
		pthread_mutex_lock(&mutex_tsp);
		(*cuts)++ ;
		pthread_mutex_unlock(&mutex_tsp);
		return;
	}


	/* un rayon de coupure � 15, pour ne pas lancer la programmation
	   lin�aire pour les petits arbres, plus rapide � calculer sans */
	if ((nb_towns - hops) > 22
			&& lower_bound_using_lp(path, hops, len, vpres) >= *min_loc) {
		// S�curisation de cuts
		pthread_mutex_lock(&mutex_tsp);
		(*cuts)++ ;
		pthread_mutex_unlock(&mutex_tsp);
		return;
	}

	// Cas ou c'est une solution retenue
	if (hops == nb_towns) {
		int me = path [hops - 1];
		int dist = tsp_distance[me][0]; // retourner en 0
		if (len + dist < *min_loc){
			//Cpy du nouveau min dans la variable locale min_loc
		// S�curisation de sol_len(pas besoin) et de sol(besoin)
			pthread_mutex_lock(&mutex_tsp2);
			*sol_len = len + dist;
			*min_loc = len + dist; //minimum = len + dist;
			memcpy(sol, path, nb_towns*sizeof(int));
			if (!quiet){
				print_solution (path, len+dist);
            }
			pthread_mutex_unlock(&mutex_tsp2);
		}
	} 
	else {
		int me = path[hops - 1];        
		for (int i = 0; i < nb_towns; i++) {
			if (!present (i, hops, path, vpres)) {
				path[hops] = i;
				vpres |= (1<<i);
				int dist = tsp_distance[me][i];
				tsp (hops + 1, len + dist, vpres, path, cuts, sol, sol_len, min_loc);
				vpres &= (~(1<<i));
			}
		}
	}
}

