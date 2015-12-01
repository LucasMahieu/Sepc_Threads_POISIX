/* impl�mentation des queues de jobs, nul besoin de lire dans un premier temps */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include "tsp-types.h"
#include "tsp-job.h"

struct tsp_job {
	tsp_path_t path;
	int hops;
	int len;
	uint64_t vpres;
};

struct tsp_cell {
	struct tsp_job tsp_job;
	struct tsp_cell *next;
};

void init_queue (struct tsp_queue *q) {
	q->first = 0;
	q->last = 0;
	q->end = 0;
	q->nbmax = 0;
	q->nb = 0;
}

int empty_queue (struct tsp_queue *q) {
	return ((q->first == 0) && (q->end == 1));
}

void add_job (struct tsp_queue *q, tsp_path_t p, int hops, int len, uint64_t vpres) {
	struct tsp_cell *ptr;

	ptr = malloc (sizeof (*ptr));
	if (!ptr) {
		printf ("L'allocation a echoue, recursion infinie ?\n");
		exit (1);
	}
	ptr->next = 0;
	ptr->tsp_job.len = len;
	ptr->tsp_job.hops = hops;
	ptr->tsp_job.vpres = vpres;
	memcpy (ptr->tsp_job.path, p, hops * sizeof(p[0]));

	if (q->first == 0) {
		q->first = q->last = ptr;
	} else {
		q->last->next = ptr;
		q->last = ptr;
	}
	q->nbmax ++;
	q->nb ++;
}


// Fonction de r�cup�ration de job dans la liste � faire en exclu
// Pas de soucis de perte de temps car aucun calcule n'est fait ici
// On r�cup�re juste les info sur un jobs et on d�cale le suivant
pthread_mutex_t mutex_get_job;
pthread_mutex_init(&mutex_get_job,NULL);
int get_job (struct tsp_queue *q, tsp_path_t p, int *hops, int *len, uint64_t *vpres) {

	pthread_mutex_lock(&mutex_get_job);
	struct tsp_cell *ptr;

	if (q->first == 0) {
		return 0;
	}

	ptr = q->first;

	q->first = ptr->next;
	if (q->first == 0) {
		q->last = 0;
	}

	*len = ptr->tsp_job.len;
	*hops = ptr->tsp_job.hops;
	*vpres = ptr->tsp_job.vpres;
	memcpy (p, ptr->tsp_job.path, *hops * sizeof(p[0]));

	free (ptr);

	q->nb --;
	if (affiche_progress)
		printf("<!- %d / %d %% ->\n",q->nb, q->nbmax);

	pthread_mutex_unlock(&mutex_get_job);
	return 1;
} 

void no_more_jobs (struct tsp_queue *q) {
	q->end = 1;
}











