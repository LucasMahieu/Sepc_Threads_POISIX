#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <sys/time.h>
#include <time.h>
#include <assert.h>
#include <complex.h>
#include <stdbool.h>
#include <unistd.h>
#include <pthread.h>
#include "tsp-types.h"
#include "tsp-job.h"
#include "tsp-genmap.h"
#include "tsp-print.h"
#include "tsp-tsp.h"
#include "tsp-lp.h"
#include "tsp-hkbound.h"


/* macro de mesure de temps, retourne une valeur en nanosecondes */
#define TIME_DIFF(t1, t2) \
	((t2.tv_sec - t1.tv_sec) * 1000000000ll + (long long int) (t2.tv_nsec - t1.tv_nsec))

/* tableau des distances */
tsp_distance_matrix_t tsp_distance ={};

/** Param�tres **/

/* nombre de villes */
int nb_towns=10;
/* graine */
long int myseed= 0;
/* nombre de threads */
int nb_threads=1;
/* Nombre de thread en court */
int current_nb_thread = 0;

/* affichage SVG */
bool affiche_sol= false;
bool affiche_progress=false;
bool quiet=false;

/* Structure � passer en argument du thread */
typedef struct args {
	struct tsp_queue* q;
	tsp_path_t solution;
	long long int *cuts;
	tsp_path_t sol;
	int* sol_len;
	pthread_mutex_t *mutex;
} args;

// mutex
extern pthread_mutex_t mutex_get_job;
extern pthread_mutex_t mutex_get_min;
extern pthread_mutex_t mutex_tsp;
extern pthread_mutex_t mutex_tsp2;


// Pour pouvoir avoir un �uivalant de clock_gettime sur mac.
// Pour mac activer le define 
// Pour linux commenter le define
//#define __MACH__
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif


void current_utc_time(struct timespec *ts);
void current_utc_time(struct timespec *ts) {

#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
	clock_serv_t cclock;
	mach_timespec_t mts;
	host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
	clock_get_time(cclock, &mts);
	mach_port_deallocate(mach_task_self(), cclock);
	ts->tv_sec = mts.tv_sec;
	ts->tv_nsec = mts.tv_nsec;
#else
	clock_gettime(CLOCK_REALTIME, ts);
#endif
}

static void generate_tsp_jobs (struct tsp_queue *q, int hops, int len, uint64_t vpres, tsp_path_t path, long long int *cuts, tsp_path_t sol, int *sol_len, int depth)
{
	if (len >= minimum) {
		(*cuts)++ ;
		return;
	}

	if (hops == depth) {
		/* On enregistre du travail � faire plus tard... */
		add_job (q, path, hops, len, vpres);
	} else {
		int me = path [hops - 1];        
		for (int i = 0; i < nb_towns; i++) {
			if (!present (i, hops, path, vpres)) {
				path[hops] = i;
				vpres |= (1<<i);
				int dist = tsp_distance[me][i];
				generate_tsp_jobs (q, hops + 1, len + dist, vpres, path, cuts, sol, sol_len, depth);
				vpres &= (~(1<<i));
			}
		}
	}
}

static void usage(const char *name) {
	fprintf (stderr, "Usage: %s [-s] <ncities> <seed> <nthreads>\n", name);
	exit (-1);
}
/////////////////////////////////////
//Fonction ex�cuter par les threads
void* work(void* arg){
    /*
	struct tsp_queue* q = *(args *)arg.q;
	tsp_path_t solution = *(args *)arg.solution;
	uint64_t* vpres = *(args *)arg.vpres;
	long long int* cuts = *(args *)arg.cuts;
	tsp_path_t sol = *(args *)arg.sol;
	int* sol_len = *(args *)arg.sol_len;
    */
	int nb_jobs_finished = 0;
    args A = *(args*)arg;
	int min_loc = getMin();
	tsp_path_t sol_loc;
	uint64_t vpres_loc=1;
	memcpy(sol_local, A.solution,nb_thread*sizeof(int));
	while (!empty_queue (&q)) {
	int hops = 0, len = 0;
	
	// DEBUT DE ZONE A DISTIBUER
	// ////////////////////////////
	// A faire en exlu, fait
	pthread_mutex_lock(&mutex_get_job);
	get_job (A.q, sol_loc, &hops, &len, &vpres_loc);
	pthread_mutex_unlock(&mutex_get_job);

	// le noeud est moins bon que la solution courante
	if (min < INT_MAX
		&& (nb_towns - hops) > 10
		&& ( (lower_bound_using_hk(sol_loc, hops, len, vpres)) >= min_loc ||
		(lower_bound_using_lp(sol_loc, hops, len, vpres_loc)) >= min_loc)
	){ return NULL; }

	// Si le noeud est meilleur que la sol courante 
	// Faire attention au ressource, cuts et sol et sol_len proteg�
	tsp(hops, len, vpres_loc, sol_loc, A.cuts, A.sol, A.sol_len, &min_loc);

	if( nb_jobs_finished == REFRESH_MIN_RATE ){
		nb_jobs_finished = 0;
		setMin(min_loc);
		int m getMin(); 
		if(m<min_loc) min_loc = m; 
	}
	nb_jobs_finished ++;
	}
	setMin(min_loc);
	return NULL;
// FIN DE ZONE A DISTIBUER 
///////////////////////////////
}


int main (int argc, char **argv)
{
	unsigned long long perf;
	tsp_path_t path;
	tsp_path_t sol;
	int sol_len;
	long long int cuts = 0;
	struct tsp_queue q;
	struct timespec t1, t2;
	// Compteur de boucle threads
	int i=0;

	/* lire les arguments */
	int opt;
	while ((opt = getopt(argc, argv, "spq")) != -1) {
		switch (opt) {
			case 's':
				affiche_sol = true;
				break;
			case 'p':
				affiche_progress = true;
				break;
			case 'q':
				quiet = true;
				break;
			default:
				usage(argv[0]);
				break;
		}
	}

	if (optind != argc-3)
		usage(argv[0]);

	nb_towns = atoi(argv[optind]);
	myseed = atol(argv[optind+1]);
	nb_threads = atoi(argv[optind+2]);
	assert(nb_towns > 0);
	assert(nb_threads > 0);

	// Creation des threads :
	pthread_t* thread_tid;
	void *status;
    /* C�ation et init des mutex*/
	pthread_mutex_t mutex_thread;
    if(pthread_mutex_init(&mutex_thread,NULL)){
        printf("errorinit  mutex_thread");
        return EXIT_FAILURE;
    }
    if(pthread_mutex_init(&mutex_get_job,NULL)){
        printf("error mutex init get_job");
        return EXIT_FAILURE;
    }

   if(pthread_mutex_init(&mutex_tsp,NULL)){
        printf("error mutex init tsp");
        return EXIT_FAILURE;
   }
    if(pthread_mutex_init(&mutex_get_min,NULL)){
        printf("error mutex init get_min");
        return EXIT_FAILURE;
    }
    if(pthread_mutex_init(&mutex_tsp2,NULL)){
        printf("error mutex init tsp2");
        return EXIT_FAILURE;
    }
	
	if(( thread_tid = (pthread_t*)calloc(nb_threads,sizeof(*thread_tid)) )==NULL){
		printf("calloc thread_tid error\n");
		return -1;
	}
	if(( status = (void*)calloc(nb_threads,sizeof(*status)) )==NULL){
		printf("calloc status error\n");
		return -1;
	}
	args* arguments;
	if(( arguments = (args*)calloc(nb_threads,sizeof(*arguments)) )==NULL){
		printf("calloc  args error\n");
		return -1;
	}
	minimum = INT_MAX;

	/* generer la carte et la matrice de distance */
	if (! quiet)
		fprintf (stderr, "ncities = %3d\n", nb_towns);
	genmap ();

	init_queue (&q);

	//clock_gettime(CLOCK_REALTIME, &t1);
	current_utc_time(&t1);

	memset (path, -1, MAX_TOWNS * sizeof (int));
	path[0] = 0;

	/* mettre les travaux dans la file d'attente */
	generate_tsp_jobs (&q, 1, 0, vpres, path, &cuts, sol, &sol_len, 3);
	no_more_jobs (&q);

	/* calculer chacun des travaux */
	//tsp-path_t est un tableau de distance
	tsp_path_t solution;
	memset (solution, -1, MAX_TOWNS * sizeof (int));
	solution[0] = 0;
	// DEBUT DE ZONE A DISTIBUER
	// ////////////////////////////
	for (i=0; i<nb_thread; i++){
		arguments[i].q = &q;
		memcpy(arguments[i].solution,solution, MAX_TOWNS*sizeof(int));
		arguments[i].cuts = &cuts;
		memcpy(arguments[i].sol, sol, MAX_TOWNS*sizeof(int));
		//arguments[i].sol_len = &sol_len;
		arguments[i].mutex = &mutex_thread;

		pthread_create(&thread_tid[i],NULL,work, &(arguments[i]) );

	}
	for(i=O; i<nb_thread; i++){
		pthread_join(thread_tid+i,status+i);
	}
	// FIN DE ZONE A DISTIBUER 
	///////////////////////////////

	//    clock_gettime(CLOCK_REALTIME, &t2);
	current_utc_time(&t2);

	sol_len = getMin();
	if (affiche_sol)
		print_solution_svg (sol, sol_len);

	perf = TIME_DIFF (t1,t2);
	printf("<!-- # = %d seed = %ld len = %d threads = %d time = %lld.%03lld ms ( %lld coupures ) -->\n",
			nb_towns, myseed, sol_len, nb_threads,
			perf/1000000ll, perf%1000000ll, cuts);

	free(thread_tid);
	free(arguments);
	return 0 ;
}
