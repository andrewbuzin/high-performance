#include <iostream>
#include <cstdlib>
#include <sstream>
#include <stdio.h>
#include <pthread.h>
#include <unistd.h>
#include <queue>
#include <cmath>

#define THREADS 3
#define EPS 0.01

using namespace std;

pthread_t threads[THREADS];     // threads

pthread_cond_t cond;            // condition variable for the pool
pthread_mutex_t mutex;          // mutex lock for the pool

pthread_mutex_t idleMutex;      // mutex lock for IDLE_THREADS and RESULT variable
int IDLE_THREADS = THREADS;

struct task {                   // a structure that is being passed to the pool
    double a;
    double b;
    double integral;
};

queue<struct task> pool;                                    // working pool

void* integrate(void *_arg);                                // a thread

double rule(double a, double b);                            // integration rule

double F (double x);                                        // the function that we want to integrate

double A = 1, B = 10;                                       // the interval over which we're going to integrate the function

double RESULT = 0.0;                                        // the result of the integration

void addToPool (double _a, double _b, double _integral);    // add a new "task" into the pool

//************************************************************************************

int main() {

    addToPool(A, B, rule(A, B));            // initialize the pool with the first task

    pthread_cond_init(&cond, NULL);         // initialize the conditional variable and both locks
    pthread_mutex_init(&mutex, NULL);
    pthread_mutex_init(&idleMutex, NULL);

    /*
    * Start the threads
    */
    for(int i = 0; i < THREADS; i++) {
		int status = pthread_create(&threads[i], NULL, integrate, NULL);
		if(status != 0) {
			printf("Failed to start a thread at step %d.\n", i);
            exit(EXIT_FAILURE);
		}
	}

//************************************************************************************

    /*
    * Listen for all the threads to go idle after emptying the pool
    */
    while (true) {
        //sleep(1);

        pthread_mutex_lock(&idleMutex);

        if (IDLE_THREADS == THREADS && pool.empty()){
            printf("***\n");
            break;
        }

        pthread_mutex_unlock(&idleMutex);
    }

    /*
    * Stop all the threads
    */
    for(int i = 0; i < THREADS; i++) {
        pthread_cancel(threads[i]);
    }

    printf("Result: %f.\n", RESULT);    // print the result

    pthread_mutex_destroy(&mutex);      // clear both locks and the conitional variable
    pthread_cond_destroy(&cond);
    pthread_mutex_destroy(&idleMutex);

    exit(EXIT_SUCCESS);
}

//************************************************************************************

/*
* A function that we're going to integrate
*/
double F (double x) {
    return sin(x) * log(x);
}

/*
* Add a new task to the pool
* The task structure consists of the interval [a, b] and the integration result over that interval
*/
void addToPool (double _a, double _b, double _integral) {

    int status;
    task newTask;
    newTask.a = _a;
    newTask.b = _b;
    newTask.integral = _integral;

    /*
    * Locking the mutex, adding a task to a pool queue, then releasing the mutex
    */
    status = pthread_mutex_lock(&mutex);
    if(status != 0)
    {
        printf("Failed to lock the mutex at function addToPool().\n");
        exit(EXIT_FAILURE);
    }

    pool.push(newTask);

    status = pthread_mutex_unlock(&mutex);
    if(status != 0)
    {
        printf("Failed to unlock the mutex at function addToPool().\n");
        exit(EXIT_FAILURE);
    }

    pthread_cond_signal(&cond); // send a signal to idle threads
}

//************************************************************************************

void* integrate(void *_arg) {
    int status;

    while (true) {
        status = pthread_mutex_lock(&mutex);
        if(status != 0)
        {
            printf("Thread failed to lock the mutex.\n");
            exit(EXIT_FAILURE);
        }

        while (pool.empty()) {
            pthread_cond_wait(&cond, &mutex);                   // waiting for a new task to appear in the pool
        }

        printf("thread %ld going busy\n", pthread_self());

        pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);   // forbid terminating this thread

        pthread_mutex_lock(&idleMutex);                         // the thread is no longer idle
        IDLE_THREADS--;
        pthread_mutex_unlock(&idleMutex);

        task _task = pool.front();
        pool.pop();

        status = pthread_mutex_unlock(&mutex);                  // release the mutex and start processing the task
        if(status != 0)
        {
            printf("Thread failed to lock the mutex.\n");
            exit(EXIT_FAILURE);
        }

        double a = _task.a;
        double b = _task.b;                                     // parsing the task
        double integral = _task.integral;

        double c = (a + b) / 2;

        double newIntegral = rule(a, c) + rule(c, b);           // calculate the integral over [a, c] and [c, b]

        pthread_mutex_lock(&idleMutex);                         // lock the second mutex in order to safely work with IDLE_THREADS and RESULT

        if (abs(newIntegral - integral) < EPS) {
            RESULT += newIntegral;
        }
        else {
            addToPool(a, c, rule(a, c));                        // add new tasks to the pool
            addToPool(c, b, rule(c, b));
        }

        IDLE_THREADS++;                                         // the task is idle again
        printf("thread %ld going idle\n", pthread_self());
        pthread_mutex_unlock(&idleMutex);

        pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);    // allow terminating the thread
    }

    pthread_exit(NULL);
}

/*
* The trapezoid rule
*/
double rule(double a, double b) {
    //sleep(1);
    return (b - a) * ( (F(a) + F(b)) / 2 );
}
