#include <iostream>
#include <stdio.h>
#include <sys/sem.h>
#include <unistd.h>
#include <cstdlib>
#include <sys/shm.h>
#include <string.h>
#include <semaphore.h>

using namespace std;

#define ARRAY_SIZE 100000           // size of the shared array
#define FORKS 4
#define ARRAY_SLICE 25000

struct data {                       // container for an array that we're going to share
	double numbers[ARRAY_SIZE];
    double sum;
};

struct sembuf semBusy = {0, -1, 0}; // sembuf structure for locking a semaphore
struct sembuf semIdle = {0, 1, 0};  // sembuf structure for unlocking a semaphore

int main() {

//******************************************************************************
    /*
    * Here we make a semaphore. Goes as follows:
    * 1. Generate a key using a path and a number
    * 2. Make said semaphore with semget
    * The semaphore is set to 0 ("resource unavailible") by default.
    */
    key_t semKey;
    int semId;

    semKey = ftok ("/dev/null", 777);

    semId = semget(semKey, 1, 0666 | IPC_CREAT);
    if (semId == -1) {
        printf("Failed to make a semaphore.\n");
        exit(EXIT_FAILURE);
    }

//******************************************************************************
    /*
    * Making a shared memory segment the same way as a semaphore.
    */
	key_t shmKey;
	int shmId;

	shmKey = ftok ("/dev/null", 888);
	shmId = shmget(shmKey, sizeof(struct data), 0666 | IPC_CREAT);
	if (shmId == -1) {
        printf("Failed to grab shared memory.\n");
        exit(EXIT_FAILURE);
    }

//******************************************************************************

    struct data* sharedData;                            // variable that's gonna be used as a linking point for shared memory

	sharedData = (struct data*)shmat(shmId, NULL, 0);   // linking the shared memory with shmat
	if (sharedData < 0) {
		printf("Failed to link variables to shared memory.\n");
		exit(EXIT_FAILURE);
	}

//******************************************************************************
/*
    int fd[2];

    if (pipe(fd)) {
        printf("Failed to allocate a pipe.\n");
        exit(EXIT_FAILURE);
    }
*/
//******************************************************************************
    
    int status;

    for (int i = 0; i < ARRAY_SIZE; i++) {
        sharedData->numbers[i] = 1.0;
    }
    sharedData->sum = 0.0;

    status = semop(semId, &semIdle, 1);     // flip the semaphore to indicate that the memory is availible
    if (status == -1) {
        printf("Failed to unlock the critical section.\n");
        exit(EXIT_FAILURE);
    }

    pid_t pid;

    for (int i = 0; i < FORKS; i++) {

        pid = fork();

        if (pid > (pid_t) 0) {

            /* This is the parent process. */
/*
            close (fd[0]);

            int bytes = write(fd[1], &i, sizeof(int));
            if (bytes < sizeof(int)) {
                printf("Pipe overflow during the step %d.\n", i);
                exit(EXIT_FAILURE);
            }
*/       
        }

        else if (pid < (pid_t) 0) {

            /* The fork failed. */
            printf("Fork failed during the step %d!\n", i);
            exit(EXIT_FAILURE);
        }
        else {

            /* This is the child process. */
/*      
            close (fd[1]);

            int forkNumber;
            read(fd[0], &forkNumber, sizeof(int));
*/
            double tmp = 0.0;

            for (int j = i * ARRAY_SLICE; j < (i + 1) * ARRAY_SLICE; j++) {
                tmp += sharedData->numbers[j];
            }

            status = semop(semId, &semBusy, 1);     // attempt to flip the semaphore to indicate that the process uses the shared memory segment
            if (status == -1) {
                printf("Failed to lock the critical section.\n");
                exit(EXIT_FAILURE);
            }

            sharedData->sum += tmp;
            
            status = semop(semId, &semIdle, 1);     // flip the semaphore to indicate that the memory is availible
            if (status == -1) {
                printf("Failed to unlock the critical section.\n");
                exit(EXIT_FAILURE);
            }

            exit(EXIT_SUCCESS);
        }
    }

    sleep(1);

    status = semop(semId, &semBusy, 1);     // try to flip the semaphore to make sure the memory not in use
    if (status == -1) {
        printf("Failed to lock the critical section.\n");
        exit(EXIT_FAILURE);
    }

    printf("Sum equals %f\n", sharedData->sum);

    shmctl(shmId, IPC_RMID, NULL);

	exit(EXIT_SUCCESS);
}
