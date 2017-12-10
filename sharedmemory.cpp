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

struct data {                       // container for an array that we're going to share
	double numbers[ARRAY_SIZE];
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
    * The child process will inherit both keys and links (obtained with semget and shmget)
    * to the semaphore as well as the shared memory.
    */
    pid_t pid = fork();

	if (pid < 0) {
		printf("Fork failed.\n");
        exit(EXIT_FAILURE);
	}

	else if (pid == 0) {                        // child process

        int status;

        status = semop(semId, &semBusy, 1);     // attempt to flip the semaphore to indicate that the process uses the shared memory segment
        if (status == -1) {                     // if semaphore is currently set to 0, wait until it's 1
            printf("Failed to lock the critical section.\n");
            exit(EXIT_FAILURE);
        }

        double sum = 0.0;

        for (int i = 0; i < ARRAY_SIZE; i++) {
            sum += sharedData->numbers[i];
        }

        printf("sum equals %f\n", sum);

        status = semop(semId, &semIdle, 1);     // flip the semaphore again to indicate that the memory is no longer in use
        if (status == -1) {
            printf("Failed to unlock the critical section.\n");
            exit(EXIT_FAILURE);
        }

		exit(EXIT_SUCCESS);
	}

	else {                                      // parent process

        int status;

        for (int i = 0; i < ARRAY_SIZE; i++) {
            sharedData->numbers[i] = 1.0;
        }

        status = semop(semId, &semIdle, 1);     // flip the semaphore to indicate that the memory is availible
        if (status == -1) {
            printf("Failed to unlock the critical section.\n");
            exit(EXIT_FAILURE);
        }

        sleep(2);

        status = semop(semId, &semBusy, 1);     // try to flip the semaphore to make sure the memory not in use
        if (status == -1) {
            printf("Failed to unlock the critical section.\n");
            exit(EXIT_FAILURE);
        }

        shmctl(shmId, IPC_RMID, NULL);

		exit(EXIT_SUCCESS);
	}


	return 0;
}
