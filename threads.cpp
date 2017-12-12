#include <iostream>
#include <cstdlib>
#include <sstream>
#include <stdio.h>
#include <pthread.h>
#include <unistd.h>

using namespace std;

void* factorial(void *_n) {     // function that is being passed to each new thread

	int n = * (int *)_n;
	int f = 1;

    for (int i = 1; i <= n; i++)
        f *= i;

    printf("i'm thread %d, my factorial is %d\n", n, f);

    pthread_exit(NULL);         // sending exit status back to the main thread
}

int main(int argc, char* argv[]) {

    if (argc != 2) {
        printf("Invalid number of arguments.\n");
        exit(EXIT_FAILURE);
    }

    istringstream ss(argv[1]);

    int n;
    if (!(ss >> n))
        printf("Invalid number argument.\n");

//*******************************************************************************

	pthread_t tid[n];
	int numbers[n];

	for(int i = 0; i < n; i++) {
		numbers[i] = i;

		/*
		* Making a thread that runs factorial(), and passing a pointer to numbers[i] to it.
		*/
		int status = pthread_create(&tid[i], NULL, factorial, (void*)&numbers[i]);
		if(status != 0) {
			printf("Failed to start a thread at step %d.\n", i);
            exit(EXIT_FAILURE);
		}

		/*
		* Waiting for the thread to finish and catching its exit status
		*/
		status = pthread_join(tid[i], NULL);

		if(status != 0)
		{
			printf("Join failed at step %d.\n", i);
			exit(EXIT_FAILURE);
		}
	}

    sleep(2);
    printf("All done.\n");
    return (EXIT_SUCCESS);

}
