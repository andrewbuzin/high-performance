#include <iostream>
#include <sstream>
#include <stdio.h>
#include <unistd.h>
#include  <sys/types.h>
#include <sys/wait.h>
#include <cstdlib>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 2) {
        printf("Invalid number of arguments.\n");
        exit(EXIT_FAILURE);
    }

    istringstream ss(argv[1]);

    int n;
    if (!(ss >> n))
        printf("Invalid number argument.\n");

    int fd[2];
    pid_t pid;

    if (pipe(fd)) {
        printf("Failed to allocate a pipe.\n");
        exit(EXIT_FAILURE);
    }

//*****************************************************************

    for (int i = 0; i < n; i++) {

        pid = fork();

        if (pid > (pid_t) 0) {

            /* This is the parent process.
                Close other end first. */

            close (fd[0]);

            pid_t thisPid = getpid();

            int bytes = write(fd[1], &thisPid, sizeof(thisPid));
            if (bytes < sizeof(thisPid)) {
                printf("Pipe overflow during the step %d.\n", i);
                exit(EXIT_FAILURE);
            }

            printf("dad sent %d\n", thisPid);
        }

        else if (pid < (pid_t) 0) {

            /* The fork failed. */
            printf("Fork failed during the step %d!\n", i);
            exit(EXIT_FAILURE);
        }
        else {

            /* This is the child process.
              Close other end first. */

            close (fd[1]);

            pid_t receivedPid;
            read(fd[0], &receivedPid, sizeof(receivedPid));

            pid_t parentPid = getppid();

            if (receivedPid == parentPid) {
                printf("i'm kid number %d, all good\n", i);
            }
            else {
                printf("i'm kid number %d, my dad's %d, what i read is %d\n", i, parentPid, receivedPid);
            }

            exit(EXIT_SUCCESS);
        }
    }

    sleep(1);
    return 0;
}
